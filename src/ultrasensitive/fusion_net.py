"""
Multi-Modal Fusion Network -- combines DNA + Methylation + Protein + Clinical signals.

Ported from ``cancer-screening/src/ml_models/fusion_net.py`` (SHA 2fe61590).

Uses cross-modal attention to learn interactions between different data modalities,
producing a unified cancer probability score and tissue-of-origin prediction.

NOTES:
  - This is a 4-modality FusionNet. DeepCatch's existing
    ``src/multimodal_fusion/advanced_fusion.py`` already implements
    ``CrossAttentionFusion`` + ``EarlyLateFusion`` + ``GCNTissueOfOrigin``
    in PyTorch (post-biomedical-review fixes).
  - Recommended migration path: use ``deepcatch/src/multimodal_fusion/`` as the
    production-grade fusion layer; this port remains as a documented
    reference implementation for the 4-modality (genomic+methyl+protein+clinical)
    pathway that ``advanced_fusion.py`` does not yet cover.
  - Deprecation notice: see ``git log`` -- the upstream impl is kept here for
    binary reproducibility with the CANCER-SEEK manuscript; after v2.4.0 it will
    route through ``deepcatch/src/multimodal_fusion/``.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class MultiModalDataset(Dataset):
    def __init__(
        self,
        genomic,
        methylation,
        protein,
        clinical,
        labels=None,
        tissues=None,
    ):
        self.genomic = torch.FloatTensor(genomic)
        self.methylation = torch.FloatTensor(methylation)
        self.protein = torch.FloatTensor(protein)
        self.clinical = torch.FloatTensor(clinical)
        self.labels = torch.FloatTensor(labels) if labels is not None else None
        self.tissues = torch.LongTensor(tissues) if tissues is not None else None

    def __len__(self):
        return len(self.genomic)

    def __getitem__(self, idx):
        item = {
            'genomic': self.genomic[idx],
            'methylation': self.methylation[idx],
            'protein': self.protein[idx],
            'clinical': self.clinical[idx],
        }
        if self.labels is not None:
            item['label'] = self.labels[idx]
        if self.tissues is not None:
            item['tissue'] = self.tissues[idx]
        return item


class CrossModalAttention(nn.Module):
    """Cross-attention layer: each modality attends to all other modalities."""

    def __init__(self, d_model: int, n_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        self.cross_attn = nn.MultiheadAttention(
            embed_dim=d_model,
            num_heads=n_heads,
            dropout=dropout,
            batch_first=True,
        )
        self.norm = nn.LayerNorm(d_model)
        self.ffn = nn.Sequential(
            nn.Linear(d_model, d_model * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 2, d_model),
        )

    def forward(self, x1: torch.Tensor, x2: torch.Tensor) -> torch.Tensor:
        attn_out, _ = self.cross_attn(x1, x2, x2)
        x = self.norm(x1 + attn_out)
        x = self.norm(x + self.ffn(x))
        return x


class ModalityEncoder(nn.Module):
    """Encode each modality to a fixed-size latent vector."""

    def __init__(self, input_dim: int, latent_dim: int = 128, dropout: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, latent_dim * 2),
            nn.LayerNorm(latent_dim * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(latent_dim * 2, latent_dim),
            nn.LayerNorm(latent_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class FusionNet(nn.Module):
    """
    Multi-modal fusion network with cross-attention (4-modality: g+m+p+c).

    Routes through ``deepcatch/src/multimodal_fusion/`` for future work;
    see module docstring for deprecation timeline.
    """

    def __init__(
        self,
        genomic_dim: int = 200,
        methylation_dim: int = 500,
        protein_dim: int = 20,
        clinical_dim: int = 16,
        latent_dim: int = 128,
        n_tissues: int = 25,
        dropout: float = 0.1,
    ):
        super().__init__()

        self.genomic_encoder = ModalityEncoder(genomic_dim, latent_dim, dropout)
        self.methylation_encoder = ModalityEncoder(methylation_dim, latent_dim, dropout)
        self.protein_encoder = ModalityEncoder(protein_dim, latent_dim, dropout)
        self.clinical_encoder = ModalityEncoder(clinical_dim, latent_dim, dropout)

        self.cross_attn_gm = CrossModalAttention(latent_dim, dropout=dropout)
        self.cross_attn_gp = CrossModalAttention(latent_dim, dropout=dropout)
        self.cross_attn_gc = CrossModalAttention(latent_dim, dropout=dropout)
        self.cross_attn_mp = CrossModalAttention(latent_dim, dropout=dropout)

        self.modality_gate = nn.Sequential(
            nn.Linear(latent_dim * 4, latent_dim * 4),
            nn.Sigmoid(),
        )

        self.shared = nn.Sequential(
            nn.Linear(latent_dim * 4, latent_dim * 2),
            nn.LayerNorm(latent_dim * 2),
            nn.GELU(),
            nn.Dropout(dropout),
        )

        self.cancer_head = nn.Linear(latent_dim * 2, 1)
        self.tissue_head = nn.Linear(latent_dim * 2, n_tissues)
        self.tf_head = nn.Sequential(
            nn.Linear(latent_dim * 2, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, 1),
            nn.Sigmoid(),
        )

    def forward(self, genomic, methylation, protein, clinical) -> dict:
        g = self.genomic_encoder(genomic)
        m = self.methylation_encoder(methylation)
        p = self.protein_encoder(protein)
        c = self.clinical_encoder(clinical)

        g = self.cross_attn_gm(g.unsqueeze(1), m.unsqueeze(1)).squeeze(1)
        g = self.cross_attn_gp(g.unsqueeze(1), p.unsqueeze(1)).squeeze(1)
        g = self.cross_attn_gc(g.unsqueeze(1), c.unsqueeze(1)).squeeze(1)

        all_z = torch.cat([g, m, p, c], dim=-1)
        gates = self.modality_gate(all_z)

        fused = torch.cat([g, m, p, c], dim=-1) * gates
        shared = self.shared(fused)

        return {
            'cancer_prob': torch.sigmoid(self.cancer_head(shared).squeeze(-1)),
            'tissue_logits': self.tissue_head(shared),
            'tumor_fraction': self.tf_head(shared).squeeze(-1),
        }


class FusionNetTrainer:
    """Minimal training loop for the 4-modality FusionNet."""

    def __init__(self, model: FusionNet, device: str = 'cuda', lr: float = 1e-4):
        self.model = model.to(device)
        self.device = device
        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=lr, weight_decay=1e-5,
        )
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            model.parameters(), T_max=100,
        )
        self.cancer_criterion = nn.BCEWithLogitsLoss(
            pos_weight=torch.tensor([50.0]).to(device),
        )
        self.tissue_criterion = nn.CrossEntropyLoss(label_smoothing=0.05)

    def train_step(self, batch) -> dict:
        self.model.train()
        self.optimizer.zero_grad()

        out = self.model(
            batch['genomic'].to(self.device),
            batch['methylation'].to(self.device),
            batch['protein'].to(self.device),
            batch['clinical'].to(self.device),
        )

        cancer_loss = self.cancer_criterion(
            out['cancer_prob'], batch['label'].to(self.device),
        )
        tissue_loss = self.tissue_criterion(
            out['tissue_logits'], batch['tissue'].to(self.device),
        )

        loss = 1.0 * cancer_loss + 0.3 * tissue_loss
        loss.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
        self.optimizer.step()

        return {
            'loss': loss.item(),
            'cancer_loss': cancer_loss.item(),
            'tissue_loss': tissue_loss.item(),
        }
