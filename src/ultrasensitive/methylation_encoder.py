"""
Deep Methylation Encoder -- Transformer-based methylation pattern analysis.

Ported from ``cancer-screening/src/ml_models/methylation_encoder.py`` (SHA ac70301).

Learns tissue-of-origin and cancer detection from 450K/850K CpG site beta values.
Uses self-attention to capture long-range methylation co-dependencies.

Architecture: 6-layer Transformer Encoder
  - Pre-training: Masked Methylation Prediction (BERT-style)
  - Fine-tuning: Cancer classification + tissue of origin + tumor fraction regression

NOTE: Reference architecture; integrates into ``deepcatch/src/foundation/``.
Phase 4 of ``deepcatch/METHYLATION_PROJECT.md`` is the validation gate.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pandas as pd
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class MethylationDataset(Dataset):
    """Dataset for methylation beta values."""

    def __init__(
        self,
        beta_values: np.ndarray,
        labels: Optional[np.ndarray] = None,
        tissue_labels: Optional[np.ndarray] = None,
    ):
        self.beta = torch.FloatTensor(beta_values)
        self.labels = torch.FloatTensor(labels) if labels is not None else None
        self.tissue_labels = torch.LongTensor(tissue_labels) if tissue_labels is not None else None

    def __len__(self) -> int:
        return len(self.beta)

    def __getitem__(self, idx: int):
        item = {'beta': self.beta[idx]}
        if self.labels is not None:
            item['label'] = self.labels[idx]
        if self.tissue_labels is not None:
            item['tissue'] = self.tissue_labels[idx]
        return item


class PositionalEncoding(nn.Module):
    """Sinusoidal positional encoding for genomic position awareness."""

    def __init__(self, d_model: int, max_len: int = 5000, dropout: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-np.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.pe[:x.size(1)]
        return self.dropout(x)


class DeepMethylationEncoder(nn.Module):
    """
    Transformer-based encoder for methylation data.

    Input: (batch, n_cpgs) beta values -> embedded -> Transformer -> predictions

    Heads:
      1. Cancer classification (binary)
      2. Tissue of origin (multi-class, 25 cancer types)
      3. Tumor fraction regression (0-100%)
    """

    def __init__(
        self,
        n_cpgs: int = 10000,
        d_model: int = 256,
        n_heads: int = 8,
        n_layers: int = 6,
        d_ff: int = 512,
        dropout: float = 0.1,
        n_tissues: int = 25,
        masked_prediction: bool = True,
    ):
        super().__init__()
        self.masked_prediction = masked_prediction

        self.beta_embedding = nn.Sequential(
            nn.Linear(1, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
            nn.Dropout(dropout),
        )

        self.pos_encoding = PositionalEncoding(d_model, max_len=n_cpgs + 1, dropout=dropout)
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model))

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_ff,
            dropout=dropout,
            activation='gelu',
            batch_first=True,
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        self.cancer_head = nn.Sequential(
            nn.Linear(d_model, d_model // 2),
            nn.LayerNorm(d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, 1),
        )

        self.tissue_head = nn.Sequential(
            nn.Linear(d_model, d_model // 2),
            nn.LayerNorm(d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, n_tissues),
        )

        self.tf_head = nn.Sequential(
            nn.Linear(d_model, d_model // 2),
            nn.LayerNorm(d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, 1),
            nn.Sigmoid(),
        )

        if masked_prediction:
            self.mask_head = nn.Sequential(
                nn.Linear(d_model, d_model // 2),
                nn.GELU(),
                nn.Linear(d_model // 2, 1),
                nn.Sigmoid(),
            )

        self.d_model = d_model

    def forward(
        self,
        beta: torch.Tensor,
        mask_pct: float = 0.0,
    ) -> dict:
        batch_size, n_cpgs = beta.shape

        x = beta.unsqueeze(-1)
        x = self.beta_embedding(x)

        cls_tokens = self.cls_token.expand(batch_size, -1, -1)
        x = torch.cat([cls_tokens, x], dim=1)

        x = self.pos_encoding(x)

        if self.masked_prediction and mask_pct > 0 and self.training:
            mask = torch.rand_like(beta) < mask_pct
            masked_beta = beta.clone()
            masked_beta[mask] = 0.5
            x_masked = torch.cat([cls_tokens, masked_beta.unsqueeze(-1)], dim=1)
            x_masked = self.beta_embedding(x_masked)
            x_masked = self.pos_encoding(x_masked)
            mask_pred = self.transformer(x_masked)
        else:
            mask_pred = None

        encoded = self.transformer(x)
        cls_repr = encoded[:, 0]

        cancer_logits = self.cancer_head(cls_repr).squeeze(-1)
        cancer_prob = torch.sigmoid(cancer_logits)

        tissue_logits = self.tissue_head(cls_repr)
        tumor_fraction = self.tf_head(cls_repr).squeeze(-1)

        return {
            'cancer_prob': cancer_prob,
            'cancer_logits': cancer_logits,
            'tissue_logits': tissue_logits,
            'tumor_fraction': tumor_fraction,
            'mask_pred': mask_pred,
            'cls_repr': cls_repr,
        }


class MethylationTrainer:
    """Training loop for the Deep Methylation Encoder."""

    def __init__(
        self,
        model: DeepMethylationEncoder,
        lr: float = 1e-4,
        weight_decay: float = 1e-5,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu',
    ):
        self.model = model.to(device)
        self.device = device
        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=lr, weight_decay=weight_decay
        )
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.model.parameters(), T_max=100
        )

        self.cancer_criterion = nn.BCEWithLogitsLoss(
            pos_weight=torch.tensor([50.0]).to(device)
        )
        self.tissue_criterion = nn.CrossEntropyLoss(label_smoothing=0.05)
        self.mse_criterion = nn.MSELoss()

        logger.info(f"MethylationTrainer initialized on device: {device}")

    def train_epoch(self, dataloader: DataLoader, mask_pct: float = 0.15) -> dict:
        self.model.train()
        total_loss = 0.0
        n_batches = 0

        for batch in dataloader:
            beta = batch['beta'].to(self.device)

            self.optimizer.zero_grad()
            outputs = self.model(beta, mask_pct=mask_pct)

            cancer_loss = self.cancer_criterion(
                outputs['cancer_logits'],
                batch['label'].to(self.device)
            )
            tissue_loss = self.tissue_criterion(
                outputs['tissue_logits'],
                batch['tissue'].to(self.device)
            )
            tf_loss = self.mse_criterion(
                outputs['tumor_fraction'],
                batch.get(
                    'tumor_fraction',
                    torch.zeros_like(outputs['tumor_fraction']),
                ).to(self.device)
            )

            mask_loss = 0.0
            if outputs['mask_pred'] is not None:
                mask_loss = F.binary_cross_entropy(
                    outputs['mask_pred'][:, 1:],
                    beta[:, :outputs['mask_pred'].shape[1]-1],
                )

            loss = (
                1.0 * cancer_loss +
                0.3 * tissue_loss +
                0.1 * tf_loss +
                0.1 * mask_loss
            )

            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        self.scheduler.step()
        return {'loss': total_loss / n_batches}

    def evaluate(self, dataloader: DataLoader) -> dict:
        self.model.eval()
        all_preds = []
        all_labels = []

        with torch.no_grad():
            for batch in dataloader:
                beta = batch['beta'].to(self.device)
                outputs = self.model(beta)
                probs = outputs['cancer_prob'].cpu().numpy()
                all_preds.extend(probs)
                if 'label' in batch:
                    all_labels.extend(batch['label'].numpy())

        from sklearn.metrics import roc_auc_score, average_precision_score
        if len(all_labels) > 0 and len(set(all_labels)) > 1:
            auc = float(roc_auc_score(all_labels, all_preds))
            ap = float(average_precision_score(all_labels, all_preds))
        else:
            auc, ap = 0.0, 0.0

        return {'auc_roc': auc, 'auc_pr': ap}

    def save(self, path: str) -> None:
        torch.save({
            'model_state': self.model.state_dict(),
            'optimizer_state': self.optimizer.state_dict(),
            'scheduler_state': self.scheduler.state_dict(),
        }, path)
        logger.info(f"Checkpoint saved to {path}")
