#!/usr/bin/env python3
"""
Foundation Model Architecture
==============================

MultiModalEncoder — Joint embedding via per-modality linear projections
and a 4-layer TransformerEncoder (~3-5M params).

PretrainHead — Masked modality prediction (reconstruct masked modalities
from context).

ContrastiveHead — Cross-modal contrastive loss (pull same-sample
modalities together, push different samples apart).
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import FoundationConfig, MODALITY_DIMS, MODALITY_NAMES


class LinearProjection(nn.Module):
    """Project a single modality's features into the joint embedding space."""

    def __init__(self, input_dim: int, embed_dim: int, dropout: float = 0.1):
        super().__init__()
        self.proj = nn.Linear(input_dim, embed_dim)
        self.norm = nn.LayerNorm(embed_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.norm(self.proj(x)))


class MultiModalEncoder(nn.Module):
    """
    Multi-modal joint encoder with per-modality projections and a
    shared TransformerEncoder backbone.

    Architecture:
        Input: Dict[str, Tensor] with 6 modality keys
        → LinearProjection (each modality → embed_dim)
        → [CLS] token prepended
        → 4-layer TransformerEncoder
        → Output: (batch, n_modalities, embed_dim) joint embedding

    Parameters
    ----------
    config : FoundationConfig
        Model configuration.
    modality_dims : dict, optional
        Override modality dimensions (default uses MODALITY_DIMS).
    """

    def __init__(
        self,
        config: FoundationConfig,
        modality_dims: Optional[Dict[str, int]] = None,
    ):
        super().__init__()
        self.config = config
        self.embed_dim = config.embed_dim

        # Determine modality dimensions
        self.modality_dims = modality_dims or MODALITY_DIMS
        self.modality_names = list(self.modality_dims.keys())
        self.n_modalities = len(self.modality_names)

        # Per-modality projections
        self.projections = nn.ModuleDict({
            name: LinearProjection(dim, self.embed_dim, config.dropout)
            for name, dim in self.modality_dims.items()
        })

        # Modality type embeddings (learned)
        self.modality_embed = nn.Parameter(
            torch.randn(self.n_modalities, self.embed_dim) * 0.02
        )

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=self.embed_dim,
            nhead=config.n_heads,
            dim_feedforward=config.ff_dim,
            dropout=config.dropout,
            activation="gelu",
            batch_first=True,
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer,
            num_layers=config.n_layers,
        )

        # Output normalization
        self.out_norm = nn.LayerNorm(self.embed_dim)

        self._init_weights()

    def _init_weights(self):
        """Initialize weights with small random values for stable training."""
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p, gain=0.5)

    def forward(
        self,
        modalities: Dict[str, torch.Tensor],
        mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Forward pass producing joint embeddings for all modalities.

        Parameters
        ----------
        modalities : dict[str, Tensor]
            Keys: modality names. Values: (batch, dim_i) tensors.
        mask : (batch, n_modalities) bool tensor, optional
            True = masked (not used as input). Masked modalities
            get zero input.

        Returns
        -------
        joint_embedding : (batch, n_modalities, embed_dim)
            Joint embedding for each modality per sample.
        """
        batch_size = None
        projected = []

        for i, name in enumerate(self.modality_names):
            if name in modalities:
                x = modalities[name]
                if x.dim() == 1:
                    x = x.unsqueeze(0)
                if batch_size is None:
                    batch_size = x.shape[0]
            else:
                # Missing modality — use zeros
                dim = self.modality_dims[name]
                if batch_size is None:
                    batch_size = 1
                x = torch.zeros(batch_size, dim, device=self.modality_embed.device)

            # Ensure float32
            if x.dtype != torch.float32:
                x = x.to(torch.float32)

            # Project
            emb = self.projections[name](x)

            # Apply mask if provided
            if mask is not None:
                is_masked = mask[:, i].unsqueeze(-1).float()  # (B, 1)
                emb = emb * (1.0 - is_masked)

            projected.append(emb)

        # Stack: (B, n_modalities, embed_dim)
        projected = torch.stack(projected, dim=1)

        # Add modality type embeddings
        mod_emb = self.modality_embed.unsqueeze(0).expand(batch_size, -1, -1)
        projected = projected + mod_emb

        # Transformer encoding (modalities attend to each other)
        encoded = self.transformer(projected)

        # Normalize output
        encoded = self.out_norm(encoded)

        return encoded

    def encode_single(
        self,
        modalities: Dict[str, torch.Tensor],
    ) -> torch.Tensor:
        """
        Encode modalities and return mean-pooled global representation.

        Parameters
        ----------
        modalities : dict[str, Tensor]

        Returns
        -------
        global_embedding : (batch, embed_dim)
            Mean-pooled across modalities.
        """
        joint = self.forward(modalities)  # (B, n_modalities, embed_dim)
        return joint.mean(dim=1)

    @property
    def num_params(self) -> int:
        """Total number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def estimate_params(self) -> int:
        """Estimate parameter count (same as num_params)."""
        return self.num_params


class PretrainHead(nn.Module):
    """
    Masked Modality Prediction Head.

    Given the joint embedding of unmasked modalities, predict the
    features of masked modalities.

    Architecture:
        joint_embedding[B, unmasked] → mean-pool → MLP → reconstructed modality features
    """

    def __init__(self, config: FoundationConfig):
        super().__init__()
        self.config = config
        self.embed_dim = config.embed_dim
        self.modality_dims = MODALITY_DIMS

        # Decoder: reconstruct each modality from joint embedding
        self.decoders = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(self.embed_dim, self.embed_dim),
                nn.GELU(),
                nn.Dropout(config.dropout),
                nn.Linear(self.embed_dim, dim),
            )
            for name, dim in self.modality_dims.items()
        })

    def forward(
        self,
        joint_embedding: torch.Tensor,
        mask: torch.Tensor,
    ) -> Dict[str, torch.Tensor]:
        """
        Predict masked modality features from joint embeddings.

        Parameters
        ----------
        joint_embedding : (batch, n_modalities, embed_dim)
            Joint embeddings from MultiModalEncoder.
        mask : (batch, n_modalities) bool tensor
            True = masked modality.

        Returns
        -------
        reconstructed : dict[str, Tensor]
            Reconstructed features for each modality.
        """
        batch_size = joint_embedding.shape[0]
        recon = {}

        # For each masked position, use context (mean of unmasked modalities)
        # to predict the masked modality's features
        for i, name in enumerate(self.modality_dims.keys()):
            is_masked = mask[:, i]  # (B,)

            # Context: mean of all modality embeddings (both masked/unmasked)
            # The encoder has already zeroed out masked inputs before transformer
            context = joint_embedding.mean(dim=1)  # (B, embed_dim)

            # Decode
            predicted = self.decoders[name](context)
            recon[name] = predicted

        return recon

    def compute_loss(
        self,
        reconstructed: Dict[str, torch.Tensor],
        modalities: Dict[str, torch.Tensor],
        mask: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute masked reconstruction loss (MSE on masked modalities only).

        Parameters
        ----------
        reconstructed : dict[str, Tensor]
            Predicted features.
        modalities : dict[str, Tensor]
            Ground truth features.
        mask : (batch, n_modalities) bool tensor

        Returns
        -------
        loss : scalar tensor
            Mean squared error.
        """
        total_loss = 0.0
        n_masked = 0

        for i, name in enumerate(self.modality_dims.keys()):
            if name not in modalities:
                continue
            target = modalities[name]
            if target.dim() == 1:
                target = target.unsqueeze(0)
            target = target.to(torch.float32)

            # Only compute loss on masked positions
            is_masked = mask[:, i]
            if is_masked.any():
                pred = reconstructed[name][is_masked]
                tgt = target[is_masked]
                total_loss += F.mse_loss(pred, tgt)
                n_masked += 1

        if n_masked == 0:
            return torch.tensor(0.0, device=mask.device)

        return total_loss / n_masked


class ContrastiveHead(nn.Module):
    """
    Cross-Modal Contrastive Learning Head.

    Pulls embeddings of different modalities from the same sample closer
    together in the joint space, while pushing embeddings from different
    samples apart.

    Uses InfoNCE-style loss with temperature scaling.
    """

    def __init__(self, config: FoundationConfig):
        super().__init__()
        self.config = config
        self.embed_dim = config.embed_dim
        self.temperature = config.temperature

        # Projection head for contrastive space
        self.projector = nn.Sequential(
            nn.Linear(self.embed_dim, self.embed_dim),
            nn.GELU(),
            nn.Linear(self.embed_dim, self.embed_dim),
        )

    def forward(
        self,
        joint_embedding: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute contrastive loss across modalities.

        Strategy:
        - For each pair of modalities from the SAME sample, treat them
          as positive pairs.
        - All other combinations are negatives.

        Parameters
        ----------
        joint_embedding : (batch, n_modalities, embed_dim)
            Joint embeddings from MultiModalEncoder.

        Returns
        -------
        loss : scalar tensor
            InfoNCE contrastive loss.
        """
        B, N, D = joint_embedding.shape

        if N < 2:
            return torch.tensor(0.0, device=joint_embedding.device)

        # Project to contrastive space
        z = self.projector(joint_embedding)  # (B, N, D)

        # Normalize
        z = F.normalize(z, dim=-1)

        # Reshape: (B * N, D)
        z_flat = z.reshape(B * N, D)

        # Cosine similarity matrix: (B*N, B*N)
        sim = z_flat @ z_flat.T  # (B*N, B*N)

        # Temperature scaling
        sim = sim / self.temperature

        # Build labels: for each (sample, modality) the positives are
        # other modalities from the SAME sample
        # sample_id = floor(flat_idx / N)
        # For index i, positives are same sample different modality
        sample_ids = torch.arange(B, device=z.device).unsqueeze(1).expand(B, N).reshape(-1)

        # Positive mask: same sample, different modality
        pos_mask = sample_ids.unsqueeze(0) == sample_ids.unsqueeze(1)  # (B*N, B*N)
        # Remove self-comparisons
        self_mask = torch.eye(B * N, dtype=torch.bool, device=z.device)
        pos_mask = pos_mask & ~self_mask

        # InfoNCE loss
        exp_sim = torch.exp(sim)

        # For each anchor, compute log-likelihood
        pos_sum = (exp_sim * pos_mask.float()).sum(dim=1)  # (B*N,)
        all_sum = exp_sim.sum(dim=1)  # (B*N,) — includes self

        # Only compute for anchors that have at least one positive
        loss_per_anchor = -torch.log(pos_sum / (all_sum + 1e-8) + 1e-8)

        # Mask: only anchors with ≥1 positive
        has_pos = pos_mask.any(dim=1)
        if has_pos.any():
            loss = loss_per_anchor[has_pos].mean()
        else:
            loss = torch.tensor(0.0, device=z.device)

        return loss
