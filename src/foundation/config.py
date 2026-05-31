#!/usr/bin/env python3
"""
Foundation Model Configuration
===============================

Hyperparameter management for the DeepCatch foundation model.
Provides default, prototype, and production configuration presets.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

# ── Modality dimension definitions ──────────────────────────────

MODALITY_DIMS: Dict[str, int] = {
    "frag_basic": 4,       # MFR, FSI, CAFF, FEM aggregated
    "frag_enhanced": 44,   # DELFI + MFS + nucleosome + 5-mer features
    "cnv": 6,              # chromosomal instability features
    "sero": 4,             # PG-I, PG-II, G-17, H. pylori
    "gnn": 1,              # field_defect_score from GATv2
    "tissue": 24,          # tissue deconvolution features
}

MODALITY_NAMES: Tuple[str, ...] = tuple(MODALITY_DIMS.keys())

N_MODALITIES: int = len(MODALITY_DIMS)


@dataclass
class FoundationConfig:
    """
    Configuration for the DeepCatch Foundation Model.

    Parameters
    ----------
    embed_dim : int
        Joint embedding dimension (all modalities projected here).
    n_modalities : int
        Number of input modalities (default 6).
    n_heads : int
        Number of attention heads in TransformerEncoder.
    n_layers : int
        Number of TransformerEncoder layers.
    ff_dim : int
        Feed-forward dimension in transformer layers.
    dropout : float
        Dropout rate in transformer and projections.
    mask_ratio : float
        Fraction of modalities to mask during pre-training.
    temperature : float
        Temperature for contrastive loss.
    lambda_mask : float
        Weight for masked prediction loss.
    lambda_contrast : float
        Weight for contrastive loss.
    pretrain_lr : float
        Learning rate for pre-training.
    finetune_lr : float
        Learning rate for downstream fine-tuning.
    batch_size : int
        Batch size for training.
    n_epochs : int
        Number of pre-training epochs.
    contrastive_margin : float
        Margin for contrastive loss (push negatives beyond this).
    seed : int
        Random seed for reproducibility.
    device : str
        Device for computation ('cpu', 'cuda', 'mps').

    Notes
    -----
    Total params ≈ 3-5M with default settings, optimized for
    laptop/Apple Silicon inference.
    """

    # ── Architecture ──
    embed_dim: int = 128
    n_modalities: int = 6
    n_heads: int = 4
    n_layers: int = 4
    ff_dim: int = 256
    dropout: float = 0.1

    # ── Pre-training ──
    mask_ratio: float = 0.3
    temperature: float = 0.1
    lambda_mask: float = 1.0
    lambda_contrast: float = 0.5
    pretrain_lr: float = 1e-4
    finetune_lr: float = 1e-5
    batch_size: int = 32
    n_epochs: int = 100
    contrastive_margin: float = 0.5

    # ── General ──
    seed: int = 42
    device: str = "cpu"

    def __post_init__(self):
        if self.embed_dim % self.n_heads != 0:
            raise ValueError(
                f"embed_dim ({self.embed_dim}) must be divisible by "
                f"n_heads ({self.n_heads})"
            )
        if not (0 < self.mask_ratio < 1):
            raise ValueError(f"mask_ratio must be in (0, 1), got {self.mask_ratio}")

    @property
    def head_dim(self) -> int:
        """Dimension per attention head."""
        return self.embed_dim // self.n_heads

    def to_dict(self) -> Dict:
        """Serialize config to dictionary."""
        return {
            "embed_dim": self.embed_dim,
            "n_modalities": self.n_modalities,
            "n_heads": self.n_heads,
            "n_layers": self.n_layers,
            "ff_dim": self.ff_dim,
            "dropout": self.dropout,
            "mask_ratio": self.mask_ratio,
            "temperature": self.temperature,
            "lambda_mask": self.lambda_mask,
            "lambda_contrast": self.lambda_contrast,
            "pretrain_lr": self.pretrain_lr,
            "finetune_lr": self.finetune_lr,
            "batch_size": self.batch_size,
            "n_epochs": self.n_epochs,
            "contrastive_margin": self.contrastive_margin,
            "seed": self.seed,
            "device": self.device,
        }

    @staticmethod
    def from_dict(d: Dict) -> "FoundationConfig":
        """Deserialize config from dictionary."""
        return FoundationConfig(**d)


# ── Presets ─────────────────────────────────────────────────────

DEFAULT_CONFIG = FoundationConfig()

PROTOTYPE_CONFIG = FoundationConfig(
    embed_dim=64,
    n_heads=2,
    n_layers=2,
    ff_dim=128,
    batch_size=16,
    n_epochs=10,
    dropout=0.2,
)

PRODUCTION_CONFIG = FoundationConfig(
    embed_dim=128,
    n_heads=4,
    n_layers=4,
    ff_dim=256,
    batch_size=64,
    n_epochs=200,
    dropout=0.1,
    pretrain_lr=1e-4,
    finetune_lr=1e-5,
)
