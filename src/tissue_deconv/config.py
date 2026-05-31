#!/usr/bin/env python3
"""
Tissue Deconvolution Configuration
====================================

Dataclass-based hyperparameter management for the cfSort-style
tissue deconvolution module.

Reference
---------
Li et al. (2023) PNAS — Comprehensive tissue methylation atlas
with 29 tissues × 521 samples, DNN-based deconvolution.
Detection limit: 0.1% tissue fraction at 20× coverage.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Dict

# ── 29 Tissues from cfSort reference atlas ──────────────────────

SUPPORTED_TISSUES: List[str] = [
    # Blood / immune
    "Whole Blood",
    "PBMC",
    "CD4+ T Cell",
    "CD8+ T Cell",
    "CD19+ B Cell",
    "CD56+ NK Cell",
    "Monocyte",
    "Neutrophil",
    "Eosinophil",
    # Solid organs — gastrointestinal
    "Liver",
    "Pancreas",
    "Stomach",
    "Colon",
    "Small Intestine",
    "Esophagus",
    # Solid organs — other
    "Lung",
    "Heart",
    "Kidney",
    "Bladder",
    "Prostate",
    "Breast",
    "Thyroid",
    "Adrenal Gland",
    "Ovary",
    # CNS
    "Brain (Cortex)",
    "Brain (Cerebellum)",
    # Other
    "Skeletal Muscle",
    "Adipose Tissue",
    "Placenta",
]

# Tissues most relevant to cancer detection
CANCER_RELEVANT_TISSUES: List[str] = [
    "Colon",
    "Stomach",
    "Pancreas",
    "Liver",
    "Esophagus",
    "Lung",
    "Breast",
    "Prostate",
    "Ovary",
    "Bladder",
    "Thyroid",
    "Kidney",
    "Brain (Cortex)",
]

# Blood-derived cfDNA background (normally dominant)
BLOOD_DERIVED_TISSUES: List[str] = [
    "Whole Blood",
    "PBMC",
    "Neutrophil",
    "Monocyte",
    "CD4+ T Cell",
    "CD8+ T Cell",
    "CD19+ B Cell",
    "CD56+ NK Cell",
    "Eosinophil",
]


@dataclass
class TissueDeconvConfig:
    """
    Configuration for cfSort-style tissue deconvolution DNN.

    Parameters
    ----------
    n_tissues : int
        Number of tissue types to deconvolve (default 29 from cfSort).
    n_cpg_features : int
        Number of tissue-discriminative CpG markers to use.
        cfSort uses ~1000 after feature selection.
    hidden_dims : list of int
        Hidden layer dimensions for the DNN.
    dropout : float
        Dropout rate after each hidden layer.
    chunk_size : int
        Processing chunk size for large methylome data (memory efficiency).
    detection_limit : float
        Minimum detectable tissue fraction (0.001 = 0.1%).
    n_ensemble : int
        Number of models in ensemble for robust prediction.
    learning_rate : float
        Adam optimizer learning rate.
    weight_decay : float
        L2 regularization strength.
    n_epochs : int
        Number of training epochs.
    batch_size : int
        Training batch size.
    patience : int
        Early stopping patience in epochs.
    lambda_sparsity : float
        Weight for sparsity penalty (most tissues should be ~0).
    lambda_smooth : float
        Weight for smoothness / entropy regularization.
    device : str
        'auto' (probe), 'mps', 'cuda', or 'cpu'.
    checkpoint_dir : str
        Directory for model checkpoints.
    """

    # ── Tissue atlas ──
    n_tissues: int = 29
    n_cpg_features: int = 1000

    # ── DNN architecture ──
    hidden_dims: List[int] = field(default_factory=lambda: [256, 128, 64])
    dropout: float = 0.3

    # ── Processing ──
    chunk_size: int = 10_000
    detection_limit: float = 0.001  # 0.1% tissue fraction

    # ── Ensemble ──
    n_ensemble: int = 3

    # ── Training ──
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    n_epochs: int = 100
    batch_size: int = 64
    patience: int = 15
    lambda_sparsity: float = 0.01
    lambda_smooth: float = 0.001

    # ── Device ──
    device: str = "auto"

    # ── Checkpointing ──
    checkpoint_dir: str = "checkpoints/tissue_deconv"

    def __post_init__(self):
        """Validate configuration."""
        if self.n_tissues != len(SUPPORTED_TISSUES):
            raise ValueError(
                f"n_tissues ({self.n_tissues}) must match "
                f"SUPPORTED_TISSUES length ({len(SUPPORTED_TISSUES)})"
            )
        if self.n_cpg_features < 10:
            raise ValueError("n_cpg_features must be ≥ 10")
        if self.hidden_dims is None or len(self.hidden_dims) < 1:
            raise ValueError("hidden_dims must be a non-empty list")
        if self.n_ensemble < 1:
            raise ValueError("n_ensemble must be ≥ 1")
        if self.detection_limit <= 0:
            raise ValueError("detection_limit must be > 0")
        if self.device == "auto":
            self.device = self._detect_device()

    @staticmethod
    def _detect_device() -> str:
        """Auto-detect best available PyTorch device."""
        try:
            import torch
            if torch.cuda.is_available():
                return "cuda"
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
        except ImportError:
            pass
        return "cpu"

    @property
    def estimated_params(self) -> int:
        """Estimate total trainable parameters for one model."""
        total = 0
        prev = self.n_cpg_features

        for h in self.hidden_dims:
            # Linear: weights + bias
            total += prev * h + h
            # BatchNorm: 2 params per unit (gamma, beta)
            total += 2 * h
            prev = h

        # Output layer: last hidden → n_tissues
        total += prev * self.n_tissues + self.n_tissues

        return total

    def to_dict(self) -> Dict:
        """Serialize config to dictionary."""
        return {
            "n_tissues": self.n_tissues,
            "n_cpg_features": self.n_cpg_features,
            "hidden_dims": self.hidden_dims,
            "dropout": self.dropout,
            "chunk_size": self.chunk_size,
            "detection_limit": self.detection_limit,
            "n_ensemble": self.n_ensemble,
            "learning_rate": self.learning_rate,
            "weight_decay": self.weight_decay,
            "n_epochs": self.n_epochs,
            "batch_size": self.batch_size,
            "patience": self.patience,
            "lambda_sparsity": self.lambda_sparsity,
            "device": self.device,
        }


# ── Pre-built defaults ──────────────────────────────────────────

DEFAULT_CONFIG = TissueDeconvConfig()

PROTOTYPE_CONFIG = TissueDeconvConfig(
    n_cpg_features=100,
    hidden_dims=[64, 32, 16],
    n_ensemble=2,
    n_epochs=10,
    dropout=0.2,
)
