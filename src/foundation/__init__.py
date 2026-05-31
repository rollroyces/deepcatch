#!/usr/bin/env python3
"""
DeepCatch Foundation Model — Self-Supervised Multi-Modal cfDNA Pre-training
=============================================================================

A lightweight (~3-5M params) multi-modal foundation model for cfDNA
cancer screening that learns universal joint embeddings from 6 modalities:

    1. Fragmentomics basic (4D)
    2. Enhanced Fragmentomics (44D)
    3. CNV (6D)
    4. Serological (4D)
    5. GNN Methylation field_defect_score (1D)
    6. Tissue Deconvolution (24D)

Pre-training is self-supervised (no labels needed) using:
  - Masked modality prediction
  - Cross-modal contrastive learning

Downstream fine-tuning supports cancer detection, TOO prediction,
and healthy aging screening with a drop-in API compatible with
the existing CrossAttentionFusion.
"""

from .config import (
    FoundationConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    PRODUCTION_CONFIG,
    MODALITY_DIMS,
)

from .model import (
    MultiModalEncoder,
    PretrainHead,
    ContrastiveHead,
)

from .data import (
    MultiModalDataGenerator,
    MODALITY_NAMES,
    MODALITY_DIMS as _MODALITY_DIMS,
)

from .pretrain import FoundationPretrainer

from .downstream import FoundationDownstream, FoundationCompatibilityWrapper

from .synthetic_benchmark import (
    FoundationBenchmark,
    run_benchmark,
)

__all__ = [
    # Config
    "FoundationConfig",
    "DEFAULT_CONFIG",
    "PROTOTYPE_CONFIG",
    "PRODUCTION_CONFIG",
    "MODALITY_DIMS",
    # Model
    "MultiModalEncoder",
    "PretrainHead",
    "ContrastiveHead",
    # Data
    "MultiModalDataGenerator",
    "MODALITY_NAMES",
    # Training
    "FoundationPretrainer",
    # Downstream
    "FoundationDownstream",
    "FoundationCompatibilityWrapper",
    # Benchmark
    "FoundationBenchmark",
    "run_benchmark",
]
