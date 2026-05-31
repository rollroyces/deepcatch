"""
DeepCatch Tissue Deconvolution — cfSort-Style DNN Head
========================================================

Stage 1 Modality (Modality 6): Tissue-of-origin deconvolution from
cfDNA methylation beta values using a cfSort-inspired DNN architecture.

Overview
--------

This module adds tissue deconvolution as a 6th modality in the DeepCatch
Stage 1 (Capture) pipeline. The DNN head predicts per-tissue cell death
proportions from cfDNA methylation data, enabling:

1. **Tissue-of-origin identification**: Which organ is shedding cfDNA?
2. **Tumor signal localization**: Elevated fractions in cancer-relevant
   tissues may indicate tumor-derived cfDNA.
3. **Multi-cancer signal**: Two tissues with similarly elevated fractions
   can suggest co-pathology or metastatic spread.

Architecture
------------

- **TissueAtlas**: Stores 29-tissue reference methylation profiles
  (synthetic or real cfSort data).
- **TissueDeconvolutionModel**: Lightweight DNN (~500K params) with
  [256, 128, 64] hidden layers, BatchNorm, ReLU, and Dropout.
- **TissueDeconvolutionEnsemble**: 3 models with different seeds,
  averaging predictions for robustness.
- **TissueDeconvTrainer**: Training on synthetic mixtures with
  KL divergence + L1 sparsity + entropy regularization.
- **TissueDeconvolutionFeatures**: 24 scalar features extracted
  from tissue fractions for fusion.
- **DEConvIntegration**: Adapter for CrossAttentionFusion compatible
  with the existing 5-modality pipeline.

Fallback Design
---------------

The module is designed to work **without** real reference data:
- Synthetic reference profiles are generated on-the-fly
- Training works entirely on synthetic mixtures
- Feature extraction returns zeros when no model is available
- All public APIs handle missing dependencies gracefully

Usage
-----

.. code-block:: python

    from src.tissue_deconv import DEConvIntegration, TissueDeconvConfig

    # With pre-trained model
    deconv = DEConvIntegration(checkpoint="checkpoints/deconv.pt")
    scores = deconv.to_modality(sample)

    # Train from scratch on synthetic data
    deconv = DEConvIntegration()
    deconv.fit_synthetic(n_samples=2000)
    fractions = deconv.predict_tissue_fractions(methylation)

    # Extract features for fusion
    features = deconv.extract_all(sample)

References
----------
.. [1] Li, S. et al. (2023). "Comprehensive tissue methylation atlas
       and DNN-based deconvolution of cfDNA." PNAS.
.. [2] Moss, J. et al. (2018). "Comprehensive human cell-type
       methylation atlas reveals origins of circulating cell-free
       DNA in health and disease." Nature Communications.
"""

from .config import (
    TissueDeconvConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    SUPPORTED_TISSUES,
    CANCER_RELEVANT_TISSUES,
    BLOOD_DERIVED_TISSUES,
)
from .tissue_atlas import TissueAtlas
from .model import (
    TissueDeconvolutionModel,
    TissueDeconvolutionEnsemble,
    DeconvLoss,
)
from .trainer import TissueDeconvTrainer
from .tissue_features import TissueDeconvolutionFeatures
from .integration import DEConvIntegration

__all__ = [
    # Config
    "TissueDeconvConfig",
    "DEFAULT_CONFIG",
    "PROTOTYPE_CONFIG",
    "SUPPORTED_TISSUES",
    "CANCER_RELEVANT_TISSUES",
    "BLOOD_DERIVED_TISSUES",
    # Atlas
    "TissueAtlas",
    # Model
    "TissueDeconvolutionModel",
    "TissueDeconvolutionEnsemble",
    "DeconvLoss",
    # Trainer
    "TissueDeconvTrainer",
    # Features
    "TissueDeconvolutionFeatures",
    # Integration
    "DEConvIntegration",
]
