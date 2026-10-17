"""
Ultra-sensitive screening modules ported from cancer-screening (Apr 2026).

This package consolidates the CANCER-SEEK scaffold:
  - UMI-based error correction for ultra-low VAF
  - Low-VAF variant calling with CHIP filtering
  - CNV detection (CBS segmentation)
  - Methylation beta-value processing
  - XGBoost baseline with SHAP explanations
  - β-VAE anomaly detector
  - Multi-modal FusionNet (cross-attention)
  - Methylation Transformer encoder
  - Spike-in LOD validation

All modules in this package are optional dependencies — they do not change
DeepCatch's existing 189-test fast-path suite. They are guarded by
``pip install 'deepcatch[ultra]'`` extras.

References:
  - ARCHITECTURE: ../../docs/CANCER_SEEK_ARCHITECTURE.md
  - ROADMAP:      ../../docs/ROADMAP_CANCER_SEEK.md
  - DATA SCIENCE: ../../docs/legacy/DATA_SCIENCE_ML.md
  - Source repo:  https://github.com/rollroyces/cancer-screening (now archived)
"""

from .umi import UMICorrector, UMIFamily, UMIRead, DuplexConsensus
from .cnv_analyzer import CNVAnalyzer
from .variant_calling import LowVAFVariantFilter, VariantCaller, VariantCall
from .methylation import MethylationProcessor
from .spike_in_simulator import SpikeInSimulator, SpikeInConfig, run_lod_validation

__all__ = [
    "UMICorrector",
    "UMIFamily",
    "UMIRead",
    "DuplexConsensus",
    "CNVAnalyzer",
    "LowVAFVariantFilter",
    "VariantCaller",
    "VariantCall",
    "MethylationProcessor",
    "SpikeInSimulator",
    "SpikeInConfig",
    "run_lod_validation",
]

__version__ = "2.3.0"
