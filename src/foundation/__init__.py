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

Lazy re-exports
---------------
The submodules below pull in torch (``model``, ``pretrain``,
``downstream``). Importing them eagerly at package-import time made
``import src.foundation`` fail in any environment without torch --
including the "Core tests" CI job, which installs only
numpy/scipy/scikit-learn/pandas/matplotlib/seaborn/pytest. That
aborted the whole pytest session at collection, so a single optional
dependency silenced every test in the job.

PEP 562 module-level ``__getattr__`` keeps the public API identical
(``from src.foundation import FoundationDownstream`` still works) while
importing each submodule only when one of its names is actually
requested. ``src.foundation.regularized_fusion`` is torch-free and now
importable in a torch-free environment.
"""
from __future__ import annotations

import importlib

# Public name -> submodule that defines it.
_LAZY_EXPORTS = {
    # Config
    "FoundationConfig": ".config",
    "DEFAULT_CONFIG": ".config",
    "PROTOTYPE_CONFIG": ".config",
    "PRODUCTION_CONFIG": ".config",
    "MODALITY_DIMS": ".config",
    # Model
    "MultiModalEncoder": ".model",
    "PretrainHead": ".model",
    "ContrastiveHead": ".model",
    # Data
    "MultiModalDataGenerator": ".data",
    "MODALITY_NAMES": ".data",
    "_MODALITY_DIMS": ".data",
    # Training
    "FoundationPretrainer": ".pretrain",
    # Downstream
    "FoundationDownstream": ".downstream",
    "FoundationCompatibilityWrapper": ".downstream",
    # Benchmark
    "FoundationBenchmark": ".synthetic_benchmark",
    "run_benchmark": ".synthetic_benchmark",
}

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


def __getattr__(name: str):
    """Resolve a public name to its submodule on first access."""
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(
            f"module 'src.foundation' has no attribute {name!r}")
    module = importlib.import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value  # cache, so this runs once per name
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY_EXPORTS))
