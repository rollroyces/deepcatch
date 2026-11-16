"""
DeepCatch Methylation GNN — Epigenetic Field Defect Detection
==============================================================

Graph Neural Network branch for Stage 1 (Capture) of the DeepCatch CET pipeline.
Detects pre-cancer epigenetic field defects through methylation network reconstruction:

1. **Graph Builder** — Constructs heterogeneous methylation graphs from cfDNA +
   reference data (Hi-C 3D contacts, co-methylation patterns, co-fragmentation).

2. **GNN Model** — GATv2-based graph attention with dual heads:
   - *Reconstruction Decoder*: self-supervised masked node prediction
   - *Anomaly Head*: per-node and graph-level anomaly scoring

3. **Trainer** — Three-phase training: masked pretraining → joint reconstruction
   + anomaly → fine-tuning with labels.

4. **Inference** — Lightweight inference pipeline producing a scalar
   ``field_defect_score`` for fusion with fragmentomics, CNV, and serological
   modalities.

.. rubric:: References

.. [1] Rao, S.S. et al. (2014) Cell 159:1665-1680. Hi-C chromatin contact maps.
.. [2] Brody, S. et al. (2022) "How Attentive are Graph Attention Networks?"
       ICLR 2022. GATv2 theoretical foundations.
.. [3] Velickovic, P. et al. (2018) ICLR. Graph Attention Networks.
.. [4] Schlichtkrull, M. et al. (2018) ESWC. Relational GCN for link prediction.
.. [5] Hu, W. et al. (2020) NeurIPS. Self-supervised graph pretraining.
"""

from __future__ import annotations

import importlib

from .config import GNNConfig, DEFAULT_GNN_CONFIG
from .graph_builder import (
    RegulatoryGraphBuilder,
    GenomicRegion,
    REGION_TYPES,
)
from .integration import (
    MethylationBranchAdapter,
    ModularArmsBuilder,
    extend_fusion_with_gnn,
)
from .data import (
    ReferenceDataCatalog,
    download_ucsc_cpg_islands,
    download_encode_hic,
    download_gencode_promoters,
    download_fantom5_enhancers,
    preprocess_methylation_betas,
)
from .real_data import (
    RealDataError,
    build_real_cohort_graphs,
    load_finaleme_cohort,
    load_real_regions,
    per_sample_methylation_data,
)

# Names that need torch (or torch_geometric) are resolved lazily.
#
# ``gnn_model``, ``gnn_trainer`` and ``gnn_inference`` subclass
# ``torch.nn.Module`` at import time, so importing them eagerly made
# ``import src.methylation_gnn`` fail in any environment without torch --
# including the "Core tests" CI job, which installs only
# numpy/scipy/scikit-learn/pandas/matplotlib/seaborn/pytest. That
# aborted the whole pytest session at collection, so one optional
# dependency silenced every test in the job.
#
# PEP 562 module-level ``__getattr__`` keeps the public API identical
# (``from src.methylation_gnn import MethylationGNN`` still works) while
# importing each submodule only when one of its names is actually
# requested. ``src.methylation_gnn.real_data`` is torch-free and is
# importable without torch.

# Public name -> submodule that defines it.
_LAZY_EXPORTS = {
    "MethylationGNN": ".gnn_model",
    "GNNTrainer": ".gnn_trainer",
    "GNNTrainerPhase": ".gnn_trainer",
    "GNNInference": ".gnn_inference",
    "MethylationGNNPredictor": ".gnn_inference",
}


def __getattr__(name: str):
    """Resolve a torch-dependent public name to its submodule on first access.

    Raises the import error from the submodule itself when torch is
    missing, so the message names the missing dependency rather than
    failing with "NoneType has no attribute Module" from deep inside a
    class definition.
    """
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(
            f"module 'src.methylation_gnn' has no attribute {name!r}")
    module = importlib.import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value  # cache, so this runs once per name
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY_EXPORTS))


__all__ = [
    # Config
    "GNNConfig",
    "DEFAULT_GNN_CONFIG",
    # Graph builder
    "RegulatoryGraphBuilder",
    "GenomicRegion",
    "REGION_TYPES",
    # GNN model
    "MethylationGNN",
    # Trainer
    "GNNTrainer",
    "GNNTrainerPhase",
    # Inference
    "GNNInference",
    "MethylationGNNPredictor",
    # Reference data
    "ReferenceDataCatalog",
    "download_ucsc_cpg_islands",
    "download_encode_hic",
    "download_gencode_promoters",
    "download_fantom5_enhancers",
    "preprocess_methylation_betas",
    # Integration
    "MethylationBranchAdapter",
    "ModularArmsBuilder",
    "extend_fusion_with_gnn",
    # Real (non-synthetic) cohorts
    "RealDataError",
    "load_finaleme_cohort",
    "per_sample_methylation_data",
    "load_real_regions",
    "build_real_cohort_graphs",
]
