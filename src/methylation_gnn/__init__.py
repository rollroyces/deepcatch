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

from .config import GNNConfig, DEFAULT_GNN_CONFIG
from .graph_builder import (
    RegulatoryGraphBuilder,
    GenomicRegion,
    REGION_TYPES,
)
from .gnn_model import MethylationGNN
from .gnn_trainer import GNNTrainer, GNNTrainerPhase
from .gnn_inference import GNNInference, MethylationGNNPredictor
from .data import (
    ReferenceDataCatalog,
    download_ucsc_cpg_islands,
    download_encode_hic,
    download_gencode_promoters,
    download_fantom5_enhancers,
    preprocess_methylation_betas,
)
from .integration import (
    MethylationBranchAdapter,
    ModularArmsBuilder,
    extend_fusion_with_gnn,
)

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
]
