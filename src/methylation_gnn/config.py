#!/usr/bin/env python3
"""
GNN Methylation Network Configuration
======================================

Dataclass-based hyperparameter management for the GNN methylation branch.
All config values have sensible defaults validated through the implementation
plan's architecture design (GATv2Conv backbone, dual head, 3-layer depth).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Union, Literal

import numpy as np

try:
    import torch

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None


# ── Node Feature Specification ─────────────────────────────────

# Feature names and descriptions for documentation and serialization
NODE_FEATURE_SPEC: Dict[str, Dict[str, str]] = {
    "mean_methylation": {
        "desc": "Mean methylation β-value across CpGs in region",
        "range": "[0, 1]",
        "source": "cfDNA or TCGA methylation array",
    },
    "methylation_entropy": {
        "desc": "Shannon entropy of per-CpG methylation pattern",
        "range": "[0, log2(n_cpg)]",
        "source": "Computed from β-value distribution",
    },
    "methylation_variance": {
        "desc": "Variance across CpGs in region (epiallelic heterogeneity)",
        "range": "[0, 0.25]",
        "source": "Computed from β-value distribution",
    },
    "cpg_density": {
        "desc": "Number of CpG dinucleotides per 100 bp",
        "range": "[0, 50+]",
        "source": "Reference genome (hg38)",
    },
    "cpg_obs_exp": {
        "desc": "Observed/expected CpG ratio (normalized for GC)",
        "range": "[0, ~2.0]",
        "source": "Reference genome (hg38)",
    },
    "gc_content": {
        "desc": "Fraction of G+C bases in region",
        "range": "[0, 1]",
        "source": "Reference genome (hg38)",
    },
    "coverage_depth": {
        "desc": "Normalized cfDNA read coverage at region",
        "range": "[0, ∞)",
        "source": "cfDNA BAM/coverage profile",
    },
    "fragment_size_mean": {
        "desc": "Mean fragment length at region (bp)",
        "range": "[90, 250]",
        "source": "cfDNA fragment length profile",
    },
    "fragment_short_frac": {
        "desc": "Fraction of fragments < 150 bp at region",
        "range": "[0, 1]",
        "source": "cfDNA fragment length profile",
    },
    "end_motif_diversity": {
        "desc": "4-mer motif diversity score at region boundaries",
        "range": "[0, 1]",
        "source": "cfDNA fragment end motifs",
    },
    # Chromatin features (from public reference data, optional)
    "dnase_signal": {
        "desc": "Mean DNase-seq signal (chromatin accessibility)",
        "range": "[0, ∞) log-normalized",
        "source": "ENCODE/Roadmap DNase-seq",
    },
    "h3k4me3_signal": {
        "desc": "Active promoter histone mark",
        "range": "[0, ∞) log-normalized",
        "source": "ENCODE/Roadmap ChIP-seq",
    },
    "h3k27ac_signal": {
        "desc": "Active enhancer histone mark",
        "range": "[0, ∞) log-normalized",
        "source": "ENCODE/Roadmap ChIP-seq",
    },
    "h3k27me3_signal": {
        "desc": "Polycomb repressive mark",
        "range": "[0, ∞) log-normalized",
        "source": "ENCODE/Roadmap ChIP-seq",
    },
    "h3k9me3_signal": {
        "desc": "Constitutive heterochromatin mark",
        "range": "[0, ∞) log-normalized",
        "source": "ENCODE/Roadmap ChIP-seq",
    },
    # One-hot region type (5 categories)
    "region_type_cpg": {"desc": "Is CpG island type", "range": "{0,1}"},
    "region_type_enhancer": {"desc": "Is enhancer type", "range": "{0,1}"},
    "region_type_promoter": {"desc": "Is promoter type", "range": "{0,1}"},
    "region_type_ctcf": {"desc": "Is CTCF binding site", "range": "{0,1}"},
    "region_type_dhs": {"desc": "Is DNase hypersensitive site", "range": "{0,1}"},
}


@dataclass
class GNNConfig:
    """
    Configuration for GNN methylation network reconstruction.

    All hyperparameters have been chosen based on:
    - GATv2 theoretical advantages over GCN/GAT for biological graphs
      (Brody et al. 2022, ICLR)
    - 3-layer architecture for sufficient receptive field across ~50K nodes
    - Self-supervised masked node prediction established for molecular graphs
      (Hu et al. 2020, NeurIPS)

    Parameters
    ----------
    n_nodes : int
        Number of genomic regulatory regions. 50K covers major CpG islands,
        enhancers, and promoters in hg38. Increase to 100K for CTCF sites.
    edge_k : int
        Maximum edges per node (k-NN truncation for scalability). 20 edges
        gives sufficient neighborhood without quadratic blowup.
    reference_genome : str
        Genome build ('hg38' or 'hg19').
    n_node_features : int
        Raw features per node (methylation β, CpG density, GC, coverage,
        fragment length, chromatin marks, region type one-hot).
    hidden_dims : list of int
        Per-layer hidden dimensions for message passing.
    n_attention_heads : int
        Number of GATv2 attention heads per layer.
    n_edge_types : int
        Number of distinct edge relation types (physical_interaction,
        co_methylation, co_fragmentation, genomic_proximity,
        regulatory_domain).
    dropout : float
        Dropout rate applied after each GNN layer.
    decoder_hidden : list of int
        Hidden dimensions for the reconstruction decoder MLP.
    anomaly_hidden : list of int
        Hidden dimensions for the anomaly scoring MLP.
    learning_rate : float
        AdamW learning rate (1e-3 works well for GATv2 on molecular graphs).
    weight_decay : float
        AdamW weight decay for regularization.
    n_epochs_pretrain : int
        Epochs for Phase 1: self-supervised masked node prediction.
    n_epochs_finetune : int
        Epochs for Phase 2: joint reconstruction + anomaly training.
    mask_ratio_pretrain : float
        Fraction of node features masked during pretraining.
    lambda_anomaly : float
        Weight for anomaly loss in joint training (up-weight if
        anomaly head dominates).
    lambda_temporal : float
        Weight for temporal consistency loss (optional, only if
        longitudinal data available).
    device : str
        'auto' (probe), 'cuda', 'mps', or 'cpu'.
    checkpoint_dir : str
        Directory for model checkpoints during training.
    patience : int
        Early stopping patience in epochs (validation loss).
    log_every : int
        Log training metrics every N epochs.
    """

    # ── Graph construction ──
    n_nodes: int = 50_000
    edge_k: int = 20
    reference_genome: str = "hg38"
    alpha_hic_weight: float = 0.7  # Hi-C vs co-fragmentation blend

    # ── Node features ──
    n_node_features: int = 20

    # ── GNN architecture ──
    hidden_dims: List[int] = field(default_factory=lambda: [64, 128, 256])
    n_attention_heads: int = 4
    n_edge_types: int = 5
    dropout: float = 0.2
    decoder_hidden: List[int] = field(default_factory=lambda: [128, 64])
    anomaly_hidden: List[int] = field(default_factory=lambda: [128, 64])

    # ── Training ──
    batch_size: int = 1  # one graph per sample → graph-level batch
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    n_epochs_pretrain: int = 100
    n_epochs_finetune: int = 50
    mask_ratio_pretrain: float = 0.30
    lambda_anomaly: float = 0.5
    lambda_temporal: float = 0.1

    # ── Device ──
    device: str = "auto"

    # ── Logging / checkpointing ──
    checkpoint_dir: str = "checkpoints/methylation_gnn"
    patience: int = 15
    log_every: int = 10

    def __post_init__(self):
        """Validate and set derived parameters."""
        if self.hidden_dims is None or len(self.hidden_dims) < 1:
            raise ValueError("hidden_dims must be a non-empty list")
        if self.n_node_features <= 0:
            raise ValueError("n_node_features must be positive")
        if self.n_edge_types < 1:
            raise ValueError("n_edge_types must be at least 1")
        if self.device == "auto":
            self.device = self._detect_device()

    @staticmethod
    def _detect_device() -> str:
        """Auto-detect best available PyTorch device."""
        if not _HAS_TORCH:
            return "cpu"
        try:
            if torch.cuda.is_available():
                return "cuda"
        except Exception:
            pass
        try:
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
        except Exception:
            pass
        return "cpu"

    @property
    def n_layers(self) -> int:
        """Number of GNN message-passing layers."""
        return len(self.hidden_dims)

    @property
    def final_dim(self) -> int:
        """Dimension after all GNN layers (used by decoder/anomaly heads)."""
        return self.hidden_dims[-1]

    @property
    def estimated_params(self) -> int:
        """
        Estimate total trainable parameters.
        Rough calculation; actual count depends on PyG implementation details.
        Target is ~2M parameters (fits comfortably on any device).
        """
        d_in = self.n_node_features
        total = 0

        # Input projection: Linear(d_in → hidden_dims[0])
        total += d_in * self.hidden_dims[0] + self.hidden_dims[0]

        # GNN layers: GATv2Conv per edge type per layer
        for i in range(len(self.hidden_dims) - 1):
            d_curr = self.hidden_dims[i]
            d_next = self.hidden_dims[i + 1]
            d_per_head = d_next // self.n_attention_heads
            # GATv2: W_src + W_dst + a_src + a_dst
            per_type = (
                d_curr * d_per_head * self.n_attention_heads  # src
                + d_curr * d_per_head * self.n_attention_heads  # dst
                + d_per_head * self.n_attention_heads  # att_src
                + d_per_head * self.n_attention_heads  # att_dst
            )
            total += self.n_edge_types * per_type
            # BatchNorm
            total += 2 * d_next

        # Reconstruction decoder
        d_last = self.hidden_dims[-1]
        decoder_layers = [d_last] + self.decoder_hidden + [self.n_node_features]
        for i in range(len(decoder_layers) - 1):
            total += decoder_layers[i] * decoder_layers[i + 1]
            total += decoder_layers[i + 1]

        # Anomaly head
        anom_layers = [d_last] + self.anomaly_hidden + [1]
        for i in range(len(anom_layers) - 1):
            total += anom_layers[i] * anom_layers[i + 1]
            total += anom_layers[i + 1]

        # Mask token
        total += self.n_node_features

        return total

    def to_dict(self) -> dict:
        """Serialize config to dictionary (for logging/checkpoints)."""
        return {
            "n_nodes": self.n_nodes,
            "edge_k": self.edge_k,
            "n_node_features": self.n_node_features,
            "hidden_dims": self.hidden_dims,
            "n_attention_heads": self.n_attention_heads,
            "n_edge_types": self.n_edge_types,
            "dropout": self.dropout,
            "decoder_hidden": self.decoder_hidden,
            "anomaly_hidden": self.anomaly_hidden,
            "learning_rate": self.learning_rate,
            "weight_decay": self.weight_decay,
            "n_epochs_pretrain": self.n_epochs_pretrain,
            "n_epochs_finetune": self.n_epochs_finetune,
            "mask_ratio_pretrain": self.mask_ratio_pretrain,
            "lambda_anomaly": self.lambda_anomaly,
            "device": self.device,
        }


# ── Pre-built defaults ──────────────────────────────────────────

# Standard configuration for 50K-node methylation graphs
DEFAULT_GNN_CONFIG = GNNConfig()

# Smaller config for rapid prototyping / CI tests
PROTOTYPE_GNN_CONFIG = GNNConfig(
    n_nodes=5_000,
    edge_k=10,
    hidden_dims=[32, 64, 128],
    n_attention_heads=2,
    n_epochs_pretrain=20,
    n_epochs_finetune=10,
)

# Large config for full-scale 100K-node production runs
PRODUCTION_GNN_CONFIG = GNNConfig(
    n_nodes=100_000,
    edge_k=30,
    hidden_dims=[128, 256, 512],
    n_attention_heads=8,
    n_epochs_pretrain=200,
    n_epochs_finetune=100,
)
