#!/usr/bin/env python3
"""
Methylation GNN — GATv2-Based Field Defect Detector
=====================================================

Graph Neural Network for detecting pre-cancer epigenetic field defects
through methylation network reconstruction.

Architecture
------------

.. image:: https://via.placeholder.com/800x400?text=GNN+Architecture+Diagram

::

    Node Features (N × 20)
         │
         ▼
    Feature Projection: Linear(20→64) + ReLU + LayerNorm
         │
         ▼
    ┌─────────────────────────────────────────────┐
    │  GATv2Conv Layer 1: 64 → 128 (4 heads)      │
    │  + BatchNorm + ReLU + Dropout(0.2)           │
    ├─────────────────────────────────────────────┤
    │  GATv2Conv Layer 2: 128 → 256 (4 heads)     │
    │  + BatchNorm + ReLU + Dropout(0.2)           │
    ├─────────────────────────────────────────────┤
    │  GATv2Conv Layer 3: 256 → 256 (4 heads)     │
    │  + BatchNorm + ReLU + Dropout(0.2)           │
    └──────────────────┬──────────────────────────┘
                       │
         ┌─────────────┴────────────┐
         ▼                          ▼
    Reconstruction Decoder    Anomaly Head
    MLP(256→128→64→20)       MLP(256→128→64→1)
         │                          │
         ▼                          ▼
    x̂ (reconstructed)       a_i ∈ [0,1] per node
         │                          │
         └──────────┬───────────────┘
                    ▼
    field_defect_score = weighted_reconstruction_error

Key design decisions
--------------------

- **GATv2Conv** over GCN/GAT: dynamic attention computes attention weights
  *after* the linear transformation, enabling more expressive attention
  patterns. This is critical for biological graphs where edge importance
  varies by node state (e.g., a promoter is more important when methylated).

- **Dual head**: reconstruction decoder (self-supervised) learns the normal
  methylation network; anomaly head learns which reconstruction errors are
  diagnostically relevant.

- **Heterogeneous edge support**: model accepts ``edge_type`` tensor for
  type-specific message passing (via separate GATv2Conv per type, aggregated
  via mean pooling). Falls back gracefully to homogeneous GATv2Conv when
  ``n_edge_types=1``.

References
----------
.. [1] Brody, S. et al. (2022). "How Attentive are Graph Attention Networks?"
       ICLR 2022.
.. [2] Velickovic, P. et al. (2018). "Graph Attention Networks." ICLR.
.. [3] Hu, W. et al. (2020). "Strategies for Pre-training Graph Neural
       Networks." ICLR.
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.nn import GATv2Conv, BatchNorm, global_mean_pool

    _HAS_PYG = True
except ImportError:
    _HAS_PYG = False
    torch = None
    nn = None
    F = None
    GATv2Conv = None
    BatchNorm = None
    global_mean_pool = None

logger = logging.getLogger(__name__)


# ── Device helper ───────────────────────────────────────────────

def _get_device() -> str:
    """Return the best available PyTorch device."""
    if torch is None:
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


# ── Sub-modules ─────────────────────────────────────────────────

class ReconstructionDecoder(nn.Module):
    """
    Decoder that reconstructs original node features from GNN embeddings.

    Maps from the final GNN hidden dimension back to the original
    feature space. Used in self-supervised pretraining for the
    masked node prediction task.
    """

    def __init__(
        self,
        in_dim: int,
        hidden_dims: List[int],
        out_dim: int,
        dropout: float = 0.2,
    ):
        super().__init__()
        layers = []
        prev_dim = in_dim
        for h_dim in hidden_dims:
            layers.append(nn.Linear(prev_dim, h_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            prev_dim = h_dim
        layers.append(nn.Linear(prev_dim, out_dim))
        self.mlp = nn.Sequential(*layers)

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Parameters
        ----------
        embeddings : (N, d_final)
            Node embeddings from the last GNN layer.

        Returns
        -------
        reconstructed : (N, d_features)
            Reconstructed node features.
        """
        return self.mlp(embeddings)


class AnomalyHead(nn.Module):
    """
    MLP head that scores each node's reconstruction error significance.

    Learns to assign higher anomaly scores to nodes whose reconstruction
    errors are diagnostically informative (i.e., their dysregulation
    indicates a field defect). Output is in [0, 1] via sigmoid.
    """

    def __init__(
        self,
        in_dim: int,
        hidden_dims: List[int],
        dropout: float = 0.2,
    ):
        super().__init__()
        layers = []
        prev_dim = in_dim
        for h_dim in hidden_dims:
            layers.append(nn.Linear(prev_dim, h_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            prev_dim = h_dim
        layers.append(nn.Linear(prev_dim, 1))
        layers.append(nn.Sigmoid())
        self.mlp = nn.Sequential(*layers)

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Parameters
        ----------
        embeddings : (N, d_final)
            Node embeddings from the last GNN layer.

        Returns
        -------
        scores : (N, 1)
            Per-node anomaly scores in [0, 1].
        """
        return self.mlp(embeddings)


# ── Main GNN Model ──────────────────────────────────────────────

class MethylationGNN(nn.Module):
    """
    GATv2-based Graph Neural Network for methylation field defect detection.

    Three-layer message passing with:
    - Learnable mask token for self-supervised pretraining
    - Reconstruction decoder (predicts original node features from context)
    - Anomaly head (per-node anomaly scores → graph-level field defect score)

    Parameters
    ----------
    n_node_features : int
        Number of raw features per node (default 20).
    hidden_dims : list of int
        Per-layer hidden dimensions. Length determines number of GNN layers.
        Default: [64, 128, 256] → 3 layers.
    n_edge_types : int
        Number of distinct edge types. If 1, uses standard homogeneous
        GATv2Conv. If > 1, uses separate GATv2Conv per edge type with
        mean aggregation (simplified HAN).
    n_attention_heads : int
        Number of attention heads per GATv2Conv layer.
    decoder_hidden : list of int
        Hidden dimensions for reconstruction decoder MLP.
    anomaly_hidden : list of int
        Hidden dimensions for anomaly scoring MLP.
    dropout : float
        Dropout rate (applied after each GNN layer activation).
    """

    def __init__(
        self,
        n_node_features: int = 20,
        hidden_dims: Optional[List[int]] = None,
        n_edge_types: int = 5,
        n_attention_heads: int = 4,
        decoder_hidden: Optional[List[int]] = None,
        anomaly_hidden: Optional[List[int]] = None,
        dropout: float = 0.2,
    ):
        if not _HAS_PYG:
            raise ImportError(
                "PyTorch Geometric is required for MethylationGNN. "
                "Install with: pip install torch_geometric"
            )

        super().__init__()

        if hidden_dims is None:
            hidden_dims = [64, 128, 256]
        if decoder_hidden is None:
            decoder_hidden = [128, 64]
        if anomaly_hidden is None:
            anomaly_hidden = [128, 64]

        self.n_node_features = n_node_features
        self.hidden_dims = hidden_dims
        self.n_edge_types = n_edge_types
        self.n_attention_heads = n_attention_heads
        self.dropout_rate = dropout

        # ── Input projection ──
        self.input_proj = nn.Sequential(
            nn.Linear(n_node_features, hidden_dims[0]),
            nn.ReLU(),
            nn.LayerNorm(hidden_dims[0]),
        )

        # ── GNN layers ──
        # For heterogeneous edges (n_edge_types > 1): build a ModuleList
        # of GATv2Conv per type per layer, then aggregate.
        # For homogeneous: single GATv2Conv per layer.
        self.n_layers = len(hidden_dims) - 1  # message passing layers
        self.gnn_layers = nn.ModuleList()
        self.layer_norms = nn.ModuleList()

        in_dim = hidden_dims[0]
        for layer_idx in range(self.n_layers):
            out_dim = hidden_dims[layer_idx + 1]
            # Each edge type gets its own GATv2Conv
            type_convs = nn.ModuleList([
                GATv2Conv(
                    in_dim,
                    out_dim // n_attention_heads,
                    heads=n_attention_heads,
                    dropout=dropout,
                    add_self_loops=True,
                    concat=True,  # concat heads → out_dim
                )
                for _ in range(self.n_edge_types)
            ])
            self.gnn_layers.append(type_convs)
            self.layer_norms.append(nn.LayerNorm(out_dim))
            in_dim = out_dim

        final_dim = hidden_dims[-1]

        # ── Heads ──
        self.decoder = ReconstructionDecoder(
            final_dim, decoder_hidden, n_node_features, dropout
        )
        self.anomaly_head = AnomalyHead(final_dim, anomaly_hidden, dropout)

        # ── Learnable mask token (for pretraining) ──
        self.mask_token = nn.Parameter(
            torch.randn(n_node_features) * 0.02
        )

        logger.info(
            "MethylationGNN: %d features → %s → decoder → anomaly head (%d params)",
            n_node_features,
            " → ".join(str(d) for d in hidden_dims),
            self.num_parameters,
        )

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_type: Optional[torch.Tensor] = None,
        return_attention: bool = False,
    ) -> Dict[str, torch.Tensor]:
        """
        Forward pass.

        Parameters
        ----------
        x : (N, d_features)
            Node feature matrix.
        edge_index : (2, E)
            Sparse adjacency (COO format).
        edge_type : (E,) or None
            Edge type IDs (0..n_edge_types-1). If None, uses homogeneous
            message passing (all edges treated as type 0).
        return_attention : bool
            If True, return attention weights (for interpretability).

        Returns
        -------
        dict with keys:
            - reconstructed : (N, d_features) reconstructed node features
            - anomaly_scores : (N, 1) per-node anomaly scores
            - node_embeddings : (N, d_final) GNN embeddings
            - graph_anomaly : (1,) aggregated graph-level anomaly score
            - attention_weights : list of tensor (optional)
        """
        # Input projection
        h = self.input_proj(x)
        attention_weights = []

        # Message passing
        for layer_idx, (type_convs, norm) in enumerate(
            zip(self.gnn_layers, self.layer_norms)
        ):
            if self.n_edge_types == 1 or edge_type is None:
                # Homogeneous: use first conv
                conv = type_convs[0]
                h_new = conv(h, edge_index)
            else:
                # Heterogeneous: mask edges by type and aggregate
                type_outputs = []
                for et in range(self.n_edge_types):
                    mask = edge_type == et
                    if mask.any():
                        ei = edge_index[:, mask]
                        if ei.numel() > 0:
                            out_et = type_convs[et](h, ei)
                            type_outputs.append(out_et)
                if type_outputs:
                    h_new = torch.stack(type_outputs).mean(dim=0)
                else:
                    # Fallback: identity (no edges of any type)
                    h_new = h

            # Residual + normalization
            if h.shape == h_new.shape:
                h = h + h_new
            else:
                h = h_new
            h = norm(h)
            h = F.relu(h)
            h = F.dropout(h, p=self.dropout_rate, training=self.training)

        # Decode
        x_reconstructed = self.decoder(h)
        anomaly_scores = self.anomaly_head(h)

        # Graph-level anomaly score (mean over nodes)
        graph_anomaly = anomaly_scores.mean(dim=0)  # (1,)

        output = {
            "reconstructed": x_reconstructed,
            "anomaly_scores": anomaly_scores,
            "node_embeddings": h,
            "graph_anomaly": graph_anomaly,
        }

        if return_attention:
            output["attention_weights"] = attention_weights

        return output

    def pretrain_forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_type: Optional[torch.Tensor] = None,
        mask_ratio: float = 0.30,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Masked node prediction forward pass (Phase 1 self-supervised).

        Randomly masks a fraction of node features and asks the model
        to reconstruct them from neighbour context. This is the core
        self-supervision signal that teaches the model what "normal"
        methylation patterns look like.

        Parameters
        ----------
        x : (N, d_features)
            Original node features.
        edge_index : (2, E)
            Graph edges.
        edge_type : (E,) or None
        mask_ratio : float
            Fraction of nodes to mask (0.30 = 30%).

        Returns
        -------
        reconstructed : (N, d_features)
            Full reconstruction (use with mask for loss).
        mask : (N,) bool
            Boolean mask of masked nodes.
        x_original : (N, d_features)
            Original features (for loss computation).
        anomaly_scores : (N, 1)
            Per-node anomaly scores.
        """
        x_original = x.clone()
        n_nodes = x.shape[0]
        device = x.device

        # Create random mask
        mask = torch.rand(n_nodes, device=device) < mask_ratio

        # Replace masked features with learnable mask token
        x_masked = x.clone()
        x_masked[mask] = self.mask_token.to(device)

        # Forward pass
        output = self.forward(x_masked, edge_index, edge_type)

        return (
            output["reconstructed"],
            mask,
            x_original,
            output["anomaly_scores"],
        )

    def predict(self, graph: "Data") -> float:
        """
        Produce field defect score for a single graph.

        Convenience method wrapping ``.eval()`` prediction to produce
        a scalar score in [0, 1].

        Parameters
        ----------
        graph : torch_geometric.data.Data
            Input graph with x, edge_index, edge_type.

        Returns
        -------
        field_defect_score : float
            Aggregated anomaly score. Higher → more likely field defect.
        """
        was_training = self.training
        self.eval()

        device = next(self.parameters()).device
        with torch.no_grad():
            x = graph.x.to(device)
            ei = graph.edge_index.to(device)
            et = graph.edge_type.to(device) if hasattr(graph, "edge_type") and graph.edge_type is not None else None

            output = self.forward(x, ei, et)
            score = float(output["graph_anomaly"].squeeze().cpu())

        if was_training:
            self.train()
        return score

    @property
    def num_parameters(self) -> int:
        """Total number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    @property
    def device_name(self) -> str:
        """Device the model is currently on."""
        try:
            return str(next(self.parameters()).device)
        except StopIteration:
            return "cpu"


# ── Loss function ───────────────────────────────────────────────

class FieldDefectLoss(nn.Module):
    """
    Combined loss for methylation GNN field defect detection.

    L_total = L_recon + λ_anomaly * L_anomaly + λ_temp * L_temporal

    Where:
    - L_recon: MSE between original and reconstructed node features
      (computed on masked nodes during pretraining, all nodes during finetune)
    - L_anomaly: BCE between graph-level anomaly score and cancer label
    - L_temporal: L1 smoothness of anomaly scores across timepoints
      (only used if longitudinal data available)
    """

    def __init__(
        self,
        lambda_anomaly: float = 0.5,
        lambda_temporal: float = 0.1,
    ):
        super().__init__()
        self.lambda_anomaly = lambda_anomaly
        self.lambda_temporal = lambda_temporal

    def forward(
        self,
        x_original: torch.Tensor,
        x_reconstructed: torch.Tensor,
        anomaly_scores: torch.Tensor,
        labels: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        prev_anomaly_scores: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, Dict[str, float]]:
        """
        Compute combined loss.

        Parameters
        ----------
        x_original : (N, d) or (B, N, d)
            Original node features.
        x_reconstructed : (N, d) or (B, N, d)
            Reconstructed node features.
        anomaly_scores : (N, 1) or (B*N, 1)
            Per-node anomaly scores.
        labels : (batch_size, 1) or (batch_size,)
            Cancer labels (0/1).
        mask : (N,) or (B, N) bool or None
            If provided, reconstruction loss is computed only on masked nodes
            (pretraining mode). If None, computed on all nodes (finetune mode).
        prev_anomaly_scores : (N, 1) or None
            Anomaly scores from previous timepoint (for temporal smoothness).

        Returns
        -------
        loss_total : scalar tensor
        metrics : dict
            Individual loss components (detached, for logging).
        """
        # Reconstruction loss
        if mask is not None and mask.any():
            loss_recon = F.mse_loss(
                x_reconstructed[mask], x_original[mask]
            )
        else:
            loss_recon = F.mse_loss(x_reconstructed, x_original)

        # Graph-level anomaly → BCE
        graph_anomaly = anomaly_scores.mean(dim=0)  # scalar
        # Ensure correct shapes for BCE
        ano_input = graph_anomaly.view(1)
        lab_input = labels.float().view(1)
        loss_anomaly = F.binary_cross_entropy(ano_input, lab_input)

        # Temporal consistency (optional)
        loss_temporal = torch.tensor(0.0, device=x_original.device)
        if prev_anomaly_scores is not None:
            loss_temporal = F.l1_loss(anomaly_scores, prev_anomaly_scores)

        # Total
        loss_total = (
            loss_recon
            + self.lambda_anomaly * loss_anomaly
            + self.lambda_temporal * loss_temporal
        )

        metrics = {
            "loss_recon": float(loss_recon.detach().cpu()),
            "loss_anomaly": float(loss_anomaly.detach().cpu()),
            "loss_temporal": float(loss_temporal.detach().cpu()),
            "loss_total": float(loss_total.detach().cpu()),
        }

        return loss_total, metrics


# ── Utility: Field defect score from reconstruction error ───────

def compute_field_defect_score(
    x_original: torch.Tensor,
    x_reconstructed: torch.Tensor,
    anomaly_scores: torch.Tensor,
) -> torch.Tensor:
    """
    Combine reconstruction error and anomaly scores into field defect score.

    s_field = Σ_i (x_i - x̂_i)² · a_i / Σ_i a_i

    This weights reconstruction error by anomaly importance, so regions
    that the anomaly head has learned are diagnostically relevant contribute
    more to the final score.

    Parameters
    ----------
    x_original : (N, d)
        Original node features.
    x_reconstructed : (N, d)
        Reconstructed node features.
    anomaly_scores : (N, 1)
        Per-node anomaly scores.

    Returns
    -------
    score : scalar tensor
        Field defect score (higher → more anomalous).
    """
    recon_error = ((x_original - x_reconstructed) ** 2).mean(dim=1, keepdim=True)
    weighted = recon_error * anomaly_scores
    score = weighted.sum() / (anomaly_scores.sum() + 1e-8)
    return score.squeeze()
