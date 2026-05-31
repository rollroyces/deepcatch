#!/usr/bin/env python3
"""
GNN Inference — Lightweight Prediction Pipeline
=================================================

Wraps a trained MethylationGNN model for production inference.
Supports single-sample and batch prediction with optional
interpretability (node-level anomaly maps, reconstruction error heatmaps).

Designed for:
- Integration with DeepCatch Stage 1 fusion layer
- API serving (via the inference server)
- Research analysis (interpretability modes)

Example
-------

.. code-block:: python

    from src.methylation_gnn import GNNInference, MethylationGNNPredictor
    from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder

    # Option 1: Low-level inference
    inf = GNNInference.load("checkpoints/finetune_best.pt")
    score, details = inf.predict(sample_graph, return_details=True)
    print(f"Field defect score: {score:.4f}")

    # Option 2: High-level predictor (graph builder + model)
    predictor = MethylationGNNPredictor(
        builder=RegulatoryGraphBuilder(n_nodes=5000),
        checkpoint="checkpoints/finetune_best.pt",
    )
    result = predictor.predict_sample(
        sample_name="sample_001",
        methylation_data={"beta_values": betas},
    )
"""

from __future__ import annotations

import logging
import os
from copy import deepcopy
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import torch
    import torch.nn.functional as F

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None

try:
    from torch_geometric.data import Data

    _HAS_PYG = True
except ImportError:
    _HAS_PYG = False
    Data = None

from .config import GNNConfig
from .gnn_model import MethylationGNN, compute_field_defect_score

logger = logging.getLogger(__name__)


class GNNInference:
    """
    Wraps a trained MethylationGNN for inference-only usage.

    Handles device placement, model loading, and prediction with
    optional interpretability outputs.

    Parameters
    ----------
    model : MethylationGNN
        Trained model (already loaded with weights).
    config : GNNConfig or None
        Configuration (defaults to model's config if available).
    device : str or None
        Device override. If None, uses model's current device.
    """

    def __init__(
        self,
        model: MethylationGNN,
        config: Optional[GNNConfig] = None,
        device: Optional[str] = None,
    ):
        if not _HAS_PYG:
            raise ImportError(
                "PyTorch Geometric is required for GNNInference. "
                "Install with: pip install torch_geometric"
            )

        self.config = config or GNNConfig()
        self.device = torch.device(device or self.config.device)
        self.model = model.to(self.device)
        self.model.eval()

        logger.info(
            "GNNInference ready on %s | %d params",
            self.device, self.model.num_parameters,
        )

    @classmethod
    def load(
        cls,
        checkpoint_path: str,
        device: Optional[str] = None,
        **model_kwargs,
    ) -> "GNNInference":
        """
        Load a trained model from checkpoint.

        Parameters
        ----------
        checkpoint_path : str
            Path to a .pt checkpoint saved by GNNTrainer.
        device : str or None
            Device to load model on. Auto-detected if None.
        **model_kwargs
            Override model architecture parameters (if checkpoint
            doesn't contain config).

        Returns
        -------
        GNNInference
            Initialized inference wrapper.

        Raises
        ------
        FileNotFoundError
            If checkpoint doesn't exist.
        """
        if not os.path.exists(checkpoint_path):
            raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

        if device is None:
            device = GNNConfig._detect_device()

        checkpoint = torch.load(
            checkpoint_path, map_location=device, weights_only=False
        )

        # Try to extract config from checkpoint
        if "config" in checkpoint:
            cfg_dict = checkpoint["config"]
            config = GNNConfig(**{
                k: v for k, v in cfg_dict.items()
                if k in GNNConfig.__dataclass_fields__
            })
            logger.info("Loaded config from checkpoint")
        else:
            config = GNNConfig()
            logger.warning("No config in checkpoint; using defaults")

        # Build model
        model = MethylationGNN(
            n_node_features=config.n_node_features,
            hidden_dims=config.hidden_dims,
            n_edge_types=config.n_edge_types,
            n_attention_heads=config.n_attention_heads,
            decoder_hidden=config.decoder_hidden,
            anomaly_hidden=config.anomaly_hidden,
            dropout=config.dropout,
        )

        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(torch.device(device))

        logger.info(
            "Loaded model from %s (epoch=%d)",
            checkpoint_path,
            checkpoint.get("epoch", 0),
        )

        return cls(model, config, device)

    def _to_device(self, graph: "Data") -> "Data":
        """Move graph tensors to inference device."""
        g = graph.clone()
        for key, val in g.items():
            if isinstance(val, torch.Tensor):
                g[key] = val.to(self.device)
        return g

    @torch.no_grad()
    def predict(
        self,
        graph: "Data",
        return_details: bool = False,
    ) -> Union[float, Dict[str, Any]]:
        """
        Predict field defect score for a single graph.

        Parameters
        ----------
        graph : Data
            PyG graph (from RegulatoryGraphBuilder.build_graph).
        return_details : bool
            If True, return dict with full prediction details.

        Returns
        -------
        float or dict
            Scalar field_defect_score ∈ [0, 1] if return_details=False.
            Full dict with node-level scores if return_details=True.
        """
        graph = self._to_device(graph)

        output = self.model(
            graph.x, graph.edge_index,
            edge_type=getattr(graph, "edge_type", None),
        )

        anomaly_score = float(output["graph_anomaly"].squeeze().cpu())
        field_score = float(
            compute_field_defect_score(
                graph.x,
                output["reconstructed"],
                output["anomaly_scores"],
            ).cpu()
        )

        if not return_details:
            return field_score

        # Full details
        recon_error = ((graph.x - output["reconstructed"]) ** 2).mean(dim=1)
        node_scores = output["anomaly_scores"].squeeze(-1)

        return {
            "field_defect_score": field_score,
            "graph_anomaly_score": anomaly_score,
            "node_anomaly_scores": node_scores.cpu().numpy(),
            "node_recon_errors": recon_error.cpu().numpy(),
            "reconstructed_features": output["reconstructed"].cpu().numpy(),
            "node_embeddings": output["node_embeddings"].cpu().numpy(),
            "top_anomalous_nodes": torch.topk(node_scores, min(10, len(node_scores))).indices.cpu().numpy().tolist(),
        }

    def predict_batch(
        self,
        graphs: List["Data"],
        return_details: bool = False,
    ) -> Union[List[float], List[Dict[str, Any]]]:
        """
        Predict field defect scores for multiple graphs.

        Parameters
        ----------
        graphs : list of Data
            List of PyG graphs.
        return_details : bool
            If True, return list of detail dicts.

        Returns
        -------
        list of float or list of dict
        """
        results = []
        for graph in graphs:
            result = self.predict(graph, return_details=return_details)
            results.append(result)
        return results


class MethylationGNNPredictor:
    """
    High-level predictor combining graph builder + trained model.

    This is the main entry point for integration with the DeepCatch
    CET pipeline. It takes raw methylation data, builds a graph,
    runs GNN inference, and returns a field defect score ready for
    fusion with fragmentomics, CNV, and serological modalities.

    Parameters
    ----------
    builder : RegulatoryGraphBuilder
        Configured graph builder (loaded with reference regions).
    checkpoint : str
        Path to trained model checkpoint.
    device : str or None
        Inference device (auto-detect if None).
    """

    def __init__(
        self,
        builder: "RegulatoryGraphBuilder",  # type: ignore
        checkpoint: str,
        device: Optional[str] = None,
    ):
        from .graph_builder import RegulatoryGraphBuilder as _RB

        self.builder: _RB = builder
        self.inference = GNNInference.load(checkpoint, device=device)

    def predict_sample(
        self,
        sample_name: str,
        methylation_data: Dict[str, np.ndarray],
        cfDNA_coverage: Optional[np.ndarray] = None,
        chromatin_data: Optional[Dict[str, np.ndarray]] = None,
        fragmentomics_data: Optional[Dict[str, np.ndarray]] = None,
        return_details: bool = False,
    ) -> Union[float, Dict[str, Any]]:
        """
        End-to-end prediction for one sample.

        Parameters
        ----------
        sample_name : str
            Sample identifier.
        methylation_data : dict
            Methylation features per region (minimum: 'beta_values').
        cfDNA_coverage : (n_nodes,) array or None
        chromatin_data : dict or None
        fragmentomics_data : dict or None
        return_details : bool
            Return full prediction details or just score.

        Returns
        -------
        float or dict
        """
        graph = self.builder.build_graph(
            sample_name=sample_name,
            methylation_data=methylation_data,
            cfDNA_coverage=cfDNA_coverage,
            chromatin_data=chromatin_data,
            fragmentomics_data=fragmentomics_data,
        )

        return self.inference.predict(graph, return_details=return_details)

    def predict_batch(
        self,
        samples: List[Dict[str, Any]],
        return_details: bool = False,
    ) -> Union[List[float], List[Dict[str, Any]]]:
        """
        End-to-end prediction for multiple samples.

        Parameters
        ----------
        samples : list of dict
            Each dict contains sample parameters (sample_name,
            methylation_data, etc.).
        return_details : bool

        Returns
        -------
        list of float or list of dict
        """
        results = []
        for i, sample in enumerate(samples):
            if isinstance(sample, dict) and "sample_name" not in sample:
                sample["sample_name"] = f"sample_{i}"
            result = self.predict_sample(**sample, return_details=return_details)
            results.append(result)
        return results
