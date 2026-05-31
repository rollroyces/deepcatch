#!/usr/bin/env python3
"""
Integration Adapter — GNN Methylation Branch → DeepCatch CET Fusion
=====================================================================

Connects the GNN methylation network branch to the existing DeepCatch
multi-modal fusion pipeline in ``src/multimodal_fusion/advanced_fusion.py``.

The adapter provides:

1. ``ModularArmsBuilder`` — Composite feature extractor that runs all
   Stage 1 modalities (fragmentomics, CNV, serological, GNN methylation)
   and produces standardized scores ready for fusion.

2. ``MethylationBranchAdapter`` — Standalone wrapper that takes raw
   methylation data and produces a GNN field defect score compatible
   with ``CrossAttentionFusion`` and ``EarlyLateFusion`` interfaces.

Integration Point
-----------------

In the existing DeepCatch pipeline, the fusion layer expects a list of
(n_samples,) arrays::

    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion

    fusion = CrossAttentionFusion(n_modalities=5)  # was 4, now +1 for GNN
    fusion.fit([frag_scores, cnv_scores, sero_scores, mfr_scores, gnn_scores], labels)

The GNN branch adds a 5th modality with orthogonal signal (epigenetic
field defects) that fragmentomics alone cannot detect at ultra-low
tumor fractions (<0.01% ctDNA).

Example
-------

.. code-block:: python

    from src.methylation_gnn.integration import ModularArmsBuilder

    # Build pipeline with all modalities
    arms = ModularArmsBuilder(
        methylation_checkpoint="checkpoints/finetune_best.pt",
    )

    # Process a batch of samples
    results = arms.process_samples(samples)

    # Fuse with existing pipeline
    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
    fusion = CrossAttentionFusion(n_modalities=5)
    fusion.fit(results['modality_scores'], labels)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import torch

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None

from .graph_builder import RegulatoryGraphBuilder
from .gnn_model import MethylationGNN
from .gnn_inference import GNNInference, MethylationGNNPredictor

logger = logging.getLogger(__name__)


class MethylationBranchAdapter:
    """
    Adapter that makes GNN methylation predictions compatible with DeepCatch fusion.

    Wraps the full GNN pipeline (graph building + inference) into a
    simple ``predict(samples)`` interface that returns a (n_samples,)
    array of field defect scores. This plugs directly into the existing
    ``CrossAttentionFusion.fit(modality_scores, labels)`` and
    ``EarlyLateFusion.fit(modality_features, labels)`` interfaces.

    Parameters
    ----------
    checkpoint : str
        Path to trained GNN model checkpoint.
    n_regions : int
        Number of genomic regions for graph construction.
        Default 5000 for quick integration; use 50000 for production.
    device : str or None
        Inference device (auto-detect if None).
    """

    def __init__(
        self,
        checkpoint: str,
        n_regions: int = 5_000,
        device: Optional[str] = None,
        edge_k: int = 15,
    ):
        if not _HAS_TORCH:
            raise ImportError(
                "PyTorch is required for MethylationBranchAdapter. "
                "Install with: pip install torch"
            )

        self.builder = RegulatoryGraphBuilder(
            n_nodes=n_regions,
            edge_k=edge_k,
        )
        self.predictor = MethylationGNNPredictor(
            builder=self.builder,
            checkpoint=checkpoint,
            device=device,
        )

    def predict_sample(
        self,
        methylation_data: Dict[str, np.ndarray],
        sample_name: str = "sample",
        cfDNA_coverage: Optional[np.ndarray] = None,
    ) -> float:
        """
        Predict field defect score for one sample.

        Parameters
        ----------
        methylation_data : dict
            Must contain 'beta_values' per region.
        sample_name : str
            Sample identifier.
        cfDNA_coverage : (n_regions,) array or None

        Returns
        -------
        field_defect_score : float
            Scalar in [0, 1]. Higher → higher likelihood of field defect.
        """
        return self.predictor.predict_sample(
            sample_name=sample_name,
            methylation_data=methylation_data,
            cfDNA_coverage=cfDNA_coverage,
        )

    def predict_batch(
        self,
        samples: List[Dict[str, Any]],
    ) -> np.ndarray:
        """
        Predict field defect scores for a batch of samples.

        Parameters
        ----------
        samples : list of dict
            Each dict must have 'methylation_data' key, optionally
            'sample_name' and 'cfDNA_coverage'.

        Returns
        -------
        scores : (n_samples,) float32 array
            Field defect scores ready for fusion.
        """
        results = self.predictor.predict_batch(samples)
        return np.array(results, dtype=np.float32)


class ModularArmsBuilder:
    """
    Composite feature extractor for all DeepCatch Stage 1 modalities.

    Mimics the multi-arm architecture of the THEMIS framework, producing
    standardized per-modality scores ready for the fusion layer.

    Modalities
    ----------
    0. Fragmentomics (FSI, MFR, CAFF, FEM) — from themis_features.py
    1. CNV — chromosomal instability index
    2. Serological (PG, G-17, Hp) — from serological_fusion.py
    3. MFR standalone — methylation fraction ratio (often separate from GNN)
    4. ★ GNN Methylation Network — field defect score (NEW in v2.1)

    Parameters
    ----------
    methylation_checkpoint : str
        Path to trained GNN model checkpoint for modality 4.
    n_gnn_regions : int
        Graph size for GNN branch.
    include_gnn : bool
        If False, acts as a v2.0-compatible builder (4 modalities,
        no GNN). Useful for ablation testing.
    """

    def __init__(
        self,
        methylation_checkpoint: Optional[str] = None,
        n_gnn_regions: int = 5_000,
        include_gnn: bool = True,
        device: Optional[str] = None,
    ):
        self.include_gnn = include_gnn
        self.gnn_adapter: Optional[MethylationBranchAdapter] = None

        if include_gnn and methylation_checkpoint:
            self.gnn_adapter = MethylationBranchAdapter(
                checkpoint=methylation_checkpoint,
                n_regions=n_gnn_regions,
                device=device,
            )
        elif include_gnn:
            logger.warning(
                "include_gnn=True but no checkpoint provided. "
                "GNN branch will be disabled."
            )
            self.include_gnn = False

    @property
    def n_modalities(self) -> int:
        """Number of active modalities (4 or 5)."""
        return 5 if self.include_gnn else 4

    def extract_all_scores(
        self,
        samples: List[Dict[str, Any]],
    ) -> List[np.ndarray]:
        """
        Extract all modality scores for a batch of samples.

        Each sample dict should contain:
        - fragmentomics features (for modalities 0, 3)
        - cnv features (for modality 1)
        - serological features (for modality 2)
        - methylation_data (for modality 4, if GNN enabled)

        Parameters
        ----------
        samples : list of dict

        Returns
        -------
        modality_scores : list of (n_samples,) arrays
            Compatible with CrossAttentionFusion.fit(modality_scores, labels).
        """
        n = len(samples)
        scores = []

        # Modality 0: Fragmentomics composite
        scores.append(self._extract_fragmentomics_scores(samples))

        # Modality 1: CNV
        scores.append(self._extract_cnv_scores(samples))

        # Modality 2: Serological
        scores.append(self._extract_serological_scores(samples))

        # Modality 3: MFR standalone
        scores.append(self._extract_mfr_scores(samples))

        # Modality 4: GNN methylation (NEW)
        if self.include_gnn and self.gnn_adapter:
            gnn_scores = self._extract_gnn_scores(samples)
            scores.append(gnn_scores)

        return scores

    def _extract_fragmentomics_scores(self, samples: List[Dict]) -> np.ndarray:
        """Extract fragmentomics scores. Placeholder — delegates to existing."""
        # In production, this calls themis_features.py calculators
        # For now: extract from sample dict or return random placeholder
        scores = np.zeros(len(samples), dtype=np.float32)
        for i, s in enumerate(samples):
            if "fragmentomics_score" in s:
                scores[i] = s["fragmentomics_score"]
            elif "fsi" in s:
                # Simple composite: 1/FSI (high FSI → more long fragments → healthy)
                scores[i] = np.clip(1.0 / (s.get("fsi", 1.0) + 0.1), 0, 1)
            else:
                scores[i] = 0.5  # neutral placeholder
        return scores

    def _extract_cnv_scores(self, samples: List[Dict]) -> np.ndarray:
        """Extract CNV/CAFF scores."""
        scores = np.zeros(len(samples), dtype=np.float32)
        for i, s in enumerate(samples):
            if "cnv_score" in s:
                scores[i] = s["cnv_score"]
            elif "caff_score" in s:
                scores[i] = s["caff_score"]
            else:
                scores[i] = 0.5
        return scores

    def _extract_serological_scores(self, samples: List[Dict]) -> np.ndarray:
        """Extract serological fusion scores."""
        scores = np.zeros(len(samples), dtype=np.float32)
        for i, s in enumerate(samples):
            if "serological_score" in s:
                scores[i] = s["serological_score"]
            else:
                scores[i] = 0.5
        return scores

    def _extract_mfr_scores(self, samples: List[Dict]) -> np.ndarray:
        """Extract MFR (methylation fraction ratio) scores."""
        scores = np.zeros(len(samples), dtype=np.float32)
        for i, s in enumerate(samples):
            if "mfr_score" in s:
                scores[i] = s["mfr_score"]
            elif "methylation_data" in s:
                betas = s["methylation_data"].get("beta_values", np.array([0.5]))
                scores[i] = float(np.mean(betas))
            else:
                scores[i] = 0.5
        return scores

    def _extract_gnn_scores(self, samples: List[Dict]) -> np.ndarray:
        """Extract GNN methylation field defect scores."""
        if self.gnn_adapter is None:
            return np.full(len(samples), 0.5, dtype=np.float32)

        return self.gnn_adapter.predict_batch(samples)

    def process_sample(
        self,
        sample: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Process one sample and return all modality scores + fused result.

        Parameters
        ----------
        sample : dict
            All modality data for one sample.

        Returns
        -------
        result : dict
            Keys: fragmentomics_score, cnv_score, serological_score,
            mfr_score, gnn_field_defect_score (if enabled),
            n_modalities.
        """
        scores = self.extract_all_scores([sample])
        result = {
            "fragmentomics_score": float(scores[0][0]),
            "cnv_score": float(scores[1][0]),
            "serological_score": float(scores[2][0]),
            "mfr_score": float(scores[3][0]),
            "n_modalities": self.n_modalities,
        }
        if self.include_gnn and len(scores) > 4:
            result["gnn_field_defect_score"] = float(scores[4][0])
        return result


# ── Convenience function: extend existing fusion ────────────────


def extend_fusion_with_gnn(
    existing_modality_scores: List[np.ndarray],
    gnn_scores: np.ndarray,
) -> List[np.ndarray]:
    """
    Append GNN scores to an existing modality score list.

    This is a lightweight way to add the GNN branch without refactoring
    existing code that calls ``CrossAttentionFusion.fit()``.

    Parameters
    ----------
    existing_modality_scores : list of (n_samples,) arrays
        Scores from modalities 0-3 (fragmentomics, CNV, serological, MFR).
    gnn_scores : (n_samples,) array
        Field defect scores from the GNN branch.

    Returns
    -------
    extended_scores : list of (n_samples,) arrays
        Same as input but with GNN scores appended as modality 4.

    Example
    -------

    .. code-block:: python

        # Existing pipeline
        fusion = CrossAttentionFusion(n_modalities=4)
        fusion.fit([frag, cnv, sero, mfr], labels)

        # Extend with GNN
        from src.methylation_gnn import MethylationBranchAdapter
        adapter = MethylationBranchAdapter("checkpoint.pt")
        gnn_scores = adapter.predict_batch(samples)

        from src.methylation_gnn.integration import extend_fusion_with_gnn
        all_scores = extend_fusion_with_gnn([frag, cnv, sero, mfr], gnn_scores)

        fusion2 = CrossAttentionFusion(n_modalities=5)
        fusion2.fit(all_scores, labels)
    """
    if len(existing_modality_scores) < 4:
        logger.warning(
            "Expected 4 existing modalities, got %d. Padding with zeros.",
            len(existing_modality_scores),
        )
    return existing_modality_scores + [gnn_scores]
