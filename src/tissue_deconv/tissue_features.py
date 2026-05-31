#!/usr/bin/env python3
"""
Tissue Deconvolution Feature Extraction
=========================================

Extracts fixed-size feature vectors from tissue deconvolution results
for fusion with other DeepCatch modalities in Stage 1 (Capture).

Features capture:
- Tissue composition of cfDNA (inferred origin of cell-free DNA)
- Dominant tissue signals (which organ is shedding)
- Tissue diversity metrics (entropy, concentration)
- Cancer-relevant tissue indices
- Ratio-based diagnostics (two-tissue comparisons)

Works even without a trained model: returns zero-filled features
as a graceful fallback.

Reference
---------
Li et al. (2023) PNAS — cfSort tissue deconvolution framework.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .config import (
    SUPPORTED_TISSUES,
    CANCER_RELEVANT_TISSUES,
    BLOOD_DERIVED_TISSUES,
    TissueDeconvConfig,
    DEFAULT_CONFIG,
)

logger = logging.getLogger(__name__)


class TissueDeconvolutionFeatures:
    """
    Extract interpretable scalar features from tissue deconvolution results.

    Designed to plug into CrossAttentionFusion and EarlyLateFusion
    as additional modality features (~20-30 scalars per sample).

    Parameters
    ----------
    config : TissueDeconvConfig
        Configuration for feature extraction thresholds.
    """

    def __init__(self, config: Optional[TissueDeconvConfig] = None):
        self.config = config if config is not None else DEFAULT_CONFIG
        self.detection_limit = self.config.detection_limit

        # Pre-compute tissue indices for fast lookups
        self._tissue_idx: Dict[str, int] = {
            t: i for i, t in enumerate(SUPPORTED_TISSUES)
        }
        self._blood_idx: List[int] = [
            self._tissue_idx[t] for t in BLOOD_DERIVED_TISSUES
            if t in self._tissue_idx
        ]
        self._cancer_idx: List[int] = [
            self._tissue_idx[t] for t in CANCER_RELEVANT_TISSUES
            if t in self._tissue_idx
        ]

    @property
    def n_features(self) -> int:
        """Total number of scalar features extracted per sample."""
        return 24

    def feature_names(self) -> List[str]:
        """Names of all extracted features in order."""
        return [
            # Top tissue signals (5)
            "tissue_top1_fraction",
            "tissue_top2_fraction",
            "tissue_top3_fraction",
            "tissue_top4_fraction",
            "tissue_top5_fraction",
            # Tissue identity (5)  —  which tissues
            "tissue_top1_is_blood",
            "tissue_top1_is_cancer",
            "tissue_top2_is_blood",
            "tissue_top2_is_cancer",
            "tissue_top3_is_blood",
            # Diversity metrics (4)
            "tissue_entropy",
            "tissue_n_active",
            "tissue_diversity_score",
            "tissue_concentration_gini",
            # Ratio diagnostics (4)
            "tissue_top2_ratio",
            "tissue_cancer_blood_ratio",
            "tissue_max_cancer_fraction",
            "tissue_gastrointestinal_sum",
            # Organ-specific (6)
            "tissue_liver_fraction",
            "tissue_lung_fraction",
            "tissue_colon_fraction",
            "tissue_pancreas_fraction",
            "tissue_breast_fraction",
            "tissue_prostate_fraction",
        ]

    # ── Main extraction ──────────────────────────────────────────

    def extract(
        self,
        tissue_fractions: np.ndarray,
    ) -> Dict[str, float]:
        """
        Extract features from tissue deconvolution fractions.

        Parameters
        ----------
        tissue_fractions : (n_tissues,) array
            Predicted tissue fractions from deconvolution model.
            Should sum to ~1.0.

        Returns
        -------
        features : dict
            Dictionary of named scalar features for fusion.
        """
        if tissue_fractions.ndim != 1:
            raise ValueError(
                f"Expected 1D array (n_tissues,), got shape {tissue_fractions.shape}"
            )
        if len(tissue_fractions) != len(SUPPORTED_TISSUES):
            raise ValueError(
                f"Expected {len(SUPPORTED_TISSUES)} tissues, "
                f"got {len(tissue_fractions)}"
            )

        frac = tissue_fractions.copy()

        # Sort tissues by fraction (descending)
        sorted_idx = np.argsort(frac)[::-1]
        sorted_frac = frac[sorted_idx]

        features: Dict[str, float] = {}

        # ── Top tissue fractions ─────────────────────────────────
        for rank in range(5):
            key = f"tissue_top{rank+1}_fraction"
            features[key] = float(sorted_frac[rank]) if rank < len(sorted_frac) else 0.0

        # ── Top tissue identity flags ────────────────────────────
        for rank in range(3):
            if rank < len(sorted_idx):
                tissue_idx = int(sorted_idx[rank])
                tissue = SUPPORTED_TISSUES[tissue_idx]
                features[f"tissue_top{rank+1}_is_blood"] = float(
                    tissue_idx in self._blood_idx
                )
                features[f"tissue_top{rank+1}_is_cancer"] = float(
                    tissue_idx in self._cancer_idx
                )
            else:
                features[f"tissue_top{rank+1}_is_blood"] = 0.0
                features[f"tissue_top{rank+1}_is_cancer"] = 0.0

        # ── Diversity metrics ────────────────────────────────────
        eps = 1e-7
        entropy = -np.sum(frac * np.log(frac + eps))
        max_entropy = np.log(len(SUPPORTED_TISSUES))
        n_active = int(np.sum(frac > self.detection_limit))

        features["tissue_entropy"] = float(entropy)
        features["tissue_n_active"] = float(n_active)
        features["tissue_diversity_score"] = float(1.0 - entropy / max_entropy)

        # Gini coefficient of tissue concentration
        sorted_asc = np.sort(frac)
        n = len(frac)
        index = np.arange(1, n + 1)
        gini = (2 * np.sum(index * sorted_asc)) / (n * np.sum(sorted_asc)) - (n + 1) / n
        features["tissue_concentration_gini"] = float(np.clip(gini, 0, 1))

        # ── Ratio diagnostics ────────────────────────────────────
        top1 = float(sorted_frac[0]) if len(sorted_frac) > 0 else 0.0
        top2 = float(sorted_frac[1]) if len(sorted_frac) > 1 else 0.0

        # Top2/Top1 ratio: when close to 1, two tissues have similar
        # death rates → possible co-pathology or specific organ damage
        features["tissue_top2_ratio"] = float(top2 / (top1 + eps))

        # Cancer vs blood ratio: elevated when tumor-derived cfDNA
        # dominates over normal blood cell turnover
        cancer_sum = float(np.sum(frac[self._cancer_idx]))
        blood_sum = float(np.sum(frac[self._blood_idx]))
        features["tissue_cancer_blood_ratio"] = float(cancer_sum / (blood_sum + eps))
        features["tissue_max_cancer_fraction"] = float(
            np.max(frac[self._cancer_idx]) if self._cancer_idx else 0.0
        )

        # Gastrointestinal sum (clinically relevant for GI cancers)
        gi_tissues = ["Colon", "Stomach", "Pancreas", "Liver", "Esophagus", "Small Intestine"]
        gi_sum = sum(
            frac[self._tissue_idx[t]]
            for t in gi_tissues
            if t in self._tissue_idx
        )
        features["tissue_gastrointestinal_sum"] = float(gi_sum)

        # ── Organ-specific fractions ─────────────────────────────
        organ_map = {
            "tissue_liver_fraction": "Liver",
            "tissue_lung_fraction": "Lung",
            "tissue_colon_fraction": "Colon",
            "tissue_pancreas_fraction": "Pancreas",
            "tissue_breast_fraction": "Breast",
            "tissue_prostate_fraction": "Prostate",
        }
        for key, tissue in organ_map.items():
            if tissue in self._tissue_idx:
                features[key] = float(frac[self._tissue_idx[tissue]])
            else:
                features[key] = 0.0

        return features

    def extract_batch(
        self,
        tissue_fractions: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """
        Extract features for a batch of samples.

        Parameters
        ----------
        tissue_fractions : (n_samples, n_tissues) array
            Predicted tissue fractions for multiple samples.

        Returns
        -------
        features : dict
            {feature_name: (n_samples,) array}.
        """
        if tissue_fractions.ndim != 2:
            raise ValueError(
                f"Expected 2D array (n_samples, n_tissues), "
                f"got shape {tissue_fractions.shape}"
            )

        n_samples = tissue_fractions.shape[0]
        feat_lists: Dict[str, List[float]] = {
            name: [] for name in self.feature_names()
        }

        for i in range(n_samples):
            sample_feats = self.extract(tissue_fractions[i])
            for name in self.feature_names():
                feat_lists[name].append(sample_feats[name])

        return {
            name: np.array(values, dtype=np.float32)
            for name, values in feat_lists.items()
        }

    def extract_from_mixture(
        self,
        methylation_data: np.ndarray,
        fractions: Optional[np.ndarray] = None,
    ) -> Dict[str, float]:
        """
        Extract features directly from methylation data + optional
        pre-computed fractions.

        If fractions is None, returns zero-filled features (graceful
        fallback when no deconvolution model is available).

        Parameters
        ----------
        methylation_data : (n_cpg_features,) array
            Methylation beta values.
        fractions : (n_tissues,) array, optional
            Pre-computed tissue fractions.

        Returns
        -------
        features : dict
        """
        if fractions is not None:
            return self.extract(fractions)

        # Graceful fallback: zero features
        logger.debug(
            "No tissue fractions provided; returning zero features "
            "(deconvolution model not available)"
        )
        return {name: 0.0 for name in self.feature_names()}

    # ── Summary statistics (for reporting) ───────────────────────

    def top_tissues(
        self,
        tissue_fractions: np.ndarray,
        n: int = 5,
        min_fraction: float = 0.001,
    ) -> List[Tuple[str, float]]:
        """
        Get top N tissues with their fractions.

        Parameters
        ----------
        tissue_fractions : (n_tissues,) array
        n : int
            Number of top tissues.
        min_fraction : float
            Minimum fraction to include.

        Returns
        -------
        top_tissues : list of (tissue_name, fraction) tuples
        """
        sorted_idx = np.argsort(tissue_fractions)[::-1]
        result = []
        for idx in sorted_idx[:n]:
            frac = float(tissue_fractions[idx])
            if frac >= min_fraction:
                result.append((SUPPORTED_TISSUES[idx], frac))
        return result

    def tissue_distribution_summary(
        self,
        tissue_fractions: np.ndarray,
    ) -> str:
        """
        Human-readable summary of tissue deconvolution results.

        Parameters
        ----------
        tissue_fractions : (n_tissues,) array

        Returns
        -------
        summary : str
        """
        blood_sum = float(np.sum(tissue_fractions[self._blood_idx]))
        cancer_sum = float(np.sum(tissue_fractions[self._cancer_idx]))
        top_tissues = self.top_tissues(tissue_fractions, n=5)

        lines = [
            f"Blood-derived: {blood_sum:.1%}",
            f"Cancer-relevant: {cancer_sum:.1%}",
            "Top tissues:",
        ]
        for tissue, frac in top_tissues:
            lines.append(f"  {tissue}: {frac:.3%}")

        return "\n".join(lines)
