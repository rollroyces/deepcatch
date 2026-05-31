#!/usr/bin/env python3
"""
DEConv Integration — Tissue Deconvolution → DeepCatch CET Fusion
===================================================================

Connects the tissue deconvolution module (cfSort-style DNN head) to the
existing DeepCatch multi-modal fusion pipeline in
``src/multimodal_fusion/advanced_fusion.py``.

The adapter provides:

1. ``DEConvIntegration`` — Wraps the full tissue deconvolution pipeline
   (atlas + model + features) into a single interface compatible with
   ``CrossAttentionFusion.fit(modality_scores, labels)``.

2. ``extract_all(sample)`` — Extracts all tissue deconvolution features
   for a single sample in a dictionary format.

3. ``to_modality(sample)`` — Produces a single scalar score (max abnormal
   tissue fraction) for simple fusion scenarios.

Integration Point
-----------------

In the existing DeepCatch pipeline with 5 modalities (fragmentomics, CNV,
serological, MFR, GNN methylation), tissue deconvolution adds a **6th
modality** from cfDNA methylation:

.. code-block:: python

    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
    from src.tissue_deconv.integration import DEConvIntegration

    deconv = DEConvIntegration()
    deconv.fit_synthetic(n_samples=2000)

    # Now 6 modalities
    fusion = CrossAttentionFusion(n_modalities=6)
    fusion.fit(
        [frag, cnv, sero, mfr, gnn, deconv_scores],
        labels,
    )

Example
-------

.. code-block:: python

    from src.tissue_deconv.integration import DEConvIntegration

    integ = DEConvIntegration(
        checkpoint="checkpoints/tissue_deconv/ensemble.pt",
    )

    # Process one sample
    tissue_scores = integ.predict_tissue_fractions(methylation_data)
    features = integ.extract_all(sample)
    modality_score = integ.to_modality(sample)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .config import (
    TissueDeconvConfig,
    DEFAULT_CONFIG,
    SUPPORTED_TISSUES,
    CANCER_RELEVANT_TISSUES,
    BLOOD_DERIVED_TISSUES,
)
from .tissue_atlas import TissueAtlas
from .tissue_features import TissueDeconvolutionFeatures

logger = logging.getLogger(__name__)


class DEConvIntegration:
    """
    Integration adapter for tissue deconvolution in DeepCatch Stage 1.

    Wraps the tissue atlas, deconvolution model, and feature extraction
    into a unified interface compatible with the fusion layer.

    Supports two modes:
    - **Real mode**: With trained model checkpoint, full deconvolution.
    - **Fallback mode**: Without model, returns zeros/neutral scores
      (graceful degradation).

    Parameters
    ----------
    checkpoint : str, optional
        Path to trained TissueDeconvTrainer checkpoint.
    config : TissueDeconvConfig, optional
        Configuration (auto-detects from checkpoint if loaded).
    use_ensemble : bool
        If True, expect ensemble checkpoint; else single model.
    device : str, optional
        Device override ('cpu', 'mps', 'cuda', 'auto').
    """

    def __init__(
        self,
        checkpoint: Optional[str] = None,
        config: Optional[TissueDeconvConfig] = None,
        use_ensemble: bool = True,
        device: Optional[str] = None,
    ):
        self.config = config if config is not None else DEFAULT_CONFIG
        self.device_str = device or self.config.device
        self.use_ensemble = use_ensemble

        # Submodules
        self.atlas = TissueAtlas(
            n_cpg_features=self.config.n_cpg_features,
        )
        self.feature_extractor = TissueDeconvolutionFeatures(self.config)

        # Deconvolution model (lazy-loaded from checkpoint or trained)
        self._trainer = None
        self._ensemble = None
        self._has_model = False

        # Load atlas (or generate synthetic)
        self.atlas.load_reference()

        # Load checkpoint if provided
        if checkpoint is not None:
            self._load_checkpoint(checkpoint)

    # ── Model loading ────────────────────────────────────────────

    def _load_checkpoint(self, checkpoint: str):
        """Load a trained model checkpoint."""
        try:
            import torch
        except ImportError:
            logger.warning(
                "PyTorch not available; cannot load checkpoint. "
                "Falling back to zero-feature mode."
            )
            return

        try:
            from .trainer import TissueDeconvTrainer

            self._trainer = TissueDeconvTrainer(
                config=self.config,
                device=self.device_str,
            )
            self._trainer.load_checkpoint(checkpoint)
            self._has_model = True
            logger.info(
                "Loaded tissue deconvolution model from %s", checkpoint
            )
        except (FileNotFoundError, RuntimeError, KeyError) as e:
            logger.warning(
                "Failed to load checkpoint '%s': %s. "
                "Falling back to zero-feature mode.",
                checkpoint, e,
            )
            self._has_model = False

    def fit_synthetic(
        self,
        n_samples: int = 2000,
        n_epochs: int = 50,
        verbose: bool = False,
    ) -> None:
        """
        Train the deconvolution model on synthetic data.

        This creates a fully functional model without requiring
        real cfSort reference data.

        Parameters
        ----------
        n_samples : int
            Number of synthetic training samples.
        n_epochs : int
            Training epochs.
        verbose : bool
            Print progress.
        """
        try:
            import torch
            from .trainer import TissueDeconvTrainer
        except ImportError:
            logger.warning("PyTorch not available; cannot train model.")
            return

        self._trainer = TissueDeconvTrainer(
            config=self.config,
            device=self.device_str,
        )
        self._trainer.atlas = self.atlas
        self._ensemble = self._trainer.quick_train(
            n_samples=n_samples,
            n_epochs=n_epochs,
            verbose=verbose,
        )
        self._has_model = True

    # ── Deconvolution prediction ─────────────────────────────────

    def predict_tissue_fractions(
        self,
        methylation_data: np.ndarray,
    ) -> np.ndarray:
        """
        Predict tissue fractions from methylation beta values.

        Parameters
        ----------
        methylation_data : (n_cpg_features,) or (n_samples, n_cpg_features) array
            Methylation beta values at tissue-discriminative CpGs.

        Returns
        -------
        fractions : (n_tissues,) or (n_samples, n_tissues) float32 array
            Predicted tissue fractions.
        """
        if not self._has_model or self._trainer is None:
            # Fallback: return uniform distribution
            logger.debug(
                "No model loaded; returning uniform tissue fractions."
            )
            if methylation_data.ndim == 1:
                return np.ones(len(SUPPORTED_TISSUES)) / len(SUPPORTED_TISSUES)
            return np.ones((len(methylation_data), len(SUPPORTED_TISSUES))) / len(SUPPORTED_TISSUES)

        result = self._trainer.predict_batch(methylation_data)
        # If input was 1D, return 1D
        if methylation_data.ndim == 1:
            return result[0]
        return result

    def predict_sample(
        self,
        methylation_data: np.ndarray,
        return_top: int = 5,
    ) -> Dict[str, float]:
        """
        Predict tissue fractions for a single sample.

        Parameters
        ----------
        methylation_data : (n_cpg_features,) array
        return_top : int
            Number of top tissues to return.

        Returns
        -------
        fractions : dict
            {tissue_name: fraction} with metadata keys.
        """
        if not self._has_model or self._trainer is None:
            return {
                "_entropy": float(np.log(len(SUPPORTED_TISSUES))),
                "_n_active": len(SUPPORTED_TISSUES),
                "_top_tissue": SUPPORTED_TISSUES[0],
                "_top_fraction": 1.0 / len(SUPPORTED_TISSUES),
            }

        return self._trainer.predict(methylation_data)

    # ── Feature extraction for fusion ────────────────────────────

    def extract_all(
        self,
        sample: Dict[str, Any],
    ) -> Dict[str, float]:
        """
        Extract all tissue deconvolution features for one sample.

        The sample dict may contain:
        - 'methylation_data': dict with 'beta_values' or
          'tissue_cpgs' keys
        - 'tissue_fractions': pre-computed fractions (if available)
        - 'tissue_deconv_score': pre-computed score (if available)

        Parameters
        ----------
        sample : dict
            Sample data dictionary.

        Returns
        -------
        features : dict
            24+ scalar features ready for fusion.
        """
        # Case 1: Pre-computed tissue fractions available
        if "tissue_fractions" in sample:
            frac = sample["tissue_fractions"]
            if isinstance(frac, np.ndarray) and len(frac) == len(SUPPORTED_TISSUES):
                return self.feature_extractor.extract(frac)

        # Case 2: Methylation data available → run deconvolution
        methylation = self._extract_methylation_cpgs(sample)
        if methylation is not None:
            fractions = self.predict_tissue_fractions(methylation)
            return self.feature_extractor.extract(fractions)

        # Case 3: Pre-computed score only
        if "tissue_deconv_score" in sample:
            score = float(sample["tissue_deconv_score"])
            return {
                name: score if name == "tissue_max_cancer_fraction" else 0.0
                for name in self.feature_extractor.feature_names()
            }

        # Case 4: Nothing available → zeros
        logger.debug("No tissue data in sample; returning zero features.")
        return {name: 0.0 for name in self.feature_extractor.feature_names()}

    def extract_batch(
        self,
        samples: List[Dict[str, Any]],
    ) -> Dict[str, np.ndarray]:
        """
        Extract tissue deconvolution features for a batch.

        Parameters
        ----------
        samples : list of dict

        Returns
        -------
        features : dict
            {feature_name: (n_samples,) array}
        """
        features_list = [self.extract_all(s) for s in samples]
        all_names = self.feature_extractor.feature_names()

        return {
            name: np.array([f[name] for f in features_list], dtype=np.float32)
            for name in all_names
        }

    def _extract_methylation_cpgs(
        self,
        sample: Dict[str, Any],
    ) -> Optional[np.ndarray]:
        """
        Extract tissue-discriminative CpG beta values from sample.

        Returns (n_cpg_features,) array or None.
        """
        meth = sample.get("methylation_data")
        if meth is None:
            return None

        if isinstance(meth, dict):
            # Try 'tissue_cpgs' key first (specific to deconvolution)
            cpg_vals = meth.get("tissue_cpgs")
            if cpg_vals is not None and isinstance(cpg_vals, np.ndarray):
                if len(cpg_vals) == self.config.n_cpg_features:
                    return cpg_vals.astype(np.float32)
                # Truncate or pad
                if len(cpg_vals) < self.config.n_cpg_features:
                    padded = np.zeros(self.config.n_cpg_features, dtype=np.float32)
                    padded[:len(cpg_vals)] = cpg_vals
                    return padded
                return cpg_vals[:self.config.n_cpg_features].astype(np.float32)

            # Try 'beta_values' and subsample
            beta = meth.get("beta_values")
            if beta is not None and isinstance(beta, np.ndarray):
                flat = beta.flatten()
                if len(flat) >= self.config.n_cpg_features:
                    # Use first n_cpg_features (in practice, select
                    # tissue-discriminative CpGs via marker selection)
                    return flat[:self.config.n_cpg_features].astype(np.float32)
                padded = np.zeros(self.config.n_cpg_features, dtype=np.float32)
                padded[:len(flat)] = flat
                return padded

        if isinstance(meth, np.ndarray):
            flat = meth.flatten()
            if len(flat) >= self.config.n_cpg_features:
                return flat[:self.config.n_cpg_features].astype(np.float32)
            padded = np.zeros(self.config.n_cpg_features, dtype=np.float32)
            padded[:len(flat)] = flat
            return padded

        return None

    # ── Single modality score (for CrossAttentionFusion) ────────

    def to_modality(
        self,
        sample: Dict[str, Any],
    ) -> float:
        """
        Compute a single scalar score for the tissue deconvolution
        modality, compatible with CrossAttentionFusion.

        The score is:
        - max cancer-relevant tissue fraction (when model available)
        - 0.0 (when no model / fallback mode)

        This reflects the primary clinical signal: elevated fraction
        from a cancer-relevant tissue suggests tumor shedding.

        Parameters
        ----------
        sample : dict
            Sample data with methylation info.

        Returns
        -------
        score : float
            Scalar in [0, 1] for fusion.
        """
        if not self._has_model:
            return 0.0

        methylation = self._extract_methylation_cpgs(sample)
        if methylation is None:
            return 0.0

        fractions = self.predict_tissue_fractions(methylation)

        cancer_idx = [
            SUPPORTED_TISSUES.index(t)
            for t in CANCER_RELEVANT_TISSUES
            if t in SUPPORTED_TISSUES
        ]
        cancer_frac = float(np.sum(fractions[cancer_idx]))

        return cancer_frac

    def to_modality_batch(
        self,
        samples: List[Dict[str, Any]],
    ) -> np.ndarray:
        """
        Compute modality scores for a batch.

        Parameters
        ----------
        samples : list of dict

        Returns
        -------
        scores : (n_samples,) float32 array
        """
        return np.array(
            [self.to_modality(s) for s in samples],
            dtype=np.float32,
        )

    # ── Serialization ────────────────────────────────────────────

    def save_checkpoint(self, path: str) -> None:
        """
        Save the deconvolution model checkpoint.

        Parameters
        ----------
        path : str
            Output checkpoint path.
        """
        if self._trainer is not None:
            self._trainer.save_checkpoint(path)
        else:
            logger.warning("No trainer to save; skipping checkpoint.")

    # ── Properties ───────────────────────────────────────────────

    @property
    def has_model(self) -> bool:
        """Whether a trained deconvolution model is loaded."""
        return self._has_model

    @property
    def n_tissues(self) -> int:
        """Number of tissues in the atlas."""
        return self.config.n_tissues

    @property
    def n_features(self) -> int:
        """Number of scalar features per sample."""
        return self.feature_extractor.n_features
