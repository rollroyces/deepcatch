#!/usr/bin/env python3
"""
Priming Integration — Priming Agent Adapter for DeepCatch Fusion
==================================================================

Connects the priming agents module to the DeepCatch multi-modal fusion
pipeline. Provides an adapter compatible with ``CrossAttentionFusion``
that adds priming-aware signal processing and prediction.

Integration Point
-----------------

In the DeepCatch CET pipeline, priming adds a **7th modality**:

.. code-block:: python

    from src.priming.integration import PrimingIntegration
    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion

    priming = PrimingIntegration()

    # Extract priming modality scores
    priming_scores = [priming.to_modality(s) for s in samples]

    # Now 7 modalities
    fusion = CrossAttentionFusion(n_modalities=7)
    fusion.fit(
        [frag, cnv, sero, mfr, gnn, deconv_scores, priming_scores],
        labels,
    )
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Union

import numpy as np

from .config import PrimingConfig, DEFAULT_CONFIG, AGENT_NAMES, AGENT_SPECS
from .pharmacokinetics import PKModel, OptimalDosingSchedule
from .signal_processing import PostPrimingDenoiser, SignalEnhancer, BaselineCorrector
from .response_predictor import PrimingResponsePredictor, PatientStratifier

logger = logging.getLogger(__name__)


class PrimingIntegration:
    """Adapter for integrating priming agents into DeepCatch fusion.

    Wraps PK/PD simulation, response prediction, signal processing,
    and patient stratification into a single interface compatible
    with the multi-modal fusion pipeline.

    Parameters
    ----------
    config : PrimingConfig, optional
        Configuration object.
    checkpoint : str, optional
        Path to predictor checkpoint.
    """

    def __init__(
        self,
        config: Optional[PrimingConfig] = None,
        checkpoint: Optional[str] = None,
    ):
        self.config = config or DEFAULT_CONFIG
        self.pk_model = PKModel()
        self.dosing = OptimalDosingSchedule(config=self.config)
        self.denoiser = PostPrimingDenoiser(window=self.config.denoising_window)
        self.enhancer = SignalEnhancer(boost_factor=self.config.detection_boost_factor)
        self.corrector = BaselineCorrector()

        try:
            self.predictor = PrimingResponsePredictor(config=self.config)
        except ImportError:
            logger.warning("PyTorch not available; response predictor disabled")
            self.predictor = None

        self.stratifier = PatientStratifier(config=self.config)
        if self.predictor is not None:
            self.stratifier.set_predictor(self.predictor)

        if checkpoint is not None:
            self.load_checkpoint(checkpoint)

    def process_signal(self, sample: dict) -> dict:
        """Process a single sample through the full priming pipeline.

        Parameters
        ----------
        sample : dict with keys:
            features : np.ndarray — raw cfDNA features
            patient_data : dict — patient profile
            priming_timing : float, optional — when priming was given (samples)
            pre_priming_features : np.ndarray, optional — pre-priming baseline
            agent : str, optional — administered agent type

        Returns
        -------
        dict with:
            priming_score : float (0-1, benefit score)
            boost_factor : float (predicted ctDNA increase)
            signal_quality : float (SNR after processing)
            denoised_features : ndarray
            enhanced_features : ndarray
            optimal_agent : str
            stratification : str
            recommended_dose_mg : float
            recommended_timing_h : float
            toxicity_risk : float
        """
        features = sample.get("features")
        patient_data = sample.get("patient_data", {})
        priming_timing = sample.get("priming_timing", None)
        pre_priming = sample.get("pre_priming_features", None)
        agent = sample.get("agent", None)

        if features is None:
            logger.warning("No features provided; returning zero-result")
            return self._empty_result()

        features = np.asarray(features, dtype=np.float64)

        # If no patient data, can't compute PK → minimal processing
        if not patient_data:
            denoised = self.denoiser.denoise(features, priming_timing)
            return {
                "priming_score": 0.0,
                "boost_factor": 1.0,
                "signal_quality": 1.0,
                "denoised_features": denoised,
                "enhanced_features": denoised,
                "optimal_agent": "none",
                "stratification": "poor_candidate",
                "recommended_dose_mg": 0.0,
                "recommended_timing_h": 0.0,
                "toxicity_risk": 1.0,
            }

        # Find best agent if not specified
        if agent is None:
            agent = self._find_best_agent(patient_data)

        # Compute optimal dosing
        dosing_result = self.dosing.compute(agent, patient_data)

        # Denoise signal
        denoised = self.denoiser.denoise(features, priming_timing)

        # Baseline correction if pre-priming available
        if pre_priming is not None:
            denoised = self.corrector.correct(denoised, pre_priming)

        # Estimate signal quality
        noise_level = self.denoiser.estimate_background_noise(denoised)
        signal_power = np.var(denoised) + 1e-10
        signal_quality = signal_power / (noise_level + 1e-10)

        # Enhance signal
        predicted_boost = dosing_result["predicted_ctDNA_boost"]
        enhanced = self.enhancer.enhance(
            features=denoised,
            signal_to_noise=signal_quality,
            priming_boost=predicted_boost,
        )

        # Stratify patient
        stratification = self.stratifier.stratify(
            patient_data=patient_data,
            agent=agent,
            predicted_boost=predicted_boost,
            predicted_toxicity=dosing_result["predicted_toxicity_risk"],
        )

        # Priming score: 0 = no benefit, 1 = maximal benefit
        strat_scores = {
            "ideal_candidate": 0.9,
            "moderate_candidate": 0.5,
            "poor_candidate": 0.1,
        }
        base_score = strat_scores.get(stratification, 0.1)
        boost_normalized = min(predicted_boost / self.config.detection_boost_factor, 1.0)
        priming_score = float(base_score * 0.5 + boost_normalized * 0.5)
        priming_score = float(np.clip(priming_score, 0.0, 1.0))

        return {
            "priming_score": priming_score,
            "boost_factor": float(predicted_boost),
            "signal_quality": float(signal_quality),
            "denoised_features": denoised,
            "enhanced_features": enhanced,
            "optimal_agent": agent,
            "stratification": stratification,
            "recommended_dose_mg": dosing_result["optimal_dose_mg"],
            "recommended_timing_h": dosing_result["optimal_timing_h"],
            "toxicity_risk": dosing_result["predicted_toxicity_risk"],
        }

    def to_modality(self, sample: dict) -> float:
        """Produce a single scalar modality score for fusion.

        Parameters
        ----------
        sample : dict
            Sample data (see process_signal).

        Returns
        -------
        float : Estimated ctDNA boost from optimal priming.
            0 = no benefit, 1 = 10x+ boost.
        """
        result = self.process_signal(sample)
        return result["priming_score"]

    def extract_all(self, sample: dict) -> Dict[str, float]:
        """Extract all priming-related features for fusion.

        Parameters
        ----------
        sample : dict
            Sample data.

        Returns
        -------
        dict[str, float] : ~15 features for multimodal fusion.
        """
        result = self.process_signal(sample)
        patient_data = sample.get("patient_data", {})

        features = {
            "priming_score": result["priming_score"],
            "priming_boost_factor": result["boost_factor"],
            "priming_signal_quality": result["signal_quality"],
            "priming_toxicity_risk": result["toxicity_risk"],
            "priming_recommended_dose": float(result["recommended_dose_mg"]),
            "priming_recommended_timing": float(result["recommended_timing_h"]),
            "priming_liver_function": float(patient_data.get("liver_function", 1.0)),
            "priming_renal_function": float(patient_data.get("renal_function", 1.0)),
            "priming_tumor_fraction": float(patient_data.get("tumor_fraction", 0.0)),
            "priming_baseline_cfdna": float(patient_data.get("baseline_cfdna_ng_ml", 20.0)),
            "priming_agent_half_life": float(
                AGENT_SPECS.get(result.get("optimal_agent", "liposome"), {}).get("half_life_hours", 18.0)
            ),
            "priming_ideal_candidate": 1.0 if result.get("stratification") == "ideal_candidate" else 0.0,
            "priming_moderate_candidate": 1.0 if result.get("stratification") == "moderate_candidate" else 0.0,
            "priming_poor_candidate": 1.0 if result.get("stratification") == "poor_candidate" else 0.0,
            "priming_denoising_quality": float(
                self.denoiser.estimate_background_noise(
                    result.get("denoised_features", np.zeros(1))
                )
            ),
        }

        return features

    def process_batch(self, samples: List[dict]) -> List[dict]:
        """Process a batch of samples.

        Parameters
        ----------
        samples : list[dict]
            List of sample dictionaries.

        Returns
        -------
        list[dict] : Processed results.
        """
        return [self.process_signal(s) for s in samples]

    def _find_best_agent(self, patient_data: dict) -> str:
        """Find the best priming agent for a patient.

        Uses PK model + stratifier to rank agents.

        Parameters
        ----------
        patient_data : dict
            Patient profile.

        Returns
        -------
        str : Best agent name.
        """
        best_agent = "liposome"  # default: longest half-life, good safety
        best_score = -np.inf

        for agent in AGENT_NAMES:
            try:
                dosing = self.dosing.compute(agent, patient_data)
                boost = dosing["predicted_ctDNA_boost"]
                tox = dosing["predicted_toxicity_risk"]
                score = boost * (1.0 - 0.5 * tox)
                if score > best_score:
                    best_score = score
                    best_agent = agent
            except Exception as e:
                logger.debug(f"Agent {agent} failed: {e}")

        return best_agent

    def load_checkpoint(self, path: str) -> None:
        """Load predictor checkpoint.

        Parameters
        ----------
        path : str
            Path to checkpoint file.
        """
        try:
            import torch

            if self.predictor is not None:
                state = torch.load(path, map_location=self.config.device, weights_only=True)
                self.predictor.load_state_dict(state)
                self.predictor.eval()
                logger.info(f"Loaded predictor checkpoint from {path}")
        except ImportError:
            logger.warning("PyTorch not available; cannot load checkpoint")
        except FileNotFoundError:
            logger.error(f"Checkpoint not found: {path}")
        except Exception as e:
            logger.error(f"Failed to load checkpoint: {e}")

    def save_checkpoint(self, path: str) -> None:
        """Save predictor checkpoint.

        Parameters
        ----------
        path : str
            Output path.
        """
        if self.predictor is not None:
            try:
                import torch

                torch.save(self.predictor.state_dict(), path)
                logger.info(f"Saved predictor checkpoint to {path}")
            except Exception as e:
                logger.error(f"Failed to save checkpoint: {e}")

    def _empty_result(self) -> dict:
        """Return an empty/zero result for missing data."""
        zero_features = np.zeros(1)
        return {
            "priming_score": 0.0,
            "boost_factor": 1.0,
            "signal_quality": 0.0,
            "denoised_features": zero_features,
            "enhanced_features": zero_features,
            "optimal_agent": "none",
            "stratification": "poor_candidate",
            "recommended_dose_mg": 0.0,
            "recommended_timing_h": 0.0,
            "toxicity_risk": 1.0,
        }
