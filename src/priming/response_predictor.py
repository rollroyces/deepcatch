#!/usr/bin/env python3
"""
Priming Response Predictor
============================

Lightweight neural network (~5K params) that predicts patient-specific
response to priming agents: ctDNA boost factor, time to peak, and toxicity risk.

Also provides rule-based + learned patient stratification to identify
ideal candidates for priming-enhanced liquid biopsy.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None
    nn = None
    F = None

from .config import (
    PrimingConfig,
    DEFAULT_CONFIG,
    AGENT_SPECS,
    AGENT_NAMES,
    N_PATIENT_FEATURES,
)

logger = logging.getLogger(__name__)


# ── Torch Model ──────────────────────────────────────────────────


class PrimingResponsePredictor(nn.Module):
    """MLP-based predictor for patient-specific priming response.

    Input: [patient_features (20) + agent_PK_params (5)] = 25 features
    Architecture: 25 → 64 → 32 → 3 (boost_factor, time_to_peak, toxicity_risk)
    ~5K trainable parameters.

    Parameters
    ----------
    config : PrimingConfig
        Configuration object.
    """

    def __init__(self, config: Optional[PrimingConfig] = None):
        if not _HAS_TORCH:
            raise ImportError("PyTorch is required for PrimingResponsePredictor")

        super().__init__()
        self.config = config or DEFAULT_CONFIG

        input_dim = self.config.n_patient_features + self.config.n_agent_params
        hidden = self.config.mlp_hidden
        output_dim = self.config.predictor_output

        layers = []
        prev_dim = input_dim
        for h in hidden:
            layers.append(nn.Linear(prev_dim, h))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.1))
            prev_dim = h
        layers.append(nn.Linear(prev_dim, output_dim))

        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass.

        Parameters
        ----------
        x : torch.Tensor, shape (batch, 25)
            Concatenated patient features + agent PK params.

        Returns
        -------
        torch.Tensor, shape (batch, 3)
            [boost_factor, time_to_peak, toxicity_risk]
        """
        out = self.net(x)
        # Apply output transformations:
        # boost_factor: softplus for positive values (add small epsilon to avoid exactly 0)
        # time_to_peak: softplus for positive hours
        # toxicity_risk: sigmoid for [0, 1]
        boost = F.softplus(out[:, 0:1]) + 1.0  # minimum boost ~2.0x (softplus(0)+1≈1.69)
        time_to_peak = F.softplus(out[:, 1:2]) + 1.0  # minimum time ~1.7h
        toxicity = torch.sigmoid(out[:, 2:3])
        return torch.cat([boost, time_to_peak, toxicity], dim=1)

    def predict(self, patient_data: dict, agent: str) -> dict:
        """Predict priming response for a single patient-agent pair.

        Parameters
        ----------
        patient_data : dict
            Patient feature dictionary. Keys should match PATIENT_FEATURE_SPEC.
        agent : str
            Agent type name.

        Returns
        -------
        dict with:
            predicted_boost_factor : float
            predicted_time_to_peak_h : float
            predicted_toxicity_risk : float
            agent : str
        """
        self.eval()

        patient_vec = self._encode_patient(patient_data)
        agent_vec = self._encode_agent(agent)
        x = np.concatenate([patient_vec, agent_vec]).astype(np.float32)
        x_tensor = torch.from_numpy(x).unsqueeze(0)

        with torch.no_grad():
            output = self.forward(x_tensor).squeeze(0).numpy()

        return {
            "predicted_boost_factor": round(float(output[0]), 3),
            "predicted_time_to_peak_h": round(float(output[1]), 3),
            "predicted_toxicity_risk": round(float(output[2]), 3),
            "agent": agent,
        }

    def _encode_patient(self, patient_data: dict) -> np.ndarray:
        """Encode patient dict → fixed-length feature vector."""
        vec = np.zeros(self.config.n_patient_features, dtype=np.float32)
        feature_map = {
            "age": 0,
            "weight_kg": 1,
            "bmi": 2,
            "liver_function": 3,
            "renal_function": 4,
            "tumor_type": 5,
            "tumor_stage": 6,
            "tumor_fraction": 7,
            "baseline_cfdna_ng_ml": 8,
            "albumin_g_L": 9,
            "bilirubin_umol_L": 10,
            "creatinine_umol_L": 11,
            "alt_IU_L": 12,
            "ast_IU_L": 13,
            "platelet_count": 14,
            "neutrophil_count": 15,
            "hemoglobin_g_L": 16,
            "crp_mg_L": 17,
            "prior_treatment": 18,
            "performance_status": 19,
        }
        for key, idx in feature_map.items():
            if key in patient_data and patient_data[key] is not None:
                vec[idx] = float(patient_data[key])
        return vec

    def _encode_agent(self, agent: str) -> np.ndarray:
        """Encode agent → PK parameter vector."""
        if agent not in AGENT_NAMES:
            raise ValueError(f"Unknown agent: {agent}. Choose from {AGENT_NAMES}")

        spec = AGENT_SPECS[agent]
        # Normalize each PK param to [0, 1] range
        vec = np.array([
            spec["half_life_hours"] / 24.0,          # normalize to ~24h
            spec["volume_of_distribution_L"] / 20.0,  # normalize to ~20L
            spec["clearance_rate_L_h"] / 5.0,         # normalize to ~5 L/h
            spec["bioavailability"],                  # already [0,1]
            spec["protein_binding"],                  # already [0,1]
        ], dtype=np.float32)
        return vec

    def estimated_params(self) -> int:
        """Return estimated number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ── Patient Stratifier ───────────────────────────────────────────


class PatientStratifier:
    """Rule-based + learned patient stratification for priming candidacy.

    Categories:
    - "ideal_candidate": High expected benefit, low risk
    - "moderate_candidate": Moderate benefit, acceptable risk
    - "poor_candidate": Low benefit or high risk

    Decision logic combines:
    - Tumor type/stage eligibility
    - Liver/renal function
    - Predicted ctDNA boost from response predictor
    - Toxicity risk
    """

    def __init__(self, config: Optional[PrimingConfig] = None):
        self.config = config or DEFAULT_CONFIG
        self.predictor: Optional[PrimingResponsePredictor] = None

    def set_predictor(self, predictor: PrimingResponsePredictor) -> None:
        """Attach a trained response predictor for learned thresholds."""
        self.predictor = predictor

    def stratify(
        self,
        patient_data: dict,
        agent: Optional[str] = None,
        predicted_boost: Optional[float] = None,
        predicted_toxicity: Optional[float] = None,
    ) -> str:
        """Stratify a patient into a candidacy category.

        Parameters
        ----------
        patient_data : dict
            Patient features.
        agent : str, optional
            Agent type. If provided and predictor is attached, predictions
            will be computed automatically.
        predicted_boost : float, optional
            Pre-computed boost prediction.
        predicted_toxicity : float, optional
            Pre-computed toxicity risk.

        Returns
        -------
        str : "ideal_candidate", "moderate_candidate", or "poor_candidate"
        """
        thresh = self.config.stratification_thresholds

        # Compute predictions if not provided
        if predicted_boost is None and self.predictor is not None and agent:
            preds = self.predictor.predict(patient_data, agent)
            predicted_boost = preds["predicted_boost_factor"]
            predicted_toxicity = preds["predicted_toxicity_risk"]

        # Rule-based checks
        liver_fn = patient_data.get("liver_function", 1.0)
        renal_fn = patient_data.get("renal_function", 1.0)
        tumor_type = patient_data.get("tumor_type", 0)
        tumor_stage = patient_data.get("tumor_stage", 1)
        performance = patient_data.get("performance_status", 0)
        prior_tx = patient_data.get("prior_treatment", 0)

        # Contraindications
        severe_organ = (liver_fn < 0.3) or (renal_fn < 0.3)
        poor_performance = performance >= 3
        contraindicated = severe_organ or poor_performance

        if contraindicated:
            return "poor_candidate"

        # Use predictions if available
        if predicted_boost is not None and predicted_toxicity is not None:
            if predicted_boost >= thresh["ideal_boost_min"] and predicted_toxicity <= thresh["ideal_risk_max"]:
                return "ideal_candidate"
            elif predicted_boost <= thresh["poor_boost_max"] or predicted_toxicity >= thresh["poor_risk_min"]:
                return "poor_candidate"
            else:
                return "moderate_candidate"

        # Fallback: rule-based stratification without ML predictions
        score = 0.0

        # Tumor type scoring (solid tumors benefit more)
        tumor_type_scores = {0: 0.8, 1: 0.7, 2: 0.6, 3: 0.5, 4: 0.4}  # lung, colorectal, breast, pancreatic, other
        score += tumor_type_scores.get(int(tumor_type), 0.5) * 0.3

        # Stage scoring (earlier stage = higher need)
        stage_scores = {1: 0.9, 2: 0.7, 3: 0.5, 4: 0.3, 0: 0.1}
        score += stage_scores.get(int(tumor_stage), 0.5) * 0.2

        # Organ function
        score += liver_fn * 0.25
        score += renal_fn * 0.15

        # Performance status
        score += (1.0 - performance / 5.0) * 0.1

        # Thresholds
        if score > 0.65:
            return "ideal_candidate"
        elif score > 0.35:
            return "moderate_candidate"
        else:
            return "poor_candidate"

    def batch_stratify(
        self, patients: List[dict], agent: Optional[str] = "liposome"
    ) -> List[str]:
        """Stratify a batch of patients.

        Parameters
        ----------
        patients : list[dict]
            List of patient data dictionaries.
        agent : str, optional
            Agent to use for predictions.

        Returns
        -------
        list[str] : Stratification categories.
        """
        results = []
        for patient in patients:
            try:
                category = self.stratify(patient, agent=agent)
            except Exception as e:
                logger.warning(f"Stratification failed: {e}")
                category = "poor_candidate"
            results.append(category)
        return results
