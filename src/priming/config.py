#!/usr/bin/env python3
"""
Priming Agents Configuration
=============================

Dataclass-based hyperparameter management for the priming agents module.
Models the interaction between cfDNA priming agents (Amplifyer Bio, Martin-Alonso
et al. 2024, Science) and DeepCatch's multi-modal detection pipeline.

All config values have defaults grounded in published literature:
- Martin-Alonso et al. (2024) Science: priming agents increase ctDNA >10x
- Wagner et al. (2023) Nat Biotechnol: ctDNA half-life ~30-120 min
- Thierry et al. (2016) Cancer Metastasis Rev: cfDNA clearance dynamics
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ── Agent PK Parameters (Literature-based) ───────────────────────
# Half-lives and properties from published PK studies for each priming agent type

AGENT_SPECS: Dict[str, dict] = {
    "scFv": {
        "description": "Single-chain variable fragment — anti-cfDNA antibody fragment",
        "half_life_hours": 2.5,
        "molecular_weight_kDa": 27.0,
        "volume_of_distribution_L": 15.0,
        "clearance_rate_L_h": 4.2,
        "bioavailability": 0.65,
        "protein_binding": 0.12,
        "reference": "Bates et al. (2019) — scFv PK in primate models",
    },
    "liposome": {
        "description": "Liposomal nanoparticle — PEGylated liposome for cfDNA scavenging",
        "half_life_hours": 18.0,
        "molecular_weight_kDa": 100_000.0,
        "volume_of_distribution_L": 3.8,
        "clearance_rate_L_h": 0.15,
        "bioavailability": 0.90,
        "protein_binding": 0.85,
        "reference": "Gabizon et al. (1994) — Stealth liposome PK",
    },
    "nanoparticle": {
        "description": "Polymeric nanoparticle — PLGA-based cfDNA binding particle",
        "half_life_hours": 6.0,
        "molecular_weight_kDa": 40_000.0,
        "volume_of_distribution_L": 8.2,
        "clearance_rate_L_h": 0.95,
        "bioavailability": 0.75,
        "protein_binding": 0.40,
        "reference": "Alexis et al. (2008) — PLGA nanoparticle PK",
    },
    "polymeric_micelle": {
        "description": "Polymeric micelle — amphiphilic block copolymer for cfDNA capture",
        "half_life_hours": 8.0,
        "molecular_weight_kDa": 15_000.0,
        "volume_of_distribution_L": 10.0,
        "clearance_rate_L_h": 0.87,
        "bioavailability": 0.70,
        "protein_binding": 0.95,
        "reference": "Yokoyama et al. (2014) — Polymeric micelle PK",
    },
    "dendrimer": {
        "description": "PAMAM dendrimer — highly branched cfDNA-binding macromolecule",
        "half_life_hours": 4.0,
        "molecular_weight_kDa": 14.0,
        "volume_of_distribution_L": 12.0,
        "clearance_rate_L_h": 2.1,
        "bioavailability": 0.55,
        "protein_binding": 0.30,
        "reference": "Malik et al. (2000) — Dendrimer PK in rodent models",
    },
}

AGENT_NAMES: Tuple[str, ...] = ("scFv", "liposome", "nanoparticle", "polymeric_micelle", "dendrimer")


# ── Patient Feature Specification ─────────────────────────────────

PATIENT_FEATURE_SPEC: Dict[str, Dict[str, str]] = {
    "age": {
        "desc": "Patient age in years",
        "range": "[18, 100]",
        "feature_idx": 0,
    },
    "weight_kg": {
        "desc": "Patient weight in kg",
        "range": "[30, 200]",
        "feature_idx": 1,
    },
    "bmi": {
        "desc": "Body mass index",
        "range": "[15, 60]",
        "feature_idx": 2,
    },
    "liver_function": {
        "desc": "Liver function score (0=impaired, 1=normal)",
        "range": "[0, 1]",
        "feature_idx": 3,
    },
    "renal_function": {
        "desc": "Renal function score (0=impaired, 1=normal)",
        "range": "[0, 1]",
        "feature_idx": 4,
    },
    "tumor_type": {
        "desc": "Encoded tumor type (0=lung, 1=colorectal, 2=breast, 3=pancreatic, 4=other)",
        "range": "{0,1,2,3,4}",
        "feature_idx": 5,
    },
    "tumor_stage": {
        "desc": "Cancer stage (0=none, 1=I, 2=II, 3=III, 4=IV)",
        "range": "{0,1,2,3,4}",
        "feature_idx": 6,
    },
    "tumor_fraction": {
        "desc": "Estimated tumor fraction in cfDNA",
        "range": "[0, 1]",
        "feature_idx": 7,
    },
    "baseline_cfdna_ng_ml": {
        "desc": "Baseline cfDNA concentration (ng/mL)",
        "range": "[1, 500]",
        "feature_idx": 8,
    },
    "albumin_g_L": {
        "desc": "Serum albumin (g/L)",
        "range": "[20, 55]",
        "feature_idx": 9,
    },
    "bilirubin_umol_L": {
        "desc": "Total bilirubin (umol/L)",
        "range": "[2, 50]",
        "feature_idx": 10,
    },
    "creatinine_umol_L": {
        "desc": "Serum creatinine (umol/L)",
        "range": "[30, 500]",
        "feature_idx": 11,
    },
    "alt_IU_L": {
        "desc": "Alanine aminotransferase (IU/L)",
        "range": "[5, 200]",
        "feature_idx": 12,
    },
    "ast_IU_L": {
        "desc": "Aspartate aminotransferase (IU/L)",
        "range": ["5", "200"],
        "feature_idx": 13,
    },
    "platelet_count": {
        "desc": "Platelet count (10^9/L)",
        "range": "[50, 600]",
        "feature_idx": 14,
    },
    "neutrophil_count": {
        "desc": "Neutrophil count (10^9/L)",
        "range": "[1, 20]",
        "feature_idx": 15,
    },
    "hemoglobin_g_L": {
        "desc": "Hemoglobin (g/L)",
        "range": "[70, 180]",
        "feature_idx": 16,
    },
    "crp_mg_L": {
        "desc": "C-reactive protein (mg/L)",
        "range": "[0, 200]",
        "feature_idx": 17,
    },
    "prior_treatment": {
        "desc": "Prior systemic treatment (0=none, 1=chemo, 2=immuno, 3=targeted)",
        "range": "{0,1,2,3}",
        "feature_idx": 18,
    },
    "performance_status": {
        "desc": "ECOG performance status (0-5)",
        "range": "{0,1,2,3,4,5}",
        "feature_idx": 19,
    },
}

N_PATIENT_FEATURES: int = len(PATIENT_FEATURE_SPEC)  # 20


@dataclass
class PrimingConfig:
    """Configuration for Priming Agents AI Module.

    Attributes
    ----------
    n_agents : int
        Number of priming agent types (scFv, liposome, nanoparticle,
        polymeric micelle, dendrimer).
    n_patient_features : int
        Number of patient-level features for response prediction.
    latent_dim : int
        Dimension of latent PK/PD representation.
    denoising_window : int
        Number of sequential samples for denoising moving average.
    detection_boost_factor : float
        Maximum expected ctDNA concentration increase from priming (>10x per
        Martin-Alonso et al. 2024).
    mlp_hidden : list[int]
        Hidden layer sizes for response predictor MLP.
    predictor_output : int
        Number of output dimensions (boost_factor, time_to_peak, toxicity_risk).
    stratification_thresholds : dict
        Cutoffs for patient stratification.
    learning_rate : float
        Learning rate for optimizer.
    device : str
        Torch device string.
    """

    n_agents: int = 5
    n_patient_features: int = N_PATIENT_FEATURES
    latent_dim: int = 64
    denoising_window: int = 3
    detection_boost_factor: float = 10.0
    mlp_hidden: list = field(default_factory=lambda: [64, 32])
    predictor_output: int = 3
    stratification_thresholds: Dict[str, float] = field(default_factory=lambda: {
        "ideal_boost_min": 5.0,
        "ideal_risk_max": 0.3,
        "poor_boost_max": 1.5,
        "poor_risk_min": 0.6,
    })
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    device: str = "cpu"
    seed: int = 42

    @property
    def n_agent_params(self) -> int:
        """Number of PK parameters per agent (half_life, Vd, CL, F, PB)."""
        return 5

    @property
    def input_dim(self) -> int:
        """Input dimension = patient features + agent PK params."""
        return self.n_patient_features + self.n_agent_params

    def to_dict(self) -> dict:
        """Serialize config to dictionary."""
        return {
            "n_agents": self.n_agents,
            "n_patient_features": self.n_patient_features,
            "latent_dim": self.latent_dim,
            "denoising_window": self.denoising_window,
            "detection_boost_factor": self.detection_boost_factor,
            "mlp_hidden": self.mlp_hidden,
            "predictor_output": self.predictor_output,
            "stratification_thresholds": self.stratification_thresholds,
            "learning_rate": self.learning_rate,
            "weight_decay": self.weight_decay,
            "device": self.device,
            "seed": self.seed,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "PrimingConfig":
        """Deserialize config from dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ── Pre-built configurations ─────────────────────────────────────

# Default: balanced for research use
DEFAULT_CONFIG = PrimingConfig()

# Prototype: minimal config for CI / rapid testing
PROTOTYPE_CONFIG = PrimingConfig(
    latent_dim=16,
    mlp_hidden=[16, 8],
    denoising_window=2,
    device="cpu",
)

# Production: larger model with tuned parameters
PRODUCTION_CONFIG = PrimingConfig(
    latent_dim=128,
    mlp_hidden=[128, 64],
    denoising_window=5,
    detection_boost_factor=12.0,
    learning_rate=5e-4,
    weight_decay=1e-4,
    stratification_thresholds={
        "ideal_boost_min": 6.0,
        "ideal_risk_max": 0.25,
        "poor_boost_max": 1.2,
        "poor_risk_min": 0.65,
    },
)
