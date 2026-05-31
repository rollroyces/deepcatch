"""
DeepCatch Priming Agents — AI-Guided Priming for Enhanced Liquid Biopsy
=========================================================================

Integrates priming agent PK/PD modeling (Amplifyer Bio, Martin-Alonso et al.
2024, Science) with DeepCatch's multi-modal foundation model for improved
ctDNA detection sensitivity.

References
----------
.. [1] Martin-Alonso, C. et al. (2024) Science 383(6678):eadf2341.
       Priming agents transiently reduce cfDNA clearance.
.. [2] Cohen, J.D. et al. (2018) Science 359(6378):926-930.
       CancerSEEK multi-analyte blood test.
"""

from __future__ import annotations

from .config import (
    PrimingConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    PRODUCTION_CONFIG,
    AGENT_SPECS,
    AGENT_NAMES,
    PATIENT_FEATURE_SPEC,
    N_PATIENT_FEATURES,
)
from .pharmacokinetics import PKModel, OptimalDosingSchedule
from .response_predictor import PrimingResponsePredictor, PatientStratifier
from .signal_processing import PostPrimingDenoiser, SignalEnhancer, BaselineCorrector
from .data import (
    generate_patient_profiles,
    apply_priming_effect,
    simulate_clinical_trial,
)
from .integration import PrimingIntegration
from .whitepaper import generate_whitepaper_sections
from .pitch_deck import generate_pitch_deck

__all__ = [
    # Config
    "PrimingConfig",
    "DEFAULT_CONFIG",
    "PROTOTYPE_CONFIG",
    "PRODUCTION_CONFIG",
    "AGENT_SPECS",
    "AGENT_NAMES",
    "PATIENT_FEATURE_SPEC",
    "N_PATIENT_FEATURES",
    # PK/PD
    "PKModel",
    "OptimalDosingSchedule",
    # ML
    "PrimingResponsePredictor",
    "PatientStratifier",
    # Signal processing
    "PostPrimingDenoiser",
    "SignalEnhancer",
    "BaselineCorrector",
    # Data
    "generate_patient_profiles",
    "apply_priming_effect",
    "simulate_clinical_trial",
    # Integration
    "PrimingIntegration",
    # Documentation
    "generate_whitepaper_sections",
    "generate_pitch_deck",
]
