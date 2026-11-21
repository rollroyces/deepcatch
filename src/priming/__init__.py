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
from .signal_processing import PostPrimingDenoiser, SignalEnhancer, BaselineCorrector
from .data import (
    generate_patient_profiles,
    apply_priming_effect,
    simulate_clinical_trial,
)
from .whitepaper import generate_whitepaper_sections
from .pitch_deck import generate_pitch_deck

# `torch` is an optional dependency: it lives in the `api` extra, not in
# the runtime dependencies. Importing the neural parts at module level
# made `import src.priming` fail in a venv without it -- and a package
# that cannot be imported cannot be tested, so the collection of
# src/priming/test_integration.py was an error rather than a skip, which
# interrupted the whole run of a CPU-only install. The parts that do not
# need torch are always available; the neural ones need the extra.
try:
    from .response_predictor import (
        PrimingResponsePredictor,
        PatientStratifier,
    )
    from .integration import PrimingIntegration
except ImportError as exc:  # pragma: no cover - needs a torch-free venv
    _OPTIONAL_TORCH_IMPORT_ERROR = exc

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

# Without torch the neural names are not bound, so `from src.priming
# import *` would raise. The list is filtered to what this install can
# actually provide.
__all__ = [name for name in __all__ if name in globals()]
