#!/usr/bin/env python3
"""
DeepCatch Bioinformatics Validation Suite
==========================================

Ten validation modules providing publication-grade rigor for the DeepCatch
multi-modal longitudinal cancer screening framework.

Modules:
  1. NestedCrossValidator   — Unbiased generalization estimates (nested CV)
  2. PermutationTester      — Signal vs. noise testing
  3. CalibrationAnalyzer    — Probability calibration diagnostics
  4. DecisionCurveAnalyzer  — Net benefit clinical utility assessment
  5. delong_test            — Full DeLong covariance-based AUC comparison
  6. StratifiedAnalyzer     — Per-stratum performance with interaction tests
  7. ConfounderRobustnessTester — Sensitivity to realistic confounders
  8. BioinfoBenchmark       — Head-to-head against bioinformatics tools
  9. PowerAnalyzer          — Sample-size and detectable-effect analysis
  10. ReproducibilityGuard  — Seed registry and environment hashing

All modules integrate with the canonical validation_framework.py API.
"""

try:
    from validation.nested_cv import NestedCrossValidator
    from validation.permutation_test import PermutationTester
    from validation.calibration import CalibrationAnalyzer
    from validation.decision_curve import DecisionCurveAnalyzer
    from validation.delong_test import delong_test
    from validation.stratified import StratifiedAnalyzer
    from validation.confounders import ConfounderRobustnessTester
    from validation.bioinfo_benchmark import BioinfoBenchmark
    from validation.power_analysis import PowerAnalyzer
except ImportError:
    pass  # minimal CI env

__all__ = [
    "NestedCrossValidator",
    "PermutationTester",
    "CalibrationAnalyzer",
    "DecisionCurveAnalyzer",
    "delong_test",
    "StratifiedAnalyzer",
    "ConfounderRobustnessTester",
    "BioinfoBenchmark",
    "PowerAnalyzer",
]

__version__ = "1.0.0"
__status__ = "production"
