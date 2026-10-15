"""Real implementation for decision_threshold.py.

Original location: src/agent6-ensemble/decision_threshold.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent6-ensemble/` directory that was never committed to git.
"""

from __future__ import annotations

"""Operating-point selection for binary classifiers.

The deepcatch pipeline needs a single threshold per cancer type that
satisfies a fixed-specificity gate (no threshold optimization on the eval
cohort — see ULTRA_EARLY_READINESS.md and PIPELINE_AUDIT.md for the
fixed-specificity contract).
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class OperatingPoint:
    threshold: float
    specificity: float
    sensitivity: float


def threshold_at_specificity(scores_neg: np.ndarray,
                             target_specificity: float = 0.99) -> float:
    """Pick the cutoff ``t`` such that ~``target_specificity`` of negative
    cohort scores fall below ``t``.
    """
    if not 0.0 < target_specificity < 1.0:
        raise ValueError("target_specificity must be in (0, 1)")
    return float(np.quantile(scores_neg, target_specificity))


def sweep_specificities(scores_pos: np.ndarray, scores_neg: np.ndarray,
                        specificities: np.ndarray | None = None
                        ) -> list[OperatingPoint]:
    """Return per-specificity operating points for plotting."""
    if specificities is None:
        specificities = np.linspace(0.80, 0.999, 30)
    out = []
    for sp in specificities:
        t = threshold_at_specificity(scores_neg, float(sp))
        sens = float((scores_pos >= t).mean())
        out.append(OperatingPoint(threshold=t, specificity=float(sp), sensitivity=sens))
    return out


__all__ = ["OperatingPoint", "threshold_at_specificity", "sweep_specificities"]
