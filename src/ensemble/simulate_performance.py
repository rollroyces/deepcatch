"""Real implementation for simulate_performance.py.

Original location: src/agent6-ensemble/simulate_performance.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent6-ensemble/` directory that was never committed to git.
"""

from __future__ import annotations

"""Monte-Carlo performance simulator for the deepcatch pipeline.

Used in MAML validation: simulate N synthetic cohorts with controlled
class-prior shift, then measure ensemble sensitivity at fixed specificity.
"""
import numpy as np


def simulate_cohort(n_pos: int = 200, n_neg: int = 800,
                    seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Synthetic two-modality cohort.

    Returns (X, y) where X is shape (n_pos + n_neg, 2) and y in {0, 1}.
    """
    rng = np.random.default_rng(seed)
    n = n_pos + n_neg
    X = rng.normal(0, 1, size=(n, 2))
    # Positives shift both modalities upward.
    X[:n_pos] += rng.normal(0.8, 0.4, size=(2,))
    y = np.zeros(n, dtype=int)
    y[:n_pos] = 1
    perm = rng.permutation(n)
    return X[perm], y[perm]


def measure_at_specificity(scores_pos: np.ndarray, scores_neg: np.ndarray,
                           target_specificity: float = 0.99) -> float:
    from .decision_threshold import threshold_at_specificity
    t = threshold_at_specificity(scores_neg, target_specificity)
    return float((scores_pos >= t).mean())


__all__ = ["simulate_cohort", "measure_at_specificity"]
