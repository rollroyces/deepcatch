"""Real implementation for compare_methods.py.

Original location: src/agent3-longitudinal/compare_methods.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Compare two longitudinal detectors on a synthetic cohort."""
from dataclasses import dataclass

import numpy as np

from .simulation import simulate_trajectory
from .bayesian_changepoint import BayesianChangepointDetector


@dataclass
class ComparisonResult:
    method_a_delay: float
    method_b_delay: float
    n: int


def compare(n_patients: int = 50, seed: int = 0) -> ComparisonResult:
    rng = np.random.default_rng(seed)
    bcd = BayesianChangepointDetector()
    delays_a, delays_b = [], []
    for i in range(n_patients):
        traj = simulate_trajectory(seed=seed + i)
        cp, _ = bcd.detect(traj)
        # Method A = BCD delay from true changepoint (4)
        delays_a.append(max(0, cp - 4))
        # Method B = simple threshold-crossing at mean + 2σ
        thresh = traj[:4].mean() + 2 * traj[:4].std()
        idx_b = next((i for i, v in enumerate(traj) if v > thresh), len(traj))
        delays_b.append(max(0, idx_b - 4))
    return ComparisonResult(
        method_a_delay=float(np.mean(delays_a)),
        method_b_delay=float(np.mean(delays_b)),
        n=n_patients,
    )


__all__ = ["ComparisonResult", "compare"]
