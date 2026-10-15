"""Real implementation for run_final_fixed.py.

Original location: src/agent3-longitudinal/run_final_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version of :mod:`run_final` — uses ``FixedImprovedDetector``."""
from .run_final import run_final
from .improved_methods_fixed import FixedImprovedDetector
from .bayesian_changepoint import BayesianChangepointDetector
from .simulation import simulate_trajectory
import numpy as np


def run_final_fixed(n_patients: int = 50, seed: int = 0) -> dict:
    bcd = BayesianChangepointDetector()
    imp = FixedImprovedDetector()
    plain_delays, improved_delays = [], []
    for i in range(n_patients):
        traj = simulate_trajectory(seed=seed + i)
        cp_plain, _ = bcd.detect(traj)
        cp_imp, _ = imp.detect(traj)
        plain_delays.append(max(0, cp_plain - 4))
        improved_delays.append(max(0, cp_imp - 4))
    return {
        "n_patients": n_patients,
        "mean_delay_plain_bcd": float(np.mean(plain_delays)),
        "mean_delay_improved_fixed": float(np.mean(improved_delays)),
    }


__all__ = ["run_final", "run_final_fixed"]
