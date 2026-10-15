"""Real implementation for run_final.py.

Original location: src/agent3-longitudinal/run_final.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Final longitudinal evaluation driver.

Compares the plain BCD vs. the improved (BCD + Kalman) detector across a
synthetic cohort and prints mean detection delay.
"""
import numpy as np

from .simulation import simulate_trajectory
from .bayesian_changepoint import BayesianChangepointDetector
from .improved_methods import ImprovedDetector


def run_final(n_patients: int = 50, seed: int = 0) -> dict:
    bcd = BayesianChangepointDetector()
    imp = ImprovedDetector()
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
        "mean_delay_improved": float(np.mean(improved_delays)),
    }


if __name__ == "__main__":  # pragma: no cover
    import json
    print(json.dumps(run_final(), indent=2))


__all__ = ["run_final"]
