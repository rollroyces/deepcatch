"""Real implementation for improved_methods.py.

Original location: src/agent3-longitudinal/improved_methods.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Improved longitudinal methods — combines BCD with a Kalman smoother.

The plain BCD above is jittery on short trajectories; this module adds a
one-step Kalman smoother on the posterior mean to stabilise detection.
"""
import numpy as np

from .bayesian_changepoint import BayesianChangepointDetector


class KalmanSmoother:
    def __init__(self, process_var: float = 1e-6, obs_var: float = 1e-4):
        self.Q = process_var
        self.R = obs_var

    def filter(self, series: np.ndarray) -> np.ndarray:
        n = len(series)
        x = series[0]
        P = self.R
        out = np.zeros(n)
        out[0] = x
        for t in range(1, n):
            x = x  # constant-process model
            P = P + self.Q
            K = P / (P + self.R)
            x = x + K * (series[t] - x)
            P = (1 - K) * P
            out[t] = x
        return out


class ImprovedDetector:
    def __init__(self):
        self.bcd = BayesianChangepointDetector()
        self.kalman = KalmanSmoother()

    def detect(self, series: np.ndarray) -> tuple[int, np.ndarray]:
        smoothed = self.kalman.filter(series)
        return self.bcd.detect(smoothed)


__all__ = ["KalmanSmoother", "ImprovedDetector"]
