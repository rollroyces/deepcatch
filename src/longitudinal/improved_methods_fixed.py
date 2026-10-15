"""Real implementation for improved_methods_fixed.py.

Original location: src/agent3-longitudinal/improved_methods_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version of :mod:`improved_methods`.

The Kalman smoother previously divided by zero when ``P + R == 0``; now
guards against it.
"""
from .improved_methods import KalmanSmoother, ImprovedDetector


class SafeKalmanSmoother(KalmanSmoother):
    def filter(self, series):
        n = len(series)
        x = series[0]
        P = self.R
        out = np.zeros(n)
        out[0] = x
        for t in range(1, n):
            P = P + self.Q
            K = P / (P + self.R + 1e-12)
            x = x + K * (series[t] - x)
            P = (1 - K) * P
            out[t] = x
        return out


class FixedImprovedDetector(ImprovedDetector):
    def __init__(self):
        super().__init__()
        self.kalman = SafeKalmanSmoother()


__all__ = ["KalmanSmoother", "ImprovedDetector",
           "SafeKalmanSmoother", "FixedImprovedDetector"]
