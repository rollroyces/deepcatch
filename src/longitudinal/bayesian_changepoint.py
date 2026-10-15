"""Real implementation for bayesian_changepoint.py.

Original location: src/agent3-longitudinal/bayesian_changepoint.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Bayesian online changepoint detection (Adams & MacKay, 2007).

Used to detect treatment response / progression in serial cfDNA draws.
"""
import numpy as np


class BayesianChangepointDetector:
    def __init__(self, hazard: float = 1 / 5.0,
                 prior_mean: float = 0.005, prior_var: float = 1e-4):
        self.hazard = hazard
        self.prior_mean = prior_mean
        self.prior_var = prior_var

    def detect(self, series: np.ndarray) -> tuple[int, np.ndarray]:
        """Return (argmax, run-length posterior map).

        ``argmax`` is the time step with the highest posterior probability of
        being a changepoint.
        """
        n = len(series)
        log_run_length = np.full(n + 1, -np.inf)
        log_run_length[0] = 0.0
        maxes = np.zeros(n + 1)
        for t in range(1, n + 1):
            pred_mean = self.prior_mean
            pred_var = self.prior_var
            log_p = -0.5 * ((series[t - 1] - pred_mean) ** 2 / pred_var
                            + np.log(2 * np.pi * pred_var))
            growth = log_run_length[:t] + log_p
            maxes[t] = growth.max()
            log_run_length[1:t + 1] = growth[:t] + np.log(1 - self.hazard)
            log_run_length[t] = (log_run_length[:t] + log_p
                                 + np.log(self.hazard)).max()
        post = np.exp(log_run_length - log_run_length.max())
        post /= post.sum() + 1e-12
        return int(np.argmax(post[1:])), post[1:]


__all__ = ["BayesianChangepointDetector"]
