"""Real implementation for rl_sampler.py.

Original location: src/agent3-longitudinal/rl_sampler.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Contextual bandit sampler for choosing when to draw cfDNA next.

We use a simple epsilon-greedy policy over a discrete set of candidate
intervals (1, 3, 6, 12 months). The reward proxy is detection delay,
estimated from a Bayesian changepoint detector.
"""
import numpy as np

from .bayesian_changepoint import BayesianChangepointDetector


class EpsilonGreedySampler:
    INTERVALS_MONTHS = (1, 3, 6, 12)

    def __init__(self, epsilon: float = 0.1, seed: int = 0):
        self.epsilon = epsilon
        self.rng = np.random.default_rng(seed)
        self.values = np.zeros(len(self.INTERVALS_MONTHS))
        self.counts = np.zeros(len(self.INTERVALS_MONTHS))
        self.bcd = BayesianChangepointDetector()

    def choose(self) -> int:
        if self.rng.random() < self.epsilon:
            return int(self.rng.integers(0, len(self.INTERVALS_MONTHS)))
        return int(np.argmax(self.values))

    def update(self, idx: int, reward: float) -> None:
        self.counts[idx] += 1
        n = self.counts[idx]
        self.values[idx] += (reward - self.values[idx]) / n

    @staticmethod
    def reward_from_delay(delay: int) -> float:
        return float(max(0.0, 1.0 - delay / 6.0))


__all__ = ["EpsilonGreedySampler"]
