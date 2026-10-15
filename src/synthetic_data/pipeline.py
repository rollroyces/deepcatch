"""Real implementation for pipeline.py.

Original location: src/agent5-synthetic-data/pipeline.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""Top-level orchestration for synthetic cohort generation.

Composes the ``MultiModalDataGenerator`` with downstream feature extractors.
"""
from typing import Iterator

import numpy as np

from .data_generator import MultiModalDataGenerator, MultiModalSample


class SyntheticCohortPipeline:
    """Yields (X, y) pairs from a cancer/healthy mixture."""

    def __init__(self, n_cancer: int = 200, n_healthy: int = 800,
                 generator: MultiModalDataGenerator | None = None,
                 seed: int = 0):
        self.n_cancer = n_cancer
        self.n_healthy = n_healthy
        self.generator = generator or MultiModalDataGenerator(seed=seed)

    def __iter__(self) -> Iterator[tuple[np.ndarray, int]]:
        for _ in range(self.n_cancer):
            s = self.generator.sample_cancer()
            yield self._featurise(s), 1
        for _ in range(self.n_healthy):
            s = self.generator.sample_healthy()
            yield self._featurise(s), 0

    @staticmethod
    def _featurise(sample: MultiModalSample) -> np.ndarray:
        # Concatenate fragment-size histogram + mean methylation + mean SNV VAF
        # into a fixed-length feature vector.
        hist, _ = np.histogram(sample.fragment_sizes, bins=range(50, 500, 10))
        return np.concatenate([
            hist.astype(float),
            sample.methylation_beta.mean(axis=0),
            [sample.snv_alt_freq.mean()],
        ])


__all__ = ["SyntheticCohortPipeline"]
