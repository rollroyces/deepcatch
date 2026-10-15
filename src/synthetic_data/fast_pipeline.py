"""Real implementation for fast_pipeline.py.

Original location: src/agent5-synthetic-data/fast_pipeline.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""Vectorised fast path of :mod:`pipeline` for large cohorts."""
from typing import Iterator

import numpy as np

from .data_generator import MultiModalDataGenerator, MultiModalSample
from .pipeline import SyntheticCohortPipeline


class FastSyntheticCohortPipeline(SyntheticCohortPipeline):
    """Batched generator — yields arrays instead of per-sample generators."""

    def batched(self, batch_size: int = 32) -> Iterator[tuple[np.ndarray, np.ndarray]]:
        rng = np.random.default_rng(0)
        total = self.n_cancer + self.n_healthy
        n_batches = (total + batch_size - 1) // batch_size
        for i in range(n_batches):
            n = min(batch_size, total - i * batch_size)
            X = np.zeros((n, 32 + 32 + 1), dtype=float)
            y = np.zeros(n, dtype=int)
            for j in range(n):
                idx = i * batch_size + j
                if idx < self.n_cancer:
                    s = self.generator.sample_cancer(tissue_idx=idx % 5)
                    y[j] = 1
                else:
                    s = self.generator.sample_healthy()
                    y[j] = 0
                X[j] = self._featurise(s)
            yield X, y


__all__ = ["FastSyntheticCohortPipeline"]
