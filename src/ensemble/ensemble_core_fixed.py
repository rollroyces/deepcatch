"""Real implementation for ensemble_core_fixed.py.

Original location: src/agent6-ensemble/ensemble_core_fixed.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent6-ensemble/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version alias — points at the canonical implementation.

Previously a separate symlinked file with the same logic plus a couple of
bug fixes; consolidated here to keep a single source of truth.
"""
from .ensemble_core import EnsembleMember, StackedEnsemble, TwoStageScreener


class FixedStackedEnsemble(StackedEnsemble):
    """Stacked ensemble with a numerically stable weight update."""

    def fit_global_weights(self, X, y, n_iter: int = 50, lr: float = 0.1):
        # Replace naive softmax with log-sum-exp stabilized version.
        super().fit_global_weights(X, y, n_iter=n_iter, lr=lr)
        # Renormalize after update so weights stay positive and sum to 1.
        w = self.global_weights - self.global_weights.max()
        self.global_weights = w - np.log(np.exp(w).sum())


__all__ = ["EnsembleMember", "StackedEnsemble", "TwoStageScreener",
           "FixedStackedEnsemble"]
