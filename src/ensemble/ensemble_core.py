"""Real implementation for ensemble_core.py.

Original location: src/agent6-ensemble/ensemble_core.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent6-ensemble/` directory that was never committed to git.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np


@dataclass
class EnsembleMember:
    """One base learner in the MAML-style ensemble."""
    name: str
    predict_fn: Callable[[np.ndarray], np.ndarray]
    task_weights: dict[str, np.ndarray] = field(default_factory=dict)

    def predict(self, X: np.ndarray, task: str | None = None) -> np.ndarray:
        out = self.predict_fn(X)
        if task and task in self.task_weights:
            out = out * self.task_weights[task]
        return out


class StackedEnsemble:
    """A minimal MAML-style stacked ensemble.

    Each ``EnsembleMember`` has per-task adaptive weights. The ensemble
    combines member predictions via a softmax over learned weights. This is
    the architecture described in README §"src/ensemble/ — MAML meta-learning".
    """

    def __init__(self, members: Sequence[EnsembleMember]):
        if not members:
            raise ValueError("members must be non-empty")
        self.members: list[EnsembleMember] = list(members)
        self.global_weights: np.ndarray = np.ones(len(self.members)) / len(self.members)
        self.task_global_weights: dict[str, np.ndarray] = {}

    def fit_global_weights(self, X: np.ndarray, y: np.ndarray,
                           n_iter: int = 50, lr: float = 0.1) -> None:
        """One-step MAML-style adaptation of combination weights."""
        preds = np.stack([m.predict_fn(X) for m in self.members], axis=0)  # (M, N)
        # preds shape (M, N, C) or (M, N); flatten last dims for residual.
        if preds.ndim == 3:
            preds_flat = preds.reshape(len(self.members), X.shape[0], -1)
        else:
            preds_flat = preds[:, :, None]
        y_exp = y.reshape(-1, 1) if y.ndim == 1 else y
        if y_exp.ndim == 1:
            y_exp = y_exp[:, None]
        for _ in range(n_iter):
            w = np.exp(self.global_weights) / np.exp(self.global_weights).sum()
            combined = np.tensordot(w, preds_flat, axes=([0], [0]))  # (N, C)
            err = combined - y_exp
            grad = np.zeros(len(self.members))
            for i, p in enumerate(preds_flat):
                if p.ndim == 2:
                    grad[i] = float(np.sum(err * p))
                else:
                    grad[i] = float(np.sum(err * p.reshape(p.shape[0], -1)))
            self.global_weights -= lr * grad / max(X.shape[0], 1)

    def predict(self, X: np.ndarray, task: str | None = None) -> np.ndarray:
        preds = [m.predict(X, task) for m in self.members]
        if task and task in self.task_global_weights:
            w = self.task_global_weights[task]
        else:
            w = np.exp(self.global_weights) / np.exp(self.global_weights).sum()
        stacked = np.stack(preds, axis=0)
        return np.tensordot(w, stacked.reshape(stacked.shape[0], X.shape[0], -1),
                            axes=([0], [0])).reshape(X.shape[0], -1)

    def set_task_weights(self, task: str, weights: np.ndarray) -> None:
        if len(weights) != len(self.members):
            raise ValueError("weights length must equal member count")
        self.task_global_weights[task] = np.asarray(weights, dtype=float)


class TwoStageScreener:
    """Stage 1 ensemble → Stage 2 triage threshold.

    Stage 1 produces a continuous score. Stage 2 maps the score to a binary
    call by comparing against a fixed specificity target, learned from a
    held-out cohort.
    """

    def __init__(self, stage1: StackedEnsemble, target_specificity: float = 0.99):
        self.stage1 = stage1
        self.target_specificity = target_specificity
        self._threshold: float | None = None

    def fit_threshold(self, X: np.ndarray, y_neg: np.ndarray) -> None:
        """Pick threshold so that P(pred >= t | y=0) ≈ target_specificity."""
        scores = self.stage1.predict(X)
        if scores.ndim == 2:
            scores = scores[:, 0]
        neg_scores = self.stage1.predict(y_neg)
        if neg_scores.ndim == 2:
            neg_scores = neg_scores[:, 0]
        cutoff = np.quantile(neg_scores, self.target_specificity)
        self._threshold = float(cutoff)

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self._threshold is None:
            raise RuntimeError("call fit_threshold() before predict()")
        s = self.stage1.predict(X)
        if s.ndim == 2:
            s = s[:, 0]
        return (s >= self._threshold).astype(int)


__all__ = ["EnsembleMember", "StackedEnsemble", "TwoStageScreener"]
