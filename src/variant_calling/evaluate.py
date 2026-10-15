"""Real implementation for evaluate.py.

Original location: src/agent1-variant-calling/evaluate.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Variant-caller evaluation: precision, recall, F1 on a synthetic call set."""
from dataclasses import dataclass

import numpy as np


@dataclass
class VariantEval:
    precision: float
    recall: float
    f1: float
    n_true: int
    n_called: int


def evaluate_calls(true_pos: np.ndarray,
                   called: np.ndarray) -> VariantEval:
    """Both inputs are boolean arrays over the same variant set."""
    tp = int((true_pos & called).sum())
    fp = int((~true_pos & called).sum())
    fn = int((true_pos & ~called).sum())
    prec = tp / max(tp + fp, 1)
    rec = tp / max(tp + fn, 1)
    f1 = 2 * prec * rec / max(prec + rec, 1e-12)
    return VariantEval(precision=prec, recall=rec, f1=f1,
                       n_true=int(true_pos.sum()), n_called=int(called.sum()))


__all__ = ["VariantEval", "evaluate_calls"]
