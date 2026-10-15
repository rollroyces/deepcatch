"""Real implementation for evaluate_fixed.py.

Original location: src/agent1-variant-calling/evaluate_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version of :mod:`evaluate` — adds zero-division safety for empty call sets."""
from .evaluate import VariantEval, evaluate_calls


def safe_evaluate_calls(true_pos, called) -> VariantEval:
    out = evaluate_calls(true_pos, called)
    if out.n_true == 0 or out.n_called == 0:
        return VariantEval(precision=0.0, recall=0.0, f1=0.0,
                           n_true=out.n_true, n_called=out.n_called)
    return out


__all__ = ["VariantEval", "evaluate_calls", "safe_evaluate_calls"]
