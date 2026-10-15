"""Real implementation for baseline_callers.py.

Original location: src/agent1-variant-calling/baseline_callers.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Baseline variant callers — naive threshold-based, used as a sanity check
against the Bayesian caller.
"""
from dataclasses import dataclass

import numpy as np

from .bayesian_caller import Call


def simple_threshold_caller(ref_count: int, alt_count: int,
                            min_alt: int = 3, min_vaf: float = 0.05) -> Call | None:
    n = ref_count + alt_count
    if n == 0:
        return None
    vaf = alt_count / n
    if alt_count < min_alt or vaf < min_vaf:
        return None
    return Call(".", 0, "N", "N", ref_count, alt_count, vaf)


__all__ = ["simple_threshold_caller"]
