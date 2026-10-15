"""Real implementation for simulate.py.

Original location: src/agent1-variant-calling/simulate.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Simulate read counts at a set of candidate variant sites.

Returns a tuple ``(reads, true_pos)`` where ``reads`` is a list of
``(ref_count, alt_count)`` per site and ``true_pos`` is a boolean mask
indicating which sites are true variants.
"""
import numpy as np


def simulate_reads(n_variants: int = 100, n_sites: int = 1000,
                   coverage: int = 30, seed: int = 0) -> tuple[list, np.ndarray]:
    rng = np.random.default_rng(seed)
    variant_sites = rng.choice(n_sites, size=min(n_variants, n_sites), replace=False)
    true_pos = np.zeros(n_sites, dtype=bool)
    true_pos[variant_sites] = True
    reads = []
    for i in range(n_sites):
        if true_pos[i]:
            alt_frac = rng.uniform(0.05, 0.5)
        else:
            alt_frac = 0.0
        n = rng.poisson(coverage)
        alt = rng.binomial(n, alt_frac)
        ref = max(0, n - alt)
        reads.append((int(ref), int(alt)))
    return reads, true_pos


__all__ = ["simulate_reads"]
