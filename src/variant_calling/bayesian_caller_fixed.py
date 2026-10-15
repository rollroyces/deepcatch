"""Real implementation for bayesian_caller_fixed.py.

Original location: src/agent1-variant-calling/bayesian_caller_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version of :mod:`bayesian_caller` — uses Beta-binomial with a
slightly informative prior and numerical log-sum-exp normalisation.
"""
import numpy as np
from scipy.special import betaln

from .bayesian_caller import Call


def posterior_variant_fixed(ref_count: int, alt_count: int,
                            prior: float = 0.01) -> float:
    n = ref_count + alt_count
    if n == 0:
        return 0.0
    a_alt, b_alt = 0.5, 50.0
    a_ref, b_ref = 50.0, 0.5
    log_post_alt = (np.log(prior)
                    + betaln(ref_count + a_alt, alt_count + b_alt)
                    - betaln(a_alt, b_alt))
    log_post_ref = (np.log(1 - prior)
                    + betaln(ref_count + a_ref, alt_count + b_ref)
                    - betaln(a_ref, b_ref))
    m = max(log_post_alt, log_post_ref)
    pa = np.exp(log_post_alt - m)
    pr = np.exp(log_post_ref - m)
    return float(pa / (pa + pr))


def call_variant_fixed(chrom: str, pos: int, ref: str, alt: str,
                        ref_count: int, alt_count: int,
                        threshold: float = 0.5) -> Call:
    p = posterior_variant_fixed(ref_count, alt_count)
    return Call(chrom, pos, ref, alt, ref_count, alt_count, p)


__all__ = ["Call", "posterior_variant_fixed", "call_variant_fixed"]
