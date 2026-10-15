"""Real implementation for bayesian_caller.py.

Original location: src/agent1-variant-calling/bayesian_caller.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Naive Bayesian variant caller — likelihood ratio test on read counts.

Returns ``P(variant | data)`` using a uniform prior and Beta-binomial
likelihoods for ref/alt allele counts.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Call:
    chrom: str
    pos: int
    ref: str
    alt: str
    ref_count: int
    alt_count: int
    p_variant: float


def posterior_variant(ref_count: int, alt_count: int,
                      prior: float = 0.5,
                      alpha_ref: float = 1.0, beta_ref: float = 50.0,
                      alpha_alt: float = 1.0, beta_alt: float = 50.0) -> float:
    """Posterior probability under a Beta-binomial model with uniform prior."""
    n = ref_count + alt_count
    if n == 0:
        return 0.0
    # Likelihood ratio of alt-model vs ref-model.
    # Both are beta-binomial with the same concentration; this reduces to
    # ratio of betas at the observed counts.
    from scipy.special import betaln
    log_lr_alt = (betaln(ref_count + alpha_ref, alt_count + alpha_alt)
                  - betaln(alpha_ref, alpha_alt))
    log_lr_ref = (betaln(ref_count + alpha_ref, n - ref_count + beta_ref)
                  - betaln(alpha_ref, beta_ref))
    # Convert to posterior assuming equal priors.
    log_post_alt = np.log(prior) + log_lr_alt
    log_post_ref = np.log(1 - prior) + log_lr_ref
    m = max(log_post_alt, log_post_ref)
    p_alt = np.exp(log_post_alt - m)
    p_ref = np.exp(log_post_ref - m)
    return float(p_alt / (p_alt + p_ref))


def call_variant(chrom: str, pos: int, ref: str, alt: str,
                 ref_count: int, alt_count: int,
                 threshold: float = 0.5) -> Call:
    p = posterior_variant(ref_count, alt_count)
    return Call(chrom, pos, ref, alt, ref_count, alt_count, p)


__all__ = ["Call", "posterior_variant", "call_variant"]
