"""Real implementation for synthetic_data_fixed.py.

Original location: src/agent2-multimodal-fusion/synthetic_data_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent2-multimodal-fusion/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version alias of :mod:`synthetic_data`.

Earlier version did not clip the SNV VAF distribution, occasionally
producing biologically implausible alt-freq > 1.0. Fixed by clipping in
the upstream data generator and re-exporting here.
"""
from .synthetic_data import aligned_sample
import numpy as np


def aligned_sample_clipped(is_cancer: bool, seed: int = 0,
                            n_fragments: int = 1000):
    s = aligned_sample(is_cancer, seed=seed, n_fragments=n_fragments)
    s.snv_alt_freq = np.clip(s.snv_alt_freq, 0, 1)
    return s


__all__ = ["aligned_sample", "aligned_sample_clipped"]
