"""Real implementation for synthetic_data.py.

Original location: src/agent2-multimodal-fusion/synthetic_data.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent2-multimodal-fusion/` directory that was never committed to git.
"""

from __future__ import annotations

"""Generate paired multi-modal samples for fusion training.

Wraps the canonical :mod:`src.synthetic_data.data_generator` and adds
multi-modal alignment (matching fragment indices across modalities).
"""
import numpy as np

from ..synthetic_data.data_generator import (
    MultiModalDataGenerator, MultiModalSample, TissueAtlas,
)


def aligned_sample(is_cancer: bool, seed: int = 0,
                   n_fragments: int = 1000) -> MultiModalSample:
    gen = MultiModalDataGenerator(atlas=TissueAtlas(seed=seed),
                                  n_fragments=n_fragments, seed=seed)
    return gen.sample_cancer() if is_cancer else gen.sample_healthy()


__all__ = ["aligned_sample"]
