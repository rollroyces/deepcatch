"""Real implementation for data_generator.py.

Original location: src/agent5-synthetic-data/data_generator.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""MultiModalDataGenerator — synthesises (fragmentomics, methylation, SNV)
tuples per simulated fragment for cancer vs healthy plasma.
"""
from dataclasses import dataclass

import numpy as np

from .cfDNA_constants import size_distribution


@dataclass
class MultiModalSample:
    fragment_sizes: np.ndarray
    end_motifs: np.ndarray  # shape (n_fragments, 256)
    methylation_beta: np.ndarray  # shape (n_fragments, n_cpgs)
    snv_alt_freq: np.ndarray  # shape (n_fragments,)


class TissueAtlas:
    """Cancer vs healthy tissue-of-origin priors.

    Each tissue is parameterised by a per-cancer shift in fragment-size
    mean and a small set of methylation markers. The atlas is consumed by
    ``MultiModalDataGenerator`` to inject tissue-specific signals into
    synthetic cohorts.
    """

    def __init__(self, n_tissues: int = 5, seed: int = 0):
        self.n_tissues = n_tissues
        rng = np.random.default_rng(seed)
        # Per-tissue mean shift in bp relative to healthy.
        self.size_shift_bp = rng.normal(0, 8, size=n_tissues)
        # Per-tissue methylation shifts per CPG (small panel: 32).
        self.n_cpgs = 32
        self.methylation_shift = rng.normal(0, 0.15, size=(n_tissues, self.n_cpgs))


class MultiModalDataGenerator:
    """Produces synthetic cfDNA samples with cancer-vs-healthy signal."""

    def __init__(self, atlas: TissueAtlas | None = None,
                 n_fragments: int = 5000, seed: int = 0):
        self.atlas = atlas or TissueAtlas(seed=seed)
        self.n_fragments = n_fragments
        self.seed = seed

    def sample_healthy(self) -> MultiModalSample:
        return self._sample(is_cancer=False, tissue_idx=0)

    def sample_cancer(self, tissue_idx: int = 1) -> MultiModalSample:
        return self._sample(is_cancer=True, tissue_idx=tissue_idx % self.atlas.n_tissues)

    def _sample(self, is_cancer: bool, tissue_idx: int) -> MultiModalSample:
        rng = np.random.default_rng(self.seed + tissue_idx)
        # Fragment sizes: healthy baseline + tissue shift if cancer.
        sizes = size_distribution(self.n_fragments, seed=self.seed + tissue_idx)
        if is_cancer:
            sizes = np.clip(sizes + int(self.atlas.size_shift_bp[tissue_idx]), 50, 500)
        # End motifs: uniform 256-vector (placeholder).
        motifs = np.full((self.n_fragments, 256), 1.0 / 256)
        # Methylation: healthy ≈ 0.5, cancer shifts per-tissue.
        beta = rng.normal(0.5, 0.1, size=(self.n_fragments, self.atlas.n_cpgs))
        if is_cancer:
            beta += self.atlas.methylation_shift[tissue_idx]
        beta = np.clip(beta, 0, 1)
        # SNV alt freq: 0 for healthy, low VAF for cancer.
        snv = rng.normal(0.005 if is_cancer else 0.0, 0.002,
                         size=self.n_fragments).clip(0, 0.05) if is_cancer else np.zeros(self.n_fragments)
        return MultiModalSample(fragment_sizes=sizes,
                                end_motifs=motifs,
                                methylation_beta=beta,
                                snv_alt_freq=snv)


__all__ = ["MultiModalSample", "TissueAtlas", "MultiModalDataGenerator"]
