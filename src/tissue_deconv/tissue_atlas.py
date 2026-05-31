#!/usr/bin/env python3
"""
Tissue Methylation Atlas
==========================

Stores reference methylation profiles for 29 tissues. In the absence of
real cfSort reference data, generates synthetic reference profiles with
known tissue-specific methylation patterns for training and evaluation.

Reference
---------
Li et al. (2023) PNAS — The cfSort reference atlas contains methylation
profiles at tissue-specific CpG sites across 29 human tissues from
521 samples.

Two modes:
1. **Real mode**: Load pre-computed reference data from a path.
2. **Fallback mode**: Generate synthetic profiles with realistic
   tissue-specific differentially methylated regions (tDMRs).
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import SUPPORTED_TISSUES, CANCER_RELEVANT_TISSUES

logger = logging.getLogger(__name__)


class TissueAtlas:
    """
    Stores reference methylation profiles per tissue for deconvolution.

    Each tissue has a characteristic methylation profile at a set of
    tissue-discriminative CpG sites. These profiles serve as the
    reference basis for unmixing cfDNA methylation mixtures.

    Parameters
    ----------
    n_cpg_features : int
        Number of CpG markers for the reference profiles.
    random_seed : int
        Seed for reproducible synthetic profile generation.
    """

    def __init__(
        self,
        n_cpg_features: int = 1000,
        random_seed: int = 42,
    ):
        self.n_cpg_features = n_cpg_features
        self.random_seed = random_seed
        self.rng = np.random.RandomState(random_seed)

        # Reference profiles: {tissue_name -> np.array(n_cpg_features,)}
        self._reference: Dict[str, np.ndarray] = {}
        self._is_synthetic: bool = False

        # Marker information for each tissue
        self._marker_indices: Dict[str, np.ndarray] = {}

    @property
    def tissues(self) -> List[str]:
        """List of tissues with loaded reference profiles."""
        return list(self._reference.keys())

    @property
    def is_synthetic(self) -> bool:
        """Whether the atlas was generated synthetically."""
        return self._is_synthetic

    @property
    def is_loaded(self) -> bool:
        """Whether reference data has been loaded."""
        return len(self._reference) > 0

    # ── Reference loading ────────────────────────────────────────

    def load_reference(
        self,
        path: Optional[str] = None,
        reference_dict: Optional[Dict[str, np.ndarray]] = None,
    ) -> None:
        """
        Load pre-computed reference methylation profiles.

        Parameters
        ----------
        path : str, optional
            Path to a NumPy archive (.npz) with keys matching tissue names.
        reference_dict : dict, optional
            Dict of {tissue_name -> (n_cpg_features,) beta_values array}.
            Takes precedence over path.
        """
        if reference_dict is not None:
            self._validate_reference_dict(reference_dict)
            self._reference = dict(reference_dict)
            self._is_synthetic = False
            self._compute_markers()
            logger.info(
                "Loaded reference atlas with %d tissues from dict",
                len(self._reference),
            )
            return

        if path is not None:
            try:
                data = np.load(path, allow_pickle=False)
                ref = {t: data[t] for t in SUPPORTED_TISSUES if t in data}
                self._validate_reference_dict(ref)
                self._reference = ref
                self._is_synthetic = False
                self._compute_markers()
                logger.info(
                    "Loaded reference atlas from %s (%d tissues)",
                    path, len(self._reference),
                )
                return
            except (FileNotFoundError, OSError, KeyError) as e:
                logger.warning(
                    "Failed to load reference from %s: %s. Falling back to synthetic.",
                    path, e,
                )

        # Fallback: generate synthetic profiles
        self._generate_synthetic_reference()
        self._is_synthetic = True
        logger.info(
            "Generated synthetic reference atlas with %d tissues × %d CpGs",
            len(self._reference), self.n_cpg_features,
        )

    def _validate_reference_dict(self, ref: Dict[str, np.ndarray]) -> None:
        """Validate reference dictionary format."""
        if not ref:
            raise ValueError("Reference dictionary is empty")
        expected_shape = (self.n_cpg_features,)
        for tissue, profile in ref.items():
            if not isinstance(profile, np.ndarray):
                raise TypeError(
                    f"Profile for '{tissue}' must be np.ndarray, got {type(profile)}"
                )
            if profile.shape != expected_shape:
                raise ValueError(
                    f"Profile for '{tissue}' has shape {profile.shape}, "
                    f"expected {expected_shape}"
                )
            if np.any((profile < 0) | (profile > 1)):
                raise ValueError(
                    f"Profile for '{tissue}' contains values outside [0, 1]"
                )

    def to_npz(self, path: str) -> None:
        """
        Save reference atlas to a compressed NumPy archive.

        Parameters
        ----------
        path : str
            Output .npz file path.
        """
        if not self._reference:
            raise RuntimeError("No reference data to save. Load or generate first.")
        np.savez_compressed(path, **self._reference)

    # ── Synthetic reference generation ───────────────────────────

    def _generate_synthetic_reference(self) -> None:
        """
        Generate synthetic tissue-specific methylation profiles.

        Strategy: Each tissue gets a baseline methylation pattern plus
        tissue-specific differentially methylated regions (tDMRs).
        Blood-derived tissues share high similarity; solid organs
        have distinct patterns.

        This creates realistic synthetic data for model training
        and evaluation when real cfSort atlas is unavailable.
        """
        n_tissues = len(SUPPORTED_TISSUES)
        n_cpg = self.n_cpg_features

        # ── Global baseline: average methylome profile ───────────
        # Most CpGs have bimodal methylation (mostly 0 or 1)
        baseline = np.zeros(n_cpg)
        # ~70% CpGs are highly methylated baseline, ~30% unmethylated
        n_high = int(0.7 * n_cpg)
        high_indices = self.rng.choice(n_cpg, size=n_high, replace=False)
        baseline[high_indices] = 0.85 + self.rng.beta(2, 2, size=n_high) * 0.15
        low_mask = ~np.isin(np.arange(n_cpg), high_indices)
        baseline[low_mask] = self.rng.beta(1, 5, size=low_mask.sum()) * 0.15

        # ── Tissue-specific marker CpGs ──────────────────────────
        # Each tissue has n_markers CpGs where its profile deviates
        # from baseline
        n_markers_per_tissue = max(5, n_cpg // len(SUPPORTED_TISSUES))

        # Allocate non-overlapping marker regions when possible
        markers_per_tissue: Dict[str, np.ndarray] = {}
        used_indices: set = set()

        for tissue in SUPPORTED_TISSUES:
            available = list(set(range(n_cpg)) - used_indices)
            if len(available) >= n_markers_per_tissue:
                chosen = self.rng.choice(
                    available, size=n_markers_per_tissue, replace=False
                )
            else:
                # Fallback: allow overlap if not enough unique CpGs
                chosen = self.rng.choice(
                    n_cpg, size=n_markers_per_tissue, replace=False
                )
            markers_per_tissue[tissue] = chosen
            used_indices.update(chosen)

        # ── Build profiles ───────────────────────────────────────
        profiles: Dict[str, np.ndarray] = {}

        for tissue in SUPPORTED_TISSUES:
            profile = baseline.copy()
            markers = markers_per_tissue[tissue]

            # Tissue-specific deviation from baseline
            if tissue in BLOOD_DERIVED():
                # Blood-derived: small deviations (healthy background)
                deviation = self.rng.normal(0, 0.02, size=len(markers))
            else:
                # Solid organs: larger, tissue-specific deviations
                # Half go hypermethylated, half hypomethylated relative to baseline
                n_half = len(markers) // 2
                deviation = np.zeros(len(markers))
                deviation[:n_half] = self.rng.uniform(0.1, 0.4, size=n_half)
                deviation[n_half:] = self.rng.uniform(-0.4, -0.1, size=len(markers) - n_half)

            profile[markers] = np.clip(profile[markers] + deviation, 0.0, 1.0)

            # Add small noise for biological variability
            noise = self.rng.normal(0, 0.01, size=n_cpg)
            profile = np.clip(profile + noise, 0.0, 1.0)

            profiles[tissue] = profile.astype(np.float32)

        self._reference = profiles
        self._marker_indices = {
            t: markers_per_tissue[t] for t in SUPPORTED_TISSUES
        }

    # ── Marker CpG retrieval ─────────────────────────────────────

    def _compute_markers(self) -> None:
        """Compute tissue-discriminative marker CpG indices from profiles."""
        if not self._reference:
            return

        # Only use tissues that are in the reference
        available_tissues = [t for t in SUPPORTED_TISSUES if t in self._reference]
        if not available_tissues:
            return

        all_profiles = np.stack(
            [self._reference[t] for t in available_tissues]
        )  # (n_available, n_cpg)
        baseline = np.mean(all_profiles, axis=0)

        n_markers = max(5, self.n_cpg_features // max(1, len(available_tissues)))

        for tissue in available_tissues:
            profile = self._reference[tissue]
            # Find CpGs where this tissue deviates most from baseline
            deviation = np.abs(profile - baseline)
            top_indices = np.argsort(deviation)[-n_markers:]
            self._marker_indices[tissue] = top_indices

    def get_marker_cpgs(
        self,
        tissue: str,
        n_top: int = 100,
    ) -> np.ndarray:
        """
        Get top tissue-discriminative CpG indices for a tissue.

        Parameters
        ----------
        tissue : str
            Tissue name from SUPPORTED_TISSUES.
        n_top : int
            Number of top markers to return.

        Returns
        -------
        cpgs : (n_top,) int array
            Indices of top discriminative CpGs for this tissue.
        """
        if tissue not in self._reference:
            raise ValueError(
                f"Tissue '{tissue}' not in reference atlas. "
                f"Available: {list(self._reference.keys())}"
            )

        if tissue in self._marker_indices:
            markers = self._marker_indices[tissue]
            return markers[:n_top] if len(markers) >= n_top else markers
        else:
            # Fallback: compute on the fly
            available = list(self._reference.keys())
            all_profiles = np.stack(
                [self._reference[t] for t in available]
            )
            baseline = np.mean(all_profiles, axis=0)
            deviation = np.abs(self._reference[tissue] - baseline)
            top = np.argsort(deviation)[-min(n_top, len(deviation)):]
            return top

    def get_reference_profile(self, tissue: str) -> np.ndarray:
        """
        Get the reference methylation profile for a tissue.

        Parameters
        ----------
        tissue : str
            Tissue name.

        Returns
        -------
        profile : (n_cpg_features,) float32 array
        """
        if tissue not in self._reference:
            raise KeyError(f"Tissue '{tissue}' not found in reference atlas")
        return self._reference[tissue].copy()

    def get_all_profiles_matrix(self) -> np.ndarray:
        """
        Get all reference profiles as a matrix.

        Returns
        -------
        profiles : (n_tissues, n_cpg_features) float32 array
            Each row is a tissue reference profile.
        """
        return np.stack(
            [self._reference[t] for t in SUPPORTED_TISSUES]
        ).astype(np.float32)

    # ── In silico mixture generation ──────────────────────────────

    def generate_synthetic_mixture(
        self,
        tissue_fractions: Dict[str, float],
        noise: float = 0.01,
    ) -> np.ndarray:
        """
        Create a synthetic cfDNA methylation mixture given tissue fractions.

        The mixture is a weighted sum of tissue-specific reference profiles,
        simulating the methylation signature of cfDNA originating from
        different tissues in proportion to their cell death rates.

        Parameters
        ----------
        tissue_fractions : dict
            Mapping from tissue name to fraction (should sum to ~1.0).
        noise : float
            Gaussian noise standard deviation to add (simulates
            sequencing and biological noise).

        Returns
        -------
        mixture : (n_cpg_features,) float32 array
            Synthetic cfDNA methylation beta-values.
        """
        if not self._reference:
            self.load_reference()

        total = sum(tissue_fractions.values())
        if abs(total - 1.0) > 0.05:
            logger.warning(
                "Tissue fractions sum to %.4f (expected ~1.0). Normalizing.",
                total,
            )
            tissue_fractions = {
                k: v / total for k, v in tissue_fractions.items()
            }

        mixture = np.zeros(self.n_cpg_features, dtype=np.float32)

        for tissue, fraction in tissue_fractions.items():
            if fraction <= 0:
                continue
            if tissue not in self._reference:
                logger.warning(
                    "Tissue '%s' not in reference atlas, skipping.", tissue
                )
                continue
            mixture += fraction * self._reference[tissue]

        # Add noise
        if noise > 0:
            mixture += self.rng.normal(0, noise, size=self.n_cpg_features)

        return np.clip(mixture, 0.0, 1.0).astype(np.float32)

    def generate_training_dataset(
        self,
        n_samples: int = 1000,
        noise: float = 0.01,
        max_active_tissues: int = 5,
        blood_fraction_range: Tuple[float, float] = (0.70, 0.99),
        cancer_tissue_prob: float = 0.4,
        seed: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate a synthetic training dataset of mixtures and fractions.

        Each sample is a random mixture of tissues with:
        - Blood-derived tissues dominating (70-99% of cfDNA)
        - 0-5 solid organ tissues at trace levels
        - Some samples with elevated cancer-relevant tissue fractions

        Parameters
        ----------
        n_samples : int
            Number of synthetic mixtures to generate.
        noise : float
            Measurement noise standard deviation.
        max_active_tissues : int
            Max number of non-blood tissues active per sample.
        blood_fraction_range : tuple
            (min, max) range for total blood-derived fraction.
        cancer_tissue_prob : float
            Probability that a cancer-relevant tissue gets elevated
            fraction (simulating tumor shedding).
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        mixtures : (n_samples, n_cpg_features) float32 array
        fractions : (n_samples, n_tissues) float32 array
        """
        if not self._reference:
            self.load_reference()

        rng = np.random.RandomState(seed if seed is not None else self.random_seed)

        blood_tissues = BLOOD_DERIVED()
        solid_tissues = [
            t for t in SUPPORTED_TISSUES if t not in blood_tissues
        ]
        cancer_tissues = [t for t in CANCER_RELEVANT_TISSUES if t in solid_tissues]

        mixtures = np.zeros((n_samples, self.n_cpg_features), dtype=np.float32)
        fractions = np.zeros((n_samples, len(SUPPORTED_TISSUES)), dtype=np.float32)

        for i in range(n_samples):
            frac = np.zeros(len(SUPPORTED_TISSUES), dtype=np.float32)

            # ── Blood fraction: 70-99% ───────────────────────────
            total_blood = rng.uniform(*blood_fraction_range)
            remaining = 1.0 - total_blood

            # Distribute blood fraction across blood-derived tissues
            n_blood_active = rng.randint(3, min(len(blood_tissues) + 1, 8))
            active_blood = rng.choice(blood_tissues, size=n_blood_active, replace=False)
            blood_weights = rng.dirichlet(np.ones(n_blood_active))
            for t_name, w in zip(active_blood, blood_weights):
                idx = SUPPORTED_TISSUES.index(t_name)
                frac[idx] = w * total_blood

            # ── Solid organ fractions: ~1-30% ────────────────────
            n_solid = rng.randint(0, max_active_tissues + 1)
            if n_solid > 0 and remaining > 0.001:
                active_solid = rng.choice(solid_tissues, size=n_solid, replace=False)
                solid_weights = rng.dirichlet(np.ones(n_solid))

                for t_name, w in zip(active_solid, solid_weights):
                    # Cancer-relevant tissues may get boosted
                    if t_name in cancer_tissues and rng.random() < cancer_tissue_prob:
                        boost = rng.uniform(1.5, 4.0)
                        w *= boost

                    idx = SUPPORTED_TISSUES.index(t_name)
                    frac[idx] = w * remaining

            # Renormalize
            total = frac.sum()
            if total > 0:
                frac = frac / total

            fractions[i] = frac

            # ── Generate mixture ─────────────────────────────────
            tissue_dict = {
                SUPPORTED_TISSUES[j]: float(frac[j])
                for j in range(len(SUPPORTED_TISSUES))
                if frac[j] > 0
            }
            mixtures[i] = self.generate_synthetic_mixture(
                tissue_dict, noise=noise
            )

        return mixtures, fractions


# ── Module-level helper ──────────────────────────────────────────

_BLOOD_DERIVED_CACHE: Optional[List[str]] = None


def BLOOD_DERIVED() -> List[str]:
    """Get blood-derived tissue names (module-level singleton)."""
    from .config import BLOOD_DERIVED_TISSUES
    return BLOOD_DERIVED_TISSUES
