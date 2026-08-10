#!/usr/bin/env python3
"""
Multi-Modal Data Generation for Foundation Model Training
==========================================================

Simulates cfDNA multi-modal data for self-supervised pre-training.

Public data used for inspiration:
- TCGA methylation (Illumina 450K/EPIC arrays)
- GEO cfDNA methylation profiles (GSE122126, GSE110729)
- DELFI fragmentomics from Cristiano et al. (2019) Nature

Training data is synthetic but designed to match real cfDNA distributions:
- Correlated modalities (cancer affects multiple signals simultaneously)
- Realistic noise profiles (technical/stochastic variation in cfDNA)
- Class-imbalanced sampling (mimics real screening populations)
"""

from __future__ import annotations

import hashlib
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── Modality definitions ────────────────────────────────────────

MODALITY_DIMS: Dict[str, int] = {
    "frag_basic": 4,       # MFR, FSI, CAFF, FEM aggregated scores
    "frag_enhanced": 44,   # DELFI + MFS + nucleosome + 5-mer features
    "cnv": 6,              # chromosomal instability features
    "sero": 4,             # PG-I, PG-II, G-17, H. pylori
    "gnn": 1,              # field_defect_score from GATv2
    "tissue": 24,          # tissue deconvolution features
}

MODALITY_NAMES: Tuple[str, ...] = tuple(MODALITY_DIMS.keys())

# ── Reference ranges for realistic simulation ───────────────────

# Typical healthy ranges for each modality (mean, std)
HEALTHY_RANGES: Dict[str, Tuple[float, float]] = {
    "frag_basic": (0.3, 0.15),     # Healthy FSI ≈ 0.3
    "frag_enhanced": (0.0, 0.5),   # Z-scored enhanced features
    "cnv": (0.0, 0.3),             # CNA scores (z-scored)
    "sero": (0.5, 0.25),           # Normal serological markers
    "gnn": (0.0, 0.15),            # Low field defect score
    "tissue": (0.0, 0.2),          # Low tissue scores
}

# Cancer-related shifts for each modality (mean shift, additional std)
CANCER_SHIFTS: Dict[str, Tuple[float, float]] = {
    "frag_basic": (0.4, 0.15),
    "frag_enhanced": (0.5, 0.4),
    "cnv": (0.6, 0.3),
    "sero": (0.3, 0.2),
    "gnn": (0.5, 0.2),
    "tissue": (0.4, 0.25),
}


class MultiModalDataGenerator:
    """
    Generate synthetic multi-modal cfDNA data for foundation model
    pre-training and downstream fine-tuning.

    Modes:
    - 'pretrain': Generate unlabeled data for self-supervised training.
    - 'downstream': Generate labeled data (cancer vs healthy) for fine-tuning.
    - 'too': Generate tissue-of-origin labeled data.

    Parameters
    ----------
    seed : int
        Random seed for reproducibility.
    noise_level : float
        Amount of Gaussian noise added to features (0 = no noise,
        1 = high noise).
    modality_correlation : float
        Strength of inter-modality correlation (0 = independent,
        1 = perfectly correlated).
    """

    def __init__(
        self,
        seed: int = 42,
        noise_level: float = 0.1,
        modality_correlation: float = 0.5,
    ):
        self.rng = np.random.RandomState(seed)
        self.noise_level = noise_level
        self.modality_correlation = modality_correlation

    def _deterministic_key(self, sample_id: int, modality_name: str) -> int:
        """Generate a deterministic hash for consistent feature assignment."""
        msg = f"{sample_id}:{modality_name}".encode()
        return int(hashlib.md5(msg).hexdigest()[:8], 16)

    def generate_single_sample(
        self,
        sample_id: int,
        is_cancer: Optional[bool] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Generate multi-modal features for one sample.

        Parameters
        ----------
        sample_id : int
            Unique sample identifier (determines feature values).
        is_cancer : bool or None
            If None, randomly assign (50/50 split).

        Returns
        -------
        modalities : dict[str, ndarray]
            Features for each modality.
        """
        # Use deterministic RNG for this sample
        key = self._deterministic_key(sample_id, "flag")
        sample_rng = np.random.RandomState(key)

        # Determine label
        if is_cancer is None:
            is_cancer = sample_rng.random() > 0.7  # 30% cancer (realistic prevalence)

        modalities = {}

        for name, dim in MODALITY_DIMS.items():
            healthy_mean, healthy_std = HEALTHY_RANGES[name]
            cancer_shift_mean, cancer_shift_std = CANCER_SHIFTS[name]

            # Base healthy features
            base_key = self._deterministic_key(sample_id, name)
            feat_rng = np.random.RandomState(base_key)

            features = feat_rng.randn(dim) * healthy_std + healthy_mean

            if is_cancer:
                # Cancer shift with correlation across modalities
                # The same sample_id → consistent shift direction
                shift_rng = np.random.RandomState(
                    self._deterministic_key(sample_id, "cancer_shift")
                )
                shift_magnitude = shift_rng.randn() * cancer_shift_std + cancer_shift_mean
                # Apply shift (positive for most cancer signals)
                direction = feat_rng.choice([-1, 1], size=dim)
                features = features + direction * shift_magnitude * 0.3

            # Add modality correlation (shared component across modalities)
            if self.modality_correlation > 0:
                shared_key = self._deterministic_key(sample_id, "shared")
                shared_rng = np.random.RandomState(shared_key)
                shared_component = shared_rng.randn(dim) * self.modality_correlation * 0.1
                features = features + shared_component

            # Add technical noise
            noise = sample_rng.randn(dim) * self.noise_level * abs(healthy_std)
            features = features + noise

            # Zero out some features occasionally (missing data simulation)
            missing_mask = sample_rng.random(dim) < 0.02  # 2% missing rate
            features[missing_mask] = 0.0

            # Store
            modalities[name] = features.astype(np.float32)

        return modalities

    def generate_batch(
        self,
        batch_size: int,
        start_id: int = 0,
        is_cancer: Optional[np.ndarray] = None,
    ) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
        """
        Generate a batch of multi-modal samples.

        Parameters
        ----------
        batch_size : int
            Number of samples.
        start_id : int
            Starting sample ID.
        is_cancer : (batch_size,) bool array or None
            Labels (None = random).

        Returns
        -------
        modalities : dict[str, (batch_size, dim) ndarray]
            Batched modality features.
        labels : (batch_size,) int64 array
            1 = cancer, 0 = healthy.
        """
        batch_modalities: Dict[str, List[np.ndarray]] = {
            name: [] for name in MODALITY_NAMES
        }
        labels_list = []

        for i in range(batch_size):
            sample_id = start_id + i
            label = is_cancer[i] if is_cancer is not None else None

            sample = self.generate_single_sample(
                sample_id=sample_id,
                is_cancer=label,
            )

            for name in MODALITY_NAMES:
                batch_modalities[name].append(sample[name])

            if is_cancer is not None:
                labels_list.append(1 if is_cancer[i] else 0)

        # Stack batches
        batched = {
            name: np.stack(vals, axis=0)
            for name, vals in batch_modalities.items()
        }

        labels = np.array(labels_list, dtype=np.int64) if labels_list else np.array([])

        return batched, labels

    def generate_dataset(
        self,
        n_samples: int,
        prefix: str = "synthetic",
        cancer_prevalence: float = 0.3,
    ) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
        """
        Generate a full synthetic dataset.

        Parameters
        ----------
        n_samples : int
            Total number of samples.
        prefix : str
            Prefix for sample IDs (affects hash).
        cancer_prevalence : float
            Fraction of cancer samples.

        Returns
        -------
        modalities : dict[str, (n_samples, dim) ndarray]
        labels : (n_samples,) int64 array
        """
        n_cancer = int(n_samples * cancer_prevalence)
        is_cancer = np.zeros(n_samples, dtype=bool)
        is_cancer[:n_cancer] = True
        self.rng.shuffle(is_cancer)

        # Use a different prefix for hash so these are different from
        # pre-training data
        full_prefix = prefix
        all_modalities: Dict[str, List[np.ndarray]] = {
            name: [] for name in MODALITY_NAMES
        }

        for i in range(n_samples):
            # Deterministic sample id — built-in hash() is randomized per
            # process (PYTHONHASHSEED) and broke reproducibility across runs.
            sample_id = int(hashlib.md5(f"{full_prefix}_{i}".encode()).hexdigest()[:8], 16)
            sample = self.generate_single_sample(
                sample_id=sample_id,
                is_cancer=bool(is_cancer[i]),
            )
            for name in MODALITY_NAMES:
                all_modalities[name].append(sample[name])

        batched = {
            name: np.stack(vals, axis=0)
            for name, vals in all_modalities.items()
        }

        return batched, is_cancer.astype(np.int64)

    def generate_contrastive_pairs(
        self,
        batch_size: int,
        start_id: int = 0,
    ) -> Tuple[
        Dict[str, np.ndarray],   # anchor modalities
        Dict[str, np.ndarray],   # positive (same sample, same modalities)
        Dict[str, np.ndarray],   # negative (different sample)
    ]:
        """
        Generate contrastive pairs for pre-training.

        Returns:
            anchor_modalities: batch of samples
            positive_modalities: same samples (could add augmentations)
            negative_modalities: different (shuffled) samples
        """
        anchor, _ = self.generate_batch(batch_size, start_id)
        # Positive: same samples with slight noise augmentation
        positive = {
            name: arr + self.rng.randn(*arr.shape) * self.noise_level * 0.5
            for name, arr in anchor.items()
        }
        # Negative: different samples (shifted batch)
        neg_start = start_id + batch_size
        negative, _ = self.generate_batch(batch_size, neg_start)

        return anchor, positive, negative

    @staticmethod
    def modalities_to_array(
        modalities: Dict[str, np.ndarray],
    ) -> np.ndarray:
        """
        Convert modality dict to a flat array for simple fusion.

        Parameters
        ----------
        modalities : dict[str, ndarray]

        Returns
        -------
        stacked : (n_samples, total_dim) ndarray
        """
        arrays = [modalities[name] for name in MODALITY_NAMES]
        return np.column_stack(arrays)

    @staticmethod
    def modalities_to_scores(
        modalities: Dict[str, np.ndarray],
    ) -> List[np.ndarray]:
        """
        Convert modality dict to per-modality scalar scores
        (compatible with CrossAttentionFusion.fit interface).

        Each modality is reduced to a single scalar per sample
        (mean of features).

        Parameters
        ----------
        modalities : dict[str, ndarray]

        Returns
        -------
        scores : list of (n_samples,) arrays
        """
        scores = []
        for name in MODALITY_NAMES:
            arr = modalities[name]
            if arr.ndim == 1:
                scores.append(arr)
            else:
                scores.append(arr.mean(axis=1))
        return scores

    def augment(
        self,
        modalities: Dict[str, np.ndarray],
        strength: float = 0.1,
    ) -> Dict[str, np.ndarray]:
        """
        Apply data augmentation for contrastive learning.

        Parameters
        ----------
        modalities : dict[str, ndarray]
        strength : float
            Augmentation strength.

        Returns
        -------
        augmented : dict[str, ndarray]
        """
        augmented = {}
        for name, arr in modalities.items():
            noise = self.rng.randn(*arr.shape) * strength
            # Also randomly zero out some features
            if arr.ndim > 1:
                mask = self.rng.random(arr.shape) > 0.05  # 5% dropout
                augmented[name] = (arr * mask + noise).astype(np.float32)
            else:
                augmented[name] = (arr + noise).astype(np.float32)
        return augmented
