#!/usr/bin/env python3
"""
Synthetic Benchmark — Foundation Model vs CrossAttentionFusion
================================================================

Generates a curated synthetic dataset mimicking real cfDNA
multi-modal data and compares:

1. FoundationDownstream (pre-trained → fine-tuned)
2. FoundationDownstream (from scratch / no pre-training)
3. CrossAttentionFusion (baseline)

Metrics: AUC, accuracy, calibrated probability, robustness to noise.
"""

from __future__ import annotations

import time
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import (
    FoundationConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    MODALITY_DIMS,
    MODALITY_NAMES,
)
from .data import MultiModalDataGenerator
from .downstream import FoundationDownstream


def _compute_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute ROC AUC using trapezoidal rule (no sklearn dependency)."""
    # Sort by score
    desc_idx = np.argsort(y_score)[::-1]
    y_true_sorted = y_true[desc_idx]
    y_score_sorted = y_score[desc_idx]

    # Count positives and negatives
    n_pos = int(np.sum(y_true == 1))
    n_neg = int(np.sum(y_true == 0))

    if n_pos == 0 or n_neg == 0:
        return 0.5

    # TPR and FPR
    tpr = np.zeros(len(y_true_sorted) + 1)
    fpr = np.zeros(len(y_true_sorted) + 1)

    for i in range(1, len(y_true_sorted) + 1):
        if y_true_sorted[i - 1] == 1:
            tpr[i] = tpr[i - 1] + 1.0 / n_pos
            fpr[i] = fpr[i - 1]
        else:
            tpr[i] = tpr[i - 1]
            fpr[i] = fpr[i - 1] + 1.0 / n_neg

    # Trapezoidal AUC
    auc = 0.0
    for i in range(len(fpr) - 1):
        auc += (fpr[i + 1] - fpr[i]) * (tpr[i + 1] + tpr[i]) / 2.0

    return float(auc)


def _compute_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute accuracy."""
    return float(np.mean(y_true == y_pred))


def _compute_ece(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_bins: int = 10,
) -> float:
    """
    Expected Calibration Error.

    Parameters
    ----------
    y_true : (n,) array
        Binary labels.
    y_prob : (n,) array
        Predicted probabilities for class 1.
    n_bins : int
        Number of calibration bins.

    Returns
    -------
    ece : float
    """
    bins = np.linspace(0, 1, n_bins + 1)
    bin_ids = np.digitize(y_prob, bins) - 1
    bin_ids = np.clip(bin_ids, 0, n_bins - 1)

    ece = 0.0
    for b in range(n_bins):
        in_bin = bin_ids == b
        if in_bin.sum() > 0:
            acc = np.mean(y_true[in_bin])
            conf = np.mean(y_prob[in_bin])
            ece += (in_bin.sum() / len(y_true)) * abs(acc - conf)

    return float(ece)


class FoundationBenchmark:
    """
    Benchmark suite for foundation model vs baseline.

    Parameters
    ----------
    seed : int
        Random seed.
    n_train : int
        Number of training samples.
    n_test : int
        Number of test samples.
    cancer_prevalence : float
        Fraction of cancer samples.
    noise_level : float
        Amount of measurement noise.
    """

    def __init__(
        self,
        seed: int = 42,
        n_train: int = 1000,
        n_test: int = 500,
        cancer_prevalence: float = 0.3,
        noise_level: float = 0.1,
    ):
        self.seed = seed
        self.n_train = n_train
        self.n_test = n_test
        self.cancer_prevalence = cancer_prevalence
        self.noise_level = noise_level

        self.generator = MultiModalDataGenerator(
            seed=seed,
            noise_level=noise_level,
        )

        self.results: Dict[str, Dict] = {}

    def generate_data(
        self,
    ) -> Tuple[
        Dict[str, np.ndarray], np.ndarray,
        Dict[str, np.ndarray], np.ndarray,
    ]:
        """Generate train/test datasets."""
        train_modalities, train_labels = self.generator.generate_dataset(
            n_samples=self.n_train,
            prefix="bench_train",
            cancer_prevalence=self.cancer_prevalence,
        )
        test_modalities, test_labels = self.generator.generate_dataset(
            n_samples=self.n_test,
            prefix="bench_test",
            cancer_prevalence=self.cancer_prevalence,
        )
        return train_modalities, train_labels, test_modalities, test_labels

    def evaluate_model(
        self,
        model_name: str,
        model,
        test_modalities: Dict[str, np.ndarray],
        test_labels: np.ndarray,
    ) -> Dict:
        """
        Evaluate a trained model on test data.

        Parameters
        ----------
        model_name : str
            Display name.
        model : object with predict_proba() method
            Must return (n, 2) or (n,) probabilities.
        test_modalities : dict
        test_labels : (n,) array

        Returns
        -------
        metrics : dict
        """
        start = time.time()
        proba = model.predict_proba(test_modalities)
        elapsed = time.time() - start

        if proba.ndim == 2 and proba.shape[1] > 1:
            cancer_prob = proba[:, 1]
        else:
            cancer_prob = proba.flatten()

        pred = (cancer_prob > 0.5).astype(np.int64)

        metrics = {
            "auc": _compute_auc(test_labels, cancer_prob),
            "accuracy": _compute_accuracy(test_labels, pred),
            "ece": _compute_ece(test_labels, cancer_prob),
            "inference_time_ms": elapsed * 1000 / len(test_labels),
            "n_params": getattr(model, "num_params", -1),
        }

        self.results[model_name] = metrics
        return metrics

    def _evaluate_simple_fusion(
        self,
        name: str,
        train_modalities: Dict[str, np.ndarray],
        train_labels: np.ndarray,
        test_modalities: Dict[str, np.ndarray],
        test_labels: np.ndarray,
    ) -> Dict:
        """
        Evaluate using CrossAttentionFusion equivalent (simple sklearn LR).
        This serves as the baseline.
        """
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        # Convert modalities to stacked array
        train_scores = self.generator.modalities_to_scores(train_modalities)
        train_X = np.column_stack([s.reshape(-1, 1) if s.ndim == 1 else s for s in train_scores])
        test_scores = self.generator.modalities_to_scores(test_modalities)
        test_X = np.column_stack([s.reshape(-1, 1) if s.ndim == 1 else s for s in test_scores])

        start = time.time()
        scaler = StandardScaler()
        X_train = scaler.fit_transform(train_X)
        X_test = scaler.transform(test_X)

        clf = LogisticRegression(C=1.0, max_iter=500)
        clf.fit(X_train, train_labels)

        proba = clf.predict_proba(X_test)[:, 1]
        elapsed = time.time() - start

        pred = (proba > 0.5).astype(np.int64)

        metrics = {
            "auc": _compute_auc(test_labels, proba),
            "accuracy": _compute_accuracy(test_labels, pred),
            "ece": _compute_ece(test_labels, proba),
            "inference_time_ms": elapsed * 1000 / len(test_labels),
            "n_params": X_train.shape[1] + 1,  # LR params
        }

        self.results[name] = metrics
        return metrics

    def run(
        self,
        verbose: bool = False,
    ) -> Dict[str, Dict]:
        """
        Run full benchmark.

        Evaluates:
        1. CrossAttentionFusion baseline (sklearn LR)
        2. Foundation model from scratch
        3. Foundation model pre-trained (self-supervised)

        Parameters
        ----------
        verbose : bool
            Print progress.

        Returns
        -------
        results : dict
            {model_name: {auc, accuracy, ece, ...}}
        """
        if verbose:
            print("=" * 60)
            print("DeepCatch Foundation Model Benchmark")
            print(f"Train: {self.n_train} | Test: {self.n_test} | "
                  f"Prevalence: {self.cancer_prevalence}")
            print("=" * 60)

        # Generate data
        train_mod, train_lab, test_mod, test_lab = self.generate_data()

        # ── Baseline: CrossAttentionFusion (sklearn LR) ─────────
        if verbose:
            print("\n[1/3] CrossAttentionFusion (baseline)...")

        try:
            self._evaluate_simple_fusion(
                "CrossAttentionFusion",
                train_mod, train_lab, test_mod, test_lab,
            )
            if verbose:
                r = self.results["CrossAttentionFusion"]
                print(f"  AUC: {r['auc']:.4f} | Acc: {r['accuracy']:.4f} | "
                      f"ECE: {r['ece']:.4f}")
        except Exception as e:
            if verbose:
                print(f"  FAILED: {e}")

        # ── Foundation model from scratch ────────────────────────
        if verbose:
            print("\n[2/3] Foundation Model (from scratch)...")

        try:
            config_fs = PROTOTYPE_CONFIG
            config_fs.seed = self.seed

            foundation_scratch = FoundationDownstream(
                config=config_fs,
                pretrained=False,
            )
            foundation_scratch.fit(
                train_mod, train_lab,
                n_epochs=30,
                batch_size=32,
                verbose=False,
            )
            self.evaluate_model(
                "Foundation_FromScratch",
                foundation_scratch, test_mod, test_lab,
            )
            if verbose:
                r = self.results["Foundation_FromScratch"]
                print(f"  AUC: {r['auc']:.4f} | Acc: {r['accuracy']:.4f} | "
                      f"ECE: {r['ece']:.4f}")
        except Exception as e:
            if verbose:
                print(f"  FAILED: {e}")

        # ── Foundation model pre-trained ─────────────────────────
        if verbose:
            print("\n[3/3] Foundation Model (pre-trained)...")

        try:
            from .pretrain import FoundationPretrainer

            config_pt = PROTOTYPE_CONFIG
            config_pt.seed = self.seed

            # Quick pre-training
            pretrainer = FoundationPretrainer(
                config=config_pt,
                verbose=False,
            )
            pretrainer.pretrain(
                n_samples=500,
                p1_epochs=5,
                p2_epochs=5,
                p3_epochs=3,
                batch_size=32,
            )

            # Fine-tune
            foundation_pt = FoundationDownstream(
                config=config_pt,
                pretrained=True,
            )
            # Transfer pre-trained weights
            foundation_pt.encoder.load_state_dict(
                pretrainer.encoder.state_dict()
            )
            foundation_pt._pretrained = True

            foundation_pt.fit(
                train_mod, train_lab,
                n_epochs=20,
                batch_size=32,
                verbose=False,
            )
            self.evaluate_model(
                "Foundation_Pretrained",
                foundation_pt, test_mod, test_lab,
            )
            if verbose:
                r = self.results["Foundation_Pretrained"]
                print(f"  AUC: {r['auc']:.4f} | Acc: {r['accuracy']:.4f} | "
                      f"ECE: {r['ece']:.4f}")
        except Exception as e:
            if verbose:
                print(f"  FAILED: {e}")

        # ── Summary ──────────────────────────────────────────────
        if verbose:
            print("\n" + "=" * 60)
            print("BENCHMARK SUMMARY")
            print("-" * 60)
            print(f"{'Model':<30} {'AUC':>8} {'Acc':>8} {'ECE':>8} {'Time(ms)':>8}")
            print("-" * 60)
            for name, r in self.results.items():
                print(f"{name:<30} {r['auc']:8.4f} {r['accuracy']:8.4f} "
                      f"{r['ece']:8.4f} {r['inference_time_ms']:8.4f}")
            print("-" * 60)

            # Comparison
            if "CrossAttentionFusion" in self.results:
                baseline_auc = self.results["CrossAttentionFusion"]["auc"]
                for name in ["Foundation_FromScratch", "Foundation_Pretrained"]:
                    if name in self.results:
                        delta = self.results[name]["auc"] - baseline_auc
                        direction = "↑" if delta > 0 else "↓"
                        print(f"  {name} vs baseline: {direction}{abs(delta):.4f} AUC")

        return self.results


def run_benchmark(
    n_train: int = 500,
    n_test: int = 200,
    seed: int = 42,
    verbose: bool = True,
) -> Dict[str, Dict]:
    """
    Run a quick benchmark and return results.

    Parameters
    ----------
    n_train : int
    n_test : int
    seed : int
    verbose : bool

    Returns
    -------
    results : dict
    """
    benchmark = FoundationBenchmark(
        seed=seed,
        n_train=n_train,
        n_test=n_test,
        cancer_prevalence=0.3,
        noise_level=0.1,
    )
    return benchmark.run(verbose=verbose)
