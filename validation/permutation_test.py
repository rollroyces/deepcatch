#!/usr/bin/env python3
"""
MODULE 2: Permutation Testing
==============================

Reference: Ojala & Garriga (2010), "Permutation Tests for Studying Classifier
Performance", JMLR 11:1833-1863.

Tests the null hypothesis: "The model is fitting noise, not signal."

For each model:
  1. Train normally, record reference score.
  2. Shuffle labels, re-train, re-evaluate (repeat n_permutations times).
  3. p-value = fraction of permuted scores ≥ reference score.

The p-value answers: "If there's no signal in the data, how often would I see
a score this good by chance?" If p < 0.05 (corrected), reject the null.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from sklearn.base import BaseEstimator, clone
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    recall_score,
    roc_auc_score,
)

from validation_framework import BootstrapCI, SignificanceTester

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray
ModelFactory = Callable[[], BaseEstimator]


# ── permutation test result ─────────────────────────────────────────────────


@dataclass
class PermutationResult:
    """Result of a permutation test for one model."""
    model_name: str
    reference_score: float
    permuted_scores: List[float]
    p_value: float
    p_value_corrected: float  # Bonferroni if multi-model
    mean_permuted: float
    std_permuted: float
    ci95_permuted: Tuple[float, float]
    z_score: float  # (reference - mean_permuted) / std_permuted
    n_permutations: int
    significant: bool  # p < 0.05
    compute_time_seconds: float


class PermutationTester:
    """Full permutation testing framework for classifier validation.

    Usage:
        tester = PermutationTester(n_permutations=1000)
        result = tester.test(model_factory, X, y)

    For multiple models with correction:
        results = tester.test_multiple(
            {"Logistic": lr_fn, "RandomForest": rf_fn}, X, y
        )
    """

    def __init__(
        self,
        n_permutations: int = 1000,
        scoring: str = "roc_auc",
        random_state: int = 42,
        n_jobs: int = 1,
        verbose: bool = False,
    ):
        """Args:
            n_permutations: Number of label shuffles (≥1000 for stable p-values).
            scoring: Metric to evaluate. One of: roc_auc, accuracy, f1, sensitivity.
            random_state: Base seed for reproducibility.
            n_jobs: Parallel workers (-1 = all CPUs).
            verbose: Print progress.
        """
        self.n_permutations = n_permutations
        self.scoring = scoring.lower()
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

        self._scorer = self._resolve_scorer()
        self._supported_scorers = {"roc_auc", "accuracy", "f1", "sensitivity"}

    # ── public API ────────────────────────────────────────────────────────

    def test(
        self,
        model_factory: ModelFactory,
        X: Array,
        y: Array,
        model_name: str = "Model",
    ) -> PermutationResult:
        """Run permutation test for a single model.

        Args:
            model_factory: Zero-arg callable → fresh estimator.
            X: Feature matrix.
            y: Labels (binary).
            model_name: Label for reporting.

        Returns:
            PermutationResult with p-value and diagnostics.
        """
        X = np.asarray(X)
        y = np.asarray(y).ravel().astype(int)

        n = len(y)
        if n < 20:
            raise ValueError(f"Need ≥ 20 samples for permutation test; got {n}")

        if len(np.unique(y)) < 2:
            raise ValueError("Need both classes for permutation test")

        rng = np.random.RandomState(self.random_state)

        # ── Reference score ──
        reference_score = self._compute_score(model_factory, X, y)
        if self.verbose:
            print(f"  Reference {self.scoring}: {reference_score:.4f}")

        # ── Permutation loop ──
        t0 = time.time()
        permuted_scores: List[float] = []

        for i in range(self.n_permutations):
            y_shuffled = y[rng.permutation(n)]
            perm_score = self._compute_score(model_factory, X, y_shuffled)
            permuted_scores.append(perm_score)

            if self.verbose and (i + 1) % max(1, self.n_permutations // 10) == 0:
                elapsed = time.time() - t0
                print(
                    f"  Permutation {i + 1}/{self.n_permutations} "
                    f"({elapsed:.1f}s)"
                )

        elapsed = time.time() - t0

        perms = np.array(permuted_scores)
        mean_perm = float(np.mean(perms))
        std_perm = float(np.std(perms, ddof=1)) if len(perms) > 1 else 0.0

        # p-value: fraction of permuted scores ≥ reference (or ≤, depending on metric)
        # For all our metrics, higher is better
        p_value = float(np.mean(perms >= reference_score))

        # Handle edge case: reference better than all permutations
        if p_value == 0.0:
            p_value = 1.0 / (self.n_permutations + 1)  # continuity correction

        z_score = (
            (reference_score - mean_perm) / std_perm
            if std_perm > 0
            else float("inf")
        )

        ci95 = (
            float(np.percentile(perms, 2.5)),
            float(np.percentile(perms, 97.5)),
        )

        return PermutationResult(
            model_name=model_name,
            reference_score=reference_score,
            permuted_scores=list(perms),
            p_value=p_value,
            p_value_corrected=p_value,  # will be overwritten in test_multiple
            mean_permuted=mean_perm,
            std_permuted=std_perm,
            ci95_permuted=ci95,
            z_score=z_score,
            n_permutations=self.n_permutations,
            significant=p_value < 0.05,
            compute_time_seconds=elapsed,
        )

    def test_multiple(
        self,
        model_factories: Dict[str, ModelFactory],
        X: Array,
        y: Array,
        correction: str = "bonferroni",
    ) -> Dict[str, PermutationResult]:
        """Permutation test for multiple models with multiple-testing correction.

        Args:
            model_factories: Dict of name → model factory.
            X, y: Data.
            correction: 'bonferroni' or 'bh' (Benjamini-Hochberg).

        Returns:
            Dict of name → PermutationResult with corrected p-values.
        """
        results: Dict[str, PermutationResult] = {}
        raw_p_values: Dict[str, float] = {}

        for name, factory in model_factories.items():
            res = self.test(factory, X, y, model_name=name)
            results[name] = res
            raw_p_values[name] = res.p_value

        # Apply correction
        names = list(raw_p_values.keys())
        p_vals = np.array([raw_p_values[n] for n in names])

        if correction == "bonferroni":
            corrected = SignificanceTester.bonferroni_correct(p_vals)
        elif correction == "bh":
            corrected, rejected = SignificanceTester.benjamini_hochberg(p_vals)
        else:
            raise ValueError(f"Unknown correction '{correction}'")

        for i, name in enumerate(names):
            results[name].p_value_corrected = float(corrected[i])
            results[name].significant = float(corrected[i]) < 0.05

        return results

    # ── scoring helpers ───────────────────────────────────────────────────

    def _compute_score(
        self, model_factory: ModelFactory, X: Array, y: Array
    ) -> float:
        """Train model and compute metric with internal hold-out split."""
        n = len(y)
        # Simple 80/20 single-split evaluation for speed
        rng = np.random.RandomState(self.random_state)
        idx = rng.permutation(n)
        n_train = int(0.8 * n)
        train_idx, test_idx = idx[:n_train], idx[n_train:]

        if len(np.unique(y[train_idx])) < 2 or len(np.unique(y[test_idx])) < 2:
            # Not enough class diversity — use full data (acceptable for permutation)
            return self._scorer(y, y)  # degenerate

        model = model_factory()
        model.fit(X[train_idx], y[train_idx])

        y_pred = model.predict(X[test_idx])

        if hasattr(model, "predict_proba"):
            y_score = model.predict_proba(X[test_idx])[:, 1]
        elif hasattr(model, "decision_function"):
            y_score = model.decision_function(X[test_idx])
        else:
            y_score = None

        return self._scorer(y[test_idx], y_pred, y_score)

    def _resolve_scorer(self):
        """Return a scorer fn: (y_true, y_pred, y_score) → float."""
        scorers = {
            "roc_auc": lambda yt, yp, ys: (
                roc_auc_score(yt, ys)
                if ys is not None and len(np.unique(yt)) > 1
                else float("nan")
            ),
            "accuracy": lambda yt, yp, ys: accuracy_score(yt, yp),
            "f1": lambda yt, yp, ys: f1_score(
                yt, yp, average="binary", zero_division=0
            ),
            "sensitivity": lambda yt, yp, ys: recall_score(
                yt, yp, average="binary", zero_division=0
            ),
        }
        if self.scoring not in scorers:
            raise ValueError(
                f"Unknown scoring '{self.scoring}'. "
                f"Available: {list(scorers.keys())}"
            )
        return scorers[self.scoring]

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: PermutationResult) -> str:
        """Publication-ready string report."""
        sig = "SIGNIFICANT ✓" if result.significant else "NOT SIGNIFICANT ✗"
        return (
            f"══ Permutation Test: {result.model_name} ══\n"
            f"  Reference {result.n_permutations} permutations\n"
            f"  Reference score:            {result.reference_score:.4f}\n"
            f"  Permuted mean ± std:        {result.mean_permuted:.4f} ± {result.std_permuted:.4f}\n"
            f"  Permuted CI95:              [{result.ci95_permuted[0]:.4f}, {result.ci95_permuted[1]:.4f}]\n"
            f"  z-score:                    {result.z_score:.2f}\n"
            f"  p-value (raw):              {result.p_value:.6f}\n"
            f"  p-value (corrected):        {result.p_value_corrected:.6f}\n"
            f"  Verdict:                    {sig}\n"
            f"  Time:                       {result.compute_time_seconds:.1f}s"
        )


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 2: PermutationTester — self-test")
    print("=" * 60)

    # Synthetic data with real signal
    rng = np.random.RandomState(42)
    n = 300
    X_test = rng.randn(n, 8)
    y_test = np.zeros(n, dtype=int)
    pos_idx = rng.choice(n, size=n // 3, replace=False)
    X_test[pos_idx, 0] += 2.0  # strong signal
    y_test[pos_idx] = 1

    from sklearn.linear_model import LogisticRegression

    pt = PermutationTester(n_permutations=200, scoring="roc_auc", verbose=True)

    result = pt.test(
        model_factory=lambda: LogisticRegression(solver="liblinear"),
        X=X_test,
        y=y_test,
        model_name="LogisticRegression",
    )

    print(pt.report(result))

    # Also test on pure noise to confirm non-significant
    print("\n── Negative control (pure noise) ──")
    y_noise = rng.randint(0, 2, size=n)  # random labels
    noise_result = pt.test(
        model_factory=lambda: LogisticRegression(solver="liblinear"),
        X=X_test,
        y=y_noise,
        model_name="LogisticRegression (noise)",
    )
    print(pt.report(noise_result))

    if not noise_result.significant:
        print("\n✓ Correctly identified pure noise as non-significant")
    else:
        print("\n⚠ Pure noise flagged as significant — check data generation")

    print("\nSelf-test complete.")
    sys.exit(0)
