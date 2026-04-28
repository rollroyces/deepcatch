#!/usr/bin/env python3
"""
MODULE 1: Nested Cross-Validation
==================================

Reference: Cawley & Talbot (2010), "On Over-fitting in Model Selection and
Subsequent Selection Bias in Performance Evaluation", JMLR 11:2079-2107.

Principle: The OUTER loop estimates generalization performance on data NEVER
touched during hyperparameter tuning. The INNER loop selects hyperparameters.
Without nesting, you leak test information into model selection — producing
optimistically biased performance estimates.

This is THE gold standard for unbiased generalization error estimation.
"""

from __future__ import annotations

import time
import warnings
from dataclasses import dataclass, field
from itertools import product
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
from scipy import stats
from sklearn.base import BaseEstimator, clone
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import ParameterGrid, StratifiedKFold

from validation_framework import BootstrapCI

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray
ParamDict = Dict[str, Sequence[Any]]
ScoreFn = Callable[[Array, Array], float]
ModelFactory = Callable[[], BaseEstimator]


# ── nested cross-validator ──────────────────────────────────────────────────


@dataclass
class OuterFoldResult:
    """Results for a single outer fold."""
    fold: int
    best_params: Dict[str, Any]
    inner_cv_mean: float
    inner_cv_std: float
    outer_score: float
    n_train: int
    n_test: int


@dataclass
class NestedCVResult:
    """Aggregated nested cross-validation results."""
    outer_folds: List[OuterFoldResult]
    mean_outer_score: float
    std_outer_score: float
    ci95_outer: Tuple[float, float]
    all_inner_scores: List[float]
    mean_inner_score: float
    std_inner_score: float
    # selection stability: how often was each param chosen?
    param_selection_counts: Dict[str, int]
    compute_time_seconds: float


class NestedCrossValidator:
    """Unbiased performance estimates via nested K-fold cross-validation.

    OUTER loop (5-fold): estimate generalization performance.
    INNER loop (3-fold): hyperparameter tuning.

    Architecture
    ───────────────────────────────────────────
      ┌───────────────────────────────────┐
      │             Data                   │
      └──────────┬────────────┬───────────┘
           Fold 1 (held out) │  …, Fold 5
                             │
                 ┌───────────┴───────────┐
                 │   Inner CV on K-1     │
                 │   folds → tune params │
                 │   → best config       │
                 └───────────────────────┘
                             │
                 Evaluate on held-out fold
    ───────────────────────────────────────────
    """

    def __init__(
        self,
        n_outer: int = 5,
        n_inner: int = 3,
        scoring: str = "roc_auc",
        random_state: int = 42,
        n_jobs: int = 1,
        verbose: bool = False,
    ):
        """Args:
            n_outer: Number of outer CV folds (default 5).
            n_inner: Number of inner CV folds (default 3).
            scoring: Metric for model selection (default 'roc_auc').
                One of: roc_auc, accuracy, f1, average_precision, sensitivity, specificity.
            random_state: Random seed for reproducibility.
            n_jobs: Parallelism (-1 = all CPUs). Set 1 for deterministic runs.
            verbose: Print per-fold progress.
        """
        self.n_outer = n_outer
        self.n_inner = n_inner
        self.scoring = scoring
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

        # ── scoring dispatch ──
        self._scorer_fn: ScoreFn = self._resolve_scorer(scoring)

    # ── public API ────────────────────────────────────────────────────────

    def validate(
        self,
        model_factory: ModelFactory,
        param_grid: ParamDict,
        X: Array,
        y: Array,
        stratify: Optional[Array] = None,
        sample_weight: Optional[Array] = None,
    ) -> NestedCVResult:
        """Run full nested cross-validation.

        Args:
            model_factory: Zero-arg callable returning a fresh estimator.
            param_grid: Dict of param_name → list_of_values.
            X: Feature matrix (n_samples, n_features) or DataFrame.
            y: Labels (n_samples,).
            stratify: Optional array for stratification (defaults to y).
            sample_weight: Optional sample weights.

        Returns:
            NestedCVResult with unbiased performance estimates.
        """
        X = np.asarray(X)
        y = np.asarray(y)

        if len(y) < self.n_outer * 2:
            raise ValueError(
                f"Need at least {self.n_outer * 2} samples for "
                f"{self.n_outer}-fold nested CV, got {len(y)}"
            )

        if stratify is None:
            stratify = y

        # Pre-compute parameter combinations
        param_list = list(ParameterGrid(param_grid))
        if self.verbose:
            print(
                f"Nested CV: {self.n_outer} outer × {self.n_inner} inner "
                f"folds, {len(param_list)} param combos, {len(np.unique(y))} classes"
            )

        t0 = time.time()

        outer_skf = StratifiedKFold(
            n_splits=self.n_outer,
            shuffle=True,
            random_state=self.random_state,
        )

        outer_fold_results: List[OuterFoldResult] = []
        all_inner_best_scores: List[float] = []
        param_selection_counts: Dict[str, int] = {}

        for fold_idx, (train_idx, test_idx) in enumerate(
            outer_skf.split(X, stratify)
        ):
            if self.verbose:
                print(f"  Outer fold {fold_idx + 1}/{self.n_outer} ", end="")

            X_train, y_train = X[train_idx], y[train_idx]
            X_test, y_test = X[test_idx], y[test_idx]

            # Inner CV: tune on train data only
            inner_skf = StratifiedKFold(
                n_splits=self.n_inner,
                shuffle=True,
                random_state=self.random_state + fold_idx * 10,
            )

            best_score = -np.inf
            best_params: Dict[str, Any] = {}
            best_inner_std = 0.0
            inner_scores_for_best: List[float] = []

            for params in param_list:
                fold_scores: List[float] = []

                for inner_train_idx, inner_val_idx in inner_skf.split(
                    X_train, y_train
                ):
                    model = model_factory()
                    model.set_params(**params)
                    model.fit(
                        X_train[inner_train_idx],
                        y_train[inner_train_idx],
                    )

                    y_val = y_train[inner_val_idx]
                    X_val = X_train[inner_val_idx]

                    if hasattr(model, "predict_proba"):
                        y_score = model.predict_proba(X_val)[:, 1]
                    else:
                        y_score = model.decision_function(X_val)

                    y_pred = model.predict(X_val)
                    score = self._scorer_fn(y_val, y_pred, y_score)
                    fold_scores.append(score)

                mean_score = float(np.mean(fold_scores))
                if mean_score > best_score:
                    best_score = mean_score
                    best_params = dict(params)
                    best_inner_std = float(np.std(fold_scores))
                    inner_scores_for_best = fold_scores

            # Re-train best model on ALL inner data, evaluate on outer test
            model = model_factory()
            model.set_params(**best_params)
            model.fit(X_train, y_train)

            if hasattr(model, "predict_proba"):
                y_score_test = model.predict_proba(X_test)[:, 1]
            else:
                y_score_test = model.decision_function(X_test)

            y_pred_test = model.predict(X_test)
            outer_score = self._scorer_fn(y_test, y_pred_test, y_score_test)

            # Track param selection
            param_key = str(sorted(best_params.items()))
            param_selection_counts[param_key] = (
                param_selection_counts.get(param_key, 0) + 1
            )

            outer_fold_results.append(
                OuterFoldResult(
                    fold=fold_idx + 1,
                    best_params=best_params,
                    inner_cv_mean=best_score,
                    inner_cv_std=best_inner_std,
                    outer_score=outer_score,
                    n_train=len(y_train),
                    n_test=len(y_test),
                )
            )
            all_inner_best_scores.append(best_score)

            if self.verbose:
                print(
                    f"inner={best_score:.4f}±{best_inner_std:.4f}  "
                    f"outer={outer_score:.4f}"
                )

        outer_scores = np.array([f.outer_score for f in outer_fold_results])
        mean_outer = float(np.mean(outer_scores))
        std_outer = float(np.std(outer_scores, ddof=1))

        # CI via t-distribution on outer fold scores
        if len(outer_scores) > 1:
            se = std_outer / np.sqrt(len(outer_scores))
            t_crit = float(stats.t.ppf(0.975, df=len(outer_scores) - 1))
            ci95 = (mean_outer - t_crit * se, mean_outer + t_crit * se)
        else:
            ci95 = (mean_outer, mean_outer)

        elapsed = time.time() - t0

        return NestedCVResult(
            outer_folds=outer_fold_results,
            mean_outer_score=mean_outer,
            std_outer_score=std_outer,
            ci95_outer=ci95,
            all_inner_scores=all_inner_best_scores,
            mean_inner_score=float(np.mean(all_inner_best_scores)),
            std_inner_score=float(np.std(all_inner_best_scores)),
            param_selection_counts=param_selection_counts,
            compute_time_seconds=elapsed,
        )

    def compare_models(
        self,
        model_factories: Dict[str, ModelFactory],
        param_grids: Dict[str, ParamDict],
        X: Array,
        y: Array,
    ) -> Dict[str, NestedCVResult]:
        """Run nested CV for multiple models for head-to-head comparison.

        Args:
            model_factories: Dict of name → model factory.
            param_grids: Dict of name → param grid.
            X, y: Data.

        Returns:
            Dict of name → NestedCVResult.
        """
        results: Dict[str, NestedCVResult] = {}
        for name in model_factories:
            if self.verbose:
                print(f"\n── {name} ──")
            results[name] = self.validate(
                model_factories[name],
                param_grids.get(name, {}),
                X, y,
            )
        return results

    # ── scoring ───────────────────────────────────────────────────────────

    @staticmethod
    def _resolve_scorer(scoring: str) -> ScoreFn:
        """Resolve a scoring string to a callable that takes (y_true, y_pred, y_score)."""
        scorers: Dict[str, ScoreFn] = {
            "roc_auc": lambda yt, yp, ys: (
                roc_auc_score(yt, ys) if ys is not None else float("nan")
            ),
            "accuracy": lambda yt, yp, ys: accuracy_score(yt, yp),
            "f1": lambda yt, yp, ys: f1_score(yt, yp, average="binary", zero_division=0),
            "average_precision": lambda yt, yp, ys: (
                average_precision_score(yt, ys) if ys is not None else float("nan")
            ),
            "sensitivity": lambda yt, yp, ys: recall_score(yt, yp, average="binary", zero_division=0),
            "specificity": lambda yt, yp, ys: (
                float(np.sum((yt == 0) & (yp == 0))) / float(np.sum(yt == 0))
                if np.sum(yt == 0) > 0 else float("nan")
            ),
            "precision": lambda yt, yp, ys: precision_score(yt, yp, average="binary", zero_division=0),
        }

        scoring_lower = scoring.lower()
        if scoring_lower not in scorers:
            raise ValueError(
                f"Unknown scoring '{scoring}'. Available: {list(scorers.keys())}"
            )
        return scorers[scoring_lower]

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: NestedCVResult, model_name: str = "Model") -> str:
        """Format a nested CV result as a publication-ready string."""
        lines = [
            f"══ Nested Cross-Validation: {model_name} ══",
            f"  Outer CV ({len(result.outer_folds)}-fold):",
            f"    Mean: {result.mean_outer_score:.4f}",
            f"    Std:  {result.std_outer_score:.4f}",
            f"    CI95: [{result.ci95_outer[0]:.4f}, {result.ci95_outer[1]:.4f}]",
            f"  Inner CV (best per outer fold):",
            f"    Mean: {result.mean_inner_score:.4f}",
            f"    Std:  {result.std_inner_score:.4f}",
            f"  Optimism gap (inner - outer): "
            f"{result.mean_inner_score - result.mean_outer_score:.4f}",
            f"  Parameter stability:",
        ]
        for param_key, count in sorted(
            result.param_selection_counts.items(),
            key=lambda x: -x[1],
        ):
            lines.append(f"    {count}/{len(result.outer_folds)}: {param_key}")
        lines.append(f"  Time: {result.compute_time_seconds:.1f}s")
        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 1: NestedCrossValidator — self-test")
    print("=" * 60)

    # Synthetic dataset: 500 samples, 10 features, ~30% positive class
    rng = np.random.RandomState(42)
    n = 500
    X_test = rng.randn(n, 10)
    # Add some signal in feature 0 for positive class
    y_test = np.zeros(n, dtype=int)
    pos_idx = rng.choice(n, size=int(n * 0.3), replace=False)
    X_test[pos_idx, 0] += 1.5
    y_test[pos_idx] = 1

    from sklearn.linear_model import LogisticRegression

    ncv = NestedCrossValidator(
        n_outer=5, n_inner=3, scoring="roc_auc", verbose=True
    )

    param_grid = {"C": [0.01, 0.1, 1.0, 10.0]}

    result = ncv.validate(
        model_factory=lambda: LogisticRegression(solver="liblinear", max_iter=500),
        param_grid=param_grid,
        X=X_test,
        y=y_test,
    )

    print(ncv.report(result))

    # Check that outer score is not optimistically biased relative to inner
    optimism_gap = result.mean_inner_score - result.mean_outer_score
    if optimism_gap > 0.05:
        print(f"\n⚠ Optimism gap {optimism_gap:.4f} — nested CV important!")
    else:
        print(f"\n✓ Optimism gap {optimism_gap:.4f} — small, tuning not overfitting")

    print("\nSelf-test complete.")
    sys.exit(0)
