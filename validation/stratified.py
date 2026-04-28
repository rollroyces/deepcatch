#!/usr/bin/env python3
"""
MODULE 6: Stratified Performance Analysis
==========================================

Reference: Kent et al. (2010), "Assessing and Reporting Heterogeneity in
Treatment Effects in Clinical Trials: A Proposal", Trials 11:85.

Performance broken down by clinically relevant strata:
  - Cancer type (LUAD, COADREAD, BRCA, PRAD, LGG, etc.)
  - Cancer stage (I-IV if available)
  - ctDNA fraction level (0.001%-0.01%, 0.01%-0.1%, 0.1%+)
  - Patient age group (18-40, 41-60, 61-80)
  - Fragment size profile cluster

Reports per-stratum AUC, sensitivity, specificity with CIs.
Tests for interaction effects (does performance differ significantly
between strata?).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)

from validation_framework import BootstrapCI

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray


# ── stratum result ──────────────────────────────────────────────────────────


@dataclass
class StratumResult:
    """Performance metrics for one stratum."""
    stratum_name: str
    n_samples: int
    n_pos: int
    n_neg: int
    prevalence: float
    auc: float
    auc_ci95: Tuple[float, float]
    sensitivity: float
    sensitivity_ci95: Tuple[float, float]
    specificity: float
    specificity_ci95: Tuple[float, float]
    f1: float
    f1_ci95: Tuple[float, float]


@dataclass
class InteractionTest:
    """Result of testing whether performance differs between strata."""
    metric: str
    strata_pair: Tuple[str, str]
    delta: float
    p_value: float
    p_value_corrected: float
    significant: bool


@dataclass
class StratifiedResult:
    """Full stratified performance analysis."""
    strata: List[StratumResult]
    overall: StratumResult
    interaction_tests: List[InteractionTest]
    n_strata: int
    total_samples: int


class StratifiedAnalyzer:
    """Performance analysis stratified by clinically relevant subgroups.

    Usage:
        analyzer = StratifiedAnalyzer()
        result = analyzer.analyze(
            y_true, y_score,
            strata_labels,  # array of stratum assignments
            threshold=0.5,
        )
    """

    def __init__(
        self,
        min_stratum_size: int = 10,
        n_bootstrap: int = 1000,
        ci: float = 0.95,
        random_state: int = 42,
    ):
        """Args:
            min_stratum_size: Minimum samples per stratum (skipped if fewer).
            n_bootstrap: Bootstrap replicates for CIs.
            ci: Confidence level.
            random_state: Seed.
        """
        self.min_stratum_size = min_stratum_size
        self.n_bootstrap = n_bootstrap
        self.ci = ci
        self.random_state = random_state

    # ── public API ────────────────────────────────────────────────────────

    def analyze(
        self,
        y_true: Array,
        y_score: Array,
        strata_labels: Array,
        threshold: float = 0.5,
        stratum_names: Optional[Dict[Any, str]] = None,
    ) -> StratifiedResult:
        """Full stratified performance analysis.

        Args:
            y_true: Ground-truth labels, shape (n_samples,).
            y_score: Predicted probabilities in [0, 1], shape (n_samples,).
            strata_labels: Stratum assignment per sample, shape (n_samples,).
                Can be int, string, or categorical.
            threshold: Classification threshold for binary metrics.
            stratum_names: Optional dict mapping stratum value → display name.

        Returns:
            StratifiedResult with per-stratum metrics and interaction tests.
        """
        y_true = np.asarray(y_true).ravel().astype(int)
        y_score = np.asarray(y_score).ravel().astype(float)
        strata_labels = np.asarray(strata_labels).ravel()

        self._validate(y_true, y_score, strata_labels)

        unique_strata = np.unique(strata_labels)
        strata_results: List[StratumResult] = []

        for stratum_val in unique_strata:
            mask = strata_labels == stratum_val
            n = int(np.sum(mask))

            if n < self.min_stratum_size:
                continue

            name = (
                stratum_names.get(stratum_val, str(stratum_val))
                if stratum_names
                else str(stratum_val)
            )

            yt = y_true[mask]
            ys = y_score[mask]

            if len(np.unique(yt)) < 2:
                # Single-class stratum — report N/A for AUC
                strata_results.append(StratumResult(
                    stratum_name=name,
                    n_samples=n,
                    n_pos=int(np.sum(yt == 1)),
                    n_neg=int(np.sum(yt == 0)),
                    prevalence=float(np.mean(yt)),
                    auc=float("nan"),
                    auc_ci95=(float("nan"), float("nan")),
                    sensitivity=float("nan"),
                    sensitivity_ci95=(float("nan"), float("nan")),
                    specificity=float("nan"),
                    specificity_ci95=(float("nan"), float("nan")),
                    f1=float("nan"),
                    f1_ci95=(float("nan"), float("nan")),
                ))
                continue

            y_pred = (ys >= threshold).astype(int)

            # Bootstrap CIs
            ci_record = self._bootstrap_metrics(yt, ys, n_bootstrap=self.n_bootstrap)
            sr = StratumResult(
                stratum_name=name,
                n_samples=n,
                n_pos=int(np.sum(yt == 1)),
                n_neg=int(np.sum(yt == 0)),
                prevalence=float(np.mean(yt)),
                auc=ci_record["auc"]["point"],
                auc_ci95=ci_record["auc"]["ci"],
                sensitivity=ci_record["sensitivity"]["point"],
                sensitivity_ci95=ci_record["sensitivity"]["ci"],
                specificity=ci_record["specificity"]["point"],
                specificity_ci95=ci_record["specificity"]["ci"],
                f1=ci_record["f1"]["point"],
                f1_ci95=ci_record["f1"]["ci"],
            )
            strata_results.append(sr)

        # Overall
        y_pred_all = (y_score >= threshold).astype(int)
        ci_all = self._bootstrap_metrics(y_true, y_score, n_bootstrap=self.n_bootstrap)
        overall = StratumResult(
            stratum_name="Overall",
            n_samples=len(y_true),
            n_pos=int(np.sum(y_true == 1)),
            n_neg=int(np.sum(y_true == 0)),
            prevalence=float(np.mean(y_true)),
            auc=ci_all["auc"]["point"],
            auc_ci95=ci_all["auc"]["ci"],
            sensitivity=ci_all["sensitivity"]["point"],
            sensitivity_ci95=ci_all["sensitivity"]["ci"],
            specificity=ci_all["specificity"]["point"],
            specificity_ci95=ci_all["specificity"]["ci"],
            f1=ci_all["f1"]["point"],
            f1_ci95=ci_all["f1"]["ci"],
        )

        # Interaction tests: does AUC differ between strata?
        interaction_tests = self._test_interactions(
            y_true, y_score, strata_labels, strata_results
        )

        return StratifiedResult(
            strata=strata_results,
            overall=overall,
            interaction_tests=interaction_tests,
            n_strata=len(strata_results),
            total_samples=len(y_true),
        )

    def analyze_by_cancer_type(
        self,
        y_true: Array,
        y_score: Array,
        cancer_types: Array,
        threshold: float = 0.5,
    ) -> StratifiedResult:
        """Convenience: stratify by cancer type."""
        return self.analyze(y_true, y_score, cancer_types, threshold)

    def analyze_by_ctdna_fraction(
        self,
        y_true: Array,
        y_score: Array,
        ctdna_fractions: Array,
        threshold: float = 0.5,
        bins: Optional[List[float]] = None,
    ) -> StratifiedResult:
        """Stratify by ctDNA fraction, auto-binning continuous values.

        Args:
            ctdna_fractions: Continuous ctDNA fraction values.
            bins: Optional list of bin edges. Default: [0, 0.0001, 0.001, 0.01, inf].
        """
        if bins is None:
            bins = [0, 0.0001, 0.001, 0.01, np.inf]

        ctdna = np.asarray(ctdna_fractions).ravel()
        strata: list[str] = []
        stratum_names: Dict[str, str] = {}

        for i in range(len(bins) - 1):
            lo, hi = bins[i], bins[i + 1]
            label = f"ctDNA [{lo:.4f}%-{hi:.4f}%)"
            stratum_names[label] = label

        # Categorize
        categories = np.empty(len(ctdna), dtype=object)
        for i in range(len(bins) - 1):
            lo, hi = bins[i], bins[i + 1]
            mask = (ctdna >= lo) & (ctdna < hi)
            label = f"ctDNA [{lo:.4f}%-{hi:.4f}%)"
            categories[mask] = label

        if np.sum(np.isnan(categories)) > 0:
            categories[np.isnan(categories)] = "ctDNA [unknown]"

        return self.analyze(y_true, y_score, categories, threshold, stratum_names)

    # ── interaction testing ───────────────────────────────────────────────

    def _test_interactions(
        self,
        y_true: Array,
        y_score: Array,
        strata_labels: Array,
        strata_results: List[StratumResult],
    ) -> List[InteractionTest]:
        """Test whether AUC differs significantly between strata pairs.

        Uses DeLong-style bootstrap test for paired strata.
        """
        if len(strata_results) < 2:
            return []

        tests: List[InteractionTest] = []
        p_values: List[float] = []

        # Get valid strata (with meaningful AUC)
        valid_strata = [
            s for s in strata_results
            if not np.isnan(s.auc) and s.n_samples >= self.min_stratum_size
        ]

        for i in range(len(valid_strata)):
            for j in range(i + 1, len(valid_strata)):
                sa, sb = valid_strata[i], valid_strata[j]

                # Bootstrap test of difference
                mask_a = strata_labels == sa.stratum_name
                mask_b = strata_labels == sb.stratum_name

                yt_a, ys_a = y_true[mask_a], y_score[mask_a]
                yt_b, ys_b = y_true[mask_b], y_score[mask_b]

                # Compare AUC using bootstrap
                auc_diff_stats = self._bootstrap_auc_diff(
                    yt_a, ys_a, yt_b, ys_b
                )

                p_values.append(auc_diff_stats["p_value"])
                tests.append(InteractionTest(
                    metric="AUC",
                    strata_pair=(sa.stratum_name, sb.stratum_name),
                    delta=auc_diff_stats["delta_mean"],
                    p_value=auc_diff_stats["p_value"],
                    p_value_corrected=auc_diff_stats["p_value"],  # updated below
                    significant=auc_diff_stats["p_value"] < 0.05,
                ))

        # Bonferroni correction
        if p_values:
            from validation_framework import SignificanceTester
            corrected = SignificanceTester.bonferroni_correct(np.array(p_values))

            for idx, test in enumerate(tests):
                test.p_value_corrected = float(corrected[idx])
                test.significant = float(corrected[idx]) < 0.05

        return tests

    def _bootstrap_auc_diff(
        self,
        y_true_a: Array,
        y_score_a: Array,
        y_true_b: Array,
        y_score_b: Array,
        n_bootstrap: int = 1000,
    ) -> Dict[str, float]:
        """Bootstrap test for AUC difference between two independent strata."""
        rng = np.random.RandomState(self.random_state)

        n_a = len(y_true_a)
        n_b = len(y_true_b)

        diffs = np.zeros(n_bootstrap)

        for k in range(n_bootstrap):
            idx_a = rng.choice(n_a, size=n_a, replace=True)
            idx_b = rng.choice(n_b, size=n_b, replace=True)

            try:
                auc_a = roc_auc_score(y_true_a[idx_a], y_score_a[idx_a])
                auc_b = roc_auc_score(y_true_b[idx_b], y_score_b[idx_b])
                diffs[k] = auc_b - auc_a
            except ValueError:
                diffs[k] = 0.0

        mean_diff = float(np.mean(diffs))
        p_val = float(np.mean(np.abs(diffs - mean_diff) >= np.abs(mean_diff))) \
            if abs(mean_diff) > 1e-12 else 1.0

        return {
            "delta_mean": mean_diff,
            "delta_std": float(np.std(diffs)),
            "p_value": p_val,
            "ci_lower": float(np.percentile(diffs, 2.5)),
            "ci_upper": float(np.percentile(diffs, 97.5)),
        }

    # ── bootstrap metrics ─────────────────────────────────────────────────

    def _bootstrap_metrics(
        self,
        y_true: Array,
        y_score: Array,
        threshold: float = 0.5,
        n_bootstrap: int = 1000,
    ) -> Dict[str, Dict]:
        """Compute bootstrap CIs for all key metrics."""

        y_true = np.asarray(y_true).ravel().astype(int)
        y_score = np.asarray(y_score).ravel().astype(float)

        n = len(y_true)
        rng = np.random.RandomState(self.random_state)

        aucs = np.zeros(n_bootstrap)
        sens = np.zeros(n_bootstrap)
        specs = np.zeros(n_bootstrap)
        f1s = np.zeros(n_bootstrap)

        for k in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            yt = y_true[idx]
            ys = y_score[idx]
            yp = (ys >= threshold).astype(int)

            try:
                aucs[k] = roc_auc_score(yt, ys) if len(np.unique(yt)) > 1 else float("nan")
            except ValueError:
                aucs[k] = float("nan")

            cm = confusion_matrix(yt, yp, labels=[0, 1])
            if cm.shape == (2, 2):
                tn, fp, fn, tp = cm.ravel()
                sens[k] = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                specs[k] = tn / (tn + fp) if (tn + fp) > 0 else 0.0
                f1s[k] = (2 * tp) / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0
            else:
                sens[k] = specs[k] = f1s[k] = float("nan")

        alpha_2 = (1 - self.ci) / 2

        def _ci(arr):
            valid = arr[~np.isnan(arr)]
            if len(valid) < 2:
                return (float("nan"), float("nan"))
            return (
                float(np.percentile(valid, 100 * alpha_2)),
                float(np.percentile(valid, 100 * (1 - alpha_2))),
            )

        return {
            "auc": {"ci": _ci(aucs), "point": float(np.nanmean(aucs))},
            "sensitivity": {"ci": _ci(sens), "point": float(np.nanmean(sens))},
            "specificity": {"ci": _ci(specs), "point": float(np.nanmean(specs))},
            "f1": {"ci": _ci(f1s), "point": float(np.nanmean(f1s))},
        }

    # ── validation ────────────────────────────────────────────────────────

    @staticmethod
    def _validate(y_true: Array, y_score: Array, strata: Array) -> None:
        if len(y_true) != len(strata) or len(y_score) != len(strata):
            raise ValueError(
                f"Length mismatch: y_true={len(y_true)}, "
                f"y_score={len(y_score)}, strata={len(strata)}"
            )
        if len(np.unique(y_true)) < 2:
            raise ValueError("y_true must contain both classes")
        if len(np.unique(strata)) < 2:
            raise ValueError("Need ≥ 2 strata for stratified analysis")

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: StratifiedResult) -> str:
        """Publication-ready stratified analysis report."""
        lines = [
            "══ Stratified Performance Analysis ══",
            f"  Total samples: {result.total_samples}",
            f"  Number of strata: {result.n_strata}",
            f"  Minimum stratum size: {min(s.n_samples for s in result.strata)}",
            "",
            f"  {'Stratum':>30s} {'N':>6s} {'Prev':>7s} {'AUC':>7s} {'[CI95]':>20s} {'Sens':>7s} {'Spec':>7s}",
            f"  {'─'*30} {'─'*6} {'─'*7} {'─'*7} {'─'*20} {'─'*7} {'─'*7}",
        ]

        # Overall first
        o = result.overall
        lines.append(
            f"  {'OVERALL':>30s} {o.n_samples:6d} {o.prevalence:7.3f} "
            f"{o.auc:7.4f} [{o.auc_ci95[0]:.4f}, {o.auc_ci95[1]:.4f}] "
            f"{o.sensitivity:7.4f} {o.specificity:7.4f}"
        )
        lines.append("")

        for s in result.strata:
            auc_str = f"{s.auc:.4f}" if not np.isnan(s.auc) else "   N/A"
            ci_str = (
                f"[{s.auc_ci95[0]:.4f}, {s.auc_ci95[1]:.4f}]"
                if not np.isnan(s.auc_ci95[0])
                else "       N/A"
            )
            sens_str = f"{s.sensitivity:.4f}" if not np.isnan(s.sensitivity) else "   N/A"
            spec_str = f"{s.specificity:.4f}" if not np.isnan(s.specificity) else "   N/A"

            lines.append(
                f"  {s.stratum_name:>30s} {s.n_samples:6d} {s.prevalence:7.3f} "
                f"{auc_str:>7s} {ci_str:>20s} {sens_str:>7s} {spec_str:>7s}"
            )

        # Interaction tests
        if result.interaction_tests:
            lines.append(f"\n  Interaction Tests (corrected for multiple comparisons):")
            lines.append(f"  {'Stratum A':>25s} vs {'Stratum B':<25s} {'Δ AUC':>7s} {'p':>7s} {'Verdict':>10s}")

            sig_count = 0
            for t in result.interaction_tests:
                sig = "SIG" if t.significant else "ns"
                if t.significant:
                    sig_count += 1
                a_name = t.strata_pair[0][:25]
                b_name = t.strata_pair[1][:25]
                lines.append(
                    f"  {a_name:>25s} vs {b_name:<25s} "
                    f"{t.delta:+7.4f} {t.p_value_corrected:7.4f} {sig:>10s}"
                )

            if sig_count > 0:
                lines.append(
                    f"\n  {sig_count}/{len(result.interaction_tests)} "
                    f"interactions significant — performance varies by stratum."
                )
            else:
                lines.append(
                    f"\n  No significant interactions — performance consistent across strata."
                )

        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 6: StratifiedAnalyzer — self-test")
    print("=" * 60)

    rng = np.random.RandomState(42)
    n = 600

    # Three cancer types
    cancer_types = np.array(
        ["LUAD"] * 200 + ["COADREAD"] * 200 + ["BRCA"] * 200
    )

    # True labels with type-specific prevalence
    y_true = np.zeros(n, dtype=int)
    for i, ct in enumerate(cancer_types):
        if ct == "LUAD":
            p = 0.40
        elif ct == "COADREAD":
            p = 0.25
        else:
            p = 0.35
        y_true[i] = int(rng.rand() < p)

    # Scores: LUAD easier to detect, COADREAD harder
    y_score = np.zeros(n)
    for i in range(n):
        base = 0.3 + 0.4 * y_true[i]
        if cancer_types[i] == "LUAD":
            base += 0.1
        elif cancer_types[i] == "COADREAD":
            base -= 0.05
        y_score[i] = np.clip(base + 0.1 * rng.randn(), 0, 1)

    analyzer = StratifiedAnalyzer(min_stratum_size=50)
    result = analyzer.analyze(y_true, y_score, cancer_types)

    print(analyzer.report(result))

    print("\nSelf-test complete.")
    sys.exit(0)
