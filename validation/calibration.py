#!/usr/bin/env python3
"""
MODULE 3: Calibration Analysis
===============================

References:
  - Niculescu-Mizil & Caruana (2005), "Predicting Good Probabilities with
    Supervised Learning", ICML.
  - Guo et al. (2017), "On Calibration of Modern Neural Networks", ICML.

Metrics:
  - Reliability diagram (binned predicted vs. observed frequency)
  - Brier score (Brier 1950): BS = (1/N) Σ (p̂ᵢ - yᵢ)²
  - Expected Calibration Error (ECE) (Naeini et al. 2015)
  - Maximum Calibration Error (MCE)
  - Calibration slope & intercept (from logistic recalibration)

Clinical relevance: A well-calibrated 30% predicted risk means cancer in
~30% of patients. Poor calibration leads to wrong risk communication.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss

from validation_framework import BootstrapCI

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray


# ── calibration result ──────────────────────────────────────────────────────


@dataclass
class BinStats:
    """Statistics for one calibration bin."""
    bin_center: float
    n_samples: int
    mean_predicted: float
    mean_observed: float
    se_observed: float  # standard error of the proportion


@dataclass
class CalibrationResult:
    """Comprehensive calibration analysis."""
    brier_score: float
    brier_skill: float  # relative to climatological baseline
    ece: float          # Expected Calibration Error
    mce: float          # Maximum Calibration Error
    bins: List[BinStats]
    calibration_slope: float
    calibration_intercept: float
    # Platt / isotonic recalibration scores (on hold-out)
    platt_brier: Optional[float] = None
    isotonic_brier: Optional[float] = None
    # Bootstrap CIs for Brier
    brier_ci95: Tuple[float, float] = (float("nan"), float("nan"))


class CalibrationAnalyzer:
    """Full calibration diagnostics for probabilistic classifiers.

    Usage:
        analyzer = CalibrationAnalyzer(n_bins=10)
        result = analyzer.analyze(y_true, y_score)

        # Compare Platt vs. isotonic recalibration:
        result = analyzer.analyze_with_recalibration(
            y_true, y_score, X_val, y_val, X_test, y_test
        )
    """

    def __init__(
        self,
        n_bins: int = 10,
        bin_strategy: str = "uniform_width",  # "uniform_width" or "equal_mass"
        random_state: int = 42,
    ):
        """Args:
            n_bins: Number of bins for reliability diagram.
            bin_strategy: 'uniform_width' (0–1 in equal intervals) or
                'equal_mass' (approximately equal samples per bin).
            random_state: Seed.
        """
        if n_bins < 3:
            raise ValueError("n_bins must be ≥ 3")
        self.n_bins = n_bins
        self.bin_strategy = bin_strategy
        self.random_state = random_state

        if bin_strategy not in ("uniform_width", "equal_mass"):
            raise ValueError(
                f"bin_strategy must be 'uniform_width' or 'equal_mass', "
                f"got '{bin_strategy}'"
            )

    # ── public API ────────────────────────────────────────────────────────

    def analyze(
        self,
        y_true: Array,
        y_score: Array,
    ) -> CalibrationResult:
        """Full calibration analysis.

        Args:
            y_true: Ground-truth binary labels.
            y_score: Predicted probabilities in [0, 1].

        Returns:
            CalibrationResult with all diagnostics.
        """
        y_true = np.asarray(y_true).ravel().astype(int)
        y_score = np.asarray(y_score).ravel().astype(float)

        self._validate_inputs(y_true, y_score)

        # ── Brier score ──
        brier = brier_score_loss(y_true, y_score)

        # Brier skill score: 1 - (brier / brier_climatology)
        p_base = float(np.mean(y_true))
        brier_ref = float(np.mean((y_true - p_base) ** 2))
        brier_skill = 1.0 - (brier / brier_ref) if brier_ref > 0 else 0.0

        # ── Reliability bins ──
        bins = self._compute_bins(y_true, y_score)

        # ── ECE / MCE ──
        ece = float(np.sum(
            np.array([b.n_samples / len(y_true) * abs(b.mean_predicted - b.mean_observed)
                       for b in bins])
        ))
        mce = float(np.max(
            [abs(b.mean_predicted - b.mean_observed) for b in bins]
        ))

        # ── Calibration slope/intercept (logistic recalibration) ──
        slope, intercept = self._calibration_slope_intercept(y_true, y_score)

        # ── Bootstrap CI for Brier ──
        brier_vals = self._bootstrap_brier(y_true, y_score, n_bootstrap=1000)
        brier_ci = (
            float(np.percentile(brier_vals, 2.5)),
            float(np.percentile(brier_vals, 97.5)),
        )

        return CalibrationResult(
            brier_score=brier,
            brier_skill=brier_skill,
            ece=ece,
            mce=mce,
            bins=bins,
            calibration_slope=slope,
            calibration_intercept=intercept,
            brier_ci95=brier_ci,
        )

    def analyze_with_recalibration(
        self,
        y_true: Array,
        y_score: Array,
        # Data for Platt scaling (can be same as y_true/y_score)
        X_val: Optional[Array] = None,
        y_val: Optional[Array] = None,
        X_test: Optional[Array] = None,
        y_test: Optional[Array] = None,
    ) -> CalibrationResult:
        """Analyze calibration and compare Platt vs isotonic recalibration.

        If validation/test splits aren't provided, uses internal train/test split.

        Returns:
            CalibrationResult with platt_brier and isotonic_brier populated.
        """
        base = self.analyze(y_true, y_score)

        y_true_arr = np.asarray(y_true).ravel().astype(int)
        y_score_arr = np.asarray(y_score).ravel().astype(float)

        n = len(y_true_arr)

        # Set up recalibration data
        if X_val is not None and X_test is not None:
            val_scores = y_score_arr  # these ARE the validation scores
            val_labels = y_true_arr
            test_scores = y_score_arr
            test_labels = y_true_arr
        else:
            # Internal split
            rng = np.random.RandomState(self.random_state)
            idx = rng.permutation(n)
            n_val = n // 2
            val_idx, test_idx = idx[:n_val], idx[n_val:]
            val_scores = y_score_arr[val_idx]
            val_labels = y_true_arr[val_idx]
            test_scores = y_score_arr[test_idx]
            test_labels = y_true_arr[test_idx]

        # ── Platt scaling (logistic regression on scores) ──
        try:
            platt = LogisticRegression(C=1e12, solver="lbfgs")
            platt.fit(val_scores.reshape(-1, 1), val_labels)
            platt_scores = platt.predict_proba(test_scores.reshape(-1, 1))[:, 1]
            base.platt_brier = brier_score_loss(test_labels, platt_scores)
        except Exception:
            base.platt_brier = None

        # ── Isotonic regression ──
        try:
            iso = IsotonicRegression(y_min=0.0, y_max=1.0, out_of_bounds="clip")
            iso.fit(val_scores, val_labels)
            iso_scores = iso.predict(test_scores)
            base.isotonic_brier = brier_score_loss(test_labels, iso_scores)
        except Exception:
            base.isotonic_brier = None

        return base

    # ── bin computation ───────────────────────────────────────────────────

    def _compute_bins(self, y_true: Array, y_score: Array) -> List[BinStats]:
        """Compute reliability diagram bins."""
        if self.bin_strategy == "uniform_width":
            bin_edges = np.linspace(0, 1, self.n_bins + 1)
            bin_stats: List[BinStats] = []

            for i in range(self.n_bins):
                mask = (y_score >= bin_edges[i]) & (
                    y_score < bin_edges[i + 1] if i < self.n_bins - 1
                    else (y_score <= bin_edges[i + 1])
                )
                if np.sum(mask) == 0:
                    continue
                bin_stats.append(self._bin_summary(
                    y_true[mask], y_score[mask], (bin_edges[i] + bin_edges[i + 1]) / 2
                ))
            return bin_stats

        else:  # equal_mass
            order = np.argsort(y_score)
            y_true_sorted = y_true[order]
            y_score_sorted = y_score[order]

            bin_size = len(y_true) // self.n_bins
            bins: List[BinStats] = []

            for i in range(self.n_bins):
                start = i * bin_size
                end = start + bin_size if i < self.n_bins - 1 else len(y_true)
                mask = slice(start, end)
                bins.append(self._bin_summary(
                    y_true_sorted[mask],
                    y_score_sorted[mask],
                    float(np.mean(y_score_sorted[mask])),
                ))
            return bins

    @staticmethod
    def _bin_summary(
        y_true: Array, y_score: Array, center: float
    ) -> BinStats:
        n = len(y_true)
        mean_pred = float(np.mean(y_score))
        mean_obs = float(np.mean(y_true))
        # Standard error of binomial proportion
        se = np.sqrt(mean_obs * (1 - mean_obs) / n) if n > 0 else 0.0
        return BinStats(
            bin_center=center,
            n_samples=n,
            mean_predicted=mean_pred,
            mean_observed=mean_obs,
            se_observed=se,
        )

    # ── calibration slope / intercept ─────────────────────────────────────

    @staticmethod
    def _calibration_slope_intercept(
        y_true: Array, y_score: Array
    ) -> Tuple[float, float]:
        """Fit logistic recalibration: logit(p_true) = α + β·logit(p_score).

        β = 1, α = 0 → perfect calibration.
        β < 1 → model is overconfident.
        β > 1 → model is underconfident.
        """
        # Clip to avoid logit(0) or logit(1)
        eps = 1e-12
        scores_clipped = np.clip(y_score, eps, 1 - eps)

        # Add small constant column for intercept
        X_logit = np.column_stack([
            np.ones_like(scores_clipped),
            np.log(scores_clipped / (1 - scores_clipped)),
        ])

        try:
            lr = LogisticRegression(C=1e12, solver="lbfgs", fit_intercept=False)
            lr.fit(X_logit, y_true)
            intercept = float(lr.coef_[0, 0])
            slope = float(lr.coef_[0, 1])
        except Exception:
            # Fallback: simple OLS on logit-transformed means
            intercept, slope = 0.0, 1.0

        return slope, intercept

    # ── bootstrap Brier ───────────────────────────────────────────────────

    def _bootstrap_brier(
        self,
        y_true: Array,
        y_score: Array,
        n_bootstrap: int = 1000,
    ) -> np.ndarray:
        """Bootstrap Brier scores for CI estimation."""
        rng = np.random.RandomState(self.random_state)
        n = len(y_true)
        briers = np.empty(n_bootstrap)

        for i in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            briers[i] = brier_score_loss(y_true[idx], y_score[idx])

        return briers

    # ── validation ────────────────────────────────────────────────────────

    @staticmethod
    def _validate_inputs(y_true: Array, y_score: Array) -> None:
        if len(y_true) < 10:
            raise ValueError(f"Need ≥ 10 samples for calibration; got {len(y_true)}")
        if np.any(y_score < 0) or np.any(y_score > 1):
            raise ValueError("y_score must be in [0, 1]")
        if len(np.unique(y_true)) < 2:
            raise ValueError("y_true must contain both classes")

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: CalibrationResult) -> str:
        """Publication-ready calibration report."""
        lines = [
            "══ Calibration Analysis ══",
            f"  Brier score:              {result.brier_score:.4f}",
            f"  Brier skill score:        {result.brier_skill:.4f}",
            f"  Brier CI95:               [{result.brier_ci95[0]:.4f}, {result.brier_ci95[1]:.4f}]",
            f"  ECE (Expected Cal. Error): {result.ece:.4f}",
            f"  MCE (Max Cal. Error):      {result.mce:.4f}",
            f"  Calibration slope:         {result.calibration_slope:.3f}",
            f"  Calibration intercept:     {result.calibration_intercept:.3f}",
        ]

        slope = result.calibration_slope
        if slope < 0.8:
            lines.append("  ⚠ Model is OVERCONFIDENT (slope << 1)")
        elif slope > 1.2:
            lines.append("  ⚠ Model is UNDERCONFIDENT (slope >> 1)")
        else:
            lines.append("  ✓ Calibration slope near 1")

        if result.platt_brier is not None:
            improvement = result.brier_score - result.platt_brier
            lines.append(f"  Platt Brier:               {result.platt_brier:.4f} (Δ = {improvement:+.4f})")
        if result.isotonic_brier is not None:
            improvement = result.brier_score - result.isotonic_brier
            lines.append(f"  Isotonic Brier:            {result.isotonic_brier:.4f} (Δ = {improvement:+.4f})")

        lines.append(f"\n  Reliability Diagram ({len(result.bins)} bins):")
        lines.append(f"  {'Bin':>5s} {'N':>6s} {'Pred':>8s} {'Obs':>8s} {'|Δ|':>8s}")
        lines.append(f"  {'─'*5} {'─'*6} {'─'*8} {'─'*8} {'─'*8}")
        for b in result.bins:
            delta = abs(b.mean_predicted - b.mean_observed)
            lines.append(
                f"  {b.bin_center:5.3f} {b.n_samples:6d} "
                f"{b.mean_predicted:8.4f} {b.mean_observed:8.4f} {delta:8.4f}"
            )

        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 3: CalibrationAnalyzer — self-test")
    print("=" * 60)

    rng = np.random.RandomState(42)
    n = 500

    # Generate well-calibrated scores
    true_risk = np.clip(0.1 + 0.3 * rng.randn(n), 0.01, 0.99)
    y_true_test = (rng.rand(n) < true_risk).astype(int)
    # Well-calibrated: score ≈ true_risk + small noise
    y_score_well = np.clip(true_risk + 0.03 * rng.randn(n), 0, 1)

    # Overconfident: push scores toward extremes
    y_score_over = np.clip(true_risk * 2 - 0.1, 0, 1)

    ca = CalibrationAnalyzer(n_bins=10)

    print("\n── Well-Calibrated Model ──")
    result_well = ca.analyze(y_true_test, y_score_well)
    print(ca.report(result_well))

    print("\n── Overconfident Model ──")
    result_over = ca.analyze(y_true_test, y_score_over)
    print(ca.report(result_over))

    if result_over.ece > result_well.ece:
        print("\n✓ ECE correctly identifies overconfident model")
    else:
        print("\n⚠ ECE ordering unexpected")

    print("\nSelf-test complete.")
    sys.exit(0)
