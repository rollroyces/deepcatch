#!/usr/bin/env python3
"""
MODULE 4: Decision Curve Analysis
==================================

Reference: Vickers & Elkin (2006), "Decision Curve Analysis: A Novel Method
for Evaluating Prediction Models", Medical Decision Making 26:565-574.

Decision curve analysis evaluates clinical utility — not just statistical
accuracy. It answers: "At what range of risk thresholds does this model
provide net benefit over treating all or none?"

Net Benefit:
  NB = (TP - w·FP) / N
  where:
    w = p_t / (1 - p_t)  (odds at threshold p_t)
    p_t = probability threshold for action

Models compared:
  - Treat ALL:   NB = prevalence - w·(1 - prevalence)
  - Treat NONE:  NB = 0
  - Our MODEL:   NB computed from predictions

A model provides clinical value when its NB > max(Treat All, Treat None).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

from validation_framework import BootstrapCI

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray


# ── decision curve results ──────────────────────────────────────────────────


@dataclass
class ThresholdResult:
    """Net benefit at a single threshold."""
    threshold: float
    nb_model: float
    nb_treat_all: float
    nb_treat_none: float = 0.0
    nb_ci95: Tuple[float, float] = (float("nan"), float("nan"))


@dataclass
class DecisionCurveResult:
    """Full decision curve analysis."""
    thresholds: List[ThresholdResult]
    # Range where model provides net benefit
    useful_range_start: float
    useful_range_end: float
    # Interventions avoided per 100 patients (at optimal threshold)
    interventions_avoided: float
    interventions_avoided_ci95: Tuple[float, float]
    # Test tradeoff: 1 / (NB · (p_t / (1-p_t)))
    test_tradeoff_best: Optional[float] = None
    # Multi-model comparison
    model_comparisons: Optional[Dict[str, List[ThresholdResult]]] = None


class DecisionCurveAnalyzer:
    """Decision curve analysis for clinical prediction models.

    Usage:
        dca = DecisionCurveAnalyzer()
        result = dca.analyze(y_true, y_score)

        # Compare multiple models:
        result = dca.compare_models(
            y_true,
            {"Model A": scores_a, "Model B": scores_b},
        )
    """

    def __init__(
        self,
        thresholds: Optional[Array] = None,
        n_thresholds: int = 100,
        random_state: int = 42,
    ):
        """Args:
            thresholds: Array of threshold probabilities. If None, generates
                101 evenly spaced thresholds from 0.01 to 0.99.
            n_thresholds: Number of thresholds if auto-generating.
            random_state: Seed for bootstrap CIs.
        """
        self.thresholds = (
            thresholds
            if thresholds is not None
            else np.linspace(0.01, 0.99, n_thresholds)
        )
        self.random_state = random_state

    # ── public API ────────────────────────────────────────────────────────

    def analyze(
        self,
        y_true: Array,
        y_score: Array,
        model_name: str = "Model",
    ) -> DecisionCurveResult:
        """Compute decision curve for a single model.

        Args:
            y_true: Binary ground-truth labels.
            y_score: Predicted probabilities in [0, 1].
            model_name: Label for reporting.

        Returns:
            DecisionCurveResult with net benefit across all thresholds.
        """
        y_true = np.asarray(y_true).ravel().astype(int)
        y_score = np.asarray(y_score).ravel().astype(float)

        self._validate(y_true, y_score)

        prevalence = float(np.mean(y_true))
        n = len(y_true)

        threshold_results: List[ThresholdResult] = []
        useful_start = 1.0
        useful_end = 0.0

        for pt in self.thresholds:
            w = pt / (1 - pt) if pt < 1.0 else float("inf")

            # Model net benefit
            y_pred = (y_score >= pt).astype(int)
            tp = float(np.sum((y_pred == 1) & (y_true == 1)))
            fp = float(np.sum((y_pred == 1) & (y_true == 0)))
            nb_model = (tp - w * fp) / n

            # Treat all
            nb_all = prevalence - w * (1 - prevalence)

            # Bootstrap CI for model NB
            nb_ci = self._bootstrap_nb(y_true, y_score, pt, n_bootstrap=1000)

            tr = ThresholdResult(
                threshold=float(pt),
                nb_model=nb_model,
                nb_treat_all=nb_all,
                nb_treat_none=0.0,
                nb_ci95=nb_ci,
            )
            threshold_results.append(tr)

            # Track useful range: model NB > max(treat all, treat none)
            if nb_model > max(nb_all, 0.0):
                useful_start = min(useful_start, float(pt))
                useful_end = max(useful_end, float(pt))

        # Interventions avoided (at the midpoint of useful range or at
        # the point of max net benefit)
        if useful_start <= useful_end:
            best_pt = (useful_start + useful_end) / 2
        else:
            # Find threshold maximizing NB above treat-all
            nbs = np.array([t.nb_model - max(t.nb_treat_all, 0) for t in threshold_results])
            best_idx = int(np.argmax(nbs))
            best_pt = float(self.thresholds[best_idx])
            useful_start = best_pt
            useful_end = best_pt

        ia, ia_ci = self._interventions_avoided(
            y_true, y_score, best_pt, threshold_results
        )

        # Test tradeoff: 1 / (NB_model_at_pt * w)
        w_best = best_pt / (1 - best_pt) if best_pt < 1.0 else float("inf")
        best_tr = self._closest_threshold(best_pt, threshold_results)
        if best_tr and best_tr.nb_model > 0:
            test_tradeoff = 1 / (best_tr.nb_model * w_best) if w_best > 0 else float("inf")
        else:
            test_tradeoff = None

        return DecisionCurveResult(
            thresholds=threshold_results,
            useful_range_start=useful_start,
            useful_range_end=useful_end,
            interventions_avoided=ia,
            interventions_avoided_ci95=ia_ci,
            test_tradeoff_best=test_tradeoff,
        )

    def compare_models(
        self,
        y_true: Array,
        model_scores: Dict[str, Array],
    ) -> DecisionCurveResult:
        """Compare multiple models on the same decision curve.

        Args:
            y_true: Ground-truth labels.
            model_scores: Dict of name → predicted probabilities.

        Returns:
            DecisionCurveResult with model_comparisons populated.
        """
        y_true = np.asarray(y_true).ravel().astype(int)
        n = len(y_true)
        prevalence = float(np.mean(y_true))

        # Run analysis for first model as primary
        first_name = list(model_scores.keys())[0]
        result = self.analyze(y_true, model_scores[first_name], first_name)

        # Compute NB for other models at each threshold
        comparisons: Dict[str, List[ThresholdResult]] = {}
        for name, scores in model_scores.items():
            scores_arr = np.asarray(scores).ravel().astype(float)
            tr_list: List[ThresholdResult] = []

            for pt in self.thresholds:
                w = pt / (1 - pt) if pt < 1 else float("inf")
                y_pred = (scores_arr >= pt).astype(int)
                tp = float(np.sum((y_pred == 1) & (y_true == 1)))
                fp = float(np.sum((y_pred == 1) & (y_true == 0)))
                nb = (tp - w * fp) / n

                tr_list.append(ThresholdResult(
                    threshold=float(pt),
                    nb_model=nb,
                    nb_treat_all=prevalence - w * (1 - prevalence),
                ))

            comparisons[name] = tr_list

        result.model_comparisons = comparisons
        return result

    # ── helpers ───────────────────────────────────────────────────────────

    @staticmethod
    def _bootstrap_nb(
        y_true: Array,
        y_score: Array,
        threshold: float,
        n_bootstrap: int = 1000,
    ) -> Tuple[float, float]:
        """Bootstrap CI for net benefit at a given threshold."""
        n = len(y_true)
        rng = np.random.RandomState(42)
        nbs = np.empty(n_bootstrap)

        w = threshold / (1 - threshold) if threshold < 1.0 else float("inf")

        for i in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            yt = y_true[idx]
            ys = y_score[idx]
            yp = (ys >= threshold).astype(int)
            tp = float(np.sum((yp == 1) & (yt == 1)))
            fp = float(np.sum((yp == 1) & (yt == 0)))
            nbs[i] = (tp - w * fp) / n

        return (
            float(np.percentile(nbs, 2.5)),
            float(np.percentile(nbs, 97.5)),
        )

    def _interventions_avoided(
        self,
        y_true: Array,
        y_score: Array,
        pt: float,
        threshold_results: List[ThresholdResult],
    ) -> Tuple[float, Tuple[float, float]]:
        """Compute interventions avoided per 100 patients.

        Interventions avoided = (NB_model - NB_treat_all) × 100
        at the given threshold.
        """
        tr = self._closest_threshold(pt, threshold_results)
        if tr is None:
            return 0.0, (0.0, 0.0)

        ia = (tr.nb_model - tr.nb_treat_all) * 100

        # Bootstrap CI
        n = len(y_true)
        rng = np.random.RandomState(self.random_state)
        w = pt / (1 - pt) if pt < 1.0 else float("inf")
        ias = np.empty(1000)

        for i in range(1000):
            idx = rng.choice(n, size=n, replace=True)
            yt = y_true[idx]
            ys = y_score[idx]
            yp = (ys >= pt).astype(int)
            tp = float(np.sum((yp == 1) & (yt == 1)))
            fp = float(np.sum((yp == 1) & (yt == 0)))
            nb_m = (tp - w * fp) / n
            prev = float(np.mean(yt))
            nb_a = prev - w * (1 - prev)
            ias[i] = (nb_m - nb_a) * 100

        ia_ci = (
            float(np.percentile(ias, 2.5)),
            float(np.percentile(ias, 97.5)),
        )
        return ia, ia_ci

    @staticmethod
    def _closest_threshold(
        pt: float, results: List[ThresholdResult]
    ) -> Optional[ThresholdResult]:
        """Find the ThresholdResult closest to the given pt."""
        best = None
        best_dist = float("inf")
        for t in results:
            d = abs(t.threshold - pt)
            if d < best_dist:
                best_dist = d
                best = t
        return best

    @staticmethod
    def _validate(y_true: Array, y_score: Array) -> None:
        if len(y_true) < 20:
            raise ValueError(f"Need ≥ 20 samples; got {len(y_true)}")
        if len(np.unique(y_true)) < 2:
            raise ValueError("y_true must contain both classes")
        if np.any(y_score < 0) or np.any(y_score > 1):
            raise ValueError("y_score must be in [0, 1]")

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: DecisionCurveResult) -> str:
        """Publication-ready decision curve report."""
        lines = [
            "══ Decision Curve Analysis ══",
            f"  Useful threshold range:  [{result.useful_range_start:.3f}, {result.useful_range_end:.3f}]",
            f"  Interventions avoided:   {result.interventions_avoided:.1f} per 100 patients",
            f"      CI95:                 [{result.interventions_avoided_ci95[0]:.1f}, {result.interventions_avoided_ci95[1]:.1f}]",
        ]

        if result.test_tradeoff_best is not None:
            tt = result.test_tradeoff_best
            if tt < float("inf"):
                lines.append(f"  Test tradeoff:           1 test per ~{tt:.0f} patients to find one case")
            else:
                lines.append(f"  Test tradeoff:           not estimable at useful threshold")

        # Summary table of key thresholds
        lines.append(f"\n  Net Benefit at Key Thresholds:")
        key_pts = [0.01, 0.05, 0.10, 0.20, 0.30, 0.50]
        lines.append(f"  {'pt':>6s}  {'NB Model':>10s}  {'NB All':>10s}  {'NB None':>10s}  {'Δ vs Best':>10s}")
        lines.append(f"  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}")
        for kp in key_pts:
            tr = DecisionCurveAnalyzer._closest_threshold(kp, result.thresholds)
            if tr:
                nb_best = max(tr.nb_treat_all, 0.0)
                lines.append(
                    f"  {tr.threshold:6.4f}  {tr.nb_model:10.4f}  "
                    f"{tr.nb_treat_all:10.4f}  {0.0:10.4f}  "
                    f"{tr.nb_model - nb_best:10.4f}"
                )

        # Multi-model summary
        if result.model_comparisons:
            lines.append(f"\n  Model Comparison (at pt=0.10):")
            lines.append(f"  {'Model':>20s}  {'NB':>10s}")
            lines.append(f"  {'─'*20}  {'─'*10}")
            for name, tr_list in result.model_comparisons.items():
                tr = DecisionCurveAnalyzer._closest_threshold(0.10, tr_list)
                if tr:
                    lines.append(f"  {name:>20s}  {tr.nb_model:10.4f}")

        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 4: DecisionCurveAnalyzer — self-test")
    print("=" * 60)

    rng = np.random.RandomState(42)
    n = 500
    prevalence = 0.30

    y_true = (rng.rand(n) < prevalence).astype(int)

    # Good model: scores elevated for true positives
    y_score_good = np.clip(
        0.1 + 0.0 * rng.randn(n) + 0.6 * y_true + 0.1 * rng.randn(n), 0, 1
    )

    # Random model
    y_score_random = rng.rand(n)

    dca = DecisionCurveAnalyzer(thresholds=np.linspace(0.01, 0.50, 50))
    result = dca.analyze(y_true, y_score_good, model_name="Good Model")
    print(dca.report(result))

    print("\n── Model Comparison ──")
    comp = dca.compare_models(
        y_true,
        {"Good Model": y_score_good, "Random": y_score_random},
    )
    print(dca.report(comp))

    if result.useful_range_start <= result.useful_range_end:
        print(f"\n✓ Model provides net benefit in [{result.useful_range_start:.3f}, {result.useful_range_end:.3f}]")
    else:
        print("\n⚠ Model does not provide net benefit over treat-all/none")

    print("\nSelf-test complete.")
    sys.exit(0)
