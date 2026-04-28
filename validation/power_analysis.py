#!/usr/bin/env python3
"""
MODULE 9: Sample Size & Power Analysis
=======================================

Reference: Cohen (1988), "Statistical Power Analysis for the Behavioral
Sciences", 2nd ed. Lawrence Erlbaum.

For each experiment, answers two critical questions:

1. Given observed effect size, what sample size is needed for 80% power at α=0.05?
2. Given planned sample size, what minimum effect size is detectable?

This is essential for:
  - Justifying that experiments are adequately powered
  - Honestly reporting when findings need larger samples
  - Planning future validation studies
  - Meeting reviewer expectations at top journals
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from scipy.optimize import brentq

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray


# ── power analysis result ───────────────────────────────────────────────────


@dataclass
class PowerCurvePoint:
    """One point on a power curve."""
    sample_size: int
    power: float
    detectable_effect: float


@dataclass
class ExperimentPower:
    """Power analysis for one experiment/metric."""
    experiment_name: str
    metric: str
    observed_effect_size: float  # e.g., Cohen's d
    observed_effect_se: float   # standard error of observed effect
    observed_n: int              # actual sample size used
    observed_power: float        # post-hoc power at α=0.05
    required_n_for_80pct: int    # sample size for 80% power
    required_n_for_90pct: int    # sample size for 90% power
    minimum_detectable_effect: float  # at current n, 80% power
    power_curve: List[PowerCurvePoint]
    # Interpretation
    adequately_powered: bool  # observed_power >= 0.80
    recommendation: str


@dataclass
class PowerAnalysisResult:
    """Full power analysis across experiments."""
    experiments: List[ExperimentPower]
    summary: Dict[str, Any]
    n_experiments: int
    n_adequately_powered: int
    n_underpowered: int


class PowerAnalyzer:
    """Sample size and statistical power analysis for model validation.

    Usage:
        pa = PowerAnalyzer()
        result = pa.analyze(experiments)
        result = pa.analyze_from_auc(auc, n_samples, positives_ratio)
    """

    def __init__(
        self,
        alpha: float = 0.05,
        target_power: float = 0.80,
        min_n: int = 10,
        max_n: int = 100000,
        n_points: int = 50,
    ):
        """Args:
            alpha: Significance level (default 0.05).
            target_power: Target statistical power (default 0.80).
            min_n, max_n: Range for sample size search.
            n_points: Number of points on power curve.
        """
        self.alpha = alpha
        self.target_power = target_power
        self.min_n = min_n
        self.max_n = max_n
        self.n_points = n_points

    # ── main API ──────────────────────────────────────────────────────────

    def analyze(
        self,
        experiments: List[Dict[str, Any]],
    ) -> PowerAnalysisResult:
        """Full power analysis for multiple experiments.

        Each experiment dict should contain:
          - 'name': str
          - 'metric': str (e.g., 'AUC', 'sensitivity')
          - 'observed_effect': float (effect size or AUC)
          - 'observed_se': float (standard error)
          - 'n_samples': int
          - Optional: 'n_pos': int, 'n_neg': int (for AUC-specific calculations)

        Returns:
            PowerAnalysisResult with per-experiment power estimates.
        """
        results: List[ExperimentPower] = []

        for exp in experiments:
            name = exp.get("name", "Unknown")
            metric = exp.get("metric", "effect")
            observed_effect = float(exp["observed_effect"])
            observed_se = float(exp.get("observed_se", 0.05))
            n_samples = int(exp["n_samples"])

            # Convert AUC to Cohen's d for power calculations
            if metric.lower() in ("auc", "roc_auc", "auroc"):
                cohens_d = self._auc_to_cohens_d(observed_effect)
                # For AUC power, effective n depends on positive/negative split
                n_pos = exp.get("n_pos", n_samples // 2)
                n_neg = exp.get("n_neg", n_samples // 2)
                effective_n = self._effective_n_auc(n_pos, n_neg)
            else:
                cohens_d = observed_effect
                effective_n = n_samples

            # Observed (post-hoc) power
            observed_power = self._compute_power(
                cohens_d, effective_n, self.alpha
            )

            # Required sample size for various power levels
            required_n_80 = self._required_n(
                cohens_d, alpha=self.alpha, power=0.80
            )
            required_n_90 = self._required_n(
                cohens_d, alpha=self.alpha, power=0.90
            )

            # Minimum detectable effect at current n
            mde = self._minimum_detectable_effect(
                effective_n, alpha=self.alpha, power=self.target_power
            )

            # Power curve
            n_range = np.logspace(
                np.log10(max(self.min_n, effective_n * 0.1)),
                np.log10(min(self.max_n, max(required_n_90 * 2, effective_n * 5))),
                self.n_points,
            ).astype(int)
            n_range = np.unique(np.clip(n_range, self.min_n, self.max_n))

            power_curve: List[PowerCurvePoint] = []
            for n in n_range:
                power = self._compute_power(cohens_d, n, self.alpha)
                mde_at_n = self._minimum_detectable_effect(
                    n, alpha=self.alpha, power=self.target_power
                )
                power_curve.append(PowerCurvePoint(
                    sample_size=int(n),
                    power=power,
                    detectable_effect=mde_at_n,
                ))

            # Recommendation
            adequately_powered = observed_power >= self.target_power
            if adequately_powered:
                recommendation = (
                    f"Adequately powered (power={observed_power:.2f}). "
                    f"No additional samples needed for {metric}."
                )
            elif required_n_80 <= self.max_n:
                shortage = required_n_80 - effective_n
                recommendation = (
                    f"UNDERPOWERED (power={observed_power:.2f}). "
                    f"Need ~{shortage} more samples for 80% power. "
                    f"Current findings should be interpreted cautiously."
                )
            else:
                recommendation = (
                    f"SEVERELY UNDERDEVELOPED (power={observed_power:.2f}). "
                    f"Effect too small to be confirmed with feasible sample sizes. "
                    f"Consider whether this metric is clinically meaningful."
                )

            results.append(ExperimentPower(
                experiment_name=name,
                metric=metric,
                observed_effect_size=cohens_d,
                observed_effect_se=observed_se,
                observed_n=effective_n,
                observed_power=observed_power,
                required_n_for_80pct=required_n_80,
                required_n_for_90pct=required_n_90,
                minimum_detectable_effect=mde,
                power_curve=power_curve,
                adequately_powered=adequately_powered,
                recommendation=recommendation,
            ))

        # Summary
        n_adequate = sum(1 for r in results if r.adequately_powered)
        n_under = sum(1 for r in results if not r.adequately_powered)

        summary = {
            "n_experiments": len(results),
            "n_adequately_powered": n_adequate,
            "n_underpowered": n_under,
            "pct_adequately_powered": 100 * n_adequate / len(results) if results else 0,
            "largest_required_n": max((r.required_n_for_80pct for r in results), default=0),
        }

        return PowerAnalysisResult(
            experiments=results,
            summary=summary,
            n_experiments=len(results),
            n_adequately_powered=n_adequate,
            n_underpowered=n_under,
        )

    def analyze_from_auc(
        self,
        auc: float,
        n_samples: int,
        pos_ratio: Optional[float] = None,
        experiment_name: str = "AUC Analysis",
        auc_se: float = 0.03,
    ) -> ExperimentPower:
        """Convenience: power analysis from AUC estimate.

        Args:
            auc: Observed AUC value.
            n_samples: Total sample size.
            pos_ratio: Fraction of positive samples (default: infer).
            experiment_name: Label.
            auc_se: Standard error of AUC estimate.

        Returns:
            ExperimentPower.
        """
        if pos_ratio is None:
            pos_ratio = 0.5  # assumption

        n_pos = int(n_samples * pos_ratio)
        n_neg = n_samples - n_pos

        return self.analyze([{
            "name": experiment_name,
            "metric": "AUC",
            "observed_effect": auc,
            "observed_se": auc_se,
            "n_samples": n_samples,
            "n_pos": n_pos,
            "n_neg": n_neg,
        }]).experiments[0]

    def analyze_from_sensitivity(
        self,
        sensitivity: float,
        specificity: float,
        n_pos: int,
        n_neg: int,
        experiment_name: str = "Sensitivity Analysis",
    ) -> ExperimentPower:
        """Power analysis for sensitivity, accounting for test characteristics.

        For sensitivity: effective n = n_pos
        Effect size: sensitivity - null_hypothesis_sensitivity
        """
        # Sensitivity better than random (50%) or better than a threshold
        null_sens = 0.5
        # Convert sensitivity difference to Cohen's h (arcsine transform)
        cohens_h = 2 * (np.arcsin(np.sqrt(sensitivity)) - np.arcsin(np.sqrt(null_sens)))

        return self.analyze([{
            "name": experiment_name,
            "metric": "Sensitivity",
            "observed_effect": abs(cohens_h),
            "observed_se": np.sqrt(sensitivity * (1 - sensitivity) / n_pos),
            "n_samples": n_pos,
        }]).experiments[0]

    # ── core power calculations ───────────────────────────────────────────

    def _compute_power(self, effect_size: float, n: int, alpha: float) -> float:
        """Compute statistical power for a two-sample test.

        Using non-central t-distribution.
        """
        if n < 2 or effect_size <= 0 or np.isnan(effect_size):
            return 0.0

        df = 2 * n - 2  # two independent samples
        # Non-centrality parameter
        ncp = effect_size * np.sqrt(n / 2)
        # Critical t-value
        t_crit = stats.t.ppf(1 - alpha / 2, df)
        # Power = P(|t| > t_crit | non-centrality = ncp)
        power = 1 - stats.nct.cdf(t_crit, df, ncp) + stats.nct.cdf(-t_crit, df, ncp)
        return float(min(max(power, 0.0), 1.0))

    def _required_n(
        self,
        effect_size: float,
        alpha: float = 0.05,
        power: float = 0.80,
    ) -> int:
        """Find minimum sample size per group for desired power.

        Uses numerical root-finding on the power function.
        """
        if effect_size <= 0 or np.isnan(effect_size):
            return self.max_n

        def power_diff(n):
            if n < 2:
                return -power
            return self._compute_power(effect_size, int(n), alpha) - power

        # Search for root
        try:
            n_required = brentq(
                power_diff,
                self.min_n,
                self.max_n,
                maxiter=100,
            )
            return max(self.min_n, int(np.ceil(n_required)))
        except (ValueError, RuntimeError):
            # Effect too small — can't achieve target power within max_n
            return self.max_n

    def _minimum_detectable_effect(
        self,
        n: int,
        alpha: float = 0.05,
        power: float = 0.80,
    ) -> float:
        """Minimum detectable effect size (Cohen's d) at given n and power."""
        if n < 2:
            return float("inf")

        def effect_for_power(d):
            return self._compute_power(d, n, alpha) - power

        try:
            mde = brentq(effect_for_power, 0.001, 5.0, maxiter=100)
            return float(mde)
        except (ValueError, RuntimeError):
            return 5.0  # very large effect needed

    # ── conversions ───────────────────────────────────────────────────────

    @staticmethod
    def _auc_to_cohens_d(auc: float) -> float:
        """Convert AUC to Cohen's d.

        AUC = Φ(d / √2)  where Φ is the standard normal CDF.
        Therefore: d = √2 · Φ⁻¹(AUC)

        Reference: Ruscio (2008), "A Probability-Based Measure of Effect
        Size: Robustness to Base Rates and Other Factors", Psych Methods.
        """
        # Clamp AUC to avoid numerical issues at boundaries
        auc = max(0.5001, min(0.9999, auc))
        return np.sqrt(2) * stats.norm.ppf(auc)

    @staticmethod
    def _cohens_d_to_auc(d: float) -> float:
        """Convert Cohen's d to AUC."""
        return float(stats.norm.cdf(d / np.sqrt(2)))

    @staticmethod
    def _effective_n_auc(n_pos: int, n_neg: int) -> int:
        """Effective sample size for AUC comparison.

        For AUC, the effective sample size for power is approximately
        the harmonic mean of positives and negatives.

        Reference: Hanley & McNeil (1982), Radiology 143:29-36.
        """
        if n_pos < 1 or n_neg < 1:
            return 0
        return int(2 / (1.0 / n_pos + 1.0 / n_neg))

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: PowerAnalysisResult) -> str:
        """Publication-ready power analysis report."""
        lines = [
            "══ Power & Sample Size Analysis ══",
            "",
            f"  Experiments analyzed: {result.n_experiments}",
            f"  Adequately powered (≥80%): {result.n_adequately_powered}",
            f"  Underpowered: {result.n_underpowered}",
            f"  Overall powered: {result.summary['pct_adequately_powered']:.0f}%",
            "",
            f"  {'Experiment':>30s} {'Effect':>8s} {'n':>6s} {'Power':>7s} {'n(80%)':>8s} {'n(90%)':>8s} {'MDE':>7s}",
            f"  {'─'*30} {'─'*8} {'─'*6} {'─'*7} {'─'*8} {'─'*8} {'─'*7}",
        ]

        for exp in result.experiments:
            status = "✓" if exp.adequately_powered else "✗"
            n80 = str(exp.required_n_for_80pct) if exp.required_n_for_80pct < 100000 else ">100K"
            n90 = str(exp.required_n_for_90pct) if exp.required_n_for_90pct < 100000 else ">100K"
            lines.append(
                f"  {exp.experiment_name:>30s} {exp.observed_effect_size:+8.3f} "
                f"{exp.observed_n:6d} {exp.observed_power:7.3f} "
                f"{n80:>8s} {n90:>8s} "
                f"{exp.minimum_detectable_effect:7.3f} {status}"
            )

        lines.append(f"\n  Power Curves (key experiments):")
        for exp in result.experiments[:6]:  # Show top 6
            lines.append(f"\n  ── {exp.experiment_name} ({exp.metric}) ──")
            lines.append(f"  Observed effect: {exp.observed_effect_size:.3f}, n={exp.observed_n}, power={exp.observed_power:.2f}")
            lines.append(f"  Required n for 80% power: {exp.required_n_for_80pct}")
            lines.append(f"  Required n for 90% power: {exp.required_n_for_90pct}")

            # Power curve summary
            if exp.power_curve:
                # Show key points
                key_powers = [0.50, 0.80, 0.90, 0.95]
                points_shown = 0
                for kp in key_powers:
                    for pt in exp.power_curve:
                        if abs(pt.power - kp) < 0.02:
                            lines.append(f"    n={pt.sample_size:6d} → power={pt.power:.2f}")
                            points_shown += 1
                            break
                if points_shown == 0:
                    # Show range
                    lo = exp.power_curve[0]
                    hi = exp.power_curve[-1]
                    lines.append(f"    Range: n={lo.sample_size} (power={lo.power:.2f}) to n={hi.sample_size} (power={hi.power:.2f})")

            lines.append(f"  {exp.recommendation}")

        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 9: PowerAnalyzer — self-test")
    print("=" * 60)

    pa = PowerAnalyzer(alpha=0.05, target_power=0.80)

    experiments = [
        {
            "name": "CET Detection (AUC=0.95)",
            "metric": "AUC",
            "observed_effect": 0.95,
            "observed_se": 0.02,
            "n_samples": 200,
            "n_pos": 60,
            "n_neg": 140,
        },
        {
            "name": "GNN Fusion (AUC=0.75)",
            "metric": "AUC",
            "observed_effect": 0.75,
            "observed_se": 0.04,
            "n_samples": 600,
            "n_pos": 180,
            "n_neg": 420,
        },
        {
            "name": "Contrastive Learner (AUC=0.68)",
            "metric": "AUC",
            "observed_effect": 0.68,
            "observed_se": 0.05,
            "n_samples": 100,
            "n_pos": 30,
            "n_neg": 70,
        },
        {
            "name": "Bayesian Caller (AUC=0.88)",
            "metric": "AUC",
            "observed_effect": 0.88,
            "observed_se": 0.03,
            "n_samples": 500,
            "n_pos": 150,
            "n_neg": 350,
        },
    ]

    result = pa.analyze(experiments)
    print(pa.report(result))

    # Verify: larger AUC should need fewer samples
    cet = result.experiments[0]
    gnn = result.experiments[1]
    if cet.required_n_for_80pct < gnn.required_n_for_80pct:
        print(f"\n✓ CET (AUC=0.95) needs fewer samples than GNN (AUC=0.75)")
    else:
        print(f"\n⚠ Sample size ordering unexpected")

    print("\nSelf-test complete.")
    sys.exit(0)
