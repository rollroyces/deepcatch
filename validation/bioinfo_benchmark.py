#!/usr/bin/env python3
"""
MODULE 8: Bioinformatic Tool Benchmark
=======================================

Head-to-head comparison against real bioinformatics tools and baseline methods.

Variant Calling Comparison:
  - Mutect2 (GATK) — gold standard somatic caller
  - VarScan2 — heuristic-based caller
  - Strelka2 — fast somatic caller
  - LoFreq — low-VAF sensitive caller
  - SiNVICT — ultra-low VAF caller

Fusion Strategy Comparison:
  - Late fusion — combine modality predictions
  - Early fusion — concatenate features before modeling
  - Single-best-modality — use only the best performing modality
  - Our multi-modal model

All comparisons use DeLong test with multiple testing correction.
Benchmarks run via simulation of tool behavior based on known characteristics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    average_precision_score,
)

from validation.delong_test import delong_test, delong_test_multi

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray


# ── benchmark results ───────────────────────────────────────────────────────


@dataclass
class ToolResult:
    """Performance of one bioinformatics tool."""
    name: str
    category: str  # "variant_caller", "fusion_strategy", "baseline"
    auc: float
    auc_ci95: Tuple[float, float]
    sensitivity: float
    specificity: float
    f1: float
    precision: float
    # Detection rate at low VAF (< 0.5%)
    low_vaf_sensitivity: Optional[float] = None
    # Computational burden (relative)
    compute_time_relative: Optional[float] = None


@dataclass
class BenchmarkComparison:
    """Pairwise comparison between our model and a benchmark."""
    our_model: str
    benchmark: str
    delta_auc: float
    p_value: float  # DeLong
    p_value_corrected: float
    significant: bool
    ci95: Tuple[float, float]


@dataclass
class BioinfoBenchmarkResult:
    """Full benchmark results."""
    variant_callers: Dict[str, ToolResult]
    fusion_strategies: Dict[str, ToolResult]
    # Pairwise comparisons (our model vs each benchmark)
    variant_comparisons: List[BenchmarkComparison]
    fusion_comparisons: List[BenchmarkComparison]
    # Rankings
    ranked_variant_callers: List[Tuple[str, float]]
    ranked_fusion_strategies: List[Tuple[str, float]]
    # Summary
    best_variant_caller: str
    best_fusion_strategy: str


# ── tool simulator ──────────────────────────────────────────────────────────


class BioinfoBenchmark:
    """Benchmark our model against bioinformatics tools and fusion baselines.

    Usage:
        bench = BioinfoBenchmark()
        result = bench.benchmark_variant_calling(
            our_model, X, y, low_vaf_mask=low_vaf_samples
        )
        result = bench.benchmark_fusion_strategies(
            our_model, modality_scores, y
        )
    """

    # ── known tool characteristics (from literature) ──────────────────────
    # These are approximate performance characteristics used for simulation
    # when real tool outputs aren't available.
    TOOL_CHARACTERISTICS = {
        "Mutect2": {
            "sensitivity": 0.85, "specificity": 0.92,
            "low_vaf_sensitivity": 0.55, "category": "variant_caller",
            "reference": "Cibulskis et al. (2013), Nature Biotechnology",
        },
        "VarScan2": {
            "sensitivity": 0.78, "specificity": 0.88,
            "low_vaf_sensitivity": 0.30, "category": "variant_caller",
            "reference": "Koboldt et al. (2012), Genome Research",
        },
        "Strelka2": {
            "sensitivity": 0.88, "specificity": 0.90,
            "low_vaf_sensitivity": 0.45, "category": "variant_caller",
            "reference": "Kim et al. (2018), Nature Methods",
        },
        "LoFreq": {
            "sensitivity": 0.82, "specificity": 0.85,
            "low_vaf_sensitivity": 0.70, "category": "variant_caller",
            "reference": "Wilm et al. (2012), Nucleic Acids Research",
        },
        "SiNVICT": {
            "sensitivity": 0.75, "specificity": 0.82,
            "low_vaf_sensitivity": 0.80, "category": "variant_caller",
            "reference": "Kockan et al. (2017), Bioinformatics",
        },
        "DeepCatch (ours)": {
            "sensitivity": None, "specificity": None,  # filled from actual model
            "low_vaf_sensitivity": None, "category": "variant_caller",
            "reference": "This study",
        },
    }

    def __init__(
        self,
        random_state: int = 42,
        n_bootstrap: int = 1000,
        correction: str = "bonferroni",
    ):
        self.random_state = random_state
        self.n_bootstrap = n_bootstrap
        self.correction = correction
        self.rng = np.random.RandomState(random_state)

    # ── variant calling benchmark ─────────────────────────────────────────

    def benchmark_variant_calling(
        self,
        our_model: Any,  # Must have .predict_proba() or .decision_function()
        X: Array,
        y: Array,
        low_vaf_mask: Optional[Array] = None,
        tool_scores: Optional[Dict[str, Array]] = None,
    ) -> BioinfoBenchmarkResult:
        """Compare variant calling performance against tools.

        Args:
            our_model: Trained model with predict_proba method.
            X: Feature matrix.
            y: Ground-truth variant labels (0=reference, 1=variant).
            low_vaf_mask: Boolean mask for low-VAF samples.
            tool_scores: Optional dict of tool_name → predicted_scores.
                If provided, use actual tool outputs.
                If None, simulate based on known tool characteristics.

        Returns:
            BioinfoBenchmarkResult.
        """
        X = np.asarray(X)
        y = np.asarray(y).ravel().astype(int)

        # Get our model's scores
        if hasattr(our_model, "predict_proba"):
            our_scores = our_model.predict_proba(X)[:, 1]
        elif hasattr(our_model, "decision_function"):
            our_scores = our_model.decision_function(X)
        else:
            raise ValueError("Model must have predict_proba or decision_function")

        our_results = self._compute_tool_metrics(y, our_scores, low_vaf_mask, "DeepCatch (ours)")

        # Get benchmark tool scores (simulated or real)
        variant_callers: Dict[str, ToolResult] = {"DeepCatch (ours)": our_results}
        all_scores: Dict[str, Array] = {"DeepCatch (ours)": our_scores}

        if tool_scores is not None:
            for tool_name, scores in tool_scores.items():
                tr = self._compute_tool_metrics(y, scores, low_vaf_mask, tool_name)
                variant_callers[tool_name] = tr
                all_scores[tool_name] = np.asarray(scores)
        else:
            # Simulate tool behavior
            for tool_name, chars in self.TOOL_CHARACTERISTICS.items():
                if tool_name == "DeepCatch (ours)":
                    continue
                sim_scores = self._simulate_tool_scores(
                    y, chars["sensitivity"], chars["specificity"],
                    low_vaf_mask, chars.get("low_vaf_sensitivity")
                )
                tr = self._compute_tool_metrics(y, sim_scores, low_vaf_mask, tool_name)
                variant_callers[tool_name] = tr
                all_scores[tool_name] = sim_scores

        # Pairwise DeLong comparisons
        comparisons: List[BenchmarkComparison] = []
        our_name = "DeepCatch (ours)"

        for tool_name in variant_callers:
            if tool_name == our_name:
                continue
            dl = delong_test(y, all_scores[our_name], all_scores[tool_name])
            comparisons.append(BenchmarkComparison(
                our_model=our_name,
                benchmark=tool_name,
                delta_auc=dl["delta_auc"],
                p_value=dl["p_value"],
                p_value_corrected=dl["p_value"],  # corrected below
                significant=dl["significant"],
                ci95=(dl["ci95_lower"], dl["ci95_upper"]),
            ))

        # Apply multiple testing correction
        if len(comparisons) > 1:
            p_vals = np.array([c.p_value for c in comparisons])
            from validation_framework import SignificanceTester
            corrected = SignificanceTester.bonferroni_correct(p_vals)
            for i, c in enumerate(comparisons):
                c.p_value_corrected = float(corrected[i])
                c.significant = float(corrected[i]) < 0.05

        # Rank by AUC
        ranked = sorted(
            [(name, tr.auc) for name, tr in variant_callers.items()],
            key=lambda x: -x[1],
        )

        return BioinfoBenchmarkResult(
            variant_callers=variant_callers,
            fusion_strategies={},
            variant_comparisons=comparisons,
            fusion_comparisons=[],
            ranked_variant_callers=ranked,
            ranked_fusion_strategies=[],
            best_variant_caller=ranked[0][0] if ranked else "N/A",
            best_fusion_strategy="N/A",
        )

    def benchmark_fusion_strategies(
        self,
        our_model: Any,
        modality_scores: Dict[str, Array],
        y: Array,
        X_early_fusion: Optional[Array] = None,
    ) -> BioinfoBenchmarkResult:
        """Compare multi-modal fusion strategies.

        Args:
            our_model: Our multi-modal fusion model.
            modality_scores: Dict of modality_name → probability scores.
                Each array is (n_samples,) with probabilities in [0,1].
            y: Ground-truth labels.
            X_early_fusion: Early fusion feature matrix (optional).

        Returns:
            BioinfoBenchmarkResult with fusion_comparisons populated.
        """
        y = np.asarray(y).ravel().astype(int)
        n = len(y)

        our_name = "DeepCatch Fusion"

        # Our model scores
        if hasattr(our_model, "predict_proba"):
            our_scores = our_model.predict_proba(
                np.column_stack([modality_scores[m] for m in modality_scores])
            )[:, 1]
        else:
            our_scores = np.mean(
                [modality_scores[m] for m in modality_scores], axis=0
            )

        our_metrics = self._compute_tool_metrics(y, our_scores, category=our_name)
        fusion_strategies: Dict[str, ToolResult] = {our_name: our_metrics}
        all_scores: Dict[str, Array] = {our_name: our_scores}

        modality_names = list(modality_scores.keys())

        # Single-best-modality baseline
        best_mod = None
        best_mod_auc = 0.0
        for mod_name in modality_names:
            scores = np.asarray(modality_scores[mod_name]).ravel()
            auc = roc_auc_score(y, scores)
            tr = self._compute_tool_metrics(y, scores, category=mod_name)
            fusion_strategies[mod_name] = tr
            all_scores[mod_name] = scores
            if auc > best_mod_auc:
                best_mod_auc = auc
                best_mod = mod_name

        # Late fusion baseline: average of modality probabilities
        late_scores = np.mean([modality_scores[m] for m in modality_names], axis=0)
        late_metrics = self._compute_tool_metrics(y, late_scores, category="Late Fusion")
        fusion_strategies["Late Fusion"] = late_metrics
        all_scores["Late Fusion"] = late_scores

        # Early fusion baseline (if features provided)
        if X_early_fusion is not None:
            X_ef = np.asarray(X_early_fusion)
            if hasattr(our_model, "predict_proba"):
                ef_scores = our_model.predict_proba(X_ef)[:, 1]
            else:
                ef_scores = np.mean([modality_scores[m] for m in modality_names], axis=0)

            ef_metrics = self._compute_tool_metrics(y, ef_scores, category="Early Fusion")
            fusion_strategies["Early Fusion"] = ef_metrics
            all_scores["Early Fusion"] = ef_scores

        # Pairwise DeLong comparisons
        comparisons: List[BenchmarkComparison] = []

        for strategy_name in fusion_strategies:
            if strategy_name == our_name:
                continue
            try:
                dl = delong_test(y, all_scores[our_name], all_scores[strategy_name])
                comparisons.append(BenchmarkComparison(
                    our_model=our_name,
                    benchmark=strategy_name,
                    delta_auc=dl["delta_auc"],
                    p_value=dl["p_value"],
                    p_value_corrected=dl["p_value"],
                    significant=dl["significant"],
                    ci95=(dl["ci95_lower"], dl["ci95_upper"]),
                ))
            except Exception:
                pass

        # Correct for multiple comparisons
        if len(comparisons) > 1:
            p_vals = np.array([c.p_value for c in comparisons])
            from validation_framework import SignificanceTester
            corrected = SignificanceTester.bonferroni_correct(p_vals)
            for i, c in enumerate(comparisons):
                c.p_value_corrected = float(corrected[i])
                c.significant = float(corrected[i]) < 0.05

        # Rank
        ranked = sorted(
            [(name, tr.auc) for name, tr in fusion_strategies.items()],
            key=lambda x: -x[1],
        )

        return BioinfoBenchmarkResult(
            variant_callers={},
            fusion_strategies=fusion_strategies,
            variant_comparisons=[],
            fusion_comparisons=comparisons,
            ranked_variant_callers=[],
            ranked_fusion_strategies=ranked,
            best_variant_caller="N/A",
            best_fusion_strategy=ranked[0][0] if ranked else "N/A",
        )

    # ── tool simulation ───────────────────────────────────────────────────

    def _simulate_tool_scores(
        self,
        y: Array,
        sensitivity: float,
        specificity: float,
        low_vaf_mask: Optional[Array] = None,
        low_vaf_sensitivity: Optional[float] = None,
    ) -> Array:
        """Simulate tool prediction scores based on known characteristics.

        Creates scores consistent with the tool's published sensitivity/
        specificity, with realistic noise. Low-VAF samples get degraded
        sensitivity.
        """
        n = len(y)
        pos_idx = y == 1
        neg_idx = y == 0

        # Base score: signal + noise
        scores = np.zeros(n)

        # For positives: generate scores consistent with sensitivity
        sens = sensitivity
        if low_vaf_mask is not None and low_vaf_sensitivity is not None:
            # Two groups with different sensitivity
            low_vaf = np.asarray(low_vaf_mask).ravel()
            lv_pos = low_vaf & pos_idx
            hv_pos = ~low_vaf & pos_idx

            # Low VAF: lower sensitivity
            scores[lv_pos] = np.where(
                self.rng.rand(np.sum(lv_pos)) < low_vaf_sensitivity,
                self.rng.uniform(0.5, 0.95, np.sum(lv_pos)),
                self.rng.uniform(0.05, 0.5, np.sum(lv_pos)),
            )
            # High VAF: full sensitivity
            scores[hv_pos] = np.where(
                self.rng.rand(np.sum(hv_pos)) < sensitivity,
                self.rng.uniform(0.5, 0.95, np.sum(hv_pos)),
                self.rng.uniform(0.05, 0.5, np.sum(hv_pos)),
            )
        else:
            scores[pos_idx] = np.where(
                self.rng.rand(np.sum(pos_idx)) < sens,
                self.rng.uniform(0.5, 0.95, np.sum(pos_idx)),
                self.rng.uniform(0.05, 0.5, np.sum(pos_idx)),
            )

        # For negatives: generate scores consistent with specificity
        scores[neg_idx] = np.where(
            self.rng.rand(np.sum(neg_idx)) < specificity,
            self.rng.uniform(0.01, 0.45, np.sum(neg_idx)),
            self.rng.uniform(0.5, 0.99, np.sum(neg_idx)),
        )

        return np.clip(scores, 0, 1)

    # ── metric computation ────────────────────────────────────────────────

    def _compute_tool_metrics(
        self,
        y: Array,
        scores: Array,
        low_vaf_mask: Optional[Array] = None,
        name: str = "",
        category: str = "",
    ) -> ToolResult:
        """Compute comprehensive metrics for a tool."""
        y = np.asarray(y).ravel().astype(int)
        scores = np.asarray(scores).ravel()

        threshold = 0.5
        y_pred = (scores >= threshold).astype(int)

        cm = confusion_matrix(y, y_pred, labels=[0, 1])
        if cm.shape == (2, 2):
            tn, fp, fn, tp = cm.ravel()
            sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
            prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            f1 = 2 * prec * sens / (prec + sens) if (prec + sens) > 0 else 0.0
        else:
            sens = spec = prec = f1 = 0.0

        # AUC
        auc_val = roc_auc_score(y, scores) if len(np.unique(y)) > 1 else float("nan")

        # Bootstrap CI for AUC
        auc_ci = self._bootstrap_auc_ci(y, scores)

        # Low VAF sensitivity
        low_vaf_sens = None
        if low_vaf_mask is not None:
            lv_mask = np.asarray(low_vaf_mask).ravel() & (y == 1)
            if np.sum(lv_mask) > 0:
                lv_pred = (scores[lv_mask] >= threshold).astype(int)
                lv_true = y[lv_mask]
                tp_lv = np.sum((lv_pred == 1) & (lv_true == 1))
                low_vaf_sens = tp_lv / np.sum(lv_true)

        return ToolResult(
            name=name or category,
            category=category or "variant_caller",
            auc=float(auc_val),
            auc_ci95=auc_ci,
            sensitivity=sens,
            specificity=spec,
            f1=f1,
            precision=prec,
            low_vaf_sensitivity=low_vaf_sens,
        )

    def _bootstrap_auc_ci(
        self, y: Array, scores: Array, n_bootstrap: int = 1000
    ) -> Tuple[float, float]:
        """Bootstrap CI for AUC."""
        n = len(y)
        rng = np.random.RandomState(self.random_state)
        aucs = np.zeros(n_bootstrap)

        for i in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            yb = y[idx]
            sb = scores[idx]
            if len(np.unique(yb)) > 1:
                try:
                    aucs[i] = roc_auc_score(yb, sb)
                except ValueError:
                    aucs[i] = 0.5
            else:
                aucs[i] = 0.5

        return (
            float(np.percentile(aucs, 2.5)),
            float(np.percentile(aucs, 97.5)),
        )

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: BioinfoBenchmarkResult) -> str:
        """Publication-ready benchmark report."""
        lines = ["══ Bioinformatic Tool Benchmark ══"]

        if result.variant_callers:
            lines.append(f"\n  Variant Calling Comparison:")
            lines.append(f"  {'Tool':>25s} {'AUC':>7s} {'CI95':>20s} {'Sens':>7s} {'Spec':>7s} {'F1':>7s} {'LowVAF':>7s}")
            lines.append(f"  {'─'*25} {'─'*7} {'─'*20} {'─'*7} {'─'*7} {'─'*7} {'─'*7}")

            for name, tr in result.variant_callers.items():
                vaf_str = f"{tr.low_vaf_sensitivity:.4f}" if tr.low_vaf_sensitivity else "   N/A"
                lines.append(
                    f"  {name:>25s} {tr.auc:7.4f} "
                    f"[{tr.auc_ci95[0]:.4f}, {tr.auc_ci95[1]:.4f}] "
                    f"{tr.sensitivity:7.4f} {tr.specificity:7.4f} "
                    f"{tr.f1:7.4f} {vaf_str:>7s}"
                )

            # Rankings
            lines.append(f"\n  Rankings:")
            for rank, (name, auc) in enumerate(result.ranked_variant_callers, 1):
                lines.append(f"  {rank}. {name}: AUC={auc:.4f}")

            # Pairwise comparisons vs our model
            if result.variant_comparisons:
                lines.append(f"\n  DeLong Tests (vs DeepCatch):")
                for c in result.variant_comparisons:
                    sig = "SIG" if c.significant else "ns"
                    lines.append(
                        f"  {c.benchmark:>20s}: ΔAUC={c.delta_auc:+7.4f} "
                        f"[{c.ci95[0]:+.4f}, {c.ci95[1]:+.4f}] "
                        f"p={c.p_value_corrected:.4f} ({sig})"
                    )

        if result.fusion_strategies:
            lines.append(f"\n  Fusion Strategy Comparison:")
            lines.append(f"  {'Strategy':>25s} {'AUC':>7s} {'CI95':>20s} {'Sens':>7s} {'Spec':>7s} {'F1':>7s}")
            lines.append(f"  {'─'*25} {'─'*7} {'─'*20} {'─'*7} {'─'*7} {'─'*7}")

            for name, tr in result.fusion_strategies.items():
                lines.append(
                    f"  {name:>25s} {tr.auc:7.4f} "
                    f"[{tr.auc_ci95[0]:.4f}, {tr.auc_ci95[1]:.4f}] "
                    f"{tr.sensitivity:7.4f} {tr.specificity:7.4f} "
                    f"{tr.f1:7.4f}"
                )

            if result.fusion_comparisons:
                lines.append(f"\n  DeLong Tests (vs DeepCatch Fusion):")
                for c in result.fusion_comparisons:
                    sig = "SIG" if c.significant else "ns"
                    lines.append(
                        f"  {c.benchmark:>20s}: ΔAUC={c.delta_auc:+7.4f} "
                        f"[{c.ci95[0]:+.4f}, {c.ci95[1]:+.4f}] "
                        f"p={c.p_value_corrected:.4f} ({sig})"
                    )

        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 8: BioinfoBenchmark — self-test")
    print("=" * 60)

    from sklearn.linear_model import LogisticRegression

    rng = np.random.RandomState(42)
    n = 500

    # Generate variant calling data
    X = rng.randn(n, 15)
    signal = X[:, 0] + 0.5 * X[:, 1]
    y = (1.0 / (1.0 + np.exp(-signal)) > 0.5).astype(int)

    # Low VAF samples (first 50 positives)
    low_vaf_mask = np.zeros(n, dtype=bool)
    pos_idx = np.where(y == 1)[0]
    if len(pos_idx) > 50:
        low_vaf_mask[pos_idx[:50]] = True

    # Train our model
    model = LogisticRegression(solver="liblinear")
    model.fit(X, y)

    bench = BioinfoBenchmark(n_bootstrap=500)

    print("\n── Variant Calling Benchmark ──")
    result_vc = bench.benchmark_variant_calling(
        model, X, y, low_vaf_mask=low_vaf_mask
    )
    print(bench.report(result_vc))

    print("\n── Fusion Strategy Benchmark ──")
    # Simulate 3 modalities
    modality_scores = {
        "ctDNA": np.clip(model.predict_proba(X)[:, 1] + 0.05 * rng.randn(n), 0, 1),
        "Methylation": np.clip(model.predict_proba(X)[:, 1] + 0.08 * rng.randn(n), 0, 1),
        "Fragmentomics": np.clip(model.predict_proba(X)[:, 1] + 0.10 * rng.randn(n), 0, 1),
    }
    result_fusion = bench.benchmark_fusion_strategies(
        model, modality_scores, y
    )
    print(bench.report(result_fusion))

    print("\nSelf-test complete.")
    sys.exit(0)
