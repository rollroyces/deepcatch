#!/usr/bin/env python3
"""
MODULE 7: Confounder Robustness Suite
======================================

Reference: Lipsitch, Tchetgen & Cohen (2010), "Negative Controls: A Tool for
Detecting Confounding and Bias in Observational Studies", Epidemiology
21:383-388.

Tests model performance under realistic confounders known to affect ctDNA-based
cancer screening:

1. CHIP (Clonal Hematopoiesis of Indeterminate Potential)
   → Age-related somatic mutations in blood, not cancer.
   → Simulated as background "mutation noise" scaling with age.

2. Batch Effects
   → Systematic technical variation between sequencing runs.
   → Simulated as additive/multiplicative shifts per batch.

3. Inflammatory Conditions
   → Transient cfDNA spikes from non-malignant sources.
   → Simulated as temporary elevation of background signal.

4. Variable Blood Volume
   → ±50% variation in input plasma volume.
   → Simulated as scaling of observed fragment count.

5. Sequencing Depth Variation
   → 10K-100K× coverage variation between samples.
   → Simulated as Poisson thinning of features.

6. Library Preparation Variability
   → GC bias and fragment size selection biases.
   → Simulated as systematic feature shifts based on GC content.

Reports: performance degradation vs confounder strength.
Identifies which confounders are most damaging.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray
ModelFactory = Callable[[], Any]


# ── confounder definitions ──────────────────────────────────────────────────


@dataclass
class ConfounderConfig:
    """Configuration for a single confounder scenario."""
    name: str
    description: str
    # Strength levels to test
    strengths: List[float]
    strength_labels: List[str]
    # Whether this confounder is expected to affect all samples
    affects_all: bool = True
    # Whether the effect is stronger on positives or negatives
    bias_toward: str = "both"  # "positive", "negative", "both"


@dataclass
class ConfounderResult:
    """Results for one confounder at one strength level."""
    confounder: str
    strength: float
    strength_label: str
    baseline_auc: float
    confounded_auc: float
    auc_degradation: float  # baseline - confounded
    auc_degradation_pct: float  # percentage degradation
    # Bootstrap CI for degradation
    degradation_ci95: Tuple[float, float]
    # Whether degradation is statistically significant
    p_value: float
    significant: bool


@dataclass
class ConfounderSuiteResult:
    """Full confounder robustness results."""
    individual_results: List[ConfounderResult]
    summary_table: Dict[str, List[ConfounderResult]]  # confounder → results
    # Rank confounders by maximum degradation
    ranked_by_impact: List[Tuple[str, float, float]]  # (name, max_degradation, pct)
    n_confounders: int
    n_samples: int


# ── confounder robustness tester ────────────────────────────────────────────


class ConfounderRobustnessTester:
    """Tests model robustness to realistic clinical confounders.

    Usage:
        tester = ConfounderRobustnessTester()
        results = tester.test_all(model_factory, X, y,
                                   patient_ages=ages,
                                   gc_content=gc_content)
    """

    def __init__(
        self,
        n_bootstrap: int = 500,
        random_state: int = 42,
        verbose: bool = False,
    ):
        """Args:
            n_bootstrap: Bootstrap replicates for degradation CIs.
            random_state: Seed.
            verbose: Print progress.
        """
        self.n_bootstrap = n_bootstrap
        self.random_state = random_state
        self.verbose = verbose
        self.rng = np.random.RandomState(random_state)

    # ── main testing API ──────────────────────────────────────────────────

    def test_all(
        self,
        model_factory: ModelFactory,
        X: Array,
        y: Array,
        patient_ages: Optional[Array] = None,
        batch_labels: Optional[Array] = None,
        gc_content: Optional[Array] = None,
        fragment_sizes: Optional[Array] = None,
        custom_confounders: Optional[List[ConfounderConfig]] = None,
    ) -> ConfounderSuiteResult:
        """Run all confounder robustness tests.

        Args:
            model_factory: Zero-arg callable returning a trained model (or
                a train function: (X, y) → model with .predict_proba).
                If the factory returns an untrained model, it will be fit
                on X, y before each test.
            X: Feature matrix, shape (n_samples, n_features).
            y: Ground-truth labels, shape (n_samples,).
            patient_ages: Optional ages for CHIP simulation.
            batch_labels: Optional batch assignments.
            gc_content: Optional per-feature GC content for library bias.
            fragment_sizes: Optional fragment size profiles.
            custom_confounders: Additional confounders to test.

        Returns:
            ConfounderSuiteResult with full degradation analysis.
        """
        X = np.asarray(X)
        y = np.asarray(y).ravel().astype(int)

        # Build default confounder list
        confounders = self._build_default_confounders(
            X, patient_ages, batch_labels, gc_content, fragment_sizes
        )
        if custom_confounders:
            confounders.extend(custom_confounders)

        # Get baseline performance
        baseline_auc = self._evaluate_auc(model_factory, X, y)

        if self.verbose:
            print(f"Baseline AUC: {baseline_auc:.4f}")
            print(f"Testing {len(confounders)} confounders...")

        all_results: List[ConfounderResult] = []
        summary: Dict[str, List[ConfounderResult]] = {}

        for conf in confounders:
            conf_results: List[ConfounderResult] = []

            for strength, label in zip(conf.strengths, conf.strength_labels):
                if self.verbose:
                    print(f"  {conf.name} @ {label}...", end=" ")

                # Apply confounder
                X_perturbed = self._apply_confounder(
                    conf, X, y, strength, patient_ages, batch_labels,
                    gc_content, fragment_sizes
                )

                # Evaluate
                conf_auc = self._evaluate_auc(model_factory, X_perturbed, y)
                degradation = baseline_auc - conf_auc
                degradation_pct = (
                    100 * degradation / baseline_auc if baseline_auc > 0 else 0.0
                )

                # Bootstrap CI for degradation
                deg_ci, p_val = self._bootstrap_degradation(
                    model_factory, X, X_perturbed, y, baseline_auc, conf_auc
                )

                cr = ConfounderResult(
                    confounder=conf.name,
                    strength=strength,
                    strength_label=label,
                    baseline_auc=baseline_auc,
                    confounded_auc=conf_auc,
                    auc_degradation=degradation,
                    auc_degradation_pct=degradation_pct,
                    degradation_ci95=deg_ci,
                    p_value=p_val,
                    significant=p_val < 0.05,
                )
                conf_results.append(cr)
                all_results.append(cr)

                if self.verbose:
                    sig = "*" if cr.significant else " "
                    print(f"AUC={conf_auc:.4f} Δ={degradation:+.4f}{sig}")

            summary[conf.name] = conf_results

        # Rank by max degradation
        ranked: List[Tuple[str, float, float]] = []
        for name, results in summary.items():
            max_deg = max(r.auc_degradation for r in results)
            max_pct = max(r.auc_degradation_pct for r in results)
            ranked.append((name, max_deg, max_pct))
        ranked.sort(key=lambda x: -x[1])

        return ConfounderSuiteResult(
            individual_results=all_results,
            summary_table=summary,
            ranked_by_impact=ranked,
            n_confounders=len(confounders),
            n_samples=len(y),
        )

    # ── confounder application ────────────────────────────────────────────

    def _build_default_confounders(
        self,
        X: Array,
        patient_ages: Optional[Array],
        batch_labels: Optional[Array],
        gc_content: Optional[Array],
        fragment_sizes: Optional[Array],
    ) -> List[ConfounderConfig]:
        """Build the default set of bioinformatic confounders."""
        confounders: List[ConfounderConfig] = []

        # 1. CHIP (Clonal Hematopoiesis)
        if patient_ages is not None:
            confounders.append(ConfounderConfig(
                name="CHIP (Clonal Hematopoiesis)",
                description="Age-related background mutations in blood",
                strengths=[0.0, 0.5, 1.0, 1.5, 2.0],
                strength_labels=["None", "Mild", "Moderate", "Strong", "Extreme"],
                bias_toward="negative",  # CHIP primarily affects healthy individuals
            ))
        else:
            # Use age-independent CHIP simulation
            confounders.append(ConfounderConfig(
                name="CHIP (simulated, age-independent)",
                description="Random background mutation noise",
                strengths=[0.0, 0.3, 0.6, 1.0],
                strength_labels=["None", "Low", "Medium", "High"],
                bias_toward="negative",
            ))

        # 2. Batch effects
        if batch_labels is not None:
            confounders.append(ConfounderConfig(
                name="Batch Effects",
                description="Systematic shifts between sequencing runs",
                strengths=[0.0, 0.2, 0.5, 1.0],
                strength_labels=["None", "Small", "Medium", "Large"],
            ))
        else:
            confounders.append(ConfounderConfig(
                name="Batch Effects (simulated)",
                description="Per-sample random shifts",
                strengths=[0.0, 0.2, 0.5, 1.0],
                strength_labels=["None", "Small", "Medium", "Large"],
            ))

        # 3. Inflammatory conditions
        confounders.append(ConfounderConfig(
            name="Inflammatory Conditions",
            description="Transient cfDNA spikes from non-malignant sources",
            strengths=[0.0, 0.2, 0.5, 1.0],
            strength_labels=["None", "Mild", "Moderate", "Severe"],
            affects_all=False,  # Only affects subset of patients
            bias_toward="negative",  # Inflammatory spikes mimic cancer signal
        ))

        # 4. Variable blood volume (±50%)
        confounders.append(ConfounderConfig(
            name="Blood Volume Variation",
            description="Scaling of observed fragment counts",
            strengths=[1.0, 0.75, 0.5, 0.25],
            strength_labels=["100%", "75%", "50%", "25%"],
            bias_toward="both",
        ))

        # 5. Sequencing depth variation
        confounders.append(ConfounderConfig(
            name="Sequencing Depth",
            description="Coverage variation (Poisson thinning)",
            strengths=[1.0, 0.7, 0.4, 0.1],
            strength_labels=["100%", "70%", "40%", "10%"],
            bias_toward="both",
        ))

        # 6. Library preparation GC bias
        if gc_content is not None:
            confounders.append(ConfounderConfig(
                name="Library GC Bias",
                description="GC-content-dependent amplification bias",
                strengths=[0.0, 0.05, 0.10, 0.20],
                strength_labels=["None", "Mild", "Moderate", "Severe"],
            ))
        else:
            confounders.append(ConfounderConfig(
                name="Library GC Bias (simulated)",
                description="Random per-feature amplification bias",
                strengths=[0.0, 0.05, 0.10, 0.20],
                strength_labels=["None", "Mild", "Moderate", "Severe"],
            ))

        return confounders

    def _apply_confounder(
        self,
        conf: ConfounderConfig,
        X: Array,
        y: Array,
        strength: float,
        patient_ages: Optional[Array],
        batch_labels: Optional[Array],
        gc_content: Optional[Array],
        fragment_sizes: Optional[Array],
    ) -> Array:
        """Apply a confounder to the feature matrix.

        Returns:
            Perturbed copy of X.
        """
        X_pert = X.copy().astype(float)
        n, p = X_pert.shape

        if "CHIP" in conf.name:
            # Age-dependent background mutation noise
            if patient_ages is not None and len(patient_ages) == n:
                # Normalize ages to [0, 1] for scaling
                age_norm = patient_ages / 80.0
                noise = strength * age_norm[:, np.newaxis] * self.rng.randn(n, p) * 0.1
                # More noise on negatives (healthy patients with CHIP mimic cancer)
                neg_mask = y == 0
                X_pert[neg_mask] += noise[neg_mask] * 1.5
                X_pert[y == 1] += noise[y == 1] * 0.5
            else:
                # Flat noise, biased toward negatives
                noise = strength * 0.1 * self.rng.randn(n, p)
                X_pert[y == 0] += noise[y == 0] * 2.0
                X_pert[y == 1] += noise[y == 1] * 0.5

        elif "Batch" in conf.name:
            if batch_labels is not None and len(batch_labels) == n:
                batches = np.unique(batch_labels)
                for batch in batches:
                    mask = batch_labels == batch
                    shift = strength * self.rng.randn(p) * np.std(X, axis=0)
                    X_pert[mask] += shift[np.newaxis, :]
            else:
                # Random per-sample shifts
                X_pert += strength * 0.3 * self.rng.randn(n, p) * np.std(X, axis=0)[np.newaxis, :]

        elif "Inflammatory" in conf.name:
            # Only affect a subset of patients (healthy ones get cfDNA spikes)
            n_affected = int(0.3 * n)
            affected_idx = self.rng.choice(
                np.where(y == 0)[0], size=min(n_affected, np.sum(y == 0)), replace=False
            )
            if len(affected_idx) > 0:
                # Add signal across ~20% of features to simulate cfDNA spike
                n_feat_affected = max(1, int(0.2 * p))
                feat_idx = self.rng.choice(p, size=n_feat_affected, replace=False)
                X_pert[np.ix_(affected_idx, feat_idx)] += strength * 1.5 * np.std(X[:, feat_idx])

        elif "Blood Volume" in conf.name:
            # Scale ALL features by volume factor (strength IS the volume fraction)
            X_pert = X_pert * strength

        elif "Sequencing Depth" in conf.name:
            # Poisson thinning: keep only a fraction of "counts"
            mask = self.rng.rand(n, p) < strength
            # Zero out features that fall below detection threshold
            X_pert[~mask] = 0.0

        elif "GC Bias" in conf.name:
            if gc_content is not None and len(gc_content) == p:
                # Feature-dependent amplification bias based on GC content
                gc_norm = (gc_content - np.mean(gc_content)) / (np.std(gc_content) + 1e-10)
                bias = 1.0 + strength * gc_norm[np.newaxis, :]
                X_pert = X_pert * bias
            else:
                # Simulate GC bias with random per-feature scaling
                bias = 1.0 + strength * self.rng.randn(p) * 0.5
                X_pert = X_pert * bias[np.newaxis, :]

        return X_pert

    # ── evaluation ────────────────────────────────────────────────────────

    def _evaluate_auc(
        self,
        model_factory: ModelFactory,
        X: Array,
        y: Array,
    ) -> float:
        """Evaluate AUC with internal 3-fold CV for robustness."""
        from sklearn.model_selection import StratifiedKFold

        y = np.asarray(y).ravel().astype(int)

        if len(np.unique(y)) < 2:
            return float("nan")

        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=self.random_state)
        aucs: List[float] = []

        for train_idx, test_idx in skf.split(X, y):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                continue

            model = model_factory()
            model.fit(X_train, y_train)

            if hasattr(model, "predict_proba"):
                y_score = model.predict_proba(X_test)[:, 1]
            elif hasattr(model, "decision_function"):
                y_score = model.decision_function(X_test)
            else:
                y_score = model.predict(X_test).astype(float)

            try:
                aucs.append(roc_auc_score(y_test, y_score))
            except ValueError:
                continue

        return float(np.mean(aucs)) if aucs else float("nan")

    def _bootstrap_degradation(
        self,
        model_factory: ModelFactory,
        X_original: Array,
        X_perturbed: Array,
        y: Array,
        baseline_auc: float,
        confounded_auc: float,
    ) -> Tuple[Tuple[float, float], float]:
        """Bootstrap CI for AUC degradation and p-value."""
        n = len(y)
        rng = np.random.RandomState(self.random_state)
        degradations = np.zeros(self.n_bootstrap)

        for i in range(self.n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            y_boot = y[idx]
            if len(np.unique(y_boot)) < 2:
                degradations[i] = 0.0
                continue

            # Train on bootstrapped original data
            model = model_factory()
            model.fit(X_original[idx], y_boot)

            if hasattr(model, "predict_proba"):
                score_orig = model.predict_proba(X_original[idx])[:, 1]
                score_pert = model.predict_proba(X_perturbed[idx])[:, 1]
            else:
                score_orig = model.decision_function(X_original[idx])
                score_pert = model.decision_function(X_perturbed[idx])

            try:
                auc_o = roc_auc_score(y_boot, score_orig)
                auc_p = roc_auc_score(y_boot, score_pert)
                degradations[i] = auc_o - auc_p
            except ValueError:
                degradations[i] = 0.0

        ci = (
            float(np.percentile(degradations, 2.5)),
            float(np.percentile(degradations, 97.5)),
        )

        # Two-sided p-value: fraction of bootstrap degradations ≤ 0
        obs_deg = baseline_auc - confounded_auc
        p_val = float(np.mean(degradations <= 0)) if obs_deg > 0 else 1.0

        return ci, p_val

    # ── reporting ─────────────────────────────────────────────────────────

    @staticmethod
    def report(result: ConfounderSuiteResult) -> str:
        """Publication-ready confounder robustness report."""
        lines = [
            "══ Confounder Robustness Analysis ══",
            f"  Samples: {result.n_samples}",
            f"  Confounders tested: {result.n_confounders}",
            "",
            "  Impact Ranking (by max AUC degradation):",
        ]

        for rank, (name, deg, pct) in enumerate(result.ranked_by_impact, 1):
            severity = (
                "🔴 Critical" if pct > 5
                else "🟡 Moderate" if pct > 2
                else "🟢 Minor" if pct > 0.5
                else "⚪ Negligible"
            )
            lines.append(f"  {rank}. {name}: ΔAUC={deg:.4f} ({pct:.1f}%) — {severity}")

        # Detailed per-confounder breakdown
        for name, results in result.summary_table.items():
            lines.append(f"\n  ── {name} ──")
            lines.append(f"  {'Strength':>12s} {'AUC':>8s} {'Δ AUC':>8s} {'Δ %':>7s} {'CI95':>20s} {'p':>7s}")
            lines.append(f"  {'─'*12} {'─'*8} {'─'*8} {'─'*7} {'─'*20} {'─'*7}")

            for r in results:
                sig = "*" if r.significant else " "
                lines.append(
                    f"  {r.strength_label:>12s} {r.confounded_auc:8.4f} "
                    f"{r.auc_degradation:+8.4f} {r.auc_degradation_pct:+6.1f}% "
                    f"[{r.degradation_ci95[0]:.4f}, {r.degradation_ci95[1]:.4f}] "
                    f"{r.p_value:7.4f}{sig}"
                )

        # Summary
        critical = sum(1 for r in result.individual_results if r.auc_degradation_pct > 5)
        moderate = sum(1 for r in result.individual_results if 2 < r.auc_degradation_pct <= 5)
        minor = sum(1 for r in result.individual_results if 0.5 < r.auc_degradation_pct <= 2)

        lines.append(f"\n  Summary: {critical} critical, {moderate} moderate, {minor} minor degradations")
        lines.append(f"  * = statistically significant (p < 0.05)")

        return "\n".join(lines)


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 7: ConfounderRobustnessTester — self-test")
    print("=" * 60)

    from sklearn.linear_model import LogisticRegression

    rng = np.random.RandomState(42)
    n = 400
    p = 20

    # Generate data with real signal
    X = rng.randn(n, p)
    signal = X[:, 0] + 0.5 * X[:, 1] - 0.3 * X[:, 2]
    y = (1.0 / (1.0 + np.exp(-signal)) > 0.5).astype(int)

    # Patient ages for CHIP
    ages = rng.randint(30, 80, size=n).astype(float)

    # Batch labels
    batch_labels = np.array(["Batch_A"] * 200 + ["Batch_B"] * 200)

    # GC content
    gc_content = 0.3 + 0.4 * rng.rand(p)

    tester = ConfounderRobustnessTester(verbose=True, n_bootstrap=200)

    def make_model():
        return LogisticRegression(solver="liblinear", C=1.0)

    result = tester.test_all(
        make_model, X, y,
        patient_ages=ages,
        batch_labels=batch_labels,
        gc_content=gc_content,
    )

    print("\n" + tester.report(result))
    print("\nSelf-test complete.")
    sys.exit(0)
