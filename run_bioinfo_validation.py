#!/usr/bin/env python3
"""
=====================================================================
DeepCatch Bioinformatics Validation Suite — Unified Runner
=====================================================================

Runs ALL validation modules in sequence:
  [1/10]  Nested Cross-Validation
  [2/10]  Permutation Testing
  [3/10]  Calibration Analysis
  [4/10]  Decision Curve Analysis
  [5/10]  DeLong Statistical Tests
  [6/10]  Stratified Performance Analysis
  [7/10]  Confounder Robustness Suite
  [8/10]  Bioinformatic Tool Benchmark
  [9/10]  Sample Size & Power Analysis
  [10/10] Reproducibility Verification

Outputs a comprehensive report:
  results/BIOINFORMATICS_VALIDATION_REPORT.md

Usage:
  python run_bioinfo_validation.py
  python run_bioinfo_validation.py --quick    # Fast mode (fewer iterations)
  python run_bioinfo_validation.py --data PATH  # Custom data directory
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

# ── Ensure we can import from the project root ──────────────────────────
sys.path.insert(0, str(Path(__file__).parent))

from validation_framework import (
    DataSplitter,
    BootstrapCI,
    ThresholdCalibrator,
    KFoldValidator,
    SignificanceTester,
    ClaimValidator,
)

from validation.nested_cv import NestedCrossValidator
from validation.permutation_test import PermutationTester
from validation.calibration import CalibrationAnalyzer
from validation.decision_curve import DecisionCurveAnalyzer
from validation.delong_test import delong_test, delong_test_multi, report as delong_report
from validation.stratified import StratifiedAnalyzer
from validation.confounders import ConfounderRobustnessTester
from validation.bioinfo_benchmark import BioinfoBenchmark
from validation.power_analysis import PowerAnalyzer

# ── constants ───────────────────────────────────────────────────────────────
RESULTS_DIR = Path(__file__).parent / "results"
REPORT_PATH = RESULTS_DIR / "BIOINFORMATICS_VALIDATION_REPORT.md"
SEED_REGISTRY_PATH = Path(__file__).parent / "reproducibility" / "seed_registry.json"

Array = np.ndarray


# ═════════════════════════════════════════════════════════════════════════════
# Synthetic Data Generator (for demonstration/CI testing)
# ═════════════════════════════════════════════════════════════════════════════

def generate_synthetic_data(
    n_samples: int = 1000,
    n_features: int = 20,
    seed: int = 42,
) -> Tuple[Array, Array, Dict[str, Any]]:
    """Generate synthetic cancer screening data for validation testing.

    Creates realistic synthetic data with:
    - Signal-carrying features
    - Noise features
    - Stratified cancer types
    - Age distribution
    - Batch labels

    Used when real data is not available (e.g., CI/testing environments).
    """
    rng = np.random.RandomState(seed)

    X = np.zeros((n_samples, n_features))
    y = np.zeros(n_samples, dtype=int)

    # Feature 0-2: informative cancer signal
    signal = rng.randn(n_samples)

    # 30% prevalence
    pos_idx = rng.choice(n_samples, size=int(n_samples * 0.3), replace=False)
    y[pos_idx] = 1

    # Cancer patients get elevated signal features
    X[pos_idx, 0] = 1.5 + 0.5 * rng.randn(len(pos_idx))
    X[pos_idx, 1] = 0.8 + 0.6 * rng.randn(len(pos_idx))
    X[pos_idx, 2] = 0.3 + 0.4 * rng.randn(len(pos_idx))

    # Healthy patients
    neg_idx = np.where(y == 0)[0]
    X[neg_idx, 0] = 0.2 * rng.randn(len(neg_idx))
    X[neg_idx, 1] = 0.1 * rng.randn(len(neg_idx))
    X[neg_idx, 2] = -0.1 + 0.3 * rng.randn(len(neg_idx))

    # Feature 3-19: noise features
    X[:, 3:] = rng.randn(n_samples, n_features - 3) * 0.5

    # Metadata
    ages = rng.randint(30, 85, size=n_samples).astype(float)
    cancer_types = np.full(n_samples, "Healthy", dtype=object)
    cancer_types[pos_idx] = rng.choice(
        ["LUAD", "COADREAD", "BRCA", "PRAD", "LGG"],
        size=len(pos_idx), replace=True,
    )
    batch_labels = np.array(["Batch_A"] * (n_samples // 2) + ["Batch_B"] * (n_samples - n_samples // 2))

    # ctDNA fraction for positive patients
    ctdna_fractions = np.zeros(n_samples)
    ctdna_fractions[pos_idx] = 10 ** rng.uniform(-4, -0.5, size=len(pos_idx))

    # GC content (per feature)
    gc_content = 0.3 + 0.4 * rng.rand(n_features)

    # Low VAF mask (bottom 40% of positive ctdna fractions)
    low_vaf_mask = np.zeros(n_samples, dtype=bool)
    if np.sum(y == 1) > 0:
        pos_ctdna = ctdna_fractions[pos_idx]
        low_vaf_thresh = np.percentile(pos_ctdna, 40)
        low_vaf_pos = pos_idx[pos_ctdna <= low_vaf_thresh]
        low_vaf_mask[low_vaf_pos] = True

    metadata = {
        "ages": ages,
        "cancer_types": cancer_types,
        "batch_labels": batch_labels,
        "ctdna_fractions": ctdna_fractions,
        "gc_content": gc_content,
        "low_vaf_mask": low_vaf_mask,
        "n_features": n_features,
        "n_samples": n_samples,
        "prevalence": float(np.mean(y)),
    }

    return X, y, metadata


# ═════════════════════════════════════════════════════════════════════════════
# Module wrappers
# ═════════════════════════════════════════════════════════════════════════════

def run_module_1(X: Array, y: Array, quick: bool = False) -> str:
    """Nested Cross-Validation."""
    from sklearn.linear_model import LogisticRegression

    ncv = NestedCrossValidator(
        n_outer=3 if quick else 5,
        n_inner=2 if quick else 3,
        scoring="roc_auc",
        verbose=True,
    )
    param_grid = {"C": [0.01, 0.1, 1.0, 10.0]}

    result = ncv.validate(
        model_factory=lambda: LogisticRegression(solver="liblinear", max_iter=500),
        param_grid=param_grid,
        X=X, y=y,
    )
    return ncv.report(result, model_name="Logistic Regression (ctDNA)")


def run_module_2(X: Array, y: Array, quick: bool = False) -> str:
    """Permutation Testing."""
    from sklearn.linear_model import LogisticRegression

    tester = PermutationTester(
        n_permutations=100 if quick else 500,
        scoring="roc_auc",
        verbose=True,
    )
    result = tester.test(
        model_factory=lambda: LogisticRegression(solver="liblinear"),
        X=X, y=y,
        model_name="Logistic Regression (ctDNA)",
    )
    return tester.report(result)


def run_module_3(X: Array, y: Array, quick: bool = False) -> str:
    """Calibration Analysis."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split

    # Train model and get scores
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    model = LogisticRegression(solver="liblinear", max_iter=500)
    model.fit(X_train, y_train)
    y_score = model.predict_proba(X_val)[:, 1]

    analyzer = CalibrationAnalyzer(n_bins=10)
    result = analyzer.analyze_with_recalibration(
        y_val, y_score,
        X_val=y_score.reshape(-1, 1), y_val=y_val,
        X_test=y_score.reshape(-1, 1), y_test=y_val,
    )
    return analyzer.report(result)


def run_module_4(X: Array, y: Array, quick: bool = False) -> str:
    """Decision Curve Analysis."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split

    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    model = LogisticRegression(solver="liblinear", max_iter=500)
    model.fit(X_train, y_train)
    y_score = model.predict_proba(X_val)[:, 1]

    dca = DecisionCurveAnalyzer(
        thresholds=np.linspace(0.01, 0.50, 50) if quick else np.linspace(0.01, 0.50, 100)
    )
    result = dca.analyze(y_val, y_score, model_name="ctDNA Screening Model")
    return dca.report(result)


def run_module_5(X: Array, y: Array, quick: bool = False) -> str:
    """DeLong Statistical Tests."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )

    model_a = LogisticRegression(solver="liblinear", max_iter=500)
    model_b = RandomForestClassifier(n_estimators=50, random_state=42)
    model_a.fit(X_train, y_train)
    model_b.fit(X_train, y_train)

    scores_a = model_a.predict_proba(X_test)[:, 1]
    scores_b = model_b.predict_proba(X_test)[:, 1]

    result = delong_test(y_test, scores_a, scores_b)
    return delong_report(result)


def run_module_6(X: Array, y: Array, metadata: Dict, quick: bool = False) -> str:
    """Stratified Performance Analysis."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    model = LogisticRegression(solver="liblinear", max_iter=500)
    model.fit(X_train, y_train)
    y_score_test = model.predict_proba(X_test)[:, 1]

    analyzer = StratifiedAnalyzer(min_stratum_size=20, n_bootstrap=200 if quick else 500)
    result = analyzer.analyze(
        y_test, y_score_test,
        strata_labels=metadata["cancer_types"][len(y_train):],
    )
    return analyzer.report(result)


def run_module_7(X: Array, y: Array, metadata: Dict, quick: bool = False) -> str:
    """Confounder Robustness Suite."""
    from sklearn.linear_model import LogisticRegression

    tester = ConfounderRobustnessTester(
        n_bootstrap=100 if quick else 300,
        verbose=True,
    )
    result = tester.test_all(
        model_factory=lambda: LogisticRegression(solver="liblinear", max_iter=200),
        X=X, y=y,
        patient_ages=metadata.get("ages"),
        batch_labels=metadata.get("batch_labels"),
        gc_content=metadata.get("gc_content"),
    )
    return tester.report(result)


def run_module_8(X: Array, y: Array, metadata: Dict, quick: bool = False) -> str:
    """Bioinformatic Tool Benchmark."""
    from sklearn.linear_model import LogisticRegression

    model = LogisticRegression(solver="liblinear", max_iter=500)
    model.fit(X, y)

    bench = BioinfoBenchmark(n_bootstrap=200 if quick else 500)
    result = bench.benchmark_variant_calling(
        model, X, y,
        low_vaf_mask=metadata.get("low_vaf_mask"),
    )
    return bench.report(result)


def run_module_9(quick: bool = False) -> str:
    """Sample Size & Power Analysis."""
    pa = PowerAnalyzer(alpha=0.05, target_power=0.80)

    experiments = [
        {
            "name": "CET Detection (AUC≈0.95)",
            "metric": "AUC",
            "observed_effect": 0.95,
            "observed_se": 0.02,
            "n_samples": 200,
            "n_pos": 60,
            "n_neg": 140,
        },
        {
            "name": "GNN Fusion (AUC≈0.75)",
            "metric": "AUC",
            "observed_effect": 0.75,
            "observed_se": 0.04,
            "n_samples": 600,
            "n_pos": 180,
            "n_neg": 420,
        },
        {
            "name": "Bayesian Caller (AUC≈0.88)",
            "metric": "AUC",
            "observed_effect": 0.88,
            "observed_se": 0.03,
            "n_samples": 500,
            "n_pos": 150,
            "n_neg": 350,
        },
        {
            "name": "Contrastive Learner (AUC≈0.68)",
            "metric": "AUC",
            "observed_effect": 0.68,
            "observed_se": 0.05,
            "n_samples": 100,
            "n_pos": 30,
            "n_neg": 70,
        },
        {
            "name": "Multi-Modal Ensemble (AUC≈0.85)",
            "metric": "AUC",
            "observed_effect": 0.85,
            "observed_se": 0.03,
            "n_samples": 800,
            "n_pos": 240,
            "n_neg": 560,
        },
    ]
    result = pa.analyze(experiments)
    return pa.report(result)


def run_module_10() -> str:
    """Reproducibility Verification."""
    lines = ["══ Reproducibility Verification ══", ""]

    # System info
    lines.append(f"  Python:     {platform.python_version()}")
    lines.append(f"  Platform:   {platform.platform()}")
    lines.append(f"  Timestamp:  {datetime.now(timezone.utc).isoformat()}")

    # Package versions
    lines.append(f"\n  Package Versions:")
    for pkg in ["numpy", "scipy", "sklearn", "pandas"]:
        try:
            mod = __import__(pkg)
            ver = getattr(mod, "__version__", "unknown")
            lines.append(f"    {pkg}: {ver}")
        except ImportError:
            lines.append(f"    {pkg}: NOT INSTALLED")

    # Seed registry
    lines.append(f"\n  Seed Registry:")
    if SEED_REGISTRY_PATH.exists():
        with open(SEED_REGISTRY_PATH) as f:
            registry = json.load(f)
        n_seeds = sum(
            isinstance(v, dict) and "seed" in v
            for section in registry.get("seeds", {}).values()
            for v in (section.values() if isinstance(section, dict) else [section])
        )
        lines.append(f"    Registry: {SEED_REGISTRY_PATH}")
        lines.append(f"    Sections: {len(registry.get('seeds', {}))}")
        lines.append(f"    Master seed: {registry['conventions']['master_seed']}")
    else:
        lines.append(f"    ⚠ NOT FOUND: {SEED_REGISTRY_PATH}")

    # File hashes
    lines.append(f"\n  Source File SHA-256 Hashes:")
    source_files = [
        "validation_framework.py",
        "validation/__init__.py",
        "validation/nested_cv.py",
        "validation/permutation_test.py",
        "validation/calibration.py",
        "validation/decision_curve.py",
        "validation/delong_test.py",
        "validation/stratified.py",
        "validation/confounders.py",
        "validation/bioinfo_benchmark.py",
        "validation/power_analysis.py",
        "reproducibility/seed_registry.json",
    ]
    base = Path(__file__).parent
    for sf in source_files:
        path = base / sf
        if path.exists():
            with open(path, "rb") as f:
                sha = hashlib.sha256(f.read()).hexdigest()
            lines.append(f"    {sf:45s} {sha}")
        else:
            lines.append(f"    {sf:45s} MISSING")

    return "\n".join(lines)


# ═════════════════════════════════════════════════════════════════════════════
# Main runner
# ═════════════════════════════════════════════════════════════════════════════

MODULES = [
    ("Nested Cross-Validation", run_module_1),
    ("Permutation Testing", run_module_2),
    ("Calibration Analysis", run_module_3),
    ("Decision Curve Analysis", run_module_4),
    ("DeLong Statistical Tests", run_module_5),
    ("Stratified Performance Analysis", run_module_6),
    ("Confounder Robustness Suite", run_module_7),
    ("Bioinformatic Tool Benchmark", run_module_8),
    ("Sample Size & Power Analysis", run_module_9),
    ("Reproducibility Verification", run_module_10),
]


def main():
    parser = argparse.ArgumentParser(
        description="DeepCatch Bioinformatics Validation Suite"
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Fast mode (reduced iterations for CI/testing)",
    )
    parser.add_argument(
        "--data", type=str, default=None,
        help="Path to data directory (uses synthetic data if not provided)",
    )
    parser.add_argument(
        "--output", type=str, default=str(REPORT_PATH),
        help="Output report path",
    )
    parser.add_argument(
        "--skip", type=int, nargs="*", default=[],
        help="Module numbers to skip (1-10)",
    )
    parser.add_argument(
        "--only", type=int, nargs="*", default=None,
        help="Only run specified module numbers (1-10)",
    )
    args = parser.parse_args()

    quick = args.quick or os.environ.get("DEEPCATCH_QUICK_MODE") == "1"

    # ── Load or generate data ──────────────────────────────────────────
    if args.data and Path(args.data).exists():
        print(f"Loading data from {args.data}...")
        # Real data loading would go here
        raise NotImplementedError("Real data loading not yet implemented")
    else:
        print("Generating synthetic data (n=1000)...")
        X, y, metadata = generate_synthetic_data(n_samples=1000)

    print(f"Data: {X.shape[0]} samples, {X.shape[1]} features, "
          f"prevalence={metadata['prevalence']:.2f}")

    # ── Run modules ────────────────────────────────────────────────────
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report_lines: List[str] = []
    timings: Dict[str, float] = {}
    failures: List[str] = []

    total_start = time.time()

    for i, (name, fn) in enumerate(MODULES, 1):
        # Skip logic
        if args.only is not None and i not in args.only:
            continue
        if i in args.skip:
            report_lines.append(f"\n## [{i}/10] {name} — SKIPPED\n")
            continue

        print(f"\n{'='*70}")
        print(f"[{i}/{len(MODULES)}] {name}")
        print(f"{'='*70}")

        t0 = time.time()

        try:
            if i == 6 or i == 7 or i == 8:
                output = fn(X, y, metadata, quick=quick)
            elif i == 9 or i == 10:
                output = fn(quick=quick) if i == 9 else fn()
            else:
                output = fn(X, y, quick=quick)
            elapsed = time.time() - t0
            timings[name] = elapsed

            report_lines.append(f"\n## [{i}/10] {name}\n")
            report_lines.append(f"_Runtime: {elapsed:.1f}s_\n")
            report_lines.append("```")
            report_lines.append(output)
            report_lines.append("```\n")

            print(f"✓ Completed in {elapsed:.1f}s")

        except Exception as e:
            elapsed = time.time() - t0
            failures.append(name)
            report_lines.append(f"\n## [{i}/10] {name} — FAILED\n")
            report_lines.append(f"_Runtime: {elapsed:.1f}s_\n")
            report_lines.append(f"**Error:** {type(e).__name__}: {e}\n")
            print(f"✗ FAILED: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()

    total_elapsed = time.time() - total_start

    # ── Build report ────────────────────────────────────────────────────
    header = f"""# DeepCatch Bioinformatics Validation Report

**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}
**Mode:** {'Quick (reduced iterations)' if quick else 'Full'}
**Data:** {'Synthetic (n=1000)' if not args.data else f'Real data from {args.data}'}
**Total Runtime:** {total_elapsed:.1f}s

---

## Executive Summary

| # | Module | Status | Runtime |
|---|--------|--------|---------|
"""

    for i, (name, _) in enumerate(MODULES, 1):
        status = "✓ PASSED" if name not in failures else "✗ FAILED"
        t = timings.get(name, 0)
        header += f"| {i} | {name} | {status} | {t:.1f}s |\n"

    header += f"\n**Passed:** {len(MODULES) - len(failures)}/{len(MODULES)}\n"
    if failures:
        header += f"**Failed:** {', '.join(failures)}\n"

    header += "\n---\n"

    full_report = header + "\n".join(report_lines)

    # ── Write report ────────────────────────────────────────────────────
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(full_report)

    print(f"\n{'='*70}")
    print(f"Report written to: {output_path}")
    print(f"Total runtime: {total_elapsed:.1f}s")
    print(f"Passed: {len(MODULES) - len(failures)}/{len(MODULES)}")
    if failures:
        print(f"Failed: {', '.join(failures)}")
        return 1
    print(f"{'='*70}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
