#!/usr/bin/env python3
"""
run_jiang_analysis.py — P0-C1
==============================

Standalone CLI orchestrator for the full Jiang 4-mer motif CET analysis
pipeline.  From flat frequency files to publication-ready results in one
command.

Pipeline (10 steps)
-------------------
1.  Load data via ``FrequencyDataset``
2.  CET analysis — Mann-Whitney U per motif
3.  Compute effect sizes — Cliff's delta
4.  FDR correction — Benjamini-Hochberg
5.  Rank motifs by composite score
6.  Logistic regression fusion on top-k motifs
7.  Cross-validation with AUC
8.  (Optional) Nested cross-validation
9.  Generate plots
10. Output summary report (Markdown)

Usage
-----
    python run_jiang_analysis.py -i jiang_4mer.xlsx -o results/jiang/
    python run_jiang_analysis.py -i data.csv -l labels.csv --top-k 30 --plot
    python run_jiang_analysis.py -i frequencies.npz --nested-cv --optimal-k

Dependencies
------------
    numpy, scipy, pandas, scikit-learn, matplotlib
    src/clinical/frequency_input.py      (FrequencyDataset, PlotGenerator)
    src/clinical/cet_cross_validator.py  (assumed available; graceful fallback)
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# ── Add repo root to path ──────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.clinical.frequency_input import FrequencyDataset, PlotGenerator
from src.clinical.clinical_interpretation import ClinicalReportGenerator

logger = logging.getLogger("jiang_analysis")

# ── import cet_cross_validator (graceful fallback) ─────────────────────────
try:
    from src.clinical.cet_cross_validator import NestedCETValidator, compute_cliffs_delta, MotifRanker
    _HAS_CET_CV = True
except ImportError:
    _HAS_CET_CV = False
    logger.debug("cet_cross_validator not available — nested CV disabled")


# ═══════════════════════════════════════════════════════════════════════════
# Core analysis functions
# ═══════════════════════════════════════════════════════════════════════════

def compute_cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta — non-parametric effect size.

    δ ∈ [−1, 1]; |δ| ≥ 0.147 = small, ≥ 0.33 = medium, ≥ 0.474 = large.

    References
    ----------
    Cliff N (1993) Psychol Bull 114:494-509
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return 0.0
    # Count pairwise comparisons
    gt = sum(int(xi > yj) for xi in x for yj in y)
    lt = sum(int(xi < yj) for xi in x for yj in y)
    return (gt - lt) / (nx * ny)


def _fast_cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Vectorised Cliff's delta for moderate-sized arrays."""
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return 0.0
    # Use broadcasting for speed
    if nx * ny <= 10_000_000:
        gt = np.sum(x[:, None] > y[None, :])
        lt = np.sum(x[:, None] < y[None, :])
    else:
        # Fallback to iterative for large arrays
        gt = sum(int(xi > yj) for xi in x for yj in y)
        lt = sum(int(xi < yj) for xi in x for yj in y)
    return (gt - lt) / (nx * ny)


def benjamini_hochberg_fdr(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Returns boolean mask of rejected hypotheses.
    """
    n = len(p_values)
    if n == 0:
        return np.zeros(0, dtype=bool)

    order = np.argsort(p_values)
    sorted_p = p_values[order]
    thresholds = (np.arange(1, n + 1) / n) * alpha

    below = sorted_p <= thresholds
    if np.any(below):
        k_max = np.max(np.where(below)[0])
        rejected_sorted = np.zeros(n, dtype=bool)
        rejected_sorted[:k_max + 1] = True
    else:
        rejected_sorted = np.zeros(n, dtype=bool)

    rejected = np.zeros(n, dtype=bool)
    rejected[order] = rejected_sorted
    return rejected


def composite_score(p_values: np.ndarray, effect_sizes: np.ndarray,
                    alpha: float = 0.05) -> np.ndarray:
    """Rank motifs by composite score: −log₁₀(p) × |δ|.

    Higher score = more significant *and* larger effect.
    Non-significant motifs (p ≥ alpha) get score = 0.
    """
    p_safe = np.maximum(p_values, 1e-300)
    neg_log_p = -np.log10(p_safe)
    abs_es = np.abs(effect_sizes)

    score = neg_log_p * abs_es
    score[p_values >= alpha] = 0.0
    return score


def find_optimal_k(scores: np.ndarray, max_k: int = 100) -> int:
    """Elbow method: find k where incremental gain flattens.

    Uses second derivative of the sorted score curve.
    """
    sorted_scores = np.sort(scores)[::-1]
    n = min(max_k, len(sorted_scores))
    ks = np.arange(1, n + 1)

    # Cumulative gain
    cumsum = np.cumsum(sorted_scores[:n])
    # Elbow: point of maximum curvature
    if len(ks) >= 3:
        # Second difference
        d2 = np.diff(cumsum, 2)
        elbow = int(ks[np.argmax(np.abs(d2)) + 1])
    else:
        elbow = n

    return max(5, min(elbow, max_k))


def run_cet_per_motif(X: np.ndarray, y: np.ndarray,
                      feature_names: Optional[List[str]] = None,
                      alpha: float = 0.05) -> pd.DataFrame:
    """Run Mann-Whitney U test for each motif (CET step 2).

    Returns DataFrame with columns: motif, p_value, effect_size, significant.
    """
    n_features = X.shape[1]
    if feature_names is None:
        feature_names = [f"motif_{i:04d}" for i in range(n_features)]

    pos_mask = y == 1
    neg_mask = y == 0
    pos_X = X[pos_mask]
    neg_X = X[neg_mask]

    if len(pos_X) < 3 or len(neg_X) < 3:
        raise ValueError(
            f"Need ≥ 3 samples per class for Mann-Whitney U. "
            f"Got {len(pos_X)} positive, {len(neg_X)} negative."
        )

    p_vals = np.empty(n_features)
    deltas = np.empty(n_features)

    for j in range(n_features):
        pos = pos_X[:, j]
        neg = neg_X[:, j]
        try:
            _, p = mannwhitneyu(pos, neg, alternative='two-sided')
        except ValueError:
            p = 1.0
        p_vals[j] = p
        deltas[j] = _fast_cliffs_delta(pos, neg)

    # FDR correction
    fdr_rejected = benjamini_hochberg_fdr(p_vals, alpha)
    scores = composite_score(p_vals, deltas, alpha)

    df = pd.DataFrame({
        'motif': feature_names,
        'p_value': p_vals,
        'effect_size': deltas,
        'abs_effect_size': np.abs(deltas),
        'fdr_significant': fdr_rejected,
        'composite_score': scores,
    })
    df.sort_values('composite_score', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    n_sig = int(fdr_rejected.sum())
    logger.info(
        "CET per-motif: %d/%d motifs significant at FDR %.3f",
        n_sig, n_features, alpha,
    )
    return df


def logistic_fusion_cv(X: np.ndarray, y: np.ndarray, top_k: int = 50,
                       n_folds: int = 5, seed: int = 42,
                       C: float = 10.0, select_by: str = 'p_value') -> Dict:
    """Logistic regression fusion on top-k motifs with CV AUC.

    Uses C=10.0 (weaker regularization) by default — for strong-signal
    data like Jiang 4-mer, less L2 penalty preserves more discriminative
    information and avoids over-shrinking correlated motif features.

    Parameters
    ----------
    select_by : str
        'p_value' — select top-k by Mann-Whitney U p-value (recommended)
        'variance' — select top-k by feature variance
        'composite' — use pre-computed composite score from cet_df

    Returns dict with keys: auc_mean, auc_std, auc_folds, coefs, intercept,
                            top_indices, n_top_motifs_used.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    from scipy.stats import mannwhitneyu

    # Select top-k features
    if select_by == 'variance':
        scores = np.var(X, axis=0)
        top_idx = np.argsort(-scores)[:min(top_k, X.shape[1])]
    elif select_by == 'p_value':
        p_vals = []
        for i in range(X.shape[1]):
            try:
                _, p = mannwhitneyu(X[y == 1, i], X[y == 0, i],
                                    alternative='two-sided')
                p_vals.append(p)
            except (ValueError, ZeroDivisionError):
                p_vals.append(1.0)
        top_idx = np.argsort(p_vals)[:min(top_k, X.shape[1])]
    else:
        # Fallback: variance
        scores = np.var(X, axis=0)
        top_idx = np.argsort(-scores)[:min(top_k, X.shape[1])]

    X_top = X[:, top_idx]

    if X_top.shape[1] == 0:
        return {
            'auc_mean': 0.5, 'auc_std': 0.0,
            'auc_folds': [0.5] * n_folds,
            'coefs': np.zeros(0), 'intercept': 0.0,
            'top_indices': [],
        }

    lr = LogisticRegression(
        C=C, solver='liblinear', max_iter=5000, random_state=seed,
    )

    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)

    # Per-fold AUC (explicit loop for detail, not cross_val_predict)
    fold_aucs = []
    for train_idx, test_idx in cv.split(X_top, y):
        lr_fold = LogisticRegression(
            C=C, solver='liblinear', max_iter=5000, random_state=seed,
        )
        lr_fold.fit(X_top[train_idx], y[train_idx])
        fold_pred = lr_fold.predict_proba(X_top[test_idx])[:, 1]
        try:
            fold_auc = roc_auc_score(y[test_idx], fold_pred)
        except ValueError:
            fold_auc = 0.5
        fold_aucs.append(fold_auc)

    auc_mean = float(np.mean(fold_aucs))
    auc_std = float(np.std(fold_aucs, ddof=1))

    # Fit on all data for coefficients
    lr.fit(X_top, y)

    # Full cross_val_predict for ROC curve
    y_pred_cv = cross_val_predict(lr, X_top, y, cv=cv, method='predict_proba')[:, 1]

    return {
        'auc_mean': auc_mean,
        'auc_std': auc_std,
        'auc_folds': fold_aucs,
        'coefs': lr.coef_.flatten(),
        'intercept': float(lr.intercept_[0]),
        'top_indices': top_idx.tolist(),
        'y_pred_cv': y_pred_cv,
        'n_top_motifs_used': X_top.shape[1],
    }


# ═══════════════════════════════════════════════════════════════════════════
# Report generation
# ═══════════════════════════════════════════════════════════════════════════

def generate_summary_report(
    cet_df: pd.DataFrame,
    fusion_result: Dict,
    nested_cv_result: Optional[Dict] = None,
    args: argparse.Namespace = None,
    elapsed: float = 0.0,
) -> str:
    """Generate a Markdown summary report."""
    lines = [
        "# Jiang 4-mer Motif Analysis — Summary Report",
        "",
        f"**Date**: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}",
        f"**Runtime**: {elapsed:.1f}s",
        "",
        "## Parameters",
        "",
        f"| Parameter | Value |",
        f"|-----------|-------|",
    ]
    if args:
        lines.append(f"| Input file | `{args.input}` |")
        lines.append(f"| Output dir | `{args.output}` |")
        lines.append(f"| Top-k motifs | {args.top_k} |")
        lines.append(f"| Significance α | {args.alpha} |")
        lines.append(f"| Nested CV | {'✓' if args.nested_cv else '✗'} |")
        lines.append(f"| Optimal-k | {'✓' if args.optimal_k else '✗'} |")

    lines.extend([
        "",
        "## 1. Data Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total motifs tested | {len(cet_df)} |",
        f"| Significant (FDR) | {cet_df['fdr_significant'].sum()} |",
        f"| Median p-value | {cet_df['p_value'].median():.2e} |",
        f"| Min p-value | {cet_df['p_value'].min():.2e} |",
        f"| Median |effect size| | {cet_df['abs_effect_size'].median():.4f} |",
        f"| Max |effect size| | {cet_df['abs_effect_size'].max():.4f} |",
    ])

    # Top motifs table
    lines.extend([
        "",
        "## 2. Top Significant Motifs",
        "",
        "| Rank | Motif | p-value | Effect Size (δ) | FDR Sig | Score |",
        "|------|-------|---------|-----------------|---------|-------|",
    ])
    top_n = min(20, len(cet_df))
    for i, row in cet_df.head(top_n).iterrows():
        sig_mark = '✓' if row['fdr_significant'] else ''
        lines.append(
            f"| {i + 1} | {row['motif']} | {row['p_value']:.2e} | "
            f"{row['effect_size']:+.4f} | {sig_mark} | {row['composite_score']:.2f} |"
        )

    # Fusion results
    lines.extend([
        "",
        "## 3. Logistic Regression Fusion",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Top-k motifs used | {args.top_k if args else 50} |",
        f"| CV AUC (mean ± std) | {fusion_result['auc_mean']:.4f} ± {fusion_result['auc_std']:.4f} |",
        f"| Per-fold AUCs | {[f'{a:.4f}' for a in fusion_result['auc_folds']]} |",
    ])

    if nested_cv_result is not None:
        lines.extend([
            "",
            "## 4. Nested Cross-Validation",
            "",
            f"| Metric | Value |",
            f"|--------|-------|",
            f"| Nested CV AUC (mean ± std) | {nested_cv_result.get('mean_outer_score', 'N/A')} ± {nested_cv_result.get('std_outer_score', 'N/A')} |",
            f"| Number of outer folds | {nested_cv_result.get('n_outer_folds', 'N/A')} |",
        ])

    lines.extend([
        "",
        "## 5. Interpretation",
        "",
        f"- **{cet_df['fdr_significant'].sum()} motifs** are significantly "
        f"differentially abundant between cancer and control after FDR correction.",
        f"- The logistic regression fusion achieves **AUC = {fusion_result['auc_mean']:.4f}** "
        f"in {len(fusion_result['auc_folds'])}-fold CV.",
        "",
        "---",
        "*Generated by DeepCatch Jiang Analysis Pipeline*",
    ])

    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Jiang 4-mer Motif CET Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_jiang_analysis.py -i jiang_4mer.xlsx
  python run_jiang_analysis.py -i data.csv -l labels.csv --top-k 30 --plot
  python run_jiang_analysis.py -i frequencies.npz --nested-cv --optimal-k
        """,
    )
    p.add_argument('--input', '-i', required=True,
                   help='Path to frequency data file (.xlsx, .csv, .npy, .npz)')
    p.add_argument('--labels', '-l', default=None,
                   help='Path to label file (optional; auto-detected if in same file)')
    p.add_argument('--output', '-o', default='results/jiang_reanalysis/',
                   help='Output directory (default: results/jiang_reanalysis/)')
    p.add_argument('--cancer-type', default=None,
                   help='Cancer type column name (for multi-cancer datasets)')
    p.add_argument('--control-label', default=None,
                   help='Label value for control group (auto-detected if binary)')
    p.add_argument('--top-k', type=int, default=50,
                   help='Number of top motifs for fusion (default: 50)')
    p.add_argument('--alpha', type=float, default=0.05,
                   help='Significance threshold for FDR (default: 0.05)')
    p.add_argument('--nested-cv', action='store_true',
                   help='Enable nested cross-validation')
    p.add_argument('--optimal-k', action='store_true',
                   help='Find optimal k via elbow method')
    p.add_argument('--plot', '-p', action='store_true',
                   help='Generate plots (volcano, heatmap, ROC, feature importance)')
    p.add_argument('--seed', type=int, default=42,
                   help='Random seed (default: 42)')
    p.add_argument('--report', '-r', action='store_true',
                   help='Generate clinical interpretation report (HTML + JSON)')
    p.add_argument('--lr-C', type=float, default=10.0,
                   help='LogisticRegression C (inverse reg strength, default 10.0)')
    p.add_argument('--select-by', default='p_value',
                   choices=['p_value', 'variance'],
                   help='Feature selection method (default: p_value)')
    p.add_argument('--verbose', '-v', action='store_true',
                   help='Verbose output')
    return p


def _run_all_pairwise(ds: 'FrequencyDataset', out_dir: Path,
                      args: argparse.Namespace) -> None:
    """Run binary CET analysis for each cancer type vs all controls."""
    from src.clinical.frequency_input import FrequencyDataset

    cancer_dir = out_dir / 'per_cancer'
    cancer_dir.mkdir(parents=True, exist_ok=True)

    # Get unique classes from raw labels
    if hasattr(ds, '_raw_labels_str'):
        unique_labels_arr = np.unique(ds._raw_labels_str)
    else:
        unique_labels_arr = np.unique(ds.y)
    logger.info("Running pairwise comparisons for labels: %s", unique_labels_arr)

    # Determine control label (typically 'Control')
    control_candidates = ['Control', 'control', 'Healthy', 'healthy', 'Normal', 'normal', '0', 0]
    control_label = None
    for c in control_candidates:
        if c in unique_labels_arr:
            control_label = c
            break
    if control_label is None:
        control_label = unique_labels_arr[0]
        logger.warning("No 'Control' label found, using '%s' as reference", control_label)

    all_results = {}
    for label in unique_labels_arr:
        if label == control_label:
            continue
        logger.info("─" * 40)
        logger.info("Comparison: %s vs %s", label, control_label)

        try:
            ds_copy = FrequencyDataset.__new__(FrequencyDataset)
            ds_copy.__dict__.update(ds.__dict__.copy())
            ds_copy.filter_by_label(str(label), str(control_label))
            X, y, fnames = ds_copy.X, ds_copy.y, ds_copy.feature_names

            n_pos = int(y.sum())
            n_neg = int((1 - y).sum())
            if n_pos < 3 or n_neg < 3:
                logger.warning("  Skipping: insufficient samples (%d pos, %d neg)", n_pos, n_neg)
                continue

            cet_df = run_cet_per_motif(X, y, fnames, args.alpha)
            optimal_k = args.top_k
            if args.optimal_k:
                optimal_k = find_optimal_k(cet_df['composite_score'].values, max_k=min(200, len(cet_df)))

            fusion_result = logistic_fusion_cv(X, y, top_k=min(optimal_k, X.shape[1]),
                                seed=args.seed, C=args.lr_C,
                                select_by=args.select_by)

            label_dir = cancer_dir / str(label).replace(' ', '_')
            label_dir.mkdir(parents=True, exist_ok=True)
            cet_df.to_csv(label_dir / 'cet_motif_results.csv', index=False)

            report_md = generate_summary_report(cet_df, fusion_result, None, args, 0.0)
            (label_dir / 'summary_report.md').write_text(report_md)

            all_results[str(label)] = {
                'auc_mean': fusion_result['auc_mean'],
                'auc_std': fusion_result['auc_std'],
                'n_significant': int(cet_df['fdr_significant'].sum()),
                'n_samples': len(y),
                'optimal_k': optimal_k,
            }
            logger.info("  %s vs %s: AUC=%.4f, %d sig. motifs",
                        label, control_label, fusion_result['auc_mean'],
                        int(cet_df['fdr_significant'].sum()))
        except Exception as e:
            logger.error("  Failed for %s: %s", label, e)
            import traceback
            logger.debug(traceback.format_exc())

    # Summary table
    if all_results:
        summary_lines = [
            "# Per-Cancer CET Analysis — Summary",
            "",
            "| Cancer Type | AUC | Std | Significant Motifs | Samples (case+ctrl) | Optimal k |",
            "|------------|-----|-----|-------------------|--------------------|----------|",
        ]
        for label, r in sorted(all_results.items()):
            summary_lines.append(
                f"| {label} | {r['auc_mean']:.4f} | {r['auc_std']:.4f} | "
                f"{r['n_significant']} | {r['n_samples']} | {r['optimal_k']} |"
            )
        (out_dir / 'per_cancer_summary.md').write_text('\n'.join(summary_lines))
        logger.info("\nPer-cancer summary → %s", out_dir / 'per_cancer_summary.md')


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%H:%M:%S',
    )

    t_start = time.time()

    # ── Step 0: output directory ───────────────────────────────────────
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir / 'plots'
    plot_dir.mkdir(parents=True, exist_ok=True)

    logger.info("═" * 60)
    logger.info("Jiang 4-mer Motif CET Analysis Pipeline")
    logger.info("═" * 60)
    logger.info("Input:      %s", args.input)
    logger.info("Output:     %s", out_dir)
    logger.info("Top-k:      %d", args.top_k)
    logger.info("Alpha:      %.3f", args.alpha)
    logger.info("Nested CV:  %s", args.nested_cv)
    logger.info("Plot:       %s", args.plot)

    # ── Step 1: Load data ──────────────────────────────────────────────
    logger.info("─" * 40)
    logger.info("Step 1/10: Loading data via FrequencyDataset …")

    try:
        ds = FrequencyDataset(args.input, label_path=args.labels)
        X, y, feature_names = ds.load()
    except Exception as e:
        logger.error("Failed to load data: %s", e)
        return 1

    # ── Validate ───────────────────────────────────────────────────────
    report = ds.validate()
    if not report['passed']:
        logger.warning("Data validation found issues — proceeding with caution")
    print("\n" + ds.describe())

    if report.get('n_classes', 2) > 2:
        logger.warning(
            "Detected %d classes in labels. Use --cancer-type and --control-label "
            "to select binary comparison, or run per-cancer.",
            report.get('n_classes', '?')
        )
        if args.cancer_type:
            # Filter to specific cancer type vs control
            logger.info("Filtering: %s vs %s", args.cancer_type, args.control_label or 'rest')
            ds.filter_by_label(args.cancer_type, args.control_label)
            X, y, feature_names = ds.X, ds.y, ds.feature_names
            n_pos = int(np.sum(y == 1))
            n_neg = int(np.sum(y == 0))
            logger.info("  After filtering: %d cases, %d controls (total %d samples)",
                        n_pos, n_neg, len(y))
        else:
            # Auto-run all comparisons
            logger.info("No --cancer-type specified. Running all pairwise comparisons…")
            _run_all_pairwise(ds, out_dir, args)
            logger.info("═" * 60)
            logger.info("✅ All pairwise comparisons complete in %.1fs", time.time() - t_start)
            logger.info("═" * 60)
            return 0

    # ── Step 2: CET per-motif ──────────────────────────────────────────
    logger.info("Step 2/10: CET analysis (Mann-Whitney U per motif) …")
    cet_df = run_cet_per_motif(X, y, feature_names, args.alpha)

    # ── Steps 3-4 already done inside run_cet_per_motif ────────────────
    # ── Step 5: Rank by composite score ────────────────────────────────
    logger.info("Step 5/10: Ranking motifs by composite score …")
    cet_df.to_csv(out_dir / 'cet_motif_results.csv', index=False)
    logger.info("  Saved → %s", out_dir / 'cet_motif_results.csv')

    # ── Optimal k ──────────────────────────────────────────────────────
    optimal_k = args.top_k
    if args.optimal_k:
        logger.info("Step 5b: Finding optimal k via elbow method …")
        scores = cet_df['composite_score'].values
        optimal_k = find_optimal_k(scores, max_k=min(200, len(scores)))
        logger.info("  Optimal k = %d", optimal_k)
    else:
        optimal_k = args.top_k

    # ── Step 6: Logistic regression fusion ─────────────────────────────
    logger.info(
        "Step 6/10: Logistic regression fusion on top-%d motifs …",
        optimal_k,
    )
    fusion_result = logistic_fusion_cv(X, y, top_k=optimal_k, seed=args.seed,
                                    C=args.lr_C, select_by=args.select_by)
    logger.info("  CV AUC = %.4f ± %.4f (%d motifs, C=%.1f, select=%s)",
                fusion_result['auc_mean'], fusion_result['auc_std'],
                fusion_result['n_top_motifs_used'], args.lr_C, args.select_by)

    # Save fusion coefs
    coef_df = pd.DataFrame({
        'feature': [feature_names[i] for i in fusion_result['top_indices']],
        'coefficient': fusion_result['coefs'],
    })
    coef_df.sort_values('coefficient', key=abs, ascending=False, inplace=True)
    coef_df.to_csv(out_dir / 'fusion_coefficients.csv', index=False)

    # ── Step 7: CV AUC already done ────────────────────────────────────
    logger.info("Step 7/10: Cross-validation AUC = %.4f ± %.4f (already computed)",
                fusion_result['auc_mean'], fusion_result['auc_std'])

    # ── Step 8: Nested CV (optional) ───────────────────────────────────
    nested_cv_result = None
    if args.nested_cv:
        logger.info("Step 8/10: Nested cross-validation …")
        if _HAS_CET_CV:
            try:
                cv = NestedCETValidator(
                    outer_folds=5, inner_folds=3, k_motifs=optimal_k,
                    random_state=args.seed,
                )
                ncv = cv.validate(X, y)
                outer_aucs = ncv.get('outer_aucs', [])
                mean_auc = float(np.mean(outer_aucs)) if outer_aucs else float('nan')
                std_auc = float(np.std(outer_aucs, ddof=1)) if len(outer_aucs) > 1 else 0.0
                nested_cv_result = {
                    'mean_outer_score': mean_auc,
                    'std_outer_score': std_auc,
                    'n_outer_folds': len(outer_aucs),
                }
                logger.info(
                    "  Nested CV AUC = %.4f (std=%.4f, %d folds)",
                    mean_auc, std_auc, len(outer_aucs),
                )
            except Exception as e:
                logger.warning("NestedCETValidator failed: %s — falling back to sklearn", e)
                nested_cv_result = _run_sklearn_nested_cv(X, y, optimal_k, args.seed)
        else:
            logger.info("  Nested CV: using sklearn fallback")
            nested_cv_result = _run_sklearn_nested_cv(X, y, optimal_k, args.seed)

        if nested_cv_result:
            pd.DataFrame([nested_cv_result]).to_csv(
                out_dir / 'nested_cv_results.csv', index=False,
            )
    else:
        logger.info("Step 8/10: Nested CV skipped (use --nested-cv to enable)")

    # ── Step 9: Plots ──────────────────────────────────────────────────
    if args.plot:
        logger.info("Step 9/10: Generating plots …")
        _make_plots(cet_df, X, y, feature_names, fusion_result, plot_dir, args)
    else:
        logger.info("Step 9/10: Plots skipped (use --plot to enable)")

    # ── Step 10: Summary report ────────────────────────────────────────
    logger.info("Step 10/10: Generating summary report …")
    elapsed = time.time() - t_start
    report_md = generate_summary_report(
        cet_df, fusion_result, nested_cv_result, args, elapsed,
    )
    report_path = out_dir / 'summary_report.md'
    report_path.write_text(report_md)
    logger.info("  Report saved → %s", report_path)

    # ── Step 11: Clinical report (optional) ────────────────────────────
    if args.report:
        logger.info("Step 11/11: Generating clinical interpretation report …")
        try:
            crg = ClinicalReportGenerator(
                cet_df, fusion_result,
                threshold_sens=0.70, threshold_spec=0.95,
            )

            # HTML report
            html = crg.generate_html_report()
            html_path = out_dir / 'clinical_report.html'
            html_path.write_text(html)
            logger.info("  Clinical HTML → %s", html_path)

            # JSON export
            json_path = out_dir / 'clinical_report.json'
            crg.export_json(str(json_path))
            logger.info("  Clinical JSON → %s", json_path)

            # Text briefing
            briefing = crg.generate_briefing()
            briefing_path = out_dir / 'clinical_briefing.txt'
            briefing_path.write_text(briefing)
            logger.info("  Clinical briefing → %s", briefing_path)
        except Exception as e:
            logger.error("Clinical report generation failed: %s", e)
    else:
        logger.info("Step 11/11: Clinical report skipped (use --report to enable)")

    # ── Done ───────────────────────────────────────────────────────────
    logger.info("═" * 60)
    logger.info("✅ Pipeline complete in %.1fs", elapsed)
    logger.info("   Results directory: %s", out_dir)
    logger.info("   CET results:       %s", out_dir / 'cet_motif_results.csv')
    logger.info("   Fusion coefs:      %s", out_dir / 'fusion_coefficients.csv')
    logger.info("   Summary report:    %s", report_path)
    if args.report:
        logger.info("   Clinical HTML:     %s", out_dir / 'clinical_report.html')
        logger.info("   Clinical JSON:     %s", out_dir / 'clinical_report.json')
        logger.info("   Clinical briefing: %s", out_dir / 'clinical_briefing.txt')
    logger.info("═" * 60)

    return 0


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _run_sklearn_nested_cv(X: np.ndarray, y: np.ndarray,
                           top_k: int = 50, seed: int = 42) -> Dict:
    """Fallback nested CV using sklearn directly."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_score

    var = np.var(X, axis=0)
    top_idx = np.argsort(-var)[:min(top_k, X.shape[1])]
    X_top = X[:, top_idx]

    if X_top.shape[1] == 0:
        return {'mean_outer_score': 0.5, 'ci95_outer': [0.5, 0.5],
                'optimism_gap': 0.0}

    lr = LogisticRegression(solver='liblinear', max_iter=1000, random_state=seed)
    param_grid = {'C': [0.01, 0.1, 1.0, 10.0]}

    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)

    gs = GridSearchCV(lr, param_grid, cv=inner_cv, scoring='roc_auc')
    outer_scores = cross_val_score(gs, X_top, y, cv=outer_cv, scoring='roc_auc')

    mean_outer = float(np.mean(outer_scores))
    std_outer = float(np.std(outer_scores, ddof=1))

    from scipy import stats
    if len(outer_scores) > 1:
        se = std_outer / np.sqrt(len(outer_scores))
        t_crit = stats.t.ppf(0.975, df=len(outer_scores) - 1)
        ci95 = (mean_outer - t_crit * se, mean_outer + t_crit * se)
    else:
        ci95 = (mean_outer, mean_outer)

    # Inner score approximation
    inner_score = float(np.mean([
        GridSearchCV(lr, param_grid, cv=inner_cv, scoring='roc_auc')
        .fit(X_top, y).best_score_
    ]))

    return {
        'mean_outer_score': mean_outer,
        'std_outer_score': std_outer,
        'ci95_outer': [float(ci95[0]), float(ci95[1])],
        'optimism_gap': float(inner_score - mean_outer),
        'fold_scores': [float(s) for s in outer_scores],
    }


def _make_plots(cet_df: pd.DataFrame, X: np.ndarray, y: np.ndarray,
                feature_names: List[str], fusion_result: Dict,
                plot_dir: Path, args: argparse.Namespace) -> None:
    """Generate all four standard plots."""
    pg = PlotGenerator()

    # Volcano
    try:
        pg.volcano_plot(
            cet_df['p_value'].values,
            cet_df['effect_size'].values,
            labels=cet_df['motif'].tolist(),
            title='Jiang 4-mer Motif Volcano Plot',
            save_path=str(plot_dir / 'volcano.png'),
            alpha=args.alpha,
            annotate_top=15,
        )
    except Exception as e:
        logger.warning("Volcano plot failed: %s", e)

    # Heatmap
    try:
        n_top_hm = min(30, X.shape[1])
        sample_lbl = ['Cancer' if yi else 'Control' for yi in y]
        pg.heatmap(
            X, feature_names=feature_names,
            sample_labels=sample_lbl, n_top=n_top_hm,
            save_path=str(plot_dir / 'heatmap.png'),
        )
    except Exception as e:
        logger.warning("Heatmap failed: %s", e)

    # ROC
    try:
        if 'y_pred_cv' in fusion_result and len(np.unique(y)) >= 2:
            pg.roc_curve(
                y, fusion_result['y_pred_cv'],
                title='Logistic Regression Fusion ROC',
                save_path=str(plot_dir / 'roc_curve.png'),
            )
    except Exception as e:
        logger.warning("ROC plot failed: %s", e)

    # Feature importance
    try:
        if len(fusion_result['coefs']) > 0:
            top_names = [feature_names[i] for i in fusion_result['top_indices']]
            pg.feature_importance(
                fusion_result['coefs'],
                feature_names=top_names,
                n_top=min(20, len(fusion_result['coefs'])),
                title='Top Motif Logistic Regression Coefficients',
                save_path=str(plot_dir / 'feature_importance.png'),
            )
    except Exception as e:
        logger.warning("Feature importance plot failed: %s", e)

    logger.info("  Plots saved → %s/", plot_dir)


# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    sys.exit(main())
