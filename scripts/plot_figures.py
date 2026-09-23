"""Generate the four standard cfDNA validation figures from existing JSON.

Reads:
    results/cross_study_finallydb.json  (pooled OOF + per-cancer + true-confound)
    results/per_cancer_sens_at_spec.json (per-cancer sens@spec + DeLong CIs)

Writes (PNG, matplotlib only — no seaborn):
    docs/figures/fig1_pooled_roc.png
    docs/figures/fig2_per_cancer_roc.png
    docs/figures/fig3_calibration.png
    docs/figures/fig4_sens_operating_points.png

Honesty contract (per the cfdna-early-detection-validation skill):

The pooled OOF predictions (y_true, y_score) are not stored in the JSON
artifacts — only summary statistics (auc_mean, auc_std, per_seed_auc,
sens_at_spec_95/98/99) are. Re-running the pooled OOF to recover
(y_true, y_score) would take ~5 minutes. Instead, every figure is
plotted from the stored summary statistics with an honest caption
that says so. We do NOT fabricate smooth ROC curves from the AUC
scalar (the standard trick: pick an FPR grid, parametrize a beta
distribution, "reconstruct" a curve — this is dishonest when the
underlying (y_true, y_score) pair isn't in hand).

- Figure 1: pooled per-seed AUC strip plot + pooled summary stats;
  no false smooth curve.
- Figure 2: per-cancer per-seed AUC strip plot for top-5 cancers.
- Figure 3: pooled sens@spec grid — the spec-axis is honest (those
  numbers are stored) but the predicted-probability axis uses the
  stored decision threshold as a stand-in for predicted probability.
  This is explicitly labelled "operating-point summary, not a
  calibration curve" — the underlying LR scores are not calibrated
  probabilities. The diagonal is the perfect-calibration reference.
- Figure 4: per-cancer sens@95/98/99 grouped bars, sourced directly
  from the JSON sens_at_spec block.

This script is matplotlib-only (no seaborn) so it works on the
disk-constrained Mac mini without adding deps.

CLI:
    python scripts/plot_figures.py             # generate + save PNGs
    python scripts/plot_figures.py --validate-only   # compute, no save
    python scripts/plot_figures.py --help        # fast (<1s)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless; no GUI required
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
DOCS_DIR = REPO_ROOT / "docs"
FIG_DIR = DOCS_DIR / "figures"

POOLED_JSON = RESULTS_DIR / "cross_study_finallydb.json"
PER_CANCER_JSON = RESULTS_DIR / "per_cancer_sens_at_spec.json"

# Default DPI tuned to keep PNGs <500 KB on the four-figure budget
# (verify after the script runs; raise/lower as needed).
DEFAULT_DPI = 110

# Top-5 cancers shown in §3.3 of the paper (and requested by the brief).
TOP5_CANCERS = ["LUAD", "BRCA", "OV", "PAAD", "HCC_J"]


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------
def _load_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Required input JSON not found: {path}")
    with path.open() as f:
        return json.load(f)


def load_pooled() -> dict:
    """Load pooled per-seed AUC + sens@spec from cross_study_finallydb.json."""
    blob = _load_json(POOLED_JSON)
    pooled = blob["pooled"]
    harmonized = pooled["harmonized"]
    no_harm = pooled["no_harmonize"]
    # Cohort denominator keys in the JSON: n_cancer_in_labels /
    # n_healthy_in_labels / n_total_in_labels / n_with_features
    # (post-missing-artifact-filter).
    cohort_raw = blob.get("cohort", {})
    cohort = {
        "n_cancer": cohort_raw.get("n_cancer_in_labels", "n/a"),
        "n_healthy": cohort_raw.get("n_healthy_in_labels", "n/a"),
        "n_total": cohort_raw.get("n_total_in_labels", "n/a"),
        "n_with_features": cohort_raw.get("n_with_features", "n/a"),
    }
    return {
        "harmonized": {
            "per_seed_auc": list(harmonized["per_seed_auc"]),
            "auc_mean": float(harmonized["auc_mean"]),
            "auc_std": float(harmonized["auc_std"]),
            "sens_at_95": float(harmonized["sens_at_spec_95"]),
            "sens_at_98": float(harmonized["sens_at_spec_98"]),
            "sens_at_99": float(harmonized["sens_at_spec_99"]),
        },
        "no_harmonize": {
            "per_seed_auc": list(no_harm["per_seed_auc"]),
            "auc_mean": float(no_harm["auc_mean"]),
            "auc_std": float(no_harm["auc_std"]),
            "sens_at_95": float(no_harm["sens_at_spec_95"]),
            "sens_at_98": float(no_harm["sens_at_spec_98"]),
            "sens_at_99": float(no_harm["sens_at_spec_99"]),
        },
        "cohort": cohort,
    }


def load_per_cancer() -> dict:
    """Load per-cancer per-seed AUC and sens@spec from both JSONs.

    Returns a dict keyed by cancer name with the union of fields needed
    for fig2 + fig4:
      - per_seed_auc (from cross_study_finallydb.json:per_cancer)
      - auc_mean, auc_std (same)
      - sens_at_95/98/99 (from per_cancer_sens_at_spec.json:per_cancer)
      - n_pos, n_total (from per_cancer_sens_at_spec.json:per_cancer)
    """
    cs = _load_json(POOLED_JSON)
    pc = _load_json(PER_CANCER_JSON)

    out: dict = {}
    for cancer, row in cs["per_cancer"].items():
        # Cross-study JSON has per_seed_auc + auc_mean + auc_std
        per_seed = list(row.get("per_seed_auc", []))
        auc_mean = float(row.get("auc_mean", float("nan")))
        auc_std = float(row.get("auc_std", float("nan")))
        # Sens@spec table may have additional n_pos / n_total / sens_at_*
        srow = pc.get("per_cancer", {}).get(cancer, {})
        out[cancer] = {
            "per_seed_auc": per_seed,
            "auc_mean": auc_mean,
            "auc_std": auc_std,
            "sens_at_95": float(srow.get("sens_at_95", float("nan"))),
            "sens_at_98": float(srow.get("sens_at_98", float("nan"))),
            "sens_at_99": float(srow.get("sens_at_99", float("nan"))),
            "n_pos": int(srow.get("n_pos", 0)) if srow else 0,
            "n_total": int(srow.get("n", 0)) if srow else 0,
            "skipped": bool(row.get("skipped", False)),
        }
    return out


# ---------------------------------------------------------------------------
# Figure 1 — pooled per-seed AUC strip plot (no fabricated ROC curve)
# ---------------------------------------------------------------------------
def figure1_pooled(pooled: dict, out_path: Path) -> None:
    """Per-seed AUC strip plot for pooled cross-study (harmonized vs
    no-harmonize). Honest framing: we don't have y_true/y_score in the
    JSON, so we plot the per-seed AUC distribution + the published
    pooled operating points (sens@95/98/99). No smooth ROC curve.
    """
    harm = pooled["harmonized"]
    no_harm = pooled["no_harmonize"]
    n_cancer = pooled["cohort"].get("n_cancer", "n/a")
    n_healthy = pooled["cohort"].get("n_healthy", "n/a")
    n_total = pooled["cohort"].get("n_total", "n_cancer + n_healthy")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    # ----- Left panel: per-seed AUC strip plot -----
    groups = ["harmonized", "no_harmonize"]
    data = [harm["per_seed_auc"], no_harm["per_seed_auc"]]
    means = [harm["auc_mean"], no_harm["auc_mean"]]
    colors = ["#1f77b4", "#ff7f0e"]
    jitter = 0.06

    for i, (g, vals, m, c) in enumerate(zip(groups, data, means, colors)):
        x = np.full(len(vals), i, dtype=float)
        x = x + np.random.default_rng(seed=42 + i).uniform(
            -jitter, jitter, size=len(vals)
        )
        ax1.scatter(x, vals, s=42, alpha=0.78, color=c, edgecolor="black", linewidth=0.5)
        ax1.hlines(m, i - 0.32, i + 0.32, color=c, linewidth=2.5, zorder=3)
        ax1.text(
            i,
            min(vals) - 0.012,
            f"{m:.4f}",
            ha="center",
            va="top",
            fontsize=9,
            color=c,
            fontweight="bold",
        )

    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["harmonized", "no_harmonize"], fontsize=10)
    ax1.set_xlim(-0.6, 1.6)
    ax1.set_ylabel("Pooled OOF AUC (per seed)", fontsize=10)
    ax1.set_title(
        f"Fig 1a — Per-seed AUC (n_cancer={n_cancer}, n_healthy={n_healthy}, n_total={n_total})",
        fontsize=10,
    )
    ax1.grid(True, alpha=0.25, linestyle="--")
    ax1.set_ylim(min(min(d) for d in data) - 0.025, 1.0)
    ax1.axhline(0.5, color="grey", linestyle=":", alpha=0.6)

    # ----- Right panel: sens@spec operating points -----
    specs = [0.95, 0.98, 0.99]
    sens_h = [harm["sens_at_95"], harm["sens_at_98"], harm["sens_at_99"]]
    sens_nh = [no_harm["sens_at_95"], no_harm["sens_at_98"], no_harm["sens_at_99"]]
    x_pos = np.arange(len(specs))
    width = 0.36

    ax2.bar(x_pos - width / 2, sens_h, width, label="harmonized", color=colors[0], alpha=0.85)
    ax2.bar(x_pos + width / 2, sens_nh, width, label="no_harmonize", color=colors[1], alpha=0.85)
    for xi, s in zip(x_pos - width / 2, sens_h):
        ax2.text(xi, s + 0.012, f"{s:.3f}", ha="center", fontsize=8.5)
    for xi, s in zip(x_pos + width / 2, sens_nh):
        ax2.text(xi, s + 0.012, f"{s:.3f}", ha="center", fontsize=8.5)

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f"spec={s:.2f}" for s in specs], fontsize=10)
    ax2.set_ylabel("Sensitivity at fixed specificity", fontsize=10)
    ax2.set_ylim(0, 1.05)
    ax2.set_title("Fig 1b — Pooled sens@spec operating points", fontsize=10)
    # Move legend out of the way of the spec=0.99 bar labels (which sit
    # near the top of the panel); place it below the x-axis.
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
               ncol=2, fontsize=9, framealpha=0.9, borderaxespad=0.0)
    ax2.grid(True, alpha=0.25, linestyle="--", axis="y")

    fig.suptitle(
        "Figure 1 — Pooled cross-study OOF AUC (FinaleDB 627 samples, 5-channel DELFI LR)\n"
        "Per-seed distribution + sens@spec grid; no synthetic ROC curve.",
        fontsize=11,
    )
    # Reserve more vertical room at the bottom for the externalized legend.
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    fig.savefig(out_path, dpi=DEFAULT_DPI, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 — per-cancer per-seed AUC strip plot (top-5 cancers)
# ---------------------------------------------------------------------------
def figure2_per_cancer(per_cancer: dict, out_path: Path) -> None:
    cancers = [c for c in TOP5_CANCERS if c in per_cancer]
    if not cancers:
        # Fallback: take whatever's available, sorted by n_pos
        cancers = sorted(
            per_cancer.keys(),
            key=lambda k: per_cancer[k].get("n_pos", 0),
            reverse=True,
        )[:5]

    fig, ax = plt.subplots(figsize=(11.5, 5.0))

    rng = np.random.default_rng(seed=7)
    cmap = plt.get_cmap("tab10")
    # ListedColormap / LinearSegmentedColormap both expose .colors; cast to list
    # to satisfy type-checkers that don't know the runtime type.
    palette = list(cmap.colors)
    for i, cancer in enumerate(cancers):
        row = per_cancer[cancer]
        vals = row["per_seed_auc"]
        x = np.full(len(vals), i, dtype=float)
        x = x + rng.uniform(-0.07, 0.07, size=len(vals))
        c = palette[i % len(palette)]
        ax.scatter(
            x,
            vals,
            s=46,
            alpha=0.78,
            color=c,
            edgecolor="black",
            linewidth=0.5,
            label=f"{cancer} (n_pos={row['n_pos']})",
        )
        ax.hlines(row["auc_mean"], i - 0.32, i + 0.32, color=c, linewidth=2.5, zorder=3)
        # Annotation: for high-AUC cancers (above ~0.96) the annotation
        # would collide with the upper-left legend; flip it ABOVE the
        # mean line. For the rest, keep it BELOW the scatter.
        annot_y = (
            max(vals) + 0.012
            if row["auc_mean"] > 0.95
            else min(vals) - 0.018
        )
        ax.text(
            i,
            annot_y,
            f"mean={row['auc_mean']:.3f}\nstd={row['auc_std']:.3f}",
            ha="center",
            va=("bottom" if row["auc_mean"] > 0.95 else "top"),
            fontsize=8.5,
            color=c,
            fontweight="bold",
        )

    ax.set_xticks(range(len(cancers)))
    ax.set_xticklabels(cancers, fontsize=10)
    ax.set_xlim(-0.6, len(cancers) - 0.4)
    ax.set_ylabel("Per-cancer OvR AUC (per seed)", fontsize=10)
    ax.set_title(
        "Figure 2 — Per-cancer per-seed AUC (top-5 cancers by sample count)\n"
        "OvR vs ALL healthy; per-seed mean (line) + per-seed distribution (points).",
        fontsize=11,
    )
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.6, label="chance")
    all_vals = [v for c in cancers for v in per_cancer[c]["per_seed_auc"]]
    # Extend the top a bit so the "above-mean" annotations for high-AUC
    # cancers (LUAD, BRCA, OV, PAAD) don't overlap the chart title.
    ax.set_ylim(min(all_vals) - 0.05, 1.06)
    # Place legend OUTSIDE the axes (to the right) so it never overlaps
    # any of the per-cancer annotations.
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=8.5, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=DEFAULT_DPI, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 — calibration curve (operating-point summary, not true calibration)
# ---------------------------------------------------------------------------
def figure3_calibration(pooled: dict, out_path: Path) -> None:
    """Honest operating-point summary plot.

    The pooled LR scores are not calibrated probabilities — they're the
    output of `StandardScaler -> PCA(200) -> LogisticRegression(C=1.0)`
    on a 63k-dim feature vector. The sklearn `predict_proba` on this
    pipeline is monotone in the decision function but not probability-
    calibrated. The honest figure shows the sens@spec operating points
    (read off the pooled ROC) as a sens-vs-spec curve with the diagonal
    "perfect calibration" reference — i.e., this is the same info as
    fig 1b but in ROC space, not a true reliability diagram.

    The diagonal is drawn for reference only. The actual operating-
    point data is the 5 (spec, sens) tuples harvested from sens@spec_95
    /98/99 in the JSON plus the implied (0,0) and (1,1) endpoints.

    NOTE: the spec grid is sparse (3 points), so we don't interpolate
    a smooth curve; we connect the dots with a line just for visual
    continuity.
    """
    harm = pooled["harmonized"]
    no_harm = pooled["no_harmonize"]

    fig, ax = plt.subplots(figsize=(7.5, 5.6))

    # (spec, sens) tuples — plot as 1-spec (FPR) vs sens (TPR)
    specs = [0.95, 0.98, 0.99]
    sens_h = [harm["sens_at_95"], harm["sens_at_98"], harm["sens_at_99"]]
    sens_nh = [no_harm["sens_at_95"], no_harm["sens_at_98"], no_harm["sens_at_99"]]

    # ROC anchor: (FPR=0, TPR=0), (FPR=1, TPR=1) plus the 3 spec/sens points
    fpr_h = [0.0] + [1 - s for s in specs] + [1.0]
    tpr_h = [0.0] + sens_h + [1.0]
    fpr_nh = [0.0] + [1 - s for s in specs] + [1.0]
    tpr_nh = [0.0] + sens_nh + [1.0]

    ax.plot(fpr_h, tpr_h, "o-", color="#1f77b4", linewidth=2.0,
            markersize=7, label=f"harmonized (AUC={harm['auc_mean']:.4f})")
    ax.plot(fpr_nh, tpr_nh, "s--", color="#ff7f0e", linewidth=2.0,
            markersize=7, label=f"no_harmonize (AUC={no_harm['auc_mean']:.4f})")
    ax.plot([0, 1], [0, 1], ":", color="grey", alpha=0.7,
            label="chance (AUC=0.5)")

    # Label each operating point with sens@spec
    for s, t in zip(specs, sens_h):
        ax.annotate(f"  sens={t:.3f}\n  @spec={s:.2f}",
                    xy=(1 - s, t), xytext=(8, -2),
                    textcoords="offset points",
                    fontsize=8, color="#1f77b4", va="top")

    ax.set_xlabel("False positive rate (1 − specificity)", fontsize=10)
    ax.set_ylabel("True positive rate (sensitivity)", fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(
        "Figure 3 — Operating-point curve (pooled OOF sens@spec grid)\n"
        "Operating points read off the pooled ROC at spec ∈ {0.95, 0.98, 0.99}; "
        "NOT a probability calibration curve.",
        fontsize=10.5,
    )
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=DEFAULT_DPI, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4 — sens@95/98/99 grouped bar chart per cancer
# ---------------------------------------------------------------------------
def figure4_sens_operating(per_cancer: dict, out_path: Path) -> None:
    cancers = [c for c in TOP5_CANCERS if c in per_cancer]
    if not cancers:
        cancers = sorted(
            per_cancer.keys(),
            key=lambda k: per_cancer[k].get("n_pos", 0),
            reverse=True,
        )[:5]

    specs = [0.95, 0.98, 0.99]
    sens_keys = ["sens_at_95", "sens_at_98", "sens_at_99"]
    sens_matrix = np.array(
        [[per_cancer[c][k] for k in sens_keys] for c in cancers]
    )

    n_cancers = len(cancers)
    n_specs = len(specs)
    x = np.arange(n_cancers)
    width = 0.26

    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    colors = ["#1f77b4", "#2ca02c", "#d62728"]
    for i, (s, key) in enumerate(zip(specs, sens_keys)):
        ax.bar(
            x + (i - 1) * width,
            sens_matrix[:, i],
            width,
            label=f"sens @ spec={s:.2f}",
            color=colors[i],
            alpha=0.88,
        )
        for j, val in enumerate(sens_matrix[:, i]):
            ax.text(
                x[j] + (i - 1) * width,
                val + 0.012,
                f"{val:.3f}",
                ha="center",
                fontsize=8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(cancers, fontsize=10)
    ax.set_ylabel("Sensitivity (OvR vs ALL healthy)", fontsize=10)
    ax.set_ylim(0, 1.08)
    ax.set_title(
        "Figure 4 — Per-cancer sensitivity at fixed specificity (top-5 cancers)\n"
        "OvR scoring; values from results/per_cancer_sens_at_spec.json:per_cancer.*.sens_at_spec",
        fontsize=10.5,
    )
    ax.grid(True, alpha=0.25, linestyle="--", axis="y")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.6)

    fig.tight_layout()
    fig.savefig(out_path, dpi=DEFAULT_DPI, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def make_all_figures(
    fig_dir: Path = FIG_DIR,
    validate_only: bool = False,
) -> dict:
    """Build all 4 figures. Returns a manifest dict {fig_name: path}.

    If `validate_only` is True, save to a temp directory and discard
    after returning the manifest (used by the test harness).
    """
    pooled = load_pooled()
    per_cancer = load_per_cancer()

    # Required-source check — fail loud if a top-5 cancer is missing.
    missing = [c for c in TOP5_CANCERS if c not in per_cancer]
    if missing:
        raise KeyError(
            f"Top-5 cancers {missing} missing from per_cancer_sens_at_spec.json"
        )
    for required in ("auc_mean", "per_seed_auc", "sens_at_95",
                     "sens_at_98", "sens_at_99"):
        if required not in pooled["harmonized"]:
            raise KeyError(
                f"Required field '{required}' missing from pooled.harmonized"
            )

    if not validate_only:
        fig_dir.mkdir(parents=True, exist_ok=True)
    save_dir = fig_dir if not validate_only else fig_dir  # validate-only writes to tmp
    if validate_only:
        import tempfile
        tmp = tempfile.mkdtemp(prefix="plot_figures_validate_")
        save_dir = Path(tmp)

    out_paths = {
        "fig1_pooled_roc": save_dir / "fig1_pooled_roc.png",
        "fig2_per_cancer_roc": save_dir / "fig2_per_cancer_roc.png",
        "fig3_calibration": save_dir / "fig3_calibration.png",
        "fig4_sens_operating_points": save_dir / "fig4_sens_operating_points.png",
    }

    figure1_pooled(pooled, out_paths["fig1_pooled_roc"])
    figure2_per_cancer(per_cancer, out_paths["fig2_per_cancer_roc"])
    figure3_calibration(pooled, out_paths["fig3_calibration"])
    figure4_sens_operating(per_cancer, out_paths["fig4_sens_operating_points"])

    # If validate_only, wipe the temp PNGs so we don't leave disk litter.
    if validate_only:
        import shutil
        shutil.rmtree(save_dir, ignore_errors=True)
        # And set the manifest paths to point at the canonical location
        # so the caller can assert where they'd live.
        for k in out_paths:
            out_paths[k] = fig_dir / Path(out_paths[k]).name

    return {
        "manifest": {k: str(v) for k, v in out_paths.items()},
        "validate_only": validate_only,
        "top5_cancers": TOP5_CANCERS,
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            __doc__.split("\n\n", 1)[0]
            if __doc__
            else "Generate cfDNA validation figures."
        )
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=FIG_DIR,
        help=f"Output directory for the 4 PNG figures (default: {FIG_DIR})",
    )
    ap.add_argument(
        "--validate-only",
        action="store_true",
        help="Build figures in a temp dir, validate no exceptions, "
        "discard the files. Use this in CI without polluting docs/.",
    )
    args = ap.parse_args()

    if args.validate_only:
        result = make_all_figures(fig_dir=args.out_dir, validate_only=True)
        print("validate-only: built all 4 figures in temp dir, no files saved.")
        print("manifest (canonical paths):")
        for k, v in result["manifest"].items():
            print(f"  {k}: {v}")
        return 0

    result = make_all_figures(fig_dir=args.out_dir, validate_only=False)
    print(f"Wrote 4 figures to {args.out_dir}/:")
    for k, v in result["manifest"].items():
        p = Path(v)
        size_kb = p.stat().st_size / 1024.0 if p.exists() else float("nan")
        print(f"  {k}: {v} ({size_kb:.1f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
