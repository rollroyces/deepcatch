#!/usr/bin/env python3
"""Standalone per-cancer sens@spec + PPV@prev table generator.

The standard cfDNA clinical-validation table. For each cancer type
(or pooled across all), this emits:

  - AUC with a DeLong 95% CI
  - Sensitivity at specificity ∈ {0.95, 0.98, 0.99} with DeLong 95% CIs
  - PPV at prevalence ∈ {0.001, 0.004, 0.01, 0.05, 0.10, 0.20, 0.50}
    evaluated at the spec=0.99 operating point

Usage
-----
Read a TSV with columns (sample_id, score, y, cancer_label, study):
    python scripts/per_cancer_sens_at_spec.py --scores-tsv path/to/scores.tsv \\
        --out results/per_cancer_sens_at_spec.json

Generate a synthetic fixture for smoke testing:
    python scripts/per_cancer_sens_at_spec.py --synthetic --n 200 \\
        --out results/per_cancer_sens_at_spec_synth.json

This script is deliberately framework-light — numpy + scipy only at the
math layer, pandas-free so it runs anywhere. It is the same math as
the in-repo per-cancer panel-LLR smokes (see sens_at_spec_ce.json /
sens_at_spec_sens.json shape) but standalone: it works on ANY
(panel_scores, frag_scores, cancer_labels, study_labels) inputs.
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import List, Tuple

import numpy as np

import os
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.per_cancer_sens_at_spec import (  # noqa: E402
    DEFAULT_PREVALENCES,
    DEFAULT_SPECIFICITIES,
    PPV_AT_SPEC,
    MIN_POSITIVES_FOR_CI,
    SCREENING_SPECIFICITIES,
    _spec_field_name,
    build_per_cancer_table,
)


def _spec_col_header(spec: float) -> str:
    """Compact console header for a specificity column. Uses 5 chars.

    Examples: 0.95 → ``s 95``, 0.99 → ``s 99``, 0.995 → ``s99.5``,
    0.999 → ``s99.9``.
    """
    if abs(spec - round(spec, 2)) < 1e-9:
        return f"s{int(round(spec * 100)):>4}"
    return f"@{int(round(spec * 1000)):>4}"


def _format_sens_cell(v) -> str:
    """Format a sensitivity value as a 5-char cell.

    ``None`` / NaN → ``"  NA "``. Finite float → ``"<v>"`` right-aligned.
    """
    if v is None:
        return "  NA "
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "  NA "
    import math
    if math.isnan(f) or math.isinf(f):
        return "  NA "
    return f"{f:5.3f}"


# ──────────────────────────────────────────────────────────────────────
# Synthetic fixture (for --synthetic smoke)
# ──────────────────────────────────────────────────────────────────────

def build_synthetic_fixture(
    n: int,
    seed: int,
    n_cancer_types: int = 4,
    per_cancer_positives: int = 30,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build a deterministic synthetic cfDNA scoring fixture.

    The scores are drawn from a per-cancer-type shifted normal so each
    cancer has a measurably different AUC vs the healthy baseline.

    Returns
    -------
    (y, scores, cancer_label, study_label)
        y            : 1=cancer / 0=control
        scores       : continuous "cancer-likeness" score (higher = cancer)
        cancer_label : per-sample cancer type or "HEALTHY"
        study_label  : per-sample study identifier
    """
    rng = np.random.default_rng(seed)
    # 80 healthy + n_cancer_types * per_cancer_positives
    n_healthy = n - n_cancer_types * per_cancer_positives
    if n_healthy < 20:
        raise ValueError(
            f"n={n} too small for {n_cancer_types} cancers × "
            f"{per_cancer_positives} positives (need n_healthy >= 20)."
        )
    cancers = ["BRCA", "LUAD", "CRC", "PAAD", "OV", "HCC", "ESCA"][:n_cancer_types]
    cancer_label = (
        ["HEALTHY"] * n_healthy
        + sum(([c] * per_cancer_positives for c in cancers), [])
    )
    y = np.array([0] * n_healthy + [1] * (n - n_healthy), dtype=np.int64)
    # Score: cancer mean = 0.65 + per-type bonus (0.0..0.30);
    # healthy mean = 0.30; std = 0.15.
    scores = np.empty(n, dtype=np.float64)
    for i, label in enumerate(cancer_label):
        if label == "HEALTHY":
            scores[i] = rng.normal(0.30, 0.15)
        else:
            bonus = (cancers.index(label) + 1) / (n_cancer_types + 1) * 0.30
            scores[i] = rng.normal(0.65 + bonus, 0.12)
    # Two synthetic studies to exercise n_studies provenance.
    study_label = np.array(
        ["S1"] * (n // 2) + ["S2"] * (n - n // 2), dtype=object
    )
    return y, scores, np.array(cancer_label), study_label


# ──────────────────────────────────────────────────────────────────────
# TSV I/O
# ──────────────────────────────────────────────────────────────────────

REQUIRED_TSV_COLUMNS = ("sample_id", "score", "y", "cancer_label")


def load_scores_tsv(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Load a scores TSV with required columns.

    Expected columns: sample_id, score, y, cancer_label [, study, stage].

    Returns
    -------
    (y, score, cancer_label, study_label_or_None, stage_label_or_None)
        ``stage_label_or_None`` is a per-sample stage string (e.g.
        ``"I"``, ``"II"``, ``"III"``, ``"IV"``, ``"NA"``) used by the
        stage-aware breakdown. ``None`` when the TSV has no ``stage``
        column.
    """
    import csv

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames or []
        missing = [c for c in REQUIRED_TSV_COLUMNS if c not in cols]
        if missing:
            raise ValueError(
                f"scores TSV {path!r} missing required columns: {missing}; "
                f"got {cols}."
            )
        rows = list(reader)
    if not rows:
        raise ValueError(f"scores TSV {path!r} is empty.")
    sample_ids = [r["sample_id"] for r in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError(
            f"scores TSV {path!r} has duplicate sample_id values."
        )
    y = np.array([int(r["y"]) for r in rows], dtype=np.int64)
    scores = np.array([float(r["score"]) for r in rows], dtype=np.float64)
    cancer_label = np.array([r["cancer_label"] for r in rows], dtype=object)
    study_label = (
        np.array([r["study"] for r in rows], dtype=object)
        if "study" in cols else None
    )
    stage_label = (
        np.array([r["stage"] for r in rows], dtype=object)
        if "stage" in cols else None
    )
    return y, scores, cancer_label, study_label, stage_label


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--scores-tsv",
        default=None,
        help="TSV with columns: sample_id, score, y, cancer_label [, study, stage]. "
             "Mutually exclusive with --synthetic.",
    )
    ap.add_argument(
        "--synthetic",
        action="store_true",
        help="Build a deterministic synthetic fixture (smoke mode).",
    )
    ap.add_argument(
        "--n",
        type=int,
        default=200,
        help="[synthetic] Total samples (default 200).",
    )
    ap.add_argument(
        "--seed",
        type=int,
        default=42,
        help="[synthetic] RNG seed (default 42).",
    )
    ap.add_argument(
        "--n-cancer-types",
        type=int,
        default=4,
        help="[synthetic] Number of cancer types (default 4).",
    )
    ap.add_argument(
        "--per-cancer-positives",
        type=int,
        default=30,
        help="[synthetic] Positives per cancer type (default 30).",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output JSON path. Defaults to "
             "results/per_cancer_sens_at_spec.json (or "
             "..._synth.json when --synthetic).",
    )
    ap.add_argument(
        "--no-pooled",
        action="store_true",
        help="Omit the POOLED row from the output table.",
    )
    ap.add_argument(
        "--specificities",
        type=float,
        nargs="+",
        default=list(DEFAULT_SPECIFICITIES) + list(SCREENING_SPECIFICITIES),
        help=(
            "Operating-point specificities for the sens@spec grid. "
            "Default: %(default)s. Pass --specificities 0.95 0.99 0.995 "
            "to use a custom grid (e.g. for screening-grade analyses)."
        ),
    )
    ap.add_argument(
        "--include-stage-breakdown",
        action="store_true",
        help=(
            "When the input TSV has a ``stage`` column (or when "
            "--synthetic includes stages), also emit a stage-aware "
            "breakdown (Stage I / LATE / UNKNOWN) under "
            "``stage_breakdown``. Useful for ultra-early readiness "
            "audits where Stage I detection is the headline."
        ),
    )
    args = ap.parse_args()

    if (args.scores_tsv is None) ^ args.synthetic:
        ap.error(
            "Exactly one of --scores-tsv or --synthetic is required."
        )

    stage_label = None
    if args.synthetic:
        y, scores, cancer_label, study_label = build_synthetic_fixture(
            n=args.n,
            seed=args.seed,
            n_cancer_types=args.n_cancer_types,
            per_cancer_positives=args.per_cancer_positives,
        )
        if args.include_stage_breakdown:
            # Deterministic synthetic stage assignment: every other
            # cancer sample is Stage I (so the synthetic fixture
            # exercises both Stage I and LATE buckets). Healthy
            # samples are STAGE_UNKNOWN — they are not consumed by
            # the stage breakdown.
            n_healthy = int((y == 0).sum())
            cancer_stages = []
            for _ in range(args.n_cancer_types):
                interleave = (["I", "LATE"] * ((args.per_cancer_positives + 1) // 2)
                              )[:args.per_cancer_positives]
                cancer_stages.extend(interleave)
            stage_label = np.array(
                ["UNKNOWN"] * n_healthy + cancer_stages,
                dtype=object,
            )
        default_out = _REPO_ROOT / "results" / "per_cancer_sens_at_spec_synth.json"
    else:
        y, scores, cancer_label, study_label, stage_label = load_scores_tsv(
            args.scores_tsv
        )
        if args.include_stage_breakdown and stage_label is None:
            raise ValueError(
                "--include-stage-breakdown requires the input TSV to "
                "have a 'stage' column (e.g. 'I', 'II', 'NA')."
            )
        default_out = _REPO_ROOT / "results" / "per_cancer_sens_at_spec.json"

    out_path = Path(args.out) if args.out else default_out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    table = build_per_cancer_table(
        y=y,
        s=scores,
        cancer_label=cancer_label,
        study_label=study_label,
        specificities=args.specificities,
        prevalences=DEFAULT_PREVALENCES,
        include_pooled=not args.no_pooled,
    )
    # Add a provenance footer that survives into the JSON.
    table["provenance"] = {
        "source": (
            "synthetic_fixture" if args.synthetic
            else f"scores_tsv:{args.scores_tsv}"
        ),
        "n_samples": int(len(y)),
        "n_cancer_types": int(len({c for c in cancer_label if c != "HEALTHY"})),
        "specificities": list(args.specificities),
        "prevalences": list(DEFAULT_PREVALENCES),
        "ppv_at_spec": PPV_AT_SPEC,
        "min_positives_for_ci": MIN_POSITIVES_FOR_CI,
    }
    if args.include_stage_breakdown and stage_label is not None:
        from src.per_cancer_sens_at_spec import build_stage_breakdown_table
        table["stage_breakdown"] = build_stage_breakdown_table(
            y=y, s=scores, stage_label=stage_label,
            specificities=args.specificities,
            prevalences=DEFAULT_PREVALENCES,
        )

    with open(out_path, "w") as f:
        json.dump(table, f, indent=2)

    # Print a one-screen console summary so the user can see the table
    # without opening the JSON.
    print(f"Wrote {out_path} ({len(table['per_cancer'])} cancer rows + "
          f"{'pooled' if table['pooled'] else 'no pooled'})")
    specs = list(args.specificities)
    spec_cols = " ".join(_spec_col_header(s) for s in specs)
    header = f"{'cancer':<8} {'n':>5} {'n_pos':>6} {'auc':>7} " \
             f"{'auc_ci_lo':>10} {'auc_ci_hi':>10} " + spec_cols
    print(header)
    print("-" * len(header))
    for name in sorted(table["per_cancer"]):
        r = table["per_cancer"][name]
        if r["skipped"]:
            print(f"{name:<8} {r['n']:>5} {r['n_pos']:>6}   SKIP ({r['skip_reason'][:30]})")
            continue
        sens_vals = []
        for s in specs:
            key = _spec_field_name(s)
            v = r.get(key)
            sens_vals.append(_format_sens_cell(v))
        print(
            f"{name:<8} {r['n']:>5} {r['n_pos']:>6} "
            f"{r['auc_mean']:>7.3f} {r['auc_ci'][0]:>10.3f} {r['auc_ci'][1]:>10.3f} "
            + " ".join(sens_vals)
        )
    if table["pooled"] and not table["pooled"]["skipped"]:
        r = table["pooled"]
        sens_vals = []
        for s in specs:
            key = _spec_field_name(s)
            v = r.get(key)
            sens_vals.append(_format_sens_cell(v))
        print(
            f"{'POOLED':<8} {r['n']:>5} {r['n_pos']:>6} "
            f"{r['auc_mean']:>7.3f} {r['auc_ci'][0]:>10.3f} {r['auc_ci'][1]:>10.3f} "
            + " ".join(sens_vals)
        )
        # PPV @ spec=0.99 across the prevalence grid.
        print("\nPPV @ spec=0.99 across prevalences:")
        prevs = sorted(DEFAULT_PREVALENCES)
        print("  prev:", " ".join(f"{p:>7.4f}" for p in prevs))
        print("  PPV :", " ".join(
            f"{r['ppv_at_prevalence'][f'prev_{p}']:>7.4f}"
            for p in prevs
        ))
    if "stage_breakdown" in table:
        sb = table["stage_breakdown"]
        print(
            f"\nStage breakdown (n_cancer_pos={sb['n_cancer_pos']}, "
            f"Stage I={sb['n_stage_I']}, LATE={sb['n_stage_late']}, "
            f"UNKNOWN={sb['n_stage_unknown']}):"
        )
        sens_cols = " ".join(_spec_col_header(s) for s in specs)
        stage_header = f"{'stage':<8} {'n':>5} {'n_pos':>6} {'auc':>7} " + sens_cols
        print(stage_header)
        print("-" * len(stage_header))
        for stage_name in ("I", "LATE", "UNKNOWN"):
            r = sb["per_stage"].get(stage_name)
            if r is None:
                continue
            if r.get("skipped"):
                print(
                    f"{stage_name:<8} {r.get('n', 0):>5} "
                    f"{r.get('n_pos', 0):>6}   SKIP "
                    f"({r.get('skip_reason', '')[:30]})"
                )
                continue
            sens_vals = []
            for s in specs:
                key = _spec_field_name(s)
                v = r.get(key)
                sens_vals.append(_format_sens_cell(v))
            print(
                f"{stage_name:<8} {r['n']:>5} {r['n_pos']:>6} "
                f"{r['auc_mean']:>7.3f} " + " ".join(sens_vals)
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())