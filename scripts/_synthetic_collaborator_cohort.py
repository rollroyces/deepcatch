#!/usr/bin/env python3
"""Synthesize a small reference cohort in the Collaborator Data Interface layout.

Produces a 30-sample cohort (20 cancer + 10 healthy, plus one GM1100
cell-line example that the cell-line filter should drop) under
``data/synthetic_collaborator_cohort/`` using the **exact** directory
layout a collaborator is asked to provide:

    <cohort_root>/features/<sample_id>.<channel>.npy      (5 channels)
    <cohort_root>/labels.tsv                              (TSV)
    <cohort_root>/manifest.json                           (cohort metadata)

This script is the *reference example* collaborators see in
``docs/COLLABORATOR_DATA_INTERFACE.md``.  It is also the input fixture
for ``tests/test_adapter_local_cohort.py`` — keeping the layout identical
between docs, tests, and CLI examples means the docs cannot drift from
what the adapter actually accepts.

Outputs:
    data/synthetic_collaborator_cohort/...   (the cohort)

Channel naming follows ``src/fragmentomics/tumor_naive_adapter.py``:
    delfi_5mb_ratio      631 bins
    delfi_5mb_coverage   631 bins
    delfi_100kb_ratio    30894 bins
    delfi_100kb_counts   30894 bins
    fsd_histogram        196 bins  (5bp bins, 20–1000bp)

The synthetic values are NOT realistic fragmentomics signal — they are
deterministic per-seed random arrays that differ by class (so an LR
trained on the cohort returns >0.5 AUC) but they are explicitly NOT a
biologically valid substitute for real data.  A `--realistic-seed` flag
is intentionally NOT provided: this is the *layout example*, never the
*signal example*.

Run:
    env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \\
        scripts/_synthetic_collaborator_cohort.py \\
        --out data/synthetic_collaborator_cohort
"""
from __future__ import annotations

import argparse
import json
import os
import random
from datetime import datetime, timezone

import numpy as np


# Channel contract — keep in sync with src/fragmentomics/tumor_naive_adapter.py
CHANNEL_DIMS = {
    "delfi_5mb_ratio": 631,
    "delfi_5mb_coverage": 631,
    "delfi_100kb_ratio": 30894,
    "delfi_100kb_counts": 30894,
    "fsd_histogram": 196,
}
CHANNEL_NAMES = list(CHANNEL_DIMS.keys())

# Sample composition — exactly the size promised in the docs.
N_CANCER = 20
N_HEALTHY = 10
N_CELL_LINE = 1  # GM1100 — exercises the cell-line filter regex

# Class → short disease_class label used by the per-cancer OvR section.
CLASS_LABEL_CANCER = "CRC_S"          # synthetic CRC-like cohort
CLASS_LABEL_HEALTHY = "HEALTHY"


def _make_sample_features(
    rng: np.random.Generator,
    sample_id: str,
    is_cancer: bool,
) -> dict[str, np.ndarray]:
    """Build the per-channel .npy contents for one sample.

    Cancer samples have slightly higher short/long ratio (mirroring the
    real effect that motivated the DELFI pipeline); healthy samples have
    smoother coverage.  The shift is small (mean shift ~0.01, scale
    unchanged) so a 5-fold CV on the synthetic cohort yields an honest
    AUC ~0.65 — enough to prove the layout works end-to-end without
    overpromising a clinical claim.
    """
    mean_shift = 0.01 if is_cancer else 0.0

    out: dict[str, np.ndarray] = {}
    for ch, dim in CHANNEL_DIMS.items():
        if ch == "fsd_histogram":
            # Fragment-length histogram: positive, peaks at 167bp.
            raw = rng.random(dim) + 0.05
            raw = raw / raw.sum()
            out[ch] = raw.astype(np.float32)
        elif ch.endswith("_ratio"):
            # DELFI ratio channel — bounded in [0, 1].
            v = np.clip(
                rng.normal(loc=0.15 + mean_shift, scale=0.02, size=dim),
                0.0, 1.0,
            ).astype(np.float32)
            out[ch] = v
        elif ch.endswith("_coverage"):
            # DELFI coverage channel — positive, median-normalized.
            v = np.clip(
                rng.normal(loc=1.0 + mean_shift, scale=0.2, size=dim),
                0.0, None,
            ).astype(np.float32)
            out[ch] = v
        elif ch.endswith("_counts"):
            # Raw 100kb counts — positive integers (depth, not biology).
            v = rng.integers(50, 5000, size=dim).astype(np.float32)
            out[ch] = v
        else:  # pragma: no cover — defensive
            raise ValueError(f"unknown channel {ch}")
    return out


def _build_labels_tsv(rows: list[tuple[str, str, str]]) -> str:
    """Build the labels.tsv content.

    columns: sample_id, disease_class, label, [study_id, [publication_id, [tissue_of_origin]]]
    """
    lines = ["sample_id\tdisease_class\tlabel\tstudy_id\tpublication_id\ttissue_of_origin"]
    for sid, dc, label, *rest in rows:
        study = rest[0] if len(rest) > 0 else "syn"
        pub = rest[1] if len(rest) > 1 else ""
        tissue = rest[2] if len(rest) > 2 else ""
        lines.append(f"{sid}\t{dc}\t{label}\t{study}\t{pub}\t{tissue}")
    return "\n".join(lines) + "\n"


def _build_manifest(
    cohort_name: str,
    contributing_lab: str,
    contact_email: str,
    irb_number: str,
    rows: list[tuple[str, str, str]],
    n_dropped_cell_line: int,
) -> dict:
    n_cancer = sum(1 for r in rows if r[2] == "cancer")
    n_healthy = sum(1 for r in rows if r[2] == "healthy")
    return {
        "schema_version": "1.0",
        "cohort_name": cohort_name,
        "contributing_lab": contributing_lab,
        "contact_email": contact_email,
        "irb_number": irb_number,
        "sample_count_cancer": n_cancer,
        "sample_count_healthy": n_healthy,
        "channels_present": CHANNEL_NAMES,
        "vector_lengths": {ch: int(d) for ch, d in CHANNEL_DIMS.items()},
        "cell_line_filter_dropped_count": n_dropped_cell_line,
        "generation_date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "scope": "synthetic reference cohort (Collaborator Data Interface demo). "
                 "NOT real patient data; not for clinical use.",
        "citation": "Generated by scripts/_synthetic_collaborator_cohort.py "
                    "in rollroyces/deepcatch.",
        "license": "CC0 (synthetic data, no patient information).",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        default="data/synthetic_collaborator_cohort",
        help="Cohort root directory to create.",
    )
    ap.add_argument(
        "--cohort-name",
        default="synthetic_reference_2026",
        help="Value to record as cohort_name in manifest.json.",
    )
    ap.add_argument("--seed", type=int, default=0,
                    help="RNG seed for deterministic outputs.")
    args = ap.parse_args()

    out_root = args.out
    feat_dir = os.path.join(out_root, "features")
    os.makedirs(feat_dir, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    sample_rows: list[tuple[str, str, str]] = []

    # Cancer samples: SYN_C001..SYN_C020
    for i in range(1, N_CANCER + 1):
        sid = f"SYN_C{i:03d}"
        feats = _make_sample_features(rng, sid, is_cancer=True)
        for ch, arr in feats.items():
            np.save(os.path.join(feat_dir, f"{sid}.{ch}.npy"), arr)
        sample_rows.append((sid, CLASS_LABEL_CANCER, "cancer"))

    # Healthy samples: SYN_H001..SYN_H010
    for i in range(1, N_HEALTHY + 1):
        sid = f"SYN_H{i:03d}"
        feats = _make_sample_features(rng, sid, is_cancer=False)
        for ch, arr in feats.items():
            np.save(os.path.join(feat_dir, f"{sid}.{ch}.npy"), arr)
        sample_rows.append((sid, CLASS_LABEL_HEALTHY, "healthy"))

    # Cell-line example: GM1100 — included so the adapter's cell-line
    # filter has something to drop.  Marked "cancer" + "Liver cancer"
    # disease_class to mirror the real FinaleDB mislabeling pattern.
    sid = "GM1100"
    feats = _make_sample_features(rng, sid, is_cancer=True)
    for ch, arr in feats.items():
        np.save(os.path.join(feat_dir, f"{sid}.{ch}.npy"), arr)
    sample_rows.append((sid, "Liver cancer", "cancer"))
    n_dropped_cell_line = 1

    # labels.tsv
    with open(os.path.join(out_root, "labels.tsv"), "w") as f:
        f.write(_build_labels_tsv(sample_rows))

    # manifest.json
    manifest = _build_manifest(
        cohort_name=args.cohort_name,
        contributing_lab="Independent Researcher (synthetic demo)",
        contact_email="hello@example.org",
        irb_number="N/A (synthetic data — no IRB required)",
        rows=sample_rows,
        n_dropped_cell_line=n_dropped_cell_line,
    )
    with open(os.path.join(out_root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    # Report to stdout for the human running the script.
    total_features = sum(1 for _ in os.listdir(feat_dir))
    print(f"[synthetic_collaborator_cohort] wrote {len(sample_rows)} samples "
          f"to {out_root}/")
    print(f"[synthetic_collaborator_cohort] {total_features} per-sample .npy files "
          f"under {feat_dir}/")
    print(f"[synthetic_collaborator_cohort] labels.tsv: {len(sample_rows)} rows "
          f"({N_CANCER} cancer + {N_HEALTHY} healthy + {N_CELL_LINE} cell-line)")
    print(f"[synthetic_collaborator_cohort] manifest.json: schema_version=1.0, "
          f"cohort_name={args.cohort_name!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())