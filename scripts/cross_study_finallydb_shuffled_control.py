#!/usr/bin/env python3
"""Shuffled-label null control for the cross-study FinaleDB benchmark.

This is a focused, faster re-runner for ONLY the null control — it
skips the main 5-seed x 5-fold per-cancer OvR sweep (which takes 5–10
minutes) and runs just the label-permuted pooled OOF. Same pipeline
(harmonize + PCA + LR) and same 5 seeds as `cross_study_finallydb.py`.

Mirrors the Audit-2 P0-C fix pattern in
`scripts/foundation_real_smoke.py:580-590`: label permutation once,
outside the fold loop, so every fold sees the same permutation; the
classifier is trained on shuffled labels; AUC is computed against the
shuffled labels.

The control is a sanity check, not a gate:
  shuffled_AUC < 0.55 → null control PASSED → the main 0.97 AUC is a
                          cancer signal, not a batch-leakage artifact
                          via fold structure.
  shuffled_AUC > 0.55 → WARNING printed; the run does NOT fail (a noisy
                          seed can fluctuate just above the floor;
                          such a finding would be reported honestly).

Output:
  results/cross_study_finallydb_shuffled_control.json — schema mirrors
  the top-level fields of `cross_study_finallydb.json`
  (pooled_auc_mean, pooled_auc_std, pooled_sens_at_95,
  pooled_sens_at_99, seeds, schema_version, generated_at) plus
  `control_passed: bool`.

Usage:
  env -u PYTHONPATH .venv/bin/python \\
      scripts/cross_study_finallydb_shuffled_control.py
"""
from __future__ import annotations

import argparse
import os
import sys

# Make the repo's src/ and scripts/ importable, and the cfdna pipeline.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))

from cross_study_finallydb import (  # noqa: E402
    DEFAULT_PCA_N,
    DEFAULT_SEEDS,
    PUBLICATION_REGISTRY,
    SHUFFLED_CONTROL_RUN_SEED,
    apply_cell_line_filter,
    load_labels_multiclass,
    main_with_shuffled,
)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--features-dir",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/data/features")
    ap.add_argument("--labels-multiclass",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv")
    ap.add_argument("--publications", nargs="+", default=["6", "8"],
                    help="Space-separated FinaleDB publication ids to "
                         "include (default: '6 8' = Jiang + Cristiano, "
                         "the locally-cached pair).")
    ap.add_argument("--include-snyder", action="store_true")
    ap.add_argument("--include-sun", action="store_true")
    ap.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    ap.add_argument("--pca", type=int, default=DEFAULT_PCA_N,
                    help="PCA components per fold (default 200 = same "
                         "as the main benchmark).")
    ap.add_argument("--out-json",
                    default="results/cross_study_finallydb_shuffled_control.json")
    ap.add_argument("--perm-seed", type=int,
                    default=SHUFFLED_CONTROL_RUN_SEED,
                    help="Global seed for the label permutation. "
                         f"Default {SHUFFLED_CONTROL_RUN_SEED}.")
    args = ap.parse_args()

    # Build a single args-like object that exposes both the names this
    # script accepts (--features-dir, --perm-seed, --out-json) AND the
    # names the helper inside cross_study_finallydb.py expects
    # (--shuffled-label-control-out, --shuffled-perm-seed). The helper
    # only reads `.features_dir`, `.labels_multiclass`, `.publications`,
    # `.include_snyder`, `.include_sun`, `.seeds`, `.pca`,
    # `.shuffled_label_control_out`, `.shuffled_perm_seed`.
    args.shuffled_label_control_out = args.out_json
    args.shuffled_perm_seed = args.perm_seed

    # Sanity-check: ensure args.publications is in the registry.
    known = set(PUBLICATION_REGISTRY.keys())
    unknown = {str(p) for p in args.publications} - known
    if unknown:
        ap.error(
            f"unknown publication id(s): {sorted(unknown)}. "
            f"Known: {sorted(known)}"
        )
    sys.exit(main_with_shuffled(args))


if __name__ == "__main__":
    main()
