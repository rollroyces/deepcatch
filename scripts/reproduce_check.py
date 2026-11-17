#!/usr/bin/env python3
"""Compare a freshly generated artifact directory against the committed one.

Why this file exists
--------------------
``paper/REPRODUCE.sh`` re-runs the scripts that produce every artifact the
manuscript quotes. Re-running them is not the same as *proving* the
committed numbers are reproducible: a script can succeed and still write a
different number than the one in the paper, and nobody would notice unless
the two were compared explicitly.

This script does that comparison. It takes two directories of
``results/*.json`` -- a baseline (the committed artifacts, or a snapshot
taken before the run) and the current output -- and compares every number
the documents quote, plus the supporting numbers that appear in the
manuscript tables.

The authoritative list of quoted numbers lives in ``benchmark_numbers.py``
(``QUOTED_NUMBERS`` / ``QUOTED_SPECIAL``); this script reuses it rather
than restating it, so the two can never disagree about what is quotable.

Only *numbers* are compared. The artifacts carry a ``generated_at``
timestamp and other run metadata, so comparing files byte-for-byte would
always fail even for a perfectly reproducible pipeline.

Usage
-----
    python scripts/reproduce_check.py \
        --baseline /path/to/baseline/results \
        --current  results

Exit codes
    0 = every compared number reproduced within tolerance
    1 = at least one number drifted (or was missing on one side)
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Numbers quoted in the manuscript tables that are not in the shared
# QUOTED_NUMBERS list, but that the reproduction driver regenerates and
# that the prose depends on. (artifact, path into the artifact).
EXTRA_NUMBERS = [
    # per-cancer sens@spec standalone artifact (DeLong CIs + PPV)
    ("per_cancer_sens_at_spec.json", ("pooled", "auc_mean")),
    ("per_cancer_sens_at_spec.json", ("pooled", "sens_at_spec", "0.99",
                                      "sensitivity")),
    # headline GC-corrected run: operating points + the true-confound control
    ("cross_study_finallydb_gc_corrected.json",
     ("pooled", "harmonized", "sens_at_spec_95")),
    ("cross_study_finallydb_gc_corrected.json",
     ("pooled", "harmonized", "sens_at_spec_99")),
    ("cross_study_finallydb_gc_corrected.json",
     ("true_confound_control", "cancer_jiang_healthy_cristiano",
      "harmonized", "auc_mean")),
    # legacy run
    ("cross_study_finallydb.json",
     ("pooled", "harmonized", "sens_at_spec_99")),
    ("cross_study_finallydb.json",
     ("true_confound_control", "cancer_jiang_healthy_cristiano",
      "harmonized", "auc_mean")),
    # shuffled-label null control
    ("cross_study_finallydb_shuffled_control.json",
     ("pooled", "harmonized", "auc_mean")),
    # methylation channel: the null and its two depth/coverage baselines
    ("methylation_cohort_channel.json", ("beta_summary_auc_loo",)),
    ("methylation_cohort_channel.json", ("cpgs_covered_only_auc_loo",)),
    ("methylation_cohort_channel.json", ("obs_per_cpg_only_auc_loo",)),
    # The matched-locus test, which is what the null actually rests on.
    # These three were identical on every run. (The shared-locus AUC in
    # the same block is NOT: FinaleMe's HMM training is not deterministic,
    # so it is deliberately not compared here -- see the artifact's
    # reproducibility_note.)
    ("methylation_cohort_channel.json",
     ("common_cpg_set", "per_cpg_test", "n_cpgs_tested")),
    ("methylation_cohort_channel.json",
     ("common_cpg_set", "per_cpg_test",
      "median_beta_difference_cancer_minus_control")),
    ("methylation_cohort_channel.json",
     ("common_cpg_set", "per_cpg_test", "permutation_p_value")),
    # external validation: the training-derived operating point and the
    # case-only publication's sensitivity at it.
    ("external_validation_finaledb.json",
     ("training_derived_threshold", "threshold")),
    ("external_validation_finaledb.json",
     ("per_publication_auc", "9",
      "sensitivity_at_train_threshold_99spec")),
    # The paired synthetic-plasma honest baseline behind the methods
    # paper's section 3.1. The driver's step 14 regenerates the artifact,
    # so its numbers have to be compared here like the rest -- otherwise
    # the step would write a file nobody checks.
    ("sens_at_spec_ce_n20.json", ("panel_only_auc_mean",)),
    ("sens_at_spec_ce_n20.json", ("frag_only_auc_mean",)),
    ("sens_at_spec_ce_n20.json", ("foundation_auc_mean",)),
    ("sens_at_spec_ce_n20.json", ("foundation_auc_std",)),
    ("sens_at_spec_ce_n20.json", ("lr_baseline_auc_mean",)),
    ("sens_at_spec_ce_n20.json", ("naive_avg_auc_mean",)),
    ("sens_at_spec_ce_n20.json", ("shuffled_foundation_auc_mean",)),
    ("sens_at_spec_ce_n20.json", ("delta_auc_normalized",)),
    # The loss ablation built from the two n=20 runs.
    ("sens_at_spec_ablation_n20.json", ("foundation_auc_mean_ce",)),
    ("sens_at_spec_ablation_n20.json",
     ("foundation_auc_mean_sens_at_spec",)),
    ("sens_at_spec_ablation_n20.json", ("foundation_auc_paired_t", "p_value")),
    # The readiness verdict, which ~20 documents quote. It is derived
    # from the artifacts above, so it has to be compared too -- it used
    # to be regenerated after the proof step and never checked.
    ("ultraearly_readiness.json", ("checks", 0, "value")),
    ("ultraearly_readiness.json", ("checks", 1, "value")),
    ("ultraearly_readiness.json", ("checks", 3, "value")),
]

# Numbers whose reproduction is expected to be exact to the digit rather
# than merely close. These come from fixed per-sample matrices and a
# deterministic HMM fit, so there is no seed-to-seed variation at all.
EXACT_NUMBERS = {
    "methylation_cohort_channel.json",
}

DEFAULT_TOLERANCE = 1e-6

# Numbers that are threshold statistics rather than smooth scores. A
# sensitivity-at-specificity is a step function of the score
# distribution: the operating point is fixed by the healthy tail, so a
# change of a single healthy sample moves the threshold, and a change of
# one positive sample moves the sensitivity by 1/n_pos. These are far
# more seed-sensitive than an AUC, which is why they get their own
# tolerance band. A number is treated as a threshold statistic when its
# path mentions a sens@spec field.
THRESHOLD_MARKERS = ("sens_at_spec", "sensitivity", "sens_at_95",
                     "sens_at_98", "sens_at_99")


def _is_threshold_statistic(label: str) -> bool:
    return any(marker in label for marker in THRESHOLD_MARKERS)


def _load_numbers_module():
    spec = importlib.util.spec_from_file_location(
        "benchmark_numbers", os.path.join(ROOT, "scripts",
                                          "benchmark_numbers.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_numbers = _load_numbers_module()


def dig(artifact, path):
    node = artifact
    for key in path:
        node = node[key]
    return node


def collect(results_dir: str):
    """Return {label: (value, artifact_name)} for every comparable number."""
    out = {}
    for label, value, artifact_name, _docs in _numbers.quoted_values(
            results_dir):
        out[label] = (value, artifact_name)
    for artifact_name, path in EXTRA_NUMBERS:
        artifact = _numbers.load(artifact_name, results_dir)
        if artifact is None:
            continue
        try:
            value = float(dig(artifact, path))
        except (KeyError, TypeError):
            continue
        out[f"{artifact_name}:{'.'.join(map(str, path))}"] = (value,
                                                             artifact_name)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--baseline", required=True,
                    help="directory holding the committed results/*.json")
    ap.add_argument("--current", required=True,
                    help="directory holding the freshly generated results/*.json")
    ap.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE,
                    help="absolute tolerance for smooth numbers (AUCs)")
    ap.add_argument("--threshold-tolerance", type=float, default=None,
                    help=("absolute tolerance for threshold statistics "
                          "(sens@spec and friends); defaults to --tolerance"))
    args = ap.parse_args()
    threshold_tolerance = (args.threshold_tolerance if
                           args.threshold_tolerance is not None
                           else args.tolerance)

    base = collect(args.baseline)
    cur = collect(args.current)

    if not base:
        print("[REPRODUCE-CHECK] baseline holds no comparable numbers",
              file=sys.stderr)
        return 1
    if not cur:
        print("[REPRODUCE-CHECK] current run holds no comparable numbers",
              file=sys.stderr)
        return 1

    missing = sorted(set(base) - set(cur))
    extra = sorted(set(cur) - set(base))
    rows = []
    drift = []
    for label in sorted(set(base) & set(cur)):
        base_value, artifact_name = base[label]
        cur_value, _ = cur[label]
        if artifact_name in EXACT_NUMBERS:
            tol = 0.0
        elif _is_threshold_statistic(label):
            tol = threshold_tolerance
        else:
            tol = args.tolerance
        diff = abs(cur_value - base_value)
        ok = diff <= tol
        rows.append((label, artifact_name, base_value, cur_value, diff, ok))
        if not ok:
            drift.append(label)

    print("== reproduction check (baseline vs current)")
    print(f"   {'number':<62} {'baseline':>12} {'current':>12} {'|diff|':>10}")
    for label, _artifact, base_value, cur_value, diff, ok in rows:
        flag = "ok" if ok else "DRIFT"
        print(f"   {label:<62} {base_value:>12.6f} {cur_value:>12.6f} "
              f"{diff:>10.2e}  {flag}")

    print()
    if missing:
        print(f"[REPRODUCE-CHECK] !! {len(missing)} number(s) absent from the "
              f"current run: {', '.join(missing)}")
    if extra:
        print(f"[REPRODUCE-CHECK] note: {len(extra)} number(s) present only in "
              f"the current run: {', '.join(extra)}")
    if drift:
        print(f"[REPRODUCE-CHECK] !! {len(drift)} number(s) drifted: "
              f"{', '.join(drift)}")
        return 1
    if missing:
        # A quoted number that the run did not produce is an incomplete
        # reproduction, not a successful one.
        print("[REPRODUCE-CHECK] !! the reproduction is incomplete: some "
              "quoted artifacts were not regenerated.")
        return 1

    exact = len({a for _, a, _, _, _, _ in rows} & EXACT_NUMBERS)
    print(f"[REPRODUCE-CHECK] all {len(rows)} compared number(s) reproduced "
          f"within tolerance ({exact} artifact(s) checked exactly).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
