#!/usr/bin/env python3
"""Check every number in the three candidate paragraphs against the artifacts.

What this is
------------
Three candidate paragraphs are drafted for the submission manuscript: the
PHTS screening-cohort null, the sequencing-depth control for the
held-out publications, and the SAUNA tumour-fraction gradient. Each one
is one edit away from `paper/paper.tex`, and each number in them has to
be the artifact's number and not a remembered one.

The paragraphs themselves live in the author's private working sheet
(`paper/BIORXIV_SUBMISSION.md` §8, which is gitignored), so this script
carries the numbers it checks: the check is the point, and a reader can
run it and see that every drafted figure is the committed artifact's
figure. When the paragraphs are applied, the repository's own
figure guard takes over -- this script is the check that the text is
right BEFORE it is applied.

Run from the repository root:

    python3 scripts/check_candidate_paragraph_numbers.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _load(name: str) -> dict:
    return json.loads((ROOT / "results" / f"{name}.json").read_text(
        encoding="utf-8"))


def _close(a: float, b: float, tol: float = 5e-5) -> bool:
    return abs(float(a) - float(b)) <= tol


def main() -> int:
    problems: list[str] = []
    checked = 0

    def check(label: str, got, want, tol: float = 5e-5) -> None:
        nonlocal checked
        checked += 1
        if got is None:
            problems.append(f"{label}: the artifact has no value there")
            return
        if not _close(got, want, tol):
            problems.append(f"{label}: artifact {got} vs drafted {want}")

    # ---- item 4: the PHTS null -------------------------------------------
    phts = _load("external_cohort_zenodo_phts")
    check("phts transferred auc", phts["result"]["auc"], 0.4506)
    check("phts transferred ci_lo", phts["result"]["ci_lo"], 0.3343)
    check("phts transferred ci_hi", phts["result"]["ci_hi"], 0.5669)
    check("phts within-cohort auc",
          phts["within_cohort_model"]["all_cancer_vs_no_cancer"]["auc"], 0.5465)
    check("phts within-cohort ci_lo",
          phts["within_cohort_model"]["all_cancer_vs_no_cancer"]["ci_lo"],
          0.4319)
    check("phts within-cohort ci_hi",
          phts["within_cohort_model"]["all_cancer_vs_no_cancer"]["ci_hi"],
          0.6612)
    check("phts no-self-normalisation auc",
          phts["harmonization_sensitivity"]["train_global_test_train_statistics"]
          ["auc"], 0.5388)
    check("phts no-self-normalisation ci_lo",
          phts["harmonization_sensitivity"]["train_global_test_train_statistics"]
          ["ci_lo"], 0.4232)
    check("phts no-self-normalisation ci_hi",
          phts["harmonization_sensitivity"]["train_global_test_train_statistics"]
          ["ci_hi"], 0.6544)
    short_long = phts["single_channel_readout"]["genome_short_long_ratio"]
    check("phts short/long ratio auc", short_long["post_smn"]["auc"], 0.6455)
    check("phts short/long ratio ci_lo",
          short_long["post_smn"]["ci_lo"], 0.5008)
    check("phts short/long ratio ci_hi",
          short_long["post_smn"]["ci_hi"], 0.7901)
    check("phts cohort n", phts["cohort"]["n"], 99)
    check("phts cohort n_cancer", phts["cohort"]["n_cancer"], 49)
    check("phts cohort n_healthy", phts["cohort"]["n_healthy"], 50)

    # ---- item 5: the depth control on the held-out publications ----------
    ext = _load("external_validation_finaledb")
    depth = ext.get("depth_control", {})
    for publication, depth_alone, depth_removed in (("7", 0.5342, 0.7436),
                                                   ("1", 0.5271, 0.6718)):
        check(f"depth alone pub {publication}",
              depth[publication]["auc_depth_alone"], depth_alone)
        check(f"depth removed pub {publication}",
              depth[publication]["auc_model_with_depth_removed"],
              depth_removed)

    # ---- item 3: the ccRCC cohort and its depth confound ----------------
    ccrcc = _load("external_cohort_zenodo_ccrcc")
    check("ccrcc auc", ccrcc["result"]["auc"], 0.9444)
    check("ccrcc ci_lo", ccrcc["result"]["ci_lo"], 0.8156)
    check("ccrcc ci_hi", ccrcc["result"]["ci_hi"], 1.0)
    check("ccrcc exact permutation p",
          ccrcc["result"]["exact_permutation_p"], 0.004329, 1e-5)
    check("ccrcc depth alone", ccrcc["depth_control"]["auc_depth_alone"],
          0.8333)
    check("ccrcc depth removed",
          ccrcc["depth_control"]["auc_model_with_depth_removed"], 0.9722)

    # ---- item 6: the SAUNA gradient --------------------------------------
    sauna = _load("external_cohort_zenodo_sauna")
    dose = sauna["result"]["dose_response"]
    check("sauna rho all", dose["rho_model_vs_tumour_fraction"], 0.6611)
    check("sauna rho all p", dose["permutation_p"], 0.064147, 1e-5)
    check("sauna rho crc", dose["rho_crc_only"], 0.8154)
    check("sauna rho crc p", dose["permutation_p_crc_only"], 0.035548, 1e-5)
    check("sauna heldout short/long",
          sauna["cohort"]["library"]["heldout_median_short_long_ratio"], 0.0849)
    check("sauna training short/long",
          sauna["cohort"]["library"]["training_median_short_long_ratio"],
          0.2182)
    check("sauna cohort n", sauna["cohort"]["n"], 9)
    check("sauna cohort n_crc", sauna["cohort"]["n_crc"], 7)
    check("sauna longitudinal monotone",
          float(sauna["longitudinal_series"]["short_long_ratio_monotone"]),
          1.0)
    series = sauna["longitudinal_series"]
    draws = series["draws"]
    check("sauna longitudinal short/long 1",
          draws[0]["short_long_ratio"], 0.20709, 1e-4)
    check("sauna longitudinal short/long 2",
          draws[1]["short_long_ratio"], 0.23023, 1e-4)
    check("sauna longitudinal short/long 3",
          draws[2]["short_long_ratio"], 0.26068, 1e-4)
    check("sauna longitudinal tumour fraction 1",
          draws[0]["tumour_fraction"], 0.2903, 1e-4)
    check("sauna longitudinal tumour fraction 2",
          draws[1]["tumour_fraction"], 0.3271, 1e-4)
    check("sauna longitudinal tumour fraction 3",
          draws[2]["tumour_fraction"], 0.4029, 1e-4)

    print(f"checked {checked} drafted numbers against the artifacts")
    for problem in problems:
        print(f"  MISMATCH: {problem}")
    if problems:
        print(f"{len(problems)} drafted numbers do not match the artifacts")
        return 1
    print("every drafted number is the artifact's number")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
