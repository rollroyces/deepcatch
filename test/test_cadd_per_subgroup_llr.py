"""Tests for the per-subgroup CADD panel LLR extension (DeepCatch v2.2).

Validates that:
  - cadd_per_subgroup_llr module imports cleanly
  - define_subgroups partitions the 20-patient cohort correctly by driver
    gene presence and mutation burden
  - build_uniform_weights_for_subgroup returns weight=1 for the subgroup's
    mutations and drops other patients
  - build_driver_only_weights returns weight=1 for mutations in the LUAD
    driver gene set and 0 for all others
  - run_subgroup_panel_detection preserves AUC at >= 0.5% (sanity check)
"""
import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR.parent))  # repo root (for real_tcga_validation)
sys.path.insert(0, str(SCRIPTS_DIR))

import cadd_per_subgroup_llr as cps  # noqa: E402


# A small synthetic cohort to keep tests fast.
SYNTH_COUNTS = {
    "synth_A": [
        {"gene": "TP53", "tumor_vaf": 0.3, "t_alt": 30, "t_depth": 100,
         "n_alt": 0, "n_ref": 100, "normal_error_rate": 0.001,
         "variant_class": "Missense_Mutation", "sample": "synth_A",
         "chrom": "chr17", "pos": 7674220, "ref": "G", "alt": "A"}
        for _ in range(60)
    ],
    "synth_B": [
        {"gene": "KRAS", "tumor_vaf": 0.4, "t_alt": 40, "t_depth": 100,
         "n_alt": 0, "n_ref": 100, "normal_error_rate": 0.001,
         "variant_class": "Missense_Mutation", "sample": "synth_B",
         "chrom": "chr12", "pos": 25245350, "ref": "C", "alt": "A"}
        for _ in range(80)
    ],
    "synth_C": [
        {"gene": "TTN", "tumor_vaf": 0.5, "t_alt": 50, "t_depth": 100,
         "n_alt": 0, "n_ref": 100, "normal_error_rate": 0.001,
         "variant_class": "Missense_Mutation", "sample": "synth_C",
         "chrom": "chr2", "pos": 178730000, "ref": "A", "alt": "T"}
        for _ in range(40)
    ],
}


def test_module_imports():
    assert hasattr(cps, "define_subgroups")
    assert hasattr(cps, "build_uniform_weights_for_subgroup")
    assert hasattr(cps, "build_topk_weights_for_subgroup")
    assert hasattr(cps, "build_driver_only_weights")
    assert hasattr(cps, "run_subgroup_panel_detection")


def test_define_subgroups_partitions_correctly():
    """Subgroups are defined by driver-gene presence and burden halves."""
    sg = cps.define_subgroups(SYNTH_COUNTS)
    assert set(sg.keys()) >= {
        "TP53_mutant", "TP53_wildtype", "KRAS_mutant", "KRAS_wildtype",
        "all_20_patients", "high_burden_top10", "low_burden_bottom10",
    }
    # synth_A has TP53; synth_B has KRAS; synth_C has neither
    assert sg["TP53_mutant"] == ["synth_A"]
    assert sg["TP53_wildtype"] == ["synth_B", "synth_C"]
    assert sg["KRAS_mutant"] == ["synth_B"]
    assert sg["KRAS_wildtype"] == ["synth_A", "synth_C"]
    # Mutation burden ordering: synth_B (80) > synth_A (60) > synth_C (40)
    # n=3, half=n//2=1; low_half = bottom-1 (synth_C); high_half = top-2
    # then sorted alphabetically for stable order: high = [synth_A, synth_B]
    assert sg["high_burden_top10"] == ["synth_A", "synth_B"]
    assert sg["low_burden_bottom10"] == ["synth_C"]


def test_build_uniform_weights_for_subgroup():
    sg_patients = ["synth_A", "synth_B"]
    w = cps.build_uniform_weights_for_subgroup(sg_patients, SYNTH_COUNTS)
    assert set(w.keys()) == set(sg_patients)
    assert np.all(w["synth_A"] == 1.0)
    assert np.all(w["synth_B"] == 1.0)
    assert len(w["synth_A"]) == 60
    assert len(w["synth_B"]) == 80


def test_build_driver_only_weights():
    """Driver-only weights must equal 1 for TP53/KRAS mutations and 0
    elsewhere."""
    drivers = ["TP53", "KRAS"]
    w = cps.build_driver_only_weights(SYNTH_COUNTS, drivers)
    assert np.all(w["synth_A"] == 1.0)  # all TP53
    assert np.all(w["synth_B"] == 1.0)  # all KRAS
    assert np.all(w["synth_C"] == 0.0)  # all TTN


def test_build_topk_weights_for_subgroup_drops_outside_patients():
    """Top-K weights for a subgroup must exclude patients outside the
    subgroup from the returned dict."""
    matches = {"snv_matched": [], "indel_matched": []}
    # Give every mutation in synth_A and synth_B a PHRED of 25 (so the
    # top-K filter behaves predictably).
    for p in ("synth_A", "synth_B"):
        for i, m in enumerate(SYNTH_COUNTS[p]):
            matches["snv_matched"].append({
                "sample": p, "chrom": m["chrom"], "pos": m["pos"],
                "ref": m["ref"], "alt": m["alt"], "cadd_phred": 25.0,
            })

    sg_patients = ["synth_A"]  # synth_B is intentionally excluded
    w = cps.build_topk_weights_for_subgroup(sg_patients, SYNTH_COUNTS,
                                            matches, top_k=10)
    assert set(w.keys()) == {"synth_A"}  # synth_B dropped
    assert len(w["synth_A"]) == 60
    assert int(w["synth_A"].sum()) == 10


def test_run_subgroup_panel_detection_auc_safe_tf():
    """At TF=1% with a synthetic cohort, AUC must be high and sens@99%
    defined."""
    weights = cps.build_uniform_weights_for_subgroup(
        ["synth_A", "synth_B"], SYNTH_COUNTS
    )
    res = cps.run_subgroup_panel_detection(
        SYNTH_COUNTS, weights, ["synth_A", "synth_B"],
        tumor_fractions=[0.01], seeds=[42, 123, 456],
    )
    auc_row = next(r for r in res["panel_llr_weighted"]
                   if r["metric"] == "auc")
    assert auc_row["mean"] >= 0.95, f"Expected AUC >=0.95 at TF=1%, got {auc_row['mean']}"
    sens99_row = next(r for r in res["panel_llr_weighted"]
                      if r["metric"] == "sens_at_99_spec")
    # Sens@99% must be a finite float in [0, 1]
    assert 0.0 <= sens99_row["mean"] <= 1.0


def test_run_subgroup_panel_detection_drops_zero_weight_patients():
    """If a subgroup patient has all-zero weights (e.g., driver-only panel
    on a patient with no driver mutations), the patient should be silently
    dropped and recorded in dropped_patients."""
    # synth_C has only TTN; driver-only panel yields all-zero weights.
    weights = cps.build_driver_only_weights(SYNTH_COUNTS, ["TP53", "KRAS"])
    sg = ["synth_A", "synth_C"]  # synth_C gets dropped
    res = cps.run_subgroup_panel_detection(
        SYNTH_COUNTS, weights, sg, tumor_fractions=[0.01], seeds=[42]
    )
    auc_row = next(r for r in res["panel_llr_weighted"]
                   if r["metric"] == "auc")
    assert "synth_C" in auc_row["dropped_patients"]
    assert auc_row["n_patients"] == 1
