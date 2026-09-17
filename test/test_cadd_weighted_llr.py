"""Tests for the CADD-weighted panel LLR extension (DeepCatch v2.2).

Validates that:
  - cadd_weighted_llr module imports cleanly
  - build_per_patient_weights correctly looks up scores by (patient, chrom,
    pos, ref, alt)
  - build_topk_per_patient_weights filters to top-K and to PHRED >= min_phred
  - run_weighted_panel_detection reuses the real_tcga_validation simulator and
    preserves AUC at >= 0.5% (sanity check)
"""
import json
import sys
from pathlib import Path

import numpy as np

# Set up import path for the scripts directory
SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR.parent))  # repo root (for real_tcga_validation)
sys.path.insert(0, str(SCRIPTS_DIR))

import cadd_weighted_llr as cwl  # noqa: E402


def test_module_imports():
    """The module should expose the expected public functions."""
    assert hasattr(cwl, "load_or_match_caddings")
    assert hasattr(cwl, "build_per_patient_weights")
    assert hasattr(cwl, "build_topk_per_patient_weights")
    assert hasattr(cwl, "run_weighted_panel_detection")


def test_build_per_patient_weights_median_imputation():
    """Weights for matched mutations must equal CADD PHRED; unmatched must
    fall back to the median PHRED across matched."""
    matches = {
        "snv_matched": [
            {"sample": "p1", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "T",
             "cadd_phred": 25.0},
            {"sample": "p1", "chrom": "chr1", "pos": 200, "ref": "C", "alt": "G",
             "cadd_phred": 15.0},
            {"sample": "p1", "chrom": "chr1", "pos": 300, "ref": "G", "alt": "A",
             "cadd_phred": 35.0},
        ],
        "indel_matched": [],
    }
    cohort = {
        "p1": [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},  # matched -> 25
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "G"},  # matched -> 15
            {"chrom": "chr1", "pos": 250, "ref": "T", "alt": "C"},  # unmatched -> median 25
            {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "A"},  # matched -> 35
        ]
    }
    weights = cwl.build_per_patient_weights(cohort, matches, imputation="median")
    assert "p1" in weights
    w = weights["p1"]
    assert len(w) == 4
    assert w[0] == 25.0
    assert w[1] == 15.0
    assert w[2] == 25.0  # median of [25, 15, 35] = 25
    assert w[3] == 35.0


def test_build_per_patient_weights_zero_imputation():
    """Unmatched weights must be 0 when imputation='zero'."""
    matches = {
        "snv_matched": [
            {"sample": "p1", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "T",
             "cadd_phred": 25.0},
        ],
        "indel_matched": [],
    }
    cohort = {
        "p1": [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "G"},  # unmatched
        ]
    }
    weights = cwl.build_per_patient_weights(cohort, matches, imputation="zero")
    assert weights["p1"][0] == 25.0
    assert weights["p1"][1] == 0.0


def test_build_topk_per_patient_weights():
    """Top-K filter keeps top-K by PHRED and sets others to 0."""
    matches = {
        "snv_matched": [
            {"sample": "p1", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "T",
             "cadd_phred": 10.0},
            {"sample": "p1", "chrom": "chr1", "pos": 200, "ref": "C", "alt": "G",
             "cadd_phred": 30.0},
            {"sample": "p1", "chrom": "chr1", "pos": 300, "ref": "G", "alt": "A",
             "cadd_phred": 25.0},
            {"sample": "p1", "chrom": "chr1", "pos": 400, "ref": "T", "alt": "C",
             "cadd_phred": 5.0},
        ],
        "indel_matched": [],
    }
    cohort = {
        "p1": [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},  # PHRED 10
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "G"},  # PHRED 30 (top)
            {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "A"},  # PHRED 25 (2nd)
            {"chrom": "chr1", "pos": 400, "ref": "T", "alt": "C"},  # PHRED 5
        ]
    }
    # Top-K=2: keep positions 200 and 300
    weights = cwl.build_topk_per_patient_weights(cohort, matches, top_k=2)
    assert weights["p1"][0] == 0.0   # PHRED 10 not in top 2
    assert weights["p1"][1] == 1.0   # PHRED 30
    assert weights["p1"][2] == 1.0   # PHRED 25
    assert weights["p1"][3] == 0.0   # PHRED 5


def test_build_topk_min_phred():
    """Top-K with min_phred filter excludes mutations below threshold."""
    matches = {
        "snv_matched": [
            {"sample": "p1", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "T",
             "cadd_phred": 10.0},  # below threshold 20
            {"sample": "p1", "chrom": "chr1", "pos": 200, "ref": "C", "alt": "G",
             "cadd_phred": 30.0},  # above
            {"sample": "p1", "chrom": "chr1", "pos": 300, "ref": "G", "alt": "A",
             "cadd_phred": 22.0},  # above
        ],
        "indel_matched": [],
    }
    cohort = {
        "p1": [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "G"},
            {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "A"},
        ]
    }
    # Top-K=10, threshold PHRED >= 20: keep positions 200 and 300, drop 100
    weights = cwl.build_topk_per_patient_weights(cohort, matches, top_k=10, min_phred=20.0)
    assert weights["p1"][0] == 0.0
    assert weights["p1"][1] == 1.0
    assert weights["p1"][2] == 1.0


def test_weighted_panel_detection_baseline_uniform_auc_preserved():
    """Weighted-aggregation with median imputation should preserve baseline AUC
    at non-ultra-low TFs. We run a tiny synthetic cohort."""
    from real_tcga_validation import simulate_cfdna_from_real

    cohort = {
        "synth_patient_1": [
            {"gene": "TP53", "tumor_vaf": 0.3, "t_alt": 30, "t_depth": 100,
             "n_alt": 0, "n_ref": 100, "normal_error_rate": 0.001,
             "variant_class": "Missense_Mutation", "sample": "synth_patient_1",
             "chrom": "chr17", "pos": 7674220, "ref": "G", "alt": "A"}
            for _ in range(200)
        ],
        "synth_patient_2": [
            {"gene": "KRAS", "tumor_vaf": 0.4, "t_alt": 40, "t_depth": 100,
             "n_alt": 0, "n_ref": 100, "normal_error_rate": 0.001,
             "variant_class": "Missense_Mutation", "sample": "synth_patient_2",
             "chrom": "chr12", "pos": 25245350, "ref": "C", "alt": "A"}
            for _ in range(150)
        ],
    }
    # Fake matches for half the panel positions; impute median for the other half
    matches = {"snv_matched": [], "indel_matched": []}
    for i, mut in enumerate(cohort["synth_patient_1"]):
        if i % 2 == 0:
            matches["snv_matched"].append({**mut, "cadd_phred": 25.0})
    for i, mut in enumerate(cohort["synth_patient_2"]):
        if i % 2 == 0:
            matches["snv_matched"].append({**mut, "cadd_phred": 18.0})
    weights = cwl.build_per_patient_weights(cohort, matches, imputation="median")

    # Run only at TF=1% (a "safe" TF where AUC should be ~1)
    res = cwl.run_weighted_panel_detection(
        cohort, weights, tumor_fractions=[0.01], seeds=[42, 123, 456]
    )
    auc_at_tf1 = res["panel_llr_weighted"][0]  # metric='auc' record
    assert auc_at_tf1["metric"] == "auc"
    assert auc_at_tf1["mean"] >= 0.95, f"Expected >=0.95 at TF=1%, got {auc_at_tf1['mean']}"
    # Uniform should also be high
    auc_uniform_at_tf1 = res["panel_llr_uniform"][0]
    assert auc_uniform_at_tf1["mean"] >= 0.95


def test_weighted_panel_detection_reproducibility():
    """A second call with the same seed should produce the same AUC."""
    cohort = {
        "synth_patient_3": [
            {"gene": "TP53", "tumor_vaf": 0.3, "t_alt": 30, "t_depth": 100,
             "n_alt": 0, "n_ref": 100, "normal_error_rate": 0.001,
             "variant_class": "Missense_Mutation", "sample": "synth_patient_3",
             "chrom": "chr17", "pos": 7674220, "ref": "G", "alt": "A"}
            for _ in range(100)
        ]
    }
    matches = {"snv_matched": [
        {**cohort["synth_patient_3"][0], "cadd_phred": 25.0}
    ], "indel_matched": []}
    weights = cwl.build_per_patient_weights(cohort, matches, imputation="median")
    res1 = cwl.run_weighted_panel_detection(
        cohort, weights, tumor_fractions=[0.05], seeds=[42]
    )
    res2 = cwl.run_weighted_panel_detection(
        cohort, weights, tumor_fractions=[0.05], seeds=[42]
    )
    auc1 = res1["panel_llr_weighted"][0]["per_seed"]["42"]
    auc2 = res2["panel_llr_weighted"][0]["per_seed"]["42"]
    assert abs(auc1 - auc2) < 1e-6, f"Non-reproducible: {auc1} vs {auc2}"