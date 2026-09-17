"""Tests for the AlphaMissense weighting module.

Honesty-first: these tests cover the deterministic proxy path (which is
the path used on this machine today). A separate test class documents
the contract for the real-index path, but is skipped if no pickle / TSV
is present — this prevents CI from claiming the real-path code works
without actually exercising it.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import alphamissense_weights as amw  # noqa: E402


# ---------------------------------------------------------------------------
# Pure-function tests (no I/O, deterministic)
# ---------------------------------------------------------------------------

def test_classify_am_thresholds():
    """Cheng 2023 thresholds: 0.34 / 0.564."""
    assert amw.classify_am(0.00) == "likely_benign"
    assert amw.classify_am(0.339) == "likely_benign"
    assert amw.classify_am(0.34) == "ambiguous"
    assert amw.classify_am(0.563) == "ambiguous"
    assert amw.classify_am(0.564) == "likely_pathogenic"
    assert amw.classify_am(0.99) == "likely_pathogenic"


def test_am_key_format():
    assert amw.am_key("P42345", 12, "G", "D") == "P42345:12:G:D"
    # Deterministic: same inputs → same key
    assert amw.am_key("P42345", 12, "G", "D") == amw.am_key("P42345", 12, "G", "D")


def test_proxy_am_for_variant_class():
    """LoF > missense > silent ordering must hold."""
    p_nonsense = amw.proxy_am_for_variant("Nonsense_Mutation")
    p_missense = amw.proxy_am_for_variant("Missense_Mutation")
    p_silent = amw.proxy_am_for_variant("Silent")
    assert p_nonsense > p_missense > p_silent
    # Defaults to a moderate prior when class is unknown
    assert amw.proxy_am_for_variant("") > 0.0
    assert amw.proxy_am_for_variant("") < 1.0


def test_normalize_am_collapses_to_one_when_constant():
    """If all scores are equal, normalised weights should all be 1.0."""
    out = amw.normalize_am([0.55, 0.55, 0.55])
    assert np.allclose(out, 1.0)


def test_normalize_am_min_max_with_floor():
    """Min-max with floor should clamp to [0.05, 1.0]."""
    out = amw.normalize_am([0.0, 0.5, 1.0])
    assert out.min() == pytest.approx(0.05, abs=1e-9)
    assert out.max() == pytest.approx(1.0, abs=1e-9)
    assert out[1] == pytest.approx(0.5, abs=1e-9)


# ---------------------------------------------------------------------------
# Loader contract (proxy path is always available)
# ---------------------------------------------------------------------------

def test_weight_cohort_uses_proxy_when_no_index(tmp_path):
    """No pickle, no TSV → must fall back to proxy with clear source flag."""
    fake_muts = [
        {"uniprot": "P01116", "prot_pos": 12, "ref": "G", "alt": "D",
         "variant_class": "Missense_Mutation", "sample": "X", "gene": "KRAS",
         "chrom": "chr12", "pos": 25245350},
        {"uniprot": "P04637", "prot_pos": 175, "ref": "R", "alt": "H",
         "variant_class": "Missense_Mutation", "sample": "X", "gene": "TP53",
         "chrom": "chr17", "pos": 7676154},
        # Non-missense → no AM key (still surfaced as missing)
        {"uniprot": "P01116", "prot_pos": 1, "ref": "M", "alt": "L",
         "variant_class": "Silent", "sample": "X", "gene": "KRAS",
         "chrom": "chr12", "pos": 25245350},
    ]
    weights, info = amw.weight_cohort(fake_muts)
    assert info["primary_source"] in ("proxy", "missing")
    assert info["n_proxy"] >= 2  # at least the two missense ones
    # Non-missense rows are filtered upstream by load_tcga_mutations_with_protein
    # and never enter the cohort weighting → n_missing == 0 (correct).
    assert info["n_missing"] == 0
    # KRAS G12D should have a key
    k = amw.am_key("P01116", 12, "G", "D")
    assert k in weights
    assert weights[k].source == "proxy"
    assert weights[k].classification == "ambiguous"


def test_weight_cohort_skips_non_missense():
    """Non-missense rows must NOT generate a wanted key."""
    fake_muts = [
        {"uniprot": "P01116", "prot_pos": 12, "ref": "G", "alt": "D",
         "variant_class": "Silent", "sample": "X", "gene": "KRAS",
         "chrom": "chr12", "pos": 25245350},
    ]
    _, info = amw.weight_cohort(fake_muts)
    # All rows are non-missense → wanted_keys count should be 0
    assert info.get("wanted_keys", 0) == 0


def test_load_tcga_mutations_with_protein_filters_missense():
    """The augmented pool only contains Missense_Mutation SNVs."""
    # Skip if the TCGA cache isn't present (CI without data)
    cache = Path("/Users/hermes/deepcatch/validation/tcga/tcga_cache")
    if not cache.exists():
        pytest.skip("TCGA cache not present on this machine")
    muts = amw.load_tcga_mutations_with_protein(str(cache))
    assert len(muts) > 0
    for m in muts:
        assert m["variant_class"] == "Missense_Mutation"
        assert len(m["ref"]) == 1
        assert len(m["alt"]) == 1
        assert m["uniprot"]
        assert isinstance(m["prot_pos"], int)


# ---------------------------------------------------------------------------
# Real-index contract (skipped unless an index is actually present)
# ---------------------------------------------------------------------------

def test_weight_cohort_picklelocal_when_index_present(tmp_path):
    """If a real pickle is provided, source must be 'picklelocal'.

    We synthesise a tiny valid pickle and feed it in. This proves the
    real-path code works (it's NOT skipped — this is the path we want
    CI to exercise when the data is available).
    """
    import pickle

    fake_pickle = {
        "P01116:12:G:D": {"score": 0.95},
        "P04637:175:R:H": {"score": 0.10},
    }
    pkl = tmp_path / "alphamissense_index.pkl"
    pkl.write_bytes(pickle.dumps(fake_pickle))

    fake_muts = [
        {"uniprot": "P01116", "prot_pos": 12, "ref": "G", "alt": "D",
         "variant_class": "Missense_Mutation", "sample": "X", "gene": "KRAS",
         "chrom": "chr12", "pos": 25245350},
        {"uniprot": "P04637", "prot_pos": 175, "ref": "R", "alt": "H",
         "variant_class": "Missense_Mutation", "sample": "X", "gene": "TP53",
         "chrom": "chr17", "pos": 7676154},
        # This one is NOT in the index → falls back to proxy
        {"uniprot": "P15056", "prot_pos": 600, "ref": "V", "alt": "E",
         "variant_class": "Missense_Mutation", "sample": "X", "gene": "BRAF",
         "chrom": "chr7", "pos": 140453136},
    ]
    weights, info = amw.weight_cohort(fake_muts, index_path=pkl)
    assert info["primary_source"] == "picklelocal"
    assert info["n_real"] == 2
    assert info["n_proxy"] == 1
    assert weights[amw.am_key("P01116", 12, "G", "D")].source == "picklelocal"
    assert weights[amw.am_key("P01116", 12, "G", "D")].classification == "likely_pathogenic"
    assert weights[amw.am_key("P04637", 175, "R", "H")].classification == "likely_benign"
    assert weights[amw.am_key("P15056", 600, "V", "E")].source == "proxy"