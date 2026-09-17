"""Tests for the AlphaGenome Atlas weighting module and proxy-mode runner.

Honesty-first: these tests cover the deterministic proxy path (which is
the path used on this machine today, since no AlphaGenome API key is set
and no Tabix cache exists). The proxy path produces a clearly-flagged
``avi_primary_source: "proxy"`` JSON artefact; the tests below check
that the proxy flag is present and that the proxy itself is a stable
per-variant-class lookup (NOT a fabricated AVI score).

A separate skipped-by-default class documents the contract for the
real Atlas API path so CI cannot falsely claim that path works.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import alphagenome_weights as aw  # noqa: E402


# ---------------------------------------------------------------------------
# Pure-function tests (no I/O, deterministic)
# ---------------------------------------------------------------------------

def test_variant_key_is_canonical():
    """Chrom 'chr1' and '1' should produce the same key."""
    assert aw.variant_key("1", 100, "A", "T") == aw.variant_key("chr1", 100, "A", "T")
    assert aw.variant_key("chrX", 1234, "C", "G") == "chrX:1234:C>G"


def test_proxy_avi_for_variant_class_ordering():
    """LoF > missense > silent ordering must hold for the proxy.

    This is what the proxy captures: a coarse biological ranking of
    variant classes. The numbers are NOT published AVI scores; they are
    fixed Phred-like anchors used purely to exercise the weighted
    aggregator end-to-end.
    """
    p_nonsense = aw.proxy_avi_for_variant("Nonsense_Mutation")
    p_missense = aw.proxy_avi_for_variant("Missense_Mutation")
    p_silent = aw.proxy_avi_for_variant("Silent")
    assert p_nonsense > p_missense > p_silent
    # Unknown class falls back to a moderate default (5.0)
    assert aw.proxy_avi_for_variant("Unknown_Class") == 5.0
    assert aw.proxy_avi_for_variant("") == 5.0


def test_normalize_avi_floor_and_constant_collapse():
    """Floor of 0.05 prevents zero-silencing; constant input -> all 1.0."""
    out = aw.normalize_avi([0.55, 0.55, 0.55])
    assert np.allclose(out, 1.0)
    out2 = aw.normalize_avi([0.0, 5.0, 10.0])
    assert out2.min() >= 0.05  # floor enforced
    assert out2.max() == 1.0


def test_results_json_schema_and_proxy_flag():
    """The committed proxy-run JSON must exist, parse, and be clearly proxy.

    This is the headline honesty test: the artefact must NOT look like
    real AlphaGenome data. The proxy flag is the only honest way to ship
    these numbers.
    """
    results_path = Path(__file__).resolve().parent.parent / "results" / "alphagenome_weighted_llr.json"
    assert results_path.exists(), f"missing {results_path}"
    with open(results_path) as f:
        r = json.load(f)

    # Top-level shape
    assert "metadata" in r and "results" in r
    md = r["metadata"]
    for k in ("runner", "date", "n_patients", "seeds_used",
              "tumor_fractions", "methods_compared",
              "avi_primary_source", "avi_n_real", "avi_n_proxy"):
        assert k in md, f"metadata missing key: {k}"

    # Proxy flag is THE honesty gate
    assert md["avi_primary_source"] == "proxy", (
        f"expected proxy, got {md['avi_primary_source']!r} — do not commit "
        "non-proxy numbers as if they came from this proxy runner"
    )
    assert md["avi_n_real"] == 0, (
        "this run was supposed to be proxy-only; avi_n_real must be 0"
    )
    assert md["avi_n_proxy"] > 0

    # Method coverage
    expected_methods = {
        "panel_llr_uniform", "panel_llr_avi",
        "panel_llr_topk_500", "panel_llr_topk_1000", "panel_llr_topk_2000",
    }
    assert expected_methods.issubset(set(r["results"].keys()))

    # AUC regression floor at 0.1% ctDNA (proxy run must not regress)
    for m in expected_methods:
        row = next((x for x in r["results"][m] if x["tumor_fraction"] == 0.001), None)
        assert row is not None, f"{m} missing 0.001 row"
        assert row["mean"] >= 0.92, (
            f"{m} AUC @ 0.1% = {row['mean']:.4f} regressed below the 0.92 floor"
        )


def test_tos_compliance_comments_preserved():
    """The ToS-compliance docstring must remain in alphagenome_weights.py.

    The whole reason this proxy exists is to honour the AlphaGenome
    Terms of Use (no training of other models with AVI scores, etc).
    Stripping or weakening these comments is a silent regression.
    """
    src_path = Path(__file__).resolve().parent.parent / "alphagenome_weights.py"
    src = src_path.read_text()
    # Both anchors must be present
    assert "TERMS-OF-SERVICE NOTES" in src
    assert "should not be used for the" in src
    assert "non-training clause" in src or "non-commercial" in src
    # The proxy is explicitly labelled as such
    assert "proxy" in src.lower()
    assert "PROXY_USED" in src


def test_weight_cohort_falls_back_to_proxy_without_key():
    """With no API key and no Tabix, weight_cohort must use the proxy path.

    This is the behaviour the runner relies on when ALPHAGENOME_API_KEY
    is unset. If this ever silently switches to atlas_api or tabix_local
    the JSON will lie.
    """
    fake_muts = [
        {"chrom": "1", "pos": 100, "ref": "A", "alt": "T", "variant_class": "Missense_Mutation"},
        {"chrom": "2", "pos": 200, "ref": "C", "alt": "G", "variant_class": "Nonsense_Mutation"},
        {"chrom": "X", "pos": 300, "ref": "G", "alt": "A", "variant_class": "Silent"},
    ]
    # No api_key, no tabix_path -> proxy fallback
    weights, primary_source = aw.weight_cohort(fake_muts)
    assert primary_source == "proxy"
    for vk, w in weights.items():
        assert w.source == "proxy"
        assert w.avi_score > 0.0  # the proxy table assigns positive values
    # Ordering sanity: nonsense > missense > silent within these three
    vk_mis = aw.variant_key("1", 100, "A", "T")
    vk_non = aw.variant_key("2", 200, "C", "G")
    vk_sil = aw.variant_key("X", 300, "G", "A")
    assert weights[vk_non].avi_score > weights[vk_mis].avi_score > weights[vk_sil].avi_score


# ---------------------------------------------------------------------------
# Real Atlas API path — documented but skipped (no API key in CI)
# ---------------------------------------------------------------------------

@pytest.mark.skip(reason="No ALPHAGENOME_API_KEY; real Atlas path untested in CI")
def test_fetch_avi_via_atlas_requires_key():
    """If a key is ever added, fetch_avi_via_atlas must succeed."""
    pytest.skip("requires real AlphaGenome API key")
