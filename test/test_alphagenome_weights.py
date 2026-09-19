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

    # Method coverage (5 -> 8 after the 3-fix patch: 2 baselines + 3 top-K by
    # w*LLR + 3 top-K by raw AVI)
    expected_methods = {
        "panel_llr_uniform", "panel_llr_avi",
        "panel_llr_topk_500", "panel_llr_topk_1000", "panel_llr_topk_2000",
        "panel_llr_topk_by_avi_500", "panel_llr_topk_by_avi_1000",
        "panel_llr_topk_by_avi_2000",
    }
    assert expected_methods.issubset(set(r["results"].keys())), (
        f"missing methods: {expected_methods - set(r['results'].keys())}"
    )

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


# ---------------------------------------------------------------------------
# Regression tests for the 3-fix patch (issue #1: indels; #2: top-K-by-AVI)
# ---------------------------------------------------------------------------

def test_load_tcga_mutations_keeps_indels():
    """After patch 1, indels must be retained (was 0, now 992 of 124,841).

    This is the headline honesty regression: previously the loader dropped
    Frame_Shift_Del/Ins, In_Frame_Del/Ins, Splice_Site variants with
    `len(ref) != 1`, costing 992 mutations (~0.79% of cohort).
    """
    cache_dir = Path(__file__).resolve().parent.parent / "validation" / "tcga" / "tcga_cache"
    if not cache_dir.exists():
        pytest.skip(f"{cache_dir} not present — TCGA cache not downloaded")
    muts = aw.load_tcga_mutations_with_ref_alt(str(cache_dir))
    assert len(muts) > 0
    indels = [m for m in muts if len(m["ref"]) != 1 or len(m["alt"]) != 1]
    assert len(indels) > 0, (
        "loader dropped indels again — patch 1 regressed; expected >= 992 "
        "indels in LUAD MAF cache"
    )
    # Sanity: indel fraction is small (~0.79%) but non-zero
    assert len(indels) < len(muts) * 0.05, (
        f"indel fraction {len(indels)/len(muts):.2%} suspiciously large"
    )
    # Indels should all carry a non-empty variant_class so the proxy has a target
    assert all(m.get("variant_class") for m in indels), (
        "indels without variant_class cannot receive proxy weights"
    )


def test_indel_keys_have_source_proxy_when_weighted():
    """After patch 1, indels weighted via the proxy path must be tagged
    ``source="proxy"`` (not silently missing). The proxy is honest because
    it is clearly flagged; silently dropping indels is not.
    """
    fake_muts = [
        # SNV
        {"chrom": "1", "pos": 100, "ref": "A", "alt": "T",
         "variant_class": "Missense_Mutation"},
        # Indels (the four classes the proxy table covers)
        {"chrom": "2", "pos": 200, "ref": "AT", "alt": "A",
         "variant_class": "Frame_Shift_Del"},
        {"chrom": "3", "pos": 300, "ref": "A", "alt": "ATG",
         "variant_class": "Frame_Shift_Ins"},
        {"chrom": "4", "pos": 400, "ref": "AT", "alt": "A",
         "variant_class": "In_Frame_Del"},
        {"chrom": "5", "pos": 500, "ref": "A", "alt": "ATG",
         "variant_class": "In_Frame_Ins"},
    ]
    weights, primary_source = aw.weight_cohort(fake_muts)
    assert primary_source == "proxy"
    for vk, w in weights.items():
        # variant_key format: "chrN:POS:REF>ALT"  -> 4 colon-separated fields
        # where REF>ALT may itself contain colons if REF or ALT is multi-base.
        # We split on ":" max 3 times to separate (chrom, pos, ref>alt).
        parts = vk.split(":", 3)       # ['chrN', 'POS', 'REF>ALT']
        if len(parts) != 4:
            # Indel keys may have extra colons in the allele (e.g. insertion).
            # Re-join the trailing portion as the allele field.
            allele_field = ":".join(parts[2:])
        else:
            allele_field = parts[3]
        ref, alt = allele_field.split(">", 1)
        is_indel = len(ref) != 1 or len(alt) != 1
        if is_indel:
            assert w.source == "proxy", (
                f"indel {vk} should be weighted by proxy; got source={w.source!r}"
            )
            assert w.avi_score > 0.0, (
                f"indel {vk} proxy weight should be positive; got {w.avi_score}"
            )


def test_panel_score_with_weights_by_avi_is_deterministic():
    """The new top-K-by-AVI selector must be deterministic.

    Same inputs -> same output, every time. (We compare two separate calls;
    np.argsort on identical input is guaranteed stable for our numpy version.)
    """
    # Import the runner's helper without running the whole module CLI.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from alphagenome_panel_run import _panel_score_with_weights_by_avi

    rng = np.random.default_rng(0)
    llr = rng.normal(size=200)
    avi = rng.uniform(0.05, 1.0, size=200)
    weights = avi.copy()      # weights = AVI norm in the runner
    out1 = _panel_score_with_weights_by_avi(llr, avi, weights, top_k=20)
    out2 = _panel_score_with_weights_by_avi(llr, avi, weights, top_k=20)
    assert out1 == out2, f"non-deterministic: {out1} vs {out2}"
    # Length-mismatch guard
    with pytest.raises(ValueError):
        _panel_score_with_weights_by_avi(llr, avi[:50], weights, top_k=20)


def test_run_comparison_includes_topk_by_avi_methods():
    """The methods list in run_comparison must include the new top-K-by-AVI
    family alongside the existing top-K-by-w*LLR family.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import alphagenome_panel_run as apr

    # Inspect the source to confirm the method-family strings are emitted.
    # (We don't call run_comparison directly — it requires a full cohort.)
    src = (Path(__file__).resolve().parent.parent / "alphagenome_panel_run.py").read_text()
    assert "panel_llr_topk_by_avi_" in src, (
        "panel_llr_topk_by_avi_ not present in alphagenome_panel_run.py"
    )
    # Both the joint-criterion family and the new pure-prior family must exist
    assert "panel_llr_topk_" in src
    # And the helper must exist
    assert hasattr(apr, "_panel_score_with_weights_by_avi")
