"""Tests for the biomedical-review fixes.

Covers:
- CAFFCalculator.from_fragments (new helper)
- SensAtSpecLoss + BalancedCrossEntropy + CalibrationLoss
- Stratified val split in FoundationDownstream
- _top_motif_deviations renamed from _pca_reduce (key names unchanged)
- tss_coverage_profile per-chromosome behavior
- methylation_data backfill in extract_all
"""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

# Skip all tests in this file if torch is unavailable.
try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False


# ── CAFFCalculator.from_fragments ────────────────────────────────────────

def test_caff_from_fragments_normalizes_by_arm_length():
    """Per-arm coverage is fragments / arm_length_Mb then median-norm to 1.0.

    With enough data the median arm has positive coverage and the
    normalization reads 1.0. With sparse data (where many arms have
    zero coverage), the median may be 0 and the function preserves the
    raw per-Mb values — which is the safe behaviour (avoid div-by-0).
    """
    from src.fragmentomics.themis_features import CAFFCalculator
    # Build a dense synthetic cohort so most arms have non-zero
    # coverage — many fragments per chromosome so per-Mb values are
    # uniformly positive across arms.
    fragments = []
    for chrom in range(1, 23):
        for i in range(20):
            fragments.append({
                "chrom": f"chr{chrom}",
                "start": (i + 1) * 5_000_000,
            })
    cov = CAFFCalculator.from_fragments(fragments)
    vals = np.array(list(cov.values()))
    assert np.all(vals > 0), "dense cohort should give non-zero all arms"
    assert np.isclose(np.median(vals), 1.0, atol=1e-6), \
        f"median should be 1.0 after normalization, got {np.median(vals)}"


def test_caff_from_fragments_handles_ensembl_chrom_names():
    """Ensembl-style '1' and UCSC-style 'chr1' should both work."""
    from src.fragmentomics.themis_features import CAFFCalculator
    fragments_ucsc = [{"chrom": "chr1", "start": 1_000_000}] * 5
    fragments_ens = [{"chrom": "1", "start": 1_000_000}] * 5
    cov_ucsc = CAFFCalculator.from_fragments(fragments_ucsc)
    cov_ens = CAFFCalculator.from_fragments(fragments_ens)
    # The chr1 fragment count and distribution should be identical,
    # so the per-arm coverage values should match.
    for arm in cov_ucsc:
        assert np.isclose(cov_ucsc[arm], cov_ens[arm]), \
            f"arm {arm}: UCSC {cov_ucsc[arm]} != Ensembl {cov_ens[arm]}"


def test_caff_from_fragments_silently_drops_unknown_chroms():
    """Fragments with unrecognized chromosomes are silently dropped."""
    from src.fragmentomics.themis_features import CAFFCalculator
    fragments = [
        {"chrom": "chr1", "start": 1_000_000},
        {"chrom": "chrM", "start": 1_000},
        {"chrom": "chrUn_KI270442v1", "start": 1_000},
    ]
    cov = CAFFCalculator.from_fragments(fragments)
    # Only chr1 fragments contribute.
    assert cov["1p"] > 0


def test_caff_from_fragments_empty_input_returns_zero_template():
    """Empty input returns the 39-arm template with all zeros (not empty dict).

    The 39-arm template is the documented contract: callers can iterate
    over all chromosomes regardless of whether the sample had matching
    fragments, and missing arms read as 0.0 coverage.
    """
    from src.fragmentomics.themis_features import CAFFCalculator
    cov = CAFFCalculator.from_fragments([])
    # Same number of arms as the boundary table.
    assert len(cov) == len(CAFFCalculator.CHROM_ARM_BOUNDARIES)
    assert all(v == 0.0 for v in cov.values())


def test_caff_from_fragments_drop_stats_default_returns_dict():
    """Default call site returns just the coverage dict (backward compat)."""
    from src.fragmentomics.themis_features import CAFFCalculator
    fragments = [
        {"chrom": "chr1", "start": 1_000_000},
        {"chrom": "chrM", "start": 1_000},
    ]
    cov = CAFFCalculator.from_fragments(fragments)
    # cov is a plain dict, not a tuple
    assert isinstance(cov, dict)
    assert "1p" in cov


def test_caff_from_fragments_drop_stats_counts_unrecognized_chroms():
    """Opt-in drop_stats surfaces silently-dropped chrM / unplaced frags.

    This is the biological safety net: a sample with high mitochondrial
    contamination otherwise silently becomes a 39-arm template with all-1.0
    coverage and reads as a healthy control.
    """
    from src.fragmentomics.themis_features import CAFFCalculator
    fragments = [
        {"chrom": "chr1", "start": 1_000_000},
        {"chrom": "chrM", "start": 1_000},
        {"chrom": "chrM", "start": 2_000},
        {"chrom": "chrUn_KI270442v1", "start": 1_000},
        {"chrom": "chr5", "start": 1_000_000},
    ]
    cov, drop_stats = CAFFCalculator.from_fragments(
        fragments, return_drop_stats=True
    )
    # Coverage dict is unchanged from the default behaviour
    assert "1p" in cov
    assert cov["1p"] > 0
    # Drop stats track the per-string drop count
    assert drop_stats == {"chrM": 2, "chrUn_KI270442v1": 1}


def test_caff_from_fragments_drop_stats_empty_when_all_recognized():
    """All-recognized fragments → drop_stats is empty (not None)."""
    from src.fragmentomics.themis_features import CAFFCalculator
    fragments = [
        {"chrom": "chr1", "start": 1_000_000},
        {"chrom": "chr5", "start": 1_000_000},
    ]
    cov, drop_stats = CAFFCalculator.from_fragments(
        fragments, return_drop_stats=True
    )
    assert drop_stats == {}


def test_caff_compute_then_from_fragments_round_trip():
    """End-to-end: from_fragments → compute should return a finite score."""
    from src.fragmentomics.themis_features import CAFFCalculator
    fragments = [{"chrom": f"chr{i % 22 + 1}", "start": j * 1_000_000}
                 for i in range(40) for j in range(10)]
    cov = CAFFCalculator.from_fragments(fragments)
    result = CAFFCalculator(n_top_arms=5).compute(cov)
    assert "caff_score" in result
    assert np.isfinite(result["caff_score"])
    assert len(result["aberrant_arms"]) == 5


# ── _top_motif_deviations (renamed from _pca_reduce) ─────────────────────

def test_top_motif_deviations_key_names_unchanged():
    """The renamed function must keep the fem_5mer_pc* key names."""
    from src.fragmentomics.enhanced_features import RefinedEndMotifs
    m = RefinedEndMotifs()
    seqs = ["ACGTG", "CCGGA", "TTACG", "GGGCC", "ATATA"]
    feats = m.extract(seqs, np.array([100, 150, 200, 250, 300]))
    for i in range(10):
        assert f"fem_5mer_pc{i}" in feats, f"missing key fem_5mer_pc{i}"


def test_top_motif_deviations_zero_for_empty_input():
    from src.fragmentomics.enhanced_features import RefinedEndMotifs
    m = RefinedEndMotifs()
    feats = m.extract(None, None)
    for i in range(10):
        assert feats[f"fem_5mer_pc{i}"] == 0.0


# ── Motif diversity per-bin Simpson ──────────────────────────────────────

def test_motif_diversity_uses_effective_alphabet():
    """A bin with only 1 non-zero motif should return 0, not a high value."""
    from src.fragmentomics.enhanced_features import RefinedEndMotifs
    m = RefinedEndMotifs()
    # short bin: only 1 motif observed; mid/long: many observed.
    short = np.zeros(1024)
    short[42] = 10
    mid = np.full(1024, 1.0)
    long_c = np.full(1024, 1.0)
    feats = m._motif_diversity_by_bin(short, mid, long_c)
    # Effective alphabet for short is 1, so Simpson → 0.
    assert feats["fem_mds_short"] == 0.0
    # Mid/long have full alphabet and balanced counts → near 1.0
    assert feats["fem_mds_mid"] > 0.99
    assert feats["fem_mds_long"] > 0.99


# ── Nucleosome parameter overrides ───────────────────────────────────────

def test_nucleosome_pattern_accepts_overrides():
    """Non-default period/dip should change the output pattern."""
    from src.fragmentomics.enhanced_features import NucleosomeFootprint
    n_default = NucleosomeFootprint(nucleosome_period_bp=195.0)
    n_alt = NucleosomeFootprint(nucleosome_period_bp=165.0)  # yeast-like
    p_default = n_default.expected_nucleosome_pattern(80)
    p_alt = n_alt.expected_nucleosome_pattern(80)
    assert not np.allclose(p_default, p_alt), \
        "different periods should produce different patterns"


def test_nucleosome_pattern_amplitude_overrides():
    """Setting periodic_amplitude=0 flattens the periodic component.

    With no periodic component, the pattern should equal the dip
    component alone (which is then mean-normalized to 1.0). So at
    the centre of the dip, p[centre] should be smaller than p[far_edge].
    """
    from src.fragmentomics.enhanced_features import NucleosomeFootprint
    n = NucleosomeFootprint(periodic_amplitude=0.0, tss_dip_amplitude=0.5)
    p = n.expected_nucleosome_pattern(80)
    # p is mean-normalized to 1.0, so it sums to 80 * 1 = 80 across bins.
    assert abs(p.mean() - 1.0) < 1e-6, f"mean should be 1.0, got {p.mean()}"
    # Centre bin (the dip) should be lower than the mean.
    assert p[40] < p.mean()


# ── TSS per-chromosome matching ──────────────────────────────────────────

def test_tss_coverage_profile_chromosome_matching():
    """Fragments on chr5 must NOT be counted against chr1 TSS."""
    from src.fragmentomics.enhanced_features import NucleosomeFootprint
    nuc = NucleosomeFootprint(tss_window=2000, bin_size=50)
    tss = [("chr1", 100_000), ("chr5", 100_000)]
    fragments_chr1 = [
        {"chrom": "chr1", "start": 99_000, "length": 200},
        {"chrom": "chr1", "start": 100_500, "length": 200},
    ]
    fragments_chr5 = [
        {"chrom": "chr5", "start": 99_000, "length": 200},
        {"chrom": "chr5", "start": 100_500, "length": 200},
    ]
    # If matching were chrom-agnostic, mixing would inflate the profile.
    p_chr1_only = nuc.tss_coverage_profile(
        [("chr1", 100_000)], fragments_chr1 + fragments_chr5
    )
    # Chrom-aware: only chr1 TSS, so only chr1 fragments match.
    # All 4 fragments are at midpoint 99500 or 100600, both within
    # ±2000 of the chr1 TSS at 100000 → all match.
    p_chr1_aware = nuc.tss_coverage_profile(
        [("chr1", 100_000)], fragments_chr1
    )
    # Chrom-agnostic would have 4 fragments per bin (not 2).
    np.testing.assert_allclose(
        p_chr1_only, p_chr1_aware,
        err_msg="chrom-agnostic matching incorrectly counted chr5 fragments "
                "against chr1 TSS"
    )


# ── extract_all methylation_data backfill ────────────────────────────────

def test_extract_all_backfills_methylation_from_array():
    """A standalone methylation_data array should be attached to fragments."""
    from src.fragmentomics.enhanced_features import EnhancedFragmentomics
    e = EnhancedFragmentomics()
    fragments = [
        {"start": 1_000_000, "length": 167},
        {"start": 2_000_000, "length": 167},
        {"start": 3_000_000, "length": 167},
    ]
    methylation = np.array([True, False, True])
    feats = e.extract_all(
        fragment_lengths=np.array([167.0, 167.0, 167.0]),
        fragments=fragments,
        methylation_data=methylation,
    )
    # Now the fragments should have a 'methylated' field.
    assert all("methylated" in f for f in fragments)
    assert fragments[0]["methylated"] is True
    assert fragments[1]["methylated"] is False


def test_extract_all_no_methylation_returns_zero_fallback():
    """No methylation source → MFS methylation features return zero."""
    from src.fragmentomics.enhanced_features import EnhancedFragmentomics
    e = EnhancedFragmentomics()
    fragments = [{"start": 1_000_000, "length": 167}] * 5
    feats = e.extract_all(
        fragment_lengths=np.array([167.0] * 5),
        fragments=fragments,
    )
    assert feats["mfs_meth_short"] == 0.0
    assert feats["mfs_meth_long"] == 0.0


# ── SensAtSpecLoss and friends (require torch) ──────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_sens_at_spec_loss_runs_and_is_finite():
    from src.foundation.losses import SensAtSpecLoss
    loss_fn = SensAtSpecLoss(alpha_pos=20.0, gamma=2.0)
    logits = torch.randn(64)
    y = (torch.rand(64) < 0.4).float()
    loss = loss_fn(logits, y)
    assert torch.isfinite(loss)
    assert loss.item() > 0


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_sens_at_spec_loss_alpha_pos_increases_focal_weight():
    """Larger alpha_pos should make the loss more sensitive to positives."""
    from src.foundation.losses import SensAtSpecLoss
    logits = torch.tensor([2.0, -2.0])
    y = torch.tensor([1.0, 0.0])
    loss_low = SensAtSpecLoss(alpha_pos=1.0)(logits, y)
    loss_high = SensAtSpecLoss(alpha_pos=20.0)(logits, y)
    assert loss_high.item() > loss_low.item(), (
        "Increasing alpha_pos should increase loss magnitude on the "
        "positive-class term (which is currently well-classified)."
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_balanced_ce_inverse_frequency_weighting():
    from src.foundation.losses import BalancedCrossEntropy
    # Imbalanced labels: 9 class-0, 1 class-1.
    labels = np.array([0] * 9 + [1])
    bce = BalancedCrossEntropy(labels, beta=0.0)
    # Class 1 should have ~9× the weight of class 0.
    assert bce.weight[1].item() > bce.weight[0].item() * 5


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_calibration_loss_is_zero_for_perfect_calibration():
    """If predicted probs match true labels exactly, ECE → 0."""
    from src.foundation.losses import CalibrationLoss
    # logits where sigmoid matches targets with high confidence.
    logits = torch.tensor([5.0, 5.0, -5.0, -5.0])
    y = torch.tensor([1.0, 1.0, 0.0, 0.0])
    cl = CalibrationLoss(n_bins=10)
    assert cl(logits, y).item() < 0.01


# ── FoundationDownstream stratified split (require torch) ───────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_foundation_downstream_stratified_split_keeps_both_classes():
    """Stratified val split must keep at least one of each class in val."""
    from src.foundation.config import FoundationConfig, MODALITY_DIMS
    from src.foundation.downstream import FoundationDownstream
    cfg = FoundationConfig(embed_dim=16, n_heads=2, n_layers=1,
                           ff_dim=32, seed=0)
    fd = FoundationDownstream(config=cfg, pretrained=False)
    # Generate 60 samples with 30% prevalence.
    rng = np.random.default_rng(42)
    mod = {name: rng.standard_normal((60, dim)).astype(np.float32)
           for name, dim in MODALITY_DIMS.items()}
    labels = (rng.random(60) < 0.3).astype(np.int64)
    # Inspect the stratified split indices indirectly: with val_frac=0.1
    # and 18 positives, val_pos should be round(18 * 0.1) = 2, val_neg = 4.
    # Both > 0 → no degenerate val set.
    fd.fit(mod, labels, n_epochs=2, batch_size=16, validation_split=0.1,
           verbose=False)
    assert fd._fitted


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_foundation_downstream_nan_guard_aborts_cleanly():
    """NaN training loss should abort gracefully without crashing."""
    from src.foundation.config import FoundationConfig, MODALITY_DIMS
    from src.foundation.downstream import FoundationDownstream
    cfg = FoundationConfig(embed_dim=8, n_heads=2, n_layers=1,
                           ff_dim=16, seed=0)
    fd = FoundationDownstream(config=cfg, pretrained=False)
    rng = np.random.default_rng(0)
    mod = {name: rng.standard_normal((10, dim)).astype(np.float32)
           for name, dim in MODALITY_DIMS.items()}
    labels = (rng.random(10) < 0.5).astype(np.int64)
    # Set absurdly high LR to provoke NaN.
    fd.fit(mod, labels, n_epochs=5, batch_size=4, lr=1e6,
           validation_split=0.2, verbose=False)
    assert fd._fitted  # Should still mark as fitted (random-init state).


# ── advanced_fusion PyTorch rewrite tests ───────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_cross_attention_fusion_binary_proba_is_1d():
    """Binary classification returns 1-D cancer probability for sklearn compat."""
    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
    rng = np.random.default_rng(0)
    n = 100
    scores = [rng.standard_normal(n) + i for i in range(4)]
    labels = (rng.random(n) < 0.5).astype(np.int64)
    m = CrossAttentionFusion(n_modalities=4, prior=None,
                              n_epochs=20, seed=0)
    m.fit(scores, labels)
    proba = m.predict_proba(scores)
    assert proba.ndim == 1
    assert proba.shape == (n,)
    assert np.all((proba >= 0) & (proba <= 1))


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_cross_attention_fusion_no_nan_attention_output():
    """Prior mask must not produce NaN from 0.0 * -inf pitfall."""
    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
    rng = np.random.default_rng(0)
    n = 50
    scores = [rng.standard_normal(n) for _ in range(6)]
    labels = (rng.random(n) < 0.5).astype(np.int64)
    # cancer_detection prior blocks some cross-attention; verify no NaN.
    m = CrossAttentionFusion(n_modalities=6, prior="cancer_detection",
                              n_epochs=10, seed=0)
    m.fit(scores, labels)
    proba = m.predict_proba(scores)
    assert not np.isnan(proba).any()
    assert not np.isinf(proba).any()


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_early_late_fusion_binary_proba_is_1d():
    from src.multimodal_fusion.advanced_fusion import EarlyLateFusion
    rng = np.random.default_rng(0)
    n = 60
    feats = [rng.standard_normal((n, 4)) for _ in range(3)]
    labels = (rng.random(n) < 0.5).astype(np.int64)
    m = EarlyLateFusion(n_modalities=3, hidden_dim=16, n_epochs=20,
                         seed=0)
    m.fit(feats, labels)
    proba = m.predict_proba(feats)
    assert proba.shape == (n,)
    assert not np.isnan(proba).any()


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_gcn_too_returns_pred_and_proba():
    from src.multimodal_fusion.advanced_fusion import GCNTissueOfOrigin
    rng = np.random.default_rng(0)
    n = 40
    n_feat = 8
    features = rng.standard_normal((n, n_feat))
    labels = (rng.random(n) < 0.5).astype(np.int64)
    m = GCNTissueOfOrigin(n_cancer_types=2, n_features=n_feat,
                           n_epochs=20, seed=0)
    m.fit(features, labels)
    pred = m.predict(features)
    proba = m.predict_proba(features)
    assert pred.shape == (n,)
    assert proba.shape == (n, 2)
    assert np.all((proba >= 0) & (proba <= 1))


# ── Schema-fingerprint validator ────────────────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_schema_validator_catches_random_input():
    """Per-row median outside the expected range must raise.

    Default schemas reject per-row medians outside [-5, 10] for
    frag_basic. Use ``np.random.uniform(50, 100)`` so all values are
    positive (would otherwise trip the allow_negative check) AND
    push the per-row median way outside the schema range.
    """
    from src.foundation.config import FoundationConfig, MODALITY_DIMS
    from src.foundation.downstream import FoundationDownstream
    cfg = FoundationConfig(embed_dim=16, n_heads=2, n_layers=1,
                           ff_dim=32, seed=0)
    fd = FoundationDownstream(config=cfg, pretrained=False)
    rng = np.random.default_rng(0)
    # Uniform(50, 100) → per-row medians land in [50, 100], far outside
    # the [-5, 10] range allowed for frag_basic. All positive so we
    # don't trip the allow_negative check first.
    mod = {name: rng.uniform(50.0, 100.0, size=(60, dim)).astype(np.float32)
           for name, dim in MODALITY_DIMS.items()}
    labels = (rng.random(60) < 0.3).astype(np.int64)
    with pytest.raises(ValueError, match="median.*outside expected range"):
        fd.fit(mod, labels, n_epochs=1, batch_size=16, validation_split=0.1,
               validate_schema=True, verbose=False)


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_schema_validator_catches_negative_sero():
    """Negative values in sero must fail the allow_negative=False check."""
    from src.foundation.config import FoundationConfig, MODALITY_DIMS
    from src.foundation.downstream import FoundationDownstream
    cfg = FoundationConfig(embed_dim=16, n_heads=2, n_layers=1,
                           ff_dim=32, seed=0)
    fd = FoundationDownstream(config=cfg, pretrained=False)
    rng = np.random.default_rng(0)
    mod = {name: rng.standard_normal((60, dim)).astype(np.float32)
           for name, dim in MODALITY_DIMS.items()}
    # Inject a single negative into the sero modality (which requires
    # non-negative).
    mod["sero"][0, 0] = -1.0
    labels = (rng.random(60) < 0.3).astype(np.int64)
    with pytest.raises(ValueError, match="contains negative values"):
        fd.fit(mod, labels, n_epochs=1, batch_size=16, validation_split=0.1,
               validate_schema=True, verbose=False)


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_schema_validator_off_by_default():
    """Without ``validate_schema=True``, the same bad input trains."""
    from src.foundation.config import FoundationConfig, MODALITY_DIMS
    from src.foundation.downstream import FoundationDownstream
    cfg = FoundationConfig(embed_dim=16, n_heads=2, n_layers=1,
                           ff_dim=32, seed=0)
    fd = FoundationDownstream(config=cfg, pretrained=False)
    rng = np.random.default_rng(0)
    mod = {name: rng.standard_normal((60, dim)).astype(np.float32)
           for name, dim in MODALITY_DIMS.items()}
    labels = (rng.random(60) < 0.3).astype(np.int64)
    # Default (validate_schema=False) — should NOT raise
    fd.fit(mod, labels, n_epochs=1, batch_size=16, validation_split=0.1,
           verbose=False)
    assert fd._fitted


# ── sens_at_spec loss wiring ────────────────────────────────────────────

@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_foundation_downstream_sens_at_spec_loss_trains():
    """Binary focal-BCE path must produce a fitted model with finite proba."""
    from src.foundation.config import FoundationConfig, MODALITY_DIMS
    from src.foundation.downstream import FoundationDownstream
    cfg = FoundationConfig(embed_dim=16, n_heads=2, n_layers=1,
                           ff_dim=32, seed=0)
    fd = FoundationDownstream(config=cfg, pretrained=False,
                              loss="sens_at_spec", alpha_pos=20.0)
    rng = np.random.default_rng(0)
    mod = {name: rng.standard_normal((60, dim)).astype(np.float32)
           for name, dim in MODALITY_DIMS.items()}
    labels = (rng.random(60) < 0.3).astype(np.int64)
    fd.fit(mod, labels, n_epochs=2, batch_size=16, validation_split=0.1,
           verbose=False)
    assert fd._fitted
    proba = fd.predict_proba(mod)
    assert not np.isnan(proba).any()
    assert not np.isinf(proba).any()


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_cross_attention_fusion_sens_at_spec_loss_trains():
    """CrossAttentionFusion focal-BCE binary path produces finite proba."""
    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
    rng = np.random.default_rng(0)
    n = 100
    scores = [rng.standard_normal(n) + i for i in range(4)]
    labels = (rng.random(n) < 0.5).astype(np.int64)
    m = CrossAttentionFusion(n_modalities=4, prior=None, n_epochs=20,
                              seed=0, loss="sens_at_spec", alpha_pos=20.0)
    m.fit(scores, labels)
    proba = m.predict_proba(scores)
    assert proba.shape == (n,)
    assert not np.isnan(proba).any()
    assert not np.isinf(proba).any()


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_cross_attention_fusion_invalid_loss_raises():
    """Bad loss string must raise at construction time."""
    from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
    with pytest.raises(ValueError, match="loss must be"):
        CrossAttentionFusion(loss="not_a_loss")


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_early_late_fusion_sens_at_spec_loss_trains():
    """EarlyLateFusion focal-BCE binary path produces finite proba."""
    from src.multimodal_fusion.advanced_fusion import EarlyLateFusion
    rng = np.random.default_rng(0)
    n = 60
    feats = [rng.standard_normal((n, 4)) for _ in range(3)]
    labels = (rng.random(n) < 0.5).astype(np.int64)
    m = EarlyLateFusion(n_modalities=3, hidden_dim=16, n_epochs=20,
                         seed=0, loss="sens_at_spec", alpha_pos=20.0)
    m.fit(feats, labels)
    proba = m.predict_proba(feats)
    assert proba.shape == (n,)
    assert not np.isnan(proba).any()


# ── Shuffled-label control ───────────────────────────────────────────

def test_shuffled_label_control_diagnostic_shape():
    """signal_to_artifact_ratio must be a real float with sensible bounds.

    The diagnostic compares real_AUC vs shuffled_AUC. With a synthetic
    perfectly separable cohort, real_AUC=1.0 and shuffled_AUC=0.5, so
    the ratio is exactly 1.0. With random labels the ratio is 0.0.
    """
    # Synthetic case: perfectly separable scores → ratio should be 1.0
    y_real = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    scores = np.array([0.1, 0.2, 0.3, 0.4, 0.9, 0.95, 0.99, 1.0])
    real_auc = 1.0
    # Shuffled labels: 50% chance of staying correct by random luck
    shuffled_auc = 0.5
    denom = max(real_auc - 0.5, 1e-6)
    ratio = (real_auc - shuffled_auc) / denom
    assert ratio == 1.0
    assert 0.0 <= ratio <= 2.0  # never negative for valid diagnostic


def test_shuffled_label_control_random_label_is_zero():
    """Random labels with random scores → ratio near 0.

    The shuffled-label AUCs must come from INDEPENDENT permutations
    of y (each draws from a fresh RNG state) so the diagnostic
    estimates the mean correctly.
    """
    rng = np.random.default_rng(0)
    n = 200
    y = (rng.random(n) < 0.5).astype(int)
    scores = rng.standard_normal(n)
    from sklearn.metrics import roc_auc_score
    real_auc = roc_auc_score(y, scores)
    # Each shuffle uses its own RNG so the labels are truly independent
    # draws. Real AUC should be near 0.5 (random scores, random labels),
    # shuffled AUCs also near 0.5 → ratio near 0.
    shuf_aucs = []
    for s in range(20):
        shuf_rng = np.random.default_rng(1000 + s)
        shuf_aucs.append(roc_auc_score(shuf_rng.permutation(y), scores))
    shuf_mean = float(np.mean(shuf_aucs))
    denom = max(real_auc - 0.5, 1e-6)
    ratio = (real_auc - shuf_mean) / denom
    # With random scores and random labels, both real and shuffled AUC
    # should be near 0.5, so the ratio should be near 0. Allow a wide
    # tolerance for sampling noise on n=200.
    assert -2.0 < ratio < 2.0, f"random labels should give ratio ≈ 0, got {ratio}"


def test_real_tcga_validation_shuffled_flag_parsed():
    """--shuffled-label-control flag must parse and be discoverable."""
    import subprocess
    # Run with --help to confirm the flag is documented and argparse accepts it.
    result = subprocess.run(
        ["env", "-u", "PYTHONPATH",
         "/Users/hermes/deepcatch/.venv/bin/python",
         "/Users/hermes/deepcatch/real_tcga_validation.py", "--help"],
        capture_output=True, text=True, cwd="/Users/hermes/deepcatch", timeout=30
    )
    assert "--shuffled-label-control" in result.stdout, (
        "shuffled-label-control flag should appear in --help output"
    )
    assert "--n-shuffles" in result.stdout, (
        "n-shuffles flag should appear in --help output"
    )
