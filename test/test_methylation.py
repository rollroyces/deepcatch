"""
cfrNA persistence coupling + methylation reader smoke tests.

These tests guard against regression in the v2.3.0-port of
``deepcatch-methylation`` v0.10/0.11. They verify the typed-dict contract
and the basic dose --> residual ratio -- no real patient data is touched.
"""

import pytest

from src.methylation import cfrna_coupling
from src.methylation import methylation_reader


def test_cfrna_coupling_basic_residual_ratio():
    """At sampling_time = half_life, observed/expected ratio is 1.0 for a
    correctly-dosed sample."""
    obs = cfrna_coupling.couple_cfrna_to_cdm(
        prediction={
            "predicted_half_life_hours": 12.0,
            "persistence_category": "medium",
            "target_tissue": "liver",
            "codon_optimization": "balanced",
            "feature_breakdown": {},
        },
        observed_cfRNA_concentration_pg_ml=500.0,  # half of 1000 at t=0
        sampling_time_post_dose_hours=12.0,       # exactly half_life
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    # At t = half_life, expected residual = 0.5; observed = 0.5; ratio = 1.0
    assert abs(obs.expected_residual_fraction - 0.5) < 1e-6
    assert abs(obs.observed_residual_fraction - 0.5) < 1e-6
    assert abs(obs.observed_vs_expected_ratio - 1.0) < 1e-6
    assert obs.anomaly_flag is False
    assert obs.dosing_recommendation == "maintain"


def test_cfrna_coupling_anomaly_low():
    """A ratio < ANOMALY_RATIO_LOW flags the observation as anomalous."""
    obs = cfrna_coupling.couple_cfrna_to_cdm(
        prediction={
            "predicted_half_life_hours": 12.0,
            "persistence_category": "medium",
            "target_tissue": "liver",
            "codon_optimization": "balanced",
            "feature_breakdown": {},
        },
        observed_cfRNA_concentration_pg_ml=100.0,  # 0.1 of 1000 -- much less than expected
        sampling_time_post_dose_hours=12.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.observed_residual_fraction == 0.1
    assert obs.observed_vs_expected_ratio == 0.2
    assert obs.anomaly_flag is True
    assert obs.dosing_recommendation == "increase"


def test_cfrna_coupling_hold_below_min_half_life():
    """A prediction with half_life < MIN_HALFLIFE_HOURS forces a 'hold'."""
    obs = cfrna_coupling.couple_cfrna_to_cdm(
        prediction={
            "predicted_half_life_hours": 0.5,  # below 1.0 floor
            "persistence_category": "low",
            "target_tissue": "liver",
            "codon_optimization": "balanced",
            "feature_breakdown": {},
        },
        observed_cfRNA_concentration_pg_ml=500.0,
        sampling_time_post_dose_hours=12.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.dosing_recommendation == "hold"


def test_methylation_reader_liver_synthetic():
    """``make_synthetic_liver_beta_values`` produces above-threshold marker values."""
    betas = methylation_reader.make_synthetic_liver_beta_values(tissue="liver")
    assert len(betas) >= 1
    for v in betas.values():
        assert v >= 0.5  # above the tissue-call threshold


def test_methylation_reader_build_persistence_input():
    """``build_cfrna_persistence_input`` returns a TypedDict-compatible shape."""
    inp = methylation_reader.build_cfrna_persistence_input({
        "sample_id": "test",
        "tissue": "liver",
        "beta_values": methylation_reader.make_synthetic_liver_beta_values("liver"),
        "cfrna_concentration_pg_ml": 200.0,
        "sampling_time_post_dose_hours": 12.0,
        "dose_cfrna_at_t0_pg_ml": 1000.0,
    })
    assert "predicted_half_life_hours" in inp
    assert "persistence_category" in inp
    assert inp["target_tissue"] == "liver"


def test_cfrna_coupling_invalid_category_raises():
    """An unknown ``persistence_category`` raises ``ValueError``."""
    with pytest.raises(ValueError):
        cfrna_coupling.couple_cfrna_to_cdm(
            prediction={
                "predicted_half_life_hours": 12.0,
                "persistence_category": "ultra",
                "target_tissue": "liver",
                "codon_optimization": "balanced",
                "feature_breakdown": {},
            },
            observed_cfRNA_concentration_pg_ml=500.0,
            sampling_time_post_dose_hours=12.0,
            dose_cfRNA_at_t0_pg_ml=1000.0,
        )
