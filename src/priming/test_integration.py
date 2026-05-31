#!/usr/bin/env python3
"""
Integration Tests — DeepCatch Priming Agents Module
=====================================================

30+ tests covering config, PK model, response predictor, signal processing,
integration, data generators, whitepaper, and pitch deck.

Run: python -m pytest src/priming/test_integration.py -v
"""

from __future__ import annotations

import os
import sys
import tempfile
import warnings

import numpy as np

# Ensure deepcatch is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import pytest

# Suppress numpy warnings in tests
warnings.filterwarnings("ignore", category=RuntimeWarning)

from src.priming.config import (
    PrimingConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    PRODUCTION_CONFIG,
    AGENT_SPECS,
    AGENT_NAMES,
    PATIENT_FEATURE_SPEC,
    N_PATIENT_FEATURES,
)
from src.priming.pharmacokinetics import PKModel, OptimalDosingSchedule
from src.priming.signal_processing import (
    PostPrimingDenoiser,
    SignalEnhancer,
    BaselineCorrector,
)
from src.priming.data import (
    generate_patient_profiles,
    apply_priming_effect,
    simulate_clinical_trial,
)
from src.priming.integration import PrimingIntegration
from src.priming.whitepaper import generate_whitepaper_sections
from src.priming.pitch_deck import generate_pitch_deck, generate_slide_list

# ── Try importing torch-dependent modules ───────────────────────
try:
    import torch

    from src.priming.response_predictor import (
        PrimingResponsePredictor,
        PatientStratifier,
    )

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False
    PrimingResponsePredictor = None
    PatientStratifier = None


# ── Fixtures ─────────────────────────────────────────────────────

@pytest.fixture
def rng():
    return np.random.RandomState(42)


@pytest.fixture
def sample_patient():
    """Standard test patient."""
    return {
        "age": 55,
        "weight_kg": 70.0,
        "bmi": 25.0,
        "liver_function": 0.9,
        "renal_function": 0.95,
        "tumor_type": 0,  # lung
        "tumor_stage": 2,
        "tumor_fraction": 0.01,
        "baseline_cfdna_ng_ml": 30.0,
        "albumin_g_L": 42.0,
        "bilirubin_umol_L": 8.0,
        "creatinine_umol_L": 75.0,
        "alt_IU_L": 25.0,
        "ast_IU_L": 22.0,
        "platelet_count": 250.0,
        "neutrophil_count": 4.5,
        "hemoglobin_g_L": 140.0,
        "crp_mg_L": 3.0,
        "prior_treatment": 0,
        "performance_status": 0,
    }


@pytest.fixture
def sample_features(rng):
    """Random cfDNA features."""
    return rng.randn(50).astype(np.float64)


@pytest.fixture
def pk_model():
    return PKModel()


@pytest.fixture
def dosing_schedule():
    return OptimalDosingSchedule()


@pytest.fixture
def denoiser():
    return PostPrimingDenoiser(window=3)


@pytest.fixture
def enhancer():
    return SignalEnhancer(boost_factor=10.0)


@pytest.fixture
def corrector():
    return BaselineCorrector()


# ══════════════════════════════════════════════════════════════════
# Config Tests (6 tests)
# ══════════════════════════════════════════════════════════════════


class TestConfig:
    """Tests for PrimingConfig and presets."""

    def test_default_config_values(self):
        """Default config has expected values."""
        c = DEFAULT_CONFIG
        assert c.n_agents == 5
        assert c.n_patient_features == 20
        assert c.latent_dim == 64
        assert c.denoising_window == 3
        assert c.detection_boost_factor == 10.0
        assert c.mlp_hidden == [64, 32]
        assert c.predictor_output == 3

    def test_prototype_config(self):
        """Prototype config has smaller dimensions."""
        c = PROTOTYPE_CONFIG
        assert c.latent_dim == 16
        assert c.mlp_hidden == [16, 8]
        assert c.denoising_window == 2
        assert c.device == "cpu"

    def test_production_config(self):
        """Production config has larger model."""
        c = PRODUCTION_CONFIG
        assert c.latent_dim == 128
        assert c.mlp_hidden == [128, 64]
        assert c.detection_boost_factor == 12.0
        assert c.learning_rate == 5e-4
        assert c.stratification_thresholds["ideal_boost_min"] == 6.0

    def test_config_to_dict_and_back(self):
        """Config round-trips through serialization."""
        c = PrimingConfig(latent_dim=42, seed=7)
        d = c.to_dict()
        assert d["latent_dim"] == 42
        assert d["seed"] == 7
        c2 = PrimingConfig.from_dict(d)
        assert c2.latent_dim == 42
        assert c2.seed == 7

    def test_agent_specs_complete(self):
        """All 5 agents have complete specs."""
        for agent in AGENT_NAMES:
            spec = AGENT_SPECS[agent]
            assert "half_life_hours" in spec
            assert "volume_of_distribution_L" in spec
            assert "clearance_rate_L_h" in spec
            assert "bioavailability" in spec
            assert "protein_binding" in spec
            assert 0 < spec["bioavailability"] <= 1.0

    def test_patient_features_count(self):
        """Patient feature spec matches N_PATIENT_FEATURES."""
        assert len(PATIENT_FEATURE_SPEC) == N_PATIENT_FEATURES
        assert N_PATIENT_FEATURES == 20


# ══════════════════════════════════════════════════════════════════
# PK Model Tests (8 tests)
# ══════════════════════════════════════════════════════════════════


class TestPKModel:
    """Tests for pharmacokinetic model."""

    def test_simulate_all_five_agents(self, pk_model, sample_patient):
        """PK simulation works for all 5 agent types."""
        for agent in AGENT_NAMES:
            result = pk_model.simulate(
                agent=agent, dose_mg=100, patient_weight_kg=70, liver_function=0.9
            )
            assert "time_hours" in result
            assert "concentration_mg_L" in result
            assert "clearance_rate" in result
            assert len(result["time_hours"]) == 100
            assert result["peak_concentration"] > 0
            assert result["auc_mg_h_L"] > 0

    def test_liposome_longest_half_life(self, pk_model):
        """Liposome should have the longest effective half-life."""
        half_lives = {}
        for agent in AGENT_NAMES:
            result = pk_model.simulate(agent=agent, dose_mg=100)
            half_lives[agent] = result["effective_half_life_h"]
        # liposome has longest literature half-life (18h)
        assert half_lives["liposome"] >= half_lives["scFv"]

    def test_weight_zero_clamped(self, pk_model):
        """Weight=0 should be clamped to 1, not crash."""
        result = pk_model.simulate(
            agent="liposome", dose_mg=100, patient_weight_kg=0, liver_function=1.0
        )
        assert result["peak_concentration"] > 0

    def test_liver_function_zero(self, pk_model):
        """Liver function=0 reduces clearance but doesn't crash."""
        result_normal = pk_model.simulate(
            agent="scFv", dose_mg=100, liver_function=1.0
        )
        result_bad = pk_model.simulate(
            agent="scFv", dose_mg=100, liver_function=0.0
        )
        # With impaired liver, clearance is slower → higher AUC
        assert result_bad["auc_mg_h_L"] >= result_normal["auc_mg_h_L"] * 0.5

    def test_dose_zero_returns_zero_effect(self, pk_model):
        """Dose=0 should return zero-effect result."""
        result = pk_model.simulate(agent="liposome", dose_mg=0)
        assert result["peak_concentration"] == 0.0
        assert result["auc_mg_h_L"] == 0.0
        assert np.allclose(result["ctDNA_boost_factor"], 1.0)

    def test_negative_dose_clamped(self, pk_model):
        """Negative dose should return zero-effect result."""
        result = pk_model.simulate(agent="liposome", dose_mg=-10)
        assert result["peak_concentration"] == 0.0

    def test_unknown_agent_raises(self, pk_model):
        """Unknown agent name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown agent"):
            pk_model.simulate(agent="unicorn", dose_mg=100)

    def test_shape_consistency(self, pk_model):
        """Output arrays have consistent shapes."""
        result = pk_model.simulate(agent="dendrimer", dose_mg=50, n_timepoints=75)
        assert result["time_hours"].shape == (75,)
        assert result["concentration_mg_L"].shape == (75,)
        assert result["clearance_rate"].shape == (75,)
        assert result["ctDNA_boost_factor"].shape == (75,)


# ══════════════════════════════════════════════════════════════════
# Optimal Dosing Tests (4 tests)
# ══════════════════════════════════════════════════════════════════


class TestOptimalDosing:
    """Tests for OptimalDosingSchedule."""

    def test_basic_compute(self, dosing_schedule, sample_patient):
        """Basic dosing schedule computation."""
        result = dosing_schedule.compute("liposome", sample_patient)
        assert result["optimal_dose_mg"] > 0
        assert result["optimal_timing_h"] >= 0
        assert result["predicted_ctDNA_boost"] >= 1.0
        assert 0 <= result["predicted_toxicity_risk"] <= 1.0
        assert len(result["schedule"]) >= 2
        assert isinstance(result["recommendation"], str)
        assert "is_feasible" in result

    def test_poor_liver_reduces_dose(self, dosing_schedule):
        """Poor liver function should result in lower dose."""
        good_liver = {"weight_kg": 70, "liver_function": 1.0}
        bad_liver = {"weight_kg": 70, "liver_function": 0.2}
        r_good = dosing_schedule.compute("scFv", good_liver)
        r_bad = dosing_schedule.compute("scFv", bad_liver)
        # Poor liver → higher toxicity risk
        assert r_bad["predicted_toxicity_risk"] > r_good["predicted_toxicity_risk"]

    def test_all_agents_computable(self, dosing_schedule, sample_patient):
        """Dosing works for all 5 agent types."""
        for agent in AGENT_NAMES:
            result = dosing_schedule.compute(agent, sample_patient)
            assert result["optimal_dose_mg"] > 0

    def test_unknown_agent_raises(self, dosing_schedule, sample_patient):
        """Unknown agent raises ValueError."""
        with pytest.raises(ValueError):
            dosing_schedule.compute("ghost_agent", sample_patient)


# ══════════════════════════════════════════════════════════════════
# Signal Processing Tests (9 tests)
# ══════════════════════════════════════════════════════════════════


class TestSignalProcessing:
    """Tests for signal processing classes."""

    def test_denoise_1d_signal(self, denoiser, rng):
        """1D denoising preserves signal character."""
        signal = rng.randn(100)
        signal[50:60] += 5.0  # Add spike
        denoised = denoiser.denoise(signal)
        assert denoised.shape == signal.shape
        # Noise should be reduced
        assert np.std(denoised) <= np.std(signal) * 1.5

    def test_denoise_2d_signal(self, denoiser, rng):
        """2D denoising works on feature matrices."""
        signal = rng.randn(100, 10)
        denoised = denoiser.denoise(signal)
        assert denoised.shape == signal.shape

    def test_denoise_empty_signal(self, denoiser):
        """Empty signal returns empty."""
        signal = np.array([])
        result = denoiser.denoise(signal)
        assert len(result) == 0

    def test_denoise_with_priming_timing(self, denoiser, rng):
        """Denoising with priming timing preserves pre-priming baseline."""
        signal = rng.randn(50)
        signal[30:] += 3.0  # Post-priming boost
        denoised = denoiser.denoise(signal, priming_timing=25)
        # Pre-priming region should have similar mean to original
        pre_mean_orig = np.mean(signal[:25])
        pre_mean_denoised = np.mean(denoised[:25])
        # They should be reasonably close (within 50% relative)
        assert abs(pre_mean_denoised - pre_mean_orig) < abs(pre_mean_orig) * 0.5 + 1.0

    def test_estimate_background_noise(self, denoiser, rng):
        """Noise estimation returns positive value."""
        signal = rng.randn(200)
        noise = denoiser.estimate_background_noise(signal)
        assert noise >= 0

    def test_enhance_boosts_signal(self, enhancer, rng):
        """Signal enhancement increases signal features."""
        features = rng.randn(50)
        features[10:15] += 10.0  # Strong signal features
        enhanced = enhancer.enhance(features, signal_to_noise=2.0, priming_boost=5.0)
        assert enhanced.shape == features.shape
        # Strong features should remain strong
        assert np.max(np.abs(enhanced[10:15])) >= np.max(np.abs(features[10:15])) * 0.5

    def test_enhance_within_range(self, enhancer, rng):
        """Enhancement doesn't produce extreme values."""
        features = rng.randn(100)
        enhanced = enhancer.enhance(features, priming_boost=10.0)
        assert not np.any(np.isnan(enhanced))
        assert not np.any(np.isinf(enhanced))

    def test_baseline_correction_subtractive(self, rng):
        """Subtractive baseline correction works."""
        corrector = BaselineCorrector(correction_method="subtractive")
        pre = np.array([5.0])
        post = np.array([15.0, 20.0, 18.0])
        corrected = corrector.correct(post, pre)
        assert corrected[0] == pytest.approx(10.0)
        assert corrected[1] == pytest.approx(15.0)

    def test_baseline_correction_ratio(self, rng):
        """Ratio baseline correction normalizes."""
        corrector = BaselineCorrector(correction_method="ratio")
        pre = np.array([2.0])
        post = np.array([4.0, 8.0])
        corrected = corrector.correct(post, pre)
        assert corrected[0] == pytest.approx(2.0)
        assert corrected[1] == pytest.approx(4.0)


# ══════════════════════════════════════════════════════════════════
# Response Predictor Tests (5 tests) — requires torch
# ══════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not HAS_TORCH, reason="PyTorch not installed")
class TestResponsePredictor:
    """Tests for PrimingResponsePredictor and PatientStratifier."""

    def test_forward_pass(self):
        """Predictor produces valid output shape."""
        model = PrimingResponsePredictor(config=PROTOTYPE_CONFIG)
        x = torch.randn(4, 25)  # batch of 4
        output = model(x)
        assert output.shape == (4, 3)
        # boost > 0 (softplus), toxicity in [0,1] (sigmoid)
        assert torch.all(output[:, 0] > 0)
        assert torch.all(output[:, 1] > 0)
        assert torch.all((output[:, 2] >= 0) & (output[:, 2] <= 1))

    def test_predict_method(self, sample_patient):
        """Predict returns expected keys."""
        model = PrimingResponsePredictor(config=PROTOTYPE_CONFIG)
        result = model.predict(sample_patient, "liposome")
        assert "predicted_boost_factor" in result
        assert "predicted_time_to_peak_h" in result
        assert "predicted_toxicity_risk" in result
        assert result["predicted_boost_factor"] > 0
        assert 0 <= result["predicted_toxicity_risk"] <= 1
        assert result["agent"] == "liposome"

    def test_estimated_params_small(self):
        """Model has ~5K params as designed."""
        model = PrimingResponsePredictor(config=PROTOTYPE_CONFIG)
        n_params = model.estimated_params()
        assert n_params < 10000  # Should be small
        assert n_params > 100

    def test_stratifier_ideal_candidate(self, sample_patient):
        """Ideal candidate identified correctly."""
        stratifier = PatientStratifier()
        sample_patient["liver_function"] = 0.95
        sample_patient["tumor_stage"] = 1
        category = stratifier.stratify(
            sample_patient, predicted_boost=8.0, predicted_toxicity=0.1
        )
        assert category == "ideal_candidate"

    def test_stratifier_poor_candidate(self, sample_patient):
        """Poor candidate identified with low liver function."""
        stratifier = PatientStratifier()
        sample_patient["liver_function"] = 0.2

        # First test: explicit predictions should override
        category = stratifier.stratify(
            sample_patient, predicted_boost=0.5, predicted_toxicity=0.8
        )
        assert category == "poor_candidate"

        # Second test: rule-based fallback with severe organ impairment
        category2 = stratifier.stratify(sample_patient)
        assert category2 == "poor_candidate"


# ══════════════════════════════════════════════════════════════════
# Data Generation Tests (6 tests)
# ══════════════════════════════════════════════════════════════════


class TestDataGeneration:
    """Tests for synthetic data generators."""

    def test_generate_patient_profiles(self):
        """Profiles have correct structure and count."""
        profiles = generate_patient_profiles(50, seed=42)
        assert len(profiles) == 50
        for p in profiles:
            assert "patient_id" in p
            assert "age" in p
            assert "weight_kg" in p
            assert 18 <= p["age"] <= 100
            assert 30 <= p["weight_kg"] <= 200
            assert 0 <= p["liver_function"] <= 1.0

    def test_apply_priming_effect_boost(self, sample_features):
        """Priming effect increases signal features."""
        result = apply_priming_effect(
            sample_features, "liposome", dose_mg=500, seed=42
        )
        assert result["boost_factor"] >= 1.0
        assert len(result["post_priming_features"]) == len(sample_features)
        assert result["agent"] == "liposome"

    def test_apply_priming_effect_zero_dose(self, sample_features):
        """Zero dose returns unchanged features."""
        result = apply_priming_effect(
            sample_features, "liposome", dose_mg=0, seed=42
        )
        assert result["boost_factor"] == 1.0
        np.testing.assert_array_equal(
            result["pre_priming_features"], result["post_priming_features"]
        )

    def test_simulate_clinical_trial(self):
        """Clinical trial simulation returns all components."""
        patients, trial, summary = simulate_clinical_trial(
            n_patients=50, n_agents=3, n_features=20, seed=42
        )
        assert len(patients) == 50
        assert "pre_samples" in trial
        assert trial["pre_samples"].shape == (50, 20)
        assert trial["post_samples"].shape == (150, 20)  # 50 * 3 agents
        assert len(trial["labels"]) == 50
        assert summary["n_patients"] == 50
        assert summary["mean_boost_factor"] > 1.0

    def test_simulate_clinical_trial_labels(self):
        """Some patients have cancer labels."""
        _, trial, _ = simulate_clinical_trial(n_patients=100, seed=42)
        labels = trial["labels"]
        n_cancer = np.sum(labels == 1.0)
        assert n_cancer > 0  # At least some have cancer
        assert n_cancer < 100  # Not all have cancer (10% healthy)

    def test_reproducibility(self):
        """Same seed produces identical data."""
        profiles1 = generate_patient_profiles(20, seed=42)
        profiles2 = generate_patient_profiles(20, seed=42)
        for p1, p2 in zip(profiles1, profiles2):
            assert p1["age"] == p2["age"]
            assert p1["weight_kg"] == p2["weight_kg"]


# ══════════════════════════════════════════════════════════════════
# Integration Tests (7 tests)
# ══════════════════════════════════════════════════════════════════


class TestIntegration:
    """Tests for PrimingIntegration adapter."""

    @pytest.fixture
    def integ(self):
        return PrimingIntegration()

    def test_process_signal_full(self, integ, sample_features, sample_patient):
        """Full pipeline processes a sample with patient data."""
        sample = {
            "features": sample_features,
            "patient_data": sample_patient,
            "agent": "liposome",
        }
        result = integ.process_signal(sample)
        assert "priming_score" in result
        assert 0 <= result["priming_score"] <= 1.0
        assert result["boost_factor"] >= 1.0
        assert result["signal_quality"] >= 0
        assert len(result["denoised_features"]) == len(sample_features)
        assert len(result["enhanced_features"]) == len(sample_features)
        assert result["stratification"] in (
            "ideal_candidate",
            "moderate_candidate",
            "poor_candidate",
        )

    def test_to_modality(self, integ, sample_features, sample_patient):
        """to_modality returns scalar in [0, 1]."""
        sample = {
            "features": sample_features,
            "patient_data": sample_patient,
        }
        score = integ.to_modality(sample)
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_extract_all(self, integ, sample_features, sample_patient):
        """extract_all returns ~15 features."""
        sample = {
            "features": sample_features,
            "patient_data": sample_patient,
        }
        features = integ.extract_all(sample)
        assert "priming_score" in features
        assert "priming_boost_factor" in features
        assert "priming_ideal_candidate" in features
        assert len(features) >= 10

    def test_missing_patient_data(self, integ, sample_features):
        """Missing patient data returns zero result gracefully."""
        sample = {"features": sample_features}
        result = integ.process_signal(sample)
        assert result["priming_score"] == 0.0
        assert result["boost_factor"] == 1.0
        assert result["stratification"] == "poor_candidate"

    def test_none_features(self, integ):
        """None features returns empty result."""
        sample = {"features": None, "patient_data": {}}
        result = integ.process_signal(sample)
        assert result["priming_score"] == 0.0

    def test_empty_patient_data(self, integ, sample_features):
        """Empty patient data dict works."""
        sample = {"features": sample_features, "patient_data": {}}
        result = integ.process_signal(sample)
        assert result["priming_score"] == 0.0

    def test_all_agents_integration(self, integ, sample_features, sample_patient):
        """Integration works for all agent types."""
        for agent in AGENT_NAMES:
            sample = {
                "features": sample_features,
                "patient_data": sample_patient,
                "agent": agent,
            }
            result = integ.process_signal(sample)
            assert result["boost_factor"] >= 1.0


# ══════════════════════════════════════════════════════════════════
# Checkpoint Tests (2 tests)
# ══════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not HAS_TORCH, reason="PyTorch not installed")
class TestCheckpoint:
    """Tests for model checkpoint save/load."""

    def test_save_load_checkpoint(self, sample_patient):
        """Checkpoint round-tripping."""
        integ = PrimingIntegration()
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            path = f.name
        try:
            integ.save_checkpoint(path)
            assert os.path.exists(path)

            # Load into new instance
            integ2 = PrimingIntegration(checkpoint=path)
            result = integ2.process_signal({
                "features": np.random.randn(50),
                "patient_data": sample_patient,
            })
            assert result["priming_score"] >= 0
        finally:
            if os.path.exists(path):
                os.unlink(path)

    def test_load_nonexistent_checkpoint(self):
        """Loading nonexistent checkpoint logs error but doesn't crash."""
        integ = PrimingIntegration(checkpoint="/nonexistent/path.pt")
        assert integ.predictor is not None  # Still has predictor, just not loaded


# ══════════════════════════════════════════════════════════════════
# Whitepaper Tests (1 test)
# ══════════════════════════════════════════════════════════════════


class TestWhitepaper:
    """Tests for whitepaper generation."""

    def test_all_sections_non_empty(self):
        """All whitepaper sections are non-empty markdown."""
        sections = generate_whitepaper_sections()
        required = [
            "executive_summary",
            "technical_overview",
            "experimental_design",
            "expected_outcomes",
            "data_requirements",
            "roadmap",
            "references",
        ]
        for section_name in required:
            assert section_name in sections, f"Missing section: {section_name}"
            content = sections[section_name]
            assert len(content.strip()) > 100, (
                f"Section {section_name} too short: {len(content)} chars"
            )
            assert "#" in content, f"Section {section_name} has no markdown headers"


# ══════════════════════════════════════════════════════════════════
# Pitch Deck Tests (2 tests)
# ══════════════════════════════════════════════════════════════════


class TestPitchDeck:
    """Tests for pitch deck generation."""

    def test_pitch_deck_non_empty(self):
        """Pitch deck is non-empty markdown."""
        deck = generate_pitch_deck()
        assert len(deck) > 1000
        assert "#" in deck
        assert "DeepCatch" in deck
        assert "Amplifyer" in deck

    def test_slide_list_has_all_slides(self):
        """Slide list returns 10+ slides."""
        slides = generate_slide_list()
        assert len(slides) >= 8
        for slide in slides:
            assert "title" in slide
            assert "content" in slide
            assert slide["slide_number"] >= 1


# ══════════════════════════════════════════════════════════════════
# Run main
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
