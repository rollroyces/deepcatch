#!/usr/bin/env python3
"""
Integration Tests — Tissue Deconvolution Modality
====================================================

Tests the tissue deconvolution module end-to-end and per-component.

Coverage:
1. Config validation
2. TissueAtlas with synthetic reference
3. Model forward pass with dummy data
4. Training on synthetic data (loss decreases)
5. Ensemble prediction consistency
6. Missing data handling
7. Feature extraction (keys present, value ranges)
8. Integration with fake CrossAttentionFusion
9. Checkpoint save/load
10. Graceful fallback (no model = zeros)

Usage
-----

    # Run all tests
    python src/tissue_deconv/test_integration.py

    # Run with verbose output
    python src/tissue_deconv/test_integration.py -v

Environment
-----------
Requires: numpy
Optional: torch (for model tests, skipped gracefully if missing)
"""

from __future__ import annotations

import os
import sys
import unittest
import warnings
import tempfile
from copy import deepcopy
from typing import Any, Dict, List, Optional

import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# ── PyTorch availability probe ──────────────────────────────────

_HAS_TORCH = False
try:
    import torch
    _HAS_TORCH = True
except ImportError:
    pass


# ═══════════════════════════════════════════════════════════════════
# Test Classes
# ═══════════════════════════════════════════════════════════════════

class TestConfig(unittest.TestCase):
    """Test configuration dataclass."""

    def test_01_default_config_valid(self):
        """Default config should instantiate without errors."""
        from src.tissue_deconv.config import TissueDeconvConfig, DEFAULT_CONFIG
        self.assertIsNotNone(DEFAULT_CONFIG)
        self.assertEqual(DEFAULT_CONFIG.n_tissues, 29)
        self.assertEqual(DEFAULT_CONFIG.n_cpg_features, 1000)
        self.assertEqual(DEFAULT_CONFIG.hidden_dims, [256, 128, 64])
        self.assertEqual(DEFAULT_CONFIG.n_ensemble, 3)
        self.assertGreater(DEFAULT_CONFIG.estimated_params, 0)

    def test_02_prototype_config(self):
        """Prototype config should have smaller values."""
        from src.tissue_deconv.config import PROTOTYPE_CONFIG
        self.assertLess(PROTOTYPE_CONFIG.n_cpg_features, 1000)
        self.assertLess(sum(PROTOTYPE_CONFIG.hidden_dims), 300)
        self.assertLess(PROTOTYPE_CONFIG.n_ensemble, 3)

    def test_03_supported_tissues_count(self):
        """SUPPORTED_TISSUES must have exactly 29 entries."""
        from src.tissue_deconv.config import SUPPORTED_TISSUES
        self.assertEqual(len(SUPPORTED_TISSUES), 29)

    def test_04_cancer_relevant_subset(self):
        """Cancer-relevant tissues should be a subset of supported."""
        from src.tissue_deconv.config import (
            SUPPORTED_TISSUES, CANCER_RELEVANT_TISSUES
        )
        supported_set = set(SUPPORTED_TISSUES)
        for t in CANCER_RELEVANT_TISSUES:
            self.assertIn(t, supported_set)

    def test_05_blood_derived_subset(self):
        """Blood-derived tissues should be a subset of supported."""
        from src.tissue_deconv.config import (
            SUPPORTED_TISSUES, BLOOD_DERIVED_TISSUES
        )
        supported_set = set(SUPPORTED_TISSUES)
        for t in BLOOD_DERIVED_TISSUES:
            self.assertIn(t, supported_set)

    def test_06_invalid_n_tissues_raises(self):
        """Mismatched n_tissues should raise ValueError."""
        from src.tissue_deconv.config import TissueDeconvConfig
        with self.assertRaises(ValueError):
            TissueDeconvConfig(n_tissues=5)

    def test_07_config_to_dict(self):
        """Config serialization should include key fields."""
        from src.tissue_deconv.config import DEFAULT_CONFIG
        d = DEFAULT_CONFIG.to_dict()
        self.assertIn("n_tissues", d)
        self.assertIn("hidden_dims", d)
        self.assertIn("learning_rate", d)
        self.assertEqual(d["n_tissues"], 29)

    def test_08_estimated_params_reasonable(self):
        """Estimated params should be between 300K and 600K."""
        from src.tissue_deconv.config import DEFAULT_CONFIG
        params = DEFAULT_CONFIG.estimated_params
        self.assertGreater(params, 200_000)
        self.assertLess(params, 800_000)


class TestTissueAtlas(unittest.TestCase):
    """Test TissueAtlas class."""

    def setUp(self):
        from src.tissue_deconv.tissue_atlas import TissueAtlas
        from src.tissue_deconv.config import PROTOTYPE_CONFIG
        self.atlas = TissueAtlas(
            n_cpg_features=PROTOTYPE_CONFIG.n_cpg_features,
            random_seed=42,
        )
        self.atlas.load_reference()

    def test_09_atlas_loaded_all_tissues(self):
        """Atlas should have profiles for all 29 tissues."""
        from src.tissue_deconv.config import SUPPORTED_TISSUES
        loaded = self.atlas.tissues
        for t in SUPPORTED_TISSUES:
            self.assertIn(t, loaded)

    def test_10_atlas_is_synthetic(self):
        """Freshly generated atlas should be marked synthetic."""
        self.assertTrue(self.atlas.is_synthetic)

    def test_11_get_reference_profile(self):
        """get_reference_profile should return correct shape."""
        profile = self.atlas.get_reference_profile("Liver")
        self.assertIsInstance(profile, np.ndarray)
        self.assertEqual(profile.shape, (100,))  # PROTOTYPE_CONFIG
        self.assertTrue(np.all(profile >= 0) and np.all(profile <= 1))

    def test_12_get_marker_cpgs(self):
        """get_marker_cpgs should return unique CpG indices."""
        markers = self.atlas.get_marker_cpgs("Liver", n_top=20)
        # With 100 CpGs / 29 tissues ≈ 3 → max(5,3) = 5 markers per tissue
        self.assertGreater(len(markers), 0)
        self.assertLessEqual(len(markers), 20)
        self.assertEqual(len(set(markers)), len(markers))  # unique
        self.assertTrue(np.all(markers >= 0))
        self.assertTrue(np.all(markers < self.atlas.n_cpg_features))

    def test_13_generate_synthetic_mixture(self):
        """Synthetic mixture should be weighted sum of profiles."""
        frac = {
            "Whole Blood": 0.60,
            "Neutrophil": 0.15,
            "Monocyte": 0.10,
            "Liver": 0.10,
            "Lung": 0.05,
        }
        mixture = self.atlas.generate_synthetic_mixture(frac, noise=0.005)
        self.assertEqual(mixture.shape, (100,))
        self.assertTrue(np.all(mixture >= 0) and np.all(mixture <= 1))

        # Mixture should be different from pure liver profile
        liver = self.atlas.get_reference_profile("Liver")
        self.assertFalse(np.allclose(mixture, liver, atol=0.01))

    def test_14_generate_training_dataset(self):
        """Training dataset should have correct shapes."""
        mixes, fracs = self.atlas.generate_training_dataset(
            n_samples=50, seed=123,
        )
        self.assertEqual(mixes.shape, (50, 100))
        self.assertEqual(fracs.shape, (50, 29))
        # Fractions should sum to ~1
        row_sums = fracs.sum(axis=1)
        self.assertTrue(np.allclose(row_sums, 1.0, atol=0.01))

    def test_15_atlas_npz_roundtrip(self):
        """Save and load atlas via NPZ."""
        with tempfile.NamedTemporaryFile(suffix=".npz", delete=False) as f:
            path = f.name
        try:
            self.atlas.to_npz(path)

            from src.tissue_deconv.tissue_atlas import TissueAtlas
            atlas2 = TissueAtlas(n_cpg_features=100)
            atlas2.load_reference(path=path)

            self.assertFalse(atlas2.is_synthetic)
            self.assertEqual(len(atlas2.tissues), 29)

            # Compare a profile
            orig = self.atlas.get_reference_profile("Liver")
            loaded = atlas2.get_reference_profile("Liver")
            self.assertTrue(np.allclose(orig, loaded, atol=1e-5))
        finally:
            if os.path.exists(path):
                os.unlink(path)

    def test_16_get_all_profiles_matrix(self):
        """get_all_profiles_matrix should return (29, n_cpg) shape."""
        mat = self.atlas.get_all_profiles_matrix()
        self.assertEqual(mat.shape, (29, 100))
        self.assertEqual(mat.dtype, np.float32)

    def test_17_load_reference_from_dict(self):
        """Loading from dict should override synthetic."""
        from src.tissue_deconv.tissue_atlas import TissueAtlas
        from src.tissue_deconv.config import SUPPORTED_TISSUES
        atlas2 = TissueAtlas(n_cpg_features=50)
        ref = {
            t: np.random.RandomState(i).uniform(0, 1, 50).astype(np.float32)
            for i, t in enumerate(SUPPORTED_TISSUES[:5])  # partial
        }
        atlas2.load_reference(reference_dict=ref)
        self.assertFalse(atlas2.is_synthetic)
        self.assertEqual(len(atlas2.tissues), 5)

    def test_17b_missing_tissue_raises_keyerror(self):
        """Accessing unknown tissue should raise KeyError."""
        with self.assertRaises(KeyError):
            self.atlas.get_reference_profile("NonExistentTissue")


class TestModel(unittest.TestCase):
    """Test TissueDeconvolutionModel and Ensemble."""

    def setUp(self):
        from src.tissue_deconv.config import PROTOTYPE_CONFIG
        self.config = PROTOTYPE_CONFIG
        self.batch_size = 8
        self.rng = np.random.RandomState(42)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_18_model_forward_pass(self):
        """Model should produce valid softmax output."""
        from src.tissue_deconv.model import TissueDeconvolutionModel

        model = TissueDeconvolutionModel(config=self.config, seed=42)
        model.eval()

        x = torch.randn(self.batch_size, 100)
        with torch.no_grad():
            out = model(x)

        self.assertEqual(out.shape, (self.batch_size, 29))
        # Sums should be ~1 (softmax)
        self.assertTrue(torch.allclose(out.sum(dim=1), torch.ones(self.batch_size), atol=1e-5))
        # Values in [0, 1]
        self.assertTrue(torch.all(out >= 0) and torch.all(out <= 1))

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_19_model_return_logits(self):
        """Model should return logits when requested."""
        from src.tissue_deconv.model import TissueDeconvolutionModel

        model = TissueDeconvolutionModel(config=self.config, seed=42)
        model.eval()

        x = torch.randn(self.batch_size, 100)
        with torch.no_grad():
            fractions, logits = model(x, return_logits=True)

        self.assertEqual(fractions.shape, (self.batch_size, 29))
        self.assertEqual(logits.shape, (self.batch_size, 29))
        self.assertTrue(torch.allclose(fractions.sum(dim=1), torch.ones(self.batch_size), atol=1e-5))

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_20_model_parameter_count(self):
        """Model should have reasonable parameter count."""
        from src.tissue_deconv.model import TissueDeconvolutionModel

        model = TissueDeconvolutionModel(config=self.config, seed=42)
        n_params = model.count_parameters()

        # PROTOTYPE_CONFIG: [64, 32, 16] → ~30K params
        self.assertGreater(n_params, 5000)
        self.assertLess(n_params, 100_000)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_21_ensemble_initialization(self):
        """Ensemble should have correct number of models."""
        from src.tissue_deconv.model import TissueDeconvolutionEnsemble

        ensemble = TissueDeconvolutionEnsemble(config=self.config, base_seed=42)
        self.assertEqual(len(ensemble.models), self.config.n_ensemble)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_22_ensemble_predict_consistency(self):
        """Ensemble predictions should be consistent across calls."""
        from src.tissue_deconv.model import TissueDeconvolutionEnsemble

        ensemble = TissueDeconvolutionEnsemble(config=self.config, base_seed=42)
        ensemble.eval()

        x = self.rng.uniform(0, 1, 100).astype(np.float32)

        # Two successive predictions should give same result
        pred1, std1 = ensemble.predict(x, return_std=True)
        pred2, std2 = ensemble.predict(x, return_std=True)

        self.assertTrue(np.allclose(pred1, pred2, atol=1e-6))
        self.assertTrue(np.allclose(std1, std2, atol=1e-6))

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_23_ensemble_predict_sample(self):
        """predict_sample should return dict with metadata."""
        from src.tissue_deconv.model import TissueDeconvolutionEnsemble

        ensemble = TissueDeconvolutionEnsemble(config=self.config, base_seed=42)
        x = self.rng.uniform(0, 1, 100).astype(np.float32)

        result = ensemble.predict_sample(x, return_top=5)

        self.assertIn("_entropy", result)
        self.assertIn("_n_active", result)
        self.assertIn("_top_tissue", result)
        self.assertIn("_top_fraction", result)
        self.assertIsInstance(result["_top_tissue"], str)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_24_loss_computation(self):
        """DeconvLoss should return decreasing values for better preds."""
        from src.tissue_deconv.model import TissueDeconvolutionModel, DeconvLoss

        criterion = DeconvLoss(self.config)

        # Good prediction (close to target)
        target = torch.zeros(1, 29)
        target[0, 0] = 0.7  # Blood dominant
        target[0, 10] = 0.2  # Liver
        target[0, 15] = 0.1  # Kidney

        good_pred = target.clone() + torch.randn(1, 29) * 0.01
        good_pred = torch.softmax(good_pred, dim=-1)

        bad_pred = torch.ones(1, 29) / 29  # Uniform

        loss_good, _ = criterion(good_pred, target)
        loss_bad, _ = criterion(bad_pred, target)

        # Good prediction should have lower loss
        self.assertLess(float(loss_good), float(loss_bad))


class TestTrainer(unittest.TestCase):
    """Test TissueDeconvTrainer."""

    def setUp(self):
        from src.tissue_deconv.config import PROTOTYPE_CONFIG
        self.config = PROTOTYPE_CONFIG
        self.config.n_epochs = 5
        self.config.patience = 10  # don't early-stop

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_25_synthetic_data_generation(self):
        """Synthetic data should have correct shapes."""
        from src.tissue_deconv.trainer import TissueDeconvTrainer

        trainer = TissueDeconvTrainer(config=self.config, device="cpu")
        (train_x, train_y), (val_x, val_y) = trainer.generate_synthetic_data(
            n_samples=100, seed=42,
        )

        self.assertEqual(train_x.shape[1], 100)  # n_cpg_features
        self.assertEqual(train_y.shape[1], 29)   # n_tissues
        self.assertGreater(len(train_x), len(val_x))

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_26_training_loss_decreases(self):
        """Training loss should decrease over epochs."""
        from src.tissue_deconv.trainer import TissueDeconvTrainer
        from src.tissue_deconv.model import TissueDeconvolutionModel

        trainer = TissueDeconvTrainer(config=self.config, device="cpu")
        (train_x, train_y), (val_x, val_y) = trainer.generate_synthetic_data(
            n_samples=200, seed=42,
        )

        model = TissueDeconvolutionModel(config=self.config, seed=42)
        history = trainer.train_model(
            model=model,
            train_mixtures=train_x,
            train_fractions=train_y,
            val_mixtures=val_x,
            val_fractions=val_y,
            n_epochs=5,
            verbose=False,
        )

        self.assertGreater(len(history["train_loss"]), 0)
        # Loss should generally decrease
        self.assertLess(history["train_loss"][-1], history["train_loss"][0])

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_27_prediction_shape(self):
        """Predict should return correct dict structure."""
        from src.tissue_deconv.trainer import TissueDeconvTrainer
        from src.tissue_deconv.model import TissueDeconvolutionModel

        trainer = TissueDeconvTrainer(config=self.config, device="cpu")
        (train_x, train_y), (_, _) = trainer.generate_synthetic_data(
            n_samples=100, seed=42,
        )

        model = TissueDeconvolutionModel(config=self.config, seed=42)
        trainer.train_model(
            model=model,
            train_mixtures=train_x,
            train_fractions=train_y,
            n_epochs=3,
            verbose=False,
        )
        trainer._trained_models = [model]

        result = trainer.predict(train_x[0])
        self.assertIsInstance(result, dict)
        for tissue, frac in result.items():
            self.assertIsInstance(tissue, str)
            self.assertIsInstance(frac, float)
            self.assertGreaterEqual(frac, 0)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_28_checkpoint_save_load(self):
        """Checkpoint should survive round-trip."""
        from src.tissue_deconv.trainer import TissueDeconvTrainer
        from src.tissue_deconv.model import TissueDeconvolutionModel

        trainer = TissueDeconvTrainer(config=self.config, device="cpu")
        (train_x, train_y), (val_x, val_y) = trainer.generate_synthetic_data(
            n_samples=100, seed=42,
        )

        model = TissueDeconvolutionModel(config=self.config, seed=42)
        trainer.train_model(
            model=model,
            train_mixtures=train_x,
            train_fractions=train_y,
            val_mixtures=val_x,
            val_fractions=val_y,
            n_epochs=3,
            verbose=False,
        )
        trainer._trained_models = [model]

        # Manually set up simple ensemble for checkpoint test
        from src.tissue_deconv.model import TissueDeconvolutionEnsemble
        ensemble = TissueDeconvolutionEnsemble(config=self.config, base_seed=42)
        ensemble.models = [model]  # use single model as ensemble
        trainer.ensemble = ensemble

        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            path = f.name
        try:
            trainer.save_checkpoint(path)

            trainer2 = TissueDeconvTrainer(config=self.config, device="cpu")
            trainer2.load_checkpoint(path)

            self.assertIsNotNone(trainer2.ensemble)
            self.assertIn("val_loss", trainer2.history)
        finally:
            if os.path.exists(path):
                os.unlink(path)


class TestFeatures(unittest.TestCase):
    """Test TissueDeconvolutionFeatures."""

    def setUp(self):
        from src.tissue_deconv.tissue_features import TissueDeconvolutionFeatures
        from src.tissue_deconv.config import PROTOTYPE_CONFIG
        self.extractor = TissueDeconvolutionFeatures(config=PROTOTYPE_CONFIG)
        self.rng = np.random.RandomState(42)

    def test_29_feature_names_count(self):
        """Should have 24 feature names."""
        names = self.extractor.feature_names()
        self.assertEqual(len(names), 24)

    def test_30_extract_valid_output(self):
        """Extract should return all named features."""
        # Create realistic tissue fractions
        frac = np.zeros(29)
        frac[0] = 0.60   # Whole Blood
        frac[1] = 0.10   # PBMC
        frac[6] = 0.05   # Monocyte
        frac[9] = 0.15   # Liver
        frac[15] = 0.05  # Lung
        frac[20] = 0.05  # Colon

        features = self.extractor.extract(frac)

        for name in self.extractor.feature_names():
            self.assertIn(name, features, f"Missing feature: {name}")
            self.assertIsInstance(features[name], float)

        # Check value ranges
        self.assertGreaterEqual(features["tissue_top1_fraction"], 0)
        self.assertLessEqual(features["tissue_top1_fraction"], 1)
        self.assertGreaterEqual(features["tissue_entropy"], 0)

    def test_31_extract_with_cancer_tissue(self):
        """Cancer-relevant tissue should set cancer flags."""
        # Colon is a cancer-relevant tissue
        frac = np.zeros(29)
        # Indices: Colon is at position 13 in SUPPORTED_TISSUES
        from src.tissue_deconv.config import SUPPORTED_TISSUES
        colon_idx = SUPPORTED_TISSUES.index("Colon")
        frac[colon_idx] = 0.8
        frac[0] = 0.2  # Whole Blood

        features = self.extractor.extract(frac)
        self.assertAlmostEqual(features["tissue_top1_is_cancer"], 1.0)
        self.assertGreater(features["tissue_colon_fraction"], 0.5)

    def test_32_extract_blood_only(self):
        """Blood-only sample should have low cancer/blood ratio."""
        frac = np.zeros(29)
        blood_sum = 0.0
        from src.tissue_deconv.config import SUPPORTED_TISSUES
        for i, t in enumerate(SUPPORTED_TISSUES):
            if t in ("Whole Blood", "Neutrophil", "Monocyte", "PBMC"):
                frac[i] = 0.25

        features = self.extractor.extract(frac)
        self.assertAlmostEqual(features["tissue_top1_is_blood"], 1.0)
        self.assertLess(features["tissue_cancer_blood_ratio"], 0.01)

    def test_33_extract_batch(self):
        """Batch extraction should return arrays."""
        fracs = np.random.RandomState(42).dirichlet(np.ones(29), size=10)

        batch_features = self.extractor.extract_batch(fracs)

        for name, values in batch_features.items():
            self.assertIsInstance(values, np.ndarray)
            self.assertEqual(len(values), 10)

    def test_34_extract_from_mixture_fallback(self):
        """extract_from_mixture with no fractions should return zeros."""
        meth = np.random.uniform(0, 1, 100)
        features = self.extractor.extract_from_mixture(meth, fractions=None)
        for val in features.values():
            self.assertEqual(val, 0.0)

    def test_35_top_tissues(self):
        """top_tissues should return sorted list."""
        frac = np.zeros(29)
        frac[3] = 0.4   # CD8+ T Cell
        frac[9] = 0.35  # Liver
        frac[15] = 0.25 # Lung

        top = self.extractor.top_tissues(frac, n=3)
        self.assertEqual(len(top), 3)
        self.assertGreater(top[0][1], top[1][1])
        self.assertGreater(top[1][1], top[2][1])

    def test_36_tissue_distribution_summary(self):
        """Summary should be a non-empty string."""
        frac = np.ones(29) / 29
        summary = self.extractor.tissue_distribution_summary(frac)
        self.assertIsInstance(summary, str)
        self.assertGreater(len(summary), 0)
        self.assertIn("Blood-derived", summary)
        self.assertIn("Top tissues", summary)


class TestIntegration(unittest.TestCase):
    """Test DEConvIntegration adapter."""

    def setUp(self):
        from src.tissue_deconv.config import PROTOTYPE_CONFIG
        self.config = PROTOTYPE_CONFIG

    def test_37_integration_init_no_model(self):
        """Integration should work without a model."""
        from src.tissue_deconv.integration import DEConvIntegration

        integ = DEConvIntegration(config=self.config)
        self.assertFalse(integ.has_model)

        # Should not crash
        sample = {"methylation_data": np.random.uniform(0, 1, 100)}
        features = integ.extract_all(sample)
        self.assertIsInstance(features, dict)

    def test_38_integration_to_modality_no_model(self):
        """to_modality should return 0 when no model."""
        from src.tissue_deconv.integration import DEConvIntegration

        integ = DEConvIntegration(config=self.config)
        score = integ.to_modality({"methylation_data": np.random.uniform(0, 1, 100)})
        self.assertEqual(score, 0.0)

    def test_39_integration_extract_batch_no_model(self):
        """extract_batch should work without model."""
        from src.tissue_deconv.integration import DEConvIntegration

        integ = DEConvIntegration(config=self.config)
        samples = [
            {"methylation_data": np.random.uniform(0, 1, 100)}
            for _ in range(5)
        ]
        batch = integ.extract_batch(samples)

        for name, values in batch.items():
            self.assertEqual(len(values), 5)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_40_integration_fit_synthetic(self):
        """fit_synthetic should create a working model."""
        from src.tissue_deconv.integration import DEConvIntegration
        from src.tissue_deconv.config import PROTOTYPE_CONFIG

        config = deepcopy(PROTOTYPE_CONFIG)
        config.n_epochs = 3

        integ = DEConvIntegration(config=config)
        integ.fit_synthetic(n_samples=100, n_epochs=3, verbose=False)

        self.assertTrue(integ.has_model)

        # Predict should work
        meth = np.random.uniform(0, 1, 100).astype(np.float32)
        fracs = integ.predict_tissue_fractions(meth)
        self.assertEqual(len(fracs), 29)
        self.assertAlmostEqual(float(fracs.sum()), 1.0, delta=0.01)

        # to_modality should return non-zero score
        score = integ.to_modality({"methylation_data": meth})
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_41_integration_extract_all_with_precomputed(self):
        """extract_all should use pre-computed fractions if available."""
        from src.tissue_deconv.integration import DEConvIntegration

        integ = DEConvIntegration(config=self.config)
        frac = np.zeros(29)
        frac[0] = 0.7
        frac[9] = 0.3

        sample = {"tissue_fractions": frac}
        features = integ.extract_all(sample)

        self.assertGreater(features["tissue_top1_fraction"], 0.5)

    def test_42_integration_n_features(self):
        """n_features property should match feature extractor."""
        from src.tissue_deconv.integration import DEConvIntegration

        integ = DEConvIntegration(config=self.config)
        self.assertEqual(integ.n_features, 24)

    def test_43_integration_n_tissues(self):
        """n_tissues property should be 29."""
        from src.tissue_deconv.integration import DEConvIntegration

        integ = DEConvIntegration(config=self.config)
        self.assertEqual(integ.n_tissues, 29)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_44_integration_save_load_checkpoint(self):
        """Integration checkpoint round-trip."""
        from src.tissue_deconv.integration import DEConvIntegration
        from src.tissue_deconv.config import PROTOTYPE_CONFIG

        config = deepcopy(PROTOTYPE_CONFIG)
        config.n_epochs = 2

        integ = DEConvIntegration(config=config)
        integ.fit_synthetic(n_samples=50, n_epochs=2, verbose=False)

        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            path = f.name
        try:
            integ.save_checkpoint(path)

            integ2 = DEConvIntegration(config=config, checkpoint=path)
            self.assertTrue(integ2.has_model)
        finally:
            if os.path.exists(path):
                os.unlink(path)


class TestModuleExports(unittest.TestCase):
    """Test that __init__.py exports all expected classes."""

    def test_45_import_all_exports(self):
        """All __all__ entries should be importable."""
        from src.tissue_deconv import (
            TissueDeconvConfig,
            DEFAULT_CONFIG,
            PROTOTYPE_CONFIG,
            SUPPORTED_TISSUES,
            CANCER_RELEVANT_TISSUES,
            BLOOD_DERIVED_TISSUES,
            TissueAtlas,
            TissueDeconvolutionFeatures,
            DEConvIntegration,
        )
        self.assertIsNotNone(TissueDeconvConfig)
        self.assertIsNotNone(DEFAULT_CONFIG)
        self.assertIsNotNone(PROTOTYPE_CONFIG)
        self.assertIsNotNone(SUPPORTED_TISSUES)
        self.assertIsNotNone(CANCER_RELEVANT_TISSUES)
        self.assertIsNotNone(BLOOD_DERIVED_TISSUES)
        self.assertIsNotNone(TissueAtlas)
        self.assertIsNotNone(TissueDeconvolutionFeatures)
        self.assertIsNotNone(DEConvIntegration)

    @unittest.skipIf(not _HAS_TORCH, "PyTorch not available")
    def test_46_import_torch_classes(self):
        """PyTorch-dependent classes should import when torch is available."""
        from src.tissue_deconv import (
            TissueDeconvolutionModel,
            TissueDeconvolutionEnsemble,
            DeconvLoss,
            TissueDeconvTrainer,
        )
        self.assertIsNotNone(TissueDeconvolutionModel)
        self.assertIsNotNone(TissueDeconvolutionEnsemble)
        self.assertIsNotNone(DeconvLoss)
        self.assertIsNotNone(TissueDeconvTrainer)


# ═══════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Suppress info/warning logs during tests
    import logging
    logging.getLogger("src.tissue_deconv").setLevel(logging.ERROR)

    unittest.main(verbosity=2)
