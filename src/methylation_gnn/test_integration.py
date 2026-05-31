#!/usr/bin/env python3
"""
Integration Test — GNN Methylation Network in DeepCatch CET Pipeline
======================================================================

Demonstrates how the new GNN methylation branch connects with the existing
DeepCatch Stage 1 (Capture) pipeline for multi-modal cancer detection.

This script tests:
1. Module importability with graceful fallbacks
2. Config validation
3. Graph construction (with synthetic data)
4. Model architecture (forward pass with dummy data)
5. Training loop (Phase 1 pretraining + Phase 2 finetuning)
6. Inference (score extraction)
7. Integration with fusion layer (ModularArmsBuilder + CrossAttentionFusion)
8. Co-methylation matrix building

Usage
-----

    # Run all tests (mocked, no real data needed)
    python src/methylation_gnn/test_integration.py

Environment
-----------
Requires: numpy, scipy, scikit-learn
Optional: torch (for model tests), torch_geometric (for graph tests)
"""

from __future__ import annotations

import os
import sys
import unittest
import warnings
from typing import Any, Dict, List, Optional
from copy import deepcopy

import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# ── PyTorch availability probes ─────────────────────────────────

_HAS_TORCH = False
_HAS_PYG = False
try:
    import torch
    _HAS_TORCH = True
except ImportError:
    pass
try:
    import torch_geometric
    from torch_geometric.data import Data
    _HAS_PYG = True
except ImportError:
    pass


# ── Synthetic data helpers ──────────────────────────────────────

def make_synthetic_methylation_data(
    n_regions: int = 200,
    n_samples: int = 20,
    n_cancer: int = 10,
    seed: int = 42,
) -> Dict[str, Any]:
    """Create small synthetic dataset for integration testing."""
    rng = np.random.RandomState(seed)
    n_healthy = n_samples - n_cancer

    base_state = rng.choice([0.2, 0.8], size=n_regions, p=[0.3, 0.7])

    healthy_beta = np.zeros((n_regions, n_healthy), dtype=np.float32)
    for i in range(n_healthy):
        noise = rng.normal(0, 0.05, size=n_regions)
        healthy_beta[:, i] = np.clip(base_state + noise, 0.0, 1.0)

    hyper_mask = rng.rand(n_regions) < 0.3
    cancer_beta = np.zeros((n_regions, n_cancer), dtype=np.float32)
    for i in range(n_cancer):
        noise = rng.normal(0, 0.08, size=n_regions)
        state = base_state.copy()
        state[~hyper_mask] -= 0.15
        state[hyper_mask] += 0.1
        cancer_beta[:, i] = np.clip(state + noise, 0.0, 1.0)

    beta_values = np.hstack([healthy_beta, cancer_beta]).astype(np.float32)
    labels = np.array([0] * n_healthy + [1] * n_cancer, dtype=np.int64)

    return {
        "beta_values": beta_values,
        "labels": labels,
        "cpg_density": rng.uniform(0, 30, size=n_regions).astype(np.float32),
        "gc_content": rng.uniform(0.3, 0.7, size=n_regions).astype(np.float32),
        "n_regions": n_regions,
        "n_samples": n_samples,
    }


# ── Test Cases ──────────────────────────────────────────────────

class TestConfig(unittest.TestCase):
    """Test GNNConfig validation and defaults."""

    def test_default_config(self):
        from src.methylation_gnn.config import GNNConfig, DEFAULT_GNN_CONFIG
        cfg = DEFAULT_GNN_CONFIG
        self.assertEqual(cfg.n_nodes, 50_000)
        self.assertEqual(cfg.n_node_features, 20)
        self.assertEqual(cfg.hidden_dims, [64, 128, 256])
        self.assertEqual(cfg.n_layers, 3)
        self.assertEqual(cfg.final_dim, 256)
        self.assertGreater(cfg.estimated_params, 0)
        self.assertIn(cfg.device, ["cuda", "mps", "cpu"])

    def test_prototype_config(self):
        from src.methylation_gnn.config import PROTOTYPE_GNN_CONFIG
        cfg = PROTOTYPE_GNN_CONFIG
        self.assertEqual(cfg.n_nodes, 5_000)
        self.assertEqual(cfg.n_epochs_pretrain, 20)

    def test_to_dict(self):
        from src.methylation_gnn.config import GNNConfig
        cfg = GNNConfig(n_nodes=1000)
        d = cfg.to_dict()
        self.assertEqual(d["n_nodes"], 1000)
        self.assertIn("hidden_dims", d)

    def test_invalid_config(self):
        from src.methylation_gnn.config import GNNConfig
        with self.assertRaises(ValueError):
            GNNConfig(hidden_dims=[])  # type: ignore
        with self.assertRaises(ValueError):
            GNNConfig(n_node_features=0)

    def test_node_feature_spec(self):
        from src.methylation_gnn.config import NODE_FEATURE_SPEC
        self.assertIn("mean_methylation", NODE_FEATURE_SPEC)
        self.assertIn("cpg_density", NODE_FEATURE_SPEC)
        self.assertEqual(len(NODE_FEATURE_SPEC), 20)


class TestGraphBuilder(unittest.TestCase):
    """Test graph construction (works without PyTorch)."""

    def test_no_regions_fallback(self):
        """Builder should create synthetic regions if none loaded."""
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        builder = RegulatoryGraphBuilder(n_nodes=50, edge_k=5)
        # No regions loaded → should auto-generate
        builder._ensure_regions()
        self.assertEqual(len(builder._regions), 50)

    def test_load_regions_from_list(self):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder, GenomicRegion
        regions = [
            GenomicRegion(chrom="chr1", start=i * 1000, end=i * 1000 + 500, region_type="cpg_island")
            for i in range(100)
        ]
        builder = RegulatoryGraphBuilder(n_nodes=100)
        builder.load_regions_from_list(regions)
        self.assertEqual(len(builder._regions), 100)

    def test_build_node_features(self):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        builder = RegulatoryGraphBuilder(n_nodes=20)
        builder._ensure_regions()
        meth_data = {
            "beta_values": np.random.rand(20).astype(np.float32),
            "cpg_density": np.random.rand(20).astype(np.float32) * 30,
            "gc_content": np.random.rand(20).astype(np.float32),
        }
        x = builder.build_node_features(meth_data)
        self.assertEqual(x.shape, (20, 20))
        self.assertEqual(x.dtype, np.float32)

    def test_build_node_features_with_chromatin(self):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        builder = RegulatoryGraphBuilder(n_nodes=10)
        builder._ensure_regions()
        meth_data = {"beta_values": np.random.rand(10).astype(np.float32)}
        chrom_data = {
            "dnase": np.random.rand(10).astype(np.float32),
            "h3k4me3": np.random.rand(10).astype(np.float32),
        }
        x = builder.build_node_features(meth_data, chromatin_data=chrom_data)
        self.assertEqual(x.shape, (10, 20))

    def test_build_node_features_missing_data(self):
        """Features should gracefully handle completely missing data."""
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        builder = RegulatoryGraphBuilder(n_nodes=5)
        builder._ensure_regions()
        x = builder.build_node_features({})
        self.assertEqual(x.shape, (5, 20))
        # All-zero except region type one-hot
        self.assertTrue(np.all(x[:, :15] == 0))

    def test_build_edges_fallback(self):
        """Edges should not crash when all sources are empty."""
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        builder = RegulatoryGraphBuilder(n_nodes=30, edge_k=5)
        builder._ensure_regions()
        ei, ea, et = builder.build_edges()
        self.assertGreater(ei.shape[1], 0)  # should have fallback edges
        self.assertEqual(ei.shape[0], 2)

    def test_genomic_region(self):
        from src.methylation_gnn.graph_builder import GenomicRegion
        r1 = GenomicRegion("chr1", 1000, 2000, region_type="cpg_island")
        r2 = GenomicRegion("chr1", 4000, 5000, region_type="enhancer")
        r3 = GenomicRegion("chr2", 1000, 2000, region_type="promoter")
        
        self.assertEqual(r1.length, 1000)
        self.assertEqual(r1.midpoint, 1500)
        self.assertEqual(r1.type_id, 0)
        self.assertEqual(r2.type_id, 1)
        self.assertLess(r1.distance_to(r2), 1_000_000)  # same chr, finite
        self.assertGreater(r1.distance_to(r3), 1_00_000)  # diff chr, large

    @unittest.skipIf(not _HAS_PYG, "torch_geometric not installed")
    def test_build_graph_pyg(self):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        builder = RegulatoryGraphBuilder(n_nodes=50, edge_k=5)
        builder._ensure_regions()
        meth_data = {
            "beta_values": np.random.rand(50).astype(np.float32),
        }
        cov = np.random.poisson(30, size=50).astype(float)
        graph = builder.build_graph(
            sample_name="test",
            methylation_data=meth_data,
            cfDNA_coverage=cov,
            label=1,
            cancer_type="TEST",
        )
        self.assertEqual(graph.x.shape, (50, 20))
        self.assertGreater(graph.edge_index.shape[1], 0)
        self.assertEqual(graph.y.item(), 1)
        self.assertEqual(graph.sample_name, "test")


@unittest.skipIf(not _HAS_PYG, "torch_geometric not installed")
class TestGNNModel(unittest.TestCase):
    """Test GNN model architecture (requires PyTorch, PyG)."""

    @classmethod
    def setUpClass(cls):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        data = make_synthetic_methylation_data(n_regions=50, n_samples=10, n_cancer=5)
        builder = RegulatoryGraphBuilder(n_nodes=50, edge_k=5)
        builder._ensure_regions()
        cls.test_graph = builder.build_graph(
            sample_name="test",
            methylation_data={k: data[k] for k in ["beta_values"] if k in data},
            cfDNA_coverage=np.random.poisson(30, size=50).astype(float),
            label=0,
        )
        cls.test_graphs = [
            builder.build_graph(
                sample_name=f"s{i}",
                methylation_data={
                    "beta_values": data["beta_values"][:, i].astype(np.float32),
                },
                cfDNA_coverage=np.random.poisson(30, size=50).astype(float),
                label=int(data["labels"][i]),
            )
            for i in range(10)
        ]

    def test_model_creation(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        model = MethylationGNN(
            n_node_features=20,
            hidden_dims=[32, 64, 128],
            n_edge_types=5,
            n_attention_heads=2,
        )
        self.assertGreater(model.num_parameters, 0)
        self.assertLess(model.num_parameters, 5_000_000)  # reasonable size

    def test_forward_pass(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        model = MethylationGNN(
            n_node_features=20,
            hidden_dims=[32, 64, 128],
            n_edge_types=5,
            n_attention_heads=2,
        )
        graph = self.test_graph
        output = model(graph.x, graph.edge_index, graph.edge_type)
        self.assertIn("reconstructed", output)
        self.assertIn("anomaly_scores", output)
        self.assertIn("graph_anomaly", output)
        self.assertEqual(output["reconstructed"].shape, graph.x.shape)
        self.assertEqual(output["anomaly_scores"].shape[0], graph.x.shape[0])

    def test_pretrain_forward(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        model = MethylationGNN(
            n_node_features=20,
            hidden_dims=[32, 64, 128],
            n_edge_types=5,
            n_attention_heads=2,
        )
        graph = self.test_graph
        recon, mask, x_orig, anom = model.pretrain_forward(
            graph.x, graph.edge_index, graph.edge_type, mask_ratio=0.3
        )
        self.assertEqual(recon.shape, graph.x.shape)
        self.assertGreater(mask.sum(), 0)
        self.assertLess(mask.sum(), len(mask))

    def test_predict_method(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        model = MethylationGNN(
            n_node_features=20,
            hidden_dims=[32, 64, 128],
            n_edge_types=5,
            n_attention_heads=2,
        )
        score = model.predict(self.test_graph)
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    def test_homogeneous_fallback(self):
        """Model should work with single edge type."""
        from src.methylation_gnn.gnn_model import MethylationGNN
        model = MethylationGNN(
            n_node_features=20,
            hidden_dims=[32, 64],
            n_edge_types=1,
            n_attention_heads=2,
        )
        graph = self.test_graph
        # All edges as type 0
        et = torch.zeros(graph.edge_index.shape[1], dtype=torch.long)
        output = model(graph.x, graph.edge_index, et)
        self.assertIn("graph_anomaly", output)

    def test_field_defect_loss(self):
        from src.methylation_gnn.gnn_model import FieldDefectLoss, MethylationGNN
        model = MethylationGNN(
            n_node_features=20,
            hidden_dims=[32, 64],
            n_edge_types=5,
            n_attention_heads=2,
        )
        graph = self.test_graph
        output = model(graph.x, graph.edge_index, graph.edge_type)
        
        loss_fn = FieldDefectLoss(lambda_anomaly=0.5)
        loss_total, metrics = loss_fn(
            x_original=graph.x,
            x_reconstructed=output["reconstructed"],
            anomaly_scores=output["anomaly_scores"],
            labels=graph.y,
        )
        self.assertGreater(loss_total.item(), 0)
        self.assertIn("loss_recon", metrics)
        self.assertIn("loss_anomaly", metrics)

    def test_model_small(self):
        """Smallest possible model (~100K params) for embedded/IoT."""
        from src.methylation_gnn.gnn_model import MethylationGNN
        model = MethylationGNN(
            n_node_features=10,
            hidden_dims=[16, 32],
            n_edge_types=1,
            n_attention_heads=1,
            decoder_hidden=[16],
            anomaly_hidden=[16],
        )
        self.assertGreater(model.num_parameters, 0)
        # Should handle different n_features
        x = torch.randn(20, 10)
        ei = torch.randint(0, 20, (2, 50))
        output = model(x, ei)
        self.assertEqual(output["reconstructed"].shape, (20, 10))


@unittest.skipIf(not _HAS_PYG, "torch_geometric not installed")
class TestGNNTrainer(unittest.TestCase):
    """Test training pipeline (requires PyTorch, PyG)."""

    @classmethod
    def setUpClass(cls):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        from src.methylation_gnn.config import PROTOTYPE_GNN_CONFIG
        cls.config = PROTOTYPE_GNN_CONFIG
        cls.config.n_epochs_pretrain = 3
        cls.config.n_epochs_finetune = 3
        cls.config.checkpoint_dir = "/tmp/gnn_test_checkpoints"
        os.makedirs(cls.config.checkpoint_dir, exist_ok=True)
        
        data = make_synthetic_methylation_data(n_regions=100, n_samples=20, n_cancer=10)
        builder = RegulatoryGraphBuilder(n_nodes=100, edge_k=5)
        builder._ensure_regions()
        cls.graphs = [
            builder.build_graph(
                sample_name=f"s{i}",
                methylation_data={
                    "beta_values": data["beta_values"][:, i].astype(np.float32),
                    "cpg_density": data["cpg_density"],
                    "gc_content": data["gc_content"],
                },
                cfDNA_coverage=np.random.poisson(30, size=100).astype(float),
                label=int(data["labels"][i]),
            )
            for i in range(20)
        ]
        # Split
        cls.train_graphs = cls.graphs[:12]
        cls.val_graphs = cls.graphs[12:16]
        cls.test_graphs = cls.graphs[16:]

    def test_trainer_creation(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, self.config)
        self.assertEqual(trainer.phase.value, "uninitialized")

    def test_pretrain_epoch(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, self.config)
        loss_val = trainer._pretrain_epoch(self.train_graphs[:5], mask_ratio=0.3)
        self.assertGreater(loss_val["loss"], 0)

    def test_full_training_cycle(self):
        """End-to-end pretrain + finetune cycle."""
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, self.config)
        
        # Phase 1
        result = trainer.pretrain(
            self.train_graphs, self.val_graphs, n_epochs=2
        )
        self.assertEqual(result["phase"], "pretrain")
        
        # Phase 2
        result = trainer.finetune(
            self.train_graphs, self.val_graphs, n_epochs=2
        )
        self.assertEqual(result["phase"], "finetune")
        
        # Evaluate
        metrics = trainer.evaluate(self.test_graphs)
        self.assertIn("auc", metrics)
        self.assertIn("loss_recon", metrics)

    def test_checkpoint_save_load(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, self.config)
        
        trainer._save_checkpoint("test_checkpoint.pt")
        path = os.path.join(self.config.checkpoint_dir, "test_checkpoint.pt")
        self.assertTrue(os.path.exists(path))
        
        # Load into new trainer
        model2 = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer2 = GNNTrainer(model2, self.config)
        epoch = trainer2.load_checkpoint(path)
        self.assertGreaterEqual(epoch, 0)
        
    def test_predict(self):
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, self.config)
        
        score = trainer.predict(self.test_graphs[0])
        self.assertIsInstance(score, float)
        
        details = trainer.predict(self.test_graphs[0], return_details=True)
        self.assertIsInstance(details, dict)
        self.assertIn("field_defect_score", details)
        self.assertIn("node_anomaly_scores", details)


@unittest.skipIf(not _HAS_PYG, "torch_geometric not installed")
class TestGNNInference(unittest.TestCase):
    """Test inference pipeline (requires PyTorch, PyG)."""

    @classmethod
    def setUpClass(cls):
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        from src.methylation_gnn.config import PROTOTYPE_GNN_CONFIG
        
        cls.ckpt_dir = "/tmp/gnn_test_checkpoints"
        cls.ckpt_path = os.path.join(cls.ckpt_dir, "test_inference.pt")
        os.makedirs(cls.ckpt_dir, exist_ok=True)
        
        cfg = PROTOTYPE_GNN_CONFIG
        cfg.n_epochs_pretrain = 2
        cfg.n_epochs_finetune = 2
        cfg.checkpoint_dir = cls.ckpt_dir
        
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, cfg)
        
        data = make_synthetic_methylation_data(n_regions=50, n_samples=10, n_cancer=5)
        builder = RegulatoryGraphBuilder(n_nodes=50, edge_k=5)
        builder._ensure_regions()
        graphs = [
            builder.build_graph(
                sample_name=f"s{i}",
                methylation_data={
                    "beta_values": data["beta_values"][:, i].astype(np.float32),
                },
                label=int(data["labels"][i]),
            )
            for i in range(10)
        ]
        trainer.pretrain(graphs[:6], graphs[6:8], n_epochs=2)
        trainer.finetune(graphs[:6], graphs[6:8], n_epochs=2)
        trainer._save_checkpoint("test_inference.pt")
        
        cls.test_graph = graphs[0]
        cls.builder = builder

    def test_inference_load(self):
        from src.methylation_gnn.gnn_inference import GNNInference
        inf = GNNInference.load(self.ckpt_path)
        self.assertIsNotNone(inf.model)

    def test_inference_predict(self):
        from src.methylation_gnn.gnn_inference import GNNInference
        inf = GNNInference.load(self.ckpt_path)
        score = inf.predict(self.test_graph)
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    def test_inference_predict_details(self):
        from src.methylation_gnn.gnn_inference import GNNInference
        inf = GNNInference.load(self.ckpt_path)
        details = inf.predict(self.test_graph, return_details=True)
        self.assertIn("field_defect_score", details)
        self.assertIn("node_anomaly_scores", details)
        self.assertIn("top_anomalous_nodes", details)
        self.assertIsInstance(details["top_anomalous_nodes"], list)

    def test_predictor_high_level(self):
        from src.methylation_gnn.gnn_inference import MethylationGNNPredictor
        predictor = MethylationGNNPredictor(
            builder=self.builder,
            checkpoint=self.ckpt_path,
        )
        score = predictor.predict_sample(
            sample_name="test_sample",
            methylation_data={
                "beta_values": np.random.rand(50).astype(np.float32),
            },
        )
        self.assertIsInstance(score, float)


class TestIntegration(unittest.TestCase):
    """Test integration adapter (works without PyTorch for pure data flow)."""

    def test_extend_fusion_with_gnn(self):
        from src.methylation_gnn.integration import extend_fusion_with_gnn
        existing = [
            np.random.rand(20).astype(np.float32),
            np.random.rand(20).astype(np.float32),
            np.random.rand(20).astype(np.float32),
            np.random.rand(20).astype(np.float32),
        ]
        gnn = np.random.rand(20).astype(np.float32)
        extended = extend_fusion_with_gnn(existing, gnn)
        self.assertEqual(len(extended), 5)
        self.assertEqual(len(extended[0]), 20)

    def test_modular_arms_builder_no_gnn(self):
        """Arms builder without GNN (v2.0 compatibility mode)."""
        from src.methylation_gnn.integration import ModularArmsBuilder
        arms = ModularArmsBuilder(include_gnn=False)
        self.assertEqual(arms.n_modalities, 4)
        
        sample = {
            "fsi": 2.5,
            "caff_score": 0.3,
            "serological_score": 0.6,
            "mfr_score": 0.45,
        }
        result = arms.process_sample(sample)
        self.assertIn("fragmentomics_score", result)
        self.assertEqual(result["n_modalities"], 4)

    def test_modular_arms_builder_extract_all(self):
        from src.methylation_gnn.integration import ModularArmsBuilder
        arms = ModularArmsBuilder(include_gnn=False)
        samples = [
            {"fsi": 1.5 + i * 0.5, "caff_score": 0.2, "serological_score": 0.5}
            for i in range(10)
        ]
        scores = arms.extract_all_scores(samples)
        self.assertEqual(len(scores), 4)
        for s in scores:
            self.assertEqual(len(s), 10)
            self.assertTrue(np.all(np.isfinite(s)))

    def test_extend_then_fit_fusion(self):
        """Demonstrate: get scores → extend with GNN → fit fusion."""
        from src.methylation_gnn.integration import extend_fusion_with_gnn
        try:
            from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion
        except ImportError:
            self.skipTest("multimodal_fusion not importable")

        n = 50
        rng = np.random.RandomState(42)
        frag = rng.rand(n)
        cnv = rng.rand(n)
        sero = rng.rand(n)
        mfr = rng.rand(n)
        gnn = rng.rand(n)
        labels = (rng.rand(n) > 0.5).astype(np.int64)

        # v2.0: 4 modalities
        fusion4 = CrossAttentionFusion(n_modalities=4)
        fusion4.fit([frag, cnv, sero, mfr], labels)
        probs4 = fusion4.predict_proba([frag, cnv, sero, mfr])
        self.assertEqual(len(probs4), n)

        # v2.1: 5 modalities with GNN
        all_scores = extend_fusion_with_gnn([frag, cnv, sero, mfr], gnn)
        fusion5 = CrossAttentionFusion(n_modalities=5)
        fusion5.fit(all_scores, labels)
        probs5 = fusion5.predict_proba(all_scores)
        self.assertEqual(len(probs5), n)

    @unittest.skipIf(not _HAS_PYG, "torch_geometric not installed")
    def test_branch_adapter(self):
        from src.methylation_gnn.integration import MethylationBranchAdapter
        from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
        from src.methylation_gnn.gnn_model import MethylationGNN
        from src.methylation_gnn.gnn_trainer import GNNTrainer
        from src.methylation_gnn.config import PROTOTYPE_GNN_CONFIG
        
        # Train a quick model for testing
        cfg = PROTOTYPE_GNN_CONFIG
        cfg.n_epochs_pretrain = 2
        cfg.n_epochs_finetune = 2
        cfg.checkpoint_dir = "/tmp/gnn_test_checkpoints"
        ckpt = os.path.join(cfg.checkpoint_dir, "adapter_test.pt")
        
        model = MethylationGNN(n_node_features=20, hidden_dims=[32, 64, 128], n_attention_heads=2)
        trainer = GNNTrainer(model, cfg)
        builder = RegulatoryGraphBuilder(n_nodes=50, edge_k=5)
        builder._ensure_regions()
        
        data = make_synthetic_methylation_data(n_regions=50, n_samples=10, n_cancer=5)
        graphs = [
            builder.build_graph(
                sample_name=f"s{i}",
                methylation_data={"beta_values": data["beta_values"][:, i].astype(np.float32)},
                label=int(data["labels"][i]),
            )
            for i in range(10)
        ]
        trainer.pretrain(graphs[:6], graphs[6:8], n_epochs=2)
        trainer._save_checkpoint("adapter_test.pt")
        
        adapter = MethylationBranchAdapter(ckpt, n_regions=50)
        score = adapter.predict_sample(
            {"beta_values": np.random.rand(50, 1).astype(np.float32)},
            sample_name="test",
        )
        self.assertIsInstance(score, float)
        
        # Batch
        samples = [
            {"methylation_data": {"beta_values": np.random.rand(50).astype(np.float32)}}
            for _ in range(5)
        ]
        scores = adapter.predict_batch(samples)
        self.assertEqual(scores.shape, (5,))


class TestDataUtils(unittest.TestCase):
    """Test reference data utilities (pure Python, no PyTorch needed)."""

    def test_reference_catalog(self):
        from src.methylation_gnn.data import ReferenceDataCatalog
        datasets = ReferenceDataCatalog.list_datasets()
        self.assertIn("ucsc_cpg_islands", datasets)
        self.assertIn("hic_rao2014", datasets)
        self.assertIn("tcga_methylation", datasets)

    def test_dataset_metadata(self):
        from src.methylation_gnn.data import ReferenceDataCatalog
        entry = ReferenceDataCatalog.get_entry("ucsc_cpg_islands")
        self.assertIsNotNone(entry)
        self.assertIn("CpG", entry.description)

    def test_download_urls(self):
        from src.methylation_gnn.data import (
            download_ucsc_cpg_islands,
            download_gencode_promoters,
            download_fantom5_enhancers,
        )
        url1 = download_ucsc_cpg_islands("hg38")
        self.assertTrue(url1.startswith("https://"))
        self.assertIn("cpgIslandExt", url1)
        
        url2 = download_gencode_promoters("44")
        self.assertTrue(url2.startswith("https://"))
        
        url3 = download_fantom5_enhancers()
        self.assertTrue(url3.startswith("https://"))

    def test_synthetic_data(self):
        from src.methylation_gnn.data import generate_synthetic_methylation_data
        data = generate_synthetic_methylation_data(
            n_regions=200, n_samples=30, n_cancer=15
        )
        self.assertEqual(data["beta_values"].shape, (200, 30))
        self.assertEqual(data["labels"].shape, (30,))
        self.assertEqual(np.sum(data["labels"] == 1), 15)
        self.assertEqual(data["cpg_density"].shape, (200,))
        self.assertEqual(data["gc_content"].shape, (200,))

    def test_preprocess_methylation_betas(self):
        from src.methylation_gnn.data import preprocess_methylation_betas
        betas = np.random.beta(2, 3, size=(100, 20)).astype(np.float32)
        betas[np.random.rand(100, 20) < 0.05] = np.nan  # some missing
        
        regions = [("chr1", i * 1000, i * 1000 + 500) for i in range(5)]
        positions = np.array([
            ("chr1", i * 1000 + 250) for i in range(100)
        ])
        
        result = preprocess_methylation_betas(betas, positions, regions)
        self.assertIn("beta_values", result)
        self.assertIn("methylation_variance", result)
        self.assertEqual(result["beta_values"].shape, (5, 20))

    def test_preprocess_methylation_no_regions(self):
        from src.methylation_gnn.data import preprocess_methylation_betas
        betas = np.random.beta(2, 3, size=(50, 10)).astype(np.float32)
        result = preprocess_methylation_betas(betas)
        self.assertEqual(result["beta_values"].shape, (50,))

    def test_build_comethylation_matrix(self):
        from src.methylation_gnn.data import build_comethylation_matrix
        betas = np.random.beta(2, 3, size=(30, 15)).astype(np.float32)
        corr = build_comethylation_matrix(betas, min_correlation=0.3)
        self.assertEqual(corr.shape, (30, 30))
        self.assertEqual(corr.dtype, np.float32)
        # Diagonal should be zero
        self.assertEqual(corr[0, 0], 0.0)

    def test_parse_regions_from_bed(self):
        import tempfile
        from src.methylation_gnn.data import parse_regions_from_bed
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write("chr1\t1000\t2000\tCpG_island\tGENE1\t+\n")
            f.write("chr1\t5000\t6000\tEnhancer\t.\t+\n")
            f.write("track name=test\n")
            f.write("chr2\t100\t500\tPromoter\tGENE2\t-\n")
            f.write("# comment line\n")
            bed_path = f.name
        
        try:
            regions = parse_regions_from_bed(bed_path, max_regions=10)
            self.assertEqual(len(regions), 3)
            self.assertEqual(regions[0][3], "cpg_island")
            self.assertEqual(regions[0][4], "GENE1")
            self.assertEqual(regions[1][4], None)
            self.assertEqual(regions[2][3], "promoter")
        finally:
            os.unlink(bed_path)

    def test_generate_cpg_regions(self):
        from src.methylation_gnn.data import generate_cpg_regions
        chrom_sizes = {"chr1": 100000, "chr2": 50000}
        regions = generate_cpg_regions(chrom_sizes, n_regions=10)
        self.assertGreaterEqual(len(regions), 9)  # may vary due to random sampling
        for r in regions:
            self.assertEqual(len(r), 5)  # (chrom, start, end, type, None)
            self.assertIn(r[0], ["chr1", "chr2"])

    def test_print_catalog(self):
        from src.methylation_gnn.data import ReferenceDataCatalog
        text = ReferenceDataCatalog.print_catalog()
        self.assertIn("UCSC", text)  # ucsc_cpg_islands present
        # TCGA not in print_catalog (different field)


class TestImports(unittest.TestCase):
    """Test that the module can be imported correctly."""

    def test_import_module(self):
        """from src.methylation_gnn import * should not error."""
        import src.methylation_gnn as gnn
        self.assertTrue(hasattr(gnn, "GNNConfig"))
        self.assertTrue(hasattr(gnn, "RegulatoryGraphBuilder"))
        self.assertTrue(hasattr(gnn, "MethylationGNN"))

    def test_import_all_names(self):
        import src.methylation_gnn as gnn
        for name in gnn.__all__:
            self.assertTrue(hasattr(gnn, name), f"Missing import: {name}")


# ── Main ────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("DeepCatch GNN Methylation Network — Integration Tests")
    print("=" * 65)
    print(f"PyTorch: {'✓' if _HAS_TORCH else '✗ (skip model/training tests)'}")
    print(f"PyTorch Geometric: {'✓' if _HAS_PYG else '✗ (skip graph/model tests)'}")
    print(f"numpy: {np.__version__}")
    if _HAS_TORCH:
        print(f"torch: {torch.__version__}")
    if _HAS_PYG:
        print(f"torch_geometric: {torch_geometric.__version__}")
    print("=" * 65)
    print()
    unittest.main(verbosity=2)
