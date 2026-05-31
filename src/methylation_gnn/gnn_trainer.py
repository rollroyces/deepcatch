#!/usr/bin/env python3
"""
GNN Trainer — Self-Supervised Pretraining + Supervised Fine-tuning
====================================================================

Training pipeline for the MethylationGNN with three phases:

**Phase 1: Self-Supervised Pretraining (Masked Node Prediction)**
    Train on healthy control graphs only. Randomly mask 30% of node features,
    task the model to reconstruct them from neighbor context. This teaches
    the model the "normal" methylation network topology and feature
    distribution.

**Phase 2: Joint Reconstruction + Anomaly Training**
    Train on both healthy and cancer graphs. Reconstruction loss keeps the
    network grounded in normal patterns; anomaly loss trains the anomaly
    head to distinguish cancer from healthy at the graph level.

**Phase 3 (optional): Fine-tuning with Labels**
    Train only the anomaly head on labeled data to maximize cross-cancer-type
    generalization.

Training supports:
- MPS (Apple Silicon), CUDA, and CPU backends
- Mixed precision (AMP) on CUDA for memory efficiency
- Early stopping with configurable patience
- Checkpointing with best-model tracking
- Detailed logging of per-epoch metrics

Example
-------

.. code-block:: python

    from src.methylation_gnn import MethylationGNN, GNNTrainer, GNNConfig

    config = GNNConfig()
    model = MethylationGNN(n_node_features=config.n_node_features)
    trainer = GNNTrainer(model, config)

    # Phase 1: Self-supervised pretraining on healthy controls
    trainer.pretrain(train_graphs_healthy, n_epochs=100)

    # Phase 2: Joint training on all data
    trainer.finetune(train_graphs, val_graphs, n_epochs=50)

    # Evaluate
    metrics = trainer.evaluate(test_graphs)
    print(f"Test AUC: {metrics['auc']:.3f}")
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    import torch
    import torch.nn.functional as F
    from torch.optim import AdamW
    from torch.optim.lr_scheduler import CosineAnnealingLR, ReduceLROnPlateau

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None
    AdamW = None

try:
    from torch_geometric.data import Data
    from torch_geometric.loader import DataLoader

    _HAS_PYG = True
except ImportError:
    _HAS_PYG = False
    Data = None
    DataLoader = None

try:
    from sklearn.metrics import roc_auc_score, average_precision_score

    _HAS_SKLEARN = True
except ImportError:
    _HAS_SKLEARN = False

from .config import GNNConfig
from .gnn_model import MethylationGNN, FieldDefectLoss, compute_field_defect_score

logger = logging.getLogger(__name__)


class GNNTrainerPhase(Enum):
    """Training phase enum for checkpoint metadata."""
    UNINITIALIZED = "uninitialized"
    PRETRAINING = "pretraining"  # Phase 1: masked node prediction
    JOINT_FINETUNING = "joint_finetuning"  # Phase 2: reconstruction + anomaly
    HEAD_ONLY = "head_only"  # Phase 3: anomaly head only


@dataclass
class TrainingHistory:
    """Tracks training metrics across epochs."""
    epoch: List[int] = field(default_factory=list)
    train_loss: List[float] = field(default_factory=list)
    train_loss_recon: List[float] = field(default_factory=list)
    train_loss_anomaly: List[float] = field(default_factory=list)
    val_loss: List[Optional[float]] = field(default_factory=list)
    val_auc: List[Optional[float]] = field(default_factory=list)
    val_auprc: List[Optional[float]] = field(default_factory=list)
    lr: List[float] = field(default_factory=list)

    def record(
        self,
        epoch: int,
        train_metrics: Dict[str, float],
        val_metrics: Optional[Dict[str, float]] = None,
        lr: float = 0.0,
    ):
        self.epoch.append(epoch)
        self.train_loss.append(train_metrics.get("loss_total", train_metrics.get("loss", 0)))
        self.train_loss_recon.append(train_metrics.get("loss_recon", 0))
        self.train_loss_anomaly.append(train_metrics.get("loss_anomaly", 0))
        self.val_loss.append(val_metrics.get("loss_total") if val_metrics else None)
        self.val_auc.append(val_metrics.get("auc") if val_metrics else None)
        self.val_auprc.append(val_metrics.get("auprc") if val_metrics else None)
        self.lr.append(lr)


class GNNTrainer:
    """
    Trainer for MethylationGNN with multi-phase training.

    Handles device management, optimization, scheduling, checkpointing,
    early stopping, and logging.

    Parameters
    ----------
    model : MethylationGNN
        The GNN model to train.
    config : GNNConfig
        Training configuration (hyperparameters, device, etc.).
    """

    def __init__(
        self,
        model: MethylationGNN,
        config: Optional[GNNConfig] = None,
    ):
        if not _HAS_PYG:
            raise ImportError(
                "PyTorch Geometric is required for GNNTrainer. "
                "Install with: pip install torch_geometric"
            )

        if config is None:
            config = GNNConfig()

        self.config = config
        self.device = torch.device(config.device)
        self.model = model.to(self.device)
        self.phase = GNNTrainerPhase.UNINITIALIZED
        self.history = TrainingHistory()
        self.best_val_loss = float("inf")
        self.best_val_auc = 0.0
        self.best_epoch = 0
        self.current_epoch = 0

        # Optimizer
        self.optimizer = AdamW(
            model.parameters(),
            lr=config.learning_rate,
            weight_decay=config.weight_decay,
        )

        # Scheduler (cosine annealing for pretraining, reduce-on-plateau for finetune)
        self.scheduler: Optional[Any] = None

        # Loss function
        self.loss_fn = FieldDefectLoss(
            lambda_anomaly=config.lambda_anomaly,
            lambda_temporal=config.lambda_temporal,
        )

        # Ensure checkpoint directory
        os.makedirs(config.checkpoint_dir, exist_ok=True)

        logger.info(
            "GNNTrainer initialized on %s | Model: %d params | Config: %s",
            self.device,
            model.num_parameters,
            config.device,
        )

    # ── Utility: move graph to device ───────────────────────────

    def _to_device(self, graph: "Data") -> "Data":
        """Move a PyG Data object to the trainer's device."""
        g = graph.clone()
        for key in g.keys():
            val = g[key]
            if isinstance(val, torch.Tensor):
                g[key] = val.to(self.device)
        return g

    # ── Phase 1: Self-Supervised Pretraining ────────────────────

    def pretrain(
        self,
        graphs: List["Data"],
        val_graphs: Optional[List["Data"]] = None,
        n_epochs: Optional[int] = None,
        mask_ratio: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Phase 1: Self-supervised masked node prediction.

        Trains only on reconstruction loss (MSE between original and
        reconstructed masked node features). No labels required —
        the model learns the "normal" methylation network from healthy
        control samples.

        Parameters
        ----------
        graphs : list of Data
            Training graphs (typically healthy controls only).
        val_graphs : list of Data or None
            Validation graphs for early stopping.
        n_epochs : int or None
            Number of pretraining epochs. Defaults to config value.
        mask_ratio : float or None
            Fraction of nodes to mask per graph. Defaults to config value.

        Returns
        -------
        metrics : dict
            Final epoch training metrics.
        """
        self.phase = GNNTrainerPhase.PRETRAINING
        n_epochs = n_epochs or self.config.n_epochs_pretrain
        mask_ratio = mask_ratio or self.config.mask_ratio_pretrain

        # Cosine annealing for pretraining
        self.scheduler = CosineAnnealingLR(self.optimizer, T_max=n_epochs)

        logger.info(
            "Phase 1: Pretraining | %d epochs | %d graphs | mask_ratio=%.2f",
            n_epochs, len(graphs), mask_ratio,
        )

        for epoch in range(1, n_epochs + 1):
            self.current_epoch = epoch
            train_metrics = self._pretrain_epoch(graphs, mask_ratio)
            self.scheduler.step()

            val_metrics = None
            if val_graphs:
                val_metrics = self._evaluate_reconstruction(val_graphs)

            self.history.record(
                epoch=epoch,
                train_metrics=train_metrics,
                val_metrics=val_metrics,
                lr=self.optimizer.param_groups[0]["lr"],
            )

            # Logging
            if epoch % self.config.log_every == 0 or epoch == 1:
                log_msg = (
                    f"Epoch {epoch}/{n_epochs} | "
                    f"L_rec={train_metrics['loss']:.4f}"
                )
                if val_metrics:
                    log_msg += f" | Val L_rec={val_metrics['loss']:.4f}"
                logger.info(log_msg)

            # Checkpoint best model
            if val_metrics:
                val_loss = val_metrics["loss"]
                if val_loss < self.best_val_loss:
                    self.best_val_loss = val_loss
                    self.best_epoch = epoch
                    self._save_checkpoint("pretrain_best.pt")
            elif epoch % (n_epochs // 5) == 0:
                self._save_checkpoint("pretrain_best.pt")

            # Early stopping
            if (
                self.config.patience > 0
                and val_graphs
                and epoch - self.best_epoch > self.config.patience
            ):
                logger.info("Early stopping at epoch %d (patience=%d)", epoch, self.config.patience)
                break

        # Save final pretrained model
        self._save_checkpoint("pretrain_final.pt")
        logger.info("Pretraining complete. Best val loss: %.4f at epoch %d",
                     self.best_val_loss, self.best_epoch)

        return {"phase": "pretrain", "best_val_loss": self.best_val_loss}

    def _pretrain_epoch(self, graphs: List["Data"], mask_ratio: float) -> Dict[str, float]:
        """One epoch of masked node prediction."""
        self.model.train()
        total_loss = 0.0

        for graph in graphs:
            graph = self._to_device(graph)
            reconstructed, mask, x_original, _ = self.model.pretrain_forward(
                graph.x, graph.edge_index,
                edge_type=getattr(graph, "edge_type", None),
                mask_ratio=mask_ratio,
            )

            loss = F.mse_loss(reconstructed[mask], x_original[mask])

            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            total_loss += loss.item()

        return {"loss": total_loss / max(1, len(graphs))}

    def _evaluate_reconstruction(self, graphs: List["Data"]) -> Dict[str, float]:
        """Evaluate reconstruction error on validation set."""
        self.model.eval()
        total_loss = 0.0

        with torch.no_grad():
            for graph in graphs:
                graph = self._to_device(graph)
                reconstructed, mask, x_original, _ = self.model.pretrain_forward(
                    graph.x, graph.edge_index,
                    edge_type=getattr(graph, "edge_type", None),
                    mask_ratio=self.config.mask_ratio_pretrain,
                )
                loss = F.mse_loss(reconstructed[mask], x_original[mask])
                total_loss += loss.item()

        return {"loss": total_loss / max(1, len(graphs))}

    # ── Phase 2: Joint Reconstruction + Anomaly Training ─────────

    def finetune(
        self,
        train_graphs: List["Data"],
        val_graphs: Optional[List["Data"]] = None,
        n_epochs: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Phase 2: Joint training with reconstruction + anomaly losses.

        Trains both the decoder and anomaly head using supervised labels.
        Reconstruction loss keeps the model grounded in normal methylation
        patterns; anomaly loss differentiates cancer from healthy.

        Parameters
        ----------
        train_graphs : list of Data
            Training graphs with labels (y=0 healthy, y=1 cancer).
        val_graphs : list of Data or None
            Validation graphs for early stopping and AUC monitoring.
        n_epochs : int or None
            Number of finetuning epochs. Defaults to config value.

        Returns
        -------
        metrics : dict
            Final epoch metrics including best validation AUC.
        """
        self.phase = GNNTrainerPhase.JOINT_FINETUNING
        n_epochs = n_epochs or self.config.n_epochs_finetune

        # ReduceLROnPlateau for finetune phase
        self.scheduler = ReduceLROnPlateau(
            self.optimizer, mode="min", factor=0.5, patience=5
        )

        self.best_epoch = 0

        logger.info(
            "Phase 2: Joint finetuning | %d epochs | %d train graphs | %d val graphs",
            n_epochs, len(train_graphs), len(val_graphs) if val_graphs else 0,
        )

        for epoch in range(1, n_epochs + 1):
            self.current_epoch = epoch
            train_metrics = self._finetune_epoch(train_graphs)

            val_metrics = None
            if val_graphs:
                val_metrics = self.evaluate(val_graphs)
                self.scheduler.step(val_metrics["loss_total"])
            else:
                self.scheduler.step(train_metrics["loss_total"])

            self.history.record(
                epoch=epoch,
                train_metrics=train_metrics,
                val_metrics=val_metrics,
                lr=self.optimizer.param_groups[0]["lr"],
            )

            # Logging
            if epoch % self.config.log_every == 0 or epoch == 1:
                log_msg = (
                    f"Epoch {epoch}/{n_epochs} | "
                    f"L_recon={train_metrics['loss_recon']:.3f} | "
                    f"L_ano={train_metrics['loss_anomaly']:.3f}"
                )
                if val_metrics:
                    log_msg += (
                        f" | Val AUC={val_metrics.get('auc', 0):.3f}"
                        f" | Val Loss={val_metrics['loss_total']:.3f}"
                    )
                logger.info(log_msg)

            # Track best by AUC (preferred) or loss
            if val_metrics:
                current_metric = val_metrics.get("auc", 0)
                if current_metric > self.best_val_auc:
                    self.best_val_auc = current_metric
                    self.best_epoch = epoch
                    self._save_checkpoint("finetune_best.pt")
            elif epoch % max(1, n_epochs // 5) == 0:
                self._save_checkpoint("finetune_best.pt")

            # Early stopping
            if (
                self.config.patience > 0
                and val_graphs
                and epoch - self.best_epoch > self.config.patience
            ):
                logger.info("Early stopping at epoch %d", epoch)
                break

        self._save_checkpoint("finetune_final.pt")
        logger.info("Finetuning complete. Best val AUC: %.4f at epoch %d",
                     self.best_val_auc, self.best_epoch)

        return {
            "phase": "finetune",
            "best_val_auc": self.best_val_auc,
            "best_epoch": self.best_epoch,
        }

    def _finetune_epoch(self, graphs: List["Data"]) -> Dict[str, float]:
        """One epoch of joint reconstruction + anomaly training."""
        self.model.train()
        metrics_accum = {"loss_recon": 0.0, "loss_anomaly": 0.0, "loss_total": 0.0}

        for graph in graphs:
            graph = self._to_device(graph)
            output = self.model(
                graph.x, graph.edge_index,
                edge_type=getattr(graph, "edge_type", None),
            )

            loss_total, loss_components = self.loss_fn(
                x_original=graph.x,
                x_reconstructed=output["reconstructed"],
                anomaly_scores=output["anomaly_scores"],
                labels=graph.y,
                mask=None,  # all-node reconstruction in finetune
            )

            self.optimizer.zero_grad()
            loss_total.backward()
            self.optimizer.step()

            for k in metrics_accum:
                metrics_accum[k] += loss_components[k]

        n = max(1, len(graphs))
        return {k: v / n for k, v in metrics_accum.items()}

    # ── Evaluation ───────────────────────────────────────────────

    @torch.no_grad()
    def evaluate(self, graphs: List["Data"]) -> Dict[str, float]:
        """
        Evaluate model on a set of graphs.

        Computes reconstruction error, anomaly predictions, and
        classification metrics (AUC, AUPRC).

        Parameters
        ----------
        graphs : list of Data
            Evaluation graphs with labels.

        Returns
        -------
        metrics : dict
            Keys: loss_recon, loss_anomaly, loss_total, auc, auprc,
            mean_anomaly_cancer, mean_anomaly_healthy.
        """
        self.model.eval()
        metrics_accum = {"loss_recon": 0.0, "loss_anomaly": 0.0, "loss_total": 0.0}
        all_scores: List[float] = []
        all_labels: List[int] = []
        all_recon = 0.0

        for graph in graphs:
            graph = self._to_device(graph)
            output = self.model(
                graph.x, graph.edge_index,
                edge_type=getattr(graph, "edge_type", None),
            )

            # Loss components
            loss_total, loss_components = self.loss_fn(
                x_original=graph.x,
                x_reconstructed=output["reconstructed"],
                anomaly_scores=output["anomaly_scores"],
                labels=graph.y,
            )

            for k in metrics_accum:
                metrics_accum[k] += loss_components[k]

            # Anomaly scores and labels
            graph_score = float(output["graph_anomaly"].squeeze().cpu())
            all_scores.append(graph_score)
            all_labels.append(int(graph.y.item()))

        n = max(1, len(graphs))
        result = {k: v / n for k, v in metrics_accum.items()}

        # Classification metrics
        if _HAS_SKLEARN and len(set(all_labels)) > 1:
            result["auc"] = float(roc_auc_score(all_labels, all_scores))
            result["auprc"] = float(
                average_precision_score(all_labels, all_scores)
            )
        else:
            result["auc"] = 0.5
            result["auprc"] = 0.5

        # Per-class anomaly statistics
        cancer_scores = [s for s, l in zip(all_scores, all_labels) if l == 1]
        healthy_scores = [s for s, l in zip(all_scores, all_labels) if l == 0]
        result["mean_anomaly_cancer"] = float(np.mean(cancer_scores)) if cancer_scores else 0.0
        result["mean_anomaly_healthy"] = float(np.mean(healthy_scores)) if healthy_scores else 0.0

        # Effect size
        all_scores_np = np.array(all_scores)
        result["anomaly_std"] = float(np.std(all_scores_np))
        result["anomaly_range"] = float(np.max(all_scores_np) - np.min(all_scores_np))

        return result

    # ── Prediction ──────────────────────────────────────────────

    @torch.no_grad()
    def predict(
        self, graph: "Data", return_details: bool = False
    ) -> Union[float, Dict[str, Any]]:
        """
        Predict field defect score for a single graph.

        Parameters
        ----------
        graph : Data
            Input graph.
        return_details : bool
            If True, return dict with score, node scores, and reconstruction.

        Returns
        -------
        float or dict
            Field defect score in [0, 1] (higher = more cancer-like).
        """
        self.model.eval()
        graph = self._to_device(graph)

        output = self.model(
            graph.x, graph.edge_index,
            edge_type=getattr(graph, "edge_type", None),
        )

        score = float(output["graph_anomaly"].squeeze().cpu())

        if not return_details:
            return score

        field_score = float(
            compute_field_defect_score(
                graph.x,
                output["reconstructed"],
                output["anomaly_scores"],
            )
            .cpu()
        )

        return {
            "field_defect_score": field_score,
            "graph_anomaly_score": score,
            "node_anomaly_scores": output["anomaly_scores"].squeeze(-1).cpu().numpy(),
            "reconstructed": output["reconstructed"].cpu().numpy(),
        }

    # ── Checkpointing ───────────────────────────────────────────

    def _save_checkpoint(self, filename: str) -> None:
        """Save model state and training metadata."""
        path = os.path.join(self.config.checkpoint_dir, filename)
        checkpoint = {
            "model_state_dict": self.model.state_dict(),
            "optimizer_state_dict": self.optimizer.state_dict(),
            "epoch": self.current_epoch,
            "phase": self.phase.value,
            "best_val_loss": self.best_val_loss,
            "best_val_auc": self.best_val_auc,
            "config": self.config.to_dict(),
            "history": {
                "epoch": self.history.epoch,
                "train_loss": self.history.train_loss,
                "val_loss": self.history.val_loss,
                "val_auc": self.history.val_auc,
            },
        }
        torch.save(checkpoint, path)
        logger.debug("Checkpoint saved: %s", path)

    def load_checkpoint(self, path: str, load_optimizer: bool = True) -> int:
        """
        Load model and training state from checkpoint.

        Parameters
        ----------
        path : str
            Path to checkpoint file (.pt).
        load_optimizer : bool
            If True, also restore optimizer state (for resuming training).

        Returns
        -------
        epoch : int
            The epoch the checkpoint was saved at.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"Checkpoint not found: {path}")

        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(checkpoint["model_state_dict"])

        if load_optimizer and "optimizer_state_dict" in checkpoint:
            self.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])

        self.current_epoch = checkpoint.get("epoch", 0)
        self.best_val_loss = checkpoint.get("best_val_loss", float("inf"))
        self.best_val_auc = checkpoint.get("best_val_auc", 0.0)
        phase_str = checkpoint.get("phase", "uninitialized")
        self.phase = GNNTrainerPhase(phase_str)

        logger.info(
            "Loaded checkpoint from %s (epoch %d, phase=%s, best_auc=%.4f)",
            path, self.current_epoch, phase_str, self.best_val_auc,
        )
        return self.current_epoch

    # ── History viz (lightweight) ────────────────────────────────

    def get_best_epoch(self) -> Dict[str, Any]:
        """Return metrics from the best epoch (by validation AUC)."""
        if not self.history.val_auc or all(v is None for v in self.history.val_auc):
            return {"best_epoch": None, "best_val_auc": None}

        valid_aucs = [(i, v) for i, v in enumerate(self.history.val_auc) if v is not None]
        if not valid_aucs:
            return {"best_epoch": None, "best_val_auc": None}

        best_idx, best_auc = max(valid_aucs, key=lambda x: x[1])
        return {
            "best_epoch": self.history.epoch[best_idx],
            "best_val_auc": best_auc,
            "best_val_loss": self.history.val_loss[best_idx],
            "train_loss": self.history.train_loss[best_idx],
        }
