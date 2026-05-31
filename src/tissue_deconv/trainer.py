#!/usr/bin/env python3
"""
Tissue Deconvolution Trainer
==============================

Training and evaluation utilities for the TissueDeconvolutionModel.

Handles:
- Training on synthetic (or real) cfDNA mixtures
- Evaluation metrics (Pearson r, MAE per tissue)
- Checkpoint save/load
- Synthetic training data generation
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import torch
    import torch.optim as optim
    from torch.utils.data import TensorDataset, DataLoader

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None
    optim = None
    TensorDataset = None
    DataLoader = None

from .config import (
    TissueDeconvConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    SUPPORTED_TISSUES,
)
from .model import (
    TissueDeconvolutionModel,
    TissueDeconvolutionEnsemble,
    DeconvLoss,
)
from .tissue_atlas import TissueAtlas

logger = logging.getLogger(__name__)


class TissueDeconvTrainer:
    """
    Trainer for tissue deconvolution models.

    Supports both single-model and ensemble training. Uses synthetic
    data generation from TissueAtlas when real training data is unavailable.

    Parameters
    ----------
    config : TissueDeconvConfig
        Training configuration.
    device : str, optional
        Override device (cpu, cuda, mps, auto).
    """

    def __init__(
        self,
        config: Optional[TissueDeconvConfig] = None,
        device: Optional[str] = None,
    ):
        if not _HAS_TORCH:
            raise ImportError(
                "PyTorch is required for TissueDeconvTrainer. "
                "Install with: pip install torch"
            )

        self.config = config if config is not None else DEFAULT_CONFIG
        self.device_str = device or self.config.device
        self.device = torch.device(self._resolve_device())

        # Atlas for synthetic data generation
        self.atlas = TissueAtlas(
            n_cpg_features=self.config.n_cpg_features,
            random_seed=42,
        )

        # Loss function
        self.criterion = DeconvLoss(self.config)

        # Ensemble
        self.ensemble: Optional[TissueDeconvolutionEnsemble] = None
        self._trained_models: List[TissueDeconvolutionModel] = []

        # Training history
        self.history: Dict[str, List[float]] = {
            "train_loss": [],
            "val_loss": [],
            "val_pcc": [],      # Pearson correlation
            "val_mae": [],      # Mean absolute error
        }

    def _resolve_device(self) -> str:
        """Resolve device string."""
        if self.device_str not in ("cuda", "mps", "cpu", "auto"):
            return "cpu"
        if self.device_str == "auto":
            try:
                if torch.cuda.is_available():
                    return "cuda"
                if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                    return "mps"
            except Exception:
                pass
            return "cpu"
        if self.device_str == "cuda" and not torch.cuda.is_available():
            logger.warning("CUDA requested but not available. Falling back to CPU.")
            return "cpu"
        if self.device_str == "mps" and not (
            hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
        ):
            logger.warning("MPS requested but not available. Falling back to CPU.")
            return "cpu"
        return self.device_str

    # ── Synthetic data generation ────────────────────────────────

    def generate_synthetic_data(
        self,
        n_samples: int = 2000,
        noise: float = 0.01,
        val_split: float = 0.2,
        seed: int = 42,
    ) -> Tuple[
        Tuple[np.ndarray, np.ndarray],
        Tuple[np.ndarray, np.ndarray],
    ]:
        """
        Generate synthetic training and validation data.

        Parameters
        ----------
        n_samples : int
            Total number of synthetic samples.
        noise : float
            Measurement noise level.
        val_split : float
            Fraction of data for validation.
        seed : int
            Random seed.

        Returns
        -------
        (train_mix, train_frac), (val_mix, val_frac) : tuples
        """
        self.atlas.load_reference()  # Ensure atlas is loaded

        mixtures, fractions = self.atlas.generate_training_dataset(
            n_samples=n_samples,
            noise=noise,
            seed=seed,
        )

        n_val = int(n_samples * val_split)
        indices = np.random.RandomState(seed).permutation(n_samples)

        train_idx = indices[n_val:]
        val_idx = indices[:n_val]

        return (
            (mixtures[train_idx], fractions[train_idx]),
            (mixtures[val_idx], fractions[val_idx]),
        )

    # ── Single model training ────────────────────────────────────

    def train_model(
        self,
        model: TissueDeconvolutionModel,
        train_mixtures: np.ndarray,
        train_fractions: np.ndarray,
        val_mixtures: Optional[np.ndarray] = None,
        val_fractions: Optional[np.ndarray] = None,
        n_epochs: Optional[int] = None,
        batch_size: Optional[int] = None,
        learning_rate: Optional[float] = None,
        verbose: bool = True,
    ) -> Dict[str, List[float]]:
        """
        Train a single TissueDeconvolutionModel.

        Parameters
        ----------
        model : TissueDeconvolutionModel
            Model to train.
        train_mixtures : (n_train, n_cpg_features)
            Training beta-value mixtures.
        train_fractions : (n_train, n_tissues)
            Ground truth tissue fractions.
        val_mixtures : (n_val, n_cpg_features), optional
            Validation mixtures.
        val_fractions : (n_val, n_tissues), optional
            Validation fractions.
        n_epochs : int, optional
            Override config epochs.
        batch_size : int, optional
            Override config batch size.
        learning_rate : float, optional
            Override config learning rate.
        verbose : bool
            Print progress.

        Returns
        -------
        history : dict
            Training history with loss/metrics per epoch.
        """
        n_epochs = n_epochs or self.config.n_epochs
        batch_size = batch_size or self.config.batch_size
        lr = learning_rate or self.config.learning_rate
        patience = self.config.patience

        model.to(self.device)
        model.train()

        # Data loaders
        train_dataset = TensorDataset(
            torch.as_tensor(train_mixtures, dtype=torch.float32),
            torch.as_tensor(train_fractions, dtype=torch.float32),
        )
        train_loader = DataLoader(
            train_dataset, batch_size=batch_size, shuffle=True, drop_last=False,
        )

        optimizer = optim.AdamW(
            model.parameters(), lr=lr, weight_decay=self.config.weight_decay,
        )
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='min', factor=0.5, patience=5,
        )

        history: Dict[str, List[float]] = {
            "train_loss": [],
            "val_loss": [],
            "val_pcc": [],
            "val_mae": [],
        }

        best_val_loss = float("inf")
        best_state = None
        patience_counter = 0

        for epoch in range(n_epochs):
            # ── Training ─────────────────────────────────────────
            model.train()
            epoch_losses = []

            for batch_x, batch_y in train_loader:
                batch_x = batch_x.to(self.device)
                batch_y = batch_y.to(self.device)

                optimizer.zero_grad()
                pred = model(batch_x)
                loss, components = self.criterion(pred, batch_y)
                loss.backward()
                optimizer.step()

                epoch_losses.append(float(loss.detach().cpu()))

            avg_train_loss = np.mean(epoch_losses)
            history["train_loss"].append(avg_train_loss)

            # ── Validation ───────────────────────────────────────
            if val_mixtures is not None and val_fractions is not None:
                val_loss, val_pcc, val_mae = self._evaluate_model(
                    model, val_mixtures, val_fractions, batch_size
                )
                history["val_loss"].append(val_loss)
                history["val_pcc"].append(val_pcc)
                history["val_mae"].append(val_mae)

                scheduler.step(val_loss)

                # Early stopping
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_state = {
                        k: v.cpu().clone() for k, v in model.state_dict().items()
                    }
                    patience_counter = 0
                else:
                    patience_counter += 1

                if verbose and (epoch % max(1, n_epochs // 10) == 0 or epoch == n_epochs - 1):
                    logger.info(
                        f"Epoch {epoch+1:3d}/{n_epochs} | "
                        f"train_loss: {avg_train_loss:.4f} | "
                        f"val_loss: {val_loss:.4f} | "
                        f"val_pcc: {val_pcc:.4f} | "
                        f"val_mae: {val_mae:.4f}"
                    )
            else:
                patience_counter = 0

            if patience_counter >= patience:
                logger.info(f"Early stopping at epoch {epoch+1}")
                break

        # Restore best model
        if best_state is not None:
            model.load_state_dict(best_state)

        self.history = history
        return history

    def _evaluate_model(
        self,
        model: TissueDeconvolutionModel,
        mixtures: np.ndarray,
        fractions: np.ndarray,
        batch_size: int,
    ) -> Tuple[float, float, float]:
        """Evaluate model on a dataset. Returns (loss, pcc, mae)."""
        model.eval()
        dataset = TensorDataset(
            torch.as_tensor(mixtures, dtype=torch.float32),
            torch.as_tensor(fractions, dtype=torch.float32),
        )
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

        total_loss = 0.0
        all_preds = []
        all_targets = []

        with torch.no_grad():
            for batch_x, batch_y in loader:
                batch_x = batch_x.to(self.device)
                batch_y = batch_y.to(self.device)
                pred = model(batch_x)
                loss, _ = self.criterion(pred, batch_y)
                total_loss += float(loss.cpu()) * len(batch_x)
                all_preds.append(pred.cpu().numpy())
                all_targets.append(batch_y.cpu().numpy())

        avg_loss = total_loss / len(mixtures)

        preds = np.concatenate(all_preds, axis=0)
        targets = np.concatenate(all_targets, axis=0)

        # Pearson correlation (flattened, per-sample mean)
        pccs = []
        for i in range(len(preds)):
            if np.std(preds[i]) > 1e-8 and np.std(targets[i]) > 1e-8:
                pcc = np.corrcoef(preds[i], targets[i])[0, 1]
                if not np.isnan(pcc):
                    pccs.append(pcc)
        avg_pcc = np.mean(pccs) if pccs else 0.0

        # MAE
        mae = np.mean(np.abs(preds - targets))

        return avg_loss, avg_pcc, mae

    # ── Ensemble training ────────────────────────────────────────

    def train_ensemble(
        self,
        train_mixtures: np.ndarray,
        train_fractions: np.ndarray,
        val_mixtures: Optional[np.ndarray] = None,
        val_fractions: Optional[np.ndarray] = None,
        n_epochs_per_model: Optional[int] = None,
        verbose: bool = True,
    ) -> TissueDeconvolutionEnsemble:
        """
        Train an ensemble of TissueDeconvolutionModel instances.

        Parameters
        ----------
        train_mixtures, train_fractions : arrays
            Training data.
        val_mixtures, val_fractions : arrays, optional
            Validation data.
        n_epochs_per_model : int, optional
            Epochs per ensemble member.
        verbose : bool
            Print progress.

        Returns
        -------
        ensemble : TissueDeconvolutionEnsemble
            Trained ensemble.
        """
        ensemble = TissueDeconvolutionEnsemble(
            config=self.config,
            base_seed=42,
        )
        ensemble.to(self.device)

        for i, model in enumerate(ensemble.models):
            logger.info(f"Training ensemble member {i+1}/{len(ensemble.models)}")
            self.train_model(
                model=model,
                train_mixtures=train_mixtures,
                train_fractions=train_fractions,
                val_mixtures=val_mixtures,
                val_fractions=val_fractions,
                n_epochs=n_epochs_per_model,
                verbose=verbose,
            )
            self._trained_models.append(model)

        self.ensemble = ensemble
        return ensemble

    # ── Prediction ───────────────────────────────────────────────

    def predict(
        self,
        mixture: np.ndarray,
    ) -> Dict[str, float]:
        """
        Predict tissue fractions for a single cfDNA sample.

        Uses ensemble if available, otherwise the last trained model.

        Parameters
        ----------
        mixture : (n_cpg_features,) array
            Methylation beta values.

        Returns
        -------
        fractions : dict
            {tissue_name: fraction} for non-zero tissues.
        """
        if mixture.ndim == 1:
            mixture = mixture.reshape(1, -1)

        if self.ensemble is not None:
            return self.ensemble.predict_sample(mixture)

        if self._trained_models:
            model = self._trained_models[-1]
            preds = model.predict(torch.as_tensor(mixture, dtype=torch.float32))
        else:
            # No trained model: return uniform
            preds = np.ones((1, len(SUPPORTED_TISSUES))) / len(SUPPORTED_TISSUES)

        pred = preds[0]
        top_indices = np.argsort(pred)[::-1][:5]
        result = {}
        for idx in top_indices:
            frac = float(pred[idx])
            if frac > self.config.detection_limit:
                result[SUPPORTED_TISSUES[idx]] = frac

        return result

    def predict_batch(
        self,
        mixtures: np.ndarray,
    ) -> np.ndarray:
        """
        Predict tissue fractions for a batch of samples.

        Parameters
        ----------
        mixtures : (n_samples, n_cpg_features) array

        Returns
        -------
        fractions : (n_samples, n_tissues) float32 array
        """
        if self.ensemble is not None:
            return self.ensemble.predict(mixtures)

        if self._trained_models:
            model = self._trained_models[-1]
            model.eval()
            return model.predict(
                torch.as_tensor(mixtures, dtype=torch.float32)
            )

        # No model: return uniform
        n = len(mixtures)
        return np.ones((n, len(SUPPORTED_TISSUES))) / len(SUPPORTED_TISSUES)

    # ── Checkpointing ────────────────────────────────────────────

    def save_checkpoint(self, path: str) -> None:
        """
        Save trainer state (ensemble models + config + history).

        Parameters
        ----------
        path : str
            Output .pt file path.
        """
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

        checkpoint = {
            "config": self.config.to_dict(),
            "history": self.history,
            "device": self.device_str,
        }

        if self.ensemble is not None:
            checkpoint["ensemble_state_dicts"] = self.ensemble.state_dicts()

        torch.save(checkpoint, path)
        logger.info("Checkpoint saved to %s", path)

    def load_checkpoint(self, path: str) -> None:
        """
        Load trainer state from checkpoint.

        Parameters
        ----------
        path : str
            Checkpoint .pt file path.
        """
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)

        if "ensemble_state_dicts" in checkpoint:
            n_ensemble = len(checkpoint["ensemble_state_dicts"])
            # Update config if ensemble size changed
            if n_ensemble != self.config.n_ensemble:
                logger.warning(
                    "Checkpoint has %d ensemble members, config has %d. "
                    "Using checkpoint size.",
                    n_ensemble, self.config.n_ensemble,
                )
                self.config.n_ensemble = n_ensemble

            if self.ensemble is None:
                self.ensemble = TissueDeconvolutionEnsemble(
                    config=self.config,
                    base_seed=42,
                )
            self.ensemble.to(self.device)
            self.ensemble.load_state_dicts(checkpoint["ensemble_state_dicts"])
            self._trained_models = list(self.ensemble.models)

        if "history" in checkpoint:
            self.history = checkpoint["history"]

        logger.info("Checkpoint loaded from %s", path)

    # ── Quick training (convenience method) ──────────────────────

    def quick_train(
        self,
        n_samples: int = 1000,
        n_epochs: int = 50,
        force_ensemble: bool = True,
        verbose: bool = True,
    ) -> TissueDeconvolutionEnsemble:
        """
        Quick end-to-end training on synthetic data.

        Generates synthetic data, trains ensemble, and returns
        the trained ensemble ready for inference.

        Parameters
        ----------
        n_samples : int
            Number of synthetic training samples.
        n_epochs : int
            Training epochs per model.
        force_ensemble : bool
            If True, trains full ensemble even if config.n_ensemble=1.
        verbose : bool
            Print progress.

        Returns
        -------
        ensemble : TissueDeconvolutionEnsemble
        """
        logger.info("Generating %d synthetic samples...", n_samples)
        (train_x, train_y), (val_x, val_y) = self.generate_synthetic_data(
            n_samples=n_samples,
        )

        logger.info(
            "Training ensemble (%d models) for %d epochs each...",
            self.config.n_ensemble, n_epochs,
        )

        self.train_ensemble(
            train_mixtures=train_x,
            train_fractions=train_y,
            val_mixtures=val_x,
            val_fractions=val_y,
            n_epochs_per_model=n_epochs,
            verbose=verbose,
        )

        return self.ensemble
