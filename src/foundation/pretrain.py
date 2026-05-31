#!/usr/bin/env python3
"""
Foundation Model Pre-training
==============================

Self-supervised pre-training for the DeepCatch Foundation Model.

Three-phase training strategy:

    Phase 1 — Masked Modality Prediction (MMP):
        Mask 30% of modalities randomly per sample.
        Train encoder + decoder to reconstruct masked modalities
        from unmasked context. This teaches the model to capture
        cross-modal dependencies.

    Phase 2 — Contrastive Learning:
        Pull embeddings of different modalities from the same
        sample together (positives), push different samples apart
        (negatives). This builds a discriminative joint embedding
        space without requiring labels.

    Phase 3 — Joint Training:
        Combine MMP and contrastive losses with configurable
        weighting. Fine-tunes the joint representation.

Checkpoints are saved after each phase for downstream fine-tuning.
"""

from __future__ import annotations

import logging
import os
import time
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from .config import FoundationConfig, DEFAULT_CONFIG, MODALITY_DIMS, MODALITY_NAMES
from .model import MultiModalEncoder, PretrainHead, ContrastiveHead
from .data import MultiModalDataGenerator

logger = logging.getLogger(__name__)


class FoundationPretrainer:
    """
    Self-supervised pre-trainer for the DeepCatch Foundation Model.

    Handles 3-phase training with synthetic data, checkpointing,
    and progress tracking.

    Parameters
    ----------
    config : FoundationConfig
        Model and training configuration.
    device : str, optional
        Compute device (default from config).
    verbose : bool
        Print progress during training.
    """

    def __init__(
        self,
        config: Optional[FoundationConfig] = None,
        device: Optional[str] = None,
        verbose: bool = False,
    ):
        self.config = config if config is not None else DEFAULT_CONFIG
        self.device = device or self.config.device
        self.verbose = verbose

        # Init models
        self.encoder = MultiModalEncoder(self.config)
        self.pretrain_head = PretrainHead(self.config)
        self.contrastive_head = ContrastiveHead(self.config)

        # Move to device
        self.encoder.to(self.device)
        self.pretrain_head.to(self.device)
        self.contrastive_head.to(self.device)

        # Data generator
        self.data_generator = MultiModalDataGenerator(seed=self.config.seed)

        # Training state
        self._is_pretrained = False
        self._pretrain_losses: List[float] = []
        self._phase_completed: List[str] = []

    def _generate_mask(self, batch_size: int) -> torch.Tensor:
        """
        Generate a random modality mask.

        Parameters
        ----------
        batch_size : int

        Returns
        -------
        mask : (batch, n_modalities) bool
            True = masked.
        """
        n = self.config.n_modalities
        # Ensure at least 1 modality remains unmasked
        mask = torch.rand(batch_size, n) < self.config.mask_ratio
        # For each sample, ensure not all modalities are masked
        all_masked = mask.all(dim=1)
        if all_masked.any():
            # Unmask a random modality for those samples
            for i in range(batch_size):
                if mask[i].all():
                    unmask_idx = torch.randint(0, n, (1,))
                    mask[i, unmask_idx] = False
        return mask.to(self.device)

    def _modalities_to_tensors(
        self,
        modalities: Dict[str, np.ndarray],
    ) -> Dict[str, torch.Tensor]:
        """Convert numpy modality dict to torch tensors."""
        return {
            k: torch.from_numpy(v.astype(np.float32)).to(self.device)
            for k, v in modalities.items()
        }

    def _log(self, msg: str):
        if self.verbose:
            logger.info(msg)

    # ── Phase 1: Masked Modality Prediction ─────────────────────

    def pretrain_phase1_mmp(
        self,
        n_samples: int = 5000,
        batch_size: Optional[int] = None,
        n_epochs: int = 50,
        lr: Optional[float] = None,
    ) -> List[float]:
        """
        Phase 1: Masked Modality Prediction.

        Train the encoder + pretrain head to predict masked modality
        features from unmasked context.

        Parameters
        ----------
        n_samples : int
            Number of synthetic samples to generate for training.
        batch_size : int, optional
            Batch size (default from config).
        n_epochs : int
            Number of training epochs.
        lr : float, optional
            Learning rate (default from config).

        Returns
        -------
        losses : list of float
            Epoch losses.
        """
        batch_size = batch_size or self.config.batch_size
        lr = lr or self.config.pretrain_lr

        self._log(f"Phase 1: MMP — {n_samples} samples, {n_epochs} epochs")

        # Generate synthetic data
        modalities_np, _ = self.data_generator.generate_dataset(
            n_samples=n_samples,
            prefix="phase1",
        )
        modalities = self._modalities_to_tensors(modalities_np)

        # Optimizer
        params = list(self.encoder.parameters()) + list(self.pretrain_head.parameters())
        optimizer = torch.optim.AdamW(params, lr=lr)

        losses = []
        n_batches = max(1, n_samples // batch_size)

        for epoch in range(n_epochs):
            epoch_loss = 0.0

            for batch_idx in range(n_batches):
                start = (batch_idx * batch_size) % n_samples
                end = min(start + batch_size, n_samples)
                batch_size_actual = end - start

                # Extract batch
                batch_modalities = {
                    k: v[start:end] for k, v in modalities.items()
                }

                # Generate mask
                mask = self._generate_mask(batch_size_actual)

                # Forward: encoder with mask
                joint = self.encoder(batch_modalities, mask=mask)

                # Forward: pretrain head
                reconstructed = self.pretrain_head(joint, mask)

                # Loss
                loss = self.pretrain_head.compute_loss(
                    reconstructed, batch_modalities, mask
                )
                if loss.item() == 0:
                    continue

                # Backward
                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(params, max_norm=1.0)
                optimizer.step()

                epoch_loss += loss.item()

            avg_loss = epoch_loss / max(1, n_batches)
            losses.append(avg_loss)

            if self.verbose and (epoch + 1) % 10 == 0:
                self._log(f"  Epoch {epoch+1}/{n_epochs} — MMP loss: {avg_loss:.6f}")

        self._pretrain_losses.extend(losses)
        self._phase_completed.append("phase1_mmp")
        self._is_pretrained = True

        if self.verbose:
            self._log(f"  Phase 1 complete — final MMP loss: {losses[-1]:.6f}")

        return losses

    # ── Phase 2: Contrastive Learning ───────────────────────────

    def pretrain_phase2_contrastive(
        self,
        n_samples: int = 5000,
        batch_size: Optional[int] = None,
        n_epochs: int = 50,
        lr: Optional[float] = None,
    ) -> List[float]:
        """
        Phase 2: Cross-modal contrastive learning.

        Train encoder + contrastive head to pull same-sample modalities
        together and push different-sample ones apart.

        Parameters
        ----------
        n_samples : int
            Number of synthetic samples.
        batch_size : int, optional
            Batch size.
        n_epochs : int
            Number of epochs.
        lr : float, optional
            Learning rate.

        Returns
        -------
        losses : list of float
            Epoch losses.
        """
        batch_size = batch_size or self.config.batch_size
        lr = lr or self.config.pretrain_lr

        self._log(f"Phase 2: Contrastive — {n_samples} samples, {n_epochs} epochs")

        # Generate synthetic data
        modalities_np, _ = self.data_generator.generate_dataset(
            n_samples=n_samples,
            prefix="phase2",
        )
        modalities = self._modalities_to_tensors(modalities_np)

        # Optimizer
        params = list(self.encoder.parameters()) + list(self.contrastive_head.parameters())
        optimizer = torch.optim.AdamW(params, lr=lr)

        losses = []
        n_batches = max(1, n_samples // batch_size)

        for epoch in range(n_epochs):
            epoch_loss = 0.0
            n_batches_done = 0

            for batch_idx in range(n_batches):
                start = (batch_idx * batch_size) % n_samples
                end = min(start + batch_size, n_samples)
                batch_size_actual = end - start

                # Need at least 2 samples for contrastive
                if batch_size_actual < 2:
                    continue

                batch_modalities = {
                    k: v[start:end] for k, v in modalities.items()
                }

                # Forward: encoder (no mask)
                joint = self.encoder(batch_modalities)

                # Contrastive loss
                loss = self.contrastive_head(joint)
                if loss.item() == 0:
                    continue

                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(params, max_norm=1.0)
                optimizer.step()

                epoch_loss += loss.item()
                n_batches_done += 1

            avg_loss = epoch_loss / max(1, n_batches_done)
            losses.append(avg_loss)

            if self.verbose and (epoch + 1) % 10 == 0:
                self._log(f"  Epoch {epoch+1}/{n_epochs} — Contrastive loss: {avg_loss:.6f}")

        self._pretrain_losses.extend(losses)
        self._phase_completed.append("phase2_contrastive")
        self._is_pretrained = True

        if self.verbose:
            self._log(f"  Phase 2 complete — final contrastive loss: {losses[-1]:.6f}")

        return losses

    # ── Phase 3: Joint Training ─────────────────────────────────

    def pretrain_phase3_joint(
        self,
        n_samples: int = 5000,
        batch_size: Optional[int] = None,
        n_epochs: int = 30,
        lr: Optional[float] = None,
    ) -> List[float]:
        """
        Phase 3: Joint MMP + contrastive training.

        Parameters
        ----------
        n_samples : int
            Number of synthetic samples.
        batch_size : int, optional
        n_epochs : int
        lr : float, optional

        Returns
        -------
        losses : list of float
            Epoch losses (total combined).
        """
        batch_size = batch_size or self.config.batch_size
        lr = lr or (self.config.pretrain_lr * 0.5)

        self._log(f"Phase 3: Joint — {n_samples} samples, {n_epochs} epochs")

        modalities_np, _ = self.data_generator.generate_dataset(
            n_samples=n_samples,
            prefix="phase3",
        )
        modalities = self._modalities_to_tensors(modalities_np)

        params = (
            list(self.encoder.parameters())
            + list(self.pretrain_head.parameters())
            + list(self.contrastive_head.parameters())
        )
        optimizer = torch.optim.AdamW(params, lr=lr)

        losses = []
        n_batches = max(1, n_samples // batch_size)

        for epoch in range(n_epochs):
            epoch_total_loss = 0.0
            n_batches_done = 0

            for batch_idx in range(n_batches):
                start = (batch_idx * batch_size) % n_samples
                end = min(start + batch_size, n_samples)
                batch_size_actual = end - start

                if batch_size_actual < 2:
                    continue

                batch_modalities = {
                    k: v[start:end] for k, v in modalities.items()
                }

                # Mask for MMP
                mask = self._generate_mask(batch_size_actual)

                # Forward: encoder with mask
                joint = self.encoder(batch_modalities, mask=mask)

                # MMP loss
                reconstructed = self.pretrain_head(joint, mask)
                mmp_loss = self.pretrain_head.compute_loss(
                    reconstructed, batch_modalities, mask
                )

                # Contrastive loss
                contrast_loss = self.contrastive_head(joint)

                # Combined loss
                total_loss = (
                    self.config.lambda_mask * mmp_loss
                    + self.config.lambda_contrast * contrast_loss
                )

                if total_loss.item() == 0:
                    continue

                optimizer.zero_grad()
                total_loss.backward()
                torch.nn.utils.clip_grad_norm_(params, max_norm=1.0)
                optimizer.step()

                epoch_total_loss += total_loss.item()
                n_batches_done += 1

            avg_loss = epoch_total_loss / max(1, n_batches_done)
            losses.append(avg_loss)

            if self.verbose and (epoch + 1) % 10 == 0:
                self._log(
                    f"  Epoch {epoch+1}/{n_epochs} — "
                    f"Joint loss: {avg_loss:.6f}"
                )

        self._pretrain_losses.extend(losses)
        self._phase_completed.append("phase3_joint")
        self._is_pretrained = True

        if self.verbose:
            self._log(f"  Phase 3 complete — final joint loss: {losses[-1]:.6f}")

        return losses

    # ── Full pre-training pipeline ───────────────────────────────

    def pretrain(
        self,
        n_samples: int = 5000,
        p1_epochs: int = 50,
        p2_epochs: int = 50,
        p3_epochs: int = 30,
        batch_size: Optional[int] = None,
    ) -> Dict[str, List[float]]:
        """
        Run all three pre-training phases.

        Parameters
        ----------
        n_samples : int
            Synthetic training samples per phase.
        p1_epochs : int
            Phase 1 epochs.
        p2_epochs : int
            Phase 2 epochs.
        p3_epochs : int
            Phase 3 epochs.
        batch_size : int, optional

        Returns
        -------
        losses : dict
            {'phase1': [...], 'phase2': [...], 'phase3': [...]}
        """
        self._log("=" * 50)
        self._log("DeepCatch Foundation Model Pre-training")
        self._log(f"Config: {self.config.n_layers} layers, "
                   f"{self.config.embed_dim}d embeddings, "
                   f"{self.config.n_heads} heads")
        self._log(f"Params: {self.encoder.num_params:,}")
        self._log("=" * 50)

        p1_losses = self.pretrain_phase1_mmp(
            n_samples=n_samples, n_epochs=p1_epochs, batch_size=batch_size
        )
        p2_losses = self.pretrain_phase2_contrastive(
            n_samples=n_samples, n_epochs=p2_epochs, batch_size=batch_size
        )
        p3_losses = self.pretrain_phase3_joint(
            n_samples=n_samples, n_epochs=p3_epochs, batch_size=batch_size
        )

        return {
            "phase1": p1_losses,
            "phase2": p2_losses,
            "phase3": p3_losses,
        }

    # ── Checkpointing ────────────────────────────────────────────

    def save_checkpoint(self, path: str):
        """
        Save pre-training checkpoint.

        Parameters
        ----------
        path : str
            Output file path (.pt).
        """
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        checkpoint = {
            "encoder_state_dict": self.encoder.state_dict(),
            "pretrain_head_state_dict": self.pretrain_head.state_dict(),
            "contrastive_head_state_dict": self.contrastive_head.state_dict(),
            "config": self.config.to_dict(),
            "is_pretrained": self._is_pretrained,
            "phase_completed": self._phase_completed,
            "pretrain_losses": self._pretrain_losses,
        }
        torch.save(checkpoint, path)
        self._log(f"Checkpoint saved → {path}")

    def load_checkpoint(self, path: str) -> bool:
        """
        Load pre-training checkpoint.

        Parameters
        ----------
        path : str
            Checkpoint file path.

        Returns
        -------
        success : bool
        """
        if not os.path.exists(path):
            self._log(f"Checkpoint not found: {path}")
            return False

        checkpoint = torch.load(path, map_location=self.device)

        self.encoder.load_state_dict(checkpoint["encoder_state_dict"])
        self.pretrain_head.load_state_dict(checkpoint["pretrain_head_state_dict"])
        self.contrastive_head.load_state_dict(checkpoint["contrastive_head_state_dict"])
        self._is_pretrained = checkpoint.get("is_pretrained", True)
        self._phase_completed = checkpoint.get("phase_completed", [])
        self._pretrain_losses = checkpoint.get("pretrain_losses", [])

        self._log(f"Checkpoint loaded ← {path}")
        return True

    # ── Properties ───────────────────────────────────────────────

    @property
    def is_pretrained(self) -> bool:
        return self._is_pretrained

    @property
    def pretrain_losses(self) -> List[float]:
        return self._pretrain_losses

    def get_encoder(self) -> MultiModalEncoder:
        """Return the pre-trained encoder for downstream use."""
        if not self._is_pretrained:
            logger.warning(
                "Encoder not pre-trained yet. Will train from scratch downstream."
            )
        return self.encoder
