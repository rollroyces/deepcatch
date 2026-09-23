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


# ── Module-level wrappers for per-modality standardization ────────
# Exposed so callers (e.g. ``scripts/pretrain_production_finaledb.py``)
# can fit stats once on the full cohort and re-apply them at inference
# time without instantiating a full ``FoundationPretrainer``.


def _fit_modality_stats(
    modalities: Dict[str, np.ndarray],
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Fit per-modality median and MAD/std on the TRAIN fold.

    Per-modality standardization is the standard practice for
    self-supervised multi-modal pretraining (and what the
    ``cfdna-early-detection-validation`` skill recommends: "scale
    features before pretraining"). Without it, modalities with raw
    count scales (e.g. WPS in the thousands) dominate the MSE
    reconstruction loss and the encoder never actually converges.

    Parameters
    ----------
    modalities : dict[str, (n, dim_i) ndarray]
        TRAIN-fold modalities. Caller passes only training rows so
        the stats encode the training distribution, not the test
        distribution (no leakage).

    Returns
    -------
    stats : dict[str, (median, scale)]
        ``median`` is shape ``(dim_i,)``; ``scale`` is shape
        ``(dim_i,)``. ``scale`` is the MAD scaled to a robust std
        (``1.4826 * MAD``) when the MAD is nonzero, falling back
        to the column std when MAD is zero (constant features).
    """
    stats: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for name, arr in modalities.items():
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        median = np.median(arr, axis=0)
        mad = np.median(np.abs(arr - median), axis=0)
        # MAD -> robust std (Gaussian-consistent factor 1.4826).
        scale = 1.4826 * mad
        # Fallback: if MAD is zero (degenerate / constant feature)
        # use std so we still divide by a nonzero scale. If both
        # are zero (truly constant), use 1 to leave the feature
        # untouched rather than producing NaN.
        std = arr.std(axis=0)
        scale = np.where(scale > 0, scale, std)
        scale = np.where(scale > 0, scale, 1.0)
        stats[name] = (median.astype(np.float32),
                       scale.astype(np.float32))
    return stats


def _apply_modality_standardization(
    modalities: Dict[str, np.ndarray],
    stats: Dict[str, Tuple[np.ndarray, np.ndarray]],
) -> Dict[str, np.ndarray]:
    """Apply a fitted (median, scale) per modality.

    Standardizes to ``(x - median) / scale`` per feature column.
    Pure-function: returns a new dict with float32 arrays.
    """
    standardized: Dict[str, np.ndarray] = {}
    for name, arr in modalities.items():
        median, scale = stats[name]
        # Broadcast over rows: stats are shape (dim_i,)
        standardized[name] = ((arr - median) / scale).astype(np.float32)
    return standardized


class FoundationPretrainer:
    """
    Self-supervised pre-trainer for the DeepCatch Foundation Model.

    Handles 3-phase training with checkpointing and progress tracking.

    IMPORTANT — Real-data vs synthetic (BUG FIX 2026-09-23):
        Earlier versions of this class instantiated a
        ``MultiModalDataGenerator`` unconditionally and called
        ``self.data_generator.generate_dataset(...)`` from every phase,
        which silently bypassed any real modalities the caller had
        available. The result: a "pretrained" checkpoint that had
        never seen real data — only the synthetic generator's output.

        The constructor now accepts a real ``modalities`` dict and a
        ``use_real_modalities`` flag (default ``True``). When set, the
        per-phase methods draw mini-batches from the real cohort
        rather than calling the synthetic generator. Passing
        ``use_real_modalities=False`` keeps the synthetic path so
        unit tests / CI without real data still pass.

        Per-phase ``modalities`` / ``use_real_modalities`` arguments
        override the constructor-level values, so the all-synthetic
        convenience path used by ``test_integration.py`` is
        preserved.

    Parameters
    ----------
    config : FoundationConfig
        Model and training configuration.
    device : str, optional
        Compute device (default from config).
    verbose : bool
        Print progress during training.
    modalities : dict[str, np.ndarray], optional
        Real multi-modal cohort keyed by ``MODALITY_NAMES``. Each value
        must be a 2-D array of shape ``(n_samples, dim_i)``. When
        provided together with ``use_real_modalities=True`` (default),
        phases draw mini-batches from this cohort instead of calling
        ``self.data_generator.generate_dataset``.
    use_real_modalities : bool, default ``True``
        Default behaviour for phase methods. If no ``modalities``
        were passed at construction, this flag is treated as
        ``False`` (synthetic fallback) so legacy callers don't break.
    """

    def __init__(
        self,
        config: Optional[FoundationConfig] = None,
        device: Optional[str] = None,
        verbose: bool = False,
        modalities: Optional[Dict[str, np.ndarray]] = None,
        use_real_modalities: bool = True,
    ):
        self.config = config if config is not None else DEFAULT_CONFIG
        self.device = device or self.config.device
        self.verbose = verbose

        # Real-data defaults — constructor-level. Per-phase flags can
        # override these. We track the user's explicit choice so a
        # later "wants real but has none" call site can warn them.
        self._default_modalities = modalities
        self._user_wants_real = bool(use_real_modalities)
        self._default_use_real_modalities = bool(
            use_real_modalities and modalities is not None
        )

        # Init models
        self.encoder = MultiModalEncoder(self.config)
        self.pretrain_head = PretrainHead(self.config)
        self.contrastive_head = ContrastiveHead(self.config)

        # Move to device
        self.encoder.to(self.device)
        self.pretrain_head.to(self.device)
        self.contrastive_head.to(self.device)

        # Data generator — kept for the synthetic fallback path only.
        self.data_generator = MultiModalDataGenerator(seed=self.config.seed)

        # Training state
        self._is_pretrained = False
        self._pretrain_losses: List[float] = []
        self._phase_completed: List[str] = []

    # ── Real-cohort helper ──────────────────────────────────────────

    def _resolve_modalities(
        self,
        n_samples: int,
        prefix: str,
        override_modalities: Optional[Dict[str, np.ndarray]],
        override_use_real: Optional[bool],
    ) -> Tuple[Dict[str, torch.Tensor], bool]:
        """Pick the modality source for a phase.

        Returns
        -------
        modalities : dict[str, torch.Tensor]
            Modalities moved onto the configured device.
        used_real : bool
            ``True`` if these modalities came from a real cohort
            (not the synthetic generator). Callers can use this for
            logging and for the regression test in
            ``test/test_pretrain_bug_fix.py``.
        """
        # The user's explicit intent (either per-call or at construction)
        # matters here even when the *effective* flag has been normalised
        # to False because no modalities were supplied. We need to know
        # whether they THOUGHT they were using real data so we can warn
        # them if we're silently falling back to synthetic.
        use_real_intent = (
            override_use_real
            if override_use_real is not None
            else self._user_wants_real
        )
        use_real = (
            override_use_real
            if override_use_real is not None
            else self._default_use_real_modalities
        )
        mod_dict = override_modalities or self._default_modalities

        # If the caller wanted real data but didn't supply any, fall
        # back to synthetic — this preserves the pre-fix behaviour for
        # callers that forgot to pass a cohort. Warn so they don't
        # silently train on synthetic when they thought they were
        # training on real.
        if use_real_intent and mod_dict is None:
            logger.warning(
                "use_real_modalities=True but no modalities supplied; "
                "falling back to synthetic data for this phase."
            )
            use_real = False

        if use_real and mod_dict is not None:
            return self._modalities_to_tensors(mod_dict), True

        # Synthetic fallback.
        modalities_np, _ = self.data_generator.generate_dataset(
            n_samples=n_samples,
            prefix=prefix,
        )
        return self._modalities_to_tensors(modalities_np), False

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
        modalities: Optional[Dict[str, np.ndarray]] = None,
        use_real_modalities: Optional[bool] = None,
    ) -> List[float]:
        """
        Phase 1: Masked Modality Prediction.

        Train the encoder + pretrain head to predict masked modality
        features from unmasked context.

        Parameters
        ----------
        n_samples : int
            Number of samples to draw for training. When
            ``use_real_modalities=True`` this is silently capped at
            ``len(modalities)`` because we cannot synthesize new
            real samples.
        batch_size : int, optional
            Batch size (default from config).
        n_epochs : int
            Number of training epochs.
        lr : float, optional
            Learning rate (default from config).
        modalities : dict[str, np.ndarray], optional
            Per-phase override of the constructor-level real cohort.
        use_real_modalities : bool, optional
            Per-phase override of the constructor-level flag.
            ``None`` (default) means "use the constructor default".

        Returns
        -------
        losses : list of float
            Epoch losses.
        """
        batch_size = batch_size or self.config.batch_size
        lr = lr or self.config.pretrain_lr

        modalities_t, used_real = self._resolve_modalities(
            n_samples=n_samples,
            prefix="phase1",
            override_modalities=modalities,
            override_use_real=use_real_modalities,
        )
        effective_n = next(iter(modalities_t.values())).shape[0]
        source = "real" if used_real else "synthetic"
        self._log(
            f"Phase 1: MMP — {effective_n} samples ({source}), "
            f"{n_epochs} epochs"
        )

        # Per-modality standardization (median / MAD-robust-std) on the
        # TRAIN fold. Without this the MSE on modalities with raw count
        # scales (e.g. WPS in the thousands) dominates the loss and the
        # encoder never actually converges — see the docstring on
        # ``_fit_modality_stats``. We always run it: the synthetic
        # generator's healthy-range scales differ across modalities, so
        # the same pathology can appear there too.
        modalities_np = {
            name: tensor.detach().cpu().numpy()
            for name, tensor in modalities_t.items()
        }
        stats = _fit_modality_stats(modalities_np)
        modalities_np_std = _apply_modality_standardization(
            modalities_np, stats,
        )
        modalities_t = self._modalities_to_tensors(modalities_np_std)

        # Optimizer
        params = list(self.encoder.parameters()) + list(self.pretrain_head.parameters())
        optimizer = torch.optim.AdamW(params, lr=lr)

        losses = []
        n_batches = max(1, effective_n // batch_size)

        for epoch in range(n_epochs):
            epoch_loss = 0.0

            for batch_idx in range(n_batches):
                start = (batch_idx * batch_size) % effective_n
                end = min(start + batch_size, effective_n)
                batch_size_actual = end - start

                # Extract batch
                batch_modalities = {
                    k: v[start:end] for k, v in modalities_t.items()
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
        modalities: Optional[Dict[str, np.ndarray]] = None,
        use_real_modalities: Optional[bool] = None,
    ) -> List[float]:
        """
        Phase 2: Cross-modal contrastive learning.

        Train encoder + contrastive head to pull same-sample modalities
        together and push different-sample ones apart.

        Parameters
        ----------
        n_samples : int
            Number of samples to draw for training. When
            ``use_real_modalities=True`` this is capped at
            ``len(modalities)``.
        batch_size : int, optional
            Batch size.
        n_epochs : int
            Number of epochs.
        lr : float, optional
            Learning rate.
        modalities : dict[str, np.ndarray], optional
            Per-phase override of the constructor-level real cohort.
        use_real_modalities : bool, optional
            Per-phase override of the constructor-level flag.

        Returns
        -------
        losses : list of float
            Epoch losses.
        """
        batch_size = batch_size or self.config.batch_size
        lr = lr or self.config.pretrain_lr

        modalities_t, used_real = self._resolve_modalities(
            n_samples=n_samples,
            prefix="phase2",
            override_modalities=modalities,
            override_use_real=use_real_modalities,
        )
        effective_n = next(iter(modalities_t.values())).shape[0]
        source = "real" if used_real else "synthetic"
        self._log(
            f"Phase 2: Contrastive — {effective_n} samples ({source}), "
            f"{n_epochs} epochs"
        )

        # Per-modality standardization on the TRAIN fold — same rationale
        # as Phase 1 (see _fit_modality_stats docstring). Putting all
        # modalities on a comparable scale makes the joint embedding
        # space treat them symmetrically rather than letting the
        # largest-scale modality dominate the projection norms.
        modalities_np = {
            name: tensor.detach().cpu().numpy()
            for name, tensor in modalities_t.items()
        }
        stats = _fit_modality_stats(modalities_np)
        modalities_np_std = _apply_modality_standardization(
            modalities_np, stats,
        )
        modalities_t = self._modalities_to_tensors(modalities_np_std)

        # Optimizer
        params = list(self.encoder.parameters()) + list(self.contrastive_head.parameters())
        optimizer = torch.optim.AdamW(params, lr=lr)

        losses = []
        n_batches = max(1, effective_n // batch_size)

        for epoch in range(n_epochs):
            epoch_loss = 0.0
            n_batches_done = 0

            for batch_idx in range(n_batches):
                start = (batch_idx * batch_size) % effective_n
                end = min(start + batch_size, effective_n)
                batch_size_actual = end - start

                # Need at least 2 samples for contrastive
                if batch_size_actual < 2:
                    continue

                batch_modalities = {
                    k: v[start:end] for k, v in modalities_t.items()
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
        modalities: Optional[Dict[str, np.ndarray]] = None,
        use_real_modalities: Optional[bool] = None,
    ) -> List[float]:
        """
        Phase 3: Joint MMP + contrastive training.

        Parameters
        ----------
        n_samples : int
            Number of samples to draw for training.
        batch_size : int, optional
        n_epochs : int
        lr : float, optional
        modalities : dict[str, np.ndarray], optional
            Per-phase override of the constructor-level real cohort.
        use_real_modalities : bool, optional
            Per-phase override of the constructor-level flag.

        Returns
        -------
        losses : list of float
            Epoch losses (total combined).
        """
        batch_size = batch_size or self.config.batch_size
        lr = lr or (self.config.pretrain_lr * 0.5)

        modalities_t, used_real = self._resolve_modalities(
            n_samples=n_samples,
            prefix="phase3",
            override_modalities=modalities,
            override_use_real=use_real_modalities,
        )
        effective_n = next(iter(modalities_t.values())).shape[0]
        source = "real" if used_real else "synthetic"
        self._log(
            f"Phase 3: Joint — {effective_n} samples ({source}), "
            f"{n_epochs} epochs"
        )

        # Per-modality standardization on the TRAIN fold (see
        # _fit_modality_stats docstring for the rationale). Same
        # treatment as Phase 1 / Phase 2.
        modalities_np = {
            name: tensor.detach().cpu().numpy()
            for name, tensor in modalities_t.items()
        }
        stats = _fit_modality_stats(modalities_np)
        modalities_np_std = _apply_modality_standardization(
            modalities_np, stats,
        )
        modalities_t = self._modalities_to_tensors(modalities_np_std)

        params = (
            list(self.encoder.parameters())
            + list(self.pretrain_head.parameters())
            + list(self.contrastive_head.parameters())
        )
        optimizer = torch.optim.AdamW(params, lr=lr)

        losses = []
        n_batches = max(1, effective_n // batch_size)

        for epoch in range(n_epochs):
            epoch_total_loss = 0.0
            n_batches_done = 0

            for batch_idx in range(n_batches):
                start = (batch_idx * batch_size) % effective_n
                end = min(start + batch_size, effective_n)
                batch_size_actual = end - start

                if batch_size_actual < 2:
                    continue

                batch_modalities = {
                    k: v[start:end] for k, v in modalities_t.items()
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
        modalities: Optional[Dict[str, np.ndarray]] = None,
        use_real_modalities: Optional[bool] = None,
    ) -> Dict[str, List[float]]:
        """
        Run all three pre-training phases.

        Parameters
        ----------
        n_samples : int
            Samples per phase. When ``use_real_modalities=True`` this
            is capped at ``len(modalities)``.
        p1_epochs : int
            Phase 1 epochs.
        p2_epochs : int
            Phase 2 epochs.
        p3_epochs : int
            Phase 3 epochs.
        batch_size : int, optional
        modalities : dict[str, np.ndarray], optional
            Per-call override of the real cohort. Defaults to the
            constructor-level ``modalities``.
        use_real_modalities : bool, optional
            Per-call override of the flag. ``None`` (default) means
            "use the constructor default".

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
            n_samples=n_samples, n_epochs=p1_epochs, batch_size=batch_size,
            modalities=modalities, use_real_modalities=use_real_modalities,
        )
        p2_losses = self.pretrain_phase2_contrastive(
            n_samples=n_samples, n_epochs=p2_epochs, batch_size=batch_size,
            modalities=modalities, use_real_modalities=use_real_modalities,
        )
        p3_losses = self.pretrain_phase3_joint(
            n_samples=n_samples, n_epochs=p3_epochs, batch_size=batch_size,
            modalities=modalities, use_real_modalities=use_real_modalities,
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
        self.save_checkpoint_with_stats(path, modality_stats=None)

    def save_checkpoint_with_stats(
        self,
        path: str,
        modality_stats: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None,
    ) -> None:
        """Save pre-training checkpoint + optional standardization stats.

        Parameters
        ----------
        path : str
            Output file path (.pt).
        modality_stats : dict[str, (median, scale)], optional
            Per-modality standardization stats. When supplied, the
            checkpoint also stores ``modality_stats`` so downstream
            ``FoundationDownstream`` can re-apply the same transform
            before forward passes — without this the encoder would see
            raw WPS-scale tissue values at inference time even though
            it was trained on standardized values.
        """
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        # Serialize stats with tolist() so torch.save handles JSON-like
        # native Python lists/ndarrays.
        serializable_stats = None
        if modality_stats is not None:
            serializable_stats = {
                name: {
                    "median": np.asarray(med).tolist(),
                    "scale": np.asarray(sc).tolist(),
                }
                for name, (med, sc) in modality_stats.items()
            }
        checkpoint = {
            "encoder_state_dict": self.encoder.state_dict(),
            "pretrain_head_state_dict": self.pretrain_head.state_dict(),
            "contrastive_head_state_dict": self.contrastive_head.state_dict(),
            "config": self.config.to_dict(),
            "is_pretrained": self._is_pretrained,
            "phase_completed": self._phase_completed,
            "pretrain_losses": self._pretrain_losses,
            "modality_stats": serializable_stats,
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
        # Per-modality standardization stats (added in v2). Stored as a
        # plain-dict serialization; downstream code that wants to apply
        # the same transform at inference time should read
        # ``checkpoint["modality_stats"]`` directly rather than from
        # the pretrainer instance, because the pretrainer's internal
        # numpy arrays are not preserved across processes.
        self._loaded_modality_stats = checkpoint.get("modality_stats", None)

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
