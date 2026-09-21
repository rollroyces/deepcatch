#!/usr/bin/env python3
"""
Foundation Model Downstream Fine-tuning
========================================

Drop-in replacement for CrossAttentionFusion that uses the pre-trained
foundation model encoder for cancer detection, TOO prediction, and
healthy aging screening.

API is fully compatible with CrossAttentionFusion:

.. code-block:: python

    # Before (CrossAttentionFusion)
    fusion = CrossAttentionFusion(n_modalities=6)
    fusion.fit(scores, labels)
    proba = fusion.predict_proba(scores)

    # After (FoundationDownstream, pre-trained)
    fusion = FoundationDownstream(pretrained=True)
    fusion.fit(modalities, labels)       # dict of modality arrays
    proba = fusion.predict_proba(modalities)

    # After (FoundationDownstream, from scratch)
    fusion = FoundationDownstream(pretrained=False)
    fusion.fit(modalities, labels)
    proba = fusion.predict_proba(modalities)
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import (
    FoundationConfig,
    DEFAULT_CONFIG,
    PROTOTYPE_CONFIG,
    MODALITY_DIMS,
    MODALITY_NAMES,
)
from .model import MultiModalEncoder
from .pretrain import FoundationPretrainer


# ── Modality schema fingerprint ─────────────────────────────────────
# Per-modality expected_dim + per-row median range + allow_negative.
# A modality whose per-row median falls outside the expected range is
# almost certainly garbage input (random Gaussian, all-zeros, wrong
# scale). Catching this BEFORE training prevents the silent-failure
# class where the model produces plausible numbers for the wrong
# reasons. Used by FoundationDownstream.fit(validate_schema=True).

@dataclass(frozen=True)
class ModalitySchema:
    """Schema-fingerprint for a single modality."""
    name: str
    expected_dim: int
    min_median: float
    max_median: float
    allow_negative: bool = True


# Default medians derived from MultiModalDataGenerator.HEALTHY_RANGES
# (frag_basic=0.3, frag_enhanced=0.0, cnv=0.0, sero=0.5, gnn=0.0,
# tissue=0.0) plus a wide tolerance for real-data variation. These
# ranges are intentionally generous — the goal is to catch RANDOM
# input (e.g. accidental np.random.randn) rather than gate real
# biological variation.
MODALITY_SCHEMAS: Dict[str, ModalitySchema] = {
    "frag_basic": ModalitySchema(
        "frag_basic", MODALITY_DIMS["frag_basic"],
        min_median=-5.0, max_median=10.0,
        allow_negative=False,
    ),
    "frag_enhanced": ModalitySchema(
        "frag_enhanced", MODALITY_DIMS["frag_enhanced"],
        min_median=-5.0, max_median=5.0,
        allow_negative=True,
    ),
    "cnv": ModalitySchema(
        "cnv", MODALITY_DIMS["cnv"],
        min_median=-3.0, max_median=3.0,
        allow_negative=True,
    ),
    "sero": ModalitySchema(
        "sero", MODALITY_DIMS["sero"],
        min_median=0.0, max_median=1000.0,
        allow_negative=False,
    ),
    "gnn": ModalitySchema(
        "gnn", MODALITY_DIMS["gnn"],
        min_median=-3.0, max_median=3.0,
        allow_negative=True,
    ),
    "tissue": ModalitySchema(
        "tissue", MODALITY_DIMS["tissue"],
        min_median=0.0, max_median=1.0,
        allow_negative=False,
    ),
}

logger = logging.getLogger(__name__)


class ClassificationHead(nn.Module):
    """Simple MLP classification head on top of joint embeddings."""

    def __init__(
        self,
        embed_dim: int,
        n_modalities: int,
        n_classes: int = 2,
        hidden_dim: int = 64,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.cls = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, n_classes),
        )

    def forward(self, joint_embedding):
        # Mean pool across modalities
        x = joint_embedding.mean(dim=1)  # (B, embed_dim)
        return self.cls(x)


class FoundationDownstream:
    """
    Downstream classifier using the pre-trained foundation model encoder.

    Drop-in replacement for CrossAttentionFusion with two modes:
    - **Pre-trained**: Load a pre-trained encoder checkpoint, fine-tune.
    - **From scratch**: Train encoder + classifier end-to-end.

    Parameters
    ----------
    config : FoundationConfig, optional
        Model configuration.
    pretrained : bool
        If True, use prototype (lightweight) config and expect
        pre-trained weights. If False, train from scratch.
    checkpoint_path : str, optional
        Path to pre-trained checkpoint (.pt).
    freeze_encoder : bool
        If True, freeze encoder weights during fine-tuning.
    device : str, optional
        Compute device.
    """

    def __init__(
        self,
        config: Optional[FoundationConfig] = None,
        pretrained: bool = True,
        checkpoint_path: Optional[str] = None,
        freeze_encoder: bool = False,
        device: Optional[str] = None,
        loss: str = "ce",
        alpha_pos: float = 20.0,
        gamma: float = 2.0,
    ):
        """FoundationDownstream constructor.

        Parameters
        ----------
        loss : {"ce", "sens_at_spec"}
            Loss function for binary classification. "ce" (default)
            preserves all existing benchmark numbers.
            "sens_at_spec" uses focal-modulated BCE with
            ``alpha_pos`` rebalancing for ultra-low VAF cohorts.
            Multi-class always uses CE.
        alpha_pos : float
            Positive-class weight for focal-BCE. Ignored when
            ``loss="ce"``. Default 20.0 is the documented starting
            point for 0.1% VAF cohorts.
        gamma : float
            Focal modulation exponent. Ignored when ``loss="ce"``.
        """
        if loss not in ("ce", "sens_at_spec"):
            raise ValueError(
                f"loss must be 'ce' or 'sens_at_spec', got {loss!r}"
            )
        self.loss = loss
        self.alpha_pos = alpha_pos
        self.gamma = gamma
        self.config = config if config is not None else (
            PROTOTYPE_CONFIG if pretrained else DEFAULT_CONFIG
        )
        self.device = device or self.config.device
        self.freeze_encoder = freeze_encoder
        self._pretrained = pretrained
        self._fitted = False
        self._n_classes = 2

        # Build encoder
        self.encoder = MultiModalEncoder(self.config)

        # Load pre-trained weights if available
        if pretrained and checkpoint_path and os.path.exists(checkpoint_path):
            self._load_pretrained_encoder(checkpoint_path)

        # Classification head (built after first fit or explicitly)
        self.classifier: Optional[ClassificationHead] = None

        # Optimizer state
        self._optimizer = None
        self._loss_history: List[float] = []

        # Move to device
        self.encoder.to(self.device)

        if freeze_encoder:
            for param in self.encoder.parameters():
                param.requires_grad = False

    def _load_pretrained_encoder(self, checkpoint_path: str):
        """Load pre-trained encoder weights from checkpoint."""
        try:
            checkpoint = torch.load(checkpoint_path, map_location=self.device)
            self.encoder.load_state_dict(checkpoint["encoder_state_dict"])
            logger.info(f"Loaded pre-trained encoder from {checkpoint_path}")
        except (KeyError, RuntimeError, FileNotFoundError) as e:
            logger.warning(
                f"Failed to load checkpoint {checkpoint_path}: {e}. "
                f"Will train from scratch."
            )
            self._pretrained = False

    def _build_classifier(self, n_classes: int):
        """Build classification head."""
        self._n_classes = n_classes
        self.classifier = ClassificationHead(
            embed_dim=self.config.embed_dim,
            n_modalities=self.config.n_modalities,
            n_classes=n_classes,
            dropout=self.config.dropout,
        ).to(self.device)

    def _modalities_to_tensors(
        self,
        modalities: Dict[str, np.ndarray],
    ) -> Dict[str, torch.Tensor]:
        """Convert numpy modality dict to torch tensors."""
        return {
            k: torch.from_numpy(v.astype(np.float32)).to(self.device)
            for k, v in modalities.items()
        }

    def _validate_modalities(self, modalities: Dict[str, np.ndarray]):
        """Check that input modalities match expected dimensions."""
        for name in MODALITY_NAMES:
            if name not in modalities:
                raise KeyError(
                    f"Missing modality '{name}'. "
                    f"Expected keys: {list(MODALITY_NAMES)}"
                )
            arr = modalities[name]
            expected_dim = MODALITY_DIMS[name]
            if arr.ndim == 1:
                if len(arr) != expected_dim:
                    raise ValueError(
                        f"Modality '{name}' has dim {len(arr)}, "
                        f"expected {expected_dim}"
                    )
            elif arr.ndim == 2:
                if arr.shape[1] != expected_dim:
                    raise ValueError(
                        f"Modality '{name}' has dim {arr.shape[1]}, "
                        f"expected {expected_dim}"
                    )

    def _validate_modality_schema(
        self,
        modalities: Dict[str, np.ndarray],
    ) -> List[str]:
        """Schema-fingerprint check: catch silent fallbacks to random data.

        Per-modality ``expected_dim``, ``min_median``, ``max_median``,
        and ``allow_negative`` bounds are defined in
        ``MODALITY_SCHEMAS`` (this module). A modality whose per-row
        median falls outside its expected range is almost certainly
        garbage input (random Gaussian, zeros, wrong scale). This
        catches silent fallbacks before the model trains on noise.

        Returns a list of warning strings for non-fatal anomalies
        (e.g. a sample whose median is near the boundary). Raises
        ValueError when an input clearly does not match the schema
        (e.g. expected non-negative but contains negatives).

        Only called when ``fit(..., validate_schema=True)`` is set.
        """
        warnings_list: List[str] = []
        for name in MODALITY_NAMES:
            arr = modalities[name]
            schema = MODALITY_SCHEMAS.get(name)
            if schema is None:
                # No schema registered — skip silently.
                continue
            # Compute per-row median for the 2-D batched case, else
            # median of the 1-D vector.
            if arr.ndim == 2:
                row_medians = np.median(arr, axis=1)
            else:
                row_medians = np.array([float(np.median(arr))])
            med = float(np.median(row_medians))
            if not (schema.min_median <= med <= schema.max_median):
                raise ValueError(
                    f"Modality '{name}' median {med:.4g} outside expected "
                    f"range [{schema.min_median}, {schema.max_median}]. "
                    f"This usually means the modality is random noise, "
                    f"all-zeros, or wrong-scale. Pass "
                    f"validate_schema=False to skip this check."
                )
            if not schema.allow_negative and (arr < 0).any():
                raise ValueError(
                    f"Modality '{name}' contains negative values but the "
                    f"schema requires non-negative inputs."
                )
            # Soft warning when the median is in the outer 10% of the
            # allowed range — caller may want to inspect.
            span = schema.max_median - schema.min_median
            if span > 0:
                lo_warn = schema.min_median + 0.1 * span
                hi_warn = schema.max_median - 0.1 * span
                if not (lo_warn <= med <= hi_warn):
                    warnings_list.append(
                        f"Modality '{name}' median {med:.4g} is near the "
                        f"edge of expected range "
                        f"[{schema.min_median}, {schema.max_median}]."
                    )
        return warnings_list

    def fit(
        self,
        modalities: Dict[str, np.ndarray],
        labels: np.ndarray,
        n_epochs: int = 50,
        batch_size: int = 32,
        lr: Optional[float] = None,
        validation_split: float = 0.1,
        early_stopping: bool = False,
        patience: int = 10,
        verbose: bool = False,
        validate_schema: bool = False,
    ) -> "FoundationDownstream":
        """
        Fine-tune (or train from scratch) the foundation model.

        Parameters
        ----------
        modalities : dict[str, ndarray]
            Modality features. Each value is (n_samples, dim_i).
        labels : (n_samples,) array
            Class labels (0 = healthy, 1 = cancer, etc.).
        n_epochs : int
            Number of fine-tuning epochs.
        batch_size : int
            Mini-batch size.
        lr : float, optional
            Learning rate (default from config).
        validation_split : float
            Fraction of data for validation.
        early_stopping : bool
            Whether to use early stopping.
        patience : int
            Patience for early stopping.
        verbose : bool
            Print training progress.
        validate_schema : bool
            If True, run ``_validate_modality_schema`` before training
            and raise ``ValueError`` on out-of-range per-row medians
            or negative values in non-negative modalities. Default
            False to preserve existing behavior. Recommended True
            for any clinical / production training path.

        Returns
        -------
        self
        """
        self._validate_modalities(modalities)
        if validate_schema:
            schema_warnings = self._validate_modality_schema(modalities)
            for w in schema_warnings:
                logger.warning(w)

        lr = lr or self.config.finetune_lr

        # Determine number of classes
        unique_labels = np.unique(labels)
        n_classes = len(unique_labels)
        if n_classes < 2:
            raise ValueError(f"Need at least 2 classes, got {n_classes}")

        # Build classifier if needed
        if self.classifier is None or self._n_classes != n_classes:
            self._build_classifier(n_classes)

        # Convert to tensors
        modalities_t = self._modalities_to_tensors(modalities)
        labels_t = torch.from_numpy(labels.astype(np.int64)).to(self.device)

        # Stratified train/val split that preserves class balance.
        # The previous torch.randperm split could produce val sets
        # with no positives (when validation_split * n_positives < 1),
        # in which case F.cross_entropy on the val fold returns NaN
        # and the early-stopping patience increments forever. The
        # stratified split also gives more reliable val-loss estimates
        # on imbalanced cfDNA screening cohorts.
        n_samples = int(labels_t.shape[0])
        rng = np.random.default_rng(self.config.seed)
        labels_np = labels.astype(np.int64)
        idx_pos = np.where(labels_np == 1)[0]
        idx_neg = np.where(labels_np == 0)[0]
        rng.shuffle(idx_pos)
        rng.shuffle(idx_neg)
        n_val_pos = max(0, int(round(len(idx_pos) * validation_split)))
        n_val_neg = max(0, int(round(len(idx_neg) * validation_split)))
        # Guarantee at least one of each class when both classes exist.
        if n_classes >= 2 and len(idx_pos) >= 2 and n_val_pos == 0:
            n_val_pos = 1
        if n_classes >= 2 and len(idx_neg) >= 2 and n_val_neg == 0:
            n_val_neg = 1
        val_idx_np = np.concatenate(
            [idx_pos[:n_val_pos], idx_neg[:n_val_neg]]
        )
        train_idx_np = np.concatenate(
            [idx_pos[n_val_pos:], idx_neg[n_val_neg:]]
        )
        train_idx = torch.from_numpy(train_idx_np).to(self.device)
        val_idx = torch.from_numpy(val_idx_np).to(self.device)

        # Training mode
        self.encoder.train()
        self.classifier.train()

        # Optimizer
        params = list(self.encoder.parameters()) + list(self.classifier.parameters())
        self._optimizer = torch.optim.AdamW(
            params, lr=lr, weight_decay=1e-5
        )

        # Focal-BCE: only when binary and user requested it.
        focal_bce_fn = None
        if self.loss == "sens_at_spec" and n_classes == 2:
            from src.foundation.losses import focal_binary_cross_entropy
            focal_bce_fn = focal_binary_cross_entropy

        def _compute_batch_loss(logits: "torch.Tensor",
                                batch_labels: "torch.Tensor") -> "torch.Tensor":
            if focal_bce_fn is not None and n_classes == 2:
                if logits.dim() == 2 and logits.shape[-1] == 2:
                    bin_logit = logits[:, 1] - logits[:, 0]
                else:
                    bin_logit = logits.squeeze(-1)
                tgt = batch_labels.float() if batch_labels.dtype != torch.float32 else batch_labels
                return focal_bce_fn(
                    bin_logit, tgt,
                    alpha_pos=self.alpha_pos, alpha_neg=1.0,
                    gamma=self.gamma, reduction="mean",
                )
            return F.cross_entropy(logits, batch_labels)

        self._loss_history = []
        best_val_loss = float("inf")
        best_state = None
        patience_counter = 0
        nan_guard = 0  # NaN/Inf loss counter — separate from patience

        for epoch in range(n_epochs):
            # Shuffle training indices
            shuffle_idx = train_idx[torch.randperm(len(train_idx))]

            # Mini-batch training
            epoch_loss = 0.0
            n_batches = max(1, len(shuffle_idx) // batch_size)

            for batch_start in range(0, len(shuffle_idx), batch_size):
                batch_idx = shuffle_idx[batch_start:batch_start + batch_size]

                batch_mod = {
                    k: v[batch_idx] for k, v in modalities_t.items()
                }
                batch_labels = labels_t[batch_idx]

                # Forward
                joint = self.encoder(batch_mod)  # (B, N, D)
                logits = self.classifier(joint)  # (B, n_classes)
                loss = _compute_batch_loss(logits, batch_labels)

                # NaN/Inf guard: skip this step (don't poison grads).
                # If we see too many in a row, abort training — the
                # current best state (or random init) is preserved.
                if torch.isnan(loss) or torch.isinf(loss):
                    nan_guard += 1
                    if nan_guard >= patience:
                        logger.warning(
                            "FoundationDownstream: %d consecutive NaN/Inf "
                            "training losses; aborting at epoch %d.",
                            nan_guard, epoch + 1,
                        )
                        break
                    continue

                # Backward
                self._optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(params, max_norm=1.0)
                self._optimizer.step()

                epoch_loss += loss.item()
                nan_guard = 0  # reset on a successful step

            if nan_guard >= patience:
                break

            avg_loss = epoch_loss / max(1, n_batches)
            self._loss_history.append(avg_loss)

            # Validation
            if early_stopping and len(val_idx) > 0:
                self.encoder.eval()
                self.classifier.eval()
                with torch.no_grad():
                    val_mod = {k: v[val_idx] for k, v in modalities_t.items()}
                    val_joint = self.encoder(val_mod)
                    val_logits = self.classifier(val_joint)
                    val_loss = _compute_batch_loss(val_logits, labels_t[val_idx])

                self.encoder.train()
                self.classifier.train()

                # NaN/Inf val-loss guard: don't update best_state on
                # garbage and don't increment patience counter so the
                # loop terminates normally.
                if torch.isnan(val_loss) or torch.isinf(val_loss):
                    if verbose:
                        logger.warning(
                            "FoundationDownstream: NaN/Inf val loss at "
                            "epoch %d; skipping early-stopping update.",
                            epoch + 1,
                        )
                    continue

                if val_loss < best_val_loss - 1e-4:
                    best_val_loss = val_loss.item()
                    best_state = {
                        "encoder": {k: v.cpu().clone() for k, v in self.encoder.state_dict().items()},
                        "classifier": {k: v.cpu().clone() for k, v in self.classifier.state_dict().items()},
                    }
                    patience_counter = 0
                else:
                    patience_counter += 1
                    if patience_counter >= patience:
                        if verbose:
                            logger.info(
                                f"Early stopping at epoch {epoch+1}"
                            )
                        break

            if verbose and (epoch + 1) % 10 == 0:
                logger.info(
                    f"Epoch {epoch+1}/{n_epochs} — loss: {avg_loss:.6f}"
                )

        # Restore best state
        if best_state is not None:
            self.encoder.load_state_dict(best_state["encoder"])
            self.classifier.load_state_dict(best_state["classifier"])

        self._fitted = True
        return self

    def predict_proba(
        self,
        modalities: Dict[str, np.ndarray],
    ) -> np.ndarray:
        """
        Predict class probabilities.

        Parameters
        ----------
        modalities : dict[str, ndarray]
            Modality features.

        Returns
        -------
        proba : (n_samples, n_classes) ndarray
            Softmax probabilities. For binary, proba[:, 1] is cancer probability.
        """
        if not self._fitted:
            raise RuntimeError(
                "Model not fitted. Call fit() before predict_proba()."
            )

        self._validate_modalities(modalities)
        modalities_t = self._modalities_to_tensors(modalities)

        # Determine batch size
        first_arr = list(modalities.values())[0]
        n_samples = first_arr.shape[0]

        self.encoder.eval()
        self.classifier.eval()

        all_probs = []

        with torch.no_grad():
            batch_size = 256
            for start in range(0, n_samples, batch_size):
                end = min(start + batch_size, n_samples)
                batch_mod = {
                    k: v[start:end] for k, v in modalities_t.items()
                }
                joint = self.encoder(batch_mod)
                logits = self.classifier(joint)
                probs = F.softmax(logits, dim=-1)
                all_probs.append(probs.cpu().numpy())

        return np.concatenate(all_probs, axis=0)

    def predict(self, modalities: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Predict hard class labels.

        Parameters
        ----------
        modalities : dict[str, ndarray]

        Returns
        -------
        labels : (n_samples,) int64 array
        """
        proba = self.predict_proba(modalities)
        return np.argmax(proba, axis=1)

    def encode(
        self,
        modalities: Dict[str, np.ndarray],
    ) -> np.ndarray:
        """
        Get global embedding vector for downstream analysis.

        Parameters
        ----------
        modalities : dict[str, ndarray]

        Returns
        -------
        embedding : (n_samples, embed_dim) ndarray
        """
        self._validate_modalities(modalities)
        modalities_t = self._modalities_to_tensors(modalities)

        self.encoder.eval()
        with torch.no_grad():
            global_emb = self.encoder.encode_single(modalities_t)
        return global_emb.cpu().numpy()

    def save_checkpoint(self, path: str):
        """Save fine-tuned model checkpoint."""
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        checkpoint = {
            "encoder_state_dict": self.encoder.state_dict(),
            "classifier_state_dict": self.classifier.state_dict() if self.classifier else None,
            "config": self.config.to_dict(),
            "n_classes": self._n_classes,
            "fitted": self._fitted,
            "loss_history": self._loss_history,
        }
        torch.save(checkpoint, path)

    def load_checkpoint(self, path: str) -> bool:
        """Load fine-tuned model checkpoint."""
        if not os.path.exists(path):
            return False
        checkpoint = torch.load(path, map_location=self.device)
        self.encoder.load_state_dict(checkpoint["encoder_state_dict"])
        if checkpoint.get("classifier_state_dict"):
            n_classes = checkpoint.get("n_classes", 2)
            if self.classifier is None or self._n_classes != n_classes:
                self._build_classifier(n_classes)
            self.classifier.load_state_dict(checkpoint["classifier_state_dict"])
        self._n_classes = checkpoint.get("n_classes", 2)
        self._fitted = checkpoint.get("fitted", True)
        self._loss_history = checkpoint.get("loss_history", [])
        return True

    @property
    def is_fitted(self) -> bool:
        return self._fitted

    @property
    def loss_history(self) -> List[float]:
        return self._loss_history

    @property
    def num_params(self) -> int:
        return self.encoder.num_params + (
            sum(p.numel() for p in self.classifier.parameters())
            if self.classifier else 0
        )


# ── Compatibility wrapper for CrossAttentionFusion API ──────────


class FoundationCompatibilityWrapper:
    """
    Wrapper that provides exact CrossAttentionFusion API using
    the foundation model.

    Converts per-modality scalar scores (CrossAttentionFusion API)
    to full modality features and back.

    .. code-block:: python

        # CrossAttentionFusion API
        fusion = CrossAttentionFusion(n_modalities=6)
        fusion.fit(scores_list, labels)       # list of (n,) arrays
        proba = fusion.predict_proba(scores_list)

        # Foundation model (same API)
        from src.foundation import FoundationDownstream, FoundationCompatibilityWrapper
        fusion = FoundationCompatibilityWrapper(pretrained=True)
        fusion.fit(scores_list, labels)
        proba = fusion.predict_proba(scores_list)
    """

    def __init__(
        self,
        config: Optional[FoundationConfig] = None,
        pretrained: bool = True,
        checkpoint_path: Optional[str] = None,
        device: Optional[str] = None,
    ):
        self.foundation = FoundationDownstream(
            config=config,
            pretrained=pretrained,
            checkpoint_path=checkpoint_path,
            device=device,
        )
        self.n_modalities = 6
        self._fitted = False

    def _scores_to_modalities(
        self,
        modality_scores: List[np.ndarray],
    ) -> Dict[str, np.ndarray]:
        """
        Convert list of 1D scores to full modality dict.
        Each scalar score is expanded into the modality's expected dimension
        by replicating or using as the main feature.
        """
        modalities = {}
        for i, name in enumerate(MODALITY_NAMES):
            if i < len(modality_scores):
                scores = modality_scores[i]
                dim = MODALITY_DIMS[name]
                if scores.ndim == 1:
                    # Expand scalar score to expected dimension
                    # by repeating + small noise for variation
                    expanded = np.tile(scores[:, np.newaxis], (1, dim))
                    # Add small noise so features aren't identical
                    noise = np.random.randn(*expanded.shape) * 0.01
                    modalities[name] = (expanded + noise).astype(np.float32)
                else:
                    modalities[name] = scores.astype(np.float32)
            else:
                # Padding with zeros
                modalities[name] = np.zeros((1, MODALITY_DIMS[name]), dtype=np.float32)
        return modalities

    def fit(
        self,
        modality_scores: List[np.ndarray],
        labels: np.ndarray,
        **kwargs,
    ) -> "FoundationCompatibilityWrapper":
        """Fit wrapper (CrossAttentionFusion-compatible API)."""
        modalities = self._scores_to_modalities(modality_scores)
        self.foundation.fit(modalities, labels, **kwargs)
        self._fitted = True
        return self

    def predict_proba(
        self,
        modality_scores: List[np.ndarray],
    ) -> np.ndarray:
        """Predict (CrossAttentionFusion-compatible API)."""
        if not self._fitted:
            raise RuntimeError("Not fitted")
        modalities = self._scores_to_modalities(modality_scores)
        return self.foundation.predict_proba(modalities)
