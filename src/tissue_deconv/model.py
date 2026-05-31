#!/usr/bin/env python3
"""
Tissue Deconvolution Model — cfSort-Style DNN Head
=====================================================

Deep Neural Network for tissue-of-origin deconvolution from cfDNA
methylation beta values. Architecture follows the cfSort framework
(Li et al. 2023, PNAS):

- Input: methylation beta values at ~1000 tissue-discriminative CpG sites
- Hidden: [256, 128, 64] with BatchNorm1d + ReLU + Dropout
- Output: softmax over 29 tissues → tissue fractions summing to 1

Supports:
- Single model forward pass
- Ensemble of N models with different random seeds
- MPS (Apple Silicon), CUDA, and CPU
- ~500K parameters (lightweight, fits on any device)

Architecture
------------

::

    Input (N × n_cpg_features)
         │
         ▼
    Linear(1000→256) + BatchNorm1d + ReLU + Dropout(0.3)
         │
         ▼
    Linear(256→128) + BatchNorm1d + ReLU + Dropout(0.3)
         │
         ▼
    Linear(128→64) + BatchNorm1d + ReLU + Dropout(0.3)
         │
         ▼
    Linear(64→29) + Softmax(dim=1)
         │
         ▼
    Tissue Fractions (N × 29)
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None
    nn = None
    F = None

from .config import TissueDeconvConfig, DEFAULT_CONFIG, SUPPORTED_TISSUES

logger = logging.getLogger(__name__)


# ── Loss functions ──────────────────────────────────────────────

def _soft_kl_divergence(
    pred: "torch.Tensor",
    target: "torch.Tensor",
    eps: float = 1e-7,
) -> "torch.Tensor":
    """
    Soft KL divergence loss for tissue fraction prediction.

    KL(target || pred) — asymmetric, penalizes missing active tissues
    more than predicting spurious fractions.

    Both input and target are assumed to be probability distributions
    (softmax outputs or normalized fractions).
    """
    target_safe = target + eps
    pred_safe = pred + eps
    return (target_safe * (target_safe.log() - pred_safe.log())).sum(dim=-1).mean()


def _sparsity_loss(pred: "torch.Tensor", target_sparsity: float = 0.03) -> "torch.Tensor":
    """
    L1 sparsity penalty: most tissues should have near-zero fraction.

    Encourages the model to predict sparse tissue distributions,
    matching the biological reality that cfDNA originates from
    a limited number of tissues.
    """
    return pred.abs().mean()


def _smoothness_loss(pred: "torch.Tensor") -> "torch.Tensor":
    """
    Entropy-based smoothness penalty: prevents overconfident predictions.

    High entropy → uniform distribution (undesirable).
    Low entropy → peaked distribution (may miss secondary tissues).
    Optimal: moderate entropy reflecting biological complexity.
    """
    eps = 1e-7
    entropy = -(pred * (pred + eps).log()).sum(dim=-1).mean()
    # Penalize both extremes: too low entropy (overconfident) and too high
    # Target entropy for ~3 active tissues out of 29: -ln(1/29) ≈ 3.37
    target_entropy = 2.5  # roughly 3-4 non-zero tissues
    return (entropy - target_entropy) ** 2


# ═══════════════════════════════════════════════════════════════════
# DNN Model
# ═══════════════════════════════════════════════════════════════════

class TissueDeconvolutionModel(nn.Module):
    """
    cfSort-style DNN for tissue fraction prediction from methylation.

    This is a lightweight feed-forward neural network (~500K params)
    designed to deconvolve cfDNA methylation signals into per-tissue
    cell death proportions.

    Parameters
    ----------
    config : TissueDeconvConfig
        Model configuration.
    seed : int
        Random seed for weight initialization (used in ensemble).
    """

    def __init__(
        self,
        config: Optional[TissueDeconvConfig] = None,
        seed: int = 42,
    ):
        if not _HAS_TORCH:
            raise ImportError(
                "PyTorch is required for TissueDeconvolutionModel. "
                "Install with: pip install torch"
            )
        super().__init__()

        self.config = config if config is not None else DEFAULT_CONFIG
        self.n_tissues = self.config.n_tissues
        self.n_cpg_features = self.config.n_cpg_features
        self.hidden_dims = self.config.hidden_dims
        self.dropout_rate = self.config.dropout
        self.seed = seed

        # Set seed for reproducible initialization
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)

        # ── Build layers ─────────────────────────────────────────
        layers = []
        input_dim = self.n_cpg_features

        for hidden_dim in self.hidden_dims:
            layers.extend([
                nn.Linear(input_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(inplace=True),
                nn.Dropout(self.dropout_rate),
            ])
            input_dim = hidden_dim

        self.fc_layers = nn.Sequential(*layers)
        self.output_layer = nn.Linear(input_dim, self.n_tissues)

        # Initialize weights
        self._init_weights()

    def _init_weights(self):
        """Kaiming initialization for ReLU activations."""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)

    def forward(
        self,
        x: "torch.Tensor",
        return_logits: bool = False,
    ) -> Union["torch.Tensor", Tuple["torch.Tensor", "torch.Tensor"]]:
        """
        Forward pass.

        Parameters
        ----------
        x : (batch_size, n_cpg_features) tensor
            Methylation beta values at tissue-discriminative CpGs.
        return_logits : bool
            If True, returns (fractions, logits) tuple.

        Returns
        -------
        fractions : (batch_size, n_tissues) tensor
            Predicted tissue fractions (softmax, sums to 1.0).
        logits : (batch_size, n_tissues) tensor (if return_logits=True)
            Raw logits before softmax.
        """
        features = self.fc_layers(x)
        logits = self.output_layer(features)
        fractions = F.softmax(logits, dim=-1)

        if return_logits:
            return fractions, logits
        return fractions

    def predict(
        self,
        x: "torch.Tensor",
    ) -> np.ndarray:
        """
        Predict tissue fractions and return as numpy array.

        Parameters
        ----------
        x : tensor or ndarray
            Input methylation data.

        Returns
        -------
        fractions : (batch_size, n_tissues) float32 array
        """
        self.eval()
        with torch.no_grad():
            if not isinstance(x, torch.Tensor):
                x = torch.as_tensor(x, dtype=torch.float32)
            # Ensure at least 2D for BatchNorm
            if x.ndim == 1:
                x = x.unsqueeze(0)
            x = x.to(next(self.parameters()).device)
            return self.forward(x).cpu().numpy()

    @property
    def device(self) -> torch.device:
        """Device the model is on."""
        return next(self.parameters()).device

    def count_parameters(self) -> int:
        """Count trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ═══════════════════════════════════════════════════════════════════
# Ensemble Model
# ═══════════════════════════════════════════════════════════════════

class TissueDeconvolutionEnsemble:
    """
    Ensemble of N TissueDeconvolutionModel instances.

    Each model is trained with a different random seed. Predictions
    are averaged for more robust tissue fraction estimates, reducing
    variance from random initialization.

    Parameters
    ----------
    config : TissueDeconvConfig
        Model configuration (n_ensemble controls ensemble size).
    base_seed : int
        Base seed; each ensemble member gets base_seed + i.
    """

    def __init__(
        self,
        config: Optional[TissueDeconvConfig] = None,
        base_seed: int = 42,
    ):
        if not _HAS_TORCH:
            raise ImportError(
                "PyTorch is required for TissueDeconvolutionEnsemble. "
                "Install with: pip install torch"
            )

        self.config = config if config is not None else DEFAULT_CONFIG
        self.base_seed = base_seed
        self.n_ensemble = self.config.n_ensemble

        self.models: List[TissueDeconvolutionModel] = []
        self._initialize_models()

    def _initialize_models(self):
        """Create ensemble members with different seeds."""
        self.models = [
            TissueDeconvolutionModel(
                config=self.config,
                seed=self.base_seed + i,
            )
            for i in range(self.n_ensemble)
        ]

    def to(self, device: Union[str, torch.device]):
        """Move all models to device."""
        for model in self.models:
            model.to(device)
        return self

    def train(self):
        """Set all models to training mode."""
        for model in self.models:
            model.train()

    def eval(self):
        """Set all models to evaluation mode."""
        for model in self.models:
            model.eval()

    def predict(
        self,
        x: Union[np.ndarray, "torch.Tensor"],
        return_std: bool = False,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Ensemble prediction: average across all models.

        Parameters
        ----------
        x : (batch_size, n_cpg_features)
            Input methylation beta values.
        return_std : bool
            If True, returns (mean, std) tuple per tissue.

        Returns
        -------
        fractions : (batch_size, n_tissues) float32 array
            Mean tissue fractions across ensemble.
        std : (batch_size, n_tissues) float32 array (if return_std)
            Standard deviation across ensemble members.
        """
        self.eval()
        preds = []

        if not isinstance(x, torch.Tensor):
            x = torch.as_tensor(x, dtype=torch.float32)

        for model in self.models:
            pred = model.predict(x)
            preds.append(pred)

        preds = np.stack(preds, axis=0)  # (n_ensemble, batch, n_tissues)
        mean_pred = preds.mean(axis=0)

        if return_std:
            std_pred = preds.std(axis=0)
            return mean_pred, std_pred
        return mean_pred

    def predict_sample(
        self,
        x: np.ndarray,
        return_top: int = 5,
    ) -> Dict[str, float]:
        """
        Predict tissue fractions for a single sample.

        Parameters
        ----------
        x : (n_cpg_features,) array
            Methylation beta values for one sample.
        return_top : int
            Number of top tissues to include in output.

        Returns
        -------
        fractions : dict
            {tissue_name: fraction} for top tissues + '_entropy', '_n_active'.
        """
        if x.ndim == 1:
            x = x.reshape(1, -1)

        mean_pred, std_pred = self.predict(x, return_std=True)
        fractions = mean_pred[0]  # (n_tissues,)

        # Get top tissues
        top_indices = np.argsort(fractions)[::-1][:return_top]
        result = {}
        for idx in top_indices:
            tissue = SUPPORTED_TISSUES[idx]
            frac = float(fractions[idx])
            if frac > self.config.detection_limit:
                result[tissue] = frac

        # Metadata
        eps = 1e-7
        entropy = -np.sum(fractions * np.log(fractions + eps))
        max_entropy = np.log(self.config.n_tissues)
        result["_entropy"] = float(entropy)
        result["_n_active"] = int(np.sum(fractions > self.config.detection_limit))
        result["_top_tissue"] = SUPPORTED_TISSUES[int(np.argmax(fractions))]
        result["_top_fraction"] = float(np.max(fractions))

        return result

    def state_dicts(self) -> List[Dict]:
        """Get state dicts for all ensemble members."""
        return [m.state_dict() for m in self.models]

    def load_state_dicts(self, state_dicts: List[Dict]):
        """Load state dicts for all ensemble members."""
        if len(state_dicts) != self.n_ensemble:
            raise ValueError(
                f"Expected {self.n_ensemble} state dicts, got {len(state_dicts)}"
            )
        for model, sd in zip(self.models, state_dicts):
            model.load_state_dict(sd)

    @property
    def device(self) -> torch.device:
        """Device of the first ensemble member."""
        return self.models[0].device if self.models else torch.device("cpu")


# ── Loss computer ────────────────────────────────────────────────

class DeconvLoss(nn.Module):
    """
    Composite loss for tissue deconvolution training.

    Loss = KL_divergence + λ_sparsity * L1 + λ_smooth * EntropyReg

    The KL divergence ensures predicted fractions match the true
    tissue composition. The sparsity term encourages biologically
    realistic few-tissue solutions. The smoothness term prevents
    degenerate overconfident predictions.
    """

    def __init__(self, config: Optional[TissueDeconvConfig] = None):
        super().__init__()
        self.config = config if config is not None else DEFAULT_CONFIG
        self.lambda_sparsity = self.config.lambda_sparsity
        self.lambda_smooth = self.config.lambda_smooth

    def forward(
        self,
        pred: "torch.Tensor",
        target: "torch.Tensor",
    ) -> Tuple["torch.Tensor", Dict[str, float]]:
        """
        Compute composite loss.

        Parameters
        ----------
        pred : (batch_size, n_tissues) tensor
            Predicted tissue fractions (should be softmax output).
        target : (batch_size, n_tissues) tensor
            Ground truth tissue fractions.

        Returns
        -------
        total_loss : scalar tensor
        components : dict
            Individual loss components for logging.
        """
        kl_loss = _soft_kl_divergence(pred, target)
        sp_loss = _sparsity_loss(pred)
        sm_loss = _smoothness_loss(pred)

        total = kl_loss + self.lambda_sparsity * sp_loss + self.lambda_smooth * sm_loss

        components = {
            "kl_loss": float(kl_loss.detach().cpu()),
            "sparsity_loss": float(sp_loss.detach().cpu()),
            "smoothness_loss": float(sm_loss.detach().cpu()),
            "total_loss": float(total.detach().cpu()),
        }

        return total, components
