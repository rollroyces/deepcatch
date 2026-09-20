"""
Custom PyTorch losses and biologically-motivated loss helpers for cfDNA
detection.

The default cross-entropy loss optimizes accuracy at the 0.5 decision
threshold and gives ~50/50 sensitivity at 99% specificity. For ultra-low
VAF screening cohorts, the metric that matters is **sensitivity at a
fixed specificity ≥ 99%**. The losses here re-weight the gradient to
emphasize the rare-positive operating point.

Components
----------

1. ``SensAtSpecLoss`` — focal-modulated binary cross-entropy with an
   explicit negative-class down-weight. ``alpha_pos`` controls the
   relative gradient from positive vs negative samples; values in
   ``[10, 50]`` are typical starting points for 0.1% VAF cohorts.
2. ``focal_binary_cross_entropy`` — the underlying focal-modulated BCE
   for users who want to compose it with other reductions.
3. ``balanced_cross_entropy`` — inverse-frequency weighted CE for
   multi-class without focal modulation.
4. ``calibration_loss`` — expected-calibration-error (ECE) gradient
   surrogate for joint training.

References
----------
- Lin et al. (2017). "Focal Loss for Dense Object Detection." ICCV.
- Mukhoti et al. (2020). "Calibrating Deep Neural Networks using
  Focal Loss." NeurIPS.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    _HAS_TORCH = True
except ImportError:  # pragma: no cover
    _HAS_TORCH = False
    torch = None  # type: ignore
    nn = None  # type: ignore
    F = None  # type: ignore


def focal_binary_cross_entropy(
    logits: "torch.Tensor",
    targets: "torch.Tensor",
    alpha_pos: float = 1.0,
    alpha_neg: float = 1.0,
    gamma: float = 2.0,
    reduction: str = "mean",
) -> "torch.Tensor":
    """Focal binary cross-entropy with separate per-class weights.

    Focal loss (Lin 2017) down-weights easy-to-classify samples so
    that hard, near-decision-boundary examples dominate the gradient.
    Combined with explicit class weights (``alpha_pos``, ``alpha_neg``)
    this is the recommended loss for ultra-low VAF detection where
    the model otherwise collapses to "predict all negative".

    Parameters
    ----------
    logits : (N,) or (N, 1) Tensor
        Raw logits (NOT probabilities).
    targets : (N,) Tensor of float in {0, 1}
        Binary labels.
    alpha_pos : float
        Multiplier applied to the positive-class term.
    alpha_neg : float
        Multiplier applied to the negative-class term.
    gamma : float
        Focal modulation exponent. 0 = no focal, 2 = standard focal.
    reduction : str
        "mean" or "sum" or "none".

    Returns
    -------
    loss : scalar Tensor
    """
    if logits.dim() == 2 and logits.shape[-1] == 1:
        logits = logits.squeeze(-1)
    if targets.dim() == 2 and targets.shape[-1] == 1:
        targets = targets.squeeze(-1)

    p = torch.sigmoid(logits)
    # Per-sample class weight
    alpha_t = torch.where(
        targets > 0.5,
        torch.full_like(p, alpha_pos),
        torch.full_like(p, alpha_neg),
    )
    # True-class probability for the focal modulator.
    p_t = torch.where(targets > 0.5, p, 1.0 - p)
    focal = (1.0 - p_t).clamp_min(1e-8).pow(gamma)

    bce = F.binary_cross_entropy_with_logits(
        logits, targets.float(), reduction="none"
    )
    loss = alpha_t * focal * bce
    if reduction == "mean":
        return loss.mean()
    if reduction == "sum":
        return loss.sum()
    return loss


class SensAtSpecLoss(nn.Module):
    """Binary focal loss for ultra-low VAF detection.

    Designed for cfDNA screening cohorts where the metric of interest
    is sensitivity at ≥99% specificity, not overall accuracy. The
    default ``alpha_pos`` rebalances the gradient so the model has
    incentive to push the rare-positive tail above the 99%
    specificity threshold.

    Parameters
    ----------
    alpha_pos : float
        Positive-class weight. Set ``alpha_pos = 1 / prev`` for the
        textbook inverse-prevalence weighting (e.g. ``alpha_pos = 250``
        for 0.4% screening prevalence). Values in ``[10, 50]`` are
        good starting points for 0.1% VAF analytical cohorts.
    alpha_neg : float
        Negative-class weight (usually 1.0).
    gamma : float
        Focal modulation exponent. Higher = harder focus on
        misclassified samples. 2.0 is the standard.
    label_smoothing : float, optional
        Symmetric label smoothing for the binary targets. Helps
        calibrate the model at high specificity operating points.

    Examples
    --------
    >>> loss_fn = SensAtSpecLoss(alpha_pos=20.0, gamma=2.0)
    >>> logits = torch.randn(64)  # 64 binary logits
    >>> y = torch.randint(0, 2, (64,)).float()
    >>> loss = loss_fn(logits, y)
    """

    def __init__(
        self,
        alpha_pos: float = 20.0,
        alpha_neg: float = 1.0,
        gamma: float = 2.0,
        label_smoothing: float = 0.0,
    ):
        super().__init__()
        if alpha_pos <= 0 or alpha_neg <= 0:
            raise ValueError("alpha_pos and alpha_neg must be positive")
        if gamma < 0:
            raise ValueError("gamma must be non-negative")
        if not 0.0 <= label_smoothing < 0.5:
            raise ValueError(
                "label_smoothing must be in [0, 0.5)"
            )
        self.alpha_pos = alpha_pos
        self.alpha_neg = alpha_neg
        self.gamma = gamma
        self.label_smoothing = label_smoothing

    def forward(
        self,
        logits: "torch.Tensor",
        targets: "torch.Tensor",
    ) -> "torch.Tensor":
        if self.label_smoothing > 0:
            targets = targets * (1.0 - self.label_smoothing) + 0.5 * self.label_smoothing
        return focal_binary_cross_entropy(
            logits,
            targets,
            alpha_pos=self.alpha_pos,
            alpha_neg=self.alpha_neg,
            gamma=self.gamma,
            reduction="mean",
        )


class BalancedCrossEntropy(nn.Module):
    """Inverse-frequency weighted multi-class cross-entropy.

    Use this for multi-class (tissue-of-origin, multi-cancer
    classification) where the default ``F.cross_entropy`` under-trains
    rare classes. Weights are computed from per-class inverse frequency
    with optional sqrt-dampening (Cui et al. 2019, "Class-Balanced
    Loss Based on Effective Number of Samples").

    Parameters
    ----------
    labels : array-like of int
        Per-sample class labels used to compute the weights. Pass the
        training-fold labels so weights reflect training distribution.
    beta : float
        Effective-number rebalancing exponent. 0 = inverse frequency,
        1 = sqrt-inverse (Cui 2019 default). 0.999 = near-uniform.
    """

    def __init__(
        self,
        labels,
        beta: float = 0.999,
    ):
        super().__init__()
        labels = np.asarray(labels).astype(np.int64)
        n_classes = int(labels.max()) + 1
        if n_classes < 2:
            raise ValueError("Need at least 2 classes")
        counts = np.bincount(labels, minlength=n_classes).astype(np.float64)
        # Effective number per Cui 2019:
        #   E_n = (1 - beta^n) / (1 - beta)
        # beta=0 → E_n = n (raw inverse frequency)
        # beta=1 → E_n = uniform (no rebalancing; limit case)
        if beta == 0:
            weights = 1.0 / np.maximum(counts, 1)
        else:
            e_n = (1.0 - np.power(beta, counts)) / (1.0 - beta)
            weights = 1.0 / np.maximum(e_n, 1e-12)
        # Normalize so the mean weight is 1 (preserves loss magnitude).
        weights = weights / weights.mean()
        self.register_buffer(
            "weight", torch.tensor(weights, dtype=torch.float32)
        )
        self._n_classes = n_classes

    def forward(
        self,
        logits: "torch.Tensor",
        targets: "torch.Tensor",
    ) -> "torch.Tensor":
        return F.cross_entropy(logits, targets, weight=self.weight)


class CalibrationLoss(nn.Module):
    """Expected calibration error (ECE) surrogate for training.

    The closed-form ECE is non-differentiable (bin assignment uses
    ``argmax``). This surrogate replaces the bin assignment with a
    soft assignment via softmax temperature so gradients can flow.

    Reference: Mukhoti et al. (2020), "Calibrating Deep Neural
    Networks using Focal Loss." NeurIPS.

    Parameters
    ----------
    n_bins : int
        Number of probability bins for ECE.
    temperature : float
        Soft-assignment temperature. Higher = harder bins.
    """

    def __init__(self, n_bins: int = 15, temperature: float = 10.0):
        super().__init__()
        self.n_bins = n_bins
        self.temperature = temperature

    def forward(
        self,
        logits: "torch.Tensor",
        targets: "torch.Tensor",
    ) -> "torch.Tensor":
        # Convert to 2-class softmax so we can read off class-1 probs.
        if logits.dim() == 1:
            logits = torch.stack([-logits, logits], dim=-1)
        probs = F.softmax(logits, dim=-1)[:, 1]
        targets = targets.float()
        bin_edges = torch.linspace(0, 1, self.n_bins + 1, device=probs.device)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Soft assignment: weight each (sample, bin) pair by similarity
        # of the sample's predicted prob to the bin center.
        # diff: (N, n_bins)
        diff = -self.temperature * (probs.unsqueeze(-1) - bin_centers).abs()
        soft_assignment = F.softmax(diff, dim=-1)

        # Per-bin accuracy and confidence
        bin_acc = (soft_assignment * targets.unsqueeze(-1)).sum(0)
        bin_conf = (soft_assignment * probs.unsqueeze(-1)).sum(0)
        bin_count = soft_assignment.sum(0).clamp_min(1e-8)
        ece = (bin_count / bin_count.sum() * (bin_acc - bin_conf).abs()).sum()
        return ece


__all__ = [
    "SensAtSpecLoss",
    "BalancedCrossEntropy",
    "CalibrationLoss",
    "focal_binary_cross_entropy",
]
