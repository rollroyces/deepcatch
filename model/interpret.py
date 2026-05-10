#!/usr/bin/env python3
"""
Integrated Gradients for Motif Importance Explanation
======================================================

Provides feature attribution for the MotifDiversityModel, identifying
which BPE tokens (subword motifs) most influenced a cancer prediction.

.. rubric:: Why Integrated Gradients?

Integrated Gradients [1]_ is a model-agnostic attribution method that
satisfies the **completeness** axiom: attributions sum to the difference
between the prediction at the input and at a baseline.  For cfDNA
fragmentomics, this is crucial — we need to know not just *that* a
prediction was made, but *which motifs drove it*.

.. rubric:: How It Works (DNA Context)

1. **Baseline**: A uniform frequency distribution (all motifs equally
   likely) represents the "null hypothesis" of no fragmentation bias.
2. **Path integration**: Linearly interpolate between the baseline and
   the observed frequency histogram in ``n_steps`` increments.
3. **Gradient accumulation**: At each step, compute the model gradient
   w.r.t. the input.  Sum and multiply by (input − baseline).
4. **Attribution vector**: Shape ``(vocab_size,)`` — positive values
   mean the motif *increases* cancer probability, negative values mean
   it *decreases* it.

.. rubric:: Output Format

The :func:`explain_prediction` function returns:

- **attributions**: Raw IG values per motif.
- **top_motifs**: The 10 motifs with the strongest positive influence,
  each with human-readable token string and attribution score.
- **bottom_motifs**: The 10 motifs with the strongest negative influence
  (protective / healthy signal).

.. rubric:: Performance on Apple Silicon

Each IG explanation requires ``n_steps`` forward+backward passes
(default 50).  With the default 86K-parameter model, this completes
in <100 ms on M1/M2.

.. rubric:: References

.. [1] Sundararajan, M., Taly, A., & Yan, Q. (2017).
   "Axiomatic Attribution for Deep Neural Networks." ICML.
   arXiv:1703.01365
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import torch

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None  # type: ignore[assignment]


class IntegratedGradientsExplainer:
    """
    Integrated Gradients explainer for motif frequency inputs.

    Parameters
    ----------
    model : torch.nn.Module
        A trained ``MotifDiversityModel`` (or any model with a
        ``forward(frequencies) → probabilities`` signature).
    n_steps : int
        Number of interpolation steps (default 50).  Higher = more
        accurate gradients, lower = faster.
    """

    def __init__(
        self,
        model,
        n_steps: int = 50,
    ):
        if not _HAS_TORCH:
            raise ImportError(
                "PyTorch is required for Integrated Gradients. "
                "Install: pip install torch"
            )
        self.model = model
        self.n_steps = n_steps

    def attribute(
        self,
        inputs,
        baseline=None,
        return_convergence_delta: bool = False,
    ):
        """
        Compute Integrated Gradients attributions.

        Parameters
        ----------
        inputs : Tensor of shape ``(batch, vocab_size)``
            Input frequency histogram(s).
        baseline : Tensor of shape ``(1, vocab_size)`` or ``(batch, vocab_size)``, optional
            Baseline input.  Default: uniform distribution.
        return_convergence_delta : bool
            If True, also return an approximate convergence error.

        Returns
        -------
        attributions : Tensor of shape ``(batch, vocab_size)``
        delta : Tensor of shape ``(batch,)`` (only if
            ``return_convergence_delta=True``)
        """
        if baseline is None:
            vocab_size = inputs.shape[-1]
            baseline = torch.ones(1, vocab_size) / vocab_size
            baseline = baseline.to(inputs.device).expand_as(inputs)

        if baseline.shape[0] == 1 and inputs.shape[0] > 1:
            baseline = baseline.expand_as(inputs)

        attributions = torch.zeros_like(inputs)

        for step in range(self.n_steps):
            alpha = (step + 1) / self.n_steps
            scaled = baseline + alpha * (inputs - baseline)
            scaled.requires_grad_(True)

            outputs = self.model(scaled)
            if outputs.dim() > 1:
                outputs = outputs.squeeze(-1)

            grad_outputs = torch.ones_like(outputs)
            grads = torch.autograd.grad(
                outputs=outputs,
                inputs=scaled,
                grad_outputs=grad_outputs,
                create_graph=False,
                retain_graph=False,
            )[0]

            attributions += grads.detach()

        attributions = attributions * (inputs - baseline) / self.n_steps

        if return_convergence_delta:
            with torch.no_grad():
                f_x = self.model(inputs).squeeze(-1)
                f_b = self.model(baseline).squeeze(-1)
                delta = torch.abs(
                    attributions.sum(dim=-1) - (f_x - f_b)
                )
            return attributions, delta

        return attributions

    def explain(
        self, frequencies: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        High-level explanation interface.

        Parameters
        ----------
        frequencies : np.ndarray of shape ``(vocab_size,)``
            Input frequency histogram.

        Returns
        -------
        dict with:
            - ``attributions``: raw IG values per motif.
        """
        device = next(self.model.parameters()).device
        freq_tensor = (
            torch.from_numpy(frequencies).float().unsqueeze(0).to(device)
        )
        attr = self.attribute(freq_tensor)
        return {"attributions": attr.cpu().numpy().squeeze(0)}


def explain_prediction(
    frequencies: np.ndarray,
    model,
    tokenizer=None,
    top_k: int = 10,
    n_steps: int = 50,
) -> Dict:
    """
    Explain a model prediction with Integrated Gradients.

    Parameters
    ----------
    frequencies : np.ndarray of shape ``(vocab_size,)``
        Normalized motif frequency histogram.
    model : torch.nn.Module
        Trained MotifDiversityModel.
    tokenizer : DNATokenizer, optional
        Tokenizer for decoding token IDs to human-readable strings.
        If None, tokens are reported by ID.
    top_k : int
        Number of top positive and negative motifs to return.
    n_steps : int
        IG interpolation steps.

    Returns
    -------
    dict
        {
            "attributions": [...],
            "top_motifs": [{"token": "CG", "id": 12, "attribution": 0.0342}, ...],
            "bottom_motifs": [{"token": "AT", "id": 3, "attribution": -0.0198}, ...],
            "prediction": float,
            "convergence_delta": float or None,
        }
    """
    if not _HAS_TORCH:
        warnings.warn(
            "PyTorch not installed — cannot compute Integrated Gradients. "
            "Returning empty explanation."
        )
        return {
            "attributions": [],
            "top_motifs": [],
            "bottom_motifs": [],
            "prediction": 0.0,
            "convergence_delta": None,
            "error": "PyTorch not available",
        }

    explainer = IntegratedGradientsExplainer(model, n_steps=n_steps)
    device = next(model.parameters()).device

    freq_tensor = (
        torch.from_numpy(frequencies).float().unsqueeze(0).to(device)
    )

    attributions, delta = explainer.attribute(
        freq_tensor, return_convergence_delta=True
    )
    attr_np = attributions.cpu().numpy().squeeze(0)

    # Top and bottom motifs
    sorted_idx = np.argsort(attr_np)[::-1]
    top_idx = sorted_idx[:top_k]
    bottom_idx = sorted_idx[-top_k:][::-1]

    # Decode token IDs
    id_to_token: Dict[int, str] = {}
    if tokenizer is not None:
        id_to_token = {v: k for k, v in tokenizer.vocab.items()}

    top_motifs = []
    for idx in top_idx:
        top_motifs.append(
            {
                "token": id_to_token.get(int(idx), f"<{idx}>"),
                "id": int(idx),
                "attribution": float(attr_np[idx]),
            }
        )

    bottom_motifs = []
    for idx in bottom_idx:
        bottom_motifs.append(
            {
                "token": id_to_token.get(int(idx), f"<{idx}>"),
                "id": int(idx),
                "attribution": float(attr_np[idx]),
            }
        )

    # Get prediction
    with torch.no_grad():
        pred = float(model(freq_tensor).cpu().item())

    return {
        "attributions": attr_np.tolist(),
        "top_motifs": top_motifs,
        "bottom_motifs": bottom_motifs,
        "prediction": pred,
        "convergence_delta": (
            float(delta.cpu().item()) if delta is not None else None
        ),
    }


__all__ = [
    "IntegratedGradientsExplainer",
    "explain_prediction",
]
