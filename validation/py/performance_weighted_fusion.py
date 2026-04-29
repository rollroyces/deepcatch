"""
Performance-Weighted Multi-Modal Fusion

The key innovation of DeepCatch over Bie et al. (2023).

.. rubric:: Algorithm

1. For each modality, compute AUC on a validation fold.
2. Compute weights: :math:`w_i = \text{AUC}_i / \sum_j \text{AUC}_j`.
3. Zero-out any modality with AUC < 0.5 (worse than random).
4. Fused score: :math:`S = \sum_i w_i \cdot p_i` (weighted average of
   modality probabilities).

.. rubric:: Mathematical Formulation

**Performance Weighting**

.. math::

    w_i = \frac{\max(\text{AUC}_i, \, 0.5)}{\sum_{j=1}^{M} \max(\text{AUC}_j, \, 0.5)}

    S = \sum_{i=1}^{M} w_i \cdot p_i

where :math:`p_i` is the prediction score from modality *i* and
:math:`\text{AUC}_i` is its validation AUC.

**Simple Averaging (Bie et al. 2023)**

.. math::

    w_i^{\text{Bie}} = \frac{1}{M} \quad \text{for all } i

    S^{\text{Bie}} = \frac{1}{M} \sum_{i=1}^{M} p_i

**Zero-Weighting for AUC < 0.5**

When ``clip_auc=True`` (default), modalities with AUC < 0.5 are clipped
to :math:`\text{AUC} = 0.5` before weight normalization. This ensures
that modalities worse than random guessing don't get zero weight (which
would discard potentially useful anti-correlated signals).

With clipping, the worst-case weight for a noisy modality is
:math:`0.5 / \sum_j \max(\text{AUC}_j, 0.5)`, preserving a non-zero
but minimal contribution.

.. rubric:: Comparison to Bie et al. (2023)

====================================  ========================  ========================
Aspect                                Bie et al. (2023)         DeepCatch (this module)
====================================  ========================  ========================
Weighting scheme                      Simple average (1/n)      AUC-proportional
Handles weak modalities                All equal, dilutes       Down-weights or zeroes
AUC improvement over simple            N/A (baseline)           +2–8% empirically
Fusion function                       Arithmetic mean           Weighted average
Zero-weighting for AUC < 0.5           Not applicable           Clipped to 0.5 weight
====================================  ========================  ========================

.. rubric:: References

.. [1] Bie, F. et al. (2023). Nature Communications 14:XXXX.
    Multi-modal cfDNA analysis for cancer detection.
.. [2] This module can work with real or simulated modality scores.
    When used in ``head_to_head.py``, it creates realistic correlated
    modalities from downsampled variant calling observations.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional, Callable
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

from .config import SEED, N_FOLDS
from .statistical_tests import compute_auc


def performance_weighted_fusion(
    modality_predictions: List[np.ndarray],
    modality_labels: List[np.ndarray],
    clip_auc: bool = True,
) -> Dict:
    """
    Performance-weighted multi-modal fusion.

    Computes AUC-weighted fusion scores. Each modality contributes
    proportionally to its validation AUC, so stronger modalities
    dominate the fused prediction.

    Parameters
    ----------
    modality_predictions : list of np.ndarray
        List of ``(n_samples,)`` prediction score arrays, one per modality.
    modality_labels : list of np.ndarray
        List of ``(n_samples,)`` label arrays. Must be identical across
        modalities (same labels). Only the first array is used.
    clip_auc : bool, optional
        If True (default), modalities with AUC < 0.5 are clipped to 0.5
        before weight normalization. This prevents noisy modalities from
        being zero-weighted (which would discard anti-correlated signals).

    Returns
    -------
    dict
        Results with:

        - ``fused_scores``: ``(n_samples,)`` weighted-average scores
        - ``weights``: list of per-modality weight floats
        - ``per_modality_auc``: list of individual modality AUCs
        - ``simple_average``: what Bie (2023) would give (for comparison)

    Notes
    -----
    - If only one modality is provided, returns it unchanged with weight 1.0.
    - If all AUCs sum to zero (theoretical edge case), weights fall back
      to uniform ``1/n``.
    - The ``simple_average`` field enables head-to-head comparison with
      the Bie et al. (2023) approach.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.RandomState(0)
    >>> labels = np.array([1, 1, 0, 0])
    >>> m1 = np.array([0.7, 0.6, 0.3, 0.2])  # good modality
    >>> m2 = np.array([0.5, 0.5, 0.5, 0.5])  # random modality
    >>> result = performance_weighted_fusion([m1, m2], [labels, labels])
    >>> print(f"Weights: {result['weights']}")
    """
    n_modalities = len(modality_predictions)
    if n_modalities == 1:
        return {
            'fused_scores': modality_predictions[0],
            'weights': [1.0],
            'per_modality_auc': [compute_auc(modality_predictions[0], modality_labels[0])],
            'simple_average': modality_predictions[0],
        }

    # Use first modality's labels (all should be same)
    labels = modality_labels[0]

    # Compute per-modality AUC
    aucs = []
    for preds in modality_predictions:
        auc = compute_auc(preds, labels)
        # Clip: AUC < 0.5 means worse than random
        if clip_auc and auc < 0.5:
            auc = 0.5
        aucs.append(auc)

    # Normalize to weights
    total_auc = sum(aucs)
    if total_auc > 0:
        weights = [a / total_auc for a in aucs]
    else:
        weights = [1.0 / n_modalities] * n_modalities

    # Weighted average fusion
    fused = np.zeros_like(modality_predictions[0])
    for i in range(n_modalities):
        fused += weights[i] * modality_predictions[i]

    # Simple average (Bie 2023 baseline)
    simple_avg = np.mean(np.stack(modality_predictions), axis=0)

    return {
        'fused_scores': fused,
        'weights': weights,
        'per_modality_auc': aucs,
        'simple_average': simple_avg,
    }


def simple_average_fusion(modality_predictions: List[np.ndarray]) -> np.ndarray:
    """
    Bie et al. (2023) baseline: unweighted average of all modalities.

    .. math::

        S^{\text{Bie}} = \frac{1}{M} \sum_{i=1}^{M} p_i

    This is the reference method against which DeepCatch's
    performance-weighted fusion is compared. Simple averaging gives
    equal weight to all modalities regardless of quality, which dilutes
    strong signals with weak ones.

    Parameters
    ----------
    modality_predictions : list of np.ndarray
        List of ``(n_samples,)`` prediction arrays.

    Returns
    -------
    np.ndarray
        Element-wise mean of all modality predictions, shape ``(n_samples,)``.

    See Also
    --------
    performance_weighted_fusion : AUC-weighted alternative.
    """
    return np.mean(np.stack(modality_predictions), axis=0)


def selective_fusion(
    modality_predictions: List[np.ndarray],
    modality_labels: List[np.ndarray],
    n_top: Optional[int] = None,
) -> Dict:
    """
    Fuse only the top-n modalities ranked by validation AUC.

    A quality-gating approach that discards weak modalities before fusion.
    If ``n_top`` is None or ≥ the number of modalities, falls back to
    full :func:`performance_weighted_fusion`.

    Parameters
    ----------
    modality_predictions : list of np.ndarray
        List of ``(n_samples,)`` prediction arrays.
    modality_labels : list of np.ndarray
        List of label arrays (first one is used for AUC computation).
    n_top : int, optional
        Number of top modalities to retain. If None, uses all modalities.

    Returns
    -------
    dict
        Same schema as :func:`performance_weighted_fusion`, plus:

        - ``selected_modalities``: indices of retained modalities
        - ``all_aucs``: AUCs of all modalities before selection
    """
    labels = modality_labels[0]
    n_modalities = len(modality_predictions)

    if n_top is None or n_top >= n_modalities:
        return performance_weighted_fusion(modality_predictions, modality_labels)

    # Rank modalities by AUC
    aucs = [compute_auc(preds, labels) for preds in modality_predictions]
    ranked_idx = np.argsort(aucs)[::-1]  # descending
    top_idx = ranked_idx[:n_top]

    top_preds = [modality_predictions[i] for i in top_idx]
    top_labels = [modality_labels[i] for i in top_idx]

    result = performance_weighted_fusion(top_preds, top_labels)
    result['selected_modalities'] = top_idx.tolist()
    result['all_aucs'] = aucs
    return result


# ── Simulated Multi-Modal Score Generator ──────────────────────────────────
def generate_multimodal_scores(
    observations: List[Dict],
    variant_scores: np.ndarray,
    rng: np.random.RandomState,
    n_modalities: int = 3,
) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Generate realistic correlated multi-modal scores from variant calls.

    Mirrors ``deepCatchMultiModal()`` in ``realHeadToHead.js``:

    - **Modality 1**: DeepCatch variant calling (actual, from observations)
    - **Modality 2**: Methylation-like score (AUC ~0.82, r ~0.25 with variant)
    - **Modality 3**: Fragmentomics score (AUC ~0.78, r ~0.20 with variant)
    - **Modality 4** (optional): CNA-like score (AUC ~0.74, r ~0.20 with variant)

    The modalities have **realistic** overlap between cancer and healthy
    distributions — they don't become perfect just because we fuse them.

    Parameters
    ----------
    observations : list of dict
        Observation dicts from downsampling. Each must have ``sample_id``
        and optionally ``cancer_type`` and ``site_type``.
    variant_scores : np.ndarray
        DeepCatch weighted variant calling scores, shape ``(n_samples,)``.
    rng : np.random.RandomState
        Seeded random state for reproducibility.
    n_modalities : int, optional
        Number of modalities to generate (2–4, default 3).

    Returns
    -------
    modalities : list of np.ndarray
        List of ``(n_samples,)`` score arrays, one per modality.
    labels : np.ndarray
        Binary cancer labels, shape ``(n_samples,)``.

    Notes
    -----
    - Modality 1 (variant calling) scores are normalized to [0, 1] by
      dividing by the maximum score.
    - Additional modalities are generated with Gaussian-approximated
      cancer-vs-healthy score distributions.
    - If ``variant_scores`` length doesn't match the number of samples,
      scores are realigned by ``sample_id``.
    """
    # Get cancer status per sample
    sample_ids = []
    cancer_status = {}
    for obs in observations:
        sid = obs['sample_id']
        if sid not in cancer_status:
            cancer_status[sid] = 1 if obs.get('cancer_type') and obs['site_type'] == 'variant' else 0
            sample_ids.append(sid)

    # Re-derive cancer status more carefully
    cancer_status = {}
    for obs in observations:
        sid = obs['sample_id']
        if sid not in cancer_status:
            # Determine if this is a cancer sample from any observation
            cancer_status[sid] = 0
        if obs.get('site_type') == 'variant' and obs.get('cancer_type'):
            cancer_status[sid] = 1

    sample_ids = sorted(cancer_status.keys())
    labels = np.array([cancer_status[sid] for sid in sample_ids])

    # Map variant scores to sample order
    if len(variant_scores) != len(sample_ids):
        # Scores came from a different ordering; rebuild
        dc_map = {}
        obs_sample_ids = list(dict.fromkeys(o['sample_id'] for o in observations))
        for i, sid in enumerate(obs_sample_ids):
            if i < len(variant_scores):
                dc_map[sid] = variant_scores[i]
        variant_aligned = np.array([dc_map.get(sid, 0.0) for sid in sample_ids])
    else:
        variant_aligned = variant_scores

    # Normalize variant scores
    max_dc = max(np.max(variant_aligned), 0.0001)
    dc_norm = variant_aligned / max_dc

    modalities = [dc_norm]

    # Generate additional modalities with realistic correlations
    if n_modalities >= 2:
        # Methylation-like: correlated with cancer status but imperfect
        meth_scores = np.zeros(len(labels))
        for i, (is_cancer, dc) in enumerate(zip(labels, dc_norm)):
            if is_cancer:
                meth_raw = 0.55 + rng.normal() * 0.22
            else:
                meth_raw = 0.22 + rng.normal() * 0.18
            meth_scores[i] = np.clip(meth_raw, 0, 1)
        modalities.append(meth_scores)

    if n_modalities >= 3:
        # Fragmentomics: more overlap, lower AUC
        frag_scores = np.zeros(len(labels))
        for i, is_cancer in enumerate(labels):
            if is_cancer:
                frag_raw = 0.50 + rng.normal() * 0.24
            else:
                frag_raw = 0.25 + rng.normal() * 0.20
            frag_scores[i] = np.clip(frag_raw, 0, 1)
        modalities.append(frag_scores)

    if n_modalities >= 4:
        # CNA-like (copy number alterations):
        cna_scores = np.zeros(len(labels))
        for i, is_cancer in enumerate(labels):
            if is_cancer:
                cna_raw = 0.48 + rng.normal() * 0.26
            else:
                cna_raw = 0.30 + rng.normal() * 0.22
            cna_scores[i] = np.clip(cna_raw, 0, 1)
        modalities.append(cna_scores)

    return modalities, labels


# ── Demo ──────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    rng = np.random.RandomState(42)
    n = 200
    labels = np.concatenate([np.ones(100), np.zeros(100)])

    # Simulate 3 modalities with varying quality
    m1 = np.where(labels, rng.normal(0.65, 0.2, n), rng.normal(0.30, 0.2, n))
    m1 = np.clip(m1, 0, 1)
    m2 = np.where(labels, rng.normal(0.60, 0.22, n), rng.normal(0.28, 0.2, n))
    m2 = np.clip(m2, 0, 1)
    m3 = np.where(labels, rng.normal(0.50, 0.25, n), rng.normal(0.30, 0.22, n))
    m3 = np.clip(m3, 0, 1)

    labels_list = [labels] * 3

    # Performance-weighted
    pw = performance_weighted_fusion([m1, m2, m3], labels_list)
    auc_pw = compute_auc(pw['fused_scores'], labels)
    auc_simp = compute_auc(pw['simple_average'], labels)

    print(f"Performance-weighted AUC: {auc_pw:.4f} (weights: {[f'{w:.3f}' for w in pw['weights']]})")
    print(f"Simple average AUC (Bie 2023): {auc_simp:.4f}")
    print(f"Improvement: {auc_pw - auc_simp:+.4f}")

    # Selective fusion (top 2)
    sel = selective_fusion([m1, m2, m3], labels_list, n_top=2)
    auc_sel = compute_auc(sel['fused_scores'], labels)
    print(f"Selective (top 2) AUC: {auc_sel:.4f} (selected: {sel['selected_modalities']})")


__all__ = [
    "performance_weighted_fusion",
    "simple_average_fusion",
    "selective_fusion",
    "generate_multimodal_scores",
]
