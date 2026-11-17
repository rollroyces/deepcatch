#!/usr/bin/env python3
"""Label-permutation null controls for the paired-cohort benchmarks.

WHY THIS MODULE EXISTS
----------------------
The "shuffled-label control" was implemented separately in
``scripts/foundation_real_smoke.py``,
``scripts/fusion_regularized_smoke.py`` and
``scripts/cross_study_finallydb.py``. Two of those copies had drifted
apart, and one of them was a NO-OP:

``scripts/foundation_real_smoke.py`` built its "pair-broken shuffle" as::

    pos_shuf = pos_idx.copy(); rng.shuffle(pos_shuf)
    y_shuf[pos_shuf] = 1; y_shuf[neg_shuf] = 0

A permutation of a set has the same set as its image, so that writes 1s
back to exactly the positions that already held 1s: ``y_shuf == y_true``
for every seed. The "shuffled" AUCs the script reported were therefore
the HONEST AUCs, and the control proved nothing. The fingerprint is in
the committed artifact: ``shuffled_lr_baseline_auc_mean`` was
byte-identical to ``lr_baseline_auc_mean``, and
``shuffled_naive_avg_auc_mean`` was byte-identical to
``naive_avg_auc_mean``.

A control implemented twice will drift; a control implemented once
cannot. This module is the single implementation, and every script
imports it.

THE CONSTRUCTION
----------------
``foldwise_permuted_labels`` permutes the labels WITHIN each fold's
training block and, independently, within each fold's test block. The
model is then fitted on the permuted training labels and scored against
the permuted test labels, so the whole CV pipeline sees a relabelled
cohort.

``global_permuted_labels`` permutes the label vector as a whole. It is
the right control for an unpaired cohort (the cross-study FinaleDB
benchmark), where there is no fold-local pairing to break.

WHY THE NULL IS NOT 0.5 ON A PAIRED COHORT
------------------------------------------
Every patient contributes one positive and one negative arm whose panel
scores differ by ~+127 (the positive arm is higher for 100% of
patients), so a model fitted to random labels still recovers that
within-pair axis and lands BELOW 0.5. The measured null on this cohort
is 0.42 +/- 0.12 over 50 permutations per seed -- see
``results/fusion_regularized_smoke.json`` ->
``per_seed/*/diagnostics/empirical_null``. A control below 0.5 is
therefore NOT on its own evidence of leakage here; it has to be read
against that measured null.

The measured null is a property of the COHORT, not of the estimator
under test, so it is the same for every estimator measured on it.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

# Offset applied to the caller's seed before the permutation RNG is
# created. It is part of the construction's identity: changing it
# changes every published control number, so it is a named constant
# rather than a literal scattered through three scripts.
SEED_OFFSET = 7777


def foldwise_permuted_labels(
    y: np.ndarray,
    folds: Sequence[Tuple[np.ndarray, np.ndarray]],
    seed: int,
) -> np.ndarray:
    """Permute labels independently inside each fold's train/test blocks.

    Parameters
    ----------
    y : ndarray
        The label vector.
    folds : sequence of (train_idx, test_idx) pairs
        The CV splits the control must be aligned with. The returned
        vector has the same length as ``y``; row ``i`` holds the
        permuted label of row ``i``.
    seed : int
        Deterministic seed. The permutation RNG is
        ``np.random.default_rng(seed + SEED_OFFSET)``.

    Returns
    -------
    ndarray
        A relabelled copy of ``y``. The class balance of each block is
        preserved; the association between the scores and the labels is
        destroyed.
    """
    y = np.asarray(y)
    rng = np.random.default_rng(seed + SEED_OFFSET)
    out = np.zeros_like(y)
    for tr, te in folds:
        out[tr] = rng.permutation(y[tr])
        out[te] = rng.permutation(y[te])
    return out


def global_permuted_labels(y: np.ndarray, seed: int) -> np.ndarray:
    """Permute the whole label vector.

    Used for an unpaired cohort, where there is no fold-local pairing
    to break and the only thing a null control has to destroy is the
    association between the scores and the labels. The class balance of
    the vector as a whole is preserved.

    Parameters
    ----------
    y : ndarray
        The label vector.
    seed : int
        Deterministic seed; the permutation RNG is
        ``np.random.default_rng(seed)`` -- no offset, because this
        construction predates ``SEED_OFFSET`` and its published control
        numbers depend on the exact stream.

    Returns
    -------
    ndarray
        A relabelled copy of ``y``.
    """
    return np.random.default_rng(seed).permutation(np.asarray(y))


def is_degenerate(y: np.ndarray, y_perm: np.ndarray) -> bool:
    """True when the permutation left the labels untouched.

    A control that reproduces the real labels is not a control: it
    reports the honest AUC back as if it were a null. This is the
    check that would have caught the no-op shuffle in
    ``scripts/foundation_real_smoke.py``.
    """
    return bool(np.array_equal(np.asarray(y_perm), np.asarray(y)))


def assert_not_degenerate(y: np.ndarray, y_perm: np.ndarray) -> None:
    """Raise ``RuntimeError`` when the permutation is the identity.

    Used as a guard inside the scripts so a silent regression to a
    no-op control fails loudly at run time instead of publishing a
    number that means nothing.
    """
    if is_degenerate(y, y_perm):
        raise RuntimeError(
            "label-permutation control degenerated: y_perm == y. A control "
            "that reproduces the real labels is not a control."
        )
