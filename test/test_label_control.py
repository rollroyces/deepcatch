"""Regression tests for the label-permutation null control.

THE DEFECT THIS FILE GUARDS AGAINST
-----------------------------------
``scripts/foundation_real_smoke.py`` used to build its "pair-broken
shuffle" as::

    pos_shuf = pos_idx.copy(); rng.shuffle(pos_shuf)
    y_shuf[pos_shuf] = 1; y_shuf[neg_shuf] = 0

A permutation of a set has the same set as its image, so that writes 1s
back to exactly the positions that already held 1s: ``y_shuf == y_true``
for every seed. The "shuffled" AUCs the script reported were therefore
the HONEST AUCs, not a null.

The fingerprint was visible in the committed artifact: in
``results/sens_at_spec_ce_n20.json`` the control
``shuffled_lr_baseline_auc_mean`` was byte-identical to the honest
``lr_baseline_auc_mean``, and ``shuffled_naive_avg_auc_mean`` was
byte-identical to ``naive_avg_auc_mean``. A control that reproduces the
real number exactly is not a control.

These tests are the structural guard: they check the property that was
violated, not the specific numbers, so they keep working when the
cohort or the estimator changes.

Run: python -m pytest test/test_label_control.py -v
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "scripts"))

from src.foundation.label_control import (  # noqa: E402
    assert_not_degenerate,
    foldwise_permuted_labels,
    global_permuted_labels,
    is_degenerate,
)


# --------------------------------------------------------------------------
# The permutation itself
# --------------------------------------------------------------------------
def test_foldwise_permutation_changes_the_labels():
    """The control must actually relabel the cohort.

    This is the property the no-op shuffle violated: a permutation of
    the positive index set among itself, written back to the same
    positions, leaves the label vector unchanged.
    """
    y, folds = _cv_folds()
    permuted = foldwise_permuted_labels(y, folds, seed=0)
    assert not is_degenerate(y, permuted)
    assert not np.array_equal(permuted, y)


def test_global_permutation_changes_the_labels():
    y, _ = _cv_folds()
    permuted = global_permuted_labels(y, seed=0)
    assert not is_degenerate(y, permuted)
    assert not np.array_equal(permuted, y)


def _cv_folds(n=40, n_splits=5):
    """Realistic CV folds: every fold's train and test blocks together
    cover the whole cohort, exactly as ``GroupKFold.split`` produces.

    The permutation relies on this: each fold writes a permuted copy of
    its own train block and its own test block, so with a well-formed
    split the class balance of the cohort is preserved.
    """
    y = np.array([1] * (n // 2) + [0] * (n // 2), dtype=np.int64)
    groups = np.concatenate([np.arange(n // 2), np.arange(n // 2)])
    from sklearn.model_selection import GroupKFold
    folds = list(GroupKFold(n_splits=n_splits).split(
        np.zeros(n), y, groups=groups))
    return y, folds


def test_the_permutations_preserve_the_class_balance():
    """A null control must not change how many positives there are.

    The permuted cohort has to have the same class balance as the real
    one, otherwise the control is also measuring a prevalence shift and
    its AUC is not comparable to the honest AUC.
    """
    y, folds = _cv_folds()
    for permuted in (foldwise_permuted_labels(y, folds, seed=0),
                     global_permuted_labels(y, seed=0)):
        assert int((permuted == 1).sum()) == int((y == 1).sum())
        assert int((permuted == 0).sum()) == int((y == 0).sum())


def test_the_permutation_is_deterministic():
    """Same seed, same permutation -- the artifact must be reproducible."""
    # The labels are INTERLEAVED, not block-constant: permuting a block
    # of identical labels cannot detect a change of seed, because every
    # permutation of it looks the same.
    y = np.array([1, 0] * 10, dtype=np.int64)
    folds = [(np.arange(0, 10), np.arange(10, 20))]
    a = foldwise_permuted_labels(y, folds, seed=7)
    b = foldwise_permuted_labels(y, folds, seed=7)
    assert np.array_equal(a, b)
    c = foldwise_permuted_labels(y, folds, seed=8)
    assert not np.array_equal(a, c)


def test_the_degeneracy_guard_raises_on_the_identity():
    """The guard must refuse to accept the no-op as a control."""
    y = np.array([1, 0, 1, 0], dtype=np.int64)
    with pytest.raises(RuntimeError, match="degenerated"):
        assert_not_degenerate(y, y.copy())


def test_the_degeneracy_guard_accepts_a_real_permutation():
    y = np.array([1, 0, 1, 0], dtype=np.int64)
    assert_not_degenerate(y, global_permuted_labels(y, seed=3))


# --------------------------------------------------------------------------
# The committed artifacts
# --------------------------------------------------------------------------
def _paired_cohort_artifacts():
    """The committed artifacts produced by the CURRENT smoke script.

    The n=5 pair (``sens_at_spec_ce.json`` / ``sens_at_spec_sens.json``)
    is deliberately NOT included: those were generated by the audit-1
    script under the leaky ``StratifiedKFold`` era and are kept as a
    historical record only (see ``docs/SENS_AT_SPEC_ABLATION.md``,
    "Original n=5 ablation ... THIS SECTION IS HISTORICAL AND
    SUPERSEDED"). Their control numbers belong to a superseded
    protocol and are not evidence about the current one, so asserting
    anything about them would test the wrong thing -- and regenerating
    them would rewrite a historical record.
    """
    out = []
    for name in ("sens_at_spec_ce_n20.json",
                 "sens_at_spec_sens_n20.json",
                 "sparse_aware_linear.json",
                 "sparse_aware_panel.json"):
        path = _ROOT / "results" / name
        if path.exists():
            out.append((name, json.loads(path.read_text())))
    return out


@pytest.mark.parametrize("name,artifact", _paired_cohort_artifacts())
def test_the_committed_control_is_not_the_honest_number(name, artifact):
    """The committed control must NOT equal the honest number.

    This is the assertion that would have caught the defect: with the
    no-op shuffle, ``shuffled_lr_baseline_auc_mean`` was byte-identical
    to ``lr_baseline_auc_mean`` and the control reported the honest AUC
    back as if it were a null.
    """
    honest = artifact["lr_baseline_auc_mean"]
    control = artifact["shuffled_lr_baseline_auc_mean"]
    assert control != honest, (
        f"{name}: the shuffled-label control equals the honest LR AUC "
        f"({honest}). The permutation is a no-op, so the control is not "
        f"a control."
    )

    honest_naive = artifact["naive_avg_auc_mean"]
    control_naive = artifact["shuffled_naive_avg_auc_mean"]
    assert control_naive != honest_naive, (
        f"{name}: the shuffled-label control equals the honest naive-average "
        f"AUC ({honest_naive}). The permutation is a no-op, so the control "
        f"is not a control."
    )


@pytest.mark.parametrize("name,artifact", _paired_cohort_artifacts())
def test_the_committed_control_is_in_the_plausible_null_range(name, artifact):
    """The control must land where a null lands, not where the signal is.

    On this paired design the measured null is around 0.42 (it sits
    BELOW 0.5 because a model fitted to random labels can still recover
    the within-pair axis). The honest LR AUC on this cohort is ~0.90.
    A control that lands near the honest number is not a null; a
    control that lands far outside [0, 1] is a broken artifact.
    """
    control = artifact["shuffled_lr_baseline_auc_mean"]
    assert 0.0 <= control <= 1.0, (
        f"{name}: the shuffled-label control AUC is out of range: {control}"
    )
    honest = artifact["lr_baseline_auc_mean"]
    assert abs(control - honest) > 0.05, (
        f"{name}: the shuffled-label control ({control}) is within 0.05 of "
        f"the honest AUC ({honest}); it is not controlling for anything."
    )


def test_the_committed_foundation_control_is_not_the_honest_foundation():
    """The foundation control must not reproduce the honest foundation AUC.

    The control and the honest number come from the same estimator on
    the same folds, so the only thing that can separate them is the
    label assignment. If they coincide, the permutation did nothing.
    """
    found = False
    for name, artifact in _paired_cohort_artifacts():
        honest = artifact["foundation_auc_mean"]
        control = artifact["shuffled_foundation_auc_mean"]
        assert control != honest, (
            f"{name}: the shuffled-label foundation control equals the "
            f"honest foundation AUC ({honest})."
        )
        found = True
    if not found:
        pytest.skip("no paired-cohort artifacts present")


# --------------------------------------------------------------------------
# The scripts
# --------------------------------------------------------------------------
def test_the_scripts_share_one_control_implementation():
    """The control must not be re-implemented per script.

    The defect survived because the permutation was written out
    separately in three scripts, and the copies drifted apart. This
    test fails if a script grows its own permutation helper again.
    """
    for script in ("foundation_real_smoke.py", "fusion_regularized_smoke.py",
                   "cross_study_finallydb.py"):
        path = _ROOT / "scripts" / script
        if not path.exists():
            continue
        text = path.read_text()
        assert "def _foldwise_permuted_labels" not in text, (
            f"scripts/{script} re-implements the foldwise permutation; "
            f"use src/foundation/label_control.py instead."
        )
        assert "def _global_permuted_labels" not in text, (
            f"scripts/{script} re-implements the global permutation; "
            f"use src/foundation/label_control.py instead."
        )


def test_the_no_op_shuffle_construction_is_gone():
    """The exact no-op construction must not reappear.

    ``y_shuf[pos_shuf] = 1`` where ``pos_shuf`` is a permutation of the
    positive indices writes 1s back to the positions that already held
    them. The construction is easy to reintroduce by accident, so the
    source is checked for it directly.
    """
    path = _ROOT / "scripts" / "foundation_real_smoke.py"
    text = path.read_text()
    assert "y_shuf[pos_shuf] = 1" not in text, (
        "the no-op pair-broken shuffle is back in "
        "scripts/foundation_real_smoke.py"
    )
    assert "y_shuf[neg_shuf] = 0" not in text, (
        "the no-op pair-broken shuffle is back in "
        "scripts/foundation_real_smoke.py"
    )
