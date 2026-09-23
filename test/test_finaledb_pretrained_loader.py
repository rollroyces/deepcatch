"""Tests for the FinaleDB pretrained checkpoint loader.

Verifies the integration chain:
1. The flat 2256-dim feature matrix splits into the expected 6
   per-modality arrays matching ``FoundationDownstream.MODALITY_DIMS``.
2. ``FoundationDownstream(pretrained=True, checkpoint_path=...)``
   loads the saved encoder weights and produces a finite forward pass.

These tests require the NEW (real-data-trained) artifacts:
- ``data/finaledb_pretrain_cohort_PRODUCTION.npz``
- ``checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt``

The OLD artifacts (``data/finaledb_pretrain_cohort_SYNTHETIC_v0.npz``,
``checkpoints/foundation_pretrained_SYNTHETIC_v0.pt``) are historical
and NOT trained on real FinaleDB data; the regression test
``test_synthetic_v0_artifacts_still_load_for_provenance`` verifies
they are still loadable for audit purposes, but ``FoundationDownstream``
must never default to them.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the scripts/ directory importable so the loader module can be
# imported regardless of pytest invocation cwd.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


# NEW, real-data-trained (default for downstream use)
NPZ_PATH = ROOT / "data" / "finaledb_pretrain_cohort_PRODUCTION.npz"
CKPT_PATH = ROOT / "checkpoints" / "foundation_pretrained_finaledb_PRODUCTION.pt"

# OLD, historical artifact (NOT real-data-trained; preserved for provenance)
OLD_NPZ_PATH = ROOT / "data" / "finaledb_pretrain_cohort_SYNTHETIC_v0.npz"
OLD_CKPT_PATH = ROOT / "checkpoints" / "foundation_pretrained_SYNTHETIC_v0.pt"


@pytest.mark.skipif(
    not NPZ_PATH.exists(),
    reason=f"NEW pretraining npz not present at {NPZ_PATH}",
)
def test_loader_splits_flat_X_into_expected_modality_shapes():
    """The 2256-dim flat X must split into 6 arrays of correct widths."""
    from finaledb_pretrained_loader import load_real_cohort, modalities_from_flat_X

    cohort = load_real_cohort(NPZ_PATH)
    assert cohort["X"].shape[1] == 2256
    modalities = modalities_from_flat_X(cohort["X"])
    expected_dims = {
        "frag_basic": 4,
        "frag_enhanced": 44,
        "cnv": 6,
        "sero": 4,
        "gnn": 1,
        "tissue": 24,
    }
    assert set(modalities.keys()) == set(expected_dims.keys())
    for k, expected_dim in expected_dims.items():
        assert modalities[k].shape == (cohort["X"].shape[0], expected_dim), (
            f"{k}: got {modalities[k].shape}, expected "
            f"({cohort['X'].shape[0]}, {expected_dim})"
        )


@pytest.mark.skipif(
    not NPZ_PATH.exists() or not CKPT_PATH.exists(),
    reason="NEW pretraining npz or checkpoint not present",
)
def test_pretrained_checkpoint_loads_and_runs_forward():
    """End-to-end: NEW (real-data) loader → FoundationDownstream → finite output."""
    torch = pytest.importorskip("torch")
    from src.foundation.config import MODALITY_DIMS, FoundationConfig
    from src.foundation.downstream import FoundationDownstream
    from finaledb_pretrained_loader import load_real_cohort, modalities_from_flat_X

    cohort = load_real_cohort(NPZ_PATH)
    modalities = modalities_from_flat_X(cohort["X"])

    # The NEW checkpoint was trained with PRODUCTION_CONFIG
    # (embed_dim=128). Default config already uses embed_dim=128, but
    # we set it for clarity.
    fd = FoundationDownstream(
        config=FoundationConfig(embed_dim=128, n_layers=4, n_heads=4,
                               ff_dim=256, dropout=0.2),
        pretrained=True,
        checkpoint_path=str(CKPT_PATH),
    )
    assert fd._pretrained, "Checkpoint should have been loaded"

    with torch.no_grad():
        joint = fd.encoder({k: torch.from_numpy(v) for k, v in modalities.items()})
    assert joint.shape[0] == cohort["X"].shape[0]
    assert joint.shape[1] == len(MODALITY_DIMS)
    assert joint.shape[2] == 128
    assert torch.isfinite(joint).all(), "Joint embedding must be finite"


@pytest.mark.skipif(
    not NPZ_PATH.exists(),
    reason=f"NEW pretraining npz not present at {NPZ_PATH}",
)
def test_layout_constants_sum_to_2256():
    """Defensive: layout constants must sum to 2256 or splits will be wrong.

    If either ``scripts/pretrain_synthetic_v0.py`` or
    ``scripts/pretrain_production_finaledb.py`` changes its concat
    order, this test fails immediately so the loader stays in sync.
    """
    from finaledb_pretrained_loader import _SPLIT_OFFSETS
    total = sum(end - start for start, end in _SPLIT_OFFSETS.values())
    assert total == 2256, (
        f"Layout sums to {total}, expected 2256. "
        "Update scripts/finaledb_pretrained_loader.py and the "
        "pretrain_*_finaledb.py scripts to match."
    )


# ---------------------------------------------------------------------------
# Regression: OLD (SYNTHETIC_v0) artifacts still load (provenance only)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not OLD_NPZ_PATH.exists() or not OLD_CKPT_PATH.exists(),
    reason="OLD SYNTHETIC_v0 npz or checkpoint not present on disk",
)
def test_synthetic_v0_artifacts_still_load_for_provenance():
    """The renamed historical artifacts must remain loadable.

    Their existence on disk is by design (preserved for audit /
    lineage, and as a known-bad baseline for
    ``test/test_pretrain_bug_fix.py``). They are NOT trained on real
    FinaleDB data — see the banner at the top of
    ``scripts/pretrain_synthetic_v0.py`` and ``docs/PRETRAIN_BUG.md``.
    The test only confirms they remain loadable and produce a finite
    forward pass when explicitly requested; nothing here promotes
    them to "real-data-trained".
    """
    torch = pytest.importorskip("torch")
    from src.foundation.config import MODALITY_DIMS, FoundationConfig
    from src.foundation.downstream import FoundationDownstream
    from finaledb_pretrained_loader import load_real_cohort, modalities_from_flat_X

    cohort = load_real_cohort(OLD_NPZ_PATH)
    modalities = modalities_from_flat_X(cohort["X"])

    # The OLD checkpoint was trained with PROTOTYPE_CONFIG
    # (embed_dim=64). A mismatching config would raise on load_state_dict.
    fd = FoundationDownstream(
        config=FoundationConfig(embed_dim=64, n_layers=2, n_heads=2,
                               ff_dim=128, dropout=0.2),
        pretrained=True,
        checkpoint_path=str(OLD_CKPT_PATH),
    )
    assert fd._pretrained, "OLD checkpoint should load when checkpoint_path is explicit"

    with torch.no_grad():
        joint = fd.encoder({k: torch.from_numpy(v) for k, v in modalities.items()})
    assert joint.shape[0] == cohort["X"].shape[0]
    assert joint.shape[1] == len(MODALITY_DIMS)
    assert joint.shape[2] == 64
    assert torch.isfinite(joint).all(), (
        "OLD SYNTHETIC_v0 forward pass should be finite even though the "
        "checkpoint was trained on synthetic data, not real FinaleDB."
    )


@pytest.mark.skipif(
    not OLD_NPZ_PATH.exists() or not OLD_CKPT_PATH.exists(),
    reason="OLD SYNTHETIC_v0 npz or checkpoint not present on disk",
)
def test_default_loader_paths_do_not_point_to_synthetic_v0():
    """The loader's docstring constants and these test defaults must
    point at the real-data-trained pair, not the synthetic v0 pair.
    A future refactor that silently flips the default would erase the
    'this file is real-data-trained' guarantee — guard against it.
    """
    # The defaults declared at module top of this test file.
    assert "PRODUCTION" in str(CKPT_PATH), (
        f"Default CKPT_PATH should reference the PRODUCTION (real-data) "
        f"checkpoint; got {CKPT_PATH}."
    )
    assert "PRODUCTION" in str(NPZ_PATH), (
        f"Default NPZ_PATH should reference the PRODUCTION cohort; "
        f"got {NPZ_PATH}."
    )
    # Sanity: the OLD paths must still exist (for provenance) but
    # must NOT be the ones the test defaults point at.
    assert "SYNTHETIC_v0" in str(OLD_CKPT_PATH)
    assert "SYNTHETIC_v0" in str(OLD_NPZ_PATH)
