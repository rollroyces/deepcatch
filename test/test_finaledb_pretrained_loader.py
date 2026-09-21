"""Tests for the FinaleDB pretrained checkpoint loader.

Verifies the integration chain:
1. The flat 2256-dim feature matrix splits into the expected 6
   per-modality arrays matching ``FoundationDownstream.MODALITY_DIMS``.
2. ``FoundationDownstream(pretrained=True, checkpoint_path=...)``
   loads the saved encoder weights and produces a finite forward pass.

These tests require:
- ``data/finaledb_pretrain_cohort.npz`` (committed)
- ``checkpoints/foundation_pretrained_finaledb.pt`` (gitignored; the
  pretraining pipeline produces it. If absent, tests are skipped.)
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


NPZ_PATH = ROOT / "data" / "finaledb_pretrain_cohort.npz"
CKPT_PATH = ROOT / "checkpoints" / "foundation_pretrained_finaledb.pt"


@pytest.mark.skipif(
    not NPZ_PATH.exists(),
    reason=f"Pretraining npz not present at {NPZ_PATH}",
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
    reason="Pretraining npz or checkpoint not present",
)
def test_pretrained_checkpoint_loads_and_runs_forward():
    """End-to-end: loader → FoundationDownstream(pretrained=True) → finite output."""
    torch = pytest.importorskip("torch")
    from src.foundation.config import MODALITY_DIMS, FoundationConfig
    from src.foundation.downstream import FoundationDownstream
    from finaledb_pretrained_loader import load_real_cohort, modalities_from_flat_X

    cohort = load_real_cohort(NPZ_PATH)
    modalities = modalities_from_flat_X(cohort["X"])

    # The checkpoint was trained with embed_dim=64 (PROTOTYPE_CONFIG).
    # The default config in FoundationDownstream is embed_dim=128, so we
    # must use a matching config to load the weights.
    fd = FoundationDownstream(
        config=FoundationConfig(embed_dim=64, n_layers=2, n_heads=2,
                               ff_dim=128, dropout=0.2),
        pretrained=True,
        checkpoint_path=str(CKPT_PATH),
    )
    assert fd._pretrained, "Checkpoint should have been loaded"

    with torch.no_grad():
        joint = fd.encoder({k: torch.from_numpy(v) for k, v in modalities.items()})
    assert joint.shape[0] == cohort["X"].shape[0]
    assert joint.shape[1] == len(MODALITY_DIMS)
    assert joint.shape[2] == 64
    assert torch.isfinite(joint).all(), "Joint embedding must be finite"


@pytest.mark.skipif(
    not NPZ_PATH.exists(),
    reason=f"Pretraining npz not present at {NPZ_PATH}",
)
def test_layout_constants_sum_to_2256():
    """Defensive: layout constants must sum to 2256 or splits will be wrong.

    If ``scripts/pretrain_real_finaledb.py`` changes its concat order,
    this test fails immediately so the loader stays in sync.
    """
    from finaledb_pretrained_loader import _SPLIT_OFFSETS
    total = sum(end - start for start, end in _SPLIT_OFFSETS.values())
    assert total == 2256, (
        f"Layout sums to {total}, expected 2256. "
        "Update scripts/finaledb_pretrained_loader.py and "
        "scripts/pretrain_real_finaledb.py to match."
    )
