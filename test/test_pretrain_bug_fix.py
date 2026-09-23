"""Regression tests for the real-modalities bypass bug in
``src/foundation.pretrain.FoundationPretrainer``.

Bug summary
-----------
Before the fix (pre-2026-09-23), every pretrain phase called
``self.data_generator.generate_dataset(...)`` unconditionally and
silently ignored any real modality dict the caller had supplied.
The fix added a ``use_real_modalities`` flag (default ``True``)
plus per-phase ``modalities`` overrides. When the flag is on, the
phase must draw batches from the real cohort and MUST NOT call
``self.data_generator.generate_dataset``.

These tests pin the contract:

1. With ``use_real_modalities=True`` and a real cohort supplied, the
   synthetic generator's ``generate_dataset`` is NEVER invoked for
   any of the three phases (count stays at zero).
2. With ``use_real_modalities=True`` and a real cohort supplied, the
   data that the encoder actually sees is the real cohort (not
   random hash-derived synthetic data). We assert equality against
   the input modality arrays after one epoch.
3. With ``use_real_modalities=False`` (the synthetic fallback), the
   synthetic generator IS invoked — backwards-compat for the
   CI / unit-test path that has no real cohort.
4. The constructor ``modalities=`` kwarg is honored when per-phase
   overrides are not supplied.
5. Per-phase ``modalities`` overrides beat the constructor default.

The tests run on CPU with the PROTOTYPE_CONFIG for speed.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from unittest.mock import patch as mock_patch

import numpy as np
import pytest
import torch

# Make deepcatch importable regardless of pytest invocation cwd.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.foundation.config import MODALITY_DIMS, FoundationConfig, PROTOTYPE_CONFIG
from src.foundation.data import MultiModalDataGenerator
from src.foundation.pretrain import FoundationPretrainer


def _make_real_modalities(n_samples: int = 16, seed: int = 42) -> dict:
    """Build a small real-cohort dict using a real numpy seed.

    Crucially we use a fresh ``np.random.default_rng`` (not the
    ``MultiModalDataGenerator`` so the values are not synthetic
    hash-derived and so we can detect any accidental fallback.
    """
    rng = np.random.default_rng(seed)
    return {
        name: rng.standard_normal((n_samples, dim)).astype(np.float32)
        for name, dim in MODALITY_DIMS.items()
    }


class _CallCounter:
    """Lightweight wrapper that counts calls to a bound method."""

    def __init__(self, obj, attr: str):
        self.obj = obj
        self.attr = attr
        self.calls = 0
        self.original = getattr(obj, attr)
        # Capture the call so we don't lose data when the test
        # exercises the synthetic path (the inner code calls
        # ``self.data_generator.generate_dataset``).

    def __call__(self, *args, **kwargs):
        self.calls += 1
        return self.original(*args, **kwargs)


def _attach_counter(pretrainer: FoundationPretrainer) -> _CallCounter:
    counter = _CallCounter(pretrainer.data_generator, "generate_dataset")
    pretrainer.data_generator.generate_dataset = counter
    return counter


# ── 1. Real-data path does not call the synthetic generator ──────────


def test_phase1_with_real_modalities_skips_synthetic_generator():
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=16, seed=1)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    losses = pretrainer.pretrain_phase1_mmp(
        n_samples=64, n_epochs=2, batch_size=8,
    )

    assert counter.calls == 0, (
        f"generate_dataset called {counter.calls} times during phase 1 "
        "with use_real_modalities=True; expected 0."
    )
    assert len(losses) == 2
    assert all(np.isfinite(l) for l in losses)


def test_phase2_with_real_modalities_skips_synthetic_generator():
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=16, seed=2)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    losses = pretrainer.pretrain_phase2_contrastive(
        n_samples=64, n_epochs=2, batch_size=8,
    )

    assert counter.calls == 0, (
        f"generate_dataset called {counter.calls} times during phase 2 "
        "with use_real_modalities=True; expected 0."
    )
    assert len(losses) == 2
    assert all(np.isfinite(l) for l in losses)


def test_phase3_with_real_modalities_skips_synthetic_generator():
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=16, seed=3)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    losses = pretrainer.pretrain_phase3_joint(
        n_samples=64, n_epochs=2, batch_size=8,
    )

    assert counter.calls == 0, (
        f"generate_dataset called {counter.calls} times during phase 3 "
        "with use_real_modalities=True; expected 0."
    )
    assert len(losses) == 2
    assert all(np.isfinite(l) for l in losses)


# ── 2. Encoder actually sees the real data, not synthetic ────────────


def test_encoder_input_is_real_modalities_when_flag_is_true():
    """Patches encoder.forward to capture the batch it received and
    asserts the captured batch is the real cohort after per-modality
    standardization (which the pretrainer applies in every phase —
    see ``FoundationPretrainer._fit_modality_stats``)."""
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=8, seed=4)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )

    captured = {}
    real_forward = pretrainer.encoder.forward

    def spy_forward(modalities, mask=None):
        # Capture a clone of the first sample of each modality — that's
        # enough to verify provenance. We capture before the encoder
        # applies any in-place ops.
        captured["batch0"] = {
            k: v[0].detach().cpu().clone() for k, v in modalities.items()
        }
        return real_forward(modalities, mask=mask)

    pretrainer.encoder.forward = spy_forward

    pretrainer.pretrain_phase1_mmp(
        n_samples=8, n_epochs=1, batch_size=8,
    )

    assert set(captured["batch0"].keys()) == set(real.keys()), (
        f"Encoder saw keys {set(captured['batch0'].keys())}, "
        f"expected {set(real.keys())}"
    )
    # The pretrainer now standardizes each modality (median / MAD-robust-
    # std) before the encoder. Reconstruct the expected standardized
    # first sample from the raw cohort and compare.
    from src.foundation.pretrain import (
        _fit_modality_stats, _apply_modality_standardization,
    )
    stats = _fit_modality_stats(real)
    expected_std = _apply_modality_standardization(real, stats)
    for k, expected in expected_std.items():
        got = captured["batch0"][k].numpy()
        # The encoder receives the cohort on its device. Allow float
        # comparison tolerance — values are passed through unchanged
        # after standardization.
        np.testing.assert_allclose(
            got, expected[0], rtol=1e-5, atol=1e-6,
            err_msg=(
                f"Modality {k}: encoder input does not match the "
                "standardized real cohort — pretrainer is sampling "
                "synthetic data even though use_real_modalities=True."
            ),
        )


# ── 3. Synthetic fallback still works when flag is False ────────────


def test_synthetic_path_still_calls_generator_when_flag_is_false():
    """Backwards-compat: callers that explicitly opt out
    (CI / unit tests with no real cohort) must still hit the
    synthetic generator."""
    cfg = PROTOTYPE_CONFIG
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
    )
    counter = _attach_counter(pretrainer)

    losses = pretrainer.pretrain_phase1_mmp(
        n_samples=50, n_epochs=2, batch_size=10,
    )

    assert counter.calls == 1, (
        f"Expected exactly 1 call to generate_dataset when "
        f"use_real_modalities=False (or no real cohort supplied); "
        f"got {counter.calls}."
    )
    assert len(losses) == 2


def test_synthetic_fallback_when_flag_true_but_no_modalities():
    """If the caller says use_real_modalities=True but forgets to
    pass modalities, the pretrainer must fall back to synthetic
    (and emit a warning) instead of crashing."""
    cfg = PROTOTYPE_CONFIG
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=None, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    with mock_patch(
        "src.foundation.pretrain.logger.warning"
    ) as warn_spy:
        losses = pretrainer.pretrain_phase1_mmp(
            n_samples=20, n_epochs=1, batch_size=10,
        )
    assert warn_spy.called, (
        "Expected a warning when use_real_modalities=True but no "
        "modalities were supplied."
    )
    assert counter.calls == 1
    assert len(losses) == 1


# ── 4. Constructor-level modalities honored when per-phase absent ──


def test_constructor_modalities_used_by_default_in_phase1():
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=8, seed=5)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    # No per-phase modalities / use_real_modalities overrides.
    pretrainer.pretrain_phase1_mmp(
        n_samples=64, n_epochs=1, batch_size=8,
    )
    assert counter.calls == 0


# ── 5. Per-phase override beats the constructor default ─────────────


def test_per_phase_modalities_override_constructor():
    """The per-phase ``modalities=`` argument must override the
    constructor-level default."""
    cfg = PROTOTYPE_CONFIG
    real_at_construction = _make_real_modalities(n_samples=8, seed=6)
    real_at_phase = _make_real_modalities(n_samples=8, seed=7)

    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real_at_construction, use_real_modalities=True,
    )
    captured = {}
    real_forward = pretrainer.encoder.forward

    def spy_forward(modalities, mask=None):
        captured["batch0"] = {
            k: v[0].detach().cpu().clone() for k, v in modalities.items()
        }
        return real_forward(modalities, mask=mask)

    pretrainer.encoder.forward = spy_forward

    pretrainer.pretrain_phase1_mmp(
        n_samples=8, n_epochs=1, batch_size=8,
        modalities=real_at_phase,  # per-phase override
    )

    # Verify the encoder saw the per-phase cohort (after the per-modality
    # standardization the pretrainer runs on every phase), not the
    # constructor-level one.
    from src.foundation.pretrain import (
        _fit_modality_stats, _apply_modality_standardization,
    )
    expected_stats = _fit_modality_stats(real_at_phase)
    expected_std = _apply_modality_standardization(
        real_at_phase, expected_stats,
    )
    for k, expected in expected_std.items():
        got = captured["batch0"][k].numpy()
        np.testing.assert_allclose(
            got, expected[0], rtol=1e-5, atol=1e-6,
            err_msg=f"Modality {k}: per-phase override not honored.",
        )


def test_per_phase_flag_false_overrides_constructor_true():
    """Per-phase ``use_real_modalities=False`` must override the
    constructor-level True default — this lets a caller mix-and-
    match: real-data phase 1, synthetic phase 2, etc."""
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=8, seed=8)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    pretrainer.pretrain_phase1_mmp(
        n_samples=20, n_epochs=1, batch_size=10,
        use_real_modalities=False,
    )
    assert counter.calls == 1


# ── 6. Sanity: full pipeline runs with real data ─────────────────────


def test_full_pretrain_pipeline_with_real_modalities():
    """End-to-end: all 3 phases, real cohort, no synthetic calls."""
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=12, seed=9)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    counter = _attach_counter(pretrainer)

    result = pretrainer.pretrain(
        n_samples=64,
        p1_epochs=2, p2_epochs=2, p3_epochs=1,
        batch_size=8,
    )

    assert counter.calls == 0
    assert set(result.keys()) == {"phase1", "phase2", "phase3"}
    assert len(result["phase1"]) == 2
    assert len(result["phase2"]) == 2
    assert len(result["phase3"]) == 1
    for phase_losses in result.values():
        assert all(np.isfinite(l) for l in phase_losses)
    assert pretrainer.is_pretrained


# ── 7. Encoder output shape / finiteness on real cohort ──────────────


def test_encoder_forward_on_real_modalities_is_finite():
    """Smoke test: encoder consumes the real cohort end-to-end and
    produces finite outputs in the expected shape."""
    cfg = PROTOTYPE_CONFIG
    real = _make_real_modalities(n_samples=6, seed=10)
    pretrainer = FoundationPretrainer(
        config=cfg, device="cpu", verbose=False,
        modalities=real, use_real_modalities=True,
    )
    # Quick 1-epoch training so weights are non-trivial.
    pretrainer.pretrain_phase1_mmp(
        n_samples=6, n_epochs=1, batch_size=6,
    )

    pretrainer.encoder.eval()
    with torch.no_grad():
        joint = pretrainer.encoder({
            k: torch.from_numpy(v) for k, v in real.items()
        })
    assert joint.shape == (6, len(MODALITY_DIMS), cfg.embed_dim)
    assert torch.isfinite(joint).all(), (
        "Encoder output on real cohort contains NaN/Inf after "
        "use_real_modalities=True training."
    )
