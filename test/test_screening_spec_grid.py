"""The sens@spec grid must be a flat sequence of floats.

What went wrong
---------------
``scripts/cross_study_finallydb.py`` defines the screening-grade grid as::

    SCREENING_SPEC_GRID: tuple[float, ...] = (
        tuple(DEFAULT_SPECIFICITIES) + tuple(SCREENING_SPECIFICITIES),
    )

The trailing comma makes that a **one-element tuple containing** the
grid, not the grid. With ``--include-screening`` the flag then handed
``((0.95, 0.98, 0.99, 0.995, 0.999),)`` down to the per-cancer
standalone writer, which iterated over it and passed a tuple where a
float specificity was expected. The failure surfaced deep in the stack
as::

    TypeError: unsupported operand type(s) for -: 'float' and 'tuple'

which named neither the flag nor the grid, and it broke exactly the
combination the reproduction driver needs
(``--include-screening`` + ``--out-per-cancer-json``).

These tests pin the shape of both grids and the error message the
normaliser now raises.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_cross_study():
    spec = importlib.util.spec_from_file_location(
        "cross_study_finallydb",
        _REPO_ROOT / "scripts" / "cross_study_finallydb.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_cs = _load_cross_study()


def _all_floats(grid):
    return all(isinstance(v, float) for v in grid)


def test_default_spec_grid_is_flat():
    grid = _cs.DEFAULT_SPEC_GRID
    assert _all_floats(grid), (
        f"DEFAULT_SPEC_GRID must be a flat tuple of floats, got {grid!r}"
    )
    assert grid == (0.95, 0.98, 0.99)


def test_screening_spec_grid_is_flat():
    """The regression: a trailing comma made this a nested 1-tuple."""
    grid = _cs.SCREENING_SPEC_GRID
    assert _all_floats(grid), (
        "SCREENING_SPEC_GRID must be a flat tuple of floats, got "
        f"{grid!r}. A trailing comma after the grid expression turns it "
        "into a one-element tuple holding the grid, which crashes the "
        "per-cancer standalone writer."
    )
    assert grid == (0.95, 0.98, 0.99, 0.995, 0.999)


def test_screening_grid_extends_the_default_grid():
    """The extension appends; it never replaces the legacy points."""
    assert _cs.SCREENING_SPEC_GRID[:len(_cs.DEFAULT_SPEC_GRID)] == \
        _cs.DEFAULT_SPEC_GRID


def test_normaliser_accepts_a_flat_grid():
    assert _cs._normalise_spec_grid([0.95, 0.99]) == (0.95, 0.99)
    assert _cs._normalise_spec_grid((0.95, 0.99)) == (0.95, 0.99)


def test_normaliser_falls_back_to_the_default_grid():
    assert _cs._normalise_spec_grid(None) == _cs.DEFAULT_SPEC_GRID
    assert _cs._normalise_spec_grid([]) == _cs.DEFAULT_SPEC_GRID


def test_normaliser_rejects_a_nested_grid_with_a_readable_message():
    nested = ((0.95, 0.98, 0.99),)
    with pytest.raises(TypeError) as exc:
        _cs._normalise_spec_grid(nested)
    message = str(exc.value)
    assert "flat sequence of floats" in message
    assert "nested" in message


# ── the provenance block must describe the grid that was actually used ──


def _synthetic_ovr():
    """A minimal OvR artifact map: 20 cancer + 20 healthy, one cancer."""
    import numpy as np
    rng = np.random.default_rng(0)
    y = np.array([1] * 20 + [0] * 20, dtype=np.int64)
    score = np.concatenate([
        rng.normal(1.0, 0.4, 20), rng.normal(-1.0, 0.4, 20)])
    study = np.array(["s"] * 40, dtype=object)
    return {"LUAD": {"y_true": y, "score": score, "study_sub": study,
                      "n_cancer": 20, "n_healthy": 20, "skipped": False}}


def test_provenance_reports_the_operative_grid(tmp_path):
    """Before the fix, `provenance.specificities` always said
    [0.95, 0.98, 0.99] even when the screening grid was in use, so the
    artifact misdescribed its own contents."""
    out = tmp_path / "per_cancer.json"
    payload = _cs.build_per_cancer_standalone_payload(
        _synthetic_ovr(),
        str(out),
        pooled_harmonized={"auc_mean": 0.9, "auc_std": 0.0,
                           "sens_at_spec_95": 0.8, "sens_at_spec_98": 0.7,
                           "sens_at_spec_99": 0.6},
        pooled_n_pos=20, pooled_n_neg=20,
        specificities_grid=_cs.SCREENING_SPEC_GRID,
    )
    assert payload["config"]["specificities"] == list(
        _cs.SCREENING_SPEC_GRID)
    assert payload["provenance"]["specificities"] == list(
        _cs.SCREENING_SPEC_GRID)
    row = payload["per_cancer"]["LUAD"]["sens_at_spec"]
    assert set(row) == {"0.95", "0.98", "0.99", "0.995", "0.999"}


def test_provenance_reports_the_legacy_grid_when_screening_is_off(tmp_path):
    out = tmp_path / "per_cancer.json"
    payload = _cs.build_per_cancer_standalone_payload(
        _synthetic_ovr(),
        str(out),
        pooled_harmonized={"auc_mean": 0.9, "auc_std": 0.0,
                           "sens_at_spec_95": 0.8, "sens_at_spec_98": 0.7,
                           "sens_at_spec_99": 0.6},
        pooled_n_pos=20, pooled_n_neg=20,
        specificities_grid=_cs.DEFAULT_SPEC_GRID,
    )
    assert payload["config"]["specificities"] == [0.95, 0.98, 0.99]
    assert payload["provenance"]["specificities"] == [0.95, 0.98, 0.99]
    row = payload["per_cancer"]["LUAD"]["sens_at_spec"]
    assert set(row) == {"0.95", "0.98", "0.99"}
