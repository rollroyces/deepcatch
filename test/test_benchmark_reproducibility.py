"""Tests that the benchmark pipeline is reproducible.

Why this file exists
---------------------
The committed cross-study artifacts were produced by a pipeline whose PCA
was fitted WITHOUT ``random_state``. sklearn picks the randomised solver
when ``n_components < 0.8 * min(n_samples, n_features)`` -- true for
every fold of this pipeline -- and the randomised solver is seeded from
``random_state``. The practical effect was that re-running the benchmark
produced a slightly different number every time (one publication's AUC
moved 0.9789 -> 0.9821), so "the numbers in the paper" could not be
reproduced from the code in the repository.

Two guards are installed here:

* a source-level guard that fails if any ``PCA(...)`` call in the
  repository's own ``scripts/`` and ``src/`` is fitted without a
  ``random_state`` -- this one always runs, including on CI;
* a behavioural guard that fits the benchmark twice and requires the two
  runs to agree exactly. It needs the sibling pipeline checkout, so it
  is marked and skipped where that is absent.
"""
from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

_SCAN_DIRS = ("scripts", "src")


def _is_pca_call(node: ast.Call) -> bool:
    func = node.func
    if isinstance(func, ast.Name):
        return func.id == "PCA"
    if isinstance(func, ast.Attribute):
        return func.attr == "PCA"
    return False


def _pca_calls(path: Path):
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and _is_pca_call(node):
            yield node


def test_every_pca_call_is_seeded():
    """No unseeded PCA in the pipeline code.

    An unseeded PCA makes every artifact in the repository that depends
    on it irreproducible, and the drift is small enough to pass a
    plausibility check while still changing a published number.
    """
    offenders = []
    for d in _SCAN_DIRS:
        for path in sorted((_REPO_ROOT / d).rglob("*.py")):
            for node in _pca_calls(path):
                if not any(kw.arg == "random_state" for kw in node.keywords):
                    offenders.append(
                        f"{path.relative_to(_REPO_ROOT)}:{node.lineno}")
    assert not offenders, (
        "PCA is fitted without random_state (sklearn then uses the "
        "randomised solver, so the fit is not reproducible): "
        + ", ".join(offenders))


def _load_cross_study():
    spec = importlib.util.spec_from_file_location(
        "cross_study_finallydb",
        _REPO_ROOT / "scripts" / "cross_study_finallydb.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.requires_sister_pipeline
def test_pooled_oof_is_reproducible():
    """Two runs of the benchmark's fold fit must agree exactly.

    The component count is chosen so that sklearn uses the randomised
    solver, which is exactly the case an unseeded fit gets wrong.
    """
    mod = _load_cross_study()
    rng = np.random.default_rng(0)
    n = 160
    X = rng.normal(size=(n, 60))
    y = (rng.random(n) > 0.5).astype(int)
    st = np.array(["jiang"] * (n // 2) + ["cristiano"] * (n - n // 2))

    y1, s1, a1 = mod.pooled_oof(X, y, st, [42, 13], 12, False)
    y2, s2, a2 = mod.pooled_oof(X, y, st, [42, 13], 12, False)

    assert np.array_equal(y1, y2)
    assert np.array_equal(s1, s2), (
        "the pooled out-of-fold scores differ between two identical runs: "
        f"max |delta| = {np.abs(s1 - s2).max():.3e}")
    assert a1 == a2, f"per-seed AUCs differ between runs: {a1} vs {a2}"


@pytest.mark.requires_sister_pipeline
def test_pooled_oof_differs_when_the_seed_changes():
    """Guard against the opposite failure: a seed that pins everything.

    If the PCA seed were derived from the data instead of the fold
    seed, two different fold splits would produce identical scores.
    """
    mod = _load_cross_study()
    rng = np.random.default_rng(0)
    n = 160
    X = rng.normal(size=(n, 60))
    y = (rng.random(n) > 0.5).astype(int)
    st = np.array(["jiang"] * (n // 2) + ["cristiano"] * (n - n // 2))

    _, s_a, _ = mod.pooled_oof(X, y, st, [42], 12, False)
    _, s_b, _ = mod.pooled_oof(X, y, st, [7], 12, False)
    assert not np.array_equal(s_a, s_b)
