"""Tests for scripts/plot_figures.py.

Three fast tests + smoke tests covering:
1. --help exits 0 quickly (no matplotlib work / no JSON load).
2. --validate-only runs end-to-end, asserts no exceptions, and produces
   no files in docs/figures/ (it writes to a temp dir and discards).
3. Each figure's source data table is present in the JSON the script
   reads (so a missing key would fail loudly, not silently produce a
   blank chart).
4. Smoke test: generating each figure individually produces a non-empty
   PNG <500 KB (disk-pressure budget from the brief).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "plot_figures.py"
POOLED_JSON = REPO_ROOT / "results" / "cross_study_finallydb.json"
PER_CANCER_JSON = REPO_ROOT / "results" / "per_cancer_sens_at_spec.json"
FIG_DIR = REPO_ROOT / "docs" / "figures"

# Import the module-level functions so we don't have to shell out for
# the in-process unit tests. sys.path bootstrap keeps this test runnable
# from any cwd (per the cfdna-fragmentomics skill's hardcoded-paths
# anti-pattern).
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import plot_figures as pf  # noqa: E402

TOP5_CANCERS = ["LUAD", "BRCA", "OV", "PAAD", "HCC_J"]
MAX_PNG_BYTES = 500 * 1024  # 500 KB disk-pressure cap from the brief


# ---------------------------------------------------------------------------
# Test 1: --help is fast and exits 0
# ---------------------------------------------------------------------------
def test_help_is_fast_and_exits_zero() -> None:
    """--help must complete in <10s without loading JSONs or building
    figures. Catches the 'module-level work runs before argparse' bug
    class from the cfdna-fragmentomics skill."""
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        capture_output=True,
        text=True,
        timeout=10,
        cwd=REPO_ROOT,
    )
    elapsed = time.time() - t0
    assert result.returncode == 0, f"--help exit={result.returncode}\n{result.stderr}"
    assert elapsed < 10, f"--help took {elapsed:.2f}s (should be <10s)"
    assert "Generate the four standard cfDNA validation figures" in result.stdout
    assert "--validate-only" in result.stdout
    # Must NOT have done any data loading or figure work.
    assert "Per-seed distribution" not in result.stdout


# ---------------------------------------------------------------------------
# Test 2: --validate-only runs without exception and writes no docs/figures/
# ---------------------------------------------------------------------------
def test_validate_only_runs_clean_and_writes_nothing() -> None:
    """`--validate-only` must build all 4 figures in a temp dir, succeed,
    and discard the temp files (no pollution of docs/figures/)."""
    # Pre-flight: docs/figures/ contents (snapshot for comparison).
    before = set()
    if FIG_DIR.exists():
        before = {p.name for p in FIG_DIR.iterdir() if p.is_file()}

    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--validate-only"],
        capture_output=True,
        text=True,
        timeout=120,  # matplotlib + JSON load on 627 samples; plenty of headroom
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, (
        f"--validate-only failed exit={result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "validate-only: built all 4 figures" in result.stdout

    # The validation path is supposed to use a temp dir and discard it.
    # docs/figures/ should be UNCHANGED (no new files).
    after = set()
    if FIG_DIR.exists():
        after = {p.name for p in FIG_DIR.iterdir() if p.is_file()}
    assert before == after, (
        f"--validate-only must not modify docs/figures/. "
        f"Before: {sorted(before)}\nAfter:  {sorted(after)}"
    )


# ---------------------------------------------------------------------------
# Test 3: each figure's source fields are present in the JSON
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "field, json_path",
    [
        ("pooled.harmonized.per_seed_auc", POOLED_JSON),
        ("pooled.harmonized.auc_mean", POOLED_JSON),
        ("pooled.harmonized.sens_at_spec_95", POOLED_JSON),
        ("pooled.harmonized.sens_at_spec_98", POOLED_JSON),
        ("pooled.harmonized.sens_at_spec_99", POOLED_JSON),
        ("pooled.no_harmonize.per_seed_auc", POOLED_JSON),
        ("pooled.no_harmonize.sens_at_spec_95", POOLED_JSON),
        ("cohort.n_cancer_in_labels", POOLED_JSON),
        ("cohort.n_healthy_in_labels", POOLED_JSON),
        ("cohort.n_total_in_labels", POOLED_JSON),
        ("per_cancer.LUAD.per_seed_auc", POOLED_JSON),
        ("per_cancer.HCC_J.per_seed_auc", POOLED_JSON),
        ("per_cancer.LUAD.auc_mean", POOLED_JSON),
        ("per_cancer.HCC_J.auc_std", POOLED_JSON),
        ("per_cancer.LUAD.sens_at_95", PER_CANCER_JSON),
        ("per_cancer.HCC_J.sens_at_99", PER_CANCER_JSON),
        ("per_cancer.BRCA.n_pos", PER_CANCER_JSON),
        ("per_cancer.OV.n", PER_CANCER_JSON),
    ],
)
def test_source_field_present(field: str, json_path: Path) -> None:
    """Walk a dotted path into a JSON and assert the leaf exists.

    This catches the silent-failure class where the plotting script
    references a JSON field that isn't actually there (would render a
    blank chart with KeyError swallowed somewhere downstream)."""
    with json_path.open() as f:
        blob = json.load(f)

    node = blob
    for part in field.split("."):
        assert isinstance(node, dict), (
            f"{field}: expected dict at '{part}', got {type(node).__name__}"
        )
        assert part in node, (
            f"{field}: missing key '{part}' in {json_path.name}. "
            f"Available keys at this level: {list(node.keys())[:10]}"
        )
        node = node[part]

    # Leaf must be non-None (a list, scalar, or dict).
    assert node is not None, f"{field} resolved to None"


def test_top5_cancers_all_present() -> None:
    """All five top-5 cancers referenced by fig2 / fig4 are present in
    BOTH JSONs the script consumes."""
    with POOLED_JSON.open() as f:
        cs = json.load(f)
    with PER_CANCER_JSON.open() as f:
        pc = json.load(f)
    for cancer in TOP5_CANCERS:
        assert cancer in cs["per_cancer"], (
            f"{cancer} missing from {POOLED_JSON.name}:per_cancer"
        )
        assert cancer in pc["per_cancer"], (
            f"{cancer} missing from {PER_CANCER_JSON.name}:per_cancer"
        )
        # Required sub-fields for fig2 / fig4
        for required in ("per_seed_auc", "auc_mean", "auc_std"):
            assert required in cs["per_cancer"][cancer], (
                f"{cancer} missing field '{required}' in {POOLED_JSON.name}"
            )
        for required in ("sens_at_95", "sens_at_98", "sens_at_99"):
            assert required in pc["per_cancer"][cancer], (
                f"{cancer} missing field '{required}' in {PER_CANCER_JSON.name}"
            )


# ---------------------------------------------------------------------------
# Test 4: smoke — generate each figure in a temp dir, check PNG < 500 KB
# ---------------------------------------------------------------------------
def _make_figure_in_tmp(make_fn, name: str, tmp_dir: Path) -> Path:
    """Helper: call make_fn(out_path) and return the resulting file."""
    out_path = tmp_dir / f"{name}.png"
    make_fn(out_path)
    assert out_path.exists(), f"{name}: file not created at {out_path}"
    size = out_path.stat().st_size
    assert size > 0, f"{name}: empty file"
    assert size < MAX_PNG_BYTES, (
        f"{name}: PNG is {size / 1024:.1f} KB, exceeds the {MAX_PNG_BYTES // 1024} "
        f"KB budget"
    )
    return out_path


def test_smoke_all_figures_under_size_cap() -> None:
    """Generate all 4 figures in a temp dir, assert each PNG is
    non-empty and under 500 KB. Mirrors the brief's disk-pressure cap."""
    pooled = pf.load_pooled()
    per_cancer = pf.load_per_cancer()

    with tempfile.TemporaryDirectory(prefix="plot_figures_smoke_") as tmp:
        tmp_dir = Path(tmp)
        _make_figure_in_tmp(
            lambda p: pf.figure1_pooled(pooled, p), "fig1_pooled_roc", tmp_dir
        )
        _make_figure_in_tmp(
            lambda p: pf.figure2_per_cancer(per_cancer, p),
            "fig2_per_cancer_roc",
            tmp_dir,
        )
        _make_figure_in_tmp(
            lambda p: pf.figure3_calibration(pooled, p), "fig3_calibration", tmp_dir
        )
        _make_figure_in_tmp(
            lambda p: pf.figure4_sens_operating(per_cancer, p),
            "fig4_sens_operating_points",
            tmp_dir,
        )


def test_loaders_are_idempotent() -> None:
    """Loading the same JSON twice returns equal structures (catches
    accidental in-place mutation of the parsed dict)."""
    pooled_a = pf.load_pooled()
    pooled_b = pf.load_pooled()
    assert pooled_a["harmonized"]["auc_mean"] == pooled_b["harmonized"]["auc_mean"]
    assert pooled_a["harmonized"]["per_seed_auc"] == pooled_b["harmonized"]["per_seed_auc"]
    pc_a = pf.load_per_cancer()
    pc_b = pf.load_per_cancer()
    for cancer in TOP5_CANCERS:
        assert pc_a[cancer]["per_seed_auc"] == pc_b[cancer]["per_seed_auc"]
        assert pc_a[cancer]["auc_mean"] == pc_b[cancer]["auc_mean"]
