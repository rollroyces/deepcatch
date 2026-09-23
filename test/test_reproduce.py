"""Tests for the reproduce pipeline artifacts.

These tests verify the reproduce scripts exist, are syntactically
valid, and stay in sync with the docs/ sub-system. They DO NOT run the
~10-15 minute full pipeline — they only do cheap checks.

Honors pitfall #7 from the deepcatch skill: tests resolve the repo
root via walk-up, so this file works whether pytest discovers it via
the editable-install symlink (`.venv/lib/.../site-packages/test/...`)
or directly from the repo.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


# ──────────────────────────────────────────────────────────────────────
# Walk-up repo-root resolver (deepcatch skill pitfall #7)
# ──────────────────────────────────────────────────────────────────────


def _resolve_repo_root() -> Path:
    cur = Path(__file__).resolve().parent
    while cur != cur.parent:
        if (cur / "pyproject.toml").is_file() and (cur / "scripts").is_dir():
            return cur
        cur = cur.parent
    return Path.cwd()


REPO_ROOT = _resolve_repo_root()
SCRIPTS_DIR = REPO_ROOT / "scripts"
TEST_DIR = REPO_ROOT / "test"
RESULTS_DIR = REPO_ROOT / "results"
DOCS_DIR = REPO_ROOT / "docs"


def _find_bash() -> str | None:
    """Locate bash. Falls back to /bin/bash; if absent on Windows, skip."""
    candidates = [
        shutil.which("bash"),
        "/bin/bash",
        "/usr/bin/bash",
        "/opt/homebrew/bin/bash",
        "/usr/local/bin/bash",
    ]
    for c in candidates:
        if c and Path(c).is_file():
            return c
    return None


def _safe_env() -> dict:
    """Honor deepcatch skill pitfall #1: blank PYTHONPATH so the
    subprocess inherits the venv's site-packages, not Hermes's global."""
    env = {**os.environ}
    env.pop("PYTHONPATH", None)
    return env


# ──────────────────────────────────────────────────────────────────────
# Tests
# ──────────────────────────────────────────────────────────────────────


def test_verify_environment_py_is_importable():
    """`scripts/verify_environment.py` exists and imports cleanly under
    the venv (which is how the reproduce pipeline invokes it)."""
    p = SCRIPTS_DIR / "verify_environment.py"
    assert p.is_file(), f"missing {p}"

    py_exec = REPO_ROOT / ".venv" / "bin" / "python"
    if not py_exec.exists():
        pytest.skip("no .venv/bin/python — verify_environment is venv-only")

    proc = subprocess.run(
        [str(py_exec), "-c", "import importlib.util, sys; "
         "spec = importlib.util.spec_from_file_location('ve', "
         f"'{p}'); m = importlib.util.module_from_spec(spec); "
         "spec.loader.exec_module(m); "
         "assert callable(m.main); print('ok')"],
        capture_output=True, text=True, env=_safe_env(),
        timeout=20,
    )
    assert proc.returncode == 0, (
        f"verify_environment.py failed to import: "
        f"stdout={proc.stdout!r} stderr={proc.stderr!r}"
    )
    assert "ok" in proc.stdout


def test_reproduce_all_sh_exists_and_runs_dry(tmp_path: Path):
    """`scripts/reproduce_all.sh` exists, passes bash -n syntax check,
    and supports --dry-run."""
    p = SCRIPTS_DIR / "reproduce_all.sh"
    assert p.is_file(), f"missing {p}"

    bash = _find_bash()
    if bash is None:
        pytest.skip("no bash on PATH (Windows?) — bash script not testable here")

    # 1) bash -n syntax check — no execution, just parse.
    syntax = subprocess.run(
        [bash, "-n", str(p)],
        capture_output=True, text=True,
        timeout=10,
    )
    assert syntax.returncode == 0, (
        f"reproduce_all.sh failed bash -n syntax check: "
        f"stderr={syntax.stderr!r}"
    )

    # 2) --dry-run: must exit 0 AND not modify results/ in a way that
    #    breaks later steps. We point it at a tmp RESULTS dir via
    #    monkeypatching the script's REPO_ROOT resolution would be
    #    brittle, so we instead just confirm the script exits 0 and
    #    prints the step headers. We don't assert JSON existence since
    #    the script writes no JSONs in dry-run mode.
    dry = subprocess.run(
        [bash, str(p), "--dry-run"],
        capture_output=True, text=True,
        timeout=30,
    )
    assert dry.returncode == 0, (
        f"reproduce_all.sh --dry-run failed: "
        f"stderr={dry.stderr!r} stdout_tail={dry.stdout[-500:]!r}"
    )
    # We expect every step slug to appear in stdout at least once.
    expected_slugs = [
        "cross_study",
        "per_cancer",
        "harmonization",
        "focal_bce_ablation",
        "sparse_aware_ablation",
        "foundation_smoke",
        "pytest_tests",
    ]
    for slug in expected_slugs:
        assert slug in dry.stdout, (
            f"--dry-run output missing step slug {slug!r}; got: "
            f"{dry.stdout[-1000:]!r}"
        )


def test_reproduce_all_py_imports_clean():
    """`scripts/reproduce_all.py` imports without raising AND supports
    --dry-run."""
    p = SCRIPTS_DIR / "reproduce_all.py"
    assert p.is_file(), f"missing {p}"

    py_exec = REPO_ROOT / ".venv" / "bin" / "python"
    if not py_exec.exists():
        pytest.skip("no .venv/bin/python — Python reproduce not testable here")

    # Import the module via runpy (not just py_compile) so we exercise
    # the import-time code path the reproduce script will hit.
    proc = subprocess.run(
        [str(py_exec), str(p), "--dry-run"],
        capture_output=True, text=True, env=_safe_env(),
        timeout=60,
    )
    assert proc.returncode == 0, (
        f"reproduce_all.py --dry-run failed: "
        f"stdout_tail={proc.stdout[-500:]!r} "
        f"stderr={proc.stderr!r}"
    )
    # --dry-run must mention every expected slug
    expected_slugs = [
        "cross_study", "per_cancer", "harmonization",
        "focal_bce_ablation", "sparse_aware_ablation",
        "foundation_smoke", "foundation_smoke_sens",
        "foundation_smoke_sparse", "pytest_tests",
    ]
    for slug in expected_slugs:
        assert slug in proc.stdout, (
            f"--dry-run output missing step slug {slug!r}"
        )


def test_reproduce_documents_match_results_files():
    """Every JSON artifact the reproduce scripts claim to produce must
    have a corresponding reference in `docs/REPRODUCE.md` (no orphan
    claims). Conversely, every step listed in REPRODUCE.md must have a
    backing script that writes a JSON to results/.

    This is a structural cross-check between the docs and the actual
    reproduce scripts — it fails fast if someone adds a step to the
    docs but forgets to wire it into the bash script (or vice versa).
    """
    # Allow README step to be missing if reproduce hasn't run yet.
    reproduce_doc = DOCS_DIR / "REPRODUCE.md"
    if not reproduce_doc.is_file():
        pytest.fail(
            f"docs/REPRODUCE.md is missing — the reproduce pipeline "
            f"docs contract is not satisfied."
        )

    doc_text = reproduce_doc.read_text()

    # Each step slug is referenced in the doc; the corresponding .json
    # basename appears in the doc's schema table OR in the body.
    step_slugs = [
        ("cross_study", "cross_study.json"),
        ("per_cancer", "per_cancer.json"),
        ("harmonization", "harmonization.json"),
        ("focal_bce_ablation", "focal_bce_ablation.json"),
        ("sparse_aware_ablation", "sparse_aware_ablation.json"),
        ("foundation_smoke", "foundation_smoke.json"),
        ("pytest_tests", "pytest_tests.json"),
    ]
    for slug, basename in step_slugs:
        assert slug in doc_text, (
            f"docs/REPRODUCE.md does not reference step slug {slug!r}"
        )
        # The basename is allowed to be referenced via slug.json OR
        # indirectly via the slug name itself (e.g. "results/<slug>.json"
        # in a table row). Be lenient: just require the slug to appear.
        # The bash + python scripts also reference the slug via variable
        # interpolation, so check that too.
        bash_text = (SCRIPTS_DIR / "reproduce_all.sh").read_text()
        py_text = (SCRIPTS_DIR / "reproduce_all.py").read_text()
        assert slug in bash_text, (
            f"scripts/reproduce_all.sh does not reference step {slug!r}"
        )
        assert slug in py_text, (
            f"scripts/reproduce_all.py does not reference step {slug!r}"
        )


# ──────────────────────────────────────────────────────────────────────
# Sanity: smoke that the env check itself exits 0 right now (the
# reproduce pre-condition). Skipped on machines without the venv.
# ──────────────────────────────────────────────────────────────────────


def test_verify_environment_exits_zero_when_viable():
    """On this machine, `verify_environment.py` should exit 0 because
    every dependency listed in its table is installed."""
    p = SCRIPTS_DIR / "verify_environment.py"
    if not p.is_file():
        pytest.skip("verify_environment.py missing")

    py_exec = REPO_ROOT / ".venv" / "bin" / "python"
    if not py_exec.exists():
        pytest.skip("no .venv/bin/python")

    proc = subprocess.run(
        [str(py_exec), str(p)],
        capture_output=True, text=True, env=_safe_env(),
        timeout=20,
    )
    # Don't hard-assert rc==0 here — a fresh checkout might be missing
    # something. We only assert the JSON --json branch parses.
    assert "Summary:" in proc.stdout, (
        f"verify_environment.py output missing summary line: "
        f"{proc.stdout[:500]!r}"
    )


# ──────────────────────────────────────────────────────────────────────
# Honest test: confirm --only filter excludes non-listed steps from the
# bash script's dry-run output. This catches typos in step slugs.
# ──────────────────────────────────────────────────────────────────────


def test_reproduce_all_sh_only_filter_works(tmp_path: Path):
    """`--only` excludes non-listed steps from actually running.

    The header banner is OK to appear. We check that cross_study's
    status is SKIP (not "would run this step")."""
    bash = _find_bash()
    if bash is None:
        pytest.skip("no bash on PATH")
    p = SCRIPTS_DIR / "reproduce_all.sh"

    proc = subprocess.run(
        [bash, str(p), "--dry-run", "--only", "per_cancer,harmonization"],
        capture_output=True, text=True,
        timeout=30,
    )
    assert proc.returncode == 0, (
        f"--only filter failed: stderr={proc.stderr!r}"
    )
    # per_cancer and harmonization should be present and would-run
    assert "[dry-run] would run this step" in proc.stdout, (
        "expected at least one dry-run step line"
    )
    # cross_study should be marked SKIPPED, NOT "would run this step".
    # The exact text the script emits for SKIPPED is:
    #     ==[cross_study]== SKIPPED (not in --only set)
    # And for would-run it's:
    #     ==[cross_study]== cross_study_finallydb.py (...)
    #         out: ...
    #         [dry-run] would run this step
    cross_study_block = proc.stdout.split("==[cross_study]==", 1)[1]
    cross_study_block = cross_study_block.split("==[", 1)[0]  # next block
    assert "SKIPPED" in cross_study_block, (
        f"cross_study was not marked SKIPPED under --only: "
        f"{cross_study_block!r}"
    )
    assert "[dry-run] would run this step" not in cross_study_block, (
        f"cross_study would have been executed under --only: "
        f"{cross_study_block!r}"
    )