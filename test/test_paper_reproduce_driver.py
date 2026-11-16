"""The paper reproduction driver must cover every artifact the paper quotes.

Why this file exists
--------------------
``paper/REPRODUCE.sh`` re-runs the scripts that produce the numbers in
``paper/METHODS_PAPER.md`` and ``paper/paper.tex``. It used to cover only
the two cross-study sweeps, while the paper also quotes the LOSO
transfer, the external-validation cohort, the fusion ablation, the
shuffled-label null and the methylation channel. A driver that silently
does not regenerate a quoted number is a reproduction driver that lies
about what it reproduces.

The authoritative list of quoted artifacts lives in
``scripts/benchmark_numbers.py`` (``quoted_artifacts()``). These tests
read that list rather than restating it, so the driver and the
quotable-number list cannot drift apart.

The second half of the file tests ``scripts/reproduce_check.py`` itself,
including the red case: a perturbed artifact MUST be reported as drift.
A check that never fails proves nothing.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
DRIVER = _REPO_ROOT / "paper" / "REPRODUCE.sh"
CHECK = _REPO_ROOT / "scripts" / "reproduce_check.py"


def _load_numbers_module():
    spec = importlib.util.spec_from_file_location(
        "benchmark_numbers", _REPO_ROOT / "scripts" / "benchmark_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_numbers = _load_numbers_module()


def _find_bash() -> str | None:
    for c in (shutil.which("bash"), "/bin/bash", "/usr/bin/bash"):
        if c and Path(c).is_file():
            return c
    return None


def _driver_text() -> str:
    assert DRIVER.is_file(), f"missing {DRIVER}"
    return DRIVER.read_text()


def _driver_artifacts() -> set[str]:
    """The artifact names the driver declares in its ARTIFACTS array."""
    text = _driver_text()
    match = re.search(r"ARTIFACTS=\(([^)]*)\)", text)
    assert match, "paper/REPRODUCE.sh has no ARTIFACTS=(...) array"
    return set(re.findall(r"[A-Za-z0-9_.\-]+\.json", match.group(1)))


# ── driver coverage ─────────────────────────────────────────────────────


def test_driver_is_syntactically_valid():
    bash = _find_bash()
    if bash is None:
        pytest.skip("no bash on PATH")
    proc = subprocess.run([bash, "-n", str(DRIVER)],
                          capture_output=True, text=True, timeout=20)
    assert proc.returncode == 0, (
        f"paper/REPRODUCE.sh failed bash -n: {proc.stderr!r}")


def test_driver_regenerates_every_quoted_artifact():
    """Every artifact the documents quote a number from must be in the
    driver's ARTIFACTS list, or the driver does not reproduce the paper."""
    quoted = set(_numbers.quoted_artifacts())
    covered = _driver_artifacts()
    missing = sorted(quoted - covered)
    assert not missing, (
        "paper/REPRODUCE.sh does not regenerate these quoted artifacts: "
        f"{missing}. Add them to the ARTIFACTS array and to a step that "
        f"writes them."
    )


def test_driver_writes_every_covered_artifact():
    """Each artifact in the ARTIFACTS array must actually be produced by
    one of the driver's steps (its filename must appear in a command)."""
    text = _driver_text()
    covered = _driver_artifacts()
    assert covered, "the ARTIFACTS array is empty"
    for name in sorted(covered):
        assert name in text, (
            f"{name} is listed in ARTIFACTS but no step writes it"
        )


def test_driver_runs_the_reproduction_check():
    """The driver must compare the regenerated numbers against the
    committed ones, not merely regenerate them."""
    text = _driver_text()
    assert "reproduce_check.py" in text, (
        "paper/REPRODUCE.sh does not invoke scripts/reproduce_check.py — "
        "it re-runs the scripts but never proves the numbers reproduced"
    )
    assert "--baseline" in text and "--current" in text, (
        "the reproduce_check.py call must pass both a baseline and a "
        "current results directory"
    )


def test_driver_quick_mode_writes_to_a_scratch_directory():
    """--quick must not write into the repo: a reduced-seed run would
    overwrite the published artifacts with different numbers."""
    text = _driver_text()
    assert 'RESULTS_DIR="${SCRATCH}/results"' in text, (
        "in --quick mode the driver must point RESULTS_DIR at the "
        "scratch directory, not at the repo's results/"
    )
    assert "quick mode: the repo was NOT modified" in text, (
        "the driver must state plainly that quick mode leaves the repo "
        "untouched"
    )


def test_driver_accepts_env_overrides_for_the_data_paths():
    """An independent researcher must be able to point the driver at
    their own feature cache without editing the script."""
    text = _driver_text()
    for var in ("DEEPCATCH_FEATURES_DIR", "DEEPCATCH_LABELS",
                "DEEPCATCH_EXTERNAL_DIR", "DEEPCATCH_METHYLATION_DIR"):
        assert var in text, f"{var} is not honoured by paper/REPRODUCE.sh"


# ── the driver must never write to a published path by accident ────────
#
# Several scripts default their output to a path inside the repo:
# `cross_study_finallydb.py` defaults `--out-json` to
# `results/cross_study_finallydb.json`, `--out-md` to
# `docs/CROSS_STUDY_BENCHMARK.md` and `--out-per-cancer-json` to
# `results/per_cancer_sens_at_spec.json`. A driver that invokes such a
# script without naming every one of those flags lets the script's default
# decide where the published artifact goes.
#
# That is not hypothetical. The driver ran the LEGACY (no-GC) sweep without
# `--out-per-cancer-json`, so the legacy run overwrote the committed
# GC-corrected `results/per_cancer_sens_at_spec.json` with legacy numbers
# (BRCA 0.9738 -> 0.9794) — a silent drift of an artifact the paper
# quotes. The check below makes that class of accident impossible to
# reintroduce silently.


def _invocation_flags(text: str, script_name: str):
    """Every `--flag` passed to `script_name` anywhere in the driver."""
    flags = set()
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if script_name not in line:
            continue
        block = [line]
        j = i
        while j + 1 < len(lines) and lines[j].rstrip().endswith("\\"):
            j += 1
            block.append(lines[j])
        for raw in block:
            flags.update(re.findall(r"(--[a-z0-9\-]+)", raw))
    return flags


def _published_path_defaults(script_path: Path):
    """Flags of `script_path` whose argparse default is inside the repo."""
    source = script_path.read_text()
    defaults = {}
    for match in re.finditer(
            r'add_argument\(\s*"(--[a-z0-9\-]+)"[^)]*?default="((?:results|docs)/[^"]*)"',
            source, re.S):
        defaults[match.group(1)] = match.group(2)
    return defaults


@pytest.mark.parametrize("script_rel", [
    "scripts/cross_study_finallydb.py",
    "scripts/cross_study_loso.py",
    "scripts/external_validation_finaledb.py",
    "scripts/cross_study_finallydb_shuffled_control.py",
    "scripts/methylation_cohort_channel.py",
])
def test_driver_names_every_published_output_path(script_rel):
    script_path = _REPO_ROOT / script_rel
    if not script_path.is_file():
        pytest.skip(f"{script_rel} not present")
    defaults = _published_path_defaults(script_path)
    if not defaults:
        pytest.skip(f"{script_rel} has no published-path output defaults")
    passed = _invocation_flags(_driver_text(), script_path.name)
    missing = sorted(set(defaults) - passed)
    assert not missing, (
        f"paper/REPRODUCE.sh invokes {script_rel} without naming "
        f"{missing}. Those flags default to the PUBLISHED paths "
        f"{ {f: defaults[f] for f in missing} }, so the script would "
        f"overwrite a committed artifact instead of writing where the "
        f"driver intends."
    )


# ── the check itself ────────────────────────────────────────────────────


def _run_check(baseline: Path, current: Path, tolerance: str = "1e-6"):
    env = {**os.environ}
    env.pop("PYTHONPATH", None)
    return subprocess.run(
        [sys.executable, str(CHECK),
         "--baseline", str(baseline), "--current", str(current),
         "--tolerance", tolerance],
        capture_output=True, text=True, env=env, timeout=120)


def _seed_two_dirs(tmp_path: Path):
    """Copy the committed artifacts into two temp dirs."""
    baseline = tmp_path / "baseline"
    current = tmp_path / "current"
    baseline.mkdir()
    current.mkdir()
    src = _REPO_ROOT / "results"
    for path in src.glob("*.json"):
        shutil.copy(path, baseline / path.name)
        shutil.copy(path, current / path.name)
    return baseline, current


def test_reproduce_check_accepts_identical_artifacts(tmp_path):
    baseline, current = _seed_two_dirs(tmp_path)
    proc = _run_check(baseline, current)
    assert proc.returncode == 0, (
        f"identical artifacts were reported as drifted: "
        f"{proc.stdout[-2000:]!r} {proc.stderr!r}"
    )
    assert "reproduced" in proc.stdout


def test_reproduce_check_detects_a_perturbed_number(tmp_path):
    """The red case. If this passes, the check is useless."""
    baseline, current = _seed_two_dirs(tmp_path)
    target = current / "cross_study_finallydb_gc_corrected.json"
    blob = json.loads(target.read_text())
    blob["pooled"]["harmonized"]["auc_mean"] += 0.01
    target.write_text(json.dumps(blob, indent=2))

    proc = _run_check(baseline, current)
    assert proc.returncode == 1, (
        f"a perturbed artifact was NOT reported as drift: "
        f"{proc.stdout[-2000:]!r}"
    )
    assert "DRIFT" in proc.stdout
    assert "cross_study_finallydb_gc_corrected.json" in proc.stdout


def test_reproduce_check_detects_a_missing_artifact(tmp_path):
    baseline, current = _seed_two_dirs(tmp_path)
    (current / "cross_study_loso.json").unlink()
    proc = _run_check(baseline, current)
    assert proc.returncode == 1, (
        f"a missing artifact was NOT reported: {proc.stdout[-2000:]!r}"
    )


def test_reproduce_check_rejects_an_empty_baseline(tmp_path):
    baseline = tmp_path / "empty"
    baseline.mkdir()
    current, _ = _seed_two_dirs(tmp_path)
    proc = _run_check(baseline, current)
    assert proc.returncode == 1, (
        "an empty baseline must fail the check, not pass it silently"
    )


def _run_check_with_tolerances(baseline: Path, current: Path,
                                tolerance: str, threshold_tolerance: str):
    env = {**os.environ}
    env.pop("PYTHONPATH", None)
    return subprocess.run(
        [sys.executable, str(CHECK),
         "--baseline", str(baseline), "--current", str(current),
         "--tolerance", tolerance,
         "--threshold-tolerance", threshold_tolerance],
        capture_output=True, text=True, env=env, timeout=120)


def test_threshold_statistics_get_their_own_tolerance(tmp_path):
    """A sens@spec operating point is a step function of the score
    distribution, so it is far more seed-sensitive than an AUC. It must
    be judged against the threshold tolerance, not the smooth one."""
    baseline, current = _seed_two_dirs(tmp_path)
    target = current / "cross_study_finallydb_gc_corrected.json"
    blob = json.loads(target.read_text())
    blob["pooled"]["harmonized"]["sens_at_spec_95"] += 0.02
    target.write_text(json.dumps(blob, indent=2))

    # Tight smooth tolerance, wide threshold tolerance -> the sens@spec
    # move is accepted.
    proc = _run_check_with_tolerances(baseline, current, "1e-6", "5e-2")
    assert proc.returncode == 0, (
        "a sens@spec shift inside the threshold tolerance must pass: "
        f"{proc.stdout[-1500:]!r}"
    )

    # Same shift, same narrow band for both -> it must be reported.
    proc = _run_check_with_tolerances(baseline, current, "1e-6", "1e-6")
    assert proc.returncode == 1, (
        "a sens@spec shift outside the threshold tolerance must fail: "
        f"{proc.stdout[-1500:]!r}"
    )
    assert "DRIFT" in proc.stdout


def test_an_auc_shift_is_not_absorbed_by_the_threshold_tolerance(tmp_path):
    """The wide threshold band must not become a loophole for the AUCs."""
    baseline, current = _seed_two_dirs(tmp_path)
    target = current / "cross_study_finallydb_gc_corrected.json"
    blob = json.loads(target.read_text())
    blob["pooled"]["harmonized"]["auc_mean"] += 0.02
    target.write_text(json.dumps(blob, indent=2))

    proc = _run_check_with_tolerances(baseline, current, "1e-6", "5e-2")
    assert proc.returncode == 1, (
        "an AUC shift must fail even with a wide threshold tolerance: "
        f"{proc.stdout[-1500:]!r}"
    )
    assert "DRIFT" in proc.stdout
