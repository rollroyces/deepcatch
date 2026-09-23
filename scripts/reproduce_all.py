#!/usr/bin/env python3
"""Cross-platform equivalent of `scripts/reproduce_all.sh`.

Honors the same deepcatch-skill pitfalls as the bash script:

  #1  macOS global PYTHONPATH hijack — every subprocess inherits
      `env={'PYTHONPATH': ''}` so the venv's site-packages win, not
      Hermes's hermes-agent venv.
  #4  We DO NOT use shell `... | tee log` patterns. Each subprocess is
      run with `subprocess.run` and the actual return code is what
      decides. No pipefail needed because there's no pipe.
  #7  Resolves the repo root via walk-up fallback — works whether
      invoked via the editable-install symlink under
      `.venv/lib/python*/site-packages/` or directly from the repo.

The bash script remains the source of truth for what gets run. This
file is a 1:1 cross-platform mirror for users without bash (Windows,
Alpine musl, etc.). Both produce the same `results/REPRODUCE_REPORT.md`
and the same set of `results/<slug>.json` artifacts.

Usage
-----
    env -u PYTHONPATH ./.venv/bin/python scripts/reproduce_all.py
    env -u PYTHONPATH ./.venv/bin/python scripts/reproduce_all.py --dry-run
    env -u PYTHONPATH ./.venv/bin/python scripts/reproduce_all.py --only cross_study,per_cancer
    env -u PYTHONPATH ./.venv/bin/python scripts/reproduce_all.py --skip-tests

Exit code: 0 if every non-skipped step PASSed, non-zero otherwise.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ──────────────────────────────────────────────────────────────────────
# Walk-up repo-root resolver — see deepcatch skill pitfall #7
# ──────────────────────────────────────────────────────────────────────


def _resolve_repo_root() -> Path:
    """Walk up from this file until we find the repo root.

    Anchor on sentinel files (pyproject.toml + scripts/) so we work
    whether invoked via editable-install symlink under
    `.venv/lib/python*/site-packages/scripts/reproduce_all.py` or
    directly from the repo.
    """
    cur = Path(__file__).resolve().parent
    while cur != cur.parent:
        if (cur / "pyproject.toml").is_file() and (cur / "scripts").is_dir():
            return cur
        cur = cur.parent
    return Path.cwd()


REPO_ROOT = _resolve_repo_root()
PY = REPO_ROOT / ".venv" / "bin" / "python"
RESULTS = REPO_ROOT / "results"


# ──────────────────────────────────────────────────────────────────────
# Step registry — kept in sync with reproduce_all.sh
# ──────────────────────────────────────────────────────────────────────


@dataclass
class Step:
    slug: str
    description: str
    budget_seconds: int
    # What the bash script does. We mirror that here so this file is
    # the canonical cross-platform reference.
    command: List[str] = field(default_factory=list)


def _default_steps(finallydb_dir: str, tcga_cache_dir: str, py: str) -> List[Step]:
    return [
        Step(
            slug="cross_study",
            description="cross_study_finallydb.py (FinaleDB Jiang+Cristiano)",
            budget_seconds=300,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "cross_study_finallydb.py"),
                "--features-dir", finallydb_dir,
                "--out-json", str(RESULTS / "cross_study.json"),
                "--out-md", str(REPO_ROOT / "docs" / "CROSS_STUDY_BENCHMARK.md"),
            ],
        ),
        Step(
            slug="per_cancer",
            description="per_cancer_sens_at_spec.py --synthetic",
            budget_seconds=45,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "per_cancer_sens_at_spec.py"),
                "--synthetic", "--n", "240", "--seed", "42",
                "--out", str(RESULTS / "per_cancer.json"),
            ],
        ),
        Step(
            slug="harmonization",
            description="harmonization_check.py",
            budget_seconds=15,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "harmonization_check.py"),
            ],
        ),
        Step(
            slug="focal_bce_ablation",
            description="sens_at_spec_ablation.py (CE vs focal-BCE)",
            budget_seconds=30,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "sens_at_spec_ablation.py"),
                "--ce-json", str(RESULTS / "foundation_smoke_ce.json"),
                "--sens-json", str(RESULTS / "foundation_smoke_sens.json"),
                "--out", str(RESULTS / "focal_bce_ablation.json"),
            ],
        ),
        Step(
            slug="sparse_aware_ablation",
            description="sparse_aware_ablation.py (Linear vs SparseAware)",
            budget_seconds=30,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "sparse_aware_ablation.py"),
                "--linear-json", str(RESULTS / "foundation_smoke.json"),
                "--sparse-json", str(RESULTS / "foundation_smoke_sparse.json"),
                "--out", str(RESULTS / "sparse_aware_ablation.json"),
            ],
        ),
        Step(
            slug="foundation_smoke",
            description="foundation_real_smoke.py --quick (3x: CE + focal-BCE + sparse_aware)",
            budget_seconds=240,
            command=[
                # Three sequential smoke runs so the ablation steps have
                # both inputs. The bash script does this via a single
                # shell step that internally runs three commands; we
                # expand it into three Step entries with separate slugs.
                py,
                str(REPO_ROOT / "scripts" / "foundation_real_smoke.py"),
                "--quick", "--n-patients", "20", "--seeds", "5",
                "--features-dir", tcga_cache_dir,
                "--out", str(RESULTS / "foundation_smoke_ce.json"),
            ],
        ),
        Step(
            slug="pytest_tests",
            description="pytest test/ + src/foundation/test_integration.py",
            budget_seconds=180,
            command=[
                py, "-m", "pytest",
                str(REPO_ROOT / "test"),
                str(REPO_ROOT / "src" / "foundation" / "test_integration.py"),
                "--tb=line", "--timeout=180", "--no-header", "-q",
                "--deselect",
                "test/test_cross_study_per_cancer_integration.py::test_end_to_end_script_writes_both_outputs",
            ],
        ),
    ]


# Additional "internal" smoke runs that the bash script bundles inside
# the `foundation_smoke` step. We model them as separate Steps so they
# each get their own row in the report and the ablations can find their
# expected JSONs. Their `internal_parent_slug` ties them to the
# `foundation_smoke` row for the report.

def _foundation_smoke_subs(
    py: str, tcga_cache_dir: str
) -> List[Step]:
    return [
        Step(
            slug="foundation_smoke_sens",
            description="foundation_real_smoke.py --loss sens_at_spec (input for focal-BCE ablation)",
            budget_seconds=80,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "foundation_real_smoke.py"),
                "--quick", "--n-patients", "20", "--seeds", "5",
                "--features-dir", tcga_cache_dir,
                "--loss", "sens_at_spec", "--alpha-pos", "20",
                "--out", str(RESULTS / "foundation_smoke_sens.json"),
            ],
        ),
        Step(
            slug="foundation_smoke_sparse",
            description="foundation_real_smoke.py --projection-kinds sparse_aware (input for sparse-aware ablation)",
            budget_seconds=80,
            command=[
                py,
                str(REPO_ROOT / "scripts" / "foundation_real_smoke.py"),
                "--quick", "--n-patients", "20", "--seeds", "5",
                "--features-dir", tcga_cache_dir,
                "--projection-kinds", '{"frag_basic":"sparse_aware"}',
                "--out", str(RESULTS / "foundation_smoke_sparse.json"),
            ],
        ),
    ]


# ──────────────────────────────────────────────────────────────────────
# Runner
# ──────────────────────────────────────────────────────────────────────


@dataclass
class StepResult:
    slug: str
    description: str
    rc: str  # "PASS" / "FAIL" / "SKIP" / "DRY"
    wall: float
    note: str = ""
    output_path: Optional[Path] = None


def _run_subprocess(
    cmd: List[str], cwd: Path, log_path: Path
) -> Tuple[int, str]:
    """Run `cmd` to completion, streaming to `log_path`, returning (rc, tail).

    Honors pitfall #1: PYTHONPATH explicitly blanked so the global macOS
    PYTHONPATH can't hijack the venv.

    Honors pitfall #4 implicitly: we capture stdout/stderr ourselves
    rather than piping through `tee`, so the actual subprocess rc is
    what we read. No `pipefail` shell flag needed.
    """
    env = {**os.environ, "PYTHONPATH": ""}
    with log_path.open("w") as logf:
        proc = subprocess.run(
            cmd, cwd=str(cwd), env=env,
            stdout=logf, stderr=subprocess.STDOUT,
            check=False,
        )
    rc = int(proc.returncode) if proc.returncode is not None else -1
    tail = ""
    try:
        lines = log_path.read_text(errors="replace").splitlines()[-3:]
        tail = " ".join(lines).strip()[:200]
    except Exception:
        pass
    return rc, tail


def _synthesize_pytest_json(log_path: Path, json_out: Path) -> int:
    """Parse pytest output for `N passed` and write a summary JSON.

    Returns 0 on success.
    """
    try:
        text = log_path.read_text(errors="replace")
        # pytest -q prints e.g. "138 passed in 31.20s"
        import re
        m = re.search(r"(\d+) passed", text)
        n_pass = m.group(1) if m else "0"
        payload = {
            "step": "pytest_tests",
            "tests_passed": n_pass,
            "source": "results/pytest_tests.log",
            "deselected": [
                "test_cross_study_per_cancer_integration.py::test_end_to_end_script_writes_both_outputs (5-10 min subprocess)"
            ],
        }
        json_out.write_text(json.dumps(payload, indent=2))
        return 0
    except Exception as e:
        return 1


def run_step(
    step: Step,
    only_filter: Optional[List[str]],
    dry_run: bool,
    log_dir: Path,
) -> StepResult:
    """Run one step and return a StepResult."""
    if only_filter and step.slug not in only_filter:
        return StepResult(
            slug=step.slug,
            description=step.description,
            rc="SKIP",
            wall=0.0,
            note="skipped via --only",
        )

    if dry_run:
        print(f"==[{step.slug}]== {step.description}")
        print(f"    [dry-run] would run: {' '.join(step.command)}")
        return StepResult(
            slug=step.slug,
            description=step.description,
            rc="DRY",
            wall=0.0,
            note="dry-run",
        )

    log_path = log_dir / f"{step.slug}.log"
    json_out = RESULTS / f"{step.slug}.json"
    print(f"==[{step.slug}]== {step.description}")
    print(f"    out: {json_out}")
    print(f"    log: {log_path}")

    t0 = time.monotonic()
    rc, tail_msg = _run_subprocess(step.command, REPO_ROOT, log_path)
    wall = time.monotonic() - t0

    # harmonization_check.py writes to results/harmonization_check.json
    # hardcoded; move it to our slug-named JSON. Same as bash script.
    if step.slug == "harmonization":
        hardcoded = RESULTS / "harmonization_check.json"
        if hardcoded.exists() and hardcoded != json_out:
            try:
                shutil.move(str(hardcoded), str(json_out))
            except Exception:
                pass

    # foundation_smoke's "main" JSON is the CE run. The bash script
    # copies CE → foundation_smoke.json so other tools can find a
    # single canonical smoke output.
    if step.slug == "foundation_smoke":
        if (RESULTS / "foundation_smoke_ce.json").exists():
            shutil.copy(RESULTS / "foundation_smoke_ce.json", json_out)

    # pytest gets a synthesized JSON with the pass count.
    if step.slug == "pytest_tests" and rc == 0:
        _synthesize_pytest_json(log_path, json_out)

    ok = (rc == 0) and json_out.exists()
    return StepResult(
        slug=step.slug,
        description=step.description,
        rc="PASS" if ok else "FAIL",
        wall=round(wall, 1),
        note="" if ok else f"exit={rc}, {tail_msg}",
        output_path=json_out,
    )


def _print_table(results: List[StepResult], total_wall: float) -> None:
    print()
    print(f"{'=' * 78}")
    print(f"  Step  {'Status':<10} {'Wall (s)':>8}  Description")
    print(f"{'=' * 78}")
    for r in results:
        marker = {
            "PASS": "✅", "FAIL": "❌", "SKIP": "⏭️", "DRY": "🟡"
        }.get(r.rc, "?")
        print(f"  {marker}  {r.slug:<22} {r.rc:<10} {r.wall:>7.1f}  {r.description}")
    print(f"{'=' * 78}")
    print(f"  Total wall: {total_wall:.1f}s")
    print(f"{'=' * 78}")


def _write_report(
    results: List[StepResult], total_wall: float, started: float, repo_root: Path
) -> Path:
    """Mirror reproduce_all.sh's report.md exactly so both produce the
    same artifact."""
    report = RESULTS / "REPRODUCE_REPORT.md"
    lines: List[str] = []
    lines.append("# Reproduction report\n")
    lines.append(f"- Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(started))}")
    lines.append(f"- repo_root: `{repo_root}`")
    py_ver = subprocess.run(
        [str(PY), "-V"], capture_output=True, text=True, env={"PYTHONPATH": ""},
    ).stdout.strip() if PY.exists() else "(missing)"
    lines.append(f"- python:    `{PY}` ({py_ver})")
    lines.append(f"- total wall: {total_wall:.1f}s")
    lines.append("")
    lines.append("## Step results\n")
    lines.append("| Step | Description | Wall (s) | Exit | Output |")
    lines.append("|---|---|---:|---|---|")
    for r in results:
        marker = {
            "PASS": "✅", "FAIL": "❌", "SKIP": "⏭️", "DRY": "🟡"
        }.get(r.rc, "?")
        out_name = f"{r.slug}.json"
        if r.output_path and not r.output_path.exists():
            out_name = f"{r.slug}.json (MISSING)"
        lines.append(
            f"| {marker} `{r.slug}` | {r.description} | {r.wall} | {r.rc} | "
            f"`{out_name}` |"
        )
    lines.append("")
    lines.append("## Per-step logs\n")
    lines.append("Each step writes `results/<slug>.log` with the full stdout/stderr of the script.")
    lines.append("")
    lines.append("## Honest framing\n")
    lines.append("This script reproduces the **open-data methods benchmark** numbers published in")
    lines.append("`docs/CROSS_STUDY_BENCHMARK.md`, `docs/SENS_AT_SPEC_ABLATION.md`,")
    lines.append("`docs/SPARSE_AWARE_ABLATION.md`, `docs/HARMONIZATION.md`, and the foundation")
    lines.append("smoke's `docs/FOUNDATION_REAL_SMOKE.md`-shaped numbers.")
    lines.append("")
    lines.append("It does **not** perform:")
    lines.append("- Clinical validation on an external held-out cohort")
    lines.append("- Regulatory-grade quality controls")
    lines.append("- Per-sample audit/reproducibility beyond what the per-seed JSONs capture")
    lines.append("")
    lines.append("## Re-running individual steps\n")
    lines.append("```bash")
    lines.append("# Cross-study benchmark alone")
    lines.append("env -u PYTHONPATH ./.venv/bin/python scripts/cross_study_finallydb.py \\")
    lines.append("    --features-dir $FINALEDB_FEATURES_DIR \\")
    lines.append("    --out-json results/cross_study.json \\")
    lines.append("    --out-md docs/CROSS_STUDY_BENCHMARK.md")
    lines.append("")
    lines.append("# Per-cancer sens@spec alone (synthetic fixture)")
    lines.append("env -u PYTHONPATH ./.venv/bin/python scripts/per_cancer_sens_at_spec.py \\")
    lines.append("    --synthetic --n 240 --out results/per_cancer.json")
    lines.append("```")
    lines.append("")
    report.write_text("\n".join(lines))
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--only",
        help="comma-separated subset of slugs to run",
    )
    ap.add_argument(
        "--skip-tests",
        action="store_true",
        help="skip the pytest_tests step",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print plan, do not run anything",
    )
    ap.add_argument(
        "--finallydb-dir",
        default=os.environ.get(
            "FINALEDB_FEATURES_DIR",
            "/Users/hermes/cfdna-fragmentomics-pipeline/data/features",
        ),
        help="FinaleDB feature cache dir (default: cfdna-fragmentomics-pipeline)",
    )
    ap.add_argument(
        "--tcga-cache-dir",
        default=os.environ.get(
            "TCGA_CACHE_DIR",
            str(REPO_ROOT / "validation" / "tcga" / "tcga_cache"),
        ),
        help="TCGA cache dir for foundation_real_smoke.py",
    )
    args = ap.parse_args()

    # Pre-flight
    if not PY.exists():
        print(f"❌ missing venv python at {PY}", file=sys.stderr)
        print(
            "   Create with: python3 -m venv .venv && pip install -e .",
            file=sys.stderr,
        )
        return 1

    if not Path(args.finallydb_dir).is_dir():
        print(
            f"⚠️  FinaleDB features dir not found: {args.finallydb_dir}",
            file=sys.stderr,
        )

    RESULTS.mkdir(parents=True, exist_ok=True)

    only_filter: Optional[List[str]] = None
    if args.only:
        only_filter = [s.strip() for s in args.only.split(",") if s.strip()]
    if args.skip_tests and only_filter is None:
        only_filter = None  # handled per-step

    # Pre-flight env check
    print("==[verify_environment]==")
    if args.dry_run:
        print(f"    [dry-run] would run: env -u PYTHONPATH {PY} scripts/verify_environment.py")
    else:
        subprocess.run(
            [str(PY), str(REPO_ROOT / "scripts" / "verify_environment.py")],
            env={**os.environ, "PYTHONPATH": ""},
            check=False,
        )
    print()

    # Header
    print("=" * 78)
    print("  deepcatch reproduce_all.py (cross-platform equivalent of reproduce_all.sh)")
    print(f"  repo_root : {REPO_ROOT}")
    print(f"  python    : {PY}")
    print(f"  results   : {RESULTS}")
    print(f"  features  : {args.finallydb_dir}")
    print(f"  tcga      : {args.tcga_cache_dir}")
    print(f"  platform  : {platform.system()} {platform.machine()}")
    print("=" * 78)
    print()

    steps: List[Step] = []
    for s in _default_steps(args.finallydb_dir, args.tcga_cache_dir, str(PY)):
        steps.append(s)
        if s.slug == "foundation_smoke":
            for sub in _foundation_smoke_subs(str(PY), args.tcga_cache_dir):
                steps.append(sub)

    results: List[StepResult] = []
    t_total_0 = time.monotonic()
    for step in steps:
        if args.skip_tests and step.slug == "pytest_tests":
            results.append(StepResult(
                slug=step.slug,
                description=step.description,
                rc="SKIP",
                wall=0.0,
                note="skipped via --skip-tests",
            ))
            continue
        r = run_step(step, only_filter, args.dry_run, log_dir=RESULTS)
        results.append(r)
        print(f"    wall: {r.wall}s  rc: {r.rc}\n")

    t_total_1 = time.monotonic()
    total_wall = t_total_1 - t_total_0

    _print_table(results, total_wall)

    if not args.dry_run:
        report = _write_report(results, total_wall, t_total_0, REPO_ROOT)
        print()
        print(f"  Reproduction report written: {report}")
        print(f"  Total wall: {total_wall:.1f}s")

    n_fail = sum(1 for r in results if r.rc == "FAIL")
    if n_fail > 0:
        print(f"\n❌ {n_fail} step(s) failed. See REPRODUCE_REPORT.md for details.")
        return 1
    print("\n✅ All steps passed (or were skipped / dry-run).")
    return 0


if __name__ == "__main__":
    sys.exit(main())