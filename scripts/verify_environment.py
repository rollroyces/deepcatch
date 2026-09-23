#!/usr/bin/env python3
"""Pre-flight dependency / cache check for the deepcatch reproduce pipeline.

Honors the deepcatch skill pitfalls:

  #1  macOS global PYTHONPATH hijack — we resolve `_REPO_ROOT` from a
      walk-up fallback (NOT `__file__`-relative through an editable
      install symlink) and print it so the user can verify before
      running the heavier scripts.
  #7  editable-install `__file__` paths — same walk-up fallback so this
      script works whether invoked from `pytest`, a bare
      `python scripts/verify_environment.py`, or a `cd /tmp && ...`.

Usage
-----
    env -u PYTHONPATH ./.venv/bin/python scripts/verify_environment.py

Exit code: 0 if every required dependency is present AND the
      FinaleDB feature cache directory exists; non-zero otherwise.

The script is deliberately read-only and side-effect-free — it does
NOT create directories, does NOT download data, does NOT touch
results/. It only inspects.
"""

from __future__ import annotations

import importlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ──────────────────────────────────────────────────────────────────────
# Walk-up repo-root resolver — see deepcatch skill pitfall #7
# ──────────────────────────────────────────────────────────────────────


def _resolve_repo_root() -> Path:
    """Walk up from this file until we find the deepcatch repo root.

    This avoids the pitfall where, under an editable install,
    `Path(__file__).resolve()` points at
    `.venv/lib/python3.14/site-packages/scripts/verify_environment.py`
    instead of the real `scripts/verify_environment.py`. We anchor on
    sentinel files that only the repo root contains (`pyproject.toml`
    + `scripts/` together).
    """
    cur = Path(__file__).resolve().parent
    while cur != cur.parent:
        if (cur / "pyproject.toml").is_file() and (cur / "scripts").is_dir():
            return cur
        cur = cur.parent
    # Last-ditch: fall back to CWD
    return Path.cwd()


_REPO_ROOT = _resolve_repo_root()
_VENV_PY = _REPO_ROOT / ".venv" / "bin" / "python"
_DEFAULT_FEATURES_DIR = Path(
    "/Users/hermes/cfdna-fragmentomics-pipeline/data/features"
)


# ──────────────────────────────────────────────────────────────────────
# Dependency table — required versions track what was used when the
# published numbers were generated. Anything below the required version
# is treated as a FAIL (we won't try to run the heavier reproduce
# steps on a known-broken stack).
# ──────────────────────────────────────────────────────────────────────

REQUIRED_PYTHON: Tuple[int, int] = (3, 11)

REQUIRED_DEPS: List[Dict[str, object]] = [
    {
        "name": "numpy",
        "import": "numpy",
        "min_version": "1.24.0",
        "purpose": "core numerics (cross-study benchmark, ablations)",
    },
    {
        "name": "scipy",
        "import": "scipy",
        "min_version": "1.10.0",
        "purpose": "DeLong CIs, paired t-tests",
    },
    {
        "name": "scikit-learn",
        "import": "sklearn",
        "min_version": "1.3.0",
        "purpose": "LogisticRegression / GroupKFold / StandardScaler",
    },
    {
        "name": "pandas",
        "import": "pandas",
        "min_version": "2.0.0",
        "purpose": "labels_multiclass.tsv ingestion",
    },
    {
        "name": "torch",
        "import": "torch",
        "min_version": "2.13.0",
        "purpose": "Foundation model (foundation_real_smoke)",
    },
    {
        "name": "pytest",
        "import": "pytest",
        "min_version": "7.0.0",
        "purpose": "test/ + src/foundation/test_integration.py",
    },
]


def _parse_version(s: str) -> Tuple[int, ...]:
    """Parse '2.3.5' / '1.10' / '2.13.0+cpu' into a tuple of ints."""
    out: List[int] = []
    buf = ""
    for ch in s:
        if ch.isdigit():
            buf += ch
        elif buf:
            out.append(int(buf))
            buf = ""
        if len(out) == 3:  # major, minor, micro — that's enough
            break
    if buf:
        out.append(int(buf))
    return tuple(out)


def _version_ge(found: str, required: str) -> bool:
    """Return True iff found >= required (best-effort numeric compare)."""
    try:
        return _parse_version(found) >= _parse_version(required)
    except Exception:
        # If we can't parse, fail open — let the actual reproduce step
        # complain instead of double-validating here.
        return True


# ──────────────────────────────────────────────────────────────────────
# Checks
# ──────────────────────────────────────────────────────────────────────


def check_python() -> Dict[str, object]:
    vi = sys.version_info
    found = f"{vi.major}.{vi.minor}.{vi.micro}"
    ok = (vi.major, vi.minor) >= REQUIRED_PYTHON
    return {
        "name": f"python>={'.'.join(map(str, REQUIRED_PYTHON))}",
        "found": found,
        "required": ".".join(map(str, REQUIRED_PYTHON)),
        "ok": ok,
        "note": "" if ok else "Python 3.11+ required (we ship 3.14.4)",
    }


def check_venv_python() -> Dict[str, object]:
    ok = _VENV_PY.is_file()
    found = str(_VENV_PY) if ok else "(missing)"
    return {
        "name": "deepcatch/.venv",
        "found": found,
        "required": str(_VENV_PY),
        "ok": ok,
        "note": "" if ok else "Create with `python3 -m venv .venv` then "
        "`pip install -e .` per pyproject.toml",
    }


def check_dep(dep: Dict[str, object]) -> Dict[str, object]:
    name = str(dep["name"])
    try:
        mod = importlib.import_module(str(dep["import"]))
    except Exception as e:  # pragma: no cover — defensive
        return {
            "name": name,
            "found": f"MISSING ({type(e).__name__})",
            "required": f">={dep['min_version']}",
            "ok": False,
            "note": str(dep.get("purpose", "")),
        }
    found_version = getattr(mod, "__version__", "?")
    ok = _version_ge(found_version, str(dep["min_version"]))
    return {
        "name": name,
        "found": found_version,
        "required": f">={dep['min_version']}",
        "ok": ok,
        "note": str(dep.get("purpose", "")),
    }


def check_features_dir() -> Dict[str, object]:
    """FinaleDB feature cache — required by cross_study_finallydb.py.

    Default path is the standard cfdna-fragmentomics-pipeline feature
    dump. If missing, the cross-study benchmark will exit non-zero.
    Allow override via $FINALEDB_FEATURES_DIR for non-Mac users.
    """
    override = os.environ.get("FINALEDB_FEATURES_DIR")
    path = Path(override) if override else _DEFAULT_FEATURES_DIR
    ok = path.is_dir()
    note = ""
    if not ok:
        note = (
            f"Cross-study benchmark (step 1) will fail without this. "
            f"Override via FINALEDB_FEATURES_DIR env var if you keep "
            f"FinaleDB elsewhere."
        )
    n_files = sum(1 for _ in path.glob("**/*")) if ok else 0
    return {
        "name": "FinaleDB features dir",
        "found": str(path) + (f" ({n_files} files)" if ok else ""),
        "required": str(_DEFAULT_FEATURES_DIR),
        "ok": ok,
        "note": note,
    }


def check_repo_layout() -> Dict[str, object]:
    """Confirm the script dir + test dir + results dir are present."""
    paths = [
        _REPO_ROOT / "scripts",
        _REPO_ROOT / "test",
        _REPO_ROOT / "src",
        _REPO_ROOT / "docs",
        _REPO_ROOT / "pyproject.toml",
    ]
    missing = [str(p) for p in paths if not p.exists()]
    ok = not missing
    return {
        "name": "repo layout",
        "found": "OK" if ok else f"missing: {missing}",
        "required": "scripts/ + test/ + src/ + docs/ + pyproject.toml",
        "ok": ok,
        "note": "" if ok else "Run from the deepcatch repo root",
    }


def check_results_dir_writable() -> Dict[str, object]:
    """results/ should exist + be writable (we don't write anything here,
    but the reproduce scripts will)."""
    p = _REPO_ROOT / "results"
    exists = p.exists()
    writable = exists and os.access(p, os.W_OK)
    ok = exists and writable
    return {
        "name": "results/ writable",
        "found": "OK" if ok else ("exists, not writable" if exists else "missing"),
        "required": "results/ exists + writable",
        "ok": ok,
        "note": "" if ok else "mkdir -p results && chmod u+w results",
    }


# ──────────────────────────────────────────────────────────────────────
# Report
# ──────────────────────────────────────────────────────────────────────


def _row(label: str, status: str, found: str, req: str, note: str) -> str:
    icon = "✅" if status == "PASS" else "❌"
    return (
        f"| {icon} {label} | {found} | {req} | {note} |\n"
    )


def _sep(label: str) -> str:
    # Marker row — caller may inject between header and body if desired.
    return f"| -- {label} -- |\n"


def render_markdown_table(rows: List[Dict[str, object]]) -> str:
    header = [
        "| Status | Check | Found | Required | Note |",
        "|---|---|---|---|",
        "",
    ]
    body_lines = [
        _row(
            str(r["name"]),
            "PASS" if r["ok"] else "FAIL",
            str(r["found"]),
            str(r["required"]),
            str(r["note"]),
        )
        for r in rows
    ]
    return "\n".join(header) + "\n" + "".join(body_lines)


def main() -> int:
    rows: List[Dict[str, object]] = []
    rows.append(check_python())
    rows.append(check_repo_layout())
    rows.append(check_venv_python())
    rows.append(check_results_dir_writable())
    rows.append(check_features_dir())
    for dep in REQUIRED_DEPS:
        rows.append(check_dep(dep))

    print(f"# deepcatch environment check\n")
    print(f"repo_root : {_REPO_ROOT}")
    print(f"python    : {sys.executable}")
    print(f"cwd       : {Path.cwd()}")
    print()
    print(render_markdown_table(rows))

    n_pass = sum(1 for r in rows if r["ok"])
    n_fail = sum(1 for r in rows if not r["ok"])
    print(f"\nSummary: {n_pass} pass / {n_fail} fail")

    # JSON for machine-readable usage (`./verify_environment.py --json`)
    if "--json" in sys.argv:
        sys.stdout.write(
            json.dumps(
                {
                    "repo_root": str(_REPO_ROOT),
                    "python": sys.version.split()[0],
                    "rows": rows,
                    "n_pass": n_pass,
                    "n_fail": n_fail,
                },
                indent=2,
            )
            + "\n"
        )

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())