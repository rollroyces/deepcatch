"""No module may carry a LaTeX docstring that Python mangles.

Why this file exists
--------------------
Five modules documented their maths with LaTeX in ordinary (non-raw)
docstrings. Python then interpreted the escapes: `\\frac` became a form
feed followed by "rac", `\\text` kept its backslash but `\\mid`, `\\sum`
and `\\epsilon` produced invalid-escape warnings on every import. The
rendered text of those docstrings was wrong, and Python has said for
years that such sequences "will not work in the future".

The fix is to make those docstrings raw strings, which restores the
intended LaTeX. This test keeps it fixed: every tracked module is
compiled with warnings recorded, and an invalid escape sequence is a
failure.

What this test checks, and what it deliberately does not
-------------------------------------------------------
It compiles each file and looks at the warnings, which is what the
interpreter itself objects to. It does not import the modules: importing
would need the package installed and the optional extras present, and the
objection here is lexical.

Files that cannot be compiled at all are skipped rather than failed --
a syntax error is the test suite's business, and this test is about
escapes.
"""
from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _tracked_python_files() -> list[Path]:
    try:
        proc = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), "ls-files", "*.py"],
            capture_output=True, text=True, timeout=60, check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout; the tracked tree cannot be listed")
    return [_REPO_ROOT / name for name in proc.stdout.splitlines() if name]


def test_no_module_carries_an_invalid_escape_sequence():
    """LaTeX in a docstring must be a raw string."""
    offenders: list[str] = []

    for path in _tracked_python_files():
        if not path.is_file():
            continue
        try:
            source = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                compile(source, str(path), "exec")
            except SyntaxError:
                continue
            for entry in caught:
                if "escape sequence" in str(entry.message):
                    offenders.append(
                        f"{path.relative_to(_REPO_ROOT)}:{entry.lineno}  "
                        f"{entry.message}")

    assert not offenders, (
        "a module's docstring carries an invalid escape sequence:\n  "
        + "\n  ".join(offenders[:40])
        + "\nMake the docstring a raw string (r\"\"\"): Python would "
          "otherwise interpret the LaTeX and print mangled text."
    )
