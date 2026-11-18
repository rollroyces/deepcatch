"""The quickstart notebook must still run.

Why this file exists
--------------------
The README's second demo is a recording of this notebook, and the
notebook is the one place in the repository where a reader sees the
package used the way a bioinformatician would use it.

Its histogram cell wrote to a relative path -- `figures/fsi_hist.png` --
and no directory of that name exists in the repository, so the cell
failed with a FileNotFoundError wherever the notebook was run from. The
committed outputs came from a machine where the directory happened to
exist, so nothing in CI could see it: the notebook is not a test file,
and no test ran it.

What this test checks, and what it deliberately does not
-------------------------------------------------------
Every computing cell is executed, in order, with the working directory
the notebook assumes (its own directory -- cell 1 puts `../src` on the
path from there).

A cell that renders an image through IPython is skipped when IPython is
not importable: the core CI job installs numpy, matplotlib and pytest,
not IPython, and the rendering is a display convenience rather than a
computation. The cell that failed was not one of those.
The notebook is executed in a subprocess so that its `sys.path`
manipulation and its matplotlib state cannot leak into the test run.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_NOTEBOOK = _REPO_ROOT / "notebooks" / "deepcatch_quickstart.ipynb"

_RUNNER = textwrap.dedent(
    r"""
    import json, re, sys, traceback

    nb = json.load(open(sys.argv[1]))
    namespace = {"__name__": "__main__"}
    ran = 0
    skipped = 0
    for index, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        # The display cells come in pairs: one imports IPython's Image and
        # display, the next uses them. Both are rendering, not computing.
        if "IPython" in source or re.search(r"\b(display|Image)\s*\(", source):
            skipped += 1
            continue
        try:
            exec(compile(source, "<cell %d>" % index, "exec"), namespace)
        except Exception:
            traceback.print_exc()
            sys.exit("cell %d failed" % index)
        ran += 1
    print("executed %d code cells, skipped %d rendering cells" % (ran, skipped))
    """
)


def test_the_quickstart_notebook_runs():
    """The notebook a reader is told to open must execute."""
    assert _NOTEBOOK.is_file(), "the quickstart notebook is missing"

    env = dict(os.environ, MPLBACKEND="Agg")
    try:
        proc = subprocess.run(
            [sys.executable, "-c", _RUNNER, str(_NOTEBOOK)],
            cwd=str(_NOTEBOOK.parent),
            capture_output=True, text=True, timeout=600, env=env,
        )
    except subprocess.TimeoutExpired:
        pytest.fail("the quickstart notebook did not finish in ten minutes")

    assert proc.returncode == 0, (
        "the quickstart notebook does not run:\n"
        f"{proc.stdout}\n{proc.stderr}\n"
        "The README shows this notebook working, and a reader who copies it "
        "gets the same failure."
    )
    assert "executed" in proc.stdout, (
        f"the notebook run reported nothing:\n{proc.stdout}\n{proc.stderr}")
