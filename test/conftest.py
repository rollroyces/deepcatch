"""Shared pytest configuration for the DeepCatch test suite.

Some tests drive the cross-study benchmark scripts, which call into the
sibling ``cfdna-fragmentomics-pipeline`` checkout (its 5-channel
loader and its per-study harmonizer). That checkout is not present on a
CI runner or in a fresh clone, so those tests are marked and skipped
there rather than failing on ImportError.

The same applies to the tests that drive ``scripts/reproduce_all.sh``,
which invokes the repo's ``.venv`` interpreter by path.

Tests marked this way still run -- and must still pass -- wherever the
sibling checkout / venv is available, so this narrows nothing locally.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path

# The repository's own path layout, resolved from this file: the
# sibling pipeline, the feature cache and the labels are looked up
# next to the repository, then in the home directory, with the
# environment variables overriding both.
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "scripts"))
from _paths import SIBLING_PIPELINE  # noqa: E402


import os
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

SISTER_PIPELINE = Path(os.environ.get(
    "DEEPCATCH_SISTER_PIPELINE",
    str(SIBLING_PIPELINE),
))
SISTER_PIPELINE_PRESENT = (
    SISTER_PIPELINE / "scripts" / "honest_benchmark.py"
).exists()
VENV_PRESENT = (_REPO_ROOT / ".venv" / "bin" / "python").exists()

_SISTER_REASON = (
    "requires the sibling cfdna-fragmentomics-pipeline checkout "
    "(DEEPCATCH_SISTER_PIPELINE)"
)
_VENV_REASON = (
    "requires the repo .venv (scripts/reproduce_all.sh is venv-only)"
)


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_sister_pipeline: needs the sibling pipeline checkout",
    )
    config.addinivalue_line(
        "markers",
        "requires_venv: needs the repo .venv interpreter",
    )


def pytest_collection_modifyitems(config, items):
    skip_sister = pytest.mark.skip(reason=_SISTER_REASON)
    skip_venv = pytest.mark.skip(reason=_VENV_REASON)
    for item in items:
        if ("requires_sister_pipeline" in item.keywords
                and not SISTER_PIPELINE_PRESENT):
            item.add_marker(skip_sister)
        if "requires_venv" in item.keywords and not VENV_PRESENT:
            item.add_marker(skip_venv)
