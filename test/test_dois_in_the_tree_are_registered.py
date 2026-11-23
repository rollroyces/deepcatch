"""Every DOI in the tracked tree must resolve to a record.

Why this file exists
--------------------
A DOI that does not resolve is a dead link in a document a reader follows.
This repository had three of them in its bibliographies: a DOI whose digits
were transposed (`...-21394-w` for `...-24994-w`), one naming the PNAS
paper where the manuscript described the Cancer Discovery paper, and one
belonging to a preprint that could not be located at all. Two of them
looked exactly like working DOIs.

Nothing in CI could see them: the guards that existed pinned figures and
test counts, and an identifier that resolves to nothing is neither.

What this test checks, and what it deliberately does not
-------------------------------------------------------
It runs `scripts/check_dois.py`, which extracts every DOI from the tracked
tree and asks the registry that issued it -- Crossref for the journal
literature, DataCite for the Zenodo records this repository cites -- whether
it is registered.

Two registries, because two kinds of identifier are in use: a Zenodo DOI is
not in Crossref, so checking it there would report a failure that does not
exist.

The check needs the network. When a registry cannot be reached the script
exits 2 and the test skips rather than fails: an offline machine is not a
broken DOI.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CHECKER = _REPO_ROOT / "scripts" / "check_dois.py"


def _load_checker():
    """The checker, loaded from its path."""
    spec = importlib.util.spec_from_file_location("check_dois", _CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("text", [
    "10.5281/zenodo.24453044",
    "**10.5281/zenodo.24453044**",
    "**10.5281/zenodo.24453044.**",
    "_10.5281/zenodo.24453044_",
    "The DOI is 10.5281/zenodo.24453044.",
    "[10.5281/zenodo.24453044](https://doi.org/10.5281/zenodo.24453044)",
    "\\url{https://doi.org/10.5281/zenodo.24453044}",
])
def test_markup_is_not_part_of_the_identifier(text):
    """A DOI written in bold, italics, a link or a macro is the same DOI.

    The scanner used to keep the markdown emphasis, so a DOI written in
    bold was reported as an identifier that resolves to nothing -- a
    false failure in the check that is supposed to catch dead links.
    The failure it caused was the same kind of thing the check exists
    for, so it is worth a test of its own.
    """
    found = _load_checker()._dois_in_text(text)
    assert set(found) == {"10.5281/zenodo.24453044"}, (
        f"the DOI in {text!r} was read as {found}")


def test_every_doi_in_the_tree_is_registered():
    """No DOI in the tracked tree may point at a record that does not exist."""
    assert _CHECKER.is_file(), "scripts/check_dois.py is missing"

    proc = subprocess.run(
        [sys.executable, str(_CHECKER), "--quiet", "--delay", "0.3"],
        cwd=str(_REPO_ROOT), capture_output=True, text=True, timeout=900,
    )
    if proc.returncode == 2:
        pytest.skip(f"a DOI registry was unreachable: {proc.stdout[-400:]}")

    assert proc.returncode == 0, (
        "a DOI in the tracked tree does not resolve:\n"
        f"{proc.stdout[-1500:]}\n"
        "A DOI that resolves to nothing is a dead link a reader follows."
    )
