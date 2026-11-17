"""A local arXiv source bundle must be the manuscript it claims to be.

Why this file exists
--------------------
The arXiv source archive is a build product: `scripts/make_arxiv_bundle.py`
copies `paper/paper.tex`, flattens the three figures the manuscript reaches
outside the repository for, and writes `dist/deepcatch_v2_3_0.tar.gz`. It is
not a tracked file, so nothing kept it honest.

That matters because a stale bundle is not obviously stale. It compiles, it
produces the same eight pages, and it embeds every figure -- it just carries
the manuscript from an earlier commit, with that commit's numbers. The copy
on this machine had drifted that far: it still said "883 tests" long after
the manuscript said 919, and then 924.

What this test checks, and what it deliberately does not
-------------------------------------------------------
When a local bundle exists, it must be identical to a bundle built now from
the committed manuscript. The comparison uses the builder's own `build()`
function rather than a re-implementation of its rewriting rule, so the check
cannot drift away from the code that produces the archive.

It does not police the artifact CI uploads: the paper workflow rebuilds the
bundle from the committed manuscript on every run and publishes it as
`deepcatch-arxiv-source`. Where there is no local bundle -- CI, a fresh
clone -- there is nothing to check, and the test says so instead of failing.
"""
from __future__ import annotations

import importlib.util
import tarfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_BUNDLE = _REPO_ROOT / "dist" / "deepcatch_v2_3_0.tar.gz"
_BUILDER = _REPO_ROOT / "scripts" / "make_arxiv_bundle.py"


def _builder():
    """The bundle builder, loaded from its path."""
    spec = importlib.util.spec_from_file_location("make_arxiv_bundle", _BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _bundled_source(archive: Path) -> str:
    """The manuscript source inside an archive."""
    with tarfile.open(archive, "r:gz") as tar:
        names = [name for name in tar.getnames() if name.endswith(".tex")]
        assert len(names) == 1, (
            f"the bundle must hold exactly one manuscript, found {names}")
        handle = tar.extractfile(names[0])
        assert handle is not None, f"{names[0]} is not a readable member"
        return handle.read().decode("utf-8")


def test_the_local_arxiv_bundle_is_the_current_manuscript(tmp_path):
    """A local bundle must not carry an earlier commit's manuscript."""
    if not _BUNDLE.exists():
        pytest.skip(
            "no local arXiv bundle; CI builds one from the committed "
            "manuscript and publishes it as `deepcatch-arxiv-source`")

    builder = _builder()
    fresh_archive, _mapping = builder.build(tmp_path)

    assert _bundled_source(_BUNDLE) == _bundled_source(fresh_archive), (
        "the local arXiv bundle is stale: it was built from an earlier "
        "manuscript. Rebuild it with `python3 scripts/make_arxiv_bundle.py "
        "--verify`, or download the `deepcatch-arxiv-source` artifact from "
        "the latest 'Compile Paper PDF' run."
    )


def test_the_bundled_manuscript_carries_the_committed_test_count():
    """The bundle's numbers must be the committed manuscript's numbers.

    A cheaper, independent check on the same thing: the count is the number
    most likely to drift without anyone noticing, because it changes
    whenever a test file is added.
    """
    if not _BUNDLE.exists():
        pytest.skip("no local arXiv bundle")

    builder = _builder()
    manuscript = builder._MANUSCRIPT.read_text(encoding="utf-8")
    bundled = _bundled_source(_BUNDLE)

    for count in _stated_counts(manuscript):
        assert count in bundled, (
            f"the bundle does not state {count} tests; the committed "
            f"manuscript does. The bundle is from an earlier commit.")
    for count in _stated_counts(bundled):
        assert count in manuscript, (
            f"the bundle states {count} tests, which the committed "
            f"manuscript does not.")


def _stated_counts(text: str) -> list[str]:
    """The 'N tests' figures a text states."""
    import re

    return re.findall(r"(\d{3})\s+tests\b", text)
