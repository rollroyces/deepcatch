"""A version must be pointed at by its release tag, not by a commit hash.

Why this file exists
--------------------
The model card said "DeepCatch v2.3.0 (commit `7f9613b` on `main`)" --
and `7f9613b` was 134 commits behind the release. A commit hash in a
version pointer is a promise that cannot be kept: `main` moves, the
hash stays, and a reader who checks it out gets a tree whose numbers
are not the numbers the document states. Nothing failed, because a
hash is a hash.

A tag cannot drift. The release is `v2.3.0`, the tag points at the
commit that was archived and minted into a DOI, and that pairing is
permanent.

What this test checks
---------------------
* No current-facing document may point at a version with a commit
  hash. The version must be named by its tag (and, where it helps, its
  DOI).
* The tag the documents name must exist in the repository -- a pointer
  to a tag that was never cut is the same failure in the other
  direction. This half needs the tags, which a shallow checkout does
  not have; it skips there rather than failing.

The historical records are excluded: the changelog and the dated audit
reports name old commits on purpose, because describing the old state
is the point of those files.
"""
from __future__ import annotations

import importlib.util
import io
import re
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))


def _load_numbers_module():
    spec = importlib.util.spec_from_file_location(
        "benchmark_numbers", _REPO_ROOT / "scripts" / "benchmark_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_numbers = _load_numbers_module()

# A version immediately followed by a commit hash: "2.3.0 (commit
# `7f9613b` on `main`)". The window is short so that a version and an
# unrelated commit mention further apart are not paired.
_VERSION_COMMIT = re.compile(
    r"v?\d+\.\d+\.\d+[^.\n]{0,40}?commit\s+[`'\"]?[0-9a-f]{7,40}[`'\"]?")

# The release tag a document points at: "release tag `v2.3.0`".
_TAG_MENTION = re.compile(r"release tag\s+`?(v\d+\.\d+\.\d+)`?")

_DOC_SUFFIXES = (".md", ".tex", ".html", ".cff")

# The scanner must not read itself: this file names the versions it
# looks for.
_THIS_FILE = "test/test_version_pointers_name_a_release_tag.py"


def _git(*args: str) -> str | None:
    try:
        proc = subprocess.run(["git", "-C", str(_REPO_ROOT), *args],
                              capture_output=True, text=True)
    except (OSError, ValueError):
        return None
    return proc.stdout if proc.returncode == 0 else None


def _tracked(suffixes: tuple) -> list[str]:
    out = _git("ls-files", "-z")
    if out is None:
        return []
    names = [n for n in out.split("\0") if n]
    return [n for n in names if n.lower().endswith(suffixes)]


def _text(rel: str) -> str:
    try:
        return io.open(_REPO_ROOT / rel, encoding="utf-8").read()
    except (OSError, UnicodeDecodeError):
        return ""


def test_no_document_points_at_a_version_with_a_commit_hash():
    """A commit hash in a version pointer is a promise that cannot be kept."""
    problems: list[str] = []
    for rel in _tracked(_DOC_SUFFIXES):
        if rel == _THIS_FILE or rel in _numbers.HISTORICAL_RECORDS:
            continue
        if rel.startswith("review/"):
            continue
        text = _text(rel)
        for match in _VERSION_COMMIT.finditer(text):
            problems.append(
                f"{rel}: {match.group(0)!r} -- name the release tag instead")
    assert not problems, (
        "a document points at a version with a commit hash, and a commit "
        "hash goes stale the moment the branch moves:\n  "
        + "\n  ".join(problems))


def test_every_release_tag_the_documents_name_exists():
    """A pointer to a tag that was never cut is the same failure."""
    tags = _git("tag", "-l")
    if not tags:
        pytest.skip("no tags in this checkout (a shallow clone has none)")
    known = {line.strip() for line in tags.splitlines() if line.strip()}
    mentioned: set[str] = set()
    for rel in _tracked(_DOC_SUFFIXES):
        if rel == _THIS_FILE:
            continue
        mentioned.update(_TAG_MENTION.findall(_text(rel)))
    missing = sorted(tag for tag in mentioned if tag not in known)
    assert not missing, (
        "these release tags are named in the documents but do not exist: "
        + ", ".join(missing))
