"""Every document that quotes a benchmark figure must be protected.

Why this file exists
---------------------
The number guard (``test/test_paper_numbers_match_artifacts.py``) checks
that every document NAMED in ``scripts/benchmark_numbers.py`` still quotes
the figure it was named for. It does not check the other direction: a
document that quotes a figure but was never named is a copy that can
drift silently -- the artifact is regenerated, the named documents are
updated, and this one keeps the old number.

That gap was not theoretical. ``EXTRA_DOCS`` was a dict literal, and a
dict literal with a repeated key keeps only the LAST list, so seven
figures had a second entry that silently replaced the first. Thirty-one
document/figure protections were lost that way -- including the
documents the superseded-figure guard is scoped to, so a superseded
number could have come back into them unnoticed.

So there are two checks here:

1. No Python file in the repository may hold a dict literal with a
   duplicate key: the duplicate is silently lossy, and nothing about it
   looks like a mistake when reading the source.
2. Every tracked document that quotes a pinned figure must be either
   NAMED for that figure (and so held to it) or EXCUSED with a reason
   (and so held to keeping that reason true).

The audit is deliberately narrow. A figure whose three-decimal rendering
is shared with another pinned figure is skipped: a document that writes
"0.700" cannot be told apart from a document that quotes the
sensitivity-at-99 figure, and pinning it would make the guard fail for
a reason that has nothing to do with drift.
"""
from __future__ import annotations

import ast
import collections
import importlib.util
import io
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

# The scanner must not read itself: this file names the figures it looks
# for, and it is tracked, so without the exclusion it would always find
# an unprotected copy.
_THIS_FILE = "test/test_every_figure_copy_is_protected.py"

# The document types that quote research figures.
_DOC_SUFFIXES = (".md", ".tex", ".html", ".cff")


def _git(*args: str) -> str:
    try:
        proc = subprocess.run(["git", "-C", str(_REPO_ROOT), *args],
                              capture_output=True, text=True)
    except (OSError, ValueError) as exc:
        pytest.skip(f"git is not usable here: {exc}")
    if proc.returncode != 0:
        pytest.skip(f"git is not usable here: {proc.stderr.strip()}")
    return proc.stdout


def _tracked(suffixes: tuple) -> list[str]:
    """The COMMITTED files with one of the given suffixes.

    ``git ls-files`` reports the committed tree, not the working
    directory: a file that has just been written is invisible to it
    until it is committed. That is the right behaviour for this audit
    (it must see what a reader of the repository sees) but it means a
    new document is only covered once it is committed.
    """
    names = [n for n in _git("ls-files", "-z").split("\0") if n]
    return [n for n in names if n.lower().endswith(suffixes)]


def _text(rel: str) -> str:
    try:
        return io.open(_REPO_ROOT / rel, encoding="utf-8").read()
    except (OSError, UnicodeDecodeError):
        return ""


def _entries():
    """(label, value, artifact, key, documents) for every pinned figure."""
    return list(_numbers.quoted_entries())


def _figures() -> dict:
    """{key: (label, value, named documents)} for every DISTINCT figure.

    A figure can be pinned by several entries -- a number by one, and a
    lookup by one entry per document that quotes it -- so the entries
    are merged here. Otherwise the audit would check the same figure
    once per entry that names it, and report the same document many
    times over.
    """
    figures: dict = {}
    for label, value, artifact, key, documents in _entries():
        entry = figures.setdefault(key, [label, value, set()])
        entry[2].update(documents)
    return {key: (value[0], value[1], value[2])
            for key, value in figures.items()}


def _ambiguous_renderings() -> set[str]:
    """The three-decimal renderings more than one pinned figure shares.

    A document that writes "0.700" cannot be told apart from one that
    quotes the sensitivity-at-99 figure, and a document that writes
    "0.9444" cannot be told apart from one that quotes the other figure
    with the same value. Requiring such a document to be named for
    every figure that shares the digits would fail the build for a
    reason that is not drift.

    The comparison is between DISTINCT figures: one figure can be named
    for a dozen documents, and that is not a second figure sharing the
    digits.
    """
    counts = collections.Counter(
        f"{value:.3f}" for _label, value, _docs in _figures().values())
    return {rendering for rendering, count in counts.items() if count > 1}


def test_no_python_file_holds_a_duplicate_dict_literal_key():
    """A repeated dict key silently keeps only the last value.

    Nothing about it looks like a mistake when reading the source: the
    entry reads as an ordinary key. The cost is invisible too -- the
    earlier list is simply gone. This makes the loss loud.
    """
    offenders: list[str] = []
    for rel in _tracked((".py",)):
        source = _REPO_ROOT / rel
        try:
            tree = ast.parse(source.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Dict):
                continue
            keys = []
            for key in node.keys:
                if key is None:
                    continue
                try:
                    keys.append(ast.literal_eval(key))
                except (ValueError, TypeError, SyntaxError):
                    continue
            counts = collections.Counter(
                k for k in keys if isinstance(k, (str, int, float, bool, tuple)))
            for key, count in counts.items():
                if count > 1:
                    offenders.append(
                        f"{rel}:{node.lineno} key {key!r} appears {count} times")
    assert not offenders, (
        "a dict literal with a duplicate key silently drops a value:\n  "
        + "\n  ".join(offenders))


def test_every_document_quoting_a_pinned_figure_is_protected():
    """A copy of a figure that is not named (nor excused) can drift."""
    ambiguous = _ambiguous_renderings()
    problems: list[str] = []
    for key, (label, value, documents) in _figures().items():
        if f"{value:.3f}" in ambiguous:
            continue
        excused = _numbers.EXCUSED_COPIES.get(key, {})
        for rel in _tracked(_DOC_SUFFIXES):
            if rel == _THIS_FILE or rel in documents or rel in excused:
                continue
            if rel in _numbers.HISTORICAL_RECORDS:
                continue
            if _numbers.quotes(_text(rel), value):
                problems.append(
                    f"{rel} quotes {value:.6f} ({label}) but is not named "
                    "for it and has no excuse")
    assert not problems, (
        "these documents quote a pinned figure without being protected "
        "from drift:\n  " + "\n  ".join(problems))


def test_every_excuse_is_still_needed_and_states_its_reason():
    """An excuse is a claim about a document, so it has to stay true.

    Three things can go wrong: the excuse names a figure that is not
    pinned at all, the document stops quoting the figure (so the excuse
    is stale), or the document is named for the figure after all (so the
    excuse is redundant and hides a protection).
    """
    figures = _figures()
    problems: list[str] = []
    for key, excuses in _numbers.EXCUSED_COPIES.items():
        if key not in figures:
            problems.append(f"the excuse for {key} names a figure that is "
                            "not pinned")
            continue
        label, value, documents = figures[key]
        for rel, reason in excuses.items():
            if not reason or not reason.strip():
                problems.append(f"{rel} is excused from {label} without "
                                "saying why")
                continue
            if len(reason.strip()) < 30:
                problems.append(
                    f"{rel} is excused from {label} with the reason "
                    f"{reason!r}, which does not say what the number is "
                    "instead")
                continue
            if rel in documents:
                problems.append(f"{rel} is excused from {label} but is "
                                "named for it; the excuse is redundant")
                continue
            if not _numbers.quotes(_text(rel), value):
                problems.append(f"{rel} is excused from {label} but no "
                                "longer quotes it; the excuse is stale")
    assert not problems, "\n  ".join(problems)
