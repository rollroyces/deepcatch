"""The documents must agree on how many tests the suite has.

What went wrong
---------------
The repository stated its test count in fourteen places -- the
manuscript, the methods paper, the README badge, the reviewer note, the
JOSS readiness table. Those copies drifted apart (the CI split said 395
core while the total said 581), and every time a test file was added the
number had to be chased by hand. Nothing failed when they disagreed.

This test does not pin the count -- it cannot, because the count moves
whenever a test is added. It pins the AGREEMENT: every document that
states a count must state the same one, and the per-job split must add
up. A test file added without updating the documents fails here.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Documents that state the suite's size. Each must quote the same total.
COUNT_DOCS = [
    "paper/paper.tex",
    "paper/METHODS_PAPER.md",
    "paper/PAPER.md",
    "paper.md",
    "README.md",
    "REVIEWERS.md",
    "NEXT_STEPS.md",
    "TEAM.md",
    "docs/JOURNAL_SUBMISSION.md",
    "docs/PORTFOLIO_ROADMAP.md",
    "docs/PATH_TO_IMPACT.md",
    "docs/PREPRINT_ABSTRACT.md",
]

# "658 tests", "658 collected", "Tests-658%20collected", "(658 tests)"
_COUNT_PATTERNS = (
    r"(\d{3})\s+tests\b",
    r"(\d{3})\s+collected\b",
    r"Tests-(\d{3})%20",
    r"\((\d{3})\s+tests\)",
)

# The README also describes a narrower gate ("Combined fast-path gate |
# 233 collected"), which is a different number about a different set of
# tests. Only its badge counts as the suite's size.
_DOC_PATTERNS = {
    "README.md": (r"Tests-(\d{3})%20",),
}


def _stated_counts(path: str) -> set[int]:
    text = (_REPO_ROOT / path).read_text(encoding="utf-8")
    found = set()
    for pattern in _DOC_PATTERNS.get(path, _COUNT_PATTERNS):
        found.update(int(m) for m in re.findall(pattern, text))
    return found


@pytest.mark.parametrize("path", COUNT_DOCS)
def test_document_states_the_test_count(path):
    counts = _stated_counts(path)
    assert counts, (
        f"{path} does not state the suite's test count; every document that "
        f"quotes it must be updated together"
    )


def test_all_documents_quote_the_same_test_count():
    per_doc = {path: _stated_counts(path) for path in COUNT_DOCS}
    totals = set()
    for path, counts in per_doc.items():
        totals |= counts
    assert len(totals) == 1, (
        "the documents disagree on the test count: "
        + "; ".join(f"{p} says {sorted(c)}" for p, c in per_doc.items())
    )


def test_the_ci_split_adds_up_to_the_total():
    """The two CI jobs must account for every test."""
    total = None
    splits = []
    for path in COUNT_DOCS:
        text = (_REPO_ROOT / path).read_text(encoding="utf-8")
        for match in re.finditer(
                r"(\d{3})\s+tests?[^.]*?\((\d{3})\s+core\s*\+\s*(\d{3})\s+deep",
                text):
            total = int(match.group(1))
            splits.append((int(match.group(2)), int(match.group(3))))
    if not splits:
        pytest.skip("no document states the per-job split")
    for core, deep in splits:
        assert core + deep == total, (
            f"the CI split {core} + {deep} does not add up to the stated "
            f"total {total}")
