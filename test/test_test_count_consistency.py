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
#
# `docs/PREPRINT_ABSTRACT.md` is NOT in this list, and no longer states a
# count: it is generated from the manuscript's abstract
# (`scripts/generate_preprint_abstract.py`), and the manuscript's abstract
# does not quote the suite size. It used to be here because the document
# was separately written and carried "All code (820 tests)" -- a sentence
# the paper itself does not make. The guard's rule is that documents
# which state a count must state the same one, not that every document
# must state one.
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


def test_the_readme_suite_table_matches_reality():
    """The README's per-suite test table must be the counts, not old ones.

    The README states how many tests each suite contributes and what
    the fast-path gate yields. Those numbers drifted (the biomedical-
    review suite said 24 when it had 37; the gate said "233 collected"
    when it collected 540), and nothing checked them -- the agreement
    test above only compares the documents with each other. This
    collects each suite and compares.
    """
    import subprocess
    import sys

    def collected(*args) -> int:
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", *args, "--collect-only", "-q",
             "-o", "addopts="],
            cwd=_REPO_ROOT, capture_output=True, text=True, timeout=600,
        )
        if proc.returncode != 0:
            pytest.skip(f"collection failed here: {proc.stderr[-300:]}")
        for line in proc.stdout.splitlines():
            match = re.search(r"(\d+) tests? collected", line)
            if match:
                return int(match.group(1))
        return 0

    readme = (_REPO_ROOT / "README.md").read_text(encoding="utf-8")

    suites = {
        "`test/test_publication_readiness.py`": "test/test_publication_readiness.py",
        "`test/test_biomedical_review_fixes.py`": "test/test_biomedical_review_fixes.py",
        "`test/test_foundation_smoke.py`": "test/test_foundation_smoke.py",
        "`src/foundation/test_integration.py`": "src/foundation/test_integration.py",
    }
    for label, path in suites.items():
        stated = re.search(
            r"\|\s*" + re.escape(label) + r"\s*\|\s*(\d+)\s*\|", readme)
        assert stated, f"the README no longer states a count for {label}"
        real = collected(path)
        assert int(stated.group(1)) == real, (
            f"the README says {label} has {stated.group(1)} tests; it "
            f"collects {real}")

    # The fast-path gate: the README states both the collected total and
    # the outcome of running it.
    total = collected("test/", "src/foundation/test_integration.py")
    stated_total = re.search(
        r"\*\*Combined fast-path gate\*\*\s*\|\s*\*\*(\d+)\s+collected\*\*", readme)
    assert stated_total, "the README no longer states the fast-path gate total"
    assert int(stated_total.group(1)) == total, (
        f"the README says the fast-path gate collects "
        f"{stated_total.group(1)}; it collects {total}")


def test_the_stated_count_matches_the_real_collection():
    """The stated count must be the count, not a number that used to be.

    Every document agreeing on a wrong number is still wrong, so this
    test pins the AGREEMENT against reality: it collects the suite and
    compares what the documents claim with what is actually collected.

    It also protects the other direction. The count is updated by
    replacing a string across the documents, and a raw replacement
    silently rewrites any longer number that contains it -- that is how
    the external-validation figure 0.7097 became 0.7237 in two
    successive count updates. Collection is the ground truth, so this
    test fails the moment the two diverge.
    """
    import subprocess
    import sys

    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "src", "test", "--collect-only", "-q",
         "-o", "addopts="],
        cwd=_REPO_ROOT, capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        pytest.skip(f"collection failed here: {proc.stderr[-300:]}")
    real = None
    for line in proc.stdout.splitlines():
        match = re.search(r"(\d+) tests? collected", line)
        if match:
            real = int(match.group(1))
            break
    assert real is not None, (
        "could not read the collected count from pytest's output:\n"
        f"{proc.stdout[-500:]}")

    stated = set()
    for path in COUNT_DOCS:
        stated |= _stated_counts(path)
    assert stated == {real}, (
        f"the documents state {sorted(stated)} tests but the suite collects "
        f"{real}. Update every place that states the count (and the CI "
        f"split) together.")
