#!/usr/bin/env python3
"""Measure the methods manuscript's abstract, and check the documents that state it.

Why this file exists
--------------------
PLOS Computational Biology's submission guidelines say "The Abstract should
be succinct; it must not exceed 300 words", and three documents in this
repository state how long the methods manuscript's abstract is. All three
said **414 words** -- and the real number is **436**. The figure had been
measured by hand, wrongly, and nothing ever re-measured it, so it sat in
the row that says "this would come back for revision" and under-stated the
gap by 22 words.

A number in a row about a venue's limit is exactly the kind of thing the
count guard exists to catch, and it was not caught. So it is measured here,
by the same rule every time, and the documents that state it are checked
against the measurement.

The rule
--------
The abstract is the ``## Abstract`` section of ``paper/METHODS_PAPER.md`` --
the source the manuscript's LaTeX is generated from. Its words are counted
as the text a reader sees: the markdown emphasis markers are not words,
and the structured labels ("Background.", "Methods.", "Results.",
"Conclusions.") are, because they are part of the submitted text. The
generated LaTeX renders within two words of this -- the difference is LaTeX
symbols that become a single word -- so the figure is the same either way.

Usage
-----
    python3 scripts/count_manuscript_abstract.py         # print the count
    python3 scripts/count_manuscript_abstract.py --check  # fail if a document disagrees
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT = ROOT / "paper" / "METHODS_PAPER.md"

# PLOS Computational Biology's cap, from its own submission guidelines
# (https://journals.plos.org/ploscompbiol/s/submission-guidelines, checked
# 2026-11-26): "The Abstract should be succinct; it must not exceed 300
# words." Bioinformatics Advances, the other candidate venue for this
# manuscript, states no limit.
PLOS_WORD_LIMIT = 300

# The documents that state the count, and how each of them states it.
# Every one of them is required to state it: a check that only verified a
# number *if* it were there could be defeated by deleting it.
STATED_IN = {
    "CHANGELOG.md": r"abstract measures \*\*(\d+) words\*\*",
    "docs/JOURNAL_SUBMISSION.md": r"abstract measures \*\*(\d+) words\*\*",
    "paper/METHODS_PAPER_STATUS.md": r"abstract here measures (\d+)",
}


def abstract_text(path: Path = MANUSCRIPT) -> str:
    """The manuscript's abstract, as written."""
    text = path.read_text(encoding="utf-8")
    start = text.index("## Abstract") + len("## Abstract")
    end = text.index("\n---", start)
    return text[start:end].strip()


def word_count(text: str) -> int:
    """The number of words in the abstract.

    The markdown emphasis markers are not words; everything else in the
    text is, including the structured labels.
    """
    return len(text.replace("**", " ").split())


def measured(path: Path = MANUSCRIPT) -> int:
    """The abstract's measured word count."""
    return word_count(abstract_text(path))


def stated_counts() -> dict[str, list[int]]:
    """The word counts each document states, keyed by document."""
    found: dict[str, list[int]] = {}
    for rel, pattern in STATED_IN.items():
        path = ROOT / rel
        if not path.is_file():
            found[rel] = []
            continue
        text = path.read_text(encoding="utf-8")
        found[rel] = [int(m) for m in re.findall(pattern, text)]
    return found


def check() -> list[str]:
    """The problems the documents have, if any."""
    problems: list[str] = []
    count = measured()
    for rel, values in stated_counts().items():
        if not values:
            problems.append(
                f"{rel} does not state the methods abstract's word count; "
                f"state it ({count} words) or remove it from the guard")
            continue
        for value in values:
            if value != count:
                problems.append(
                    f"{rel} states {value} words; the abstract measures "
                    f"{count}")
    return problems


def main() -> int:
    count = measured()
    print(f"methods manuscript abstract: {count} words "
          f"(PLOS Comp Bio limit {PLOS_WORD_LIMIT}, "
          f"{count - PLOS_WORD_LIMIT} over)")
    problems = check()
    if problems:
        print()
        for problem in problems:
            print(f"  !! {problem}")
        return 1
    print("every document that states it states the measured count")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
