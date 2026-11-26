#!/usr/bin/env python3
"""Measure the JOSS paper's length, section by section.

Why this file exists
--------------------
JOSS's author documentation says "The paper should be between 750-1750
words" (https://joss.readthedocs.io/en/latest/paper.html). The JOSS
manuscript in this repository is longer than that, and the number of
words it had was counted by hand -- so it was the number most likely to
be wrong, and the one a reviewer counts first.

This script measures it by one rule, the same rule every time: the
front matter (the YAML block) and the markdown markup are not prose and
are not counted; everything a reader sees in a section is. The
per-section numbers are printed as well, because that is where a cut has
to be decided.

Usage
-----
    python3 scripts/count_joss_paper_words.py           # the counts
    python3 scripts/count_joss_paper_words.py --check   # fail if over the limit
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper.md"

# JOSS's guidance, from its own author documentation (checked 2026-11-26):
# "The paper should be between 750-1750 words."
MIN_WORDS = 750
MAX_WORDS = 1750


def _strip_front_matter(text: str) -> str:
    """The paper without its YAML front matter."""
    if not text.startswith("---"):
        return text
    end = text.find("\n---", 3)
    return text[end + 4:] if end != -1 else text


def _plain(text: str) -> str:
    """The text a reader sees: no emphasis markers, no link targets."""
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = text.replace("**", "").replace("*", "").replace("`", "")
    return text


def sections(text: str) -> dict[str, str]:
    """The paper's sections, keyed by heading."""
    body = _strip_front_matter(text)
    found: dict[str, str] = {}
    current = None
    lines: list[str] = []
    for line in body.splitlines():
        match = re.match(r"^#\s+(.*?)\s*$", line)
        if match:
            if current is not None:
                found[current] = "\n".join(lines)
            current = match.group(1)
            lines = []
            continue
        lines.append(line)
    if current is not None:
        found[current] = "\n".join(lines)
    return found


def section_words(text: str) -> dict[str, int]:
    """The word count of every section, and of the whole paper."""
    counts = {
        name: len(_plain(body).split())
        for name, body in sections(text).items()
    }
    counts["(total)"] = sum(counts.values())
    # The reference list is a separate block: JOSS's own format builds it
    # automatically from the .bib file, below a "# References" heading
    # that is a placeholder rather than prose. The guidance's word count
    # is about the paper, so the prose is what is measured against it --
    # and both numbers are printed, because the difference is 228 words.
    counts["(prose)"] = sum(
        count for name, count in counts.items()
        if name not in {"(total)", "References"})
    return counts


def main() -> int:
    text = PAPER.read_text(encoding="utf-8")
    counts = section_words(text)
    total = counts.pop("(total)")
    prose = counts.pop("(prose)")
    for name, count in counts.items():
        print(f"{count:5d}  {name}")
    print(f"{prose:5d}  PROSE (everything but the reference list)")
    print(f"{total:5d}  TOTAL (prose + the reference list)")
    if prose > MAX_WORDS:
        print(f"\nJOSS's guidance: {MIN_WORDS}-{MAX_WORDS} words "
              f"({prose - MAX_WORDS} over)")
    else:
        print(f"\nJOSS's guidance: {MIN_WORDS}-{MAX_WORDS} words "
              f"({MAX_WORDS - prose} to spare)")

    if "--check" in sys.argv and not (MIN_WORDS <= prose <= MAX_WORDS):
        print(f"\n!! the paper's prose is {prose} words; JOSS asks for "
              f"{MIN_WORDS}-{MAX_WORDS}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
