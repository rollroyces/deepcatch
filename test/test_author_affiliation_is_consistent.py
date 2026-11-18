"""The three author blocks must not contradict each other.

Why this file exists
--------------------
A manuscript is submitted once, and the form asks for an affiliation; the
paper's author block and the repository's citation metadata are what a
reviewer sees side by side. `docs/JOURNAL_SUBMISSION.md` says they "must
not contradict each other on the day of upload" -- and nothing checked
them. The affiliation in this repository is genuinely undecided (the
public ORCID record carries one employment, the manuscripts declare none),
which is the author's fact to settle, not the repository's. But the
*consistency* of the three statements is the repository's business, and it
can be checked without deciding the question.

Two things are checked:

* The JOSS manuscript keeps the sections JOSS asks for. `docs/JOURNAL_SUBMISSION.md`
  claims it is "in JOSS format (summary, statement of need, features,
  quality/reproducibility, documentation, code availability, limitations,
  acknowledgements)"; a section removed in an edit would leave the claim
  standing and the submission incomplete.
* The three author blocks -- the LaTeX manuscript, the JOSS manuscript and
  `CITATION.cff` -- declare the same affiliation, and the ORCID that
  identifies the author is the same in all of them.

What this test checks, and what it deliberately does not
-------------------------------------------------------
It does not decide whether an institutional affiliation is correct. It
checks that the three files agree, so that whichever way the question is
settled, it is settled in all three places at once -- and so that a
partial edit cannot leave the paper saying one thing and the citation
metadata another.

The institution markers are a small fixed vocabulary rather than a
parser: the point is to notice when one file starts naming an employer
and the others do not, not to understand what an affiliation is.
"""
from __future__ import annotations

import re
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent

# The sections JOSS asks a submission to carry, in the order it lists them.
_JOSS_SECTIONS = (
    "Summary",
    "Statement of need",
    "Features",
    "Quality and reproducibility",
    "Documentation",
    "Code availability",
    "Limitations",
    "Acknowledgements",
)

# Words that mark an employer in an author block. A file that gains one
# of these must gain it in the other two as well.
_INSTITUTION_MARKERS = (
    "university",
    "institute",
    "hospital",
    "college",
    "school of",
    "department of",
    "centre for",
    "center for",
    "laboratory",
    "academy of",
)

_ORCID = re.compile(r"\b(\d{4}-\d{4}-\d{4}-\d{3}[\dX])\b")


def _joss_manuscript() -> str:
    return (_REPO_ROOT / "paper.md").read_text(encoding="utf-8")


def _front_matter(text: str) -> str:
    """The YAML front matter of the JOSS manuscript."""
    if not text.startswith("---"):
        return ""
    end = text.find("\n---", 3)
    return text[:end] if end != -1 else ""


def _latex_author_block() -> str:
    """The `\\author{...}` block of the LaTeX manuscript."""
    text = (_REPO_ROOT / "paper" / "paper.tex").read_text(encoding="utf-8")
    match = re.search(r"\\author\{(.*?)\}\s*\n", text, re.S)
    return match.group(1) if match else ""


def _citation_author_block() -> str:
    """The `authors:` entry of CITATION.cff."""
    text = (_REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^authors:\s*\n((?:\s+.*\n)+)", text, re.M)
    return match.group(1) if match else ""


def _markers(text: str) -> set[str]:
    """The institution markers a text names."""
    lowered = text.lower()
    return {marker for marker in _INSTITUTION_MARKERS if marker in lowered}


def test_the_joss_manuscript_keeps_the_sections_joss_asks_for():
    """The JOSS package must stay submittable."""
    headings = {
        match.group(1).strip()
        for match in re.finditer(r"^#\s+(.*)$", _joss_manuscript(), re.M)
    }
    missing = [section for section in _JOSS_SECTIONS if section not in headings]
    assert not missing, (
        "the JOSS manuscript is missing the sections JOSS asks for: "
        + ", ".join(missing)
        + "\n`docs/JOURNAL_SUBMISSION.md` states that the manuscript is in "
          "JOSS format; that claim is only true while the sections are "
          "there."
    )


def test_the_author_blocks_declare_the_same_affiliation():
    """An employer named in one author block must be named in all three."""
    blocks = {
        "paper/paper.tex": _latex_author_block(),
        "paper.md": _front_matter(_joss_manuscript()),
        "CITATION.cff": _citation_author_block(),
    }
    for name, block in blocks.items():
        assert block.strip(), f"{name} has no author block to check"

    per_file = {name: _markers(block) for name, block in blocks.items()}
    distinct = {frozenset(markers) for markers in per_file.values()}
    assert len(distinct) == 1, (
        "the author blocks disagree about the affiliation:\n"
        + "\n".join(
            f"  {name}: {sorted(markers) or 'no institution named'}"
            for name, markers in per_file.items()
        )
        + "\nWhatever the affiliation is, it has to be the same in the "
          "manuscript, the JOSS paper and the citation metadata -- a "
          "reviewer reads them together."
    )


def test_the_orcid_is_the_same_in_the_three_author_blocks():
    """The identifier that names the author cannot differ between files."""
    blocks = {
        "paper/paper.tex": _latex_author_block(),
        "paper.md": _front_matter(_joss_manuscript()),
        "CITATION.cff": _citation_author_block(),
    }
    orcids = {
        name: set(_ORCID.findall(block)) for name, block in blocks.items()
    }
    empty = [name for name, found in orcids.items() if not found]
    assert not empty, (
        "the author block of " + ", ".join(empty) + " carries no ORCID. "
        "The repository's citation metadata and the model card both carry "
        "it, and bioRxiv's form asks for it."
    )
    distinct = {frozenset(found) for found in orcids.values()}
    assert len(distinct) == 1, (
        "the three author blocks carry different ORCIDs:\n"
        + "\n".join(f"  {name}: {sorted(found)}" for name, found in orcids.items())
    )


def test_the_marker_vocabulary_notices_an_employer():
    """The marker check must see an employer when one is named."""
    assert _markers("Yu Ching Lam, University of Michigan") == {"university"}
    assert _markers("Independent Researcher") == set()
    # "centering" is not "center": the marker must not fire on a substring.
    assert _markers("\\centering") == set()
