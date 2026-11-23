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


def _markdown_author_block(path: str) -> str:
    """The author line and its footnote of a markdown manuscript."""
    text = (_REPO_ROOT / path).read_text(encoding="utf-8")
    lines = text.split("\n")
    for index, line in enumerate(lines):
        if line.startswith("**Authors:**"):
            block = [line]
            for following in lines[index + 1:index + 6]:
                if following.startswith("\u00b9") or following.startswith("1 "):
                    block.append(following)
                    break
            return "\n".join(block)
    return ""


def _authors_file_block() -> str:
    """The lead author entry of the AUTHORS file."""
    text = (_REPO_ROOT / "AUTHORS").read_text(encoding="utf-8")
    match = re.search(r"^- \*\*.+?\*\*.*$", text, re.M)
    return match.group(0) if match else ""


def _latex_author_name(block: str) -> str:
    """The name in the LaTeX author block: its first line."""
    return block.strip().split("\n")[0].strip().rstrip("\\").strip()


def _joss_author_name(block: str) -> str:
    """The name in the JOSS front matter."""
    match = re.search(r"name:\s*(.+)", block)
    return match.group(1).split(" (")[0].strip() if match else ""


def _markdown_author_name(block: str) -> str:
    """The name in a markdown manuscript's author line."""
    match = re.search(r"\*\*Authors:\*\*\s*(.+)", block)
    if not match:
        return ""
    return re.split(r"[\u00b9\u00b2]", match.group(1))[0].split(" (")[0].strip()


def _citation_author_name(block: str) -> str:
    """The name in the citation metadata's structured fields."""
    given = re.search(r'given-names:\s*"?([^"\n]+)"?', block)
    family = re.search(r'family-names:\s*"?([^"\n]+)"?', block)
    if given and family:
        return f"{given.group(1).strip()} {family.group(1).strip()}"
    return ""


def _authors_file_name(block: str) -> str:
    """The name in the AUTHORS file's lead entry."""
    match = re.search(r"^- \*\*(.+?)\*\*", block)
    return match.group(1).strip() if match else ""


# Every surface an index reads the author from, with the block that
# carries the author and the reader that finds the name in it.
#
# Google Scholar matches a paper to a profile by the NAME string, the
# affiliation and the email -- it does not read an ORCID to find an
# author -- and it follows the ORCID to link the work to the profile
# afterwards. So the surfaces that matter are the ones a preprint
# server, a journal form or an archive reads: the LaTeX manuscript, the
# JOSS paper, the two markdown manuscripts, the citation metadata
# (which is what Zenodo archives), and the author list.
_AUTHOR_SURFACES = {
    "paper/paper.tex": (_latex_author_block, _latex_author_name),
    "paper.md": (lambda: _front_matter(_joss_manuscript()), _joss_author_name),
    "paper/PAPER.md": (lambda: _markdown_author_block("paper/PAPER.md"),
                       _markdown_author_name),
    "paper/METHODS_PAPER.md": (lambda: _markdown_author_block("paper/METHODS_PAPER.md"),
                               _markdown_author_name),
    "CITATION.cff": (_citation_author_block, _citation_author_name),
    "AUTHORS": (_authors_file_block, _authors_file_name),
}


def _all_author_blocks() -> dict:
    """The author block of every surface an index reads."""
    return {name: read() for name, (read, _) in _AUTHOR_SURFACES.items()}


def test_the_orcid_is_the_same_in_every_author_block():
    """The identifier that names the author cannot differ between files."""
    blocks = _all_author_blocks()
    orcids = {
        name: set(_ORCID.findall(block)) for name, block in blocks.items()
    }
    empty = [name for name, found in orcids.items() if not found]
    assert not empty, (
        "the author block of " + ", ".join(empty) + " carries no ORCID. "
        "The repository's citation metadata and the model card both carry "
        "it, bioRxiv's form asks for it, and it is what lets an index "
        "attach the work to the author's profile."
    )
    distinct = {frozenset(found) for found in orcids.values()}
    assert len(distinct) == 1, (
        "the author blocks carry different ORCIDs:\n"
        + "\n".join(f"  {name}: {sorted(found)}" for name, found in orcids.items())
    )


def test_the_author_name_is_the_same_string_everywhere():
    """The name is what an index matches a paper to a profile by.

    A name that differs by a hyphen, an initial or a particle splits the
    author's works across two profiles, and the split is invisible until
    the works have to be merged by hand. The ORCID check above covers
    the other half -- the identifier that links them afterwards.
    """
    blocks = _all_author_blocks()
    names = {name: reader(blocks[name])
             for name, (read, reader) in _AUTHOR_SURFACES.items()}
    empty = [name for name, found in names.items() if not found]
    assert not empty, (
        "no author name could be read from " + ", ".join(empty))
    distinct = set(names.values())
    assert len(distinct) == 1, (
        "the surfaces spell the author's name differently:\n"
        + "\n".join(f"  {name}: {found!r}" for name, found in names.items())
    )


def test_the_marker_vocabulary_notices_an_employer():
    """The marker check must see an employer when one is named."""
    assert _markers("Yu Ching Lam, University of Michigan") == {"university"}
    assert _markers("Independent Researcher") == set()
    # "centering" is not "center": the marker must not fire on a substring.
    assert _markers("\\centering") == set()
