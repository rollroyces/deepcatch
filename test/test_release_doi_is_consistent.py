"""The release DOI has to land in every place that cites the release.

Why this file exists
--------------------
`ZENODO_SETUP.md` says the minted DOI goes into three places besides
`CITATION.cff`: `MODEL.md`'s BibTeX entry, the manuscript's data and code
availability section, and the methods manuscript's availability section.
Four files, one DOI, one hand -- and nothing checked that they agree.

A half-done release is worse than no release: a reader who follows the
citation metadata to a Zenodo record and a reader who follows the
manuscript get different answers, and one of them is a link that
resolves to nothing. The failure is invisible to every other guard here:
the figure guards pin numbers, the count guards pin sizes, and a DOI
that appears in one file and not another is neither.

What this test checks, and what it deliberately does not
-------------------------------------------------------
* If `CITATION.cff` declares a DOI, the model card's BibTeX entry, the
  manuscript's availability section and the methods manuscript's
  availability section must contain that same DOI.
* None of the four may contain a placeholder DOI. The repository once
  carried `10.5281/zenodo.to-be-assigned-on-first-release`; it was
  removed, and this is what keeps it from coming back.
* A DOI that appears in one of the availability statements but in
  neither the citation metadata nor the list of public archives this
  repository cites is a DOI nobody declared -- most likely a release
  that was started and not finished.

The two public archives are named explicitly, because a data archive's
DOI belongs in the availability statement without being this project's
citation.
"""
from __future__ import annotations

import re
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent

_CITATION = _REPO_ROOT / "CITATION.cff"
_MODEL = _REPO_ROOT / "MODEL.md"
_TEX = _REPO_ROOT / "paper" / "paper.tex"
_METHODS = _REPO_ROOT / "paper" / "METHODS_PAPER.md"

# A DOI-shaped string.
_DOI = re.compile(r"10\.\d{4,9}/[^\s\"'<>,\\]+")

# A DOI that is a promise rather than a record.
_PLACEHOLDER = re.compile(
    r"10\.\d{4,9}/[^\s\"'<>,\]]*"
    r"(to-be-assigned|to-be-minted|to-be-added|todo|tbd|xxx|placeholder)",
    re.I,
)

# The public archives this repository cites in its availability statements.
# They are data, not this project's citation.
_PUBLIC_ARCHIVES = frozenset({
    "10.5281/zenodo.14013719",   # FinaleDB
    "10.5281/zenodo.6914806",    # CRAG/Jiang chr22 BAMs
})


def _section(text: str, start: str, end: str) -> str:
    """The text between two markers, inclusive of the start marker."""
    i = text.find(start)
    if i == -1:
        return ""
    j = text.find(end, i + len(start))
    return text[i:j] if j != -1 else text[i:]


def _citation_doi() -> str | None:
    """The DOI the citation metadata declares, if any."""
    text = _CITATION.read_text(encoding="utf-8")
    match = re.search(r'^doi:\s*"?([^"\n]+?)"?\s*$', text, re.M)
    return match.group(1) if match else None


def _model_dois() -> list[str]:
    """The DOIs in the model card's BibTeX entry."""
    text = _MODEL.read_text(encoding="utf-8")
    block = _section(text, "```bibtex", "```\n")
    return _DOI.findall(block)


def _tex_dois() -> list[str]:
    """The DOIs in the manuscript's data and code availability section."""
    text = _TEX.read_text(encoding="utf-8")
    section = _section(text, r"\section*{Data and code availability}",
                       r"\section*{")
    return _DOI.findall(section)


def _methods_dois() -> list[str]:
    """The DOIs in the methods manuscript's availability section."""
    text = _METHODS.read_text(encoding="utf-8")
    section = _section(text, "## 5. Availability", "\n## ")
    return _DOI.findall(section)


def test_the_release_doi_reaches_every_place_that_cites_it():
    """A declared DOI must appear in the model card and both manuscripts."""
    declared = _citation_doi()
    if declared is None:
        # No DOI has been minted yet; the other guards make sure a
        # placeholder cannot stand in for one.
        return

    for name, found in (
        ("MODEL.md's BibTeX entry", _model_dois()),
        ("paper/paper.tex's availability section", _tex_dois()),
        ("paper/METHODS_PAPER.md's availability section", _methods_dois()),
    ):
        assert declared in found, (
            f"{name} does not cite the DOI that CITATION.cff declares "
            f"({declared}). A release that cites a DOI in one place and not "
            "in another sends readers to two different answers."
        )


def test_no_placeholder_doi_survives_in_the_citation_metadata():
    """A promised DOI is not a DOI."""
    for name, path in (
        ("CITATION.cff", _CITATION),
        ("MODEL.md", _MODEL),
        ("paper/paper.tex", _TEX),
        ("paper/METHODS_PAPER.md", _METHODS),
    ):
        text = path.read_text(encoding="utf-8")
        placeholders = _PLACEHOLDER.findall(text)
        assert not placeholders, (
            f"{name} still carries a placeholder DOI. A placeholder renders "
            "as a link to a record that does not exist; declare no DOI until "
            "one is minted."
        )


def test_a_doi_in_an_availability_statement_is_declared_somewhere():
    """A DOI nobody declared is a release that was not finished."""
    declared = _citation_doi()
    known = set(_PUBLIC_ARCHIVES)
    if declared:
        known.add(declared)

    for name, found in (
        ("paper/paper.tex's availability section", _tex_dois()),
        ("paper/METHODS_PAPER.md's availability section", _methods_dois()),
        ("MODEL.md's BibTeX entry", _model_dois()),
    ):
        undeclared = [doi for doi in found if doi not in known]
        assert not undeclared, (
            f"{name} cites a DOI that neither CITATION.cff nor the public "
            "archives declare: " + ", ".join(undeclared)
            + "\nEither the release is half done (add it to CITATION.cff) "
              "or the DOI is wrong."
        )
