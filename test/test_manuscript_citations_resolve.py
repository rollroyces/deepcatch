"""Every citation a manuscript makes must resolve to a reference.

Why this file exists
--------------------
`paper/METHODS_PAPER.md` carried numbered citations and no reference list.
The numbers it used were also inconsistent with each other: FinaleMe was
`[8]` in one sentence and `[11]` in four others, and the held-out
publications were `[18]`, `[19]`, `[20]` while the bibliography file
ordered them differently. A reader of a candidate PLOS/Bioinformatics
manuscript sees numbers pointing at nothing, and a reviewer who follows
one lands in the wrong paper.

The LaTeX manuscript had the opposite problem: it carries its references
inline, and its bibliography file is a different document entirely --
`paper/README.md` claimed the two were the same thing.

What this test checks, and what it deliberately does not
-------------------------------------------------------
* Every numbered citation in the methods manuscript resolves to an entry
  in its reference section, and every entry is cited at least once.
* Every `\\cite{}` key in the LaTeX manuscript resolves to a `\\bibitem`
  in that same manuscript -- the manuscript is self-contained, and the
  bib file next to it is not its bibliography.
* The compiled submission PDF carries no unresolved-citation marker. LaTeX
  prints `[?]` for a key it could not resolve, and that marker is visible
  on the page, so the PDF is checked as well as the source.

It does not check that the DOIs in the reference list resolve: that
requires the network, and a bibliography whose identifiers have drifted is
a different failure (checked when the bibliography was written).
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_METHODS = _REPO_ROOT / "paper" / "METHODS_PAPER.md"
_TEX = _REPO_ROOT / "paper" / "paper.tex"
_PDF = _REPO_ROOT / "paper" / "biorxiv_submission_v2.3.0.pdf"

# A numbered citation: [3], [3,4], [3-5].
_CITATION = re.compile(r"\[(\d+(?:\s*[–-]\s*\d+)?(?:\s*,\s*\d+)*)\]")

# A numbered reference-list entry: "12. Author ...".
_ENTRY = re.compile(r"^(\d+)\.\s+\S", re.M)


def _methods_citations() -> set[int]:
    """Every reference number the methods manuscript cites."""
    text = _METHODS.read_text(encoding="utf-8")
    # Only the body: the reference list itself carries the numbers.
    body = text.split("## 7. References")[0]
    cited: set[int] = set()
    for match in _CITATION.finditer(body):
        for part in re.split(r"\s*,\s*", match.group(1)):
            if "–" in part or "-" in part:
                low, high = re.split(r"\s*[–-]\s*", part)
                cited.update(range(int(low), int(high) + 1))
            else:
                cited.add(int(part))
    return cited


def _methods_entries() -> set[int]:
    """Every reference number the methods manuscript lists."""
    text = _METHODS.read_text(encoding="utf-8")
    if "## 7. References" not in text:
        return set()
    section = text.split("## 7. References", 1)[1]
    return {int(m.group(1)) for m in _ENTRY.finditer(section)}


def test_the_methods_manuscript_has_a_reference_list():
    """A manuscript with numbered citations needs the list they point at."""
    entries = _methods_entries()
    assert entries, (
        "paper/METHODS_PAPER.md cites numbered references but has no "
        "reference section. The numbers point at nothing."
    )


def test_every_cited_reference_resolves():
    """Every number the manuscript cites must be in the reference list."""
    cited = _methods_citations()
    entries = _methods_entries()
    assert cited, "the methods manuscript cites nothing -- the extractor is wrong"

    missing = sorted(cited - entries)
    assert not missing, (
        "the methods manuscript cites reference numbers that its reference "
        f"list does not contain: {missing}. Every number must resolve, and "
        "the numbers must be the ones the list uses."
    )


def test_every_listed_reference_is_cited():
    """A reference list entry nobody cites is clutter."""
    uncited = sorted(_methods_entries() - _methods_citations())
    assert not uncited, (
        f"the methods manuscript's reference list contains entries it never "
        f"cites: {uncited}. Either cite them or drop them."
    )


def test_the_latex_manuscript_cites_resolve_against_its_own_bibliography():
    """The LaTeX manuscript is self-contained: its keys are its bibitems."""
    tex = _TEX.read_text(encoding="utf-8")
    cited: set[str] = set()
    for match in re.finditer(r"\\cite[a-zA-Z]*\{([^}]*)\}", tex):
        cited.update(key.strip() for key in match.group(1).split(","))
    bibitems = set(re.findall(r"\\bibitem\{([^}]*)\}", tex))

    assert cited, "the LaTeX manuscript cites nothing -- the extractor is wrong"
    missing = sorted(cited - bibitems)
    assert not missing, (
        "paper/paper.tex cites keys it does not define: "
        f"{missing}. LaTeX prints '[?]' for each of them, on the page."
    )


def test_the_submission_pdf_has_no_unresolved_citation():
    """The compiled PDF must not carry a '[?]' where a reference should be."""
    pypdf = pytest.importorskip("pypdf", reason="pypdf is not installed here")
    reader = pypdf.PdfReader(str(_PDF))
    text = "".join(page.extract_text() or "" for page in reader.pages)
    markers = re.findall(r"\[\s*\?\s*\]", text)
    assert not markers, (
        f"the submission PDF carries {len(markers)} unresolved citation "
        "marker(s). LaTeX could not resolve the key; a reviewer sees '[?]'."
    )
