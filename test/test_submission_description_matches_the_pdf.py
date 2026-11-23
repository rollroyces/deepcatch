"""The documents' description of the submission PDF must be the PDF's.

Why this file exists
--------------------
Three documents describe the compiled submission: what it is, how many
pages it has and how many figures it carries. One of them said "2
figures" while the PDF embeds three -- the panel-LLR performance
figure, the pooled ROC and the sensitivity-at-specificity panel. A
count in prose about a binary artifact is exactly the kind of thing
that drifts: nothing regenerates it, and the reader has no way to tell
it is wrong without opening the PDF.

What this test checks
---------------------
* The truth, read from the artifacts: the page count from the committed
  PDF, and the figure count from the manuscript's own
  ``\\includegraphics`` calls (which is what the PDF embeds).
* Every document that states a page or figure count for the submission
  states the real ones.

The counts are checked per document rather than per sentence, so a
document that describes the submission in several places has to be
right in all of them.
"""
from __future__ import annotations

import io
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

_PDF = _REPO_ROOT / "paper" / "biorxiv_submission_v2.3.0.pdf"
_TEX = _REPO_ROOT / "paper" / "paper.tex"

# The documents that describe the compiled submission. The last one is
# gitignored (it carries the author's details), so it is only checked
# when it is present.
_DOCUMENTS = (
    "paper/README.md",
    "NEXT_STEPS.md",
    "docs/JOURNAL_SUBMISSION.md",
    "paper/BIORXIV_SUBMISSION.md",
)

_PAGES = re.compile(r"(\d+)\s+pages?\b")
_FIGURES = re.compile(r"(\d+)\s+figures?\b")
_INCLUDEGRAPHICS = re.compile(r"\\includegraphics")


def _manuscript_figures() -> int:
    """How many figures the manuscript includes."""
    text = _TEX.read_text(encoding="utf-8")
    return len(_INCLUDEGRAPHICS.findall(text))


def _pdf_pages() -> int:
    """How many pages the committed submission PDF has."""
    pypdf = pytest.importorskip("pypdf", reason="pypdf is not installed here")
    return len(pypdf.PdfReader(str(_PDF)).pages)


def _pdf_images() -> int:
    """How many images the committed submission PDF embeds."""
    pypdf = pytest.importorskip("pypdf", reason="pypdf is not installed here")
    reader = pypdf.PdfReader(str(_PDF))
    return sum(len(page.images) for page in reader.pages)


def test_the_pdf_embeds_exactly_the_figures_the_manuscript_includes():
    """The manuscript and the PDF must agree on how many figures there are."""
    manuscript = _manuscript_figures()
    embedded = _pdf_images()
    assert manuscript == embedded, (
        f"the manuscript includes {manuscript} figure(s) and the PDF embeds "
        f"{embedded}; rebuild it with:\n"
        "    cd paper && tectonic -X compile paper.tex && "
        "cp paper.pdf biorxiv_submission_v2.3.0.pdf")


def test_every_description_of_the_submission_states_the_real_counts():
    """A page or figure count in prose has to be the artifact's."""
    pages = _pdf_pages()
    figures = _manuscript_figures()
    problems: list[str] = []
    for rel in _DOCUMENTS:
        path = _REPO_ROOT / rel
        if not path.is_file():
            continue
        text = io.open(path, encoding="utf-8").read()
        for match in _PAGES.finditer(text):
            if int(match.group(1)) != pages:
                problems.append(
                    f"{rel} says {match.group(1)} pages; the PDF has {pages}")
        for match in _FIGURES.finditer(text):
            if int(match.group(1)) != figures:
                problems.append(
                    f"{rel} says {match.group(1)} figures; the manuscript "
                    f"includes {figures}")
    assert not problems, "\n  ".join(problems)
