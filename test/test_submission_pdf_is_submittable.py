"""The submission PDF must be a PDF the preprint servers accept.

Why this file exists
--------------------
`paper/biorxiv_submission_v2.3.0.pdf` is the file bioRxiv receives and
`dist/deepcatch_v2_3_0.tar.gz` is the file arXiv receives. Both are
uploaded by a human, once, and the feedback when the server does not
like them is a rejection screen or a mangled PDF -- not a test failure.

There is a separate guard for whether the PDF is the manuscript
(`test_submission_pdf_matches_the_manuscript.py`: the figures and the
numbers). This one asks the other half of the question: whether the file
itself is one a preprint server can process.

What this test checks
---------------------
1. Every font in the PDF is embedded. An unembedded font is substituted
   when the PDF is rendered or converted, and the paper comes out with
   the wrong typeface -- or with missing glyphs.
2. The PDF is not encrypted. bioRxiv's submission help lists locked or
   encrypted PDFs among the formats it does not support.
3. Its text is extractable, so the server's conversion can read it (a
   scanned page cannot be proof-read or indexed).
4. It is small. arXiv's limit is 10 MB for the output PDF or for any
   individual file in a submission; this PDF is a fraction of that.

What it does not check
----------------------
The venue's other requirements -- the abstract's length, the figures'
format, the metadata -- are checked where they live: the preprint
abstract generator's own test, the manuscript's figure guard, and the
author's local submission checklist. This file covers the properties of
the PDF file itself.

Scope
-----
`pypdf` is a test-only dependency: the CI job that compiles the paper
installs it; where it is absent the test skips rather than fails.
"""
from __future__ import annotations

from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_PDF = _REPO_ROOT / "paper" / "biorxiv_submission_v2.3.0.pdf"

# arXiv's limit: 10 MB for the output PDF or for any individual file in
# the submission.
_MAX_PDF_BYTES = 10 * 1024 * 1024


def _reader():
    """The submission PDF, read with pypdf."""
    pypdf = pytest.importorskip(
        "pypdf", reason="pypdf is a test-only dependency")
    assert _PDF.exists(), f"{_PDF} is missing"
    return pypdf.PdfReader(str(_PDF))


def _fonts(reader) -> dict[str, bool]:
    """The PDF's fonts, name -> embedded?"""
    fonts: dict[str, bool] = {}
    for page in reader.pages:
        resources = page.get("/Resources") or {}
        for _, ref in (resources.get("/Font") or {}).items():
            font = ref.get_object()
            name = str(font.get("/BaseFont"))
            descriptor = font.get("/FontDescriptor")
            embedded = False
            if descriptor is not None:
                desc = descriptor.get_object()
                embedded = any(
                    key in desc for key in ("/FontFile", "/FontFile2", "/FontFile3")
                )
            if not embedded and "/DescendantFonts" in font:
                # A Type 0 composite: the descendant fonts carry the file.
                embedded = True
            fonts.setdefault(name, embedded)
    return fonts


def test_every_font_in_the_submission_pdf_is_embedded():
    """An unembedded font is substituted when the PDF is rendered."""
    fonts = _fonts(_reader())
    assert fonts, "the PDF has no fonts; the check would be vacuous"
    missing = sorted(name for name, embedded in fonts.items() if not embedded)
    assert not missing, (
        "the submission PDF has fonts that are not embedded: "
        + ", ".join(missing)
        + ". A font the reader's machine does not have is substituted, so "
          "the paper is typeset in the wrong typeface."
    )


def test_the_submission_pdf_is_not_encrypted():
    """bioRxiv does not accept locked or encrypted PDFs."""
    reader = _reader()
    assert not reader.is_encrypted, (
        "the submission PDF is encrypted; bioRxiv's submission help lists "
        "locked or encrypted PDFs among the formats it does not support"
    )


def test_the_submission_pdf_carries_text():
    """A scanned page cannot be proof-read, indexed or converted."""
    reader = _reader()
    text = "".join((page.extract_text() or "") for page in reader.pages)
    assert len(text.strip()) > 1000, (
        f"the submission PDF carries {len(text.strip())} characters of "
        "extractable text; it looks like an image, not a typeset paper"
    )


def test_the_submission_pdf_is_within_arxivs_size_limit():
    """arXiv's limit is 10 MB for the output PDF."""
    size = _PDF.stat().st_size
    assert size <= _MAX_PDF_BYTES, (
        f"the submission PDF is {size} bytes; arXiv's limit is 10 MB for "
        "the output PDF"
    )
