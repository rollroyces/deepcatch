"""The submission PDF must be the manuscript that is committed.

Why this file exists
--------------------
`paper/biorxiv_submission_v2.3.0.pdf` is the file bioRxiv receives, and
it is a committed binary: nothing compared it with `paper/paper.tex`.
A PDF that was not rebuilt after the manuscript changed would carry
every figure the manuscript has since changed, and nobody would see it
-- the drift guard checks the .tex against the artifacts, and the PDF
would simply be an older rendering of the same prose. This repository
has already had a submission PDF whose Figure 1 plotted a different
curve from the one the text described (CHANGELOG, v2.2.0).

The check reads the PDF's own text and requires every decimal figure
the manuscript states to be in it. The list is DERIVED from the .tex
rather than written out here, so the guard cannot itself drift from
the manuscript, and it grows a check for every new figure.

`pypdf` is a test-only dependency. The CI job installs it; where it is
absent the test skips rather than fails, and says so.
"""
from __future__ import annotations

import io
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_PDF = _REPO_ROOT / "paper" / "biorxiv_submission_v2.3.0.pdf"
_TEX = _REPO_ROOT / "paper" / "paper.tex"

# A figure the manuscript states: a decimal with three or four places.
# One or two places are too generic (a threshold written 0.5 in one
# place and 0.500 in another would be a false alarm), and the figures
# that matter -- AUCs, sensitivities, CIs -- all carry three or four.
_FIGURE = re.compile(r"(?<![\d.])\d+\.\d{3,4}(?![\d])")

# A three-digit count the manuscript states: cohort sizes, the CpG
# count, the test-suite size. Two-digit and one-digit numbers are left
# out (they are too generic to be a reliable signal), and four-digit
# ones are covered by the comma-formatted ones the manuscript writes.
# A digit string that is part of a longer number is not a count: "103,863"
# contributes 103, and "50{,}000" contributes nothing -- the lookbehind
# rejects a preceding comma, a closing brace and a letter, so the two
# renderings of the same quantity agree.
_COUNT = re.compile(r"(?<![\d.,}\w])\d{3}(?![\d.])")


@pytest.fixture(scope="module")
def pdf_text() -> str:
    """The submission PDF's text, whitespace collapsed."""
    # The skip lives HERE, not at module level: a module-level
    # importorskip would make the COLLECTION count depend on whether
    # pypdf is installed, and the repository's test-count guard compares
    # the collected count against the number the documents quote.
    pypdf = pytest.importorskip("pypdf", reason="pypdf is not installed here")
    assert _PDF.is_file(), f"{_PDF} is missing; it is the submission artefact"
    reader = pypdf.PdfReader(str(_PDF))
    raw = "\n".join((page.extract_text() or "") for page in reader.pages)
    return " ".join(raw.split())


def _manuscript_figures() -> list[str]:
    """Every figure and count the manuscript states, read out of the .tex."""
    tex = _TEX.read_text(encoding="utf-8")
    body = tex[tex.index("\\begin{document}"):]
    # Strip LaTeX comments -- but not an escaped percent sign, which is
    # how the abstract writes its percentages.
    body = re.sub(r"(?<!\\)%.*", "", body)
    body = re.sub(r"\\cite\{[^}]*\}", "", body)
    return sorted(set(_FIGURE.findall(body)) | set(_COUNT.findall(body)))


def _manuscript_images() -> list[Path]:
    """The figure files the manuscript includes, read out of the .tex.

    The paths in paper.tex are relative to the paper/ directory, so they
    are resolved against it.
    """
    tex = _TEX.read_text(encoding="utf-8")
    body = tex[tex.index("\\begin{document}"):]
    body = re.sub(r"(?<!\\)%.*", "", body)
    paths = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", body)
    return [(_REPO_ROOT / "paper" / path).resolve() for path in paths]


def test_the_submission_pdf_carries_every_figure_the_manuscript_states(
        pdf_text):
    """A figure that left the manuscript must have left the PDF too."""
    figures = _manuscript_figures()
    assert figures, "no figures were read out of paper.tex -- the extractor broke"

    # A number can be split across a line in the extracted text, so the
    # whitespace-free rendering is accepted as well.
    tight = pdf_text.replace(" ", "")

    missing = [f for f in figures
               if f not in pdf_text and f not in tight]
    assert not missing, (
        f"the committed submission PDF does not carry {len(missing)} figure(s) "
        f"the manuscript states: {', '.join(missing)}. The PDF is stale -- "
        f"rebuild it with:\n"
        f"    cd paper && tectonic -X compile paper.tex && "
        f"cp paper.pdf biorxiv_submission_v2.3.0.pdf\n"
        f"(or pdftex: pdflatex paper.tex && pdflatex paper.tex)")


def test_the_submission_pdf_is_the_current_manuscript(pdf_text):
    """The PDF's title, author and closing sentence are the current ones."""
    for phrase in (
        "An Open-Source Benchmark for Panel-Based",
        "Yu Ching Lam",
        "pooling to 0.7385",
        "Clinical validation on real plasma cfDNA remains the essential next step",
    ):
        assert phrase in pdf_text, (
            f"the submission PDF does not contain {phrase!r}; it is not the "
            f"current manuscript -- rebuild it from paper/paper.tex")


def test_the_submission_pdf_has_pages():
    """A submission artefact that renders nothing is not a submission."""
    pypdf = pytest.importorskip("pypdf", reason="pypdf is not installed here")
    reader = pypdf.PdfReader(str(_PDF))
    assert len(reader.pages) >= 1


def test_the_submission_pdf_embeds_the_committed_figures():
    """The figures in the PDF must be the figures the repository holds.

    A figure that was regenerated without rebuilding the PDF leaves the
    submission showing a curve the text no longer describes -- that
    happened once in this repository (CHANGELOG, v2.2.0), and nothing
    caught it because the PDF is a binary.

    The comparison is pixel-for-pixel: LaTeX re-encodes the PNG on the
    way in, so the bytes differ, but the pixels must not.
    """
    pypdf = pytest.importorskip("pypdf", reason="pypdf is not installed here")
    image_module = pytest.importorskip("PIL.Image",
                                       reason="Pillow is not installed here")
    numpy = pytest.importorskip("numpy")

    # The figures the manuscript includes, READ OUT of paper.tex rather
    # than written out here: a hardcoded list is a second thing to
    # maintain, and a new \includegraphics would be counted but never
    # pixel-compared.
    sources = _manuscript_images()
    assert sources, "no \\includegraphics found in paper.tex -- the extractor broke"
    for source in sources:
        assert source.is_file(), f"{source} is missing; the PDF embeds it"

    reader = pypdf.PdfReader(str(_PDF))
    embedded = [image for page in reader.pages for image in page.images]
    assert len(embedded) == len(sources), (
        f"the PDF embeds {len(embedded)} image(s), the manuscript includes "
        f"{len(sources)}")

    for source in sources:
        with image_module.open(source) as handle:
            expected = numpy.asarray(handle.convert("RGB"))
        matched = False
        for image in embedded:
            with image_module.open(io.BytesIO(image.data)) as handle:
                actual = numpy.asarray(handle.convert("RGB"))
            if actual.shape != expected.shape:
                continue
            if numpy.array_equal(actual, expected):
                matched = True
                break
        assert matched, (
            f"the submission PDF does not embed {source.name} pixel for "
            f"pixel; the figure was regenerated without rebuilding the PDF")
