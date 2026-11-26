"""The paste-ready bioRxiv abstract must be the manuscript's abstract.

What went wrong
---------------
`docs/PREPRINT_ABSTRACT.md` was a separately written document. It quoted
an older abstract (panel AUC 0.922 instead of 0.921, Sens@95% 0.600
instead of 0.770, the Fisher score 0.849 instead of 0.834), it carried
none of the external-validation or leave-one-study-out results, it had
an older title, and it named "Royce et al." on a single-author paper.
`docs/PATH_TO_IMPACT.md` tells the author to submit that document, so
the text that would have gone to bioRxiv was a different paper from the
PDF that would have gone with it.

The fix is a generator (`scripts/generate_preprint_abstract.py`) that
renders the manuscript's LaTeX. These tests cover the part that a
generated-document guard cannot: that the renderer is faithful, and
that it refuses rather than guesses.
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "generate_preprint_abstract",
        _REPO_ROOT / "scripts" / "generate_preprint_abstract.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_gen = _load_module()


def _tex() -> str:
    return (_REPO_ROOT / "paper" / "paper.tex").read_text(encoding="utf-8")


# The documents that tell the author where the preprint goes.
#
# The route was settled on 2026-11-25: bioRxiv rejected this manuscript
# because it requires an institutional affiliation, and arXiv asks for
# none. `docs/JOURNAL_SUBMISSION.md` then said, in its own summary,
# "Submit the preprint to bioRxiv now" while the rows below it recorded
# the rejection -- the summary and the detail had drifted apart, and
# nothing checked it. A document that sends the author to a server that
# has already refused the paper is worse than a stale number: it costs
# him a submission to find out.
_ROUTE_DOCUMENTS = (
    "NEXT_STEPS.md",
    "ZENODO_SETUP.md",
    "docs/JOURNAL_SUBMISSION.md",
    "docs/PATH_TO_IMPACT.md",
    "docs/PORTFOLIO_ROADMAP.md",
    "docs/SCHOLAR_VISIBILITY.md",
)

# An instruction to send the manuscript to bioRxiv. The closed attempt is
# still recorded in these documents -- it is *why* the route is arXiv --
# and those records are not imperative, so the check is for the
# instruction, not for the mention.
_TO_BIORXIV = re.compile(
    r"\b(?:submit|submitting|submitted|upload|uploaded|uploading|post|"
    r"posted|posting|publish|published|publishing)"
    r"\s+(?:\w+\s+){0,3}to\s+bioRxiv\b",
    re.I,
)


def test_the_abstract_is_the_manuscripts_abstract():
    """Every headline number in the manuscript is in the rendered text."""
    abstract = _gen.latex_abstract(_tex())
    for number in ("0.921", "0.770", "0.460", "0.834", "0.831", "0.9884",
                   "0.902", "0.9444", "0.708", "0.7385", "0.6805"):
        assert number in abstract, (
            f"the rendered abstract does not carry {number}; it is not the "
            f"manuscript's abstract")


def test_the_rendered_abstract_carries_no_latex():
    """The text that gets pasted must not contain markup."""
    abstract = _gen.latex_abstract(_tex())
    for leftover in ("\\", "$", "**", "\\%", "{", "}"):
        assert leftover not in abstract, (
            f"the rendered abstract still carries {leftover!r}; bioRxiv's "
            f"abstract field is plain text")


def test_the_rendered_abstract_fits_biorxivs_limit():
    """The pasted text must fit, with the margin the repo insists on."""
    length = _gen.abstract_length(_tex())
    assert length <= 2950, (
        f"the abstract is {length} characters; bioRxiv rejects anything over "
        f"3000 and the repo keeps a 2950 margin")


def test_the_documents_state_the_measured_length_and_the_live_route():
    """Every character count and every route in a document must be real.

    The character count comes first, and it has been wrong before: the
    arXiv field's length is also quoted in prose, in documents that
    are NOT generated -- and prose is not regenerated, so the number
    drifts away from the file the paste-ready text lives in. It already
    did: the field's text grew by two characters when its paragraph
    breaks were indented the way arXiv's own text processing needs, and
    three documents kept quoting the old length while the paste-ready
    file measured the new one. Every character count a document states
    must therefore be one the generator measures, or one of the two
    limits it holds.

    The route belongs in the same place, because it drifts the same way
    and costs more than a stale number. The route was settled on
    2026-11-25: bioRxiv rejected this manuscript because it requires an
    institutional affiliation, and arXiv asks for none.
    `docs/JOURNAL_SUBMISSION.md` then said, in its own summary, "Submit
    the preprint to bioRxiv now" while the rows below it recorded the
    rejection -- the summary and the detail had drifted apart, and
    nothing checked it. A document that sends the author to a server
    that has already refused the paper is worse than a stale number: it
    costs him a submission to find out.
    """
    document = _gen.build_document()
    length = _gen.abstract_length(_tex())
    assert f"{length} characters" in document, (
        f"the generated document does not state the measured length "
        f"({length} characters)")
    # And the abstract it carries is the abstract that was measured.
    assert _gen.latex_abstract(_tex()) in document

    tex = _tex()
    allowed = {
        _gen.abstract_length(tex),              # the manuscript's abstract
        len(_gen.check_arxiv_abstract(tex)),    # arXiv's metadata field
        _gen.ARXIV_LIMIT,
        _gen.BIORXIV_LIMIT,
    }
    offenders: list[str] = []
    for path in sorted(_REPO_ROOT.rglob("*.md")):
        if not path.is_file() or ".git" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        # A document that states one of these lengths is one that talks
        # about a preprint server's field. The one document with a
        # different character count (a data-access form's field limit)
        # does not mention a server at all.
        if "arXiv" not in text and "bioRxiv" not in text:
            continue
        for value in re.findall(r"\b(\d{1,3}(?:,\d{3})|\d{4})\s+characters\b",
                                 text):
            if int(value.replace(",", "")) not in allowed:
                offenders.append(
                    f"{path.relative_to(_REPO_ROOT)} states {value} characters")
    assert not offenders, (
        "a document quotes a character count that the generator does not "
        "measure:\n  " + "\n  ".join(offenders)
        + "\nThe paste-ready text is generated and its length is measured, "
          "not written by hand: re-run the generator and quote the number "
          "it prints.")

    # The route, in the documents that tell the author where to take it.
    # The closed attempt is still recorded in them -- it is *why* the
    # route is arXiv -- and those records are not imperative, so the
    # check is for the instruction, not for the mention.
    sent_to_biorxiv: list[str] = []
    for rel in _ROUTE_DOCUMENTS:
        path = _REPO_ROOT / rel
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for match in _TO_BIORXIV.finditer(text):
            sent_to_biorxiv.append(f"{rel}: {match.group(0)!r}")
        assert "arXiv" in text, (
            f"{rel} does not name the route at all; the preprint goes to "
            "arXiv, and a document that never says where it goes is one "
            "more document to drift.")
    assert not sent_to_biorxiv, (
        "a document still sends the preprint to bioRxiv:\n  "
        + "\n  ".join(sent_to_biorxiv)
        + "\nbioRxiv rejected this manuscript on 2026-11-25 under its "
          "affiliation rule (\"your manuscript is inappropriate for bioRxiv "
          "because bioRxiv requires authors to add an institutional "
          "affiliation\"), so the route is arXiv.")

    # And the route has to be named where the reader acts on it: the
    # submission sequence and the recommendation in the venue map.
    journal = (_REPO_ROOT / "docs" / "JOURNAL_SUBMISSION.md").read_text(
        encoding="utf-8")
    assert "arXiv preprint" in journal, (
        "docs/JOURNAL_SUBMISSION.md's submission sequence does not name the "
        "route; the preprint goes to arXiv.")
    parts = journal.split("**Recommendation.**", 1)
    assert len(parts) == 2, (
        "docs/JOURNAL_SUBMISSION.md has no Recommendation paragraph; the "
        "route has to be stated somewhere a reader acts on.")
    assert "arXiv" in parts[1].split("\n\n", 1)[0], (
        "docs/JOURNAL_SUBMISSION.md's recommendation does not name the "
        "route; the preprint goes to arXiv.")


def test_the_title_and_the_citation_file_agree():
    """Two copies of the paper's title must not drift apart."""
    title = _gen.latex_title(_tex())
    cff = (_REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    stated = _gen.citation_field(cff, "title")
    assert _gen.plain(stated) == title


def test_an_unknown_command_is_refused_rather_than_dropped():
    """A construct the renderer does not know must fail the build.

    An abstract that silently loses a symbol is worse than one that
    stops the build, so the converter's vocabulary is closed.
    """
    with pytest.raises(_gen.UnsupportedConstruct):
        _gen.render_text(r"a \unknowncommand{x} b")
    with pytest.raises(_gen.UnsupportedConstruct):
        _gen.render_math(r"\unknownsymbol")
    with pytest.raises(_gen.UnsupportedConstruct):
        _gen.render_text(r"unbalanced $ sign")


def test_the_math_is_rendered_as_its_characters():
    """The formulae the abstract uses come out as readable symbols."""
    assert _gen.render_math(r"\sim 10\times") == "~ 10×"
    assert _gen.render_math(r"\le 10^{-4}") == "≤ 10⁻⁴"
    assert _gen.render_math(r"\ge 50{,}000\times") == "≥ 50,000×"


def test_the_check_mode_reports_a_stale_document(tmp_path):
    """`--check` fails when the committed document is not what we produce."""
    import subprocess
    import sys

    stale = tmp_path / "PREPRINT_ABSTRACT.md"
    stale.write_text("stale\n", encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, "scripts/generate_preprint_abstract.py",
         "--check", "--out", str(stale)],
        cwd=str(_REPO_ROOT), capture_output=True, text=True, timeout=120)
    assert proc.returncode == 1, (
        f"--check accepted a stale document:\n{proc.stdout}{proc.stderr}")


def test_the_check_mode_accepts_the_committed_document():
    """The committed document is what the generator produces."""
    import subprocess
    import sys

    proc = subprocess.run(
        [sys.executable, "scripts/generate_preprint_abstract.py", "--check"],
        cwd=str(_REPO_ROOT), capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, (
        f"the committed document has drifted from paper/paper.tex:\n"
        f"{proc.stdout}{proc.stderr}")
