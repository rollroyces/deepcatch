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


def test_the_document_states_the_length_it_measured():
    """The number the document prints is the number it has."""
    document = _gen.build_document()
    length = _gen.abstract_length(_tex())
    assert f"{length} characters" in document, (
        f"the generated document does not state the measured length "
        f"({length} characters)")
    # And the abstract it carries is the abstract that was measured.
    assert _gen.latex_abstract(_tex()) in document


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
