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

And a third, which belongs to the same question -- what the repository's
manuscripts say about the work: every manuscript that is offered for
submission declares its AI assistance, names the tools, and says that the
author reviewed what they produced. arXiv's moderation policy asks for the
use of generative AI to be reported *in the work itself*, so a disclosure
that lives in a file beside the manuscript does not satisfy it, and one
that names no tool cannot be checked.

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

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# The sections JOSS asks a submission to carry, in the order it lists them.
#
# Four of them are new: JOSS changed its submission format in late 2025,
# and its own review bot now reports "The paper includes a State of the
# field section", "a Software design section", "a Research impact
# statement section" and "an AI usage disclosure section" as the checks
# it runs on the paper file. A paper in the old format is rejected in
# pre-review, so the list here is the format, not a style preference.
_JOSS_SECTIONS = (
    "Summary",
    "Statement of need",
    "State of the field",
    "Software design",
    "Features",
    "Quality and reproducibility",
    "Documentation",
    "Code availability",
    "Research impact statement",
    "Limitations",
    "AI usage disclosure",
    "Acknowledgements",
)

# The AI usage disclosure has to name the tools and confirm the human
# review -- JOSS's policy asks for both, and treats an incomplete
# disclosure as a research-integrity matter rather than a formatting
# one. A heading with a sentence under it does not satisfy either.
_AI_PROVIDERS = ("OpenAI", "Anthropic", "Google", "MiniMax")

# The manuscripts that carry an author block and are offered for
# submission. arXiv's moderation policy asks for the use of generative
# AI to be reported *in the work itself*, and the disclosure has to be
# in the manuscript that is uploaded -- not in a file beside it.
#
# Three of these carried the disclosure and `paper/PAPER.md` did not,
# and nothing checked any of them: an edit could have dropped it from
# the manuscript being uploaded and the submission would have gone out
# without it.
_AI_DISCLOSURE_SURFACES = (
    "paper/paper.tex",
    "paper/METHODS_PAPER.md",
    "paper/PAPER.md",
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


def _joss_word_count() -> int:
    """The JOSS paper's length, measured by the script that owns the rule.

    `scripts/count_joss_paper_words.py` counts the prose of every section
    -- the front matter and the markdown markup excluded -- and reports the
    per-section numbers as well. The venue map's figure is checked against
    it, so the two cannot disagree.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "count_joss_paper_words",
        _REPO_ROOT / "scripts" / "count_joss_paper_words.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.section_words(
        (_REPO_ROOT / "paper.md").read_text(encoding="utf-8"))["(prose)"]


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


def _paragraphs(text: str) -> list[str]:
    """The paragraphs of a text, whatever markup it is written in."""
    return [block.strip() for block in re.split(r"\n\s*\n", text) if block.strip()]


def _ai_disclosure(text: str) -> str:
    """The paragraphs that disclose the AI assistance."""
    return "\n".join(
        block for block in _paragraphs(text)
        if "generative ai" in block.lower())


def test_the_manuscripts_declare_the_sections_and_the_assistance_their_venues_ask_for():
    """What each manuscript must say about itself, it must still say."""
    text = _joss_manuscript()
    headings = {
        match.group(1).strip()
        for match in re.finditer(r"^#\s+(.*)$", text, re.M)
    }
    missing = [section for section in _JOSS_SECTIONS if section not in headings]
    assert not missing, (
        "the JOSS manuscript is missing the sections JOSS asks for: "
        + ", ".join(missing)
        + "\n`docs/JOURNAL_SUBMISSION.md` states that the manuscript is in "
          "JOSS format; that claim is only true while the sections are "
          "there."
    )

    # The JOSS paper's length is the number a reviewer counts first, and
    # it was hand-counted: the venue map said 2,509 words and 759 over,
    # on a paper that measures 2,562 and 812 over. The row is held to the
    # measurement now, so the figure in it cannot drift away from the
    # paper again.
    measured = _joss_word_count()
    venue_map = (_REPO_ROOT / "docs" / "JOURNAL_SUBMISSION.md").read_text(
        encoding="utf-8")
    # The count is written with a thousands separator ("2,562 words"), so
    # the pattern has to allow it and the digits have to be read back
    # without it.
    stated = re.search(r"the paper's prose is \*\*([\d,]+) words\*\*", venue_map)
    assert stated, (
        "the venue map no longer states the JOSS paper's length; put the "
        "measured number back (scripts/count_joss_paper_words.py)")
    assert int(stated.group(1).replace(",", "")) == measured, (
        f"the venue map says the JOSS paper's prose is {stated.group(1)} words; "
        f"paper.md measures {measured}")

    # The AI usage disclosure is the one section whose CONTENT JOSS
    # checks, not just its presence: the policy asks for the tools used,
    # the nature and scope of the assistance, and the author's
    # confirmation of review. A heading with a placeholder under it is
    # what the policy calls an incomplete disclosure.
    body = _section_body(text, "AI usage disclosure")
    assert len(body.split()) >= 40, (
        "the AI usage disclosure is a placeholder: JOSS asks for the "
        "tools used, the nature and scope of the assistance, and the "
        "author's confirmation of review, and treats an incomplete "
        "disclosure as a research-integrity matter. The section has "
        f"{len(body.split())} words."
    )
    assert any(provider in body for provider in _AI_PROVIDERS), (
        "the AI usage disclosure does not name any of "
        + ", ".join(_AI_PROVIDERS)
        + ". JOSS's AI usage policy asks for the tools and models used, "
          "and a disclosure that names none of them cannot be checked."
    )
    assert "review" in body.lower(), (
        "the AI usage disclosure does not confirm that the author "
        "reviewed and validated the AI-assisted output, which JOSS's "
        "policy requires."
    )

    # The other three manuscripts have no such section to check, so the
    # disclosure has to be somewhere in the manuscript itself -- and it
    # has to say the two things a disclosure is: which tools were used,
    # and that the author checked what they produced. arXiv's moderation
    # policy asks for the use of text-to-text generative AI to be
    # reported in the work itself, and the work is this file.
    for path in _AI_DISCLOSURE_SURFACES:
        body = _ai_disclosure(
            (_REPO_ROOT / path).read_text(encoding="utf-8"))
        assert body, (
            f"{path} does not disclose its AI assistance. arXiv's "
            "moderation policy asks for the use of generative AI to be "
            "reported in the work itself, and the work is the file that "
            "is uploaded."
        )
        assert any(provider in body for provider in _AI_PROVIDERS), (
            f"{path} discloses AI assistance without naming any of "
            + ", ".join(_AI_PROVIDERS)
            + ". A disclosure that names no tool cannot be checked."
        )
        assert "reviewed" in body.lower(), (
            f"{path} does not confirm that the author reviewed and edited "
            "the AI-assisted output, which is the part of the disclosure "
            "that makes it a disclosure rather than a mention."
        )


def _section_body(text: str, heading: str) -> str:
    """The body of a markdown section, up to the next heading."""
    match = re.search(
        r"^#\s+" + re.escape(heading) + r"\s*$(.*?)(?=^#\s|\Z)", text,
        re.M | re.S)
    return match.group(1).strip() if match else ""


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


# ─────────────────────────────────────────────────────────────────────
# The affiliation STRING, not only whether an employer is named.
#
# The check above asks "does one surface name an employer and the others
# not". It does not ask whether the surfaces declare the same
# affiliation, and they did not: the manuscript's author block says
# "Independent Researcher", `CITATION.cff` says "Independent Researcher
# (no institutional affiliation; verified email domain umich.edu)", the
# generated abstract says "Independent Researcher", and the submission
# sheet's form values said "Independent Researcher, Hong Kong" -- a
# locality that appears nowhere in the paper the form uploads.
#
# That matters because the form's Author field and the PDF's author block
# are read side by side: a reviewer sees two different affiliations for
# the same author, and `docs/SCHOLAR_VISIBILITY.md` §1 says Google
# Scholar matches a paper to a profile by the affiliation string. A
# substring difference is enough to keep the work off the profile.
#
# The comparison ignores a trailing parenthetical, because the citation
# metadata qualifies its value that way ("Independent Researcher (no
# institutional affiliation; ...)") and that qualifier is metadata about
# the declaration rather than a different declaration. It does NOT
# ignore a comma-suffix: "Independent Researcher, Hong Kong" is a
# different string from "Independent Researcher", and a locality on one
# surface only is exactly the drift this is meant to catch.
#
# The surfaces compared are the ones that DECLARE an affiliation as a
# field of their own: the manuscript's author block (so the title page of
# the uploaded PDF), the citation metadata (which is what Zenodo
# archives), the generated abstract that goes into the form, the JOSS
# paper's author field, the model card, and the submission sheet's two
# form rows. Deliberately out of scope: the roster line in `TEAM.md` and
# the docs-site bylines. Those are prose about the author rather than a
# declaration -- a byline is contact information, and it may carry a
# locality that the declared affiliation does not.
# ─────────────────────────────────────────────────────────────────────

_SUBMISSION_SHEET = "paper/BIORXIV_SUBMISSION.md"


def _core(declaration: str) -> str:
    """A declaration without its trailing parenthetical qualifier."""
    text = re.sub(r"\s*\([^()]*\)\s*$", "", declaration.strip())
    text = re.sub(r"\\[a-zA-Z]+", " ", text)
    text = text.replace("\\", " ").replace("{", " ").replace("}", " ")
    return re.sub(r"\s+", " ", text).strip().lower()


def _latex_affiliation() -> str:
    """The affiliation line of the LaTeX manuscript's author block."""
    lines = [line.strip() for line in _latex_author_block().splitlines()]
    kept = [
        line for line in lines
        if line
        and "orcid.org" not in line
        and "url{" not in line
        and "ORCID" not in line
    ]
    return " ".join(kept[1:])


def _citation_affiliation() -> str:
    """The `affiliation:` value of the citation metadata."""
    match = re.search(r'^\s*affiliation:\s*"?([^"\n]+)"?', _citation_author_block(),
                      re.M)
    return match.group(1).strip() if match else ""


def _abstract_affiliation() -> str:
    """The affiliation on the generated abstract's author line."""
    text = (_REPO_ROOT / "docs" / "PREPRINT_ABSTRACT.md").read_text(
        encoding="utf-8")
    match = re.search(r"\*\*Authors:\*\*\s*(.+)", text)
    if not match:
        return ""
    line = match.group(1)
    if "—" in line:
        line = line.split("—", 1)[1]
    line = re.split(r"\(ORCID", line)[0]
    return line.strip().rstrip("·").strip()


def _joss_affiliation() -> str:
    """The affiliation in the JOSS paper's author field."""
    match = re.search(r"name:\s*(.+)", _front_matter(_joss_manuscript()))
    if not match:
        return ""
    return match.group(1).split("(")[-1].strip().rstrip(")").strip()


def _model_card_affiliation() -> str:
    """The affiliation on the model card's author line."""
    text = (_REPO_ROOT / "MODEL.md").read_text(encoding="utf-8")
    match = re.search(r"\*\*Authors\*\*\s*\|\s*(.+)", text)
    if not match:
        return ""
    line = match.group(1)
    if "(" in line and ")" in line:
        line = line.split("(", 1)[1].rsplit(")", 1)[0]
    # The ORCID rides in the same parenthetical; it has its own check.
    return line.split(",")[0].strip()


def _submission_sheet_affiliations() -> list[tuple[str, str]]:
    """The affiliation the submission sheet tells you to paste.

    Only the form's own field rows count. The sheet also shows the
    affiliation inside code blocks (the proposed edits for the two
    endings in §7), and those are illustrations of the decision, not
    declarations of it.

    The two rows are read differently on purpose: the Author row carries
    the name with the affiliation in parentheses, while the Author
    affiliation row IS the affiliation, optionally with a trailing
    qualifier in parentheses.
    """
    path = _REPO_ROOT / _SUBMISSION_SHEET
    if not path.is_file():
        return []
    text = path.read_text(encoding="utf-8")
    found: list[tuple[str, str]] = []
    for label in ("Author", "Author affiliation"):
        match = re.search(r"\| " + label + r" \| (.+?) \|", text)
        if not match:
            continue
        value = match.group(1).strip()
        if label == "Author" and "(" in value and ")" in value:
            value = value.split("(", 1)[1].rsplit(")", 1)[0]
        found.append((label, value.strip()))
    return found


def test_the_declared_affiliation_is_the_same_string_everywhere():
    """Every surface that declares the affiliation must declare one string."""
    reference = _core(_latex_affiliation())
    assert reference, (
        "no affiliation could be read from paper/paper.tex's author block")

    declarations = {
        "paper/paper.tex": _latex_affiliation(),
        "CITATION.cff": _citation_affiliation(),
        "docs/PREPRINT_ABSTRACT.md": _abstract_affiliation(),
        "paper.md": _joss_affiliation(),
        "MODEL.md": _model_card_affiliation(),
    }
    for name, value in declarations.items():
        assert value, f"{name} declares no affiliation to compare"

    found = {name: _core(value) for name, value in declarations.items()}
    distinct = set(found.values())
    assert len(distinct) == 1, (
        "the surfaces declare different affiliations:\n"
        + "\n".join(f"  {name}: {value!r}" for name, value in found.items())
        + "\nThe manuscript's author block, the citation metadata, the "
          "abstract and the JOSS paper are read together -- and the "
          "preprint form's affiliation has to be the same string as the "
          "paper it uploads. `paper/BIORXIV_SUBMISSION.md` §7 records the "
          "one affiliation decision that is still open; when it is settled, "
          "change every surface on the same day."
    )

    # The two markdown manuscripts declare the affiliation in prose: a
    # parenthetical on the project line, or a phrase before the project
    # name. The declared string is searched for rather than parsed out of
    # the prose -- a manuscript that names a different employer, or a
    # truncated one, does not contain it.
    #
    # These two were unguarded until 2026-11-24: they are read by the
    # same reviewers as the PDF and the citation metadata, and the
    # affiliation decision moves them too, so a partial edit could have
    # left one of them naming a different institution than the paper.
    for path in ("paper/METHODS_PAPER.md", "paper/PAPER.md"):
        block = _markdown_author_block(path)
        assert block.strip(), f"{path} has no author block to check"
        assert reference.lower() in block.lower(), (
            f"{path} does not declare the affiliation the manuscript "
            f"declares:\n  manuscript: {reference!r}\n  {path}: "
            f"{block.strip()!r}\nThese two manuscripts are read side by "
            "side with the PDF and the citation metadata, and the "
            "affiliation they carry has to be the same one."
        )


def test_the_submission_sheet_declares_the_manuscripts_affiliation():
    """The form's value must not add what the paper does not declare.

    The sheet is gitignored, so it is absent from a fresh clone and from
    CI; the check says so instead of failing there.
    """
    sheet = _submission_sheet_affiliations()
    if not sheet:
        pytest.skip(f"{_SUBMISSION_SHEET} is not present")
    reference = _core(_latex_affiliation())
    problems = [
        f"{_SUBMISSION_SHEET}'s '{label}' row declares {value!r}, the "
        f"manuscript declares {reference!r}"
        for label, value in sheet if _core(value) != reference
    ]
    assert not problems, (
        "the submission sheet's form values disagree with the manuscript's "
        "author block:\n  " + "\n  ".join(problems)
        + "\nThe PDF that goes up and the form that describes it are read "
          "side by side, and Google Scholar matches a paper to a profile "
          "by the affiliation string."
    )


def test_the_core_reader_tolerates_a_qualifier_but_not_a_locality():
    """The comparison must not hide a locality on one surface only."""
    assert _core("Independent Researcher") == "independent researcher"
    assert _core(
        "Independent Researcher (no institutional affiliation; verified "
        "email domain umich.edu)") == "independent researcher"
    assert _core("\\small Independent Researcher") == "independent researcher"
    # A locality is a different string, and must not be swallowed.
    assert _core("Independent Researcher, Hong Kong") != "independent researcher"
    # Nor may a different affiliation hide behind a parenthetical.
    assert _core("University of Michigan (Independent Researcher)") != (
        "independent researcher")

