#!/usr/bin/env python3
"""Render `docs/PREPRINT_ABSTRACT.md` from the manuscript's LaTeX source.

Why this file exists
--------------------
`docs/PREPRINT_ABSTRACT.md` is the text the author pastes into the bioRxiv
submission form -- `docs/PATH_TO_IMPACT.md` says so explicitly ("Submit
abstract from docs/PREPRINT_ABSTRACT.md"). For a long time it was a
*separately written* document: it carried an older abstract (panel AUC 0.922
instead of 0.921, Sens@95% 0.600 instead of 0.770, Fisher 0.849 instead of
0.834), none of the external-validation or leave-one-study-out results, an
older title, and the author line "Royce et al." on a single-author paper.
Whoever submitted it would have pasted a different preprint than the PDF
they uploaded.

The fix is the same one the rest of this repository uses: one source of
truth, and a generator that reads it. The abstract in `paper/paper.tex` is
the text the compiled PDF carries and the text the reproducibility test
measures against bioRxiv's 3000-character limit, so it is the one thing
that should also produce the paste-ready document.

The LaTeX -> text conversion uses a CLOSED vocabulary. A construct that is
not in the table raises, rather than being dropped: an abstract that
silently loses a symbol is worse than one that fails the build.

Two documents are generated, because two preprint servers take two
different texts. bioRxiv takes the manuscript's abstract as it is
(2,940 characters, under its 3,000-character limit). arXiv takes a
metadata field, and arXiv's own help says of that field: "Keep it
short -- abstracts longer than 1920 characters will not be accepted; and
abridge your abstract if necessary"
(https://info.arxiv.org/help/prep.html, "Title and Abstract Fields").
The limit is on the FIELD -- the abstract page, the announcement and
the RSS feed -- and not on the paper: arXiv does not truncate the
manuscript. So the manuscript keeps its full abstract, and the form
field gets an abridged one.

The abridged text is written here rather than derived, because an
abridgement is an editorial judgement. What is NOT left to judgement is
that it must not invent a number: every figure it states has to appear
in the manuscript's abstract, and the build refuses if one does not.

Usage
-----
    python3 scripts/generate_preprint_abstract.py                      # write the doc
    python3 scripts/generate_preprint_abstract.py --print               # abstract only
    python3 scripts/generate_preprint_abstract.py --check               # CI mode
    python3 scripts/generate_preprint_abstract.py --variant arxiv       # arXiv's field

`--check` exits 1 when the committed document is not what the generator
produces, without writing anything.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, "paper", "paper.tex")
CITATION = os.path.join(ROOT, "CITATION.cff")
OUT = os.path.join(ROOT, "docs", "PREPRINT_ABSTRACT.md")
ARXIV_OUT = os.path.join(ROOT, "docs", "ARXIV_ABSTRACT.md")

# bioRxiv rejects an abstract over this many characters.
BIORXIV_LIMIT = 3000

# arXiv rejects a metadata abstract over this many characters. The limit
# is on the form field, not on the paper: arXiv's metadata help puts it
# under "Title and Abstract Fields", the text used on the abstract page,
# in the announcement and in the RSS feed. The manuscript keeps its own
# abstract.
ARXIV_LIMIT = 1920

# The abstract for arXiv's metadata field: the manuscript's own abstract
# with the secondary numbers dropped, because the field is 1,080
# characters shorter than the manuscript's text.
#
# It is written rather than derived -- an abridgement is an editorial
# judgement -- but every number in it is checked against the manuscript
# below, so it cannot state a figure the paper does not. The positives
# and the negatives are both kept: the held-out transfer with its
# interval, the caveat that the pooled figure is not a generalisation
# claim, the cost of holding a study out of training, and the note that
# the fused channel's mutation arm is a synthetic proxy. The keywords
# go to the form's Comments field instead, because arXiv's metadata has
# no keyword field of its own.
ARXIV_ABSTRACT = """\
Background. Molecular residual disease (MRD) detection from cell-free DNA (cfDNA) requires tracking a patient's tumor mutations in plasma at circulating tumor DNA (ctDNA) fractions as low as 0.1%. Commercial assays such as Signatera and CAPP-Seq demonstrate this clinically, but their validation data are not open.

Results. We present DeepCatch, an open-source panel-based MRD detection pipeline benchmarked against real tumor mutations. Using 20 real TCGA-LUAD patients (5,738 somatic mutations with real read counts, GDC open access) as ground truth and a context-aware error model, we simulate plasma cfDNA at five tumor fractions (10%-0.1%) and score three methods at fixed specificity, without test-set threshold optimisation. At 0.1% ctDNA, panel-level log-likelihood-ratio aggregation achieves AUC 0.921 (Sens@95% = 0.770), and an assay sweep shows duplex-UMI consensus (error <= 10^-4) or depth >= 50,000x each reach Sens@95% = 1.000. A tumor-naive fragmentomics channel (627 cross-study FinaleDB samples) averaged with a mutation channel reaches AUC 0.9884, but its mutation channel is a synthetic proxy: a complementarity result, not a performance claim. On never-trained publications, transfer reaches AUC 0.9444 (95% CI 0.8731-1.0000) on one held-out whole-genome study and 0.708 on a second, pooling to 0.7385; transfer is real but not uniform, an exome-targeted study only partly transfers, and holding a study out of training costs most of the pooled single-type AUC (0.6805). The pooled figure is a within-cohort separation number, not a generalisation claim.

Conclusions. DeepCatch provides a reproducible, honest benchmark for ultra-sensitive MRD detection on open-access data and turns the assay sweep into production specifications: duplex-UMI error suppression and ~50,000x depth suffice for AUC 1.000 at 0.1% ctDNA. Clinical validation on real plasma cfDNA remains the essential next step."""


# The math symbols the abstract may use. Anything else raises.
MATH_SYMBOLS = {
    "sim": "~",
    "times": "×",
    "le": "≤",
    "leq": "≤",
    "ge": "≥",
    "geq": "≥",
    "cdot": "·",
}

SUPERSCRIPT = {
    "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
    "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
    "-": "⁻", "+": "⁺",
}

SUBSCRIPT = {
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄",
    "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉",
}


class UnsupportedConstruct(ValueError):
    """The abstract uses LaTeX this converter does not know how to render."""


def read_group(text: str, start: int) -> tuple[str, int]:
    """Read a ``{...}`` group beginning at ``start``; return it and the index after it."""
    if start >= len(text) or text[start] != "{":
        raise UnsupportedConstruct(f"expected '{{' at offset {start}")
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:i], i + 1
    raise UnsupportedConstruct(f"unbalanced '{{' at offset {start}")


def render_math(source: str) -> str:
    """Render one math group (the text between ``$`` delimiters)."""
    source = source.replace("{,}", ",")
    out: list[str] = []
    i = 0
    while i < len(source):
        ch = source[i]
        if ch == "\\":
            match = re.match(r"\\([a-zA-Z]+)", source[i:])
            if not match:
                raise UnsupportedConstruct(
                    f"unsupported math escape at {source[i:i + 12]!r}")
            name = match.group(1)
            if name not in MATH_SYMBOLS:
                raise UnsupportedConstruct(f"unsupported math command \\{name}")
            out.append(MATH_SYMBOLS[name])
            i += 1 + len(name)
        elif ch == "^":
            body, i = read_group(source, i + 1)
            out.append(_script(body, SUPERSCRIPT, "superscript"))
        elif ch == "_":
            body, i = read_group(source, i + 1)
            out.append(_script(body, SUBSCRIPT, "subscript"))
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def _script(body: str, table: dict, what: str) -> str:
    try:
        return "".join(table[ch] for ch in body)
    except KeyError as exc:
        raise UnsupportedConstruct(
            f"unsupported {what} character {exc.args[0]!r}") from None


def render_text(source: str) -> str:
    """Render LaTeX body text."""
    out: list[str] = []
    i = 0
    while i < len(source):
        ch = source[i]
        if ch == "\\":
            # `\%` is an escaped percent sign, not a command.
            if source.startswith("\\%", i):
                out.append("%")
                i += 2
                continue
            match = re.match(r"\\([a-zA-Z]+)", source[i:])
            if not match:
                raise UnsupportedConstruct(
                    f"unsupported escape at {source[i:i + 12]!r}")
            name = match.group(1)
            if name in ("textbf", "emph"):
                inner, i = read_group(source, i + 1 + len(name))
                body = render_text(inner)
                out.append(f"**{body}**" if name == "textbf" else f"*{body}*")
            else:
                raise UnsupportedConstruct(f"unsupported command \\{name}")
        elif ch == "$":
            end = source.find("$", i + 1)
            if end == -1:
                raise UnsupportedConstruct("unbalanced '$'")
            out.append(render_math(source[i + 1:end]))
            i = end + 1
        elif source.startswith("--", i):
            out.append("–")
            i += 2
        elif source.startswith("{,}", i):
            out.append(",")
            i += 3
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def plain(text: str) -> str:
    """Strip the Markdown emphasis markers the renderer adds."""
    return text.replace("**", "").replace("*", "").strip()


def latex_title(tex: str) -> str:
    match = re.search(r"\\title\{(.*?)\}\n", tex, re.S)
    if not match:
        raise ValueError("paper.tex has no \\title")
    return plain(render_text(match.group(1)))


def latex_author(tex: str) -> tuple[str, str, str]:
    """The author name, the affiliation line under it, and the ORCID.

    The author block is LaTeX, and the ORCID is a link in it
    (`\\href{...}{...}`), so it has to be rendered rather than copied:
    a pasted abstract carrying `\\href{https://orcid.org/...}` is markup
    in a field that wants plain text. The ORCID is also kept out of the
    affiliation -- it identifies the author, it is not an employer.
    """
    match = re.search(r"\\author\{(.*?)\}\n", tex, re.S)
    if not match:
        raise ValueError("paper.tex has no \\author")
    lines = [ln.strip() for ln in match.group(1).split("\\\\") if ln.strip()]
    name = plain(render_text(lines[0]))
    affiliation = ""
    orcid = ""
    for line in lines[1:]:
        text = re.sub(r"\\small\s*", "", line).strip()
        text = re.sub(r"\\url\{([^}]*)\}", r"\1", text)
        text = re.sub(r"\\href\{([^}]*)\}\{([^}]*)\}", r"\2", text)
        text = text.replace("---", "—")
        found = re.search(r"\b(\d{4}-\d{4}-\d{4}-\d{3}[\dX])\b", text)
        if found:
            orcid = found.group(1)
            continue
        if text and "github.com" not in text:
            affiliation = text
    return name, affiliation, orcid


def citation_orcid(cff: str) -> str:
    """The ORCID the citation metadata records."""
    match = re.search(
        r"orcid:\s*\"?https?://orcid\.org/(\d{4}-\d{4}-\d{4}-\d{3}[\dX])", cff)
    return match.group(1) if match else ""


def latex_abstract(tex: str) -> str:
    """The abstract, rendered as the plain text that gets pasted to bioRxiv.

    The emphasis markers are stripped on purpose: the bioRxiv abstract
    field is plain text, so a document carrying `**AUC 0.921**` would be
    pasted with literal asterisks, and the character count would be
    measured against a string nobody submits.
    """
    match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    if not match:
        raise ValueError("paper.tex has no abstract")
    paragraphs = [p.strip() for p in match.group(1).split("\n\n") if p.strip()]
    return "\n\n".join(plain(render_text(p)) for p in paragraphs)


def abstract_length(tex: str) -> int:
    """The length of the pasted abstract, counted the way bioRxiv counts.

    This is the RAW length: the paragraph breaks are characters too,
    because they are characters in the text that gets pasted into the
    form. A count that collapsed the whitespace would under-count by
    three characters and would no longer be the number bioRxiv sees.
    """
    return len(latex_abstract(tex))


def citation_field(cff: str, key: str) -> str:
    match = re.search(rf'^{key}:\s*"?([^"\n]+)"?\s*$', cff, re.M)
    if not match:
        raise ValueError(f"CITATION.cff has no {key}")
    return match.group(1).strip()


# The superscript digits the renderer produces from the manuscript's
# LaTeX. They are folded to ASCII before the two texts are compared,
# because the manuscript writes "10^{-4}" as "10⁻⁴" and the form field
# writes "10^-4", and a comparison that cannot see that they are the
# same number would report a missing one.
_SUPERSCRIPT_DIGITS = {
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
    "⁻": "-",
}


def numeric_tokens(text: str) -> set[str]:
    """The numbers a text states, in a form two texts can be compared in.

    Thousands separators are dropped and the superscripts are folded to
    ASCII, so "5,738" and "5738" are one number, and "10⁻⁴" and "10^-4"
    are the same pair.
    """
    for unicode_digit, ascii_digit in _SUPERSCRIPT_DIGITS.items():
        text = text.replace(unicode_digit, ascii_digit)
    text = text.replace("^", "")
    return {
        match.group(0).replace(",", "")
        for match in re.finditer(r"\d+(?:,\d{3})*(?:\.\d+)?", text)
    }


def arxiv_field_text() -> str:
    """The abstract as arXiv's metadata field wants it.

    arXiv's own metadata help: "Carriage returns will be stripped unless
    they are followed by leading white spaces. So if you want a new
    paragraph or a table of contents, be sure to indent the lines after
    the carriage return." A blank line is not an indented line, so the
    paragraph breaks would be stripped -- and the parts of the abstract
    would be welded together ("...are not open.Results. We present...").
    The parts are therefore separated by a newline and two spaces, which
    is the pattern arXiv documents for exactly this purpose.
    """
    paragraphs = [p.strip() for p in ARXIV_ABSTRACT.strip().split("\n\n")
                  if p.strip()]
    return "\n  ".join(paragraphs)


def check_arxiv_abstract(tex: str) -> str:
    """The abridged abstract, after checking it against the manuscript.

    Two things can go wrong, and a word count catches neither of them.
    The text can be too long for the field, or it can state a number the
    manuscript does not -- and an abridgement that invents a figure is
    worse than one that is too long, because the field is what a reader
    sees first and the paper is what they check it against.
    """
    abstract = arxiv_field_text()
    length = len(abstract)
    if length > ARXIV_LIMIT:
        raise ValueError(
            f"the abridged abstract is {length} characters; arXiv's metadata "
            f"field rejects anything over {ARXIV_LIMIT}")
    missing = numeric_tokens(abstract) - numeric_tokens(latex_abstract(tex))
    if missing:
        raise ValueError(
            "the abridged abstract states "
            + ", ".join(sorted(missing))
            + ", which the manuscript's abstract does not contain")
    return abstract


def build_arxiv_document() -> str:
    """The document holding the text for arXiv's metadata field."""
    tex = open(TEX, encoding="utf-8").read()
    cff = open(CITATION, encoding="utf-8").read()

    title = latex_title(tex)
    name, affiliation, orcid = latex_author(tex)
    abstract = check_arxiv_abstract(tex)
    length = len(abstract)
    version = citation_field(cff, "version")
    repo = citation_field(cff, "repository-code")

    author_line = f"**Authors:** {name} — {affiliation}"
    if orcid:
        author_line += f" (ORCID {orcid})"

    return (
        "<!-- GENERATED FILE - do not hand-edit.\n"
        "     Source: paper/paper.tex (the manuscript the PDF is compiled from).\n"
        "     Regenerate: python3 scripts/generate_preprint_abstract.py "
        "--variant arxiv\n"
        "     Guarded by: test/test_generated_docs_match_their_generators.py\n"
        "     -->\n"
        "\n"
        f"# {title}\n"
        "\n"
        f"{author_line}\n"
        "\n"
        f"**Repository:** {repo} (v{version})\n"
        "\n"
        "## Abstract for arXiv's metadata field\n"
        "\n"
        "Paste this into the form's Abstract field, exactly as it stands:\n"
        "the two spaces after each line break are what keep the paragraph\n"
        "breaks from being stripped by arXiv's own text processing.\n"
        "\n"
        "```\n"
        f"{abstract}\n"
        "```\n"
        "\n"
        f"<!-- {length} characters against arXiv's {ARXIV_LIMIT}-character "
        "limit on the metadata field.\n"
        "     The paper keeps the manuscript's own abstract: arXiv's limit is "
        "on the form field,\n"
        "     not on the manuscript. -->\n"
    )


def build_document() -> str:
    tex = open(TEX, encoding="utf-8").read()
    cff = open(CITATION, encoding="utf-8").read()

    title = latex_title(tex)
    name, affiliation, orcid = latex_author(tex)
    abstract = latex_abstract(tex)
    length = abstract_length(tex)
    version = citation_field(cff, "version")
    repo = citation_field(cff, "repository-code")

    # The metadata in CITATION.cff and the manuscript's own title are two
    # copies of the same fact, and a citation file that names a different
    # paper than the PDF is the kind of thing a reviewer checks first.
    cff_title = citation_field(cff, "title")
    if plain(cff_title) != title:
        raise ValueError(
            f"CITATION.cff titles the paper {cff_title!r}, but paper.tex "
            f"titles it {title!r}; they must agree before the abstract is "
            f"published")

    # The ORCID is the same fact in two places: the author block a reader
    # sees on the first page, and the citation metadata a tool reads.
    cff_orcid = citation_orcid(cff)
    if orcid and cff_orcid and orcid != cff_orcid:
        raise ValueError(
            f"paper.tex gives the author ORCID {orcid}, but CITATION.cff "
            f"gives {cff_orcid}; they must agree before the abstract is "
            f"published")

    author_line = f"**Authors:** {name} — {affiliation}"
    if orcid:
        author_line += f" (ORCID {orcid})"

    return (
        "<!-- GENERATED FILE - do not hand-edit.\n"
        "     Source: paper/paper.tex (the manuscript the PDF is compiled from).\n"
        "     Regenerate: python3 scripts/generate_preprint_abstract.py\n"
        "     Guarded by: test/test_generated_docs_match_their_generators.py\n"
        "     -->\n"
        "\n"
        f"# {title}\n"
        "\n"
        f"{author_line}\n"
        "\n"
        f"**Repository:** {repo} (v{version})\n"
        "\n"
        "## Abstract\n"
        "\n"
        f"{abstract}\n"
        "\n"
        f"<!-- {length} characters against bioRxiv's {BIORXIV_LIMIT}-character "
        "limit;\n"
        "     this is the manuscript's own abstract. The live route is arXiv, and\n"
        "     its metadata field takes the abridged text in "
        "docs/ARXIV_ABSTRACT.md. -->\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--variant", choices=("biorxiv", "arxiv"),
                        default="biorxiv",
                        help="which server's field the abstract is written "
                             "for: bioRxiv takes the manuscript's abstract, "
                             "arXiv's metadata field is shorter")
    parser.add_argument("--out", default=None,
                        help="where to write the document")
    parser.add_argument("--print", dest="print_only", action="store_true",
                        help="print the abstract alone (paste-ready) and exit")
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if the committed document is stale")
    args = parser.parse_args()

    if args.out is None:
        args.out = ARXIV_OUT if args.variant == "arxiv" else OUT
    command = ("scripts/generate_preprint_abstract.py"
               + (" --variant arxiv" if args.variant == "arxiv" else ""))

    if args.print_only:
        tex = open(TEX, encoding="utf-8").read()
        if args.variant == "arxiv":
            print(check_arxiv_abstract(tex))
        else:
            print(latex_abstract(tex))
        return 0

    document = (build_arxiv_document() if args.variant == "arxiv"
                else build_document())

    if args.check:
        if not os.path.isfile(args.out):
            print(f"{args.out} does not exist; run {command}")
            return 1
        current = open(args.out, encoding="utf-8").read()
        if current != document:
            print(f"{args.out} has drifted from paper/paper.tex; regenerate it "
                  f"with {command}")
            return 1
        print(f"{args.out} matches paper/paper.tex")
        return 0

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(document)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
