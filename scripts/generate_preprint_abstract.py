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

Usage
-----
    python3 scripts/generate_preprint_abstract.py            # write the doc
    python3 scripts/generate_preprint_abstract.py --print     # abstract only
    python3 scripts/generate_preprint_abstract.py --check     # CI mode

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

# bioRxiv rejects an abstract over this many characters.
BIORXIV_LIMIT = 3000

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


def latex_author(tex: str) -> tuple[str, str]:
    """The author name and the affiliation line under it."""
    match = re.search(r"\\author\{(.*?)\}\n", tex, re.S)
    if not match:
        raise ValueError("paper.tex has no \\author")
    lines = [ln.strip() for ln in match.group(1).split("\\\\") if ln.strip()]
    name = plain(render_text(lines[0]))
    affiliation = ""
    for line in lines[1:]:
        text = re.sub(r"\\small\s*", "", line).strip()
        text = re.sub(r"\\url\{([^}]*)\}", r"\1", text)
        text = text.replace("---", "—")
        if text and "github.com" not in text:
            affiliation = text
    return name, affiliation


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


def build_document() -> str:
    tex = open(TEX, encoding="utf-8").read()
    cff = open(CITATION, encoding="utf-8").read()

    title = latex_title(tex)
    name, affiliation = latex_author(tex)
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

    return (
        "<!-- GENERATED FILE - do not hand-edit.\n"
        "     Source: paper/paper.tex (the manuscript the PDF is compiled from).\n"
        "     Regenerate: python3 scripts/generate_preprint_abstract.py\n"
        "     Guarded by: test/test_generated_docs_match_their_generators.py\n"
        "     -->\n"
        "\n"
        f"# {title}\n"
        "\n"
        f"**Authors:** {name} — {affiliation}\n"
        "\n"
        f"**Repository:** {repo} (v{version})\n"
        "\n"
        "## Abstract\n"
        "\n"
        f"{abstract}\n"
        "\n"
        f"<!-- {length} characters against bioRxiv's {BIORXIV_LIMIT}-character "
        "limit;\n"
        "     this is the text to paste into the submission form. -->\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", default=OUT,
                        help="where to write the document")
    parser.add_argument("--print", dest="print_only", action="store_true",
                        help="print the abstract alone (paste-ready) and exit")
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if the committed document is stale")
    args = parser.parse_args()

    if args.print_only:
        tex = open(TEX, encoding="utf-8").read()
        print(latex_abstract(tex))
        return 0

    document = build_document()

    if args.check:
        if not os.path.isfile(args.out):
            print(f"{args.out} does not exist; run "
                  f"scripts/generate_preprint_abstract.py")
            return 1
        current = open(args.out, encoding="utf-8").read()
        if current != document:
            print(f"{args.out} has drifted from paper/paper.tex; regenerate it "
                  f"with scripts/generate_preprint_abstract.py")
            return 1
        print(f"{args.out} matches paper/paper.tex")
        return 0

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(document)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
