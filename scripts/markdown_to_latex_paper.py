#!/usr/bin/env python3
"""Render the methods manuscript's markdown as a compilable LaTeX document.

Why this exists
---------------
`paper/METHODS_PAPER.md` is the methods manuscript. It is written in
markdown because that is what the repository edits, reviews and
diffs -- a manuscript that is only ever edited as LaTeX stops being
readable in a pull request. But the file that is eventually submitted
has to be LaTeX, and the two have to be the same document: a
hand-transcribed LaTeX copy of a markdown manuscript is two documents
that agree until the day somebody edits one of them.

So the LaTeX is generated, and this script is the generator. It reads
the markdown and writes `paper/METHODS_PAPER.tex`. It does not invent
prose: every sentence in the `.tex` comes from a sentence in the
`.md`. The captions of the figures are the one place where a sentence
is written here rather than translated, because the markdown names its
figures but does not caption them.

Two checks keep the pair honest, and `--check` runs both:

* The `.tex` is regenerated and compared with the committed copy, so a
  hand edit to the generated file is caught.
* The markdown's text and the LaTeX's text are compared token by
  token, so a converter bug -- a dropped sentence, a changed number --
  is caught even though both files come from the same run.

The LaTeX targets `pdflatex`, because the repository's CI compiles the
other manuscripts with it. That is why the unicode punctuation is
translated into LaTeX macros below rather than left as UTF-8: XeTeX
would take them as they are, and `pdflatex` would not.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_MD = _REPO_ROOT / "paper" / "METHODS_PAPER.md"
_DEFAULT_TEX = _REPO_ROOT / "paper" / "METHODS_PAPER.tex"


# --------------------------------------------------------------------------
# Inline translation
# --------------------------------------------------------------------------

# The characters LaTeX reserves for itself. Each becomes its own escaped
# form; none of them can survive a raw copy into the document.
_SPECIALS = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
    "<": r"\textless{}",
    ">": r"\textgreater{}",
}

# The unicode punctuation the manuscript uses. `pdflatex` with a T1 font
# has the dashes and the text symbols, but the mathematical ones are
# written as macros here so that they are set in the math font -- and so
# that the round-trip check has an exact inverse for each of them.
_UNICODE = {
    "—": "---",              # em dash
    "–": "--",               # en dash
    "−": "$-$",              # minus sign
    "§": r"\S{}",            # section sign
    "±": r"$\pm$",
    "×": r"$\times$",
    "→": r"$\rightarrow$",
    "≈": r"$\approx$",
    "≥": r"$\geq$",
    "Δ": r"$\Delta$",
    "β": r"$\beta$",
    "α": r"$\alpha$",
    "·": r"\textperiodcentered{}",
    "¹": r"\textsuperscript{1}",
    "²": r"\textsuperscript{2}",
    "³": r"\textsuperscript{3}",
    "⁻": r"\textsuperscript{-}",
}


def escape(text: str) -> str:
    """Escape LaTeX's reserved characters and translate the punctuation."""
    out = []
    for ch in text:
        if ch in _SPECIALS:
            out.append(_SPECIALS[ch])
        elif ch in _UNICODE:
            out.append(_UNICODE[ch])
        else:
            out.append(ch)
    return "".join(out)


def _skip_code_span(text: str, index: int) -> int:
    """The index just past the code span that starts at `index`.

    A code span ends at the next backtick. If there is none, the span
    runs to the end of the text -- which is what a markdown renderer
    does too, and it keeps the scanner from mistaking a stray backtick
    for emphasis.
    """
    closing = text.find("`", index + 1)
    return len(text) if closing == -1 else closing + 1


def _find_close(text: str, start: int, delimiter: str) -> int:
    """The index of the closing `delimiter`, or -1.

    Code spans are skipped while looking, so a `*` inside a code span
    cannot close an emphasis span -- which matters here because the
    manuscript writes `` `*.frag.tsv.bgz` `` three times.
    """
    index = start
    while index < len(text):
        ch = text[index]
        if ch == "`":
            index = _skip_code_span(text, index)
            continue
        if text.startswith(delimiter, index):
            return index
        index += 1
    return -1


_AUTOLINK = re.compile(r"<(https?://[^>\s]+)>")

# A code span that names a path is set with the url package's command
# rather than with \texttt, because a path is one long word to TeX and
# runs out of the margin instead of wrapping. The url package's command
# may be used inside another command's argument, which \url may not.
_PATH_SAFE = re.compile(r"^[A-Za-z0-9_.\-/*:~+ ]+$")


def _is_path(text: str) -> bool:
    """Is this code span a path that should be allowed to wrap?"""
    return "/" in text and bool(_PATH_SAFE.match(text))


# A break point inside a code span or a link label. A file path or a
# DOI is one long word to TeX, and without these it runs out of the
# margin instead of wrapping. Only the separators of a path are break
# points: a break after every hyphen would split `--motif-cohort-mode`
# into something that reads as a typo.
_BREAK_AFTER = re.compile(r"(\\_|/)")


def _breakable(text: str) -> str:
    """Let TeX break a path or an identifier where it reads naturally."""
    return _BREAK_AFTER.sub(r"\1\\allowbreak{}", text)


def inline(text: str) -> str:
    """Translate one line of markdown into LaTeX."""
    out = []
    buffer = []
    index = 0

    def flush() -> None:
        if buffer:
            # An equals sign is a natural place to break a line: the
            # manuscript writes `sens@95%=0.865` often enough that
            # keeping it whole overflows the margin.
            out.append(escape("".join(buffer)).replace(
                "=", r"=\allowbreak{}"))
            buffer.clear()

    while index < len(text):
        ch = text[index]
        if ch == "`":
            flush()
            end = _skip_code_span(text, index)
            code = text[index + 1:end - 1]
            if _is_path(code):
                out.append(r"\pathurl{" + code + "}")
            else:
                out.append(r"\texttt{" + _breakable(escape(code)) + "}")
            index = end
        elif text.startswith("**", index):
            closing = _find_close(text, index + 2, "**")
            if closing == -1:
                buffer.append(ch)
                index += 1
                continue
            flush()
            out.append(r"\textbf{" + inline(text[index + 2:closing]) + "}")
            index = closing + 2
        elif ch == "*":
            closing = _find_close(text, index + 1, "*")
            if closing == -1:
                buffer.append(ch)
                index += 1
                continue
            flush()
            out.append(r"\textit{" + inline(text[index + 1:closing]) + "}")
            index = closing + 1
        elif ch == "[":
            link = _link_at(text, index)
            if link is None:
                buffer.append(ch)
                index += 1
                continue
            label, url, end = link
            flush()
            out.append(r"\href{" + escape(url) + "}{"
                       + _breakable(inline(label)) + "}")
            index = end
        elif ch == "<":
            autolink = _AUTOLINK.match(text, index)
            if autolink is None:
                buffer.append(ch)
                index += 1
                continue
            flush()
            out.append(r"\textless{}\url{" + escape(autolink.group(1))
                       + r"}\textgreater{}")
            index = autolink.end()
        else:
            buffer.append(ch)
            index += 1
    flush()
    return "".join(out)


def _link_at(text: str, index: int):
    """(label, url, end) for the markdown link at `index`, or None.

    The URL is scanned with balanced parentheses: the manuscript's DOIs
    contain them (`10.1016/S0140-6736(97)02174-0`), and a link that
    stopped at the first `)` would be truncated.
    """
    depth = 0
    cursor = index
    while cursor < len(text):
        if text[cursor] == "[":
            depth += 1
        elif text[cursor] == "]":
            depth -= 1
            if depth == 0:
                break
        cursor += 1
    else:
        return None
    if cursor + 1 >= len(text) or text[cursor + 1] != "(":
        return None
    label = text[index + 1:cursor]

    paren_depth = 1
    url_start = cursor + 2
    cursor = url_start
    while cursor < len(text):
        if text[cursor] == "(":
            paren_depth += 1
        elif text[cursor] == ")":
            paren_depth -= 1
            if paren_depth == 0:
                break
        cursor += 1
    else:
        return None
    return label, text[url_start:cursor], cursor + 1


# --------------------------------------------------------------------------
# Block parsing
# --------------------------------------------------------------------------

_HEADING = re.compile(r"^(#{1,6})\s+(.*)$")
_HR = re.compile(r"^-{3,}\s*$")
_UL = re.compile(r"^(\s*)[-*]\s+(.*)$")
_OL = re.compile(r"^(\s*)(\d+)\.\s+(.*)$")
_QUOTE = re.compile(r"^>\s?(.*)$")
_TABLE_SEP = re.compile(r"^\|?[\s:\-|]+\|[\s:\-|]*$")
_LEADING_NUMBER = re.compile(r"^\d+(?:\.\d+)*\.?\s+")
_FIGURE_REF = re.compile(r"docs/figures/([A-Za-z0-9_.\-]+\.png)")


class Block:
    """One block of the markdown, with the text it carries."""

    def __init__(self, kind: str, **kwargs):
        self.kind = kind
        self.__dict__.update(kwargs)

    def raw_text(self) -> str:
        """The block's markdown, for looking up what it mentions."""
        if self.kind == "para":
            return self.text
        if self.kind == "quote":
            return "\n".join(self.paragraphs)
        if self.kind == "heading":
            return self.text
        return ""


def parse_blocks(lines: list[str]) -> list[Block]:
    """The markdown's blocks, in order."""
    blocks: list[Block] = []
    index = 0
    total = len(lines)

    while index < total:
        line = lines[index]
        stripped = line.strip()

        if not stripped:
            index += 1
            continue

        if _HR.match(stripped):
            blocks.append(Block("hr"))
            index += 1
            continue

        heading = _HEADING.match(line)
        if heading:
            blocks.append(Block("heading", level=len(heading.group(1)),
                                text=heading.group(2).strip()))
            index += 1
            continue

        # A table is a header row plus a separator row.
        if (stripped.startswith("|") and index + 1 < total
                and _TABLE_SEP.match(lines[index + 1].strip())
                and "-" in lines[index + 1]):
            header = _split_row(lines[index])
            aligns = _split_row(lines[index + 1])
            index += 2
            rows = []
            while index < total and lines[index].strip().startswith("|"):
                rows.append(_split_row(lines[index]))
                index += 1
            blocks.append(Block("table", header=header, aligns=aligns, rows=rows))
            continue

        if _QUOTE.match(line):
            paragraphs = []
            current = []
            while index < total and (_QUOTE.match(lines[index]) or
                                    (lines[index].strip() == "" and index + 1 < total
                                     and _QUOTE.match(lines[index + 1]))):
                match = _QUOTE.match(lines[index])
                if match is None:
                    if current:
                        paragraphs.append(" ".join(current))
                        current = []
                    index += 1
                    continue
                current.append(match.group(1).strip())
                index += 1
            if current:
                paragraphs.append(" ".join(current))
            blocks.append(Block("quote", paragraphs=paragraphs))
            continue

        if _UL.match(line) or _OL.match(line):
            items, index = _parse_list(lines, index)
            blocks.append(Block("list", items=items))
            continue

        # A paragraph runs until a blank line or the start of another block.
        paragraph = []
        while index < total:
            candidate = lines[index]
            if not candidate.strip():
                break
            if (paragraph and (_HEADING.match(candidate) or _HR.match(candidate.strip())
                               or candidate.strip().startswith("|")
                               or _QUOTE.match(candidate)
                               or _UL.match(candidate) or _OL.match(candidate))):
                break
            paragraph.append(candidate.strip())
            index += 1
        blocks.append(Block("para", text=" ".join(paragraph)))

    return blocks


def _split_row(row: str) -> list[str]:
    """The cells of a markdown table row."""
    stripped = row.strip()
    if stripped.startswith("|"):
        stripped = stripped[1:]
    if stripped.endswith("|"):
        stripped = stripped[:-1]
    return [cell.strip() for cell in stripped.split("|")]


def _parse_list(lines: list[str], index: int):
    """A list's items, with their nesting, and the index after it."""
    items = []
    total = len(lines)
    while index < total:
        line = lines[index]
        if not line.strip():
            # A blank line ends the list unless the next line continues it.
            if index + 1 < total and (_UL.match(lines[index + 1]) or _OL.match(lines[index + 1])):
                index += 1
                continue
            break

        match = _OL.match(line)
        ordered = match is not None
        if match is None:
            match = _UL.match(line)
        if match is None:
            # A continuation line belongs to the previous item.
            if items:
                items[-1]["text"].append(line.strip())
                index += 1
                continue
            break

        if ordered:
            number = int(match.group(2))
            text = match.group(3)
        else:
            number = None
            text = match.group(2)

        items.append({"indent": len(match.group(1)), "ordered": ordered,
                      "number": number, "text": [text]})
        index += 1

    # The nesting is the indentation, relative to the shallowest item.
    if items:
        base = min(item["indent"] for item in items)
        for item in items:
            item["level"] = (item["indent"] - base) // 2

    return _nest(items), index


def _nest(items: list[dict]) -> list[dict]:
    """Turn a flat list of items into a tree."""
    root: list[dict] = []
    stack: list[tuple[int, dict]] = []
    for item in items:
        node = {"ordered": item["ordered"], "number": item["number"],
                "text": " ".join(item["text"]), "children": []}
        while stack and stack[-1][0] >= item["level"]:
            stack.pop()
        if stack:
            stack[-1][1]["children"].append(node)
        else:
            root.append(node)
        stack.append((item["level"], node))
    return root


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

# The figures the manuscript names, with the captions the markdown does
# not carry. This is the only prose in this file that is not translated
# from the manuscript, and it is deliberately taken from the sentences
# that already describe each figure where such a sentence exists.
FIGURES = {
    "fig1_pooled_roc.png": (
        "Pooled out-of-fold ROC for the cross-study FinaleDB cohort "
        "(publications 6 + 8, n=627), in the GC-corrected configuration "
        "that is the headline."),
    "fig2_per_cancer_roc.png": (
        "Per-seed out-of-fold ROC for the top-5 cancers by count, scored "
        "against all healthy samples."),
    "fig3_calibration.png": (
        "Sensitivity at fixed specificity for the pooled operating points, "
        "in ROC space."),
    "fig4_sens_operating_points.png": (
        "Grouped sensitivity at 95, 98 and 99 per cent specificity, per "
        "cancer."),
    "fig5_encoder_architecture.png": (
        "The six modality slots, the per-modality projection, the "
        "[CLS] token and the four-layer transformer, the two pre-training "
        "phases that produce the checkpoint, and the 0.7 / 0.3 blend that "
        "turns the encoder into a score."),
    "fig6_validation_protocol.png": (
        "The per-patient split that replaced the leaky StratifiedKFold, "
        "the out-of-fold scores it produces, the metrics that are "
        "reported, and the three controls -- the true-confound control, "
        "the shuffled-label control, and the bit-identical sanity check."),
}


def _figure_number(name: str) -> int:
    match = re.match(r"fig(\d+)_", name)
    return int(match.group(1)) if match else 0


def _column_widths(rows: list[list[str]], count: int) -> list[float]:
    """The share of the text width each column gets.

    Proportional to the square root of its longest cell, so a column of
    prose gets most of the width and a column of numbers gets what it
    needs. The shares sum to one, and each column's width subtracts the
    padding on both sides of it, so the table is exactly as wide as the
    text.
    """
    weights = []
    for column in range(count):
        longest = 0
        for row in rows:
            if column < len(row):
                longest = max(longest, len(row[column]))
        weights.append(max(6.0, min(longest, 160) ** 0.55))
    total = sum(weights) or 1.0
    shares = [weight / total for weight in weights]
    floor = min(0.07, 1.0 / max(count, 1))
    shares = [max(share, floor) for share in shares]
    scale = 1.0 / sum(shares)
    return [share * scale for share in shares]


def _alignment(cell: str) -> str:
    """The alignment the separator row asks a column for."""
    cell = cell.strip()
    left = cell.startswith(":")
    right = cell.endswith(":")
    if left and right:
        return "c"
    if right:
        return "r"
    return "l"


def _render_table(block: Block) -> list[str]:
    header = block.header
    rows = block.rows
    count = max(len(header), max((len(row) for row in rows), default=0))
    aligns = [_alignment(block.aligns[i]) if i < len(block.aligns) else "l"
              for i in range(count)]
    shares = _column_widths([header] + rows, count)

    spec = []
    for share, align in zip(shares, aligns):
        width = r"p{\dimexpr" + f"{share:.4f}" + r"\textwidth-2\tabcolsep\relax}"
        if align == "r":
            spec.append(r">{\raggedleft\arraybackslash}" + width)
        elif align == "c":
            spec.append(r">{\centering\arraybackslash}" + width)
        else:
            # A narrow column of justified text is a column of rivers.
            spec.append(r">{\raggedright\arraybackslash}" + width)

    out = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{" + "".join(spec) + "}",
        r"\toprule",
    ]
    out.append(" & ".join(inline(cell) for cell in header) + r" \\")
    out.append(r"\midrule")
    for row in rows:
        padded = list(row) + [""] * (count - len(row))
        out.append(" & ".join(inline(cell) for cell in padded) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{tabular}")
    out.append(r"\end{table}")
    return out


def _render_list(nodes: list[dict], counter: list[int] | None = None) -> list[str]:
    """A list's LaTeX, with the numbers the markdown gave its items."""
    if counter is None:
        counter = [0]
    out: list[str] = []
    environment = "enumerate" if nodes and nodes[0]["ordered"] else "itemize"
    out.append(r"\begin{" + environment + "}")
    for node in nodes:
        if node["ordered"]:
            counter[0] += 1
            if node["number"] is not None and node["number"] != counter[0]:
                # The markdown's own numbering is not the sequential one
                # `enumerate` would produce, so it is stated explicitly.
                out.append(r"\item[" + str(node["number"]) + r"]")
            else:
                out.append(r"\item")
        else:
            out.append(r"\item")
        out.append(inline(node["text"]))
        if node["children"]:
            out.extend(_render_list(node["children"], counter))
    out.append(r"\end{" + environment + "}")
    return out


def _render_figure(name: str) -> list[str]:
    number = _figure_number(name)
    caption = FIGURES.get(name, "")
    return [
        r"\begin{center}",
        r"\includegraphics[width=0.92\textwidth]{" + name + "}",
        "",
        r"\textbf{Figure " + str(number) + ".} " + inline(caption),
        r"\end{center}",
    ]


# --------------------------------------------------------------------------
# The document
# --------------------------------------------------------------------------

PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{textcomp}
\usepackage{amsmath, amssymb}
\usepackage{booktabs}
\usepackage{array}
\usepackage{graphicx}
\usepackage[obeyspaces]{url}
\usepackage[hidelinks]{hyperref}
\DeclareUrlCommand\pathurl{\urlstyle{tt}}
\graphicspath{{../docs/figures/}{docs/figures/}}
\setlength{\tabcolsep}{3pt}
\setlength{\emergencystretch}{2em}
% The manuscript's figures are named in its own text, so the numbers
% below are the manuscript's numbers, not LaTeX's counter. The section
% numbers are LaTeX's, and the generator checks that they agree with
% the numbers the markdown wrote out by hand.
%
% There is no `inputenc` here because the generated file is ASCII: the
% generator translates every unicode character the manuscript uses, and
% refuses to write a file that still holds one.
"""


def render(markdown: str) -> str:
    """The LaTeX document for a markdown manuscript."""
    lines = markdown.split("\n")
    blocks = parse_blocks(lines)

    title = ""
    author_lines: list[str] = []
    notes: list[str] = []
    body: list[Block] = []
    in_front = True
    front_paragraphs: list[str] = []
    for block in blocks:
        if in_front and block.kind == "heading" and block.level == 2:
            in_front = False
        if in_front:
            if block.kind == "heading" and block.level == 1:
                title = block.text
            elif block.kind == "para":
                front_paragraphs.append(block.text)
        else:
            body.append(block)

    # The front matter is positional: the author line, the affiliation
    # line under it, and then the notes. The round-trip check fails
    # loudly if the manuscript's front matter stops looking like this.
    if front_paragraphs and front_paragraphs[0].startswith("**Authors:**"):
        author_lines.append(front_paragraphs[0][len("**Authors:**"):].strip())
        rest = front_paragraphs[1:]
    else:
        rest = front_paragraphs
    if rest and rest[0].startswith("¹"):
        author_lines.append(rest[0])
        rest = rest[1:]
    notes = rest

    out: list[str] = [PREAMBLE]
    out.append(r"\title{" + inline(title) + "}")
    # The author block is one line in the markdown and three on the
    # title page, because one line does not fit. The words are the
    # same; the middle dot is where the line is broken.
    author = " \\\\\n".join(
        inline(chunk)
        for line in author_lines
        for chunk in line.split(" · "))
    out.append(r"\author{" + author + "}")
    out.append(r"\date{\today}")
    out.append("")
    out.append(r"\begin{document}")
    out.append(r"\maketitle")
    if notes:
        out.append(r"{\small")
        for note in notes:
            out.append(inline(note))
            out.append("")
        out.append(r"\par}")
        out.append("")

    seen_figures: set[str] = set()

    def emit_figures(block: Block) -> None:
        """The figures a block mentions, the first time it mentions them.

        The manuscript names its figures in the prose, and two of them
        are named only in the note that lists them; the figure goes
        where the note is, and the number on it is the manuscript's
        number rather than LaTeX's counter.
        """
        for name in _FIGURE_REF.findall(block.raw_text()):
            if name not in seen_figures and name in FIGURES:
                seen_figures.add(name)
                out.extend(_render_figure(name))
                out.append("")

    in_abstract = False
    for block in body:
        if block.kind == "heading" and block.level == 2:
            heading = block.text
            if heading.lower() == "abstract":
                out.append(r"\begin{abstract}")
                in_abstract = True
                continue
            if in_abstract:
                out.append(r"\end{abstract}")
                in_abstract = False
            text = _LEADING_NUMBER.sub("", heading)
            out.append(r"\section{" + inline(text) + "}")
            continue
        if block.kind == "heading":
            level = min(block.level, 4)
            command = {3: "subsection", 4: "subsubsection"}[level]
            text = _LEADING_NUMBER.sub("", block.text)
            out.append(r"\%s{%s}" % (command, inline(text)))
            continue

        if block.kind == "hr":
            out.append(r"\par\bigskip\hrule\bigskip")
            continue

        if block.kind == "table":
            out.extend(_render_table(block))
            out.append("")
            continue

        if block.kind == "quote":
            out.append(r"\begin{quote}")
            for paragraph in block.paragraphs:
                out.append(inline(paragraph))
                out.append("")
            out.append(r"\end{quote}")
            out.append("")
            emit_figures(block)
            continue

        if block.kind == "list":
            out.extend(_render_list(block.items))
            out.append("")
            continue

        out.append(inline(block.text))
        out.append("")
        emit_figures(block)

    if in_abstract:
        out.append(r"\end{abstract}")
    out.append(r"\end{document}")
    document = "\n".join(out) + "\n"
    non_ascii = sorted({ch for ch in document if ord(ch) > 126})
    if non_ascii:
        raise ValueError(
            "the generated LaTeX still holds unicode characters that "
            f"pdflatex would choke on: {non_ascii!r}. Add them to "
            f"_UNICODE so that they are translated.")
    return document


# --------------------------------------------------------------------------
# The round-trip check
# --------------------------------------------------------------------------

_TOKEN = re.compile(r"\w+(?:[.\-/]\w+)*", re.UNICODE)

# The punctuation that separates words in one language and is a word
# character in the other. Each is turned into a separator on both sides.
_MD_SEPARATORS = "§±×→≈≥Δβα·"

# The superscripts, each turned into the digits it stands for, so that
# `10⁻³` reads the same on both sides.
_MD_SUPERSCRIPTS = {"¹": "1", "²": "2", "³": "3", "⁻": "-"}


def _strip_links(text: str) -> str:
    """Every markdown link replaced by its label.

    The links are found with the same balanced-parenthesis scan the
    generator uses, because the manuscript's DOIs contain parentheses
    and a naive scan would truncate them.
    """
    out = []
    index = 0
    while index < len(text):
        if text[index] != "[":
            out.append(text[index])
            index += 1
            continue
        link = _link_at(text, index)
        if link is None:
            out.append(text[index])
            index += 1
            continue
        label, _url, end = link
        out.append(" " + label + " ")
        index = end
    return "".join(out)


def _md_tokens(markdown: str) -> list[str]:
    """The words of the markdown, in order, with the markup removed."""
    text = markdown
    # Drop the front matter: the LaTeX title and author block carry it
    # in a different shape, and it is checked separately.
    body = text.split("## Abstract", 1)
    if len(body) == 2:
        text = body[1]

    text = re.sub(r"^\s*\|[\s:\-|]+\|[\s:\-|]*$", " ", text, flags=re.M)
    # The headings' manual numbers are LaTeX's job; that the two agree
    # is checked by `_heading_numbering`.
    text = re.sub(r"^#{1,6}\s+(?:\d+(?:\.\d+)*\.?\s+)?", "", text, flags=re.M)
    text = re.sub(r"^\s*>\s?", "", text, flags=re.M)
    text = re.sub(r"^\s*[-*]\s+", " ", text, flags=re.M)
    text = re.sub(r"^\s*\d+\.\s+", " ", text, flags=re.M)
    text = re.sub(r"^-{3,}\s*$", " ", text, flags=re.M)
    text = _strip_links(text)
    text = text.replace("`", " ")
    text = text.replace("**", " ").replace("*", " ")
    text = text.replace("|", " ")
    for ch in _MD_SEPARATORS:
        text = text.replace(ch, " ")
    for ch, digits in _MD_SUPERSCRIPTS.items():
        text = text.replace(ch, " " + digits + " ")
    return _TOKEN.findall(text)


_TEX_LITERALS = {
    r"\_": "_", r"\%": "%", r"\&": "&", r"\$": "$", r"\#": "#",
    r"\{": "{", r"\}": "}", r"\textbackslash{}": "\\",
}

_TEX_COMMANDS = (
    r"\\(?:begin|end)\{(?:itemize|enumerate|quote|table|tabular|abstract|"
    r"document|center)\}(?:\[[^\]]*\])?",
    r"\\item(?:\[[^\]]*\])?",
    r"\\(?:toprule|midrule|bottomrule)\b",
    r"\\includegraphics(?:\[[^\]]*\])?\{[^}]*\}",
    r"\\label\{[^}]*\}",
    r"\\href\{[^}]*\}\{",
    r"\\(?:par|bigskip|hrule|maketitle|noindent|centering)\b",
)


def _tex_tokens(tex: str) -> list[str]:
    """The words of the LaTeX, in order, with the markup removed."""
    # Both sides of the comparison start at the abstract: the title, the
    # author block and the notes are shaped differently by the two
    # languages, and they are checked separately.
    body = tex.split(r"\begin{abstract}", 1)
    if len(body) == 2:
        text = body[1]

    # The figures are rendered, not read: their captions are written
    # here rather than translated, so they are not part of the text.
    text = re.sub(
        r"\\begin\{center\}\s*\\includegraphics.*?\\end\{center\}",
        " ", text, flags=re.S)

    # The tabular preamble is a column specification, not text.
    text = re.sub(r"\\begin\{tabular\}\{", "@@TABULAR@@", text)
    text = re.sub(r"@@TABULAR@@(?:[^{}]|\{[^{}]*\})*\}", " ", text)

    for pattern in _TEX_COMMANDS:
        text = re.sub(pattern, " ", text)

    # A superscript keeps its digits: `10⁻³` is three tokens on both
    # sides, not two.
    text = re.sub(r"\\textsuperscript\{([^}]*)\}", r" \1 ", text)

    for literal, replacement in _TEX_LITERALS.items():
        text = text.replace(literal, replacement)
    # The break points are invisible in the text, so they are removed
    # without leaving a space behind them.
    text = text.replace(r"\allowbreak{}", "")
    text = re.sub(r"\\[A-Za-z]+", " ", text)
    text = text.replace("{", " ").replace("}", " ")
    text = text.replace("$", " ").replace("&", " ").replace("\\\\", " ")
    for ch in "±×→≈≥−§":
        text = text.replace(ch, " ")
    return _TOKEN.findall(text)


def _front_matter_checks(markdown: str, tex: str) -> list[str]:
    """The checks the token comparison cannot make about the front matter."""
    problems = []
    title = re.search(r"^#\s+(.*)$", markdown, flags=re.M)
    if title:
        wanted = title.group(1).strip()
        if wanted not in tex:
            problems.append(f"the title is not in the LaTeX: {wanted!r}")
    for wanted in ("Yu Ching Lam", "Independent Researcher",
                   "0009-0008-9113-769X",
                   "https://orcid.org/0009-0008-9113-769X"):
        if wanted not in tex:
            problems.append(f"the author block does not carry {wanted!r}")
    return problems


def _heading_numbering(markdown: str) -> list[str]:
    """Whether LaTeX's section numbers are the markdown's section numbers.

    The markdown writes its section numbers out by hand ("3.2.1 The
    estimator class is the problem") and the LaTeX leaves them to
    LaTeX's counters. Those two have to be the same numbers, or the
    prose -- which says "see §3.2.1" -- points at the wrong place in
    the compiled document.
    """
    problems = []
    section = subsection = subsubsection = 0
    for line in markdown.split("\n"):
        match = re.match(r"^(#{2,4})\s+(.*)$", line)
        if not match:
            continue
        level = len(match.group(1))
        heading = match.group(2).strip()
        if level == 2 and heading.lower() == "abstract":
            # The abstract is an environment, not a numbered section.
            continue
        if level == 2:
            section += 1
            subsection = 0
            subsubsection = 0
            expected = f"{section}"
        elif level == 3:
            subsection += 1
            subsubsection = 0
            expected = f"{section}.{subsection}"
        else:
            subsubsection += 1
            expected = f"{section}.{subsection}.{subsubsection}"

        written = re.match(r"^(\d+(?:\.\d+)*)\.?\s+", heading)
        if not written:
            problems.append(
                f"the heading {heading!r} carries no section number, but "
                f"LaTeX will number it {expected}")
        elif written.group(1) != expected:
            problems.append(
                f"the heading {heading!r} is numbered {written.group(1)} in "
                f"the markdown, but LaTeX's counters will call it {expected}")
    return problems


def check(markdown: str, tex: str) -> list[str]:
    """Everything wrong with the pair, as a list of messages."""
    problems: list[str] = []
    md_tokens = _md_tokens(markdown)
    tex_tokens = _tex_tokens(tex)
    if md_tokens != tex_tokens:
        for index, (left, right) in enumerate(zip(md_tokens, tex_tokens)):
            if left != right:
                context = " ".join(md_tokens[max(0, index - 8):index + 8])
                problems.append(
                    f"the texts diverge at token {index}: markdown says "
                    f"{left!r}, LaTeX says {right!r} (around: {context!r})")
                break
        else:
            problems.append(
                f"the texts have different lengths: markdown has "
                f"{len(md_tokens)} tokens, LaTeX has {len(tex_tokens)}; "
                f"markdown ends {md_tokens[len(tex_tokens):][:12]!r}, "
                f"LaTeX ends {tex_tokens[len(md_tokens):][:12]!r}")
    problems.extend(_heading_numbering(markdown))
    problems.extend(_front_matter_checks(markdown, tex))
    return problems


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--md", default=str(_DEFAULT_MD),
                        help="the markdown manuscript to read")
    parser.add_argument("--out", default=str(_DEFAULT_TEX),
                        help="the LaTeX manuscript to write")
    parser.add_argument("--check", action="store_true",
                        help="do not write: compare the committed LaTeX with "
                             "a fresh generation and compare the two texts")
    args = parser.parse_args(argv)

    md_path = Path(args.md)
    tex_path = Path(args.out)
    markdown = md_path.read_text(encoding="utf-8")
    generated = render(markdown)

    if args.check:
        problems = []
        if not tex_path.is_file():
            problems.append(f"{tex_path} does not exist; run the generator")
        else:
            committed = tex_path.read_text(encoding="utf-8")
            if committed != generated:
                problems.append(
                    f"{tex_path} has drifted from {md_path.name}. Regenerate "
                    f"it; the generated file must not be hand-edited, because "
                    f"the next generation erases the edit.")
        problems.extend(check(markdown, generated))
        for problem in problems:
            print(f"FAIL: {problem}", file=sys.stderr)
        if problems:
            return 1
        print(f"OK: {tex_path.name} is what {md_path.name} generates, and the "
              f"two texts agree.")
        return 0

    tex_path.parent.mkdir(parents=True, exist_ok=True)
    tex_path.write_text(generated, encoding="utf-8")
    print(f"wrote {tex_path} ({len(generated.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
