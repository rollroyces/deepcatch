"""A command a reader copies must name a file that exists.

Why this file exists
--------------------
`docs/USER_GUIDE.md` spent months describing a command line that the
repository does not contain. It documented `evaluate.py`,
`agent2-multimodal-fusion/evaluate.py`, `agent3-longitudinal/run_final_fixed.py`
and `validation_framework_tests.py`, with plausible flags, plausible output
and plausible clinical verdicts -- none of them backed by a file in the tree.
A reader who followed the guide got a usage error on the first command.

Nothing would have caught that on its own. The unit tests call functions,
not command lines; the documentation tests check numbers, not paths. The
guide was written against a layout the repository had two years ago
(`agent1-…`, `agent2-…`, `evaluate.py` at the root) and was never re-read
after the tree was reorganised under `src/`.

A pointer to a missing file is the same defect as a citation that resolves
to nothing: both send a reader confidently to a place where there is
nothing. The difference is that a broken citation is at least visible in a
review; a broken command line reads as working documentation until somebody
runs it.

What this test checks, and what it deliberately does not
-------------------------------------------------------
It checks that every script path named in a fenced code block or a markdown
link of the user-facing documents resolves to a file in this repository.

It does not check that the command's flags exist. That is a different
question (and a harder one: a flag list is a contract, not a path), but the
path is the part that silently rots -- a renamed file breaks the command
outright, while a renamed flag is caught by the command's own parser.

Scope
-----
`USAGE.md` is deliberately excluded: it is a two-repository guide, and most
of its commands run in the companion pipeline repository, not here. Checking
it against this tree would produce false positives, and a guard that cries
wolf gets deleted. The documents below describe this repository only.

The exemption list is explicit and each entry says why. An exemption that
does not explain itself is a bug waiting to happen.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# The user-facing documents that describe this repository.
_SCOPE = (
    "README.md",
    "README.zh-CN.md",
    "CONTRIBUTING.md",
    "docs/USER_GUIDE.md",
    "wiki/Home.md",
    "wiki/Getting-Started.md",
)

# Scripts that live in the companion repository (cfdna-fragmentomics-pipeline)
# and are named by a document inside an explicit companion-repo step.
_EXEMPT = {
    "run_cross_study.py": (
        "lives in the companion pipeline repository; the README names it "
        "inside the step that clones and runs that repository"
    ),
}

# A script path: a directory-ish prefix with an extension, or a bare
# filename with an extension.
_SCRIPT = re.compile(r"[A-Za-z0-9_./-]*[A-Za-z0-9_-]\.(?:py|sh)\b")

# A dotted module import (`from validation.py.config import ...`) is a
# module path, not a file path; the package directory is `validation/py/`.
_IMPORT = re.compile(r"\b(?:from|import)\s+$")

_FENCE = re.compile(r"```.*?```", re.S)
_INLINE = re.compile(r"`([^`\n]+)`")
_LINK = re.compile(r"\]\(([^)]+)\)")

# A command word immediately before the token: the token is what the reader
# runs, so it has to resolve exactly.
_COMMAND_WORD = re.compile(r"\b(?:python3?|bash|sh|zsh|\./)\s+$")


def _repo_filenames() -> set[str]:
    """Every file name in the repository (excluding .git)."""
    names: set[str] = set()
    for path in _REPO_ROOT.rglob("*"):
        if ".git" in path.parts:
            continue
        names.add(path.name)
    return names


def _code_blocks(text: str) -> list[str]:
    """Every fenced code block in `text`."""
    return _FENCE.findall(text)


def _script_tokens(text: str) -> list[tuple[str, str, bool]]:
    """(token, where, strict) triples for every script path named in `text`.

    Fenced code blocks and markdown link targets, plus inline code in prose:
    a guide that says "run `foo.py`" makes the same promise as a guide that
    shows the command.

    `strict` marks a token that is what the reader runs (it follows a
    command word). A strict token has to resolve to that exact path; a
    loose token -- a name quoted in a comment or a directory-tree diagram
    -- only has to exist somewhere in the repository.
    """
    out: list[tuple[str, str, bool]] = []

    def add(chunk: str, source: str) -> None:
        for line in chunk.splitlines():
            for match in _SCRIPT.finditer(line):
                token = match.group(0)
                if not token or token.startswith("."):
                    continue
                if any(ch in token for ch in "<>*"):
                    continue
                # A dotted module import is not a file path.
                before = line[: match.start()]
                if _IMPORT.search(before):
                    continue
                # Is the token what the reader runs, or is it a name quoted in
                # a comment or a directory-tree diagram?
                strict = bool(_COMMAND_WORD.search(before))
                out.append((token, source, strict))

    for block in _code_blocks(text):
        add(block, "code block")
    for span in _INLINE.findall(text):
        add(span, "inline code")
    for target in _LINK.findall(text):
        add(target, "link target")
    return out


def _resolves(token: str, doc: Path, strict: bool, names: set[str]) -> bool:
    """Does `token` name a file that exists?

    A strict token (one the reader runs) has to resolve to that exact path,
    relative to the repository root or to the document's own directory. A
    loose token -- a name quoted in a comment or a tree diagram -- only has
    to exist somewhere in the repository.
    """
    candidates = [
        _REPO_ROOT / token,
        doc.parent / token,
    ]
    if any(candidate.exists() for candidate in candidates):
        return True
    if strict:
        return False
    return "/" not in token and token in names


def test_documented_scripts_exist():
    """Every script a reader is told to run must exist in this repository."""
    offenders: list[str] = []
    checked = 0
    names = _repo_filenames()

    for name in _SCOPE:
        doc = _REPO_ROOT / name
        if not doc.exists():
            offenders.append(f"{name}: the document itself is missing")
            continue
        text = doc.read_text(encoding="utf-8")
        for token, where, strict in _script_tokens(text):
            if token in _EXEMPT:
                continue
            checked += 1
            if not _resolves(token, doc, strict, names):
                offenders.append(f"{name}: {token} ({where})")

    assert not offenders, (
        "a documented command names a file that does not exist:\n  "
        + "\n  ".join(offenders)
        + "\nEither point the command at the file that exists, or remove it. "
          "A reader who follows it gets a usage error instead of a result."
    )
    assert checked > 0, "the extractor found nothing -- the guard is vacuous"


def test_the_extractor_finds_a_known_script():
    """The extractor must actually see a script path, or the test above is
    vacuous.

    The fixture is assembled at run time so that this file does not become
    a second copy of the paths it is checking.
    """
    name = "run_full_validation.py"
    text = f"```bash\npython {name} --help\n```\n"
    tokens = [token for token, _, _ in _script_tokens(text)]
    assert tokens == [name]


def test_the_extractor_ignores_module_imports():
    """`from validation.py.config import ...` names a module, not a file."""
    text = "```python\nfrom validation.py.config import CTDNA_LEVELS\n```\n"
    assert _script_tokens(text) == []


def test_a_command_token_is_strict_and_a_quoted_name_is_not():
    """A token the reader runs must resolve exactly; a name quoted in a
    comment only has to exist somewhere.

    The distinction is what keeps the guard useful: the layout diagram in
    the README names `cross_study_finallydb.py` under `scripts/`, while the
    command `python agent2-multimodal-fusion/evaluate.py` names a path that
    does not exist and must be caught.
    """
    command = _script_tokens("```bash\npython agent2/evaluate.py\n```\n")
    assert command and command[0][2] is True

    quoted = _script_tokens("```\n├── scripts/   # entry points (x.py)\n```\n")
    assert quoted and quoted[0][2] is False


def test_a_missing_command_path_is_an_offender():
    """The resolution rule must reject a path that is not in the tree."""
    doc = _REPO_ROOT / "docs/USER_GUIDE.md"
    names = _repo_filenames()
    assert not _resolves("agent2-multimodal-fusion/evaluate.py", doc, True, names)
    assert _resolves("scripts/adapter_local_cohort.py", doc, True, names)
