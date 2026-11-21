"""A documented install must provide what the documented command needs.

What went wrong
---------------
`pytest` is a development dependency: `pyproject.toml` puts it in the
`dev` extra and not in the runtime `dependencies`. Every quickstart in
the repository installed the base distribution and then ran the test
suite --

    pip install -e .
    pytest -q                 # 953 tests, all pass

-- which fails on a clean virtualenv with `No module named pytest`
before a single test runs. The README's install-and-validate block, the
Chinese README, `CONTRIBUTING.md`, `docs/USER_GUIDE.md`,
`docs/REPRODUCE.md`, `USAGE.md`, the docs site's quickstart and the
install demo script all had it. The demo was worse twice over: it
installed with the ambient `pip` rather than the venv's own, and then
ran the venv's python.

The wiki made a quieter version of the same claim. It listed `torch`
and `statsmodels` under "Core Python dependencies (installed
automatically by `pip install -e .`)". Neither is a runtime dependency
-- `torch` is in the `api` extra and `statsmodels` is in no extra at
all -- so a reader who installed and expected torch had no torch.

Two numbers in the README had drifted the same way. The install block
claimed "15 passed" for the publication-readiness suite (it collects
16), and the test section claimed "702 passed, 19 skipped, 13
deselected" for the fast-path gate (it selects 755 of 768). Both were
unchecked: the count guard compares the documents with each other and
pins the suite's total, but nobody compared a claim with the command
that produces it.

What this test checks, and what it deliberately does not
-------------------------------------------------------
1. A document that runs `pytest` installs something that provides
   `pytest` first -- the `dev` extra, or `pytest` named as a package --
   or checks that `pytest` is importable and says what to install.
2. The premise of rule 1 holds in `pyproject.toml`: `pytest` is in the
   `dev` extra and is not a runtime dependency. If a future change
   moved `pytest` into the runtime dependencies, every quickstart would
   start working again and this rule would have to change with it, so
   the premise is pinned rather than assumed.
3. A document that says a list of packages is "installed automatically
   by `pip install -e .`" lists only packages that really are runtime
   dependencies.
4. The outcome a document claims for a documented test command is the
   number that command produces.

Rule 4 checks the TOTAL, not the split between passed and skipped: the
total comes from a collection (fast), while the split needs a full run
(minutes). A claim whose numbers do not add up to the collection is
rejected, and the message says to run the command for the real split.

Scope
-----
The documents below are the ones a reader installs from. `USAGE.md` is a
two-repository guide and is included for the install rule only: it
names the companion pipeline's `.[dev]` extra, which is the same shape
of extra.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# The documents a reader installs from and then runs a command in.
_DOCS = (
    "README.md",
    "README.zh-CN.md",
    "CONTRIBUTING.md",
    "USAGE.md",
    "docs/USER_GUIDE.md",
    "docs/REPRODUCE.md",
    "docs/PRETRAINING.md",
    "docs/CROSS_PLATFORM_BENCHMARK.md",
    "paper/METHODS_PAPER_STATUS.md",
    "wiki/Getting-Started.md",
    "docs/site/index.html",
    "docs/demo/terminal/install_demo.sh",
    "paper/REPRODUCE.sh",
    "scripts/reproduce_all.sh",
    # Documents that mention the test suite without telling a reader to
    # install it. They are in scope because a future edit could add a
    # command block, and the rule should catch that too.
    "RESULTS.md",
    "REVIEWERS.md",
    "NEXT_STEPS.md",
    "TEAM.md",
    "paper.md",
    "paper/METHODS_PAPER.md",
    "paper/PAPER.md",
    "paper/paper.tex",
    "docs/JOURNAL_SUBMISSION.md",
    "docs/PORTFOLIO_ROADMAP.md",
    "docs/PATH_TO_IMPACT.md",
    "ULTRA_EARLY_READINESS.md",
    "paper/REPRODUCE.md",
)

# An install command. The line has to START with the installer, so a
# message that merely mentions `pip install -e ".[dev]"` (inside an
# `echo`, say) is not an install -- it is advice, and advice that no
# reader executes cannot make a later command work.
_INSTALL = re.compile(
    r"^\s*(?:[\w./-]*/)?(?:pip|python3?\s+-m\s+pip|uv\s+pip)\s+install\b"
)
# ...and one that provides pytest.
_PROVIDES_PYTEST = re.compile(r"\bpytest\b|\.\[?[^\]\s]*dev[^\]\s]*\]?")
# A command that runs pytest: `pytest -q`, `python -m pytest ...`,
# `./.venv/bin/python -m pytest ...`. The token has to be in a command
# position -- after `-m`, or at the start of the line -- so that prose
# and a step table (`pytest_tests|45|pytest test/...`) are not mistaken
# for commands. A requirement (`pytest >= 7.0.0`) is not a command
# either.
_RUNS_PYTEST = re.compile(r"(?:\s-m\s+pytest\b)|(?:^\s*pytest\b(?!\s*[<>=~!]))")
# A check that pytest is importable.
_CHECKS_PYTEST = re.compile(r"import\s+pytest\b")

_HTML_PRE = re.compile(r"<pre[^>]*>(.*?)</pre>", re.S)
_FENCED = re.compile(r"^```([^\n]*)\n(.*?)^```", re.S | re.M)

# A fenced block whose language says it holds shell commands. A ```text
# block holds a package list, not a command: `pytest >= 7.0.0` in a
# dependency table is a requirement, not something a reader runs.
_SHELL_LANGUAGES = {
    "", "bash", "sh", "shell", "zsh", "console", "shell-session", "shellsession",
}

# The phrase that ties a package list to the base install.
_AUTO_INSTALLED = "installed automatically"


# --------------------------------------------------------------------------
# pyproject.toml
# --------------------------------------------------------------------------

def _pyproject() -> str:
    return (_REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")


def _section(text: str, header: str) -> str:
    """The body of a TOML table, up to the next table header."""
    match = re.search(
        r"^\[" + re.escape(header) + r"\]\s*$(.*?)(?=^\[|\Z)",
        text,
        re.S | re.M,
    )
    return match.group(1) if match else ""


def _string_list(body: str, key: str) -> list[str]:
    """The strings of a TOML array of strings, e.g. `dependencies = [...]`."""
    match = re.search(
        r"^" + re.escape(key) + r"\s*=\s*\[(.*?)\]", body, re.S | re.M
    )
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def _runtime_dependencies() -> list[str]:
    """The runtime `dependencies`, as declared."""
    return _string_list(_section(_pyproject(), "project"), "dependencies")


def _dev_extra() -> list[str]:
    """The `dev` extra, as declared."""
    body = _section(_pyproject(), "project.optional-dependencies")
    return _string_list(body, "dev")


def _package_name(requirement: str) -> str:
    """The distribution name of a requirement string."""
    return re.split(r"[<>=!~;\[\s]", requirement, maxsplit=1)[0].strip().lower()


# --------------------------------------------------------------------------
# documents
# --------------------------------------------------------------------------

def _doc_text(name: str) -> str:
    return (_REPO_ROOT / name).read_text(encoding="utf-8")


def _command_blocks(name: str) -> list[str]:
    """The blocks of a document that hold shell commands."""
    text = _doc_text(name)
    if name.endswith(".html"):
        return _HTML_PRE.findall(text)
    if name.endswith(".sh"):
        # A shell script is one command block, fences and all.
        return [text]
    return [
        body for language, body in _FENCED.findall(text)
        if language.strip().lower() in _SHELL_LANGUAGES
    ]


def _lines_of(name: str) -> list[str]:
    """The document's command lines, with continuation lines joined.

    A command spread over several lines with a trailing backslash is one
    command; joining it keeps the install that precedes it and the
    pytest run that follows it adjacent, which is what the rule is about.

    Comment lines are kept: the rule that checks a claimed outcome reads
    the comment under the command. The rules that look for a COMMAND skip
    them.
    """
    lines: list[str] = []
    buffer = ""
    for raw in "\n".join(_command_blocks(name)).splitlines():
        line = raw.rstrip()
        if buffer:
            buffer += " " + line.strip()
        else:
            buffer = line
        if buffer.endswith("\\"):
            buffer = buffer[:-1].rstrip()
            continue
        lines.append(buffer)
        buffer = ""
    if buffer:
        lines.append(buffer)
    return lines


def _is_comment(line: str) -> bool:
    """Is this line a comment rather than a command?"""
    return line.strip().startswith("#")


def _provides_pytest(line: str) -> bool:
    """Does this install command provide pytest?

    Either it names the `dev` extra (which carries pytest) or it names
    `pytest` as a package to install.
    """
    if not _INSTALL.match(line):
        return False
    return bool(_PROVIDES_PYTEST.search(line))


# --------------------------------------------------------------------------
# rule 1: a documented pytest run has a pytest-providing install before it
# --------------------------------------------------------------------------

def _unprepared_runs(lines: list[str]) -> list[str]:
    """The pytest runs in `lines` that nothing before them prepared.

    A document prepares pytest either by installing something that
    provides it, or by checking that pytest is importable and telling
    the reader what to install.
    """
    offenders: list[str] = []
    prepared = False
    for line in lines:
        if _is_comment(line):
            continue
        if _CHECKS_PYTEST.search(line):
            prepared = True
        if _INSTALL.match(line) and _provides_pytest(line):
            prepared = True
        if _RUNS_PYTEST.search(line) and not _INSTALL.match(line):
            if not prepared:
                offenders.append(line.strip())
            break
    return offenders


def test_every_documented_pytest_run_has_pytest_installed_first():
    """A reader who follows a quickstart must be able to run its tests.

    The rule: before the first `pytest` run in a document, the document
    must either install something that provides pytest, or check that
    pytest is importable and tell the reader what to install.
    """
    offenders: list[str] = []

    for name in _DOCS:
        for line in _unprepared_runs(_lines_of(name)):
            offenders.append(f"{name}: {line}")

    assert not offenders, (
        "a documented command runs pytest without anything installing it:\n  "
        + "\n  ".join(offenders)
        + "\n`pytest` is a development dependency: it is in the `dev` extra, "
          "not in the runtime dependencies. Install `pip install -e \".[dev]\"` "
          "first, or check that pytest is importable and say what to install."
    )


def test_the_rule_catches_the_defect_it_was_written_for():
    """The rule must reject a base install followed by a pytest run.

    Without this, the rule above could pass for the wrong reason -- an
    extractor that found nothing at all would report no offenders.
    """
    assert _unprepared_runs([
        "pip install -e .",
        "pytest -q",
    ]), "the rule accepted the defect it exists to catch"

    assert not _unprepared_runs([
        'pip install -e ".[dev]"',
        "pytest -q",
    ]), "the rule rejected a document that installs the dev extra"

    assert not _unprepared_runs([
        'if ! python -c "import pytest"; then echo "install .[dev]"; fi',
        "python -m pytest test/ -q",
    ]), "the rule rejected a document that checks pytest is importable"


def test_the_rule_finds_a_run_with_no_install():
    """The detector must see a pytest run, or the rule above is vacuous."""
    lines = _lines_of("README.md")
    runs = [line for line in lines
            if _RUNS_PYTEST.search(line) and not _INSTALL.match(line)
            and not _is_comment(line)]
    assert runs, "the extractor found no pytest run in README.md"

    prepared = any(_provides_pytest(line) for line in lines)
    assert prepared, "the extractor found no pytest-providing install"


def test_the_install_detector_ignores_advice_inside_an_echo():
    """A message that mentions an install is not an install.

    `scripts/reproduce_all.sh` prints the install command when the venv
    is missing. Treating that as an install would let the script pass
    this rule while still failing at the pytest step.
    """
    assert not _INSTALL.match('echo "   Create with: pip install -e ."')


def test_the_install_detector_recognises_the_real_forms():
    """The installer forms a reader actually uses must be recognised."""
    for command in (
        "pip install -e .",
        "pip install -e \".[dev]\"",
        "./.venv/bin/pip install -e \".[dev,api]\"",
        "python3 -m pip install -e .",
    ):
        assert _INSTALL.match(command), command

    assert _provides_pytest('pip install -e ".[dev]"')
    assert _provides_pytest('pip install -e ".[dev,api]"')
    assert _provides_pytest("pip install pytest")
    assert not _provides_pytest("pip install -e .")
    assert not _provides_pytest('pip install -e ".[api]"')


# --------------------------------------------------------------------------
# rule 2: the premise in pyproject.toml
# --------------------------------------------------------------------------

def test_pytest_is_a_development_dependency_not_a_runtime_one():
    """The premise of the rule above has to be true of `pyproject.toml`."""
    runtime = [_package_name(req) for req in _runtime_dependencies()]
    dev = [_package_name(req) for req in _dev_extra()]

    assert "pytest" not in runtime, (
        "pytest is a runtime dependency again; the quickstarts that install "
        "the base distribution would work, and this guard's rule would have "
        "to be rewritten."
    )
    assert "pytest" in dev, (
        "the `dev` extra does not carry pytest; a document that installs "
        "only the dev extra could not run the test suite."
    )


# --------------------------------------------------------------------------
# rule 3: "installed automatically" lists only runtime dependencies
# --------------------------------------------------------------------------

def _auto_install_offenders(text: str, runtime: set[str]) -> list[str]:
    """The packages a document claims the base install provides, minus
    the ones that really are runtime dependencies.
    """
    offenders: list[str] = []
    for match in re.finditer(re.escape(_AUTO_INSTALLED), text):
        tail = text[match.end():]
        block = re.search(r"```[a-z]*\n(.*?)```", tail, re.S)
        if not block:
            continue
        for line in block.group(1).splitlines():
            entry = line.split("#")[0].strip()
            if not entry:
                continue
            package = _package_name(entry.split()[0])
            if not package or package in runtime:
                continue
            offenders.append(entry)
    return offenders


def test_a_list_claimed_as_automatically_installed_is_one():
    """A document's automatic-install list must match `pyproject.toml`.

    The wiki listed `torch` and `statsmodels` under "Core Python
    dependencies (installed automatically by `pip install -e .`)". Neither
    is a runtime dependency, so the list was a promise the install did
    not keep.
    """
    runtime = {_package_name(req) for req in _runtime_dependencies()}
    offenders: list[str] = []

    for name in _DOCS:
        for entry in _auto_install_offenders(_doc_text(name), runtime):
            offenders.append(f"{name}: {entry}")

    assert not offenders, (
        "a document says these are installed automatically by "
        "`pip install -e .`, but they are not runtime dependencies:\n  "
        + "\n  ".join(offenders)
        + "\nList them as optional, with the extra that provides each."
    )


def test_the_automatic_install_check_catches_an_extra_package():
    """The check must reject a package the base install does not provide."""
    text = (
        "Core dependencies (installed automatically by `pip install -e .`):\n"
        "```text\nnumpy >= 1.24.0\ntorch >= 2.0.0\n```\n"
    )
    offenders = _auto_install_offenders(text, {"numpy"})
    assert offenders == ["torch >= 2.0.0"], offenders


def test_the_wiki_lists_the_runtime_dependencies_it_promises():
    """The wiki's own list must be the runtime list, not a subset."""
    text = _doc_text("wiki/Getting-Started.md")
    assert _AUTO_INSTALLED in text, "the wiki no longer makes this claim"

    runtime = {_package_name(req) for req in _runtime_dependencies()}
    block = re.search(
        r"```[a-z]*\n(.*?)```", text.split(_AUTO_INSTALLED, 1)[1], re.S
    )
    assert block, "the wiki no longer has the dependency list"
    listed = {
        _package_name(line.split("#")[0].strip().split()[0])
        for line in block.group(1).splitlines()
        if line.split("#")[0].strip()
    }
    missing = runtime - listed
    assert not missing, (
        f"the wiki's automatic-install list omits {sorted(missing)}; they "
        f"are runtime dependencies"
    )


# --------------------------------------------------------------------------
# rule 4: a claimed outcome is the outcome
# --------------------------------------------------------------------------

_CLAIM_SPLIT = re.compile(
    r"(\d+)\s+passed,\s*(\d+)\s+skipped,\s*(\d+)\s+deselected"
)
_CLAIM_PASSED = re.compile(r"(\d+)\s+passed\b")
_CLAIM_COLLECTED = re.compile(
    r"(\d+)\s+collected\s*\(\s*(\d+)\s+run,\s*(\d+)\s+deselected"
)

_COLLECTED = re.compile(r"(\d+)(?:/(\d+))? tests? collected")


def _collect(targets: list[str], marker: str | None) -> tuple[int, int]:
    """(total, selected) for a pytest collection over `targets`.

    `total` is every test under the targets; `selected` is the number
    the marker keeps. Without a marker they are the same number.
    """
    def run(extra: list[str]) -> tuple[int, int]:
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", *targets, "--collect-only", "-q",
             "-o", "addopts=", *extra],
            cwd=_REPO_ROOT, capture_output=True, text=True, timeout=600,
        )
        if proc.returncode != 0:
            pytest.skip(f"collection failed here: {proc.stderr[-300:]}")
        for line in proc.stdout.splitlines():
            match = _COLLECTED.search(line)
            if match:
                first = int(match.group(1))
                second = int(match.group(2)) if match.group(2) else None
                if second is None:
                    return first, first
                return second, first
        pytest.skip(f"could not read a collected count:\n{proc.stdout[-300:]}")

    total, _ = run([])
    if marker is None:
        return total, total
    _, selected = run(["-m", marker])
    return total, selected


def _targets_of(command: str) -> list[str]:
    """The repository paths a pytest command names."""
    out = []
    for token in command.split():
        token = token.strip("\"'")
        if not token or token.startswith("-"):
            continue
        if "/" in token or token.endswith(".py"):
            if (_REPO_ROOT / token).exists():
                out.append(token)
    return out


def _marker_of(command: str) -> str | None:
    """The `-m` marker a pytest command carries."""
    match = re.search(r'-m\s+"([^"]+)"', command)
    return match.group(1) if match else None


def test_a_documented_test_command_outcome_is_the_outcome():
    """A claim like "# -> 736 passed" must be what the command produces.

    The README claimed "702 passed, 19 skipped, 13 deselected" for the
    fast-path gate. The gate selects 755 of 768 collected tests, so the
    claim was 34 short and nothing noticed: the count guard compares the
    documents with each other, not with the suite.

    The check is on the total. The split between passed and skipped needs
    a full run (minutes), so a claim that adds up is accepted here even
    if its two halves are swapped; a claim that does not add up is not.
    """
    offenders: list[str] = []

    for name in ("README.md", "README.zh-CN.md"):
        lines = _lines_of(name)
        for index, line in enumerate(lines):
            if _is_comment(line) or not _RUNS_PYTEST.search(line):
                continue
            # The claim is on the next line, as a trailing comment.
            if index + 1 >= len(lines):
                continue
            claim = lines[index + 1]
            if "→" not in claim and "->" not in claim:
                continue
            targets = _targets_of(line)
            if not targets:
                continue
            total, selected = _collect(targets, _marker_of(line))

            split = _CLAIM_SPLIT.search(claim)
            collected = _CLAIM_COLLECTED.search(claim)
            passed = _CLAIM_PASSED.search(claim)
            if split:
                claim_passed, claim_skipped, claim_deselected = (
                    int(split.group(1)), int(split.group(2)),
                    int(split.group(3)),
                )
                if claim_passed + claim_skipped != selected:
                    offenders.append(
                        f"{name}: {claim.strip()} -- the command selects "
                        f"{selected} of {total} collected tests"
                    )
                elif claim_deselected != total - selected:
                    offenders.append(
                        f"{name}: {claim.strip()} -- the command deselects "
                        f"{total - selected} of {total} collected tests"
                    )
            elif collected:
                claim_total, claim_run, claim_deselected = (
                    int(collected.group(1)), int(collected.group(2)),
                    int(collected.group(3)),
                )
                if claim_total != total:
                    offenders.append(
                        f"{name}: {claim.strip()} -- the command collects "
                        f"{total}"
                    )
                elif claim_run != selected or claim_deselected != total - selected:
                    offenders.append(
                        f"{name}: {claim.strip()} -- the command runs "
                        f"{selected} and deselects {total - selected} of "
                        f"{total} collected tests"
                    )
            elif passed:
                if int(passed.group(1)) != selected:
                    offenders.append(
                        f"{name}: {claim.strip()} -- the command selects "
                        f"{selected} of {total} collected tests"
                    )

    assert not offenders, (
        "a document claims an outcome its command does not produce:\n  "
        + "\n  ".join(offenders)
        + "\nRun the command and use the numbers it prints."
    )


def test_the_claim_parser_reads_the_shapes_the_readme_uses():
    """The parser must read the claim shapes, or the rule above is vacuous."""
    assert _CLAIM_SPLIT.search(
        "# → 736 passed, 19 skipped, 13 deselected (slow real-data smoke)"
    )
    assert _CLAIM_PASSED.search("# → 16 passed")
    assert _CLAIM_COLLECTED.search("# → 228 collected (217 run, 11 deselected as slow)")


# --------------------------------------------------------------------------
# rule 5: the two READMEs state the same numbers
# --------------------------------------------------------------------------

_EN_TABLE = "**Test coverage by suite**"
_ZH_TABLE = "**按套件划分的测试覆盖**"


def _table_after(text: str, marker: str) -> list[str]:
    """The markdown table that follows `marker`."""
    index = text.find(marker)
    assert index >= 0, f"the document no longer has {marker!r}"
    rows: list[str] = []
    for line in text[index:].splitlines():
        stripped = line.strip()
        if stripped.startswith("|"):
            rows.append(stripped)
        elif rows:
            break
    return rows


def _numbers_in(rows: list[str]) -> list[int]:
    """The numbers a table states, in order."""
    out: list[int] = []
    for row in rows:
        out.extend(int(token) for token in re.findall(r"\d+", row))
    return sorted(out)


def test_the_two_readmes_state_the_same_test_numbers():
    """The Chinese README's test table must state the English one's numbers.

    The two READMEs are translations of each other and had drifted apart:
    the English badge said 953 while the Chinese one said 919, and the
    Chinese table said the gate collects 734 where the English one said
    768. The count guard compares the documents with each other, so it
    could not see a drift inside one of them -- and the Chinese README
    was not in its list at all.

    Only the numbers are compared: the labels are translated, and a
    translation is not a defect.
    """
    english = _table_after(_doc_text("README.md"), _EN_TABLE)
    chinese = _table_after(_doc_text("README.zh-CN.md"), _ZH_TABLE)

    assert english, "the English README no longer has the test table"
    assert chinese, "the Chinese README no longer has the test table"

    assert _numbers_in(english) == _numbers_in(chinese), (
        "the two READMEs state different test numbers:\n"
        f"  README.md:        {_numbers_in(english)}\n"
        f"  README.zh-CN.md:  {_numbers_in(chinese)}"
    )


def test_the_table_reader_finds_a_table():
    """The reader must find the table, or the rule above is vacuous."""
    rows = _table_after(_doc_text("README.md"), _EN_TABLE)
    assert len(rows) > 3, rows
    assert _numbers_in(rows), rows
