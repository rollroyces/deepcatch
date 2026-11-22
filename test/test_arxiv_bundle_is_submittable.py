"""The arXiv source bundle must be a bundle arXiv can actually process.

Why this file exists
--------------------
arXiv does not compile a submission with the author's toolchain: it
unpacks the archive and compiles it with its own TeX Live, with no
shell escape and no external programs. A source bundle that works
locally can still come back as "no output produced" -- the submission
never reaches the queue, and the author learns it from an error screen
rather than from a check.

The bundle this repository builds is deliberately plain: `article`,
`hyperref`, `geometry`, `amsmath`, `booktabs`, `graphicx`, and the
bibliography written inline. Nothing in it needs anything arXiv does not
have. That is a property worth keeping -- it is exactly the kind of
thing that decays silently, when a manuscript later grows a package
(`minted` needs shell escape, `fontspec` needs XeLaTeX, `biblatex` needs
biber) and the change is invisible locally.

What this test checks
---------------------
1. The manuscript does not use a package arXiv cannot run.
2. It does not ask for shell escape.
3. It contains no absolute path -- arXiv's machine has none of ours.
4. Every figure it includes is a member of the bundle.
5. The bundle holds no build leftovers (`.aux`, `.log`, ...) and no file
   exceeds arXiv's 10 MB limit.
6. If it uses a bibliography file, the compiled `.bbl` is in the bundle,
   because arXiv does not run BibTeX or biber.
7. Its document class is one arXiv accepts.

What it does not check
----------------------
Whether the manuscript compiles: `scripts/make_arxiv_bundle.py --verify`
does that, and `test_arxiv_bundle_is_current.py` checks the bundle
matches the committed manuscript. This file checks the rules arXiv
applies on top of "it compiles".

Scope
-----
The bundle is a build product, not a tracked file. Where there is no
local bundle -- CI's core job, a fresh clone -- there is nothing to
check and the test says so instead of failing; the paper workflow builds
one from the committed manuscript and publishes it as the
`deepcatch-arxiv-source` artifact.
"""
from __future__ import annotations

import re
import tarfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_BUNDLE = _REPO_ROOT / "dist" / "deepcatch_v2_3_0.tar.gz"

# arXiv's limit: 10 MB for the output PDF or for any individual file in
# the submission.
_MAX_FILE_BYTES = 10 * 1024 * 1024

# Packages arXiv's TeX Live cannot run, and why.
#
#   minted, shadings  -- need \write18 (shell escape), which arXiv refuses
#   fontspec, unicode-math, xltxtra
#                      -- need XeTeX or LuaTeX, not pdfLaTeX
#   svg               -- needs Inkscape, which arXiv does not run
#   gv                -- needs Ghostscript's ps2eps wrapper
#   fontawesome       -- not in arXiv's TeX Live
#   biblatex          -- needs biber; handled separately, because biblatex
#                        with a .bbl in the bundle is fine
_FORBIDDEN = {
    "minted": "needs \\write18 (shell escape), which arXiv refuses",
    "shadings": "is minted's shading backend",
    "fontspec": "needs XeTeX or LuaTeX, not pdfLaTeX",
    "unicode-math": "needs XeTeX or LuaTeX, not pdfLaTeX",
    "xltxtra": "needs XeTeX or LuaTeX, not pdfLaTeX",
    "svg": "needs Inkscape, which arXiv does not run",
    "gv": "needs Ghostscript's ps2eps wrapper",
    "fontawesome": "is not in arXiv's TeX Live",
}

# arXiv's own accepted document classes.
_ACCEPTED_CLASSES = {"article", "revtex", "elsarticle"}

# Files that must not travel with a submission.
_BUILD_LEFTOVERS = {
    ".aux", ".log", ".out", ".toc", ".synctex.gz", ".ds_store", ".bbl.tmp",
}


def _members() -> dict[str, bytes]:
    """The bundle's members, name -> content."""
    if not _BUNDLE.exists():
        pytest.skip(
            "no local arXiv bundle; CI builds one from the committed "
            "manuscript and publishes it as `deepcatch-arxiv-source`")
    with tarfile.open(_BUNDLE, "r:gz") as tar:
        out = {}
        for member in tar.getmembers():
            if not member.isfile():
                continue
            handle = tar.extractfile(member)
            assert handle is not None, f"{member.name} is not readable"
            out[member.name] = handle.read()
        return out


def _manuscript(members: dict[str, bytes]) -> str:
    """The bundle's manuscript source."""
    names = [name for name in members if name.endswith(".tex")]
    assert len(names) == 1, (
        f"the bundle must hold exactly one manuscript, found {names}")
    return members[names[0]].decode("utf-8")


def _packages(manuscript: str) -> list[str]:
    """The packages the manuscript loads."""
    found: list[str] = []
    for group in re.findall(r"\\usepackage(?:\[[^\]]*\])?\{([^}]*)\}", manuscript):
        found.extend(name.strip() for name in group.split(","))
    return [name for name in found if name]


def _forbidden(packages: list[str]) -> list[str]:
    """The packages arXiv cannot process, with the reason."""
    return [
        f"{name} ({_FORBIDDEN[name]})"
        for name in packages
        if name in _FORBIDDEN
    ]


def test_the_bundle_uses_no_package_arxiv_cannot_run():
    """Every package must be one arXiv's own TeX Live can load."""
    manuscript = _manuscript(_members())
    offenders = _forbidden(_packages(manuscript))
    assert not offenders, (
        "the arXiv bundle uses packages arXiv cannot process:\n  "
        + "\n  ".join(offenders)
        + "\narXiv compiles the submission with its own TeX Live, without "
          "shell escape. The bundle compiles locally, so this only shows up "
          "after the upload."
    )


def test_the_package_check_catches_a_package_arxiv_cannot_run():
    """The check must reject what arXiv rejects, and nothing else."""
    assert _forbidden(["hyperref", "geometry", "booktabs"]) == []
    assert _forbidden(["minted"]) == [
        "minted (needs \\write18 (shell escape), which arXiv refuses)"]


def test_the_package_reader_reads_a_group_and_an_option():
    """The reader must see a package however the manuscript writes it."""
    manuscript = (
        "\\documentclass[11pt]{article}\n"
        "\\usepackage[margin=1in]{geometry}\n"
        "\\usepackage{amsmath, amssymb}\n"
    )
    assert _packages(manuscript) == ["geometry", "amsmath", "amssymb"]


def test_the_bundle_asks_for_no_shell_escape():
    """`\\write18` is how a package asks for shell escape; arXiv refuses."""
    manuscript = _manuscript(_members())
    assert "\\write18" not in manuscript, (
        "the arXiv bundle asks for shell escape (\\write18); arXiv compiles "
        "submissions without it."
    )


def test_the_bundle_contains_no_absolute_path():
    """arXiv's machine has none of our paths."""
    manuscript = _manuscript(_members())
    absolute = re.findall(r"(?<![\w.])/(?:Users|home|tmp|opt|var)/[\w./-]*",
                           manuscript)
    assert not absolute, (
        "the arXiv bundle contains absolute paths: "
        + ", ".join(sorted(set(absolute)))
        + ". arXiv unpacks the archive on its own machine; a path that "
          "exists here does not exist there."
    )


def test_every_figure_the_manuscript_includes_is_in_the_bundle():
    """A figure the manuscript names but the bundle lacks is a dead build."""
    members = _members()
    manuscript = _manuscript(members)
    wanted = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", manuscript)
    assert wanted, "the manuscript includes no figure; the check is vacuous"
    missing = [name for name in wanted if name not in members]
    assert not missing, (
        "the arXiv bundle does not contain the figures the manuscript "
        f"includes: {missing}. The bundle's builder flattens the paths; a "
        "figure it does not flatten is a compile error on arXiv."
    )


def test_the_bundle_carries_no_build_leftovers():
    """arXiv asks for the sources, not the compiler's scratch files."""
    members = _members()
    leftovers = [
        name for name in members
        if any(name.lower().endswith(suffix) for suffix in _BUILD_LEFTOVERS)
    ]
    assert not leftovers, (
        "the arXiv bundle carries build leftovers: "
        + ", ".join(sorted(leftovers))
        + ". arXiv asks for the sources and the figures, not the "
          "compiler's intermediate files."
    )


def _sizes() -> list[tuple[str, int]]:
    """The bundle's members, name -> size in bytes."""
    if not _BUNDLE.exists():
        pytest.skip(
            "no local arXiv bundle; CI builds one from the committed "
            "manuscript and publishes it as `deepcatch-arxiv-source`")
    with tarfile.open(_BUNDLE, "r:gz") as tar:
        return [(member.name, member.size) for member in tar.getmembers()]


def test_no_file_in_the_bundle_exceeds_arxivs_limit():
    """arXiv's limit is 10 MB for the output PDF or for any single file."""
    too_big = [
        f"{name} ({size} bytes)"
        for name, size in _sizes()
        if size > _MAX_FILE_BYTES
    ]
    assert not too_big, (
        "the arXiv bundle holds a file above arXiv's 10 MB limit:\n  "
        + "\n  ".join(too_big)
        + "\narXiv warns above 10 MB for the output PDF or for any "
          "individual file, and rejects the upload."
    )


def test_the_bibliography_is_something_arxiv_can_use():
    """arXiv does not run BibTeX or biber.

    A manuscript that reads its references from a `.bib` needs the
    compiled `.bbl` in the bundle; one that writes them inline (as this
    manuscript does) needs nothing.
    """
    members = _members()
    manuscript = _manuscript(members)

    uses_bibtex = bool(re.search(r"\\(?:bibliography|addbibresource)\{", manuscript))
    if not uses_bibtex:
        return

    has_bbl = any(name.endswith(".bbl") for name in members)
    has_bib = any(name.endswith(".bib") for name in members)
    assert has_bbl, (
        "the arXiv bundle's manuscript reads its references from a "
        "bibliography file but the bundle carries no compiled .bbl; arXiv "
        "does not run BibTeX, so the references would not appear"
    )
    assert has_bib, (
        "the arXiv bundle's manuscript reads its references from a "
        "bibliography file but the bundle carries no .bib either"
    )


def test_the_document_class_is_one_arxiv_accepts():
    """arXiv accepts article, revtex and elsarticle."""
    manuscript = _manuscript(_members())
    match = re.search(r"\\documentclass(?:\[[^\]]*\])?\{([^}]*)\}", manuscript)
    assert match, "the manuscript has no \\documentclass"
    name = match.group(1).strip()
    assert name in _ACCEPTED_CLASSES, (
        f"the arXiv bundle's document class is {name!r}; arXiv accepts "
        f"{sorted(_ACCEPTED_CLASSES)}"
    )
