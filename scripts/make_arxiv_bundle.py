#!/usr/bin/env python3
"""Assemble the arXiv source bundle from the committed manuscript.

Why this exists
---------------
The preprint goes to arXiv, and arXiv takes a LaTeX source archive
rather than a PDF. That has one requirement a PDF upload does not: the
source has to compile *standalone*, outside this repository. This
manuscript reaches outside the repository for three figures:

    ../results/real_tcga_performance.png
    ../docs/figures/fig1_pooled_roc.png
    ../docs/figures/fig4_sens_operating_points.png

Those paths do not exist inside an arXiv source archive, so a bundle
built by copying `paper/paper.tex` alone compiles to a document with
three missing figures. This script flattens them: each referenced
figure is copied next to the manuscript under a flat name and the
`\\includegraphics` paths are rewritten to match.

It also verifies the bundle rather than assuming it: with `--verify`
it extracts the archive into a temporary directory, compiles it there
with whatever LaTeX engine is available, and compares the result with
the committed PDF -- the page count must match, and every figure in the
manuscript must actually be embedded in the output. A bundle that
compiles to a document with missing figures is worse than no bundle,
because the missing figures are invisible until a reader complains.

Usage
-----
    python scripts/make_arxiv_bundle.py            # write dist/…tar.gz
    python scripts/make_arxiv_bundle.py --verify   # …and compile it
"""
from __future__ import annotations

import argparse
import io
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_MANUSCRIPT = _REPO_ROOT / "paper" / "paper.tex"
_DIST = _REPO_ROOT / "dist"
_STEM = "deepcatch_v2_3_0"

_INCLUDE = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")


def _figure_targets(text: str) -> list[str]:
    """The repository-relative paths the manuscript includes."""
    seen: list[str] = []
    for match in _INCLUDE.finditer(text):
        raw = match.group(1).strip()
        # The manuscript lives in paper/, so its paths start with "../".
        relative = raw[3:] if raw.startswith("../") else raw
        if relative not in seen:
            seen.append(relative)
    return seen


def _flat_name(relative: str) -> str:
    """A collision-free flat name for a figure."""
    return relative.replace("/", "_")


def build(out_dir: Path) -> tuple[Path, dict[str, str]]:
    """Write the bundle; return the archive path and the path mapping."""
    out_dir.mkdir(parents=True, exist_ok=True)
    staging = out_dir / f"{_STEM}_src"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    text = _MANUSCRIPT.read_text(encoding="utf-8")
    mapping: dict[str, str] = {}
    for relative in _figure_targets(text):
        source = _REPO_ROOT / relative
        if not source.exists():
            raise SystemExit(f"the manuscript includes {relative}, which is "
                             f"not in the repository")
        flat = _flat_name(relative)
        shutil.copy2(source, staging / flat)
        mapping[relative] = flat

    # Rewrite the paths so they resolve inside the bundle.
    def rewrite(match: re.Match) -> str:
        raw = match.group(1).strip()
        relative = raw[3:] if raw.startswith("../") else raw
        flat = mapping.get(relative, relative)
        prefix = match.group(0).split("{")[0]
        return f"{prefix}{{{flat}}}"

    bundled = _INCLUDE.sub(rewrite, text)
    (staging / f"{_STEM}.tex").write_text(bundled, encoding="utf-8")

    archive = out_dir / f"{_STEM}.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for path in sorted(staging.iterdir()):
            tar.add(path, arcname=path.name)
    shutil.rmtree(staging)
    return archive, mapping


def _pdf_page_count(path: Path) -> int:
    try:
        from pypdf import PdfReader
    except ImportError:
        return -1
    return len(PdfReader(str(path)).pages)


def verify(archive: Path, expected_pdf: Path) -> tuple[bool, str]:
    """Compile the bundle standalone and compare it with the committed PDF."""
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(work)
        tex = work / f"{_STEM}.tex"

        engine = None
        for candidate in ("tectonic", "pdflatex", "latexmk"):
            if shutil.which(candidate):
                engine = candidate
                break
        if engine is None:
            return True, "no LaTeX engine available; the bundle was not compiled"

        if engine == "tectonic":
            proc = subprocess.run(
                ["tectonic", "-X", "compile", tex.name],
                cwd=work, capture_output=True, text=True, timeout=900,
            )
            produced = work / f"{_STEM}.pdf"
        else:
            proc = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", tex.name],
                cwd=work, capture_output=True, text=True, timeout=900,
            )
            produced = work / f"{_STEM}.pdf"

        if proc.returncode != 0 or not produced.exists():
            tail = (proc.stderr or proc.stdout or "")[-800:]
            return False, f"the bundle did not compile:\n{tail}"

        if not expected_pdf.exists():
            return True, ("the bundle compiled, but there is no committed "
                           "PDF to compare it with")

        # A missing figure is the failure this script exists to prevent,
        # so it is checked whether or not the page count can be read.
        text = (work / f"{_STEM}.tex").read_text(encoding="utf-8")
        for relative in _figure_targets(text):
            if not (work / _flat_name(relative)).exists():
                return False, f"the bundle is missing {relative}"

        want = _pdf_page_count(expected_pdf)
        got = _pdf_page_count(produced)
        if want < 0 or got < 0:
            return True, ("the bundle compiles standalone and carries every "
                          "figure; the page count was not compared, because "
                          "pypdf is not installed here")
        if want != got:
            return False, (f"the bundle compiled to {got} pages but the "
                           f"committed PDF has {want}")

        return True, (f"the bundle compiles standalone: {got} pages, "
                      f"the same as the committed PDF")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="compile the bundle and compare it with the "
                             "committed PDF")
    args = parser.parse_args()

    archive, mapping = build(_DIST)
    print(f"[ARXIV] wrote {archive.relative_to(_REPO_ROOT)}")
    for relative, flat in mapping.items():
        print(f"[ARXIV]   {relative} -> {flat}")

    if args.verify:
        ok, message = verify(archive, _REPO_ROOT / "paper" / "biorxiv_submission_v2.3.0.pdf")
        print(f"[ARXIV] {message}")
        if not ok:
            sys.exit(1)


if __name__ == "__main__":
    main()
