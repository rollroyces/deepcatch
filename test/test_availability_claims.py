"""What a reader is told exists must exist for them — in the manuscript
and in the reproduction driver.

What went wrong
---------------
The long-form methods manuscript's Availability section listed two
pre-trained checkpoints under "**Pretrained checkpoints:**" as though
they shipped with the repository:

    - `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt`
    - `checkpoints/foundation_pretrained_SYNTHETIC_v0.pt`

Neither file is in the repository, and neither ever was: there is no
`checkpoints/` directory at all, and the only model file the repository
carries is `model/motif_model_checkpoint.pt`. A reviewer who follows
the Availability section finds nothing there -- the same defect as a
citation that resolves to nothing, and the same defect as the Zenodo
placeholder DOI this repository already removed.

The honest statement is that the checkpoints are *regenerated* rather
than distributed, and that no reported number depends on them: the
encoder behind the §3.1 fusion comparison is trained in-run with
`pretrained=False`, so §3 reproduces from a clean checkout. This test
pins that structure, so a future edit cannot quietly turn the
Availability section back into a list of files that are not there.

The same lesson applies to the reproduction driver, which is what a
reviewer runs first: its default paths pointed at one machine's home
directory, and the message it printed when the feature cache was
missing did not say that the cache is in neither repository and has to
be built from the public API first. Those cases are checked here too.

What it does NOT do
-------------------
It does not decide whether the checkpoints *should* be committed. That
is a repository-size judgement for the author. It only requires that a
path the manuscript names either exists, or sits in a paragraph that
says it is not distributed.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_MANUSCRIPT = _REPO_ROOT / "paper" / "METHODS_PAPER.md"

# A repository-relative path to a model file, in backticks.
_MODEL_PATH = re.compile(r"`([A-Za-z0-9_./-]+\.pt)`")

# The figures the manuscript pulls in from the repository.
_INCLUDE = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")

# The phrase that makes a missing file an honest statement rather than a
# broken promise.
_NOT_DISTRIBUTED = "not distributed"


def _tracked(*patterns: str) -> list[str]:
    """The tracked files matching the given git pathspecs."""
    proc = subprocess.run(
        ["git", "-C", str(_REPO_ROOT), "ls-files", *patterns],
        capture_output=True, text=True, timeout=60, check=True,
    )
    return [name for name in proc.stdout.split("\n") if name]


def _paragraphs(text: str) -> list[str]:
    """The document's paragraphs, blank-line separated."""
    return [p for p in re.split(r"\n\s*\n", text) if p.strip()]


def test_a_named_checkpoint_either_exists_or_is_declared_not_distributed():
    """No path may be offered as available unless the file is there."""
    text = _MANUSCRIPT.read_text(encoding="utf-8")
    offenders = []
    for paragraph in _paragraphs(text):
        for path in _MODEL_PATH.findall(paragraph):
            if (_REPO_ROOT / path).exists():
                continue
            if _NOT_DISTRIBUTED in paragraph:
                continue
            offenders.append(path)
    assert not offenders, (
        "the manuscript names a model file the repository does not "
        f"contain, without saying it is not distributed: {offenders}. "
        "Either commit the file, or say in the same paragraph that it is "
        "not distributed and how to regenerate it."
    )


def test_the_way_to_regenerate_the_checkpoint_is_a_script_that_exists():
    """A disclaimer must point at a remedy that is actually there.

    "Not distributed" is only useful if the reader can regenerate the
    file. The remedy is a script in this repository, so its existence is
    checkable -- and if the script is renamed, this fails rather than
    leaving the manuscript pointing at a path that does not exist.
    """
    text = _MANUSCRIPT.read_text(encoding="utf-8")
    scripts = {
        path
        for paragraph in _paragraphs(text)
        for path in re.findall(r"`(scripts/[A-Za-z0-9_./-]+\.py)`", paragraph)
    }
    pretrain = {
        path for path in scripts if "pretrain" in path
    }
    assert pretrain, (
        "the manuscript no longer names a pre-training script, so a "
        "reader told that the checkpoints are not distributed has no "
        "way to regenerate them"
    )
    missing = sorted(p for p in pretrain if not (_REPO_ROOT / p).exists())
    assert not missing, (
        f"the manuscript points at pre-training scripts that do not "
        f"exist: {missing}"
    )


def test_no_tracked_script_points_at_one_machine_home_directory():
    """Every script a reader is told to run must work for a reader.

    Twenty-seven scripts and five shell scripts hardcoded one machine's
    home directory. Those paths do not exist for a reader -- and several
    had stopped existing for the author too, because the repository had
    been moved, so they were broken for everyone.

    The paths now resolve through `scripts/_paths.py`, which tries the
    environment variable, then the layout the documentation tells a
    reader to use (the sibling repository next to this one), then the
    author's own layout (the home directory). This test is the
    portability half of the same lesson the checkpoint guard above
    teaches: a path that exists on the author's machine is not a path a
    reader has.
    """
    offenders = {}
    for name in _tracked("*.py", "*.sh"):
        path = _REPO_ROOT / name
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        hits = []
        for match in re.finditer(
                r"/(?:Users|home)/[A-Za-z0-9_.-]+", text):
            found = match.group(0)
            # A container's home directory (`/home/node` in the Docker
            # image) is part of that image, not of a machine. Only a
            # home directory that exists HERE can be a path that fails
            # for someone else.
            if found.startswith("/home/") and not Path(found).exists():
                continue
            hits.append(found)
        if hits:
            offenders[name] = sorted(set(hits))
    assert not offenders, (
        "a tracked script points at a home directory, which does not "
        f"exist for a reader: {offenders}. Resolve the path from the "
        "repository's own location (see scripts/_paths.py) instead."
    )


def test_no_tracked_document_gives_a_reader_someone_elses_path():
    """A copy-paste command must run on the reader's machine.

    Twenty documents carried commands like
    commands that invoked one machine's virtualenv by absolute path,
    including the collaborator quick-start -- the entry point for the
    wet-lab collaborators this project is recruiting. A reader who pasted
    one got "no such file or directory" and no clue what to do instead.

    Dated records are exempt: the changelog, the audit reports and the
    generated run reports name the machine they were produced on, and
    that is the point of a record.
    """
    exempt = {
        "CHANGELOG.md",
        "AUDIT_REPORT.md",
        "AUDIT_REPORT_2.md",
        "docs/PUBLICATION_READINESS.md",
        "results/REPRODUCE_REPORT.md",
        "results/CROSS_PLATFORM_HG19_REFS_READY.md",
    }
    offenders = {}
    for name in _tracked("*.md", "*.tex", "*.cff", "*.yaml", "*.yml"):
        if name in exempt or name.startswith("review/"):
            continue
        path = _REPO_ROOT / name
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        hits = sorted({
            m.group(0) for m in re.finditer(
                r"/(?:Users|home)/[A-Za-z0-9_.-]+", text)
        })
        if hits:
            offenders[name] = hits
    assert not offenders, (
        "a tracked document hands the reader another machine's path: "
        f"{offenders}. Write the path relative to the repository (or "
        "relative to the reader's home directory) instead."
    )


def test_the_reproduce_driver_says_how_to_build_the_missing_feature_cache():
    """The missing-cache message must name the remedy, not just the knob.

    The extracted FinaleDB features are gitignored in the pipeline
    repository (~600 MB for 627 samples) and are rebuilt from the public
    API. A reader who is told only "set DEEPCATCH_FEATURES_DIR" has to
    guess where the directory comes from.
    """
    driver = (_REPO_ROOT / "paper" / "REPRODUCE.sh").read_text(encoding="utf-8")
    assert "run_cross_study.py" in driver, (
        "the driver's missing-feature message no longer says how to build "
        "the feature cache; the features are gitignored in the pipeline "
        "repository and must be rebuilt from the public FinaleDB API"
    )


def test_the_readme_says_where_the_benchmark_input_comes_from():
    """The front door must say how to get the data the numbers need.

    The README's quick start ran the cross-study benchmark and never
    said where the feature cache comes from -- and the cache is in
    *neither* repository, it is built from the public FinaleDB API. A
    reader who followed the README alone got a path that does not exist
    on their machine, and no indication of the one command that builds
    it. The README also called a virtualenv interpreter that its own
    installation section never created.

    Both READMEs now name the companion repository and the command that
    materialises the cache, and this keeps them saying so.
    """
    for doc in ("README.md", "README.zh-CN.md"):
        text = (_REPO_ROOT / doc).read_text(encoding="utf-8")
        assert "cfdna-fragmentomics-pipeline" in text, (
            f"{doc} does not name the companion repository whose feature "
            "cache every headline number is computed from"
        )
        assert "run_cross_study.py" in text, (
            f"{doc} does not name the command that builds the feature "
            "cache -- it is in neither repository, so without this line "
            "the quick start cannot work"
        )


def test_every_figure_the_manuscript_includes_exists():
    """The manuscript's figures must be in the repository.

    The arXiv source bundle is assembled from the manuscript and the
    figures it includes, with their paths flattened. A manuscript that
    includes a file the repository does not have produces a bundle that
    compiles to a document with a missing figure -- and the missing
    figure is invisible until a reader complains.
    """
    tex = (_REPO_ROOT / "paper" / "paper.tex").read_text(encoding="utf-8")
    missing = []
    for raw in _INCLUDE.findall(tex):
        stripped = raw.strip()
        relative = stripped[3:] if stripped.startswith("../") else stripped
        if not (_REPO_ROOT / relative).exists():
            missing.append(relative)
    assert not missing, (
        "the manuscript includes figures that are not in the repository: "
        f"{missing}. The arXiv source bundle would compile with those "
        "figures missing."
    )


def test_the_repository_really_has_no_checkpoint_directory_to_discover():
    """The guard above must not be vacuous because the files appeared.

    If the checkpoints are ever committed, the paragraph that says they
    are not distributed becomes wrong. This test states the current
    situation explicitly, so that day the Availability section has to be
    rewritten rather than quietly contradicting itself.
    """
    checkpoints = _REPO_ROOT / "checkpoints"
    if not checkpoints.exists():
        pytest.skip("no checkpoints/ directory in this checkout")
    shipped = sorted(p.name for p in checkpoints.glob("*.pt"))
    assert shipped, (
        "checkpoints/ exists but holds no .pt file; the Availability "
        "section's wording must match whatever is actually there"
    )
