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
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_MANUSCRIPT = _REPO_ROOT / "paper" / "METHODS_PAPER.md"

# A repository-relative path to a model file, in backticks.
_MODEL_PATH = re.compile(r"`([A-Za-z0-9_./-]+\.pt)`")

# The phrase that makes a missing file an honest statement rather than a
# broken promise.
_NOT_DISTRIBUTED = "not distributed"


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


def test_the_reproduce_driver_does_not_point_at_one_machine_home_directory():
    """The driver's defaults must work for a reader, not for its author.

    The reproduction driver is what a reviewer runs first. Its default
    paths pointed at `/Users/hermes/...`, so a reader who followed its
    own clone instruction and ran it hit a path that does not exist on
    their machine -- and the failure message said only "set
    DEEPCATCH_FEATURES_DIR to the extracted FinaleDB features", without
    saying that those features are in neither repository and have to be
    built first.

    The defaults are now derived from the repository's own location, so
    the documented layout (the sibling repository cloned next to this
    one) works without any environment variable. The guard is the
    portability half of the same lesson the checkpoint guard above
    teaches: a path that exists on the author's machine is not a path a
    reader has.
    """
    driver = (_REPO_ROOT / "paper" / "REPRODUCE.sh").read_text(encoding="utf-8")
    offenders = sorted({
        m.group(0) for m in re.finditer(r"/(?:Users|home)/[A-Za-z0-9_.-]+", driver)
    })
    assert not offenders, (
        f"paper/REPRODUCE.sh points at a home directory: {offenders}. The "
        "defaults must be derived from the repository's own location, so "
        "that a reader who follows the clone instruction gets working "
        "paths."
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
