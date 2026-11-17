"""The manuscript must not offer a file the repository does not contain.

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
