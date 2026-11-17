"""Where this repository, its sibling pipeline and its caches live.

Why this file exists
--------------------
The scripts that need the FinaleDB feature cache, the multi-class
labels, the pre-trained checkpoints or the FinaleMe installation used
to hardcode one machine's home directory. Those paths do not exist for
a reader, and several of them had stopped existing for the author too
-- the repository moved, so the stale ones were broken for everyone.

The resolution is the same everywhere:

1. the environment variable, when it is set;
2. the sibling repository cloned NEXT TO this repository, which is where
   the README, `paper/REPRODUCE.md` and `RESULTS.md` tell a reader to
   put it;
3. `~/cfdna-fragmentomics-pipeline`, the layout this repository was
   developed in.

Everything is resolved from this file's own location, so the repository
can live anywhere.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def _first_existing(*candidates: Path) -> Path:
    """The first candidate that exists, or the last one."""
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def _sibling(repo_name: str, env_var: str) -> Path:
    override = os.environ.get(env_var)
    if override:
        return Path(override).expanduser()
    return _first_existing(
        REPO_ROOT.parent / repo_name,
        Path.home() / repo_name,
    )


# The pipeline repository: the feature cache and the labels live in it.
SIBLING_PIPELINE = _sibling(
    "cfdna-fragmentomics-pipeline", "DEEPCATCH_SISTER_PIPELINE")
FEATURES_DIR = Path(
    os.environ.get("DEEPCATCH_FEATURES_DIR")
    or SIBLING_PIPELINE / "data" / "features")
LABELS_TSV = Path(
    os.environ.get("DEEPCATCH_LABELS")
    or SIBLING_PIPELINE / "labels_multiclass.tsv")

# The methylation repository: the cross-platform fusion inputs live in it.
METHYLATION_PIPELINE = _sibling(
    "deepcatch-methylation", "DEEPCATCH_METHYLATION_PIPELINE")

# FinaleMe, the tool that decodes methylation from fragment files. It is
# an installation, not a repository, so it is looked up in the home
# directory the tool was installed into.
FINALEME_HOME = Path(
    os.environ.get("FINALEME_HOME")
    or Path.home() / ".hermes" / ".local" / "finaleme")
