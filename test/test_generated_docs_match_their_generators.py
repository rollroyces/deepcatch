"""A committed generated document must still be what its generator produces.

What went wrong
---------------
`docs/PUBLICATION_READINESS.md` is produced by
`scripts/_verify_publication_readiness.py`, and for a while its "Honest
reading" said the held-out publications were "**BLOCKED** ... the API/S3
are not reachable to fetch them" -- while the same document's API verdict
read `ok`, and while `results/external_validation_finaledb.json` recorded
those publications as scored. The correcting note had been hand-edited
into a *generated* file, so the next run of the generator would have
deleted it.

Nothing checked the relationship between the two: the document could say
anything, and the only guard on its numbers compared it against the
artifacts it quoted -- which the generator writes itself, so a hand edit
that changed a number would have been caught only by luck.

This test closes that loop for every generator that can run offline from
the repository's own committed artifacts: it regenerates the document and
the artifact into a temporary directory and asserts that both committed
copies are exactly what the generator produces. A hand edit, a stale
document, an artifact regenerated without its document (or the reverse)
all fail here.

The generators that need the sibling repository's feature cache or its
labels file cannot run in CI; their documents are guarded by the number
and cohort checks in `test_paper_numbers_match_artifacts.py` and
`test_publication_readiness.py` instead.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Generators that run offline from this repository's committed artifacts.
OFFLINE_GENERATORS = [
    {
        "name": "ultraearly_readiness",
        "script": "scripts/ultraearly_readiness.py",
        "doc": "docs/ULTRA_EARLY_READINESS_AUTO.md",
        "artifact": "results/ultraearly_readiness.json",
        "args": ["--results-dir", "results"],
    },
    {
        # The paste-ready bioRxiv abstract. It was a separately written
        # document for months, and it carried an older abstract, an older
        # title and the author line "Royce et al." on a single-author
        # paper -- so the text the author submits and the PDF they upload
        # were different papers. Now the manuscript's LaTeX is the only
        # source, and this generator is what renders it.
        "name": "preprint_abstract",
        "script": "scripts/generate_preprint_abstract.py",
        "doc": "docs/PREPRINT_ABSTRACT.md",
        "artifact": None,
        "out_flag": "--out",
    },
]

# The generators that also produce a JSON artifact. A generator that only
# writes a document has no artifact to compare.
GENERATORS_WITH_ARTIFACTS = tuple(
    spec for spec in OFFLINE_GENERATORS if spec.get("artifact"))

# An ISO-8601 instant, however many fractional digits it carries.
_TIMESTAMP = re.compile(
    r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:[+-]\d{2}:\d{2}|Z)?"
)
# The generator stamps the commit it ran at; that moves with every commit
# and says nothing about the content, so it is normalised away.
_COMMIT_LINE = re.compile(r"^- \*\*Repo commit:\*\* .*$", re.MULTILINE)


def _normalise_markdown(text: str) -> str:
    text = _TIMESTAMP.sub("<TS>", text)
    return _COMMIT_LINE.sub("- **Repo commit:** <COMMIT>", text)


def _normalise_json(value):
    """Replace the generation-time fields wherever they appear."""
    if isinstance(value, dict):
        out = {}
        for key, item in value.items():
            if key in ("generated_at", "repo_commit", "checked_at"):
                out[key] = "<TS>"
            else:
                out[key] = _normalise_json(item)
        return out
    if isinstance(value, list):
        return [_normalise_json(item) for item in value]
    if isinstance(value, str):
        return _TIMESTAMP.sub("<TS>", value)
    return value


def _run_generator(spec: dict, tmp_path: Path) -> tuple[Path, Path]:
    out_json = tmp_path / "regenerated.json"
    out_md = tmp_path / "regenerated.md"
    cmd = [
        sys.executable, "-u", spec["script"],
        *spec.get("args", []),
        spec.get("out_flag", "--out-md"), str(out_md),
    ]
    if spec.get("artifact"):
        cmd += ["--out-json", str(out_json)]
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          cwd=str(_REPO_ROOT), timeout=600)
    assert proc.returncode == 0, (
        f"{spec['script']} failed (rc={proc.returncode}):\n"
        f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
    )
    return out_md, out_json


@pytest.mark.parametrize("spec", OFFLINE_GENERATORS, ids=lambda s: s["name"])
def test_the_committed_document_is_what_the_generator_produces(spec, tmp_path):
    """The committed markdown must equal a fresh generation."""
    out_md, _ = _run_generator(spec, tmp_path)
    committed = _REPO_ROOT / spec["doc"]
    assert committed.is_file(), f"{spec['doc']} is missing"

    expected = _normalise_markdown(out_md.read_text(encoding="utf-8"))
    actual = _normalise_markdown(committed.read_text(encoding="utf-8"))
    assert actual == expected, (
        f"{spec['doc']} has drifted from {spec['script']}. Regenerate it "
        f"(the generator writes it from the committed artifacts) or fix the "
        f"generator -- do not hand-edit the generated file, because the next "
        f"generation will erase the edit."
    )


@pytest.mark.parametrize("spec", GENERATORS_WITH_ARTIFACTS, ids=lambda s: s["name"])
def test_the_committed_artifact_is_what_the_generator_produces(spec, tmp_path):
    """The committed JSON artifact must equal a fresh generation too."""
    _, out_json = _run_generator(spec, tmp_path)
    committed = _REPO_ROOT / spec["artifact"]
    assert committed.is_file(), f"{spec['artifact']} is missing"

    expected = _normalise_json(json.loads(out_json.read_text(encoding="utf-8")))
    actual = _normalise_json(json.loads(committed.read_text(encoding="utf-8")))
    assert actual == expected, (
        f"{spec['artifact']} has drifted from {spec['script']}: the artifact "
        f"and the document it describes must both come from the same run."
    )
