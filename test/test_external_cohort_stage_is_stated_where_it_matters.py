"""The held-out cohort's stage must be stated wherever its AUCs are read.

Why this file exists
--------------------
The external AUCs (0.9444, 0.708, 0.7385 pooled) are the project's
headline generalisation number, and the documents that quote them did
not say what stage the cohort was at. A reader would reasonably take
0.9444 for an early-detection result. It is not one: where the public
records characterise the cohort, its cancer samples are advanced-stage
(publication 1 is Stage IV throughout; publication 9 is metastatic).

The characterisation is recorded in
``results/external_cohort_stage_annotation.json``, with the record it
came from, the sentence it says and the date it was read. That artifact
is the only place the claim lives, and these tests are what keep the
documents and the artifact from drifting apart:

* the artifact must carry a source, a quote and a retrieval date for
  every publication it describes -- a characterisation with no source
  is an assertion;
* every document that states the external numbers must state the stage
  as well, and must name the artifact, so the claim traces to something
  a reader can check;
* the statuses are pinned, so a publication cannot quietly move from
  "unresolved" to "resolved" (or the other way) without the pin being
  updated deliberately.

What this test does not do
--------------------------
It does not fetch anything. The sources are public records, read on the
date the artifact records, and the artifact is the record of that
reading -- the same contract the other published artifacts have.
"""
from __future__ import annotations

import io
import json
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

_ARTIFACT = "results/external_cohort_stage_annotation.json"

# The documents that read the external AUCs, and what each of them has to
# say about the stage. The phrase list is the characterisation itself:
# "Stage IV" for publication 1 and "metastatic" for publication 9.
STATED_IN = {
    "docs/EXTERNAL_VALIDATION.md": ["Stage IV", "metastatic"],
    "paper/paper.tex": ["Stage IV", "metastatic"],
    "paper/METHODS_PAPER.md": ["Stage IV", "metastatic"],
    "ULTRA_EARLY_READINESS.md": ["Stage IV", "metastatic"],
    "docs/PUBLIC_DATA_SOURCES.md": ["Stage IV", "metastatic"],
}

# The stage status of each publication in the held-out cohort. Publication
# 7 is "unresolved": its per-sample annotation is in the paper's
# supplementary tables, which are served behind a proof-of-work and
# reCAPTCHA gate. Pinning it keeps "unresolved" from being quietly
# turned into a stage.
EXPECTED_STATUS = {
    "1": "all_stage_iv",
    "7": "unresolved",
    "9": "metastatic",
}

_ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def _artifact() -> dict:
    path = _REPO_ROOT / _ARTIFACT
    if not path.is_file():
        pytest.skip(f"{_ARTIFACT} is not present")
    return json.loads(path.read_text(encoding="utf-8"))


def _text(rel: str) -> str:
    """The document's text, with LaTeX's escaped underscores normalised."""
    raw = io.open(_REPO_ROOT / rel, encoding="utf-8").read()
    return raw.replace("\\_", "_")


def test_the_artifact_records_a_source_and_a_quote_for_every_publication():
    """A characterisation with no source is an assertion."""
    artifact = _artifact()
    per_publication = artifact.get("per_publication")
    assert per_publication, f"{_ARTIFACT} has no per_publication entries"
    assert set(per_publication) == set(EXPECTED_STATUS), (
        "the held-out cohort's publications changed; update the pins")
    for pub, entry in per_publication.items():
        for field in ("stage_status", "stage_label", "source", "source_url",
                      "quote", "retrieved", "note"):
            value = entry.get(field)
            assert value, f"publication {pub} has no {field} in {_ARTIFACT}"
        assert _ISO_DATE.match(entry["retrieved"]), (
            f"publication {pub}'s retrieval date is not a plain date: "
            f"{entry['retrieved']!r}")
        assert "http" in entry["source_url"], (
            f"publication {pub}'s source has no URL a reader can check")
    for field in ("question", "method", "conclusion", "implication",
                  "honest_framing"):
        assert artifact.get(field), f"{_ARTIFACT} has no {field}"


def test_the_stage_status_of_each_publication_is_pinned():
    """Moving a publication from 'unresolved' to a stage is deliberate."""
    artifact = _artifact()
    found = {pub: entry.get("stage_status")
             for pub, entry in artifact["per_publication"].items()}
    assert found == EXPECTED_STATUS, (
        f"the stage statuses changed: {found}. If a collaborator resolved "
        "publication 7 from the supplementary tables, update the pins here "
        "and in the documents at the same time.")


def test_every_document_that_reads_the_external_numbers_states_the_stage():
    """The claim and the numbers must travel together."""
    problems: list[str] = []
    for doc, phrases in STATED_IN.items():
        text = _text(doc)
        if "external_cohort_stage_annotation.json" not in text:
            problems.append(f"{doc} does not name {_ARTIFACT}")
        for phrase in phrases:
            if phrase not in text:
                problems.append(
                    f"{doc} does not state {phrase!r} (the cohort's stage)")
    assert not problems, "\n  ".join(problems)


def test_no_document_claims_a_stage_the_artifact_does_not_carry():
    """Only the characterisation the artifact records may be stated.

    Publication 1 is Stage IV and publication 9 is metastatic; nothing
    else is established. A document that says the held-out cohort has
    early-stage cases, or that its stages are known, is claiming more
    than the public records support.
    """
    artifact = _artifact()
    assert artifact["per_publication"]["7"]["stage_status"] == "unresolved", (
        "publication 7 is resolved; this test's claims need revisiting")
    overclaiming = (
        "stage I sensitivity",
        "early-stage sensitivity",
        "the held-out cohort's stages are known",
        "stage labels are available",
    )
    problems: list[str] = []
    for doc in STATED_IN:
        text = _text(doc).lower()
        for phrase in overclaiming:
            if phrase in text:
                problems.append(f"{doc} claims {phrase!r}")
    assert not problems, (
        "a document claims more than the public records support:\n  "
        + "\n  ".join(problems))
