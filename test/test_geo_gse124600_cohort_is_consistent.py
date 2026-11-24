"""The GSE124600 cohort size must be the one GEO actually records.

Why this file exists
--------------------
The GEO survey said GSE124600 carries "142 cancer and 93 healthy plasma
samples" with `crc grade: Stage_I`, `Stage_II` and `Stage_III`. The
control count was wrong, and the stage set was incomplete. The series
record says plasma from 147 CRC patients and 134 controls, of which 5
patients and 2 controls failed quality control, so the deposited cohort
is 142 cancer + 132 control = 274 plasma samples -- and 5 of the cancer
samples are `Stage_IV`, which the survey had left out.

A wrong cohort size in a document a reviewer will read is not a typo: it
changes what the cohort is, and "93 healthy" against a real 132 is the
kind of error that makes a reader distrust every other number nearby.
So the measurements live in ``results/geo_gse124600_cohort.json`` --
counted from the deposited supplementary files and the GEO sample
records, not quoted from prose -- and these tests are what keep the
documents and that artifact from drifting apart again:

* the artifact's arithmetic must hold: the quality-control arithmetic,
  the stage distribution, the matrix column counts, and the coverage;
* every document that mentions the accession must state the cohort size
  the artifact records, must name the artifact so the claim is checkable,
  and must not repeat the wrong figure;
* the document check is mechanical -- it finds the documents by scanning
  the tree for the accession, so a new document that quotes GSE124600 is
  covered without anyone remembering to add it to a list.

What this test does not do
--------------------------
It does not fetch anything. The sources are public records, read on the
date the artifact records, and the artifact is the record of that
reading -- the same contract as
``results/external_cohort_stage_annotation.json``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

_ARTIFACT = "results/geo_gse124600_cohort.json"

# The cohort as GEO records it. These are pins, not documentation: if the
# artifact's numbers change, the pins change deliberately and in the same
# commit, rather than a wrong number creeping in with a file rename.
EXPECTED_PLASMA_SAMPLES = 274
EXPECTED_CRC_PATIENTS = 142
EXPECTED_CONTROLS = 132
EXPECTED_CPG_LOCI = 205250
EXPECTED_TISSUE_COLUMNS = 66

EXPECTED_STAGES = {
    "Stage_I": 30,
    "Stage_II": 65,
    "Stage_III": 42,
    "Stage_IV": 5,
}

# The stage labels a document has to spell out. A document that lists the
# cohort's stages but omits Stage_IV is describing 207 of the 274 plasma
# samples, and a reader cannot tell that from the sentence.
REQUIRED_STAGE_LABELS = ("Stage_I", "Stage_II", "Stage_III", "Stage_IV")

# The wrong figure, in the shapes it was written. A document that carries
# one of these is asserting the cohort is bigger in one arm than it is.
FORBIDDEN_FIGURES = ("93 healthy", "93 control", "93 and")

_ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

_SKIP_DIRECTORIES = {
    ".git", ".github", "node_modules", ".venv", "venv", "__pycache__",
    ".pytest_cache", ".mypy_cache", "build", "dist", ".ruff_cache",
}


def _artifact() -> dict:
    path = _REPO_ROOT / _ARTIFACT
    if not path.is_file():
        pytest.skip(f"{_ARTIFACT} is not present")
    return json.loads(path.read_text(encoding="utf-8"))


def _documents_mentioning(accession: str) -> list[Path]:
    """Every markdown document in the tree that mentions the accession."""
    found: list[Path] = []
    for path in sorted(_REPO_ROOT.rglob("*.md")):
        if any(part in _SKIP_DIRECTORIES for part in path.relative_to(_REPO_ROOT).parts):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if accession in text:
            found.append(path)
    return found


def test_the_artifact_records_where_every_number_came_from():
    """A number with no source is an assertion."""
    artifact = _artifact()
    assert artifact.get("accession") == "GSE124600"
    assert _ISO_DATE.match(artifact.get("retrieved", "")), (
        "the artifact has no plain retrieval date")
    sources = artifact.get("sources") or {}
    assert sources, f"{_ARTIFACT} has no sources"
    for name, url in sources.items():
        assert "http" in url, f"source {name!r} has no URL a reader can check"
    design = artifact.get("design") or {}
    assert design.get("quote"), "the artifact records no design quote"
    assert design.get("source"), "the artifact records no design source"
    for matrix in artifact.get("matrices") or []:
        assert matrix.get("file"), "a matrix entry has no file name"
        assert matrix.get("bytes"), f"{matrix.get('file')} has no byte count"
        assert _SHA256.match(matrix.get("sha256", "")), (
            f"{matrix.get('file')} has no sha256")
        assert matrix.get("cpg_loci"), f"{matrix.get('file')} has no locus count"
        assert matrix.get("sample_columns"), (
            f"{matrix.get('file')} has no sample column count")


def test_the_quality_control_arithmetic_holds():
    """142 + 132 = 274, and both arms come from the series record."""
    artifact = _artifact()
    design = artifact["design"]
    assert design["crc_patients_enrolled"] - design["crc_patients_excluded_by_qc"] == (
        EXPECTED_CRC_PATIENTS)
    assert design["controls_enrolled"] - design["controls_excluded_by_qc"] == (
        EXPECTED_CONTROLS)
    assert design["crc_patients_after_qc"] == EXPECTED_CRC_PATIENTS
    assert design["controls_after_qc"] == EXPECTED_CONTROLS
    assert (EXPECTED_CRC_PATIENTS + EXPECTED_CONTROLS
            == EXPECTED_PLASMA_SAMPLES)


def test_the_stage_distribution_accounts_for_every_labelled_sample():
    """Every labelled plasma sample has a stage, and the stages are named."""
    artifact = _artifact()
    stages = artifact["stage_labels"]["distribution"]
    assert stages == EXPECTED_STAGES, (
        f"the stage distribution changed: {stages}. If the labels or the "
        "counts moved, update the pins here and in the documents at the "
        "same time.")
    assert sum(stages.values()) == EXPECTED_CRC_PATIENTS
    assert artifact["stage_labels"]["on_plasma_samples"] == EXPECTED_CRC_PATIENTS
    cohort = artifact["cohort"]
    assert cohort["plasma_samples_with_a_stage_label"] == EXPECTED_CRC_PATIENTS
    assert cohort["plasma_samples_without_a_stage_label"] == EXPECTED_CONTROLS


def test_the_matrices_cover_every_plasma_sample():
    """The columns in the deposited matrices are the samples GEO records."""
    artifact = _artifact()
    matrices = artifact["matrices"]
    plasma = [m for m in matrices if m["sample_kind"] == "plasma"]
    tissue = [m for m in matrices if m["sample_kind"] == "tissue"]
    assert plasma, "no plasma matrix in the artifact"
    assert tissue, "no tissue matrix in the artifact"

    plasma_columns = sum(m["sample_columns"] for m in plasma)
    tissue_columns = sum(m["sample_columns"] for m in tissue)
    assert plasma_columns == EXPECTED_PLASMA_SAMPLES, (
        f"the plasma matrices hold {plasma_columns} columns, not "
        f"{EXPECTED_PLASMA_SAMPLES}")
    assert tissue_columns == EXPECTED_TISSUE_COLUMNS
    assert artifact["coverage"]["plasma_sample_columns_in_matrices"] == (
        plasma_columns)
    assert artifact["coverage"]["plasma_samples_in_geo"] == EXPECTED_PLASMA_SAMPLES
    assert artifact["coverage"]["tissue_sample_columns_in_matrices"] == (
        tissue_columns)
    assert artifact["cohort"]["plasma_samples"] == EXPECTED_PLASMA_SAMPLES

    for matrix in matrices:
        assert matrix["cpg_loci"] == EXPECTED_CPG_LOCI, (
            f"{matrix['file']} holds {matrix['cpg_loci']} loci, not "
            f"{EXPECTED_CPG_LOCI}")


def test_the_cohort_size_is_pinned():
    """The cohort cannot quietly change size."""
    artifact = _artifact()
    assert artifact["cohort"]["plasma_samples"] == EXPECTED_PLASMA_SAMPLES
    assert artifact["coverage"]["plasma_coverage"] == (
        f"{EXPECTED_PLASMA_SAMPLES}/{EXPECTED_PLASMA_SAMPLES}")
    assert artifact["coverage"]["cpg_loci_per_matrix"] == [EXPECTED_CPG_LOCI]


def test_every_document_that_mentions_the_accession_states_the_recorded_cohort():
    """The claim and the numbers must travel together.

    The documents are found by scanning the tree, not from a list, so a
    new document that quotes GSE124600 is covered without anyone
    remembering to register it.
    """
    documents = _documents_mentioning("GSE124600")
    assert documents, "nothing mentions GSE124600; the survey was removed?"
    problems: list[str] = []
    for path in documents:
        rel = path.relative_to(_REPO_ROOT)
        text = path.read_text(encoding="utf-8")
        if Path(_ARTIFACT).name not in text:
            problems.append(f"{rel} does not name {_ARTIFACT}")
        for number in (str(EXPECTED_CRC_PATIENTS), str(EXPECTED_CONTROLS)):
            if number not in text:
                problems.append(
                    f"{rel} does not state {number} "
                    "(the cohort's cancer / control plasma counts)")
        for label in REQUIRED_STAGE_LABELS:
            if label not in text:
                problems.append(f"{rel} does not state {label}")
        for wrong in FORBIDDEN_FIGURES:
            if wrong in text:
                problems.append(
                    f"{rel} states {wrong!r}, which is not the cohort size")
    assert not problems, "\n  ".join(problems)
