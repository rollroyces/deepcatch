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

# The same cohort, measured: what the deposited methylation matrix
# contains, and what can be read from it.
_METHYLATION_ARTIFACT = "results/methylation_cohort_geo_gse124600.json"
_METHYLATION_COVERAGE_AUC = 0.7688
_METHYLATION_STAGE_RHO = 0.0961
_METHYLATION_ZERO_FRACTION = 0.911

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


# ─────────────────────────────────────────────────────────────────────
# What the deposited methylation matrix actually contains.
#
# The survey established that GSE124600 carries stage-labelled plasma
# samples and per-sample methylation matrices. It did not open the
# matrices. Opening them answers the readiness document's question
# ("a methylation cohort would need a methylation feature extractor and
# its own model") with a measurement instead of a plan:
#
#   * the matrix is a methylated-CpG CALL matrix with normalised counts
#     (91% of the cells are 0, the maximum is 6,949.25 -- so a cell is
#     "this CpG was not called methylated in this sample", not "0%
#     methylated");
#   * the per-locus summaries are ANTI-correlated with cancer
#     (AUC 0.27-0.41) and the only class signal is the number of CpGs
#     called methylated per sample (AUC 0.7688) -- a burden/coverage
#     signal, which is what this repository's own methylation channel
#     found on FinaleMe;
#   * the called-CpG count does NOT rise with stage (Spearman rho
#     0.0961, p = 0.26 over the 142 cancer samples), so the
#     stage-labelled cohort does not deliver a stage-stratified result.
#
# The last point is the one the readiness document was waiting for, and
# it is a negative answer. It has to stay attached to the numbers.
# ─────────────────────────────────────────────────────────────────────


def _methylation_artifact() -> dict:
    path = _REPO_ROOT / _METHYLATION_ARTIFACT
    if not path.is_file():
        pytest.skip(f"{_METHYLATION_ARTIFACT} is not present")
    return json.loads(path.read_text(encoding="utf-8"))


def test_the_methylation_cohort_artifact_is_internally_consistent():
    """The measured matrix must agree with the survey's cohort."""
    artifact = _methylation_artifact()
    cohort = artifact["cohort"]
    expected_cancer = EXPECTED_PLASMA_SAMPLES - EXPECTED_CONTROLS
    assert cohort["n_cancer"] == expected_cancer
    assert cohort["n_healthy"] == EXPECTED_CONTROLS
    assert cohort["n_samples"] == EXPECTED_PLASMA_SAMPLES
    assert sum(cohort["stage_distribution"].values()) == expected_cancer
    assert set(cohort["stage_distribution"]) == set(REQUIRED_STAGE_LABELS)

    semantics = artifact["matrix_semantics"]
    zero_fraction = semantics["fraction_of_cells_exactly_zero"]
    assert zero_fraction == pytest.approx(_METHYLATION_ZERO_FRACTION, abs=1e-3)
    # A fractional methylome cannot exceed 1; a count can exceed 6,000.
    # That is what tells the two apart, and it is why the per-locus
    # summaries of this matrix are anti-correlated with cancer.
    assert semantics["max_value"] > 100, (
        "the maximum value is no longer above 100; if the matrix now holds "
        "fractions, the reading of the whole artifact changes and this "
        "test's pins have to be revisited")

    for block in artifact["summaries"].values():
        for name, value in block.items():
            if isinstance(value, dict) and "auc" in value:
                assert 0.0 <= value["auc"] <= 1.0
                assert 0.0 < value["permutation_p"] <= 1.0

    per_stage = artifact["per_stage_called_loci"]
    assert sum(block["n"] for block in per_stage.values()) == expected_cancer
    for stage, block in per_stage.items():
        assert block["n"] == artifact["cohort"]["stage_distribution"][stage]
        assert 0.0 <= block["auc_vs_healthy"] <= 1.0


def test_the_methylation_findings_are_pinned():
    """The numbers that answer the readiness gap cannot drift silently."""
    artifact = _methylation_artifact()
    coverage_auc = None
    for block in artifact["summaries"].values():
        if "called_loci_count" in block:
            coverage_auc = block["called_loci_count"]["auc"]
    assert coverage_auc == _METHYLATION_COVERAGE_AUC, (
        f"the called-locus AUC changed to {coverage_auc}. If the cohort or "
        "the analysis changed, update the pins here and in the documents "
        "at the same time.")
    rho = artifact["stage_monotonicity"]["spearman_rho_stage_vs_called_loci"]
    assert rho == _METHYLATION_STAGE_RHO, (
        f"the stage monotonicity changed to {rho}; the readiness document's "
        "answer to the stage gap moves with it")
    assert artifact["stage_monotonicity"]["p"] > 0.05, (
        "the called-CpG count now rises with stage; the negative answer "
        "about the stage gap no longer holds and the documents must say so")


def test_every_document_quoting_the_methylation_finding_names_its_artifact():
    """A finding that answers the readiness gap has to be checkable."""
    needle = str(_METHYLATION_COVERAGE_AUC)
    documents = _documents_mentioning("GSE124600")
    problems: list[str] = []
    for path in documents:
        text = path.read_text(encoding="utf-8")
        if needle not in text:
            continue
        rel = path.relative_to(_REPO_ROOT)
        if Path(_METHYLATION_ARTIFACT).name not in text:
            problems.append(f"{rel} does not name {_METHYLATION_ARTIFACT}")
        if str(_METHYLATION_STAGE_RHO) not in text:
            problems.append(
                f"{rel} does not state the stage monotonicity "
                f"({_METHYLATION_STAGE_RHO}) -- the negative answer that "
                "travels with the finding")
    assert not problems, "\n  ".join(problems)
