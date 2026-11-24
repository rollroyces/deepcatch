"""The Zenodo survey and the ccRCC cohort must agree with their artifacts.

Why this file exists
--------------------
The goal asks for public cfDNA features that need no IRB, and Zenodo had
two known records and no search. The search is now recorded
(`results/zenodo_cfdna_survey.json`), and one of the records it found --
Zenodo 21269430, per-sample hg38 fragment BED files from six ccRCC
patients and six healthy controls -- was extracted with the repository's
own extractor and scored with the same protocol as the FinaleDB
held-out cohorts (`results/external_cohort_zenodo_ccrcc.json`).

The number is the interesting one: AUC 0.9444 [0.8156, 1.0] with an
exact permutation p of 0.0043, in a cohort of a different cancer type,
a different laboratory and a different data source from the one the
manuscript's held-out cohort comes from. It is also the number most
likely to be misread, because the same cohort carries three honest
caveats -- every patient is pT3a (stage III), the assay is about 30x
shallower than the training cohort's, and the absolute threshold does
not travel (the ranking does). Those caveats belong with the number,
so these tests check that they travel with it:

* the artifacts' own arithmetic must hold;
* the cohort's size, stage and numbers are pinned, so a re-run that
  moves them has to be deliberate;
* every document that mentions the accession must state the cohort's
  numbers and its limits;
* the extractor's BED-header fix is tested both ways -- a file with a
  header and the same file without one must give the same features.

What this test does not do
--------------------------
It does not fetch anything. The sources are public records, read on the
date the artifacts record, and the artifacts are the record of that
reading.
"""
from __future__ import annotations

import io
import json
import re
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

_SURVEY_ARTIFACT = "results/zenodo_cfdna_survey.json"
_COHORT_ARTIFACT = "results/external_cohort_zenodo_ccrcc.json"

_CCRCC_DOI = "10.5281/zenodo.21269430"

# The cohort, pinned. A re-run that moves any of these has to be a
# deliberate act, not a silent drift.
EXPECTED_N_CANCER = 6
EXPECTED_N_HEALTHY = 6
EXPECTED_AUC = 0.9444
EXPECTED_CI_LO = 0.8156
EXPECTED_CI_HI = 1.0
EXPECTED_PERMUTATION_P = 0.0043
EXPECTED_PERMUTATIONS = 924
EXPECTED_STAGE = "pT3a"
EXPECTED_SPECIFICITY_AT_TRAIN_THRESHOLD = 0.5

_ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

# The limits that have to travel with the number. A document that states
# the AUC without them is overclaiming: the cohort is stage III, and the
# score scale does not travel. The threshold caveat is required as the
# word rather than as a number -- "0.5" is not a distinctive string, and
# a limit that cannot be recognised is not a check.
REQUIRED_LIMITS = (
    "pT3a",
    "0.0043",
    "threshold",
)


def _load(rel: str) -> dict:
    path = _REPO_ROOT / rel
    if not path.is_file():
        pytest.skip(f"{rel} is not present")
    return json.loads(path.read_text(encoding="utf-8"))


def _survey() -> dict:
    return _load(_SURVEY_ARTIFACT)


def _cohort() -> dict:
    return _load(_COHORT_ARTIFACT)


def _documents_mentioning(needle: str) -> list[Path]:
    found: list[Path] = []
    for path in sorted(_REPO_ROOT.rglob("*.md")):
        parts = path.relative_to(_REPO_ROOT).parts
        if any(part in {".git", "node_modules", "__pycache__", "dist",
                        "build", ".pytest_cache"} for part in parts):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if needle in text:
            found.append(path)
    return found


def test_the_survey_records_its_queries_and_the_candidates_it_read():
    """A survey with no queries is not a survey."""
    survey = _survey()
    assert _ISO_DATE.match(survey.get("surveyed", "")), (
        "the survey artifact has no plain survey date")
    method = survey.get("method") or {}
    queries = method.get("queries") or []
    assert len(queries) >= 4, (
        "the Zenodo survey ran fewer queries than it claims")
    for entry in queries:
        assert entry.get("query"), "a query entry has no query text"
        if "error" not in entry:
            assert isinstance(entry.get("dataset_hits"), int), (
                f"the query {entry['query']!r} has no hit count")
    candidates = survey.get("candidates") or []
    assert candidates, "the survey artifact records no candidate records"
    for candidate in candidates:
        assert candidate.get("doi"), "a candidate has no DOI"
        if "error" not in candidate:
            assert candidate.get("title"), (
                f"{candidate['doi']} has no title")
            assert candidate.get("publication_date"), (
                f"{candidate['doi']} has no publication date")
            assert isinstance(candidate.get("n_files"), int), (
                f"{candidate['doi']} has no file count")


def test_the_ccrcc_record_is_one_of_the_candidates():
    """The cohort that was scored must be in the survey that found it."""
    survey = _survey()
    dois = {c.get("doi") for c in survey.get("candidates") or []}
    assert _CCRCC_DOI in dois, (
        f"{_CCRCC_DOI} is not among the survey's candidates")
    candidate = next(c for c in survey["candidates"] if c["doi"] == _CCRCC_DOI)
    assert candidate["n_files"] >= 2, (
        "the ccRCC record's file list does not hold the fragment archive and "
        "the sample metadata")
    assert candidate["total_bytes"] > 0


def test_the_cohort_artifact_is_internally_consistent():
    """The numbers in the artifact must add up to each other."""
    artifact = _cohort()
    result = artifact["result"]
    cohort = artifact["cohort"]
    assert cohort["n"] == EXPECTED_N_CANCER + EXPECTED_N_HEALTHY
    assert result["n_pos"] == EXPECTED_N_CANCER
    assert result["n_neg"] == EXPECTED_N_HEALTHY
    assert result["ci_lo"] <= result["auc"] <= result["ci_hi"], (
        "the DeLong interval does not contain the point estimate")
    # The specificity at the transferred threshold is a count, not a
    # claim: it must be reproducible from the two counts.
    healthy_above = result["n_healthy_above_train_threshold"]
    assert result["specificity_at_train_threshold"] == pytest.approx(
        1.0 - healthy_above / EXPECTED_N_HEALTHY, abs=1e-4)
    assert result["n_cancer_above_train_threshold"] + healthy_above == (
        result["n_above_train_threshold"])
    # The permutation test is exact over all C(12,6) assignments.
    assert result["exact_permutations"] == EXPECTED_PERMUTATIONS
    assert 0.0 < result["exact_permutation_p"] <= 1.0
    # Every cancer sample is labelled 1 and every control 0.
    labels = cohort["labels"]
    assert sum(labels.values()) == EXPECTED_N_CANCER
    assert len(labels) == EXPECTED_N_HEALTHY + EXPECTED_N_CANCER


def test_the_cohort_numbers_are_pinned():
    """The cohort cannot quietly change size, stage or number."""
    artifact = _cohort()
    result = artifact["result"]
    cohort = artifact["cohort"]
    assert result["auc"] == EXPECTED_AUC, (
        f"the ccRCC AUC changed to {result['auc']}. If the cohort or the "
        "protocol changed, update the pins here and in the documents at the "
        "same time.")
    assert result["ci_lo"] == EXPECTED_CI_LO
    assert result["ci_hi"] == EXPECTED_CI_HI
    assert round(result["exact_permutation_p"], 4) == EXPECTED_PERMUTATION_P
    assert result["specificity_at_train_threshold"] == (
        EXPECTED_SPECIFICITY_AT_TRAIN_THRESHOLD)
    assert cohort["doi"] == _CCRCC_DOI
    # The stage of the PATIENTS. The controls carry a placeholder ("-"),
    # which is not a stage and must not be read as one.
    stages = {stage for sample, stage in cohort["stage"].items()
              if cohort["labels"][sample] == 1 and stage
              and stage != "-"}
    assert stages == {EXPECTED_STAGE}, (
        f"the cohort's stages changed to {sorted(stages)}; the limit that "
        "travels with the number has to change with it")
    assert cohort["n_cancer"] == EXPECTED_N_CANCER
    assert cohort["n_healthy"] == EXPECTED_N_HEALTHY


def test_every_document_that_mentions_the_cohort_states_its_limits():
    """The number and the caveats travel together.

    The documents are found by scanning the tree, so a document that
    quotes the cohort without its limits is covered without anyone
    remembering to register it.
    """
    documents = _documents_mentioning(_CCRCC_DOI)
    assert documents, "nothing mentions the ccRCC cohort"
    problems: list[str] = []
    for path in documents:
        rel = path.relative_to(_REPO_ROOT)
        text = path.read_text(encoding="utf-8")
        if Path(_COHORT_ARTIFACT).name not in text:
            problems.append(f"{rel} does not name {_COHORT_ARTIFACT}")
        if str(EXPECTED_AUC) not in text:
            problems.append(f"{rel} does not state the AUC")
        for limit in REQUIRED_LIMITS:
            if limit not in text:
                problems.append(
                    f"{rel} does not state {limit!r} -- one of the limits that "
                    "has to travel with the number")
    assert not problems, "\n  ".join(problems)


def test_the_extractor_reads_a_bed_file_with_a_header():
    """The deposited files are BED files, and BED files have a header."""
    import sys

    sys.path.insert(0, str(_REPO_ROOT / "src"))
    from fragmentomics import finale_features as ff

    rows = [
        "#chrom\tstart\tend\tname\tmapq\tstrand",
        "chr1\t1000\t1120\t.\t60\t+",
        "chr1\t2000\t2200\t.\t60\t-",
        "chr2\t3000\t3300\t.\t60\t+",
    ]
    body = "\n".join(rows[1:]) + "\n"

    with_header = _REPO_ROOT / "results" / "_tmp_header_test.frag.bed.gz"
    without_header = _REPO_ROOT / "results" / "_tmp_noheader_test.frag.bed.gz"
    try:
        import gzip
        with gzip.open(with_header, "wt") as handle:
            handle.write("\n".join(rows) + "\n")
        with gzip.open(without_header, "wt") as handle:
            handle.write(body)
        with_header_codes, _, _, _ = ff.read_fragments(str(with_header))
        plain_codes, _, _, _ = ff.read_fragments(str(without_header))
        assert len(with_header_codes) == 3, (
            "the header line was read as data (or a row was dropped)")
        assert list(with_header_codes) == list(plain_codes), (
            "a file with a header and the same file without one must give the "
            "same fragments")
    finally:
        for path in (with_header, without_header):
            if path.exists():
                path.unlink()


def test_the_cohort_artifact_names_its_sources():
    """The artifact must say where the data came from."""
    artifact = _cohort()
    cohort = artifact["cohort"]
    assert cohort["record"].startswith("https://zenodo.org/records/")
    design = artifact.get("design") or {}
    for field in ("question", "train", "test", "harmonization", "model",
                  "threshold"):
        assert design.get(field), f"the artifact's design has no {field}"
    limits = artifact.get("limits") or []
    assert len(limits) >= 3, (
        "the artifact records fewer limits than the cohort has")
    for limit in limits:
        assert limit.strip(), "an empty limit string is not a limit"
