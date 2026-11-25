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
_PHTS_ARTIFACT = "results/external_cohort_zenodo_phts.json"
_PHTS_PREP_ARTIFACT = "results/zenodo_phts_prep.json"

# The column names a deposit may use for its per-sample id.
_SAMPLE_COLUMNS = ("sample", "sample_id", "id")

_CCRCC_DOI = "10.5281/zenodo.21269430"
_PHTS_DOI = "10.5281/zenodo.15109495"

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

# The Cleveland Clinic PHTS screening cohort, pinned. The arm sizes
# are the depositor's own split and do not move; the numbers below
# them are the measured ones.
EXPECTED_PHTS_N = 99
EXPECTED_PHTS_N_CANCER = 49
EXPECTED_PHTS_N_HEALTHY = 50
EXPECTED_PHTS_PRE_DIAGNOSIS = 9
EXPECTED_PHTS_POST_DIAGNOSIS = 40
EXPECTED_PHTS_N_STAGED = 39
EXPECTED_PHTS_MAX_LEAD_TIME_DAYS = 7
EXPECTED_PHTS_AUC = 0.4506
EXPECTED_PHTS_CI_LO = 0.3343
EXPECTED_PHTS_CI_HI = 0.5669
EXPECTED_PHTS_PERMUTATION_P = 0.8025
EXPECTED_PHTS_PRE_AUC = 0.5
EXPECTED_PHTS_SHORT_LONG_AUC = 0.6455
EXPECTED_PHTS_ABSOLUTE_AUC = 0.5388

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

# The screening cohort's limits. "pre-diagnosis" reads as an
# early-detection arm unless the interval is stated next to it, and the
# arm's size is what makes it a measurement rather than a result.
REQUIRED_PHTS_LIMITS = (
    "lead time",
    "7 days",
    "n=9",
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


def _phts() -> dict:
    return _load(_PHTS_ARTIFACT)


def _phts_prep() -> dict:
    return _load(_PHTS_PREP_ARTIFACT)


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
            # A deposit that carries a per-sample cohort table must have
            # been read: the labels live in the table, not in the record's
            # description, and a survey that stops at the description
            # reports a labelled cohort as unlabelled.
            table = candidate.get("cohort_table")
            if candidate.get("cohort_table_files"):
                assert table is not None, (
                    f"{candidate['doi']} declares a cohort table "
                    f"({', '.join(candidate['cohort_table_files'])}) but the "
                    "survey did not read it")
            if table is not None:
                assert table.get("columns"), (
                    f"{candidate['doi']}'s cohort table has no columns")
                assert table.get("rows"), (
                    f"{candidate['doi']}'s cohort table has no rows")
                for row in table["rows"]:
                    # A per-sample table identifies its samples in a column
                    # whose name is not standardised across deposits:
                    # SAUNA calls it "sample", the Cleveland Clinic PHTS
                    # cohort calls it "id" (its values are sample ids of
                    # the form CCF00016-01-001). What matters is that a
                    # row can be tied to a sample.
                    assert any(name in row for name in _SAMPLE_COLUMNS), (
                        f"{candidate['doi']}'s cohort table has a row without "
                        "a sample column (none of "
                        f"{', '.join(_SAMPLE_COLUMNS)})")


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

    # The same arithmetic for the screening cohort: 99 samples, 49 with
    # a malignant neoplasm and 50 without.
    phts = _phts()
    phts_cohort = phts["cohort"]
    phts_result = phts["result"]
    assert phts_cohort["n"] == EXPECTED_PHTS_N
    assert phts_cohort["n_cancer"] == EXPECTED_PHTS_N_CANCER
    assert phts_cohort["n_healthy"] == EXPECTED_PHTS_N_HEALTHY
    assert phts_result["n_pos"] == EXPECTED_PHTS_N_CANCER
    assert phts_result["n_neg"] == EXPECTED_PHTS_N_HEALTHY
    assert phts_result["ci_lo"] <= phts_result["auc"] <= phts_result["ci_hi"], (
        "the DeLong interval does not contain the point estimate")
    assert len(phts_cohort["labels"]) == EXPECTED_PHTS_N
    assert sum(phts_cohort["labels"].values()) == EXPECTED_PHTS_N_CANCER
    # The arms are the depositor's own split, and they have to add up.
    arms = phts_cohort["arms"]
    assert arms["no_cancer"] == EXPECTED_PHTS_N_HEALTHY
    assert arms["pre_smn"] + arms["pre_pmn"] == EXPECTED_PHTS_PRE_DIAGNOSIS
    assert arms["post_smn"] + arms["post_pmn"] == EXPECTED_PHTS_POST_DIAGNOSIS
    assert (arms["no_cancer"] + arms["pre_smn"] + arms["post_smn"]
            + arms["pre_pmn"] + arms["post_pmn"]) == EXPECTED_PHTS_N
    # The stage counts are the depositor's: 0: 5, 1: 22, 2: 10, 3: 2.
    counts = phts_cohort["stage_counts"]
    assert sum(counts.values()) == EXPECTED_PHTS_N_STAGED
    assert counts["0"] == 5 and counts["1"] == 22
    assert counts["2"] == 10 and counts["3"] == 2


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

    # The screening cohort's numbers, pinned the same way. These are the
    # measured ones, and a re-run that moves them has to be deliberate.
    phts = _phts()
    phts_result = phts["result"]
    phts_cohort = phts["cohort"]
    assert phts_result["auc"] == EXPECTED_PHTS_AUC, (
        f"the screening cohort's AUC changed to {phts_result['auc']}. If "
        "the cohort or the protocol changed, update the pins here and in "
        "the documents at the same time.")
    assert phts_result["ci_lo"] == EXPECTED_PHTS_CI_LO
    assert phts_result["ci_hi"] == EXPECTED_PHTS_CI_HI
    assert round(phts_result["permutation_p"], 4) == EXPECTED_PHTS_PERMUTATION_P
    assert phts_result["score_saturation"]["saturated"] is True, (
        "the cohort's scores are no longer saturated; the mechanism the "
        "null is attributed to has to be re-measured")
    # The stage arms, pinned: a null in every one of them.
    arms = phts["stage_question"]["arms_against_no_cancer"]
    assert arms["pre_diagnosis"]["auc"] == EXPECTED_PHTS_PRE_AUC
    assert arms["pre_diagnosis"]["n_pos"] == EXPECTED_PHTS_PRE_DIAGNOSIS
    for name in ("stage_0", "stage_0_1", "stage_2_3", "all_staged"):
        block = arms[name]
        assert block["ci_lo"] <= 0.5 <= block["ci_hi"], (
            f"the {name} arm's interval no longer includes 0.5 "
            f"({block['ci_lo']}, {block['ci_hi']}) -- the null has moved")
    # The channel that DOES separate, and the direction the depositor's
    # own paper reports: this is what makes the null a transfer failure
    # rather than an absence of signal, so it is pinned too.
    channel = phts["single_channel_readout"]["genome_short_long_ratio"]
    assert channel["post_smn"]["auc"] == EXPECTED_PHTS_SHORT_LONG_AUC
    assert channel["post_smn"]["ci_lo"] > 0.5, (
        "the short/long ratio no longer separates the second-neoplasm arm; "
        "the mechanism the null is attributed to has to be re-measured")
    assert channel["post_smn"]["median_pos"] > channel["post_smn"]["median_neg"], (
        "the direction of the depositor's published effect has flipped")
    assert phts["published_direction_check"]["direction_reproduced"] is True
    assert phts_cohort["doi"] == _PHTS_DOI
    # The control that normalises the held-out cohort by the TRAINING
    # cohort's statistics, so its absolute position in feature space is
    # preserved. If the null is an artefact of the harmonisation, this is
    # where it would show.
    absolute = phts["harmonization_sensitivity"][
        "train_global_test_train_statistics"]
    assert absolute["auc"] == EXPECTED_PHTS_ABSOLUTE_AUC
    assert absolute["ci_lo"] <= 0.5 <= absolute["ci_hi"], (
        "the absolute-harmonisation control no longer includes 0.5 -- the "
        "null has moved and the mechanism attributed to it has to be "
        "re-measured")


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

    # The screening cohort has the same obligation, and its limits are
    # the ones a reader is most likely to skip: that the "pre-diagnosis"
    # arm was drawn within a week of the diagnosis, and that it holds
    # nine samples.
    phts_documents = _documents_mentioning(_PHTS_DOI)
    assert phts_documents, "nothing mentions the screening cohort"
    problems = []
    for path in phts_documents:
        rel = path.relative_to(_REPO_ROOT)
        text = path.read_text(encoding="utf-8")
        if Path(_PHTS_ARTIFACT).name not in text:
            problems.append(f"{rel} does not name {_PHTS_ARTIFACT}")
        for limit in REQUIRED_PHTS_LIMITS:
            if limit not in text:
                problems.append(
                    f"{rel} does not state {limit!r} -- one of the limits "
                    "that has to travel with the number")
    assert not problems, "\n  ".join(problems)


def test_the_extractor_reads_a_bed_file_with_a_header():
    """The deposited files are BED files, and BED files have a header."""
    import gzip
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

    # The Cleveland Clinic PHTS deposit ships BEDPE, not fragment BED:
    # one row per read PAIR, with the two mates' coordinates in separate
    # columns. The extractor's first three columns are the FIRST mate's
    # coordinates, so the extractor reads a BEDPE as a set of fragments
    # that are one read long. It does not fail -- it is quietly wrong.
    # The preparation script is what converts one into the other, and
    # this is the check that the conversion is the right one.
    sys.path.insert(0, str(_REPO_ROOT / "scripts"))
    from prepare_zenodo_phts_cohort import convert_member

    bedpe_rows = [
        "chr1\t1000\t1120\tchr1\t1100\t1180\tread1\t60\t+\t-\t180",
        "chr1\t5000\t5146\tchr1\t5000\t5099\tread2\t60\t-\t+\t146",
        "chr2\t2000\t2146\tchr2\t2100\t2200\tread3\t60\t+\t-\t200",
    ]
    raw = _REPO_ROOT / "results" / "_tmp_bedpe_test.frag.bed.gz"
    converted = _REPO_ROOT / "results" / "_tmp_bedpe_converted.frag.bed.gz"
    try:
        with gzip.open(raw, "wt") as handle:
            handle.write("\n".join(bedpe_rows) + "\n")
        stats: dict = {}
        with gzip.open(raw, "rt") as handle:
            convert_member(handle, str(converted), stats)
        assert stats["rows"] == 3, "a converted pair was dropped"
        assert stats["rows_dropped"] == 0

        _codes, start, end, _mapq = ff.read_fragments(str(converted))
        converted_lengths = sorted(int(e - s) for s, e in zip(start, end))
        assert converted_lengths == [146, 180, 200], (
            "the converted file does not carry the fragments' outer "
            f"coordinates (got {converted_lengths})")

        # The same file read without the conversion gives the first
        # mate's length instead of the fragment's -- which is the whole
        # reason the conversion exists.
        _codes, start, end, _mapq = ff.read_fragments(str(raw))
        raw_lengths = sorted(int(e - s) for s, e in zip(start, end))
        assert raw_lengths == [120, 146, 146], (
            "the unconverted BEDPE no longer reads as the first mate's "
            f"coordinates (got {raw_lengths}); the check that follows it "
            "would no longer mean anything")
        assert raw_lengths != converted_lengths
    finally:
        for path in (raw, converted):
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

    # The screening cohort's artifact says the same kind of things,
    # and it has to carry the depositor's own paper.
    phts = _phts()
    phts_cohort = phts["cohort"]
    assert phts_cohort["record"].startswith("https://zenodo.org/records/")
    assert phts_cohort["doi"] == _PHTS_DOI
    phts_design = phts.get("design") or {}
    for field in ("question", "train", "test", "harmonization", "model",
                  "threshold", "paper"):
        assert phts_design.get(field), (
            f"the screening cohort's design has no {field}")
    assert "xcrm" in phts_design["paper"], (
        "the screening cohort's artifact does not name the depositor's paper")
    phts_limits = phts.get("limits") or []
    assert len(phts_limits) >= 3, (
        "the screening cohort's artifact records fewer limits than the "
        "cohort has")
    for limit in phts_limits:
        assert limit.strip(), "an empty limit string is not a limit"

    # The lead time is the fact that decides what the "pre" arms are.
    # It has to be in the artifact, and the limit that says it has to
    # be in the limits -- a reader who sees "pre-diagnosis" without it
    # would read an early-detection claim into it.
    lead = phts.get("lead_time") or {}
    assert lead.get("max_lead_time_days") == EXPECTED_PHTS_MAX_LEAD_TIME_DAYS, (
        "the cohort's longest lead time changed; the limit that travels "
        "with the numbers has to change with them")
    assert lead.get("by_arm", {}).get("pre_diagnosis"), (
        "the artifact does not report the pre-diagnosis arm's intervals")
    assert any("lead-time" in limit for limit in phts_limits), (
        "the limits do not say that this is not a lead-time cohort")
    verdict = (phts.get("stage_question") or {}).get("verdict") or {}
    assert verdict.get("is_a_lead_time_cohort") is False, (
        "the verdict claims this is a lead-time cohort; the depositor's "
        "own date_diff column says otherwise")

    # The preparation report records the conversion and the genome
    # build, rather than leaving either of them implicit.
    prep = _phts_prep()
    conversion = prep.get("conversion") or {}
    assert "BEDPE" in conversion.get("from", ""), (
        "the preparation report does not say what the deposit ships")
    assert conversion.get("target_suffix"), (
        "the preparation report does not name the files it wrote")
    build = prep.get("genome_build") or {}
    assert build.get("verdict") == "hg19", (
        f"the deposit's genome build reads as {build.get('verdict')!r}; the "
        "extractor's grid is hg38 and the mismatch has to be accounted for")
    assert build.get("chromosomes_past_the_hg38_length"), (
        "the build verdict names no chromosome that runs past hg38's "
        "length, so it is not evidence of anything")
    assert build.get("depositor_statement"), (
        "the build verdict is not cross-checked against the depositor's "
        "own statement")
    dropped = build.get("hg38_grid_fraction_dropped")
    assert dropped is not None and 0.0 <= dropped < 0.05, (
        "the fraction of fragments the hg38 grid drops is not recorded, or "
        "it is not the fraction the coordinates imply")
    samples = prep.get("per_sample") or {}
    assert samples, "the preparation report records no per-sample check"
    for sample, stats in samples.items():
        assert stats["mode_fragment_length"] >= 150, (
            f"{sample}'s fragment lengths do not peak on the nucleosome "
            "peak")
        assert stats["di_nucleosome_fraction"] >= 0.008, (
            f"{sample}'s fragment lengths have no di-nucleosome shoulder")
