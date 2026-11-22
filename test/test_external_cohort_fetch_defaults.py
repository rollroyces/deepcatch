"""The fetcher's defaults must reproduce the committed external cohort.

Why this file exists
--------------------
`scripts/fetch_external_cohort.py` selects the held-out cohort from the
FinaleDB registry and downloads it. Its quota caps how many samples each
publication contributes, and a smaller cap yields a DIFFERENT cohort --
so the published numbers do not reproduce with it.

The quota's default was empty, and `--max-per-publication` defaulted to
0, so the documented command -- the one the reproduction driver tells a
reader to run -- selected nothing and stopped with "nothing selected".
The `--help` example (`9=240,7=31,1=12`) was worse: it ran, and it
selected a cohort that is not the published one (publication 1
contributed 12 samples where the committed cohort has 52).

What this test checks
---------------------
The default quota is not empty, it names every publication the published
artifact reports, and each cap is at least the number of samples that
publication contributed to the committed cohort. That is what makes
`python3 scripts/fetch_external_cohort.py` reproduce the cohort the
manuscript quotes.

The check is offline: it compares the script's defaults with the
committed artifact, so it needs neither the API nor the S3 bucket.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_FETCHER = _REPO_ROOT / "scripts" / "fetch_external_cohort.py"
_ARTIFACT = _REPO_ROOT / "results" / "external_validation_finaledb.json"


def _fetcher():
    """The fetcher script, loaded from its path."""
    spec = importlib.util.spec_from_file_location("fetch_external_cohort", _FETCHER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _default_quota() -> dict[int, int]:
    """The script's default quota, as {publication: cap}."""
    quota = _fetcher().build_parser().get_default("quota")
    out: dict[int, int] = {}
    for part in str(quota).split(","):
        part = part.strip()
        if not part:
            continue
        pub, _, n = part.partition("=")
        out[int(pub)] = int(n)
    return out


def _published_cohort() -> dict[int, int]:
    """The per-publication sample counts of the committed cohort."""
    artifact = json.loads(_ARTIFACT.read_text(encoding="utf-8"))
    return {
        int(pub): int(row["n"])
        for pub, row in artifact["per_publication_auc"].items()
        if row.get("n")
    }


def test_the_default_quota_is_not_empty():
    """An empty default quota selects nothing, and says only that."""
    quota = _default_quota()
    assert quota, (
        "the fetcher's default quota is empty, so the documented command "
        "selects no samples and stops with 'nothing selected'"
    )


def test_the_default_quota_can_reproduce_the_committed_cohort():
    """Every publication of the published cohort must be reachable."""
    quota = _default_quota()
    published = _published_cohort()
    assert published, "the artifact reports no cohort sizes"

    for pub, n in sorted(published.items()):
        assert pub in quota, (
            f"publication {pub} contributed {n} samples to the committed "
            f"cohort but is not in the default quota ({quota})"
        )
        assert quota[pub] >= n, (
            f"the default quota caps publication {pub} at {quota[pub]} "
            f"samples, but the committed cohort has {n}: the published "
            "numbers would not reproduce"
        )


def test_the_quota_reader_reads_the_script_s_own_format():
    """The guard's own parser must read the format the script writes."""
    quota = _default_quota()
    assert quota, quota
    assert all(isinstance(k, int) and isinstance(v, int)
               for k, v in quota.items()), quota


def test_the_fetcher_exposes_its_parser():
    """The script must keep the parser reachable, or the guard is vacuous."""
    parser = _fetcher().build_parser()
    assert parser.get_default("quota") is not None
    assert parser.get_default("publications") == [9, 7, 1]
