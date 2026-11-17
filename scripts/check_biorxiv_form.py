#!/usr/bin/env python3
"""Check that the bioRxiv submission form still matches the manuscript.

`paper/BIORXIV_SUBMISSION.md` is the field-by-field form Royce fills in on
bioRxiv. It states the abstract's character count and repeats the headline
numbers, so it can go stale the moment `paper/paper.tex` changes. This
script compares the two and reports any drift.

The character count uses exactly the method the regression test
(`test/test_paper_numbers_match_artifacts.py::
test_abstract_fits_the_biorxiv_character_limit`) uses, so the form and the
test can never disagree about how long the abstract is.

It is a development helper, not part of the published package.
"""
from __future__ import annotations

import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(ROOT, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
import benchmark_numbers as _numbers  # noqa: E402
import generate_preprint_abstract as _abstract  # noqa: E402

TEX = os.path.join(ROOT, "paper", "paper.tex")
FORM = os.path.join(ROOT, "paper", "BIORXIV_SUBMISSION.md")

LIMIT = 3000
MARGIN = 2950


def artifact_value(artifact: str, path) -> float | None:
    """The value an artifact holds at ``path`` (None if it does not)."""
    full = os.path.join(ROOT, "results", artifact)
    if not os.path.isfile(full):
        return None
    with open(full, encoding="utf-8") as fh:
        node = json.load(fh)
    for key in path:
        if isinstance(node, list):
            node = node[int(key)]
        elif isinstance(node, dict):
            node = node.get(key)
        else:
            return None
        if node is None:
            return None
    try:
        return float(node)
    except (TypeError, ValueError):
        return None


def matches(written: str, value: float) -> bool:
    """Does the written number equal the artifact value at its precision?

    "0.9444" and 0.94444... match; "0.944" does not, because the
    artifact was rounded to a different precision than the form used.
    """
    clean = written.replace(",", "")
    if "." not in clean:
        return abs(float(clean) - value) < 0.5
    decimals = len(clean.split(".")[1])
    return abs(round(value, decimals) - float(clean)) < 1e-12


def abstract_text(tex: str) -> str:
    """The abstract as bioRxiv would receive it.

    This is the SAME renderer the paste-ready document
    (`docs/PREPRINT_ABSTRACT.md`, written by
    `scripts/generate_preprint_abstract.py`) and the regression test use,
    so the form, that document and the test cannot disagree about how
    long the abstract is. The previous measure deleted the maths, so it
    counted a string nobody would ever paste.
    """
    return _abstract.latex_abstract(tex)


def main() -> int:
    if not os.path.isfile(FORM):
        print(f"the submission form is not present at {FORM} "
              f"(it is gitignored: it carries the author's details). "
              f"Nothing to check.")
        return 0
    tex = open(TEX, encoding="utf-8").read()
    form = open(FORM, encoding="utf-8").read()
    flat_form = re.sub(r"\s+", " ", form)

    abstract = abstract_text(tex)
    length = len(abstract)

    print(f"paper.tex abstract: {length} characters "
          f"(bioRxiv limit {LIMIT}, the test also demands <= {MARGIN})")

    problems = []
    if length > MARGIN:
        problems.append(
            f"the abstract is {length} characters -- over the {MARGIN} the "
            f"regression test insists on, and only {LIMIT - length} of "
            f"margin under bioRxiv's limit")

    # The form states the abstract's length. If it states a different one,
    # the form is stale -- and the author would copy the wrong number.
    stated = re.findall(r"abstract is (\d{3,4})[ ,]", flat_form)
    for value in stated:
        if int(value) != length:
            problems.append(
                f"the submission form says the abstract is {value} "
                f"characters; it is {length}")

    # The form's headline table repeats numbers that live in artifacts.
    # Each one is RESOLVED against the artifact it claims to come from: a
    # form that quotes a number the artifact no longer holds would put a
    # stale figure on the bioRxiv page. The previous version of this
    # script only reported which numbers the form contained, and its
    # pairs still named the superseded external AUCs (0.9915, 0.8819),
    # so a stale form would have passed unnoticed.
    for number, artifact, path in (
        ("0.9444", "external_validation_finaledb.json",
         ("per_publication_auc", "7", "auc")),
        ("0.708", "external_validation_finaledb.json",
         ("per_publication_auc", "1", "auc")),
        ("0.7385", "external_validation_finaledb.json",
         ("pooled_external_auc_matched_publications_only", "auc")),
        ("0.9825", "fusion_regularized_smoke.json",
         ("summary", "regularized_late_fusion", "auc_mean")),
        ("0.9444", "methylation_cohort_channel.json",
         ("beta_summary_auc_loo",)),
        ("103,863", "methylation_channel_diagnosis.json",
         ("matched_locus_hmm_beta", "n_cpgs_tested")),
    ):
        if number not in flat_form:
            continue
        where = f"{artifact}:{'.'.join(str(p) for p in path)}"
        value = artifact_value(artifact, path)
        if value is None:
            problems.append(
                f"the form quotes {number}, but {artifact} holds no value "
                f"at {where}")
            continue
        if not matches(number, value):
            problems.append(
                f"the form quotes {number} as {where}, but that artifact "
                f"holds {value}")
            continue
        print(f"  form quotes {number} ({where})")

    # A form that quotes a figure the v2.3.0 regeneration superseded would
    # put a stale number on the bioRxiv page. The list of superseded
    # figures lives with the other quotable numbers, so the form and the
    # drift guard cannot disagree about what is stale.
    for value in _numbers.SUPERSEDED:
        if value in flat_form:
            problems.append(
                f"the form quotes {value}, which the v2.3.0 regeneration "
                f"superseded")

    if problems:
        print()
        for p in problems:
            print(f"  !! {p}")
        return 1
    print("the submission form agrees with the current manuscript")
    return 0


if __name__ == "__main__":
    sys.exit(main())
