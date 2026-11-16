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

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, "paper", "paper.tex")
FORM = os.path.join(ROOT, "paper", "BIORXIV_SUBMISSION.md")

LIMIT = 3000
MARGIN = 2950


def abstract_text(tex: str) -> str:
    """The abstract as bioRxiv would receive it.

    Same normalisation as the regression test: drop one level of
    ``\\cmd{...}`` argument braces, then the remaining commands, then the
    math, then collapse whitespace.
    """
    match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    if not match:
        raise ValueError("paper.tex has no abstract")
    plain = re.sub(r"\\[a-zA-Z]+\{([^}]*)\}", r"\1", match.group(1))
    plain = re.sub(r"\\[a-zA-Z]+", "", plain)
    plain = re.sub(r"\$[^$]*\$", "0", plain)
    return " ".join(plain.split())


def main() -> int:
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
    # If one of them has moved, the form would put a stale number on the
    # bioRxiv page.
    for number, artifact, key in (
        ("0.9915", "external_validation_finaledb.json",
         "per_publication_auc.7.auc"),
        ("0.8819", "external_validation_finaledb.json", "pooled"),
        ("0.9825", "fusion_regularized_smoke.json", "regularized_late_fusion"),
        ("0.9444", "methylation_channel_diagnosis.json", "beta_summary_auc_loo"),
        ("103,863", "methylation_channel_diagnosis.json", "n_cpgs_tested"),
    ):
        if number in flat_form:
            print(f"  form quotes {number} ({artifact}:{key})")

    if problems:
        print()
        for p in problems:
            print(f"  !! {p}")
        return 1
    print("the submission form agrees with the current manuscript")
    return 0


if __name__ == "__main__":
    sys.exit(main())
