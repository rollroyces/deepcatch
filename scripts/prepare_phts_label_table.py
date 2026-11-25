#!/usr/bin/env python3
"""Write the Cleveland Clinic PHTS label table from the survey artifact.

Why this script exists
-----------------------
The depositor's own table (``phtscfx_meta.csv``, Zenodo 15109495) is a
clinical table, not a cohort table: its sample column is called ``id``,
its class column is ``cancer`` ("Yes"/"No"), and it carries two columns
the cohort extractor never asked for and the interpretation needs --
``date_diff`` (the days between the plasma draw and the cancer
diagnosis) and ``draw_after_smn``.

The survey artifact (``results/zenodo_cfdna_survey.json``) already holds
that table, all 99 rows and all 25 columns, because the survey reads the
depositor's table rather than the record's description. This script
turns those rows into the tab-separated form the cohort scorer takes,
without dropping a column: the point is that the arm labels stay
readable next to the dates that define them.

The output is written to ``data/external_cohort/`` next to the cohort's
features. The ``data/`` tree is ignored -- the features are far too
large to version, and the results are what the repository commits -- so
the table is derived, not stored. This script is the derivation.

Usage:
    python3 scripts/prepare_phts_label_table.py
"""
from __future__ import annotations

import argparse
import json
import os

SURVEY_ARTIFACT = "results/zenodo_cfdna_survey.json"
PHTS_DOI = "10.5281/zenodo.15109495"
DEFAULT_OUT = "data/external_cohort/phts_sample_metadata.tsv"

# The depositor's sample column, and the class column the label is
# derived from. Both are renamed; every other column is kept as
# published.
_SAMPLE_COLUMN = "id"
_CLASS_COLUMN = "cancer"


def _table(survey_path: str) -> tuple[dict, list[dict]]:
    with open(survey_path, encoding="utf-8") as handle:
        survey = json.load(handle)
    candidates = [c for c in survey.get("candidates") or []
                  if c.get("doi") == PHTS_DOI]
    if not candidates:
        raise SystemExit(f"{PHTS_DOI} is not among the survey's candidates")
    table = candidates[0].get("cohort_table") or {}
    rows = table.get("rows") or []
    if not rows:
        raise SystemExit(f"{PHTS_DOI}'s survey entry carries no cohort table")
    return table, rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--survey", default=SURVEY_ARTIFACT)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    table, rows = _table(args.survey)
    columns = list(table.get("columns") or sorted(rows[0]))
    if _SAMPLE_COLUMN not in columns:
        raise SystemExit(f"the table has no {_SAMPLE_COLUMN!r} column")

    # The output's columns: the key, the derived label, and then every
    # clinical column the depositor published, in the depositor's order.
    out_columns = (["sample_id", "label"]
                   + [c for c in columns if c != _SAMPLE_COLUMN])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    written = 0
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("\t".join(out_columns) + "\n")
        for row in sorted(rows, key=lambda r: r.get(_SAMPLE_COLUMN, "")):
            sample = (row.get(_SAMPLE_COLUMN) or "").strip()
            if not sample:
                continue
            is_cancer = (row.get(_CLASS_COLUMN) or "").strip().lower() == "yes"
            values = {"sample_id": sample, "label": "1" if is_cancer else "0"}
            for column in columns:
                if column == _SAMPLE_COLUMN:
                    continue
                values[column] = str(row.get(column, "NA"))
            handle.write("\t".join(values.get(c, "NA")
                                    for c in out_columns) + "\n")
            written += 1

    print(f"wrote {args.out}: {written} samples, {len(out_columns)} columns")
    print("      (every clinical column the depositor published is kept; "
          "date_diff and draw_after_smn among them)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
