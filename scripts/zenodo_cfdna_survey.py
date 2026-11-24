#!/usr/bin/env python3
"""Survey Zenodo for public cfDNA datasets that release processed features.

The goal asks for processed features -- fragment-length matrices, BED
files, methylation arrays -- that are public and re-fetchable without an
IRB or a data-use agreement. GEO was searched over 99 series and the
answer is recorded in `docs/PUBLIC_DATA_SOURCES.md` §5; this covers the
Zenodo side, which until now held two known records (the CRAG
re-packaging and FinaleMe) and no search.

What is recorded
----------------
Every query that was run, with the hit count it returned on the day, and
for the candidates that matter: the DOI, the date, the creators, the
files and their sizes, and the depositor's own description. The hit
counts are a property of the index on that day, not a stable number, and
they are recorded with their date rather than quoted as a fact about
Zenodo.

Usage:
    python scripts/zenodo_cfdna_survey.py --out-json results/zenodo_cfdna_survey.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

API = "https://zenodo.org/api/records"

# The queries that were run. They are written out in the artifact so the
# survey can be repeated: Zenodo's free-text search is an OR over the
# terms, which is why the hit counts are large and the candidates are
# picked by hand from the first page of each.
QUERIES = [
    "cell-free DNA fragmentomics",
    "cfDNA fragment length",
    "circulating tumor DNA fragmentome",
    "plasma fragmentome",
    "DELFI fragmentomics",
    "cfDNA plasma methylation",
    "cfDNA fragment length matrix",
    "liquid biopsy fragmentome",
]

# The records that carry cfDNA fragment or end-motif features -- the
# artifact types this repository can actually score.
CANDIDATE_DOIS = [
    "10.5281/zenodo.21269430",   # per-sample fragment BED files, ccRCC + healthy
    "10.5281/zenodo.15097509",   # SAUNA: per-sample fragment files + a cohort table
    "10.5281/zenodo.17253504",   # per-sample 2/3/4-mer end-motif frequency tables
    "10.5281/zenodo.23195882",   # LEAF-1 pretraining corpus manifest (fragments on HF)
    "10.5281/zenodo.23196827",   # chromatin landscapes + cfDNA fragmentation in CRC
    "10.5281/zenodo.15738112",   # genome-wide end-motif variations, immunotherapy
    "10.5281/zenodo.15109495",   # cfDNA fragmentomics, second malignant neoplasm risk
]

# Records whose cohort table is a separate file in the deposit. The table
# is where the per-sample labels live, so a survey that reads only the
# record's description misses them -- the SAUNA deposit carries a
# metadata.tsv with each sample's group, coverage and tumour fraction.
COHORT_TABLE_FILES = {
    "10.5281/zenodo.15097509": "metadata.tsv",
}

_USER_AGENT = ("DeepCatch-survey/2.3 (research; "
               "https://github.com/rollroyces/deepcatch)")


def _get(url: str) -> dict:
    request = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": _USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def _record(doi: str) -> dict:
    """One record, by its Zenodo id (the DOI's last segment)."""
    record_id = doi.rsplit(".", 1)[-1]
    return _get(f"{API}/{record_id}")


def _keywords(meta: dict) -> list[str]:
    raw = meta.get("keywords") or []
    out = []
    for item in raw:
        out.append(str(item.get("keyword", "")) if isinstance(item, dict)
                   else str(item))
    return [k for k in out if k]


def _strip_html(text: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", text or "")).strip()


def _cohort_table(doi: str, file_name: str) -> dict:
    """A deposit's per-sample cohort table, read from the deposit itself.

    The record's description says what the study is; the table says what
    each sample is. A survey that reads only the description reports a
    cohort as unlabelled when the labels are in the deposit's own
    metadata file.
    """
    record_id = doi.rsplit(".", 1)[-1]
    url = f"https://zenodo.org/records/{record_id}/files/{file_name}?download=1"
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8", errors="replace")
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return {"file": file_name, "rows": []}
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) != len(header):
            continue
        rows.append(dict(zip(header, parts)))
    return {"file": file_name, "columns": header, "rows": rows}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-json", default="results/zenodo_cfdna_survey.json")
    args = ap.parse_args()

    queries = []
    seen: dict[str, str] = {}
    for query in QUERIES:
        url = f"{API}?" + urllib.parse.urlencode(
            {"q": query, "size": 25, "type": "dataset"})
        try:
            payload = _get(url)
        except Exception as error:  # noqa: BLE001 - a failed query is recorded
            queries.append({"query": query, "error": str(error)})
            time.sleep(2)
            continue
        hits = payload.get("hits", {})
        queries.append({
            "query": query,
            "dataset_hits": int(hits.get("total", 0)),
            "returned": len(hits.get("hits", [])),
        })
        for hit in hits.get("hits", []):
            doi = hit.get("doi")
            if doi:
                seen.setdefault(doi, (hit.get("metadata", {})
                                      .get("title") or "").strip())
        time.sleep(1.5)

    candidates = []
    for doi in CANDIDATE_DOIS:
        try:
            record = _record(doi)
        except Exception as error:  # noqa: BLE001
            candidates.append({"doi": doi, "error": str(error)})
            continue
        meta = record.get("metadata", {})
        files = record.get("files") or []
        candidates.append({
            "doi": doi,
            "title": (meta.get("title") or "").strip(),
            "publication_date": (meta.get("publication_date") or "")[:10],
            "creators": [c.get("name") for c in (meta.get("creators") or [])],
            "description": _strip_html(meta.get("description", ""))[:1200],
            "keywords": _keywords(meta),
            "n_files": len(files),
            "total_bytes": sum(f.get("size", 0) for f in files),
            "files": [{"key": f.get("key"), "size": f.get("size")}
                      for f in files][:25],
            "related_identifiers": [
                {"relation": r.get("relation"), "identifier": r.get("identifier")}
                for r in (meta.get("related_identifiers") or [])],
            "html": record.get("links", {}).get("html"),
        })
        table_file = COHORT_TABLE_FILES.get(doi)
        if table_file:
            # The expected file name is recorded even when the fetch
            # fails, so a missing cohort table is visible in the artifact
            # rather than looking like a deposit that has none.
            candidates[-1]["cohort_table_files"] = [table_file]
            try:
                table = _cohort_table(doi, table_file)
            except Exception as error:  # noqa: BLE001
                table = {"file": table_file, "error": str(error)}
            candidates[-1]["cohort_table"] = table

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "surveyed": "2026-11-24",
        "method": {
            "api": API,
            "type_filter": "dataset",
            "page_size": 25,
            "note": "Zenodo's free-text search is an OR over the query terms, "
                    "so the hit counts are large; the candidates were picked by "
                    "hand from the first page of each query and are listed "
                    "here with their DOIs so a reader can check each one.",
            "queries": queries,
        },
        "records_seen_in_the_first_pages": len(seen),
        "candidates": candidates,
        "conclusion": {
            "fragment_files_usable_directly": [
                "10.5281/zenodo.21269430 -- per-sample hg38 fragment BED files, "
                "12 samples (6 ccRCC, 6 healthy), with a sample metadata table "
                "carrying pT stage",
                "10.5281/zenodo.15097509 -- per-sample hg38 fragment "
                "files, 9 samples at 15-33x coverage, with a cohort table "
                "giving each sample's group, coverage and tumour "
                "fraction: 7 CRC (tumour fractions 0.005-0.40) and 2 "
                "benign, one patient sampled at three timepoints "
                "(tumour fraction 0.29, 0.33, 0.40) -- a public "
                "MRD-style cohort with a graded cancer label"],
            "processed_end_motif_tables": [
                "10.5281/zenodo.17253504 -- per-sample 2/3/4-mer end-motif "
                "frequency tables for healthy, breast, metastatic colorectal "
                "and NSCLC cohorts"],
            "manifests_pointing_elsewhere": [
                "10.5281/zenodo.23195882 -- the LEAF-1 pretraining corpus; the "
                "fragments (211 GB) are on Hugging Face, this record holds the "
                "manifest and checksums"],
        },
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(f"wrote {args.out_json}: {len(queries)} queries, "
          f"{len(candidates)} candidate records")
    for entry in queries:
        print(f"  {entry.get('query')!r}: {entry.get('dataset_hits')} hits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
