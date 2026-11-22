#!/usr/bin/env python3
"""Fetch a pure-external cfDNA cohort from FinaleDB and extract features.

Why
---
The cross-study benchmark trains and tests inside the *same* merged
cohort (FinaleDB publications 6 + 8: Jiang 2015 and Cristiano 2019).
This script builds a cohort from publications that are **never** used
for training, so the model is evaluated on data it has not seen at
all -- an external validation rather than an internal one.

Public data, no IRB required
----------------------------
FinaleDB (https://pubmed.ncbi.nlm.nih.gov/33258919/) is an open
collection of ~2,500 uniformly processed cfDNA WGS datasets. The
per-sample fragment records are served from a public S3 bucket:

    entries/<seqrun_id>/hg38/<seqrun_id>.hg38.frag.tsv.bgz

Sample metadata (disease, tissue, publication) comes from the public
REST API:

    http://finaledb.research.cchmc.org/api/v1/seqrun

Labels
------
FinaleDB's ``disease`` field is the clinical diagnosis:

- a diagnosis ending in "cancer"  -> label 1 (cancer)
- an explicitly benign diagnosis  -> label 0 (non-cancer control):
  Healthy, Hepatitis B, Cirrhosis, Inflammatory bowel disease,
  Systemic lupus erythematosus, Liver transplant

Any other diagnosis is skipped rather than guessed at.

Only ``blood plasma`` samples are used: publication 3 (Burnham 2018) is
urinary cfDNA, a different analyte, and publication 4 (Rabinowitz
2019) carries no disease metadata.

Usage
-----
    python scripts/fetch_external_cohort.py \
        --out-dir data/external_cohort \
        --publications 9 7 1 \
        --max-per-disease 120 \
        --max-mb 400

Outputs
-------
    <out-dir>/features/<sample>.delfi_*.npy / .fsd.json
    <out-dir>/manifest.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

from src.fragmentomics.finale_features import extract_sample  # noqa: E402

API = "http://finaledb.research.cchmc.org/api/v1/seqrun"
S3 = "https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org"

# Publications used for TRAINING by the cross-study benchmark. Never
# fetched here -- an external cohort must not overlap the training set.
TRAINING_PUBLICATIONS = {6, 8}

# Explicit benign diagnoses -> non-cancer control (label 0).
NON_CANCER_DISEASES = {
    "Healthy",
    "Hepatitis B",
    "Cirrhosis",
    "Inflammatory bowel disease",
    "Systemic lupus erythematosus",
    "Liver transplant",
}


def label_for(disease: str | None):
    """Return 1 (cancer), 0 (benign control) or None (unusable)."""
    if not disease:
        return None
    d = disease.strip()
    if d.lower().endswith("cancer"):
        return 1
    if d in NON_CANCER_DISEASES:
        return 0
    return None


def fetch_seqruns(api: str = API, page_size: int = 100) -> list[dict]:
    """Page through the seqrun endpoint (it uses offset/total paging)."""
    rows: list[dict] = []
    offset = 0
    while True:
        with urllib.request.urlopen(
                f"{api}?page_size={page_size}&offset={offset}", timeout=90) as u:
            d = json.load(u)
        results = d.get("results", [])
        if not results:
            break
        rows.extend(results)
        total = d.get("total") or d.get("count") or 0
        offset += len(results)
        if total and len(rows) >= total:
            break
        if offset > 10000:
            break
    return rows


def frag_url(seqrun_id: int) -> str:
    return f"{S3}/entries/EE{seqrun_id}/hg38/EE{seqrun_id}.hg38.frag.tsv.bgz"


def frag_size_mb(seqrun_id: int) -> float | None:
    """HEAD the object; returns size in MB, or None if not accessible."""
    try:
        req = urllib.request.Request(frag_url(seqrun_id), method="HEAD")
        with urllib.request.urlopen(req, timeout=45) as u:
            n = u.headers.get("Content-Length")
            return int(n) / 1e6 if n else None
    except Exception:
        return None


def select_cohort(rows: list[dict], publications: set[int],
                  quotas: dict[int, int], max_mb: float) -> list[dict]:
    """Pick the external cohort: held-out publications, plasma only.

    Selection is round-robin over (label, disease) buckets inside each
    publication, up to ``max_per_publication`` samples. Round-robin
    matters: FinaleDB publication 9 holds 950 breast + 558 prostate
    cancer samples and no controls at all, while publication 1 holds
    one or two samples of each of 18 diagnoses. Taking the first N of a
    publication would therefore yield an all-cancer, single-disease
    cohort from one study and nothing from another.

    Within a bucket, samples are ordered by seqrun id, so a re-run
    selects exactly the same cohort.
    """
    buckets: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        pub = (r.get("publication") or {}).get("id")
        if pub in TRAINING_PUBLICATIONS or pub not in publications:
            continue
        s = r.get("sample") or {}
        if s.get("tissue") != "blood plasma":
            continue
        lab = label_for(s.get("disease"))
        if lab is None:
            continue
        buckets[(pub, lab, s.get("disease"))].append(r)

    for key in buckets:
        buckets[key].sort(key=lambda r: r["id"])

    selected = []
    for pub in sorted(publications):
        max_per_publication = quotas.get(pub, 0)
        if max_per_publication <= 0:
            continue
        keys = sorted(k for k in buckets if k[0] == pub)
        # Round-robin across the buckets until the publication cap is hit.
        picked: list[dict] = []
        depth = 0
        while len(picked) < max_per_publication:
            added = False
            for key in keys:
                group = buckets[key]
                if depth < len(group):
                    r = group[depth]
                    picked.append({
                        "seqrun_id": r["id"],
                        "sample": (r.get("sample") or {}).get("name",
                                                              f"EE{r['id']}"),
                        "publication": pub,
                        "disease": key[2],
                        "label": key[1],
                        "frag_num": (r.get("fragNum") or {}).get("hg38"),
                    })
                    added = True
                    if len(picked) >= max_per_publication:
                        break
            if not added:
                break
            depth += 1
        selected.extend(picked)
    return selected


def existing_stats(features_dir: str, sample: str) -> dict | None:
    """The extraction stats for a sample already on disk, or None.

    ``extract_sample`` returns these numbers; they are also recoverable
    from the two summary files it writes, so a re-run can record a
    previously fetched sample in the manifest without downloading it
    again.
    """
    delfi_path = os.path.join(features_dir, f"{sample}.delfi.json")
    fsd_path = os.path.join(features_dir, f"{sample}.fsd.json")
    if not (os.path.exists(delfi_path) and os.path.exists(fsd_path)):
        return None
    try:
        with open(delfi_path) as f:
            delfi = json.load(f)
        with open(fsd_path) as f:
            fsd = json.load(f)
    except (OSError, ValueError):
        return None
    return {
        "sample": sample,
        "genome_short_long_ratio": delfi.get("genome_short_long_ratio"),
        "mean_short_frac": delfi.get("mean_short_frac"),
        "fsd_fragment_count": fsd.get("fragment_count"),
        "fsd_short_fraction": fsd.get("short_fraction_100_150"),
        "fsd_short_long_ratio": fsd.get("short_long_ratio"),
    }


def download_and_extract(job: dict, out_dir: str, raw_dir: str,
                         max_mb: float, keep_raw: bool,
                         skip_existing: bool = False) -> dict:
    """Download one sample's fragment file, extract its features, clean up."""
    sid = job["seqrun_id"]
    sample = job["sample"]
    raw = os.path.join(raw_dir, f"{sample}.frag.tsv.bgz")
    t0 = time.time()
    try:
        size = frag_size_mb(sid)
        if size is None:
            return {**job, "status": "s3-unavailable"}
        if size > max_mb:
            return {**job, "status": "too-large", "size_mb": round(size, 1)}
        if skip_existing:
            stats = existing_stats(os.path.join(out_dir, "features"), sample)
            if stats is not None:
                # Already extracted in an earlier run: record it without
                # downloading the fragment file again.
                return {**job, **stats, "status": "ok",
                        "size_mb": round(size, 1), "seconds": 0.0,
                        "reused_existing_features": True}
        urllib.request.urlretrieve(frag_url(sid), raw)
        stats = extract_sample(raw, os.path.join(out_dir, "features"), sample)
        if not keep_raw:
            os.remove(raw)
        return {**job, "status": "ok", "size_mb": round(size, 1),
                "seconds": round(time.time() - t0, 1), **stats}
    except Exception as exc:  # network / parse failures must not kill the run
        if os.path.exists(raw):
            os.remove(raw)
        return {**job, "status": "error", "error": f"{type(exc).__name__}: {exc}"}


def build_parser() -> argparse.ArgumentParser:
    """The argument parser.

    Split out of `main()` so that a test can read the defaults without
    running the fetch: the default quota is what makes the documented
    command reproduce the committed cohort, and that is a property of
    the script, not of a network call.
    """
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="data/external_cohort")
    ap.add_argument("--publications", type=int, nargs="+", default=[9, 7, 1])
    ap.add_argument("--quota", default="9=240,7=31,1=240",
                    help="per-publication sample budget, e.g. '9=240,7=31,"
                         "1=240'. The default is the quota that reproduces "
                         "the committed cohort (results/"
                         "external_validation_finaledb.json): 240 + 31 + 60 "
                         "samples selected, 322 extracted after the "
                         "too-large or unreachable ones are dropped. A "
                         "smaller cap for a publication yields a DIFFERENT "
                         "cohort, so the published numbers do not reproduce "
                         "with it. Publications not listed are skipped.")
    ap.add_argument("--max-per-publication", type=int, default=240,
                    help="fallback budget for publications not in --quota")
    ap.add_argument("--max-mb", type=float, default=400.0,
                    help="skip samples whose fragment file exceeds this size")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        help=("reuse the features of samples already extracted into "
              "--out-dir/features instead of downloading their fragment "
              "files again. Use this to enlarge an existing cohort: the "
              "manifest still describes every sample, and only the new "
              "ones cost bandwidth."))
    return ap


def main() -> int:
    args = build_parser().parse_args()

    out_dir = os.path.abspath(args.out_dir)
    raw_dir = os.path.join(out_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "features"), exist_ok=True)

    print(f"[1/4] querying FinaleDB seqruns ...", flush=True)
    rows = fetch_seqruns()
    print(f"      {len(rows)} seqruns")

    print(f"[2/4] selecting held-out plasma cohort from publications "
          f"{sorted(args.publications)} ...", flush=True)
    quotas: dict[int, int] = {}
    for part in args.quota.split(","):
        part = part.strip()
        if not part:
            continue
        pub, _, n = part.partition("=")
        quotas[int(pub)] = int(n)
    for pub in args.publications:
        quotas.setdefault(pub, args.max_per_publication)
    print(f"      per-publication quota: {quotas}")
    cohort = select_cohort(rows, set(args.publications), quotas, args.max_mb)
    counts = Counter((c["publication"], c["disease"], c["label"]) for c in cohort)
    for (pub, dis, lab), n in sorted(counts.items()):
        print(f"      pub {pub:<3} {dis:<34} label={lab} n={n}")
    print(f"      total selected: {len(cohort)}")
    if not cohort:
        print("nothing selected")
        return 1

    print(f"[3/4] downloading + extracting ({args.workers} workers) ...",
          flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        results = list(ex.map(
            lambda j: download_and_extract(j, out_dir, raw_dir, args.max_mb,
                                           args.keep_raw, args.skip_existing),
            cohort))
    ok = [r for r in results if r["status"] == "ok"]
    print(f"      extracted {len(ok)}/{len(cohort)} in {time.time() - t0:.0f}s")
    for status, n in Counter(r["status"] for r in results).items():
        if status != "ok":
            print(f"      {status}: {n}")

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "source": "FinaleDB (public API + public S3), no IRB required",
        "training_publications_excluded": sorted(TRAINING_PUBLICATIONS),
        "publications": sorted(args.publications),
        "n_extracted": len(ok),
        "n_selected": len(cohort),
        "samples": ok,
    }
    mpath = os.path.join(out_dir, "manifest.json")
    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[4/4] wrote {mpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
