#!/usr/bin/env python3
"""Download per-case MAF for the 20 TCGA-LUAD cohort patients and extract
Uniprot + amino acid info.

This is a one-shot enrichment step. The result is saved to
``results/cohort_with_amino_acids.jsonl`` which the joint CADD+AM script
consumes directly.

Honest scope:
    20 patients only — the ones already in ``cadd_mutations_full.json``.
    Missense SNVs only (AlphaMissense is missense-only).
    We do NOT touch non-missense or non-SNV mutations.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import requests

ROOT = Path("/Users/hermes/deepcatch")
COHORT = ROOT / "results" / "cadd_mutations_full.json"
OUT_PATH = ROOT / "results" / "cohort_with_amino_acids.jsonl"
CACHE_MAF = ROOT / "validation" / "tcga" / "tcga_cache" / "downloaded_mafs"
CACHE_MAF.mkdir(parents=True, exist_ok=True)

GDC = "https://api.gdc.cancer.gov"


def _gdc_post(url: str, body: dict, timeout: int = 30) -> dict:
    r = requests.post(url, json=body, timeout=timeout)
    r.raise_for_status()
    return r.json()


def find_maf_for_case(submitter_id: str) -> Optional[dict]:
    """Return the file record for the per-case MuTect2 MAF (.maf.gz)."""
    data = _gdc_post(f"{GDC}/files", {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.submitter_id", "value": submitter_id}},
                {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}},
                {"op": "=", "content": {"field": "data_format", "value": "MAF"}},
            ],
        },
        "fields": "file_id,file_name,file_size",
        "size": 5,
    })
    hits = data.get("data", {}).get("hits", [])
    # Prefer the per-sample MuTect2 MAF (not protected/consolidated)
    for h in hits:
        if h.get("file_name", "").endswith(".maf.gz") and "wxs.MuTect2" in h.get("file_name", ""):
            return h
    return hits[0] if hits else None


def download_maf(file_id: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    url = f"{GDC}/data/{file_id}"
    r = requests.get(url, timeout=60, stream=True)
    r.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 16):
            f.write(chunk)
    return dest


_AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U",
}


def parse_maf_with_aa(path: Path) -> List[Dict]:
    rows: List[Dict] = []
    with gzip.open(path, "rt", errors="replace") as fp:
        text = fp.read()
    lines = text.splitlines()
    hdr_idx = None
    for i, ln in enumerate(lines):
        if ln.startswith("Hugo_Symbol"):
            hdr_idx = i
            break
    if hdr_idx is None:
        return rows
    header = lines[hdr_idx].split("\t")
    col = {n: header.index(n) for n in (
        "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome",
        "Start_Position", "End_Position", "Reference_Allele",
        "Tumor_Seq_Allele2", "Variant_Classification",
        "SWISSPROT", "Protein_position", "HGVSp",
        "Amino_acids", "HGVSp_Short",
    ) if n in header}
    for ln in lines[hdr_idx + 1:]:
        if not ln.strip():
            continue
        f_ = ln.split("\t")
        try:
            chrom = f_[col["Chromosome"]]
        except (KeyError, IndexError):
            continue
        vc = f_[col.get("Variant_Classification", -1)] if "Variant_Classification" in col else ""
        if vc != "Missense_Mutation":
            continue
        ref = f_[col["Reference_Allele"]]
        alt = f_[col["Tumor_Seq_Allele2"]]
        if len(ref) != 1 or len(alt) != 1:
            continue
        sw = f_[col["SWISSPROT"]].strip().split(".")[0]
        if not sw:
            continue
        pp_raw = f_[col["Protein_position"]].strip()
        try:
            pp = int(pp_raw.split("/")[0].split("-")[0].strip())
        except ValueError:
            continue
        # AA from Amino_acids "R/L" or HGVSp_Short "p.R247L"
        ref_aa = alt_aa = None
        if "Amino_acids" in col:
            aa = f_[col["Amino_acids"]].strip()
            if "/" in aa:
                a, b = aa.split("/", 1)
                if len(a) == 1 and len(b) == 1:
                    ref_aa, alt_aa = a, b
        if (ref_aa is None or alt_aa is None) and "HGVSp_Short" in col:
            m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", f_[col["HGVSp_Short"]].strip())
            if m:
                ref_aa = _AA3.get(m.group(1), m.group(1)[0])
                alt_aa = _AA3.get(m.group(3), m.group(3)[0])
        if (ref_aa is None or alt_aa is None) and "HGVSp" in col:
            m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", f_[col["HGVSp"]].strip())
            if m:
                ref_aa = _AA3.get(m.group(1), m.group(1)[0])
                alt_aa = _AA3.get(m.group(3), m.group(3)[0])
        rows.append({
            "sample": f_[col["Tumor_Sample_Barcode"]][:12],
            "chrom": chrom,
            "pos": int(f_[col["Start_Position"]]),
            "ref": ref,
            "alt": alt,
            "uniprot": sw,
            "prot_pos": pp,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "variant_class": vc,
        })
    return rows


def main():
    cohort = json.load(open(COHORT))["cohort_20_patients"]
    samples = sorted(cohort.keys())
    print(f"Cohort: {len(samples)} patients, "
          f"{sum(len(v) for v in cohort.values())} total mutations")

    all_rows = []
    for i, s in enumerate(samples):
        cached = CACHE_MAF / f"{s}.maf.gz"
        if cached.exists() and cached.stat().st_size > 1000:
            print(f"[{i+1}/{len(samples)}] {s}: using cached {cached.name}")
            rows = parse_maf_with_aa(cached)
            print(f"   → {len(rows)} missense SNVs with AA info")
            all_rows.extend(rows)
            continue
        try:
            print(f"[{i+1}/{len(samples)}] {s}: locating MAF…", end=" ", flush=True)
            rec = find_maf_for_case(s)
            if rec is None:
                print("no MAF file found")
                continue
            print(f"downloading {rec.get('file_name')}…", end=" ", flush=True)
            download_maf(rec["file_id"], cached)
            rows = parse_maf_with_aa(cached)
            print(f"→ {len(rows)} missense SNVs")
            all_rows.extend(rows)
        except Exception as e:
            print(f"error: {e}")
        time.sleep(0.2)

    with open(OUT_PATH, "w") as f:
        for r in all_rows:
            f.write(json.dumps(r) + "\n")
    print(f"\nWrote {len(all_rows)} rows to {OUT_PATH}")


if __name__ == "__main__":
    main()
