"""Parallel CADD tabix matching for large cohorts.

Uses ThreadPoolExecutor to launch many tabix calls concurrently.
"""
import json
import subprocess
import sys
import time
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/Users/hermes/deepcatch")
SNV_TSV = ROOT / "data" / "cadd" / "cadd_v1.6_gnomad_r3_snv.tsv.gz"
COHORT_PATH = ROOT / "results" / "cadd_mutations_gdc_validation.json"
OUT = ROOT / "results" / "cadd_matches_gdc_validation.json"


def lookup_one(chrom, pos, ref, alt):
    chrom_num = chrom.replace("chr", "")
    try:
        result = subprocess.run(
            ["tabix", str(SNV_TSV), f"{chrom_num}:{pos}-{pos}"],
            capture_output=True, text=True, timeout=15,
        )
    except subprocess.TimeoutExpired:
        return None
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        if fields[2] == ref and fields[3] == alt:
            try:
                return float(fields[5]), float(fields[4])
            except ValueError:
                continue
    return None


def main():
    cohort = json.load(open(COHORT_PATH))['cohort']
    all_muts = [m for v in cohort.values() for m in v]
    print(f"Total mutations: {len(all_muts)}")
    snv = [m for m in all_muts if len(m['ref']) == 1 and len(m['alt']) == 1]
    non_snv = [m for m in all_muts if not (len(m['ref']) == 1 and len(m['alt']) == 1)]
    print(f"  SNVs: {len(snv)}, non-SNV (indels): {len(non_snv)}")

    WORKERS = 32
    print(f"Matching with {WORKERS} workers...")
    s = time.time()

    matched = []
    unmatched_indices = []

    def do(i):
        return i, lookup_one(snv[i]['chrom'], snv[i]['pos'], snv[i]['ref'], snv[i]['alt'])

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(do, i) for i in range(len(snv))]
        n_done = 0
        for fut in as_completed(futs):
            i, hit = fut.result()
            if hit is not None:
                phred, raw = hit
                matched.append({**snv[i], 'cadd_phred': phred, 'cadd_raw': raw})
            else:
                unmatched_indices.append(i)
            n_done += 1
            if n_done % 5000 == 0:
                rate = n_done / (time.time() - s)
                print(f"  ...{n_done}/{len(snv)} matched={len(matched)} "
                      f"unmatched={len(unmatched_indices)} "
                      f"rate={rate:.0f}/s elapsed={time.time()-s:.1f}s")

    elapsed = time.time() - s
    print(f"\nDone in {elapsed:.1f}s")
    print(f"  Matched: {len(matched)}/{len(snv)} ({len(matched)/len(snv)*100:.1f}%)")

    # Build output
    out = {
        'snv_matched': matched,
        'snv_unmatched_snv_only': [snv[i] for i in unmatched_indices],
        'non_snv_indels': non_snv,
        'n_total': len(all_muts),
        'n_total_snv': len(snv),
        'n_matched_snv': len(matched),
        'n_unmatched_snv': len(unmatched_indices),
        'n_non_snv': len(non_snv),
        'phred_stats': (
            {'min': min(m['cadd_phred'] for m in matched),
             'median': sorted([m['cadd_phred'] for m in matched])[len(matched)//2],
             'max': max(m['cadd_phred'] for m in matched)}
            if matched else None
        ),
    }
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Saved to {OUT}")


if __name__ == '__main__':
    main()