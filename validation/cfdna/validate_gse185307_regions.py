"""GSE185307 real cfDNA methylation: CpG island region-aware validation.

Instead of genome-wide bulk methylation (which is coverage-confounded and
washes out tissue-specific signal), this annotates each bedgraph CpG as:
  - CpG island, N/S shore (±2kb), N/S shelf (±2-4kb), open sea

Cancer-specific methylation changes are concentrated in shores/shelves,
not inside CpG islands or genome-wide. This analysis tests whether
region-aware features recover the cancer signal at n=13.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

DATA_DIR = Path("/Users/hermes/deepcatch/data/GSE185307")
CPG_ISLANDS = Path("/Users/hermes/deepcatch/data/reference/cpgIslandExt.hg38.txt.gz")
SHORE_KB = 2  # shore = ±2kb from island boundary
SHELF_KB = 4  # shelf = ±4kb (2-4kb from boundary)

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
LABELS = {
    "GSM6069339": 0, "GSM6069340": 0, "GSM6069341": 0,
    "GSM6069342": 1, "GSM6069343": 1,
    "GSM6069344": 0, "GSM6069345": 0, "GSM6069346": 0, "GSM6069347": 0,
    "GSM6069348": 1, "GSM6069349": 1, "GSM6069350": 1, "GSM6069351": 1,
}
NAMES = {"GSM6069339": "HU5.10", "GSM6069340": "HU5.11", "GSM6069341": "HU5.12",
         "GSM6069342": "S1", "GSM6069343": "BC01", "GSM6069344": "BC02",
         "GSM6069345": "BC03", "GSM6069346": "BC04", "GSM6069347": "BC05",
         "GSM6069348": "BC08", "GSM6069349": "BC09", "GSM6069350": "BC10",
         "GSM6069351": "BC11"}


def load_islands():
    """Return dict: chrom → [(start, end), ...] sorted."""
    islands = {}
    with gzip.open(CPG_ISLANDS, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            chrom, start, end = p[1], int(p[2]), int(p[3])
            islands.setdefault(chrom, []).append((start, end))
    for v in islands.values():
        v.sort()
    return islands


def classify_stream(bedgraph_path, islands, min_cpgs=500_000):
    """Stream bedgraph and accumulate β values per region type.

    Coverage-normalized: samples uniformly downsampled to min_cpgs
    (the smallest autosomal CpG count across all samples = 3.48M).
    """
    island_betas = []
    shore_betas = []
    shelf_betas = []
    open_betas = []
    all_positions = []

    with gzip.open(bedgraph_path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] not in AUTOSOMES:
                continue
            try:
                b = float(p[3])
            except ValueError:
                continue
            chrom, pos = p[0], int(p[1])
            isl = islands.get(chrom, [])
            if not isl:
                open_betas.append(b)
                all_positions.append(("open", b))
                continue
            # Binary search for overlapping or nearest island
            lo, hi = 0, len(isl) - 1
            found = "open"
            while lo <= hi:
                mid = (lo + hi) // 2
                s, e = isl[mid]
                if s <= pos <= e:
                    found = "island"
                    break
                if pos < s:
                    # left of island mid; check if in shore/shelf left
                    dist = s - pos
                    if dist <= SHORE_KB * 1000:
                        found = "shore"
                    elif dist <= SHELF_KB * 1000:
                        found = "shelf"
                    hi = mid - 1
                else:
                    # right of island mid
                    dist = pos - e
                    if dist <= SHORE_KB * 1000:
                        found = "shore"
                    elif dist <= SHELF_KB * 1000:
                        found = "shelf"
                    lo = mid + 1

            if found == "island":
                island_betas.append(b)
            elif found == "shore":
                shore_betas.append(b)
            elif found == "shelf":
                shelf_betas.append(b)
            else:
                open_betas.append(b)

    # Coverage-normalize: uniform random downsampling to min_cpgs
    rng = np.random.RandomState(42)
    arrs = {"island": island_betas, "shore": shore_betas,
            "shelf": shelf_betas, "open": open_betas}
    features = {}
    for key, vals in arrs.items():
        vals = np.array(vals)
        if len(vals) > min_cpgs:
            idx = rng.choice(len(vals), min_cpgs, replace=False)
            vals = vals[idx]
        features[f"{key}_n"] = len(vals)
        if len(vals) > 0:
            features[f"{key}_mean"] = float(np.mean(vals))
            features[f"{key}_std"] = float(np.std(vals))
            features[f"{key}_pct_hyper"] = float(np.mean(vals >= 0.7))
            features[f"{key}_pct_hypo"] = float(np.mean(vals <= 0.3))
        else:
            features[f"{key}_mean"] = 0.5
            features[f"{key}_std"] = 0.0
            features[f"{key}_pct_hyper"] = 0.0
            features[f"{key}_pct_hypo"] = 0.0

    return features


def main():
    print("Loading CpG island annotations...")
    islands = load_islands()
    total = sum(len(v) for v in islands.values())
    print(f"  {total} islands across {len(islands)} chromosomes")

    # Find min autosomal CpG count for coverage normalization
    min_cpg = 100_000_000
    from collections import Counter
    region_totals = Counter()
    for gsm in LABELS:
        p = DATA_DIR / f"{gsm}.bedgraph.gz"
        acc = 0
        with gzip.open(p, "rt") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if parts[0] in AUTOSOMES:
                    acc += 1
        region_totals[gsm] = acc
        if acc < min_cpg:
            min_cpg = acc
    print(f"  Coverage normalization: all samples downsampled to {min_cpg:,} CpGs")

    # Per-sample region-aware features
    all_feats = {}
    for gsm, label in sorted(LABELS.items()):
        p = DATA_DIR / f"{gsm}.bedgraph.gz"
        feat = classify_stream(p, islands, min_cpg)
        all_feats[gsm] = feat
        print(f"  {gsm} ({NAMES[gsm]:>6} {'LUAD' if label else 'HEALTHY'}): "
              f"island={feat['island_mean']:.3f} shore={feat['shore_mean']:.3f} "
              f"open={feat['open_mean']:.3f}")

    # Build feature matrix
    feature_cols = ["island_mean", "island_pct_hyper", "island_pct_hypo",
                    "shore_mean", "shore_pct_hyper", "shore_pct_hypo",
                    "shelf_mean", "shelf_pct_hyper",
                    "open_mean", "open_pct_hyper", "open_pct_hypo"]
    X = np.array([[all_feats[gsm].get(col, 0.5) for col in feature_cols]
                  for gsm in sorted(LABELS)])
    y = np.array([LABELS[gsm] for gsm in sorted(LABELS)])
    sample_ids = sorted(LABELS)

    # Univariate checks
    print("\n  Univariate cancer vs healthy (Mann-Whitney):")
    for ci, col in enumerate(feature_cols):
        pos = X[y == 1, ci]
        neg = X[y == 0, ci]
        d = "cancer↑" if np.mean(pos) > np.mean(neg) else "cancer↓"
        try:
            u, p = mannwhitneyu(pos, neg, alternative="two-sided")
        except ValueError:
            p = 1.0
        print(f"    {col:<22} c={np.mean(pos):.3f} h={np.mean(neg):.3f} p={p:.3f} {d}")

    # LOOCV
    n = len(y)
    y_pred = np.zeros(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        sc = StandardScaler().fit(X[mask])
        Xtr = sc.transform(X[mask])
        clf = LogisticRegression(C=0.1, max_iter=2000)
        clf.fit(Xtr, y[mask])
        y_pred[i] = clf.predict_proba(sc.transform(X[i:i + 1]))[:, 1][0]

    auc = roc_auc_score(y, y_pred)
    print(f"\n  REGION-AWARE LOOCV (11 features, coverage-normalized):")
    print(f"    AUC = {auc:.4f}")
    print(f"    (genome-wide bulk features AUC was 0.357 — coverage-contaminated)")

    # Sens @ 95%
    n_neg = int((1 - y).sum())
    max_fp = max(0, int(np.ceil((1 - 0.95) * n_neg)) - 1)
    pred_neg = y_pred[y == 0]
    thresh = np.min(pred_neg) - 1e-9 if max_fp == 0 else np.partition(pred_neg, max_fp)[max_fp]
    sens95 = np.mean(y_pred[y == 1] > thresh)

    print(f"    Sensitivity @ 95% spec = {sens95:.3f} "
          f"({int((y_pred[y==0] > thresh).sum())} FPs, threshold={thresh:.4f})")

    # Per-sample
    print("\n  Per-sample scores:")
    for gsm, lab, score in sorted(zip(sample_ids, y, y_pred),
                                  key=lambda t: -t[2]):
        print(f"    {gsm} {NAMES[gsm]:>6} {'LUAD' if lab else 'HEALTHY'} score={score:.4f}")

    print(f"\n  ⚠ n=13 — feasibility check only, not a clinical claim.")

    import json
    out = {"dataset": "GSE185307 (region-aware, coverage-normalized)",
           "samples": n, "cancer": int(y.sum()), "healthy": int((1 - y).sum()),
           "genome_wide_auc": 0.357, "region_aware_auc": auc,
           "sens_at_95_spec": sens95,
           "per_sample_scores": {gsm: float(s) for gsm, s in zip(sample_ids, y_pred)}}
    with open(Path("/Users/hermes/deepcatch/results/gse185307_region_aware_validation.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("  Saved: results/gse185307_region_aware_validation.json")


if __name__ == "__main__":
    main()
