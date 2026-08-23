"""Real-data validation: cfDNA methylation cancer vs healthy (GSE185307).

13 real plasma cfDNA samples (6 lung adenocarcinoma, 7 healthy controls),
Oxford Nanopore 5mC methylation bedgraphs (genome-wide, ~4.9M CpG sites each).

Design (honest for n=13):
- Per-sample genome-wide methylation features (autosomal only, to avoid
  sex-chromosome confounds)
- LOOCV logistic regression (leave-one-out is the only honest CV at n=13)
- Report AUC, sensitivity @ 95% spec, plus univariate direction checks
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path("/Users/hermes/deepcatch/data/GSE185307")

# sample → label (1 = lung adenocarcinoma, 0 = healthy control)
LABELS = {
    "GSM6069339": 0, "GSM6069340": 0, "GSM6069341": 0,   # HU005.10-12 healthy
    "GSM6069342": 1, "GSM6069343": 1,                     # S1, BC01 LUAD
    "GSM6069344": 0, "GSM6069345": 0, "GSM6069346": 0, "GSM6069347": 0,  # BC02-05 healthy
    "GSM6069348": 1, "GSM6069349": 1, "GSM6069350": 1, "GSM6069351": 1,  # BC08-11 LUAD
}

# Chromosomes to include (autosomes only — exclude X/Y to avoid sex confound)
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def load_sample_features(path: Path) -> dict:
    """Compute per-sample methylation features from a bedgraph.gz file."""
    betas_by_chr = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, beta = parts[0], parts[3]
            if chrom not in AUTOSOMES:
                continue
            try:
                b = float(beta)
            except ValueError:
                continue
            betas_by_chr.setdefault(chrom, []).append(b)

    all_betas = np.concatenate([np.array(v) for v in betas_by_chr.values()])
    n = len(all_betas)
    if n == 0:
        return {}

    # Global methylation burden & distribution
    mean_beta = float(np.mean(all_betas))
    median_beta = float(np.median(all_betas))
    pct_hyper = float(np.mean(all_betas >= 0.7))      # hypermethylated fraction
    pct_hypo = float(np.mean(all_betas <= 0.3))       # hypomethylated fraction
    beta_std = float(np.std(all_betas))

    # Methylation-pattern entropy (5 bins) — epiallelic heterogeneity proxy
    bins = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.001])
    hist, _ = np.histogram(all_betas, bins=bins)
    p = hist / max(hist.sum(), 1)
    p = p[p > 0]
    entropy = float(-(p * np.log2(p)).sum())  # max = log2(5) ≈ 2.32

    # Per-chromosome mean (heterogeneity across chromosomes)
    chr_means = [float(np.mean(v)) for v in betas_by_chr.values()]
    chr_std = float(np.std(chr_means)) if len(chr_means) > 1 else 0.0

    return {
        "n_cpgs": n,
        "mean_beta": mean_beta,
        "median_beta": median_beta,
        "pct_hyper": pct_hyper,
        "pct_hypo": pct_hypo,
        "beta_std": beta_std,
        "entropy": entropy,
        "chr_std": chr_std,
    }


def main():
    features, labels = [], []
    sample_ids = []
    for gsm, label in sorted(LABELS.items()):
        path = DATA_DIR / f"{gsm}.bedgraph.gz"
        if not path.exists():
            print(f"  ⚠ missing {path.name}")
            continue
        feat = load_sample_features(path)
        if not feat:
            print(f"  ⚠ no autosomal data in {path.name}")
            continue
        features.append(feat)
        labels.append(label)
        sample_ids.append(gsm)
        print(f"  ✓ {gsm} (label={label}): {feat['n_cpgs']:,} CpGs  "
              f"mean β={feat['mean_beta']:.4f}  entropy={feat['entropy']:.3f}")

    df = pd.DataFrame(features)
    y = np.array(labels)
    print(f"\nSamples: {len(df)}  (cancer={y.sum()}, healthy={(1 - y).sum()})")

    # ── Univariate direction checks (Mann-Whitney) ───────────────────────
    from scipy.stats import mannwhitneyu
    print("\n  Univariate cancer vs healthy (Mann-Whitney):")
    for col in ["mean_beta", "median_beta", "pct_hyper", "pct_hypo",
                "beta_std", "entropy", "chr_std"]:
        pos = df[col][y == 1].values
        neg = df[col][y == 0].values
        try:
            u, p = mannwhitneyu(pos, neg, alternative="two-sided")
        except ValueError:
            p = 1.0
        direction = "cancer↑" if np.mean(pos) > np.mean(neg) else "cancer↓"
        print(f"    {col:<11} cancer={np.mean(pos):.4f} healthy={np.mean(neg):.4f} "
              f"p={p:.3f} {direction}")

    # ── LOOCV logistic regression ────────────────────────────────────────
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score

    X = df[["mean_beta", "median_beta", "pct_hyper", "pct_hypo",
            "beta_std", "entropy", "chr_std"]].values
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
    print(f"\n  LOOCV Logistic Regression (7 features, standardized):")
    print(f"    AUC = {auc:.4f}")
    # sens @ 95% spec: with 7 healthy, 95% spec ≈ 0 FPs allowed → 0.35 FPs floor
    n_neg = int((1 - y).sum())
    max_fp = max(0, int(np.ceil((1 - 0.95) * n_neg)) - 1)
    pred_neg = y_pred[y == 0]
    if max_fp == 0:
        threshold = np.min(pred_neg) - 1e-9  # zero FPs
    else:
        threshold = np.partition(pred_neg, max_fp)[max_fp]
    sens95 = np.mean(y_pred[y == 1] > threshold)
    print(f"    Sensitivity @ 95% spec = {sens95:.3f} "
          f"(threshold={threshold:.4f}, {int((y_pred[y==0] > threshold).sum())} FPs)")

    # Per-sample scores
    print("\n  Per-sample cancer scores (higher = cancer):")
    for gsm, lab, score in sorted(zip(sample_ids, y, y_pred), key=lambda t: -t[2]):
        print(f"    {gsm}  label={lab}  score={score:.4f}")

    # ── Honest interpretation ────────────────────────────────────────────
    print(f"\n  ⚠ n=13 — this is a feasibility check, not a clinical claim.")
    print(f"  Key finding to verify against literature: cancer cfDNA is")
    print(f"  expected to show GLOBAL HYPOMETHYLATION (↓ mean β) and/or")
    print(f"  ↑ heterogeneity. Check the univariate directions above.")

    # Save results
    out = {
        "dataset": "GSE185307 (cfDNA ONT methylation, 13 samples)",
        "samples": len(df), "cancer": int(y.sum()), "healthy": int((1 - y).sum()),
        "features": list(df.columns),
        "loocv_auc": auc, "sens_at_95_spec": sens95,
        "per_sample_scores": {gsm: float(s) for gsm, s in zip(sample_ids, y_pred)},
        "univariate": {c: {"cancer_mean": float(np.mean(df[c][y == 1])),
                           "healthy_mean": float(np.mean(df[c][y == 0]))}
                       for c in ["mean_beta", "pct_hyper", "pct_hypo", "entropy"]},
    }
    import json
    with open("/Users/hermes/deepcatch/results/gse185307_methylation_validation.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\n  Saved: results/gse185307_methylation_validation.json")


if __name__ == "__main__":
    main()
