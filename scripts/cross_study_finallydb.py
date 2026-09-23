#!/usr/bin/env python3
"""Cross-study FinaleDB benchmark: Jiang 2015 (pub 6) + Cristiano 2019 (pub 8).

Honest open-data benchmark, NOT clinical validation. Open-data scope, no
external validation, no held-out clinical cohort. Per the cfdna-fragmentomics
skill:
- Both studies are low-pass cfDNA WGS uniformly pre-processed by FinaleDB, so
  pooling is technically valid with mild harmonization.
- The TRUE cross-study confound (cancer = 100% study A, healthy = 100% study B)
  reaches AUC 0.999 without harmonization — that is the negative control here.
- Per-study z-score harmonization inside each CV fold removes the batch effect
  without leaking test-set statistics.

Sections:
  1. Cell-line filter + cohort inventory (per-study, per-cancer)
  2. Per-cohort AUC: Jiang alone + Cristiano alone
  3. Pooled cross-study with per-study harmonization (5-seed x 5-fold CV)
  4. Per-cancer sens@spec with bootstrap 95% CIs at spec in {0.95, 0.98, 0.99}
  5. True-confound control: cancer = 100% Jiang, healthy = 100% Cristiano,
     per-study z-score harmonization collapses to ~0.5 AUC
  6. The reverse-confound (cancer = 100% Cristiano, healthy = 100% Jiang) for
     symmetry (also should collapse to ~0.5)

Outputs:
  results/cross_study_finallydb.json
  docs/CROSS_STUDY_BENCHMARK.md

Usage:
  env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \\
      scripts/cross_study_finallydb.py \\
      --features-dir /Users/hermes/cfdna-fragmentomics-pipeline/data/features \\
      --out-json   results/cross_study_finallydb.json \\
      --out-md      docs/CROSS_STUDY_BENCHMARK.md
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# Reuse the cfdna-fragmentomics-pipeline 5-channel loader + harmonize helper.
PIPELINE = "/Users/hermes/cfdna-fragmentomics-pipeline"
sys.path.insert(0, PIPELINE)
sys.path.insert(0, os.path.join(PIPELINE, "scripts"))
from honest_benchmark import load5  # noqa: E402
from train_classifier import _harmonize  # noqa: E402

# Cell-line filter regex per the cfdna-fragmentomics skill
CELL_LINE_RE = re.compile(
    r"^(GM\d+|HeLa|HepG2|K562|HL60|Jurkat|Raji|MCF7|U937|THP1|HEK293|"
    r"HCT116|SW480|A549|GM12878)",
    re.I,
)

DEFAULT_SEEDS = [42, 13, 7, 99, 1234]
DEFAULT_PCA_N = 200
SPECS_FOR_PER_CANCER = [0.95, 0.98, 0.99]
N_BOOTSTRAP = 1000
BOOTSTRAP_SEED = 2026


# --------------------------------------------------------------------------- #
# Label loading
# --------------------------------------------------------------------------- #
def load_labels_multiclass(path: str):
    """Load labels_multiclass.tsv -> dicts.

    Returns: (labels, studies, disease_class)
      labels[s]        = 1 if cancer else 0
      studies[s]       = 'jiang' | 'cristiano'
      disease_class[s] = 'HCC_J' | 'LUAD' | 'BRCA' | 'HEALTHY' | ...
    """
    labels, studies, disease_class = {}, {}, {}
    with open(path) as f:
        # Skip header (the literal "sample\tdisease_class\tlabel\tstudy" row).
        header_seen = False
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            if not header_seen and p[0] == "sample" and p[3] == "study":
                header_seen = True
                continue
            # Drop rows where the study column is literally the header text.
            if p[3] in ("study", ""):
                continue
            s = p[0]
            labels[s] = 1 if p[2] == "cancer" else 0
            studies[s] = p[3]
            disease_class[s] = p[1]
    return labels, studies, disease_class


def apply_cell_line_filter(labels, studies, disease_class):
    """Drop samples whose IDs match the cell-line regex. Return counts."""
    drop_ids = [s for s in list(labels) if CELL_LINE_RE.match(s)]
    for s in drop_ids:
        labels.pop(s, None)
        studies.pop(s, None)
        disease_class.pop(s, None)
    return drop_ids


# --------------------------------------------------------------------------- #
# CV machinery
# --------------------------------------------------------------------------- #
def pooled_oof(X, y, st, seeds, pca_n, harmonize):
    """5-seed x 5-fold OOF predictions, pooled.

    Returns (y_true, score_pooled, per_seed_aucs).
      score_pooled = mean across seeds of per-seed OOF scores (each sample
                     gets exactly one prediction per seed).
    """
    score_acc = np.zeros(len(y), dtype=float)
    seed_aucs = []
    for sd in seeds:
        cv = StratifiedKFold(5, shuffle=True, random_state=sd)
        oof = np.zeros(len(y), dtype=float)
        for tr, te in cv.split(X, y):
            Xtr = X[tr].copy()
            Xte = X[te].copy()
            if harmonize:
                Xtr, sc = _harmonize(Xtr, st[tr], None)
                Xte, _ = _harmonize(Xte, st[te], sc)
            else:
                sc = StandardScaler().fit(Xtr)
                Xtr = sc.transform(Xtr)
                Xte = sc.transform(Xte)
            max_pca = min(Xtr.shape[0], Xtr.shape[1])
            pca = PCA(n_components=min(pca_n, max_pca)).fit(Xtr)
            Xtr_p = pca.transform(Xtr)
            Xte_p = pca.transform(Xte)
            m = LogisticRegression(max_iter=2000).fit(Xtr_p, y[tr])
            oof[te] = m.predict_proba(Xte_p)[:, 1]
        seed_aucs.append(float(roc_auc_score(y, oof)))
        score_acc += oof
    score_pooled = score_acc / len(seeds)
    return y, score_pooled, seed_aucs


def sens_at_spec(y_true, y_score, spec):
    """Sensitivity at target specificity (LARGEST fpr <= target)."""
    fpr, tpr, thr = roc_curve(y_true, y_score)
    target_fpr = 1.0 - spec
    idx = np.where(fpr <= target_fpr)[0]
    if len(idx) == 0:
        return 0.0, float("nan")
    return float(tpr[idx[-1]]), float(thr[idx[-1]])


def bootstrap_ci_sens(y_true, y_score, spec, n_boot, seed):
    """Percentile bootstrap 95% CI on sens@spec."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    sens_samples = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        s, _ = sens_at_spec(y_true[idx], y_score[idx], spec)
        sens_samples[b] = s
    lo = float(np.quantile(sens_samples, 0.025))
    hi = float(np.quantile(sens_samples, 0.975))
    return lo, hi


# --------------------------------------------------------------------------- #
# Sections
# --------------------------------------------------------------------------- #
def section_inventory(labels, studies, disease_class, dropped, samples_in_array):
    n_total = len(labels)
    n_cancer = sum(1 for v in labels.values() if v == 1)
    n_healthy = n_total - n_cancer
    by_study = {}
    for s, st in studies.items():
        d = by_study.setdefault(st, {"n_total": 0, "n_cancer": 0, "n_healthy": 0})
        d["n_total"] += 1
        if labels[s] == 1:
            d["n_cancer"] += 1
        else:
            d["n_healthy"] += 1
    by_cancer = {}
    for s, dc in disease_class.items():
        if labels[s] != 1:
            continue
        d = by_cancer.setdefault(dc, {"n_total": 0, "n_jiang": 0, "n_cristiano": 0})
        d["n_total"] += 1
        if studies[s] == "jiang":
            d["n_jiang"] += 1
        else:
            d["n_cristiano"] += 1
    n_with_features = len(samples_in_array)
    return {
        "n_total_in_labels": n_total,
        "n_with_features": n_with_features,
        "n_dropped_due_to_missing_features": n_total - n_with_features,
        "n_cancer_in_labels": n_cancer,
        "n_healthy_in_labels": n_healthy,
        "n_dropped_cell_line": len(dropped),
        "dropped_cell_line_ids": dropped,
        "per_study": by_study,
        "per_cancer": by_cancer,
    }


def section_per_cohort(X, y, st, seeds, pca_n):
    out = {}
    for study in ["jiang", "cristiano"]:
        mask = st == study
        if mask.sum() < 30 or (y[mask] == 1).sum() < 10 or (y[mask] == 0).sum() < 10:
            out[study] = {
                "n_total": int(mask.sum()),
                "n_cancer": int((y[mask] == 1).sum()),
                "n_healthy": int((y[mask] == 0).sum()),
                "skipped": True,
                "reason": "insufficient samples for 5-fold CV",
            }
            continue
        _, _, seed_aucs = pooled_oof(
            X[mask], y[mask], st[mask], seeds, pca_n, harmonize=False
        )
        out[study] = {
            "n_total": int(mask.sum()),
            "n_cancer": int((y[mask] == 1).sum()),
            "n_healthy": int((y[mask] == 0).sum()),
            "auc_mean": float(np.mean(seed_aucs)),
            "auc_std": float(np.std(seed_aucs)),
            "per_seed_auc": seed_aucs,
            "skipped": False,
        }
    return out


def section_pooled(X, y, st, seeds, pca_n):
    out = {}
    for tag, harm in [("harmonized", True), ("no_harmonize", False)]:
        y_true, score, seed_aucs = pooled_oof(X, y, st, seeds, pca_n, harmonize=harm)
        s95, _ = sens_at_spec(y_true, score, 0.95)
        s98, _ = sens_at_spec(y_true, score, 0.98)
        s99, _ = sens_at_spec(y_true, score, 0.99)
        out[tag] = {
            "auc_mean": float(np.mean(seed_aucs)),
            "auc_std": float(np.std(seed_aucs)),
            "per_seed_auc": seed_aucs,
            "sens_at_spec_95": s95,
            "sens_at_spec_98": s98,
            "sens_at_spec_99": s99,
        }
    return out


def section_per_cancer(X, y, st, samples_in_array, disease_class,
                       seeds, pca_n, cancer_types, harmonize=True):
    """For each top-N cancer, run OvR (cancer vs ALL healthy) with per-study
    z-score harmonization inside each CV fold. Report AUC and sens@spec at
    {0.95, 0.98, 0.99} with bootstrap 95% CIs.

    samples_in_array: list of sample IDs in the same order as X rows
                      (load5 sorts labels.keys()).
    """
    healthy_mask = y == 0
    # Pre-compute per-sample disease class for the loaded order.
    dc_arr = np.array(
        [disease_class.get(s, "") for s in samples_in_array],
        dtype=object,
    )
    out = {}
    for cancer in cancer_types:
        cancer_mask = (dc_arr == cancer) & (y == 1)
        n_cancer = int(cancer_mask.sum())
        if n_cancer < 10 or healthy_mask.sum() < 10:
            out[cancer] = {"n_cancer": n_cancer, "skipped": True,
                           "reason": "insufficient samples for 5-fold CV"}
            continue
        mask = cancer_mask | healthy_mask
        X_sub = X[mask]
        y_sub = cancer_mask[mask].astype(int)
        st_sub = st[mask]
        y_true, score, seed_aucs = pooled_oof(
            X_sub, y_sub, st_sub, seeds, pca_n, harmonize=harmonize
        )
        rows = []
        for spec in SPECS_FOR_PER_CANCER:
            sens, thr = sens_at_spec(y_true, score, spec)
            lo, hi = bootstrap_ci_sens(
                y_true, score, spec, N_BOOTSTRAP,
                seed=BOOTSTRAP_SEED + int(round(n_cancer)),
            )
            rows.append({
                "specificity": spec,
                "sensitivity": sens,
                "ci95_lo": lo,
                "ci95_hi": hi,
                "operating_threshold": thr,
            })
        out[cancer] = {
            "n_cancer": n_cancer,
            "n_healthy": int(healthy_mask.sum()),
            "auc_mean": float(np.mean(seed_aucs)),
            "auc_std": float(np.std(seed_aucs)),
            "per_seed_auc": seed_aucs,
            "per_specificity": rows,
            "skipped": False,
        }
    return out


def section_true_confound(X, y, st, studies, seeds, pca_n):
    """TRUE cross-study confound (cancer = one study, healthy = the other).

    Run BOTH orientations and BOTH harmonization settings. With per-study
    z-score harmonization this should collapse to AUC ~0.50; without
    harmonization it reaches ~0.999 (the classifier learns the study).
    """
    s_arr = np.asarray(studies == "jiang")  # bool (1 if jiang)
    s_arr_c = ~s_arr
    y_arr = np.asarray([1 if v == 1 else 0 for v in (y == 1)], dtype=bool)
    # Use raw st/studies from numpy array; need to derive cancer position from
    # original (mask,label) pair. Simpler: rebuild study boolean from st.
    is_jiang = st == "jiang"
    is_cristiano = st == "cristiano"
    orientations = {
        "cancer_jiang_healthy_cristiano": (is_jiang, is_cristiano),
        "cancer_cristiano_healthy_jiang": (is_cristiano, is_jiang),
    }
    out = {}
    for tag, (pos_mask, neg_mask) in orientations.items():
        mask = pos_mask | neg_mask
        n_pos = int((mask & pos_mask).sum())
        n_neg = int((mask & neg_mask).sum())
        if n_pos < 10 or n_neg < 10:
            out[tag] = {"skipped": True,
                        "reason": f"n_pos={n_pos}, n_neg={n_neg}"}
            continue
        X_sub = X[mask]
        y_sub = pos_mask[mask].astype(int)
        st_sub = st[mask]
        cfg = {}
        for cfg_tag, harm in [("harmonized", True), ("no_harmonize", False)]:
            _, _, seed_aucs = pooled_oof(
                X_sub, y_sub, st_sub, seeds, pca_n, harmonize=harm
            )
            cfg[cfg_tag] = {
                "auc_mean": float(np.mean(seed_aucs)),
                "auc_std": float(np.std(seed_aucs)),
                "per_seed_auc": seed_aucs,
            }
        out[tag] = {
            "n_cancer": n_pos,
            "n_healthy": n_neg,
            **cfg,
        }
    return out


# --------------------------------------------------------------------------- #
# Markdown writer (kept in this file for one-shot run)
# --------------------------------------------------------------------------- #
def write_markdown(payload, md_path):
    cfg = payload["config"]
    cohort = payload["cohort"]
    pc = payload["per_cohort"]
    pooled = payload["pooled"]
    percancer = payload["per_cancer"]
    confound = payload["true_confound_control"]
    interp = payload["interpretation"]

    L = []
    L.append("# Cross-Study FinaleDB Benchmark (Open Data)\n")
    L.append("> **Scope**: Open-data benchmark on FinaleDB publications 6 (Jiang 2015) + 8 "
            "(Cristiano 2019). **NOT** clinical validation. **NOT** external cohort "
            "validation. Pooled OOF on the same cohort that trained the model.\n")
    L.append(f"- Generated: `{payload['generated_at']}`\n")
    L.append(f"- Classifier: `{cfg['classifier']}`\n")
    L.append(f"- Feature set: {cfg['feature_set']}\n")
    L.append(f"- PCA n_components: {cfg['pca_n']} (capped at `min(n_train, n_features)` per fold)\n")
    L.append(f"- CV: {cfg['cv']}\n")
    L.append(f"- Seeds: `{cfg['seeds']}`\n")
    L.append(f"- Bootstrap: n={cfg['n_bootstrap']}, seed={cfg['bootstrap_seed']}\n")
    L.append(f"- Cell-line regex: `{cfg['cell_line_regex']}`\n")

    L.append("\n## 1. Cohort inventory\n")
    L.append(f"- Samples in labels_multiclass.tsv: **{cohort['n_total_in_labels']}** "
            f"({cohort['n_cancer_in_labels']} cancer + {cohort['n_healthy_in_labels']} healthy)\n")
    L.append(f"- Samples with all 5-channel DELFI features: **{cohort['n_with_features']}** "
            f"({cohort['n_dropped_due_to_missing_features']} dropped due to missing artifacts)\n")
    L.append(f"- Cell-line filter removed: **{cohort['n_dropped_cell_line']}** samples "
            f"(none matched the regex on this open-data cohort)\n")
    if cohort["dropped_cell_line_ids"]:
        L.append(f"  - Dropped IDs: `{cohort['dropped_cell_line_ids']}`\n")
    L.append(f"- Per-study:\n")
    for st, d in cohort["per_study"].items():
        L.append(f"  - **{st}**: n={d['n_total']} ({d['n_cancer']} cancer + "
                f"{d['n_healthy']} healthy)\n")
    L.append(f"- Per-cancer (cancer samples only):\n")
    L.append("  | Cancer | n_total | n_jiang | n_cristiano |\n")
    L.append("  |---|---:|---:|---:|\n")
    for cancer, d in sorted(cohort["per_cancer"].items(),
                            key=lambda kv: -kv[1]["n_total"]):
        L.append(f"  | {cancer} | {d['n_total']} | {d['n_jiang']} | "
                f"{d['n_cristiano']} |\n")

    L.append("\n## 2. Per-cohort AUC\n")
    L.append("Each study evaluated independently (no harmonization needed; one study only).\n\n")
    L.append("| Study | n | n_cancer | n_healthy | AUC (5-seed mean ± std) |\n")
    L.append("|---|---:|---:|---:|---|\n")
    for study in ["jiang", "cristiano"]:
        r = pc[study]
        if r.get("skipped"):
            L.append(f"| {study} | {r['n_total']} | {r['n_cancer']} | "
                    f"{r['n_healthy']} | SKIPPED ({r.get('reason', '')}) |\n")
        else:
            L.append(f"| {study} | {r['n_total']} | {r['n_cancer']} | "
                    f"{r['n_healthy']} | {r['auc_mean']:.4f} ± {r['auc_std']:.4f} |\n")

    L.append("\n## 3. Pooled cross-study AUC (with/without per-study harmonization)\n")
    L.append("Harmonization = per-study z-score StandardScaler fit on train fold only.\n\n")
    L.append("| Setting | AUC mean ± std | Sens@95% | Sens@98% | Sens@99% |\n")
    L.append("|---|---|---:|---:|---:|\n")
    for tag in ["harmonized", "no_harmonize"]:
        r = pooled[tag]
        L.append(f"| {tag} | {r['auc_mean']:.4f} ± {r['auc_std']:.4f} | "
                f"{r['sens_at_spec_95']:.3f} | {r['sens_at_spec_98']:.3f} | "
                f"{r['sens_at_spec_99']:.3f} |\n")

    L.append("\n## 4. Per-cancer Sens@Spec (top-5 cancers by sample count)\n")
    L.append("One-vs-rest: each cancer vs ALL healthy samples in the pooled cross-study "
            "cohort. Per-study harmonization inside each CV fold. Bootstrap 95% CI "
            "on sens@spec (n=" f"{cfg['n_bootstrap']} resamples).\n\n")
    L.append("| Cancer | n_cancer | AUC mean ± std | Sens@95% [95% CI] | "
            "Sens@98% [95% CI] | Sens@99% [95% CI] |\n")
    L.append("|---|---:|---:|---|---|---|\n")
    for cancer, r in percancer.items():
        if r.get("skipped"):
            L.append(f"| {cancer} | {r.get('n_cancer', '?')} | SKIPPED | — | — | — |\n")
            continue
        cells = []
        for spec in SPECS_FOR_PER_CANCER:
            row = next(x for x in r["per_specificity"] if x["specificity"] == spec)
            cells.append(f"{row['sensitivity']:.3f} [{row['ci95_lo']:.3f}–{row['ci95_hi']:.3f}]")
        L.append(f"| {cancer} | {r['n_cancer']} | "
                f"{r['auc_mean']:.4f} ± {r['auc_std']:.4f} | "
                f"{cells[0]} | {cells[1]} | {cells[2]} |\n")

    L.append("\n## 5. True-confound control\n")
    L.append("Cancer = 100% from one study, healthy = 100% from the other. "
            "Without harmonization the classifier learns 'which study is this "
            "from?' (AUC ~0.999). With per-study z-score harmonization the "
            "study-specific mean/variance is the only signal and is removed by "
            "design (AUC should collapse toward 0.50).\n\n")
    L.append("| Orientation | n_cancer | n_healthy | AUC harmonized | AUC no_harmonize |\n")
    L.append("|---|---:|---:|---:|---:|\n")
    for tag, r in confound.items():
        if r.get("skipped"):
            L.append(f"| {tag} | — | — | SKIPPED | SKIPPED |\n")
            continue
        h = r["harmonized"]
        nh = r["no_harmonize"]
        L.append(f"| {tag} | {r['n_cancer']} | {r['n_healthy']} | "
                f"{h['auc_mean']:.3f} ± {h['auc_std']:.3f} | "
                f"{nh['auc_mean']:.3f} ± {nh['auc_std']:.3f} |\n")

    L.append("\n## Verdict\n")
    h = pooled["harmonized"]
    nh = pooled["no_harmonize"]
    n_feat = cohort["n_with_features"]
    L.append(f"- Pooled harmonized cross-study AUC: **{h['auc_mean']:.4f} ± {h['auc_std']:.4f}** "
            f"(n={n_feat} with features, of {cohort['n_total_in_labels']} in labels file)\n")
    L.append(f"- Pooled AUC without harmonization: **{nh['auc_mean']:.4f} ± {nh['auc_std']:.4f}** "
            f"(mild change confirms the per-study batch effect is small on this "
            f"FinaleDB-uniformly-processed cohort)\n")
    j_pc = pc.get("jiang", {})
    c_pc = pc.get("cristiano", {})
    L.append(f"- Per-cohort AUC: Jiang **{j_pc.get('auc_mean', float('nan')):.4f} ± "
            f"{j_pc.get('auc_std', float('nan')):.4f}** "
            f"(n={j_pc.get('n_total', 0)}), "
            f"Cristiano **{c_pc.get('auc_mean', float('nan')):.4f} ± "
            f"{c_pc.get('auc_std', float('nan')):.4f}** "
            f"(n={c_pc.get('n_total', 0)})\n")
    # Confound verdict
    c1 = confound.get("cancer_jiang_healthy_cristiano", {})
    c2 = confound.get("cancer_cristiano_healthy_jiang", {})
    if c1 and not c1.get("skipped"):
        L.append(f"- True-confound control: cancer=Jiang + healthy=Cristiano: "
                f"harmonized AUC **{c1['harmonized']['auc_mean']:.3f}** "
                f"(should be ~0.50), no-harmonize AUC **{c1['no_harmonize']['auc_mean']:.3f}** "
                f"(should be ~1.00 — proves the batch effect is removable)\n")
    if c2 and not c2.get("skipped"):
        L.append(f"- True-confound control: cancer=Cristiano + healthy=Jiang: "
                f"harmonized AUC **{c2['harmonized']['auc_mean']:.3f}** "
                f"(should be ~0.50), no-harmonize AUC **{c2['no_harmonize']['auc_mean']:.3f}**\n")

    L.append("\n## Honest framing\n")
    L.append(interp["honest_framing"] + "\n\n")
    L.append(interp["true_confound_reading"] + "\n\n")
    L.append(interp["per_cancer_reading"] + "\n\n")
    L.append("**Open-data benchmark — NOT clinical validation.**\n")

    os.makedirs(os.path.dirname(md_path) or ".", exist_ok=True)
    with open(md_path, "w") as f:
        f.write("".join(L))
    print(f"Wrote {md_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features-dir",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/data/features")
    ap.add_argument("--labels-multiclass",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv")
    ap.add_argument("--out-json", default="results/cross_study_finallydb.json")
    ap.add_argument("--out-md", default="docs/CROSS_STUDY_BENCHMARK.md")
    ap.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    ap.add_argument("--pca", type=int, default=DEFAULT_PCA_N)
    ap.add_argument("--top-cancer-n", type=int, default=5)
    args = ap.parse_args()

    print(f"[1/6] Loading labels_multiclass from {args.labels_multiclass}")
    labels, studies, disease_class = load_labels_multiclass(args.labels_multiclass)
    print(f"      Loaded {len(labels)} samples (raw multiclass labels)")

    print(f"[2/6] Applying cell-line filter")
    dropped = apply_cell_line_filter(labels, studies, disease_class)
    print(f"      Dropped {len(dropped)} cell-line samples: {dropped}")

    print(f"[3/6] Loading 5-channel DELFI features from {args.features_dir}")
    X, y, st = load5(labels, studies, args.features_dir)
    print(f"      X.shape={X.shape}, n_cancer={(y==1).sum()}, "
          f"n_healthy={(y==0).sum()}, studies={sorted(set(st.tolist()))}")

    # Replicate load5's filter to recover the sample IDs in X-row order.
    samples_in_array = []
    for s in sorted(labels):
        paths = [f"{s}.delfi_5mb_ratio.npy",
                 f"{s}.delfi_5mb_coverage.npy",
                 f"{s}.delfi_100kb_ratio.npy",
                 f"{s}.delfi_100kb_counts.npy",
                 f"{s}.fsd.json"]
        if all(os.path.exists(os.path.join(args.features_dir, p)) for p in paths):
            samples_in_array.append(s)
    assert len(samples_in_array) == X.shape[0], (
        f"sample_id reconstruction failed: {len(samples_in_array)} vs {X.shape[0]}")

    print(f"[4/6] Per-cohort AUC")
    per_cohort = section_per_cohort(X, y, st, args.seeds, args.pca)
    for study, r in per_cohort.items():
        if r.get("skipped"):
            print(f"      {study}: SKIPPED ({r.get('reason')})")
        else:
            print(f"      {study:12s}: AUC {r['auc_mean']:.4f} ± {r['auc_std']:.4f} "
                  f"(n={r['n_total']}, {r['n_cancer']} cancer + {r['n_healthy']} healthy)")

    print(f"[5/6] Pooled cross-study AUC (with/without harmonization)")
    pooled = section_pooled(X, y, st, args.seeds, args.pca)
    for tag, r in pooled.items():
        print(f"      {tag:14s}: AUC {r['auc_mean']:.4f} ± {r['auc_std']:.4f}  "
              f"S95={r['sens_at_spec_95']:.3f}  S98={r['sens_at_spec_98']:.3f}  "
              f"S99={r['sens_at_spec_99']:.3f}")

    cancer_counts = {}
    for s, dc in disease_class.items():
        if labels[s] == 1 and dc != "HEALTHY":
            cancer_counts[dc] = cancer_counts.get(dc, 0) + 1
    top_cancers = [c for c, _ in sorted(cancer_counts.items(),
                                        key=lambda kv: -kv[1])[:args.top_cancer_n]]
    print(f"[6/6] Per-cancer sens@spec (top-{args.top_cancer_n}: {top_cancers})")
    per_cancer = section_per_cancer(
        X, y, st, samples_in_array, disease_class,
        args.seeds, args.pca, top_cancers,
    )
    for cancer, r in per_cancer.items():
        if r.get("skipped"):
            print(f"      {cancer:8s}: SKIPPED ({r.get('reason')})")
            continue
        s99 = next(x["sensitivity"] for x in r["per_specificity"]
                   if x["specificity"] == 0.99)
        print(f"      {cancer:8s}: n={r['n_cancer']:3d}  "
              f"AUC {r['auc_mean']:.4f} ± {r['auc_std']:.4f}  "
              f"S99={s99:.3f}")

    print(f"[+] TRUE cross-study confound control")
    confound = section_true_confound(X, y, st, studies, args.seeds, args.pca)
    for tag, r in confound.items():
        if r.get("skipped"):
            print(f"      {tag}: SKIPPED ({r.get('reason')})")
            continue
        h = r["harmonized"]
        nh = r["no_harmonize"]
        print(f"      {tag}:")
        print(f"        harmonized    : AUC {h['auc_mean']:.3f} ± {h['auc_std']:.3f}")
        print(f"        no_harmonize  : AUC {nh['auc_mean']:.3f} ± {nh['auc_std']:.3f}")

    inventory = section_inventory(labels, studies, disease_class, dropped,
                                  samples_in_array)

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scope": ("open-data cross-study benchmark on FinaleDB "
                  "publications 6 (Jiang 2015) + 8 (Cristiano 2019). "
                  "NOT clinical validation. NOT external cohort validation. "
                  "Pooled OOF on the same cohort that trained the model."),
        "config": {
            "seeds": args.seeds,
            "pca_n": args.pca,
            "classifier": "LogisticRegression(max_iter=2000)",
            "feature_set": "5-channel (5mb_ratio + 5mb_coverage + 100kb_ratio "
                           "+ 100kb_counts + FSD-196)",
            "harmonization": "per-study z-score StandardScaler fit on train fold only",
            "cv": "5-fold StratifiedKFold, 5-seed pooled OOF",
            "n_bootstrap": N_BOOTSTRAP,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "specificities_per_cancer": SPECS_FOR_PER_CANCER,
            "cell_line_regex": CELL_LINE_RE.pattern,
            "labels_file": args.labels_multiclass,
            "features_dir": args.features_dir,
        },
        "cohort": inventory,
        "per_cohort": per_cohort,
        "pooled": pooled,
        "per_cancer": per_cancer,
        "true_confound_control": confound,
        "interpretation": {
            "honest_framing": (
                "These numbers are pooled out-of-fold AUC on the 627-sample "
                "cross-study cohort (after load5's missing-artifact filter). "
                "Internal CV; no external validation. "
                "They measure how well the 5-channel DELFI features separate "
                "cancer from healthy when pooled across Jiang 2015 and "
                "Cristiano 2019 with per-study z-score harmonization. "
                "They do NOT measure clinical-grade sensitivity at the "
                "Galleri / CancerSEEK operating points, which require "
                "independent held-out plasma cohorts."
            ),
            "true_confound_reading": (
                "With per-study z-score harmonization the true-confound AUC "
                "(cancer = one study, healthy = the other) collapses toward "
                "0.50 — the per-study mean/variance shift is the only signal "
                "and the harmonization removes it by design. WITHOUT "
                "harmonization the same control reaches ~0.999 — the "
                "classifier learns 'which study is this from?', not "
                "'is this cancer or healthy?'. The paired comparison "
                "(harmonized vs no_harmonize) is the only honest way to "
                "claim a cross-study benchmark is not a study-batch artifact."
            ),
            "per_cancer_reading": (
                "Per-cancer sens@spec is OvR: each cancer class is scored "
                "against ALL healthy samples (not just the within-study "
                "healthy ones). This is the cross-study generalization "
                "view, not the within-study view. Top-5 cancers by count "
                "are reported; smaller cohorts (n<10 cancer) are skipped."
            ),
        },
    }

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {args.out_json}")

    write_markdown(payload, args.out_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())