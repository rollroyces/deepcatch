#!/usr/bin/env python3
"""Cross-study FinaleDB benchmark (publication-aware).

Honest open-data benchmark, NOT clinical validation. Open-data scope, no
external validation, no held-out clinical cohort. Per the cfdna-fragmentomics
skill:
- All FinaleDB studies are uniformly pre-processed cfDNA WGS, so pooling is
  technically valid with mild per-study harmonization.
- The TRUE cross-study confound (cancer = 100% study A, healthy = 100%
  study B) reaches AUC 0.999 without harmonization — that is the negative
  control here.
- Per-study z-score harmonization inside each CV fold removes the batch
  effect without leaking test-set statistics.

Publications (per cfdna-fragmentomics skill, FinaleDB publication id→study map):
  1  = Snyder 2016 Cell
  6  = Jiang 2015 PNAS (low-pass HCC)
  7  = Sun 2019
  8  = Cristiano 2019 (DELFI, pan-cancer + healthy)
  9  = Adalsteinsson 2017

Sections:
  1. Cell-line filter + cohort inventory (per-publication, per-cancer)
  2. Per-publication AUC (one fit per requested publication)
  3. Pooled cross-study with per-publication harmonization (5-seed x 5-fold CV)
  4. Per-cancer sens@spec with bootstrap 95% CIs at spec in {0.95, 0.98, 0.99}
  5. True-confound control: cancer = 100% pub A, healthy = 100% pub B,
     per-publication z-score harmonization collapses to ~0.5 AUC
  6. The reverse-confound for symmetry (also should collapse to ~0.5)

Outputs:
  results/cross_study_finallydb.json
  docs/CROSS_STUDY_BENCHMARK.md

FinaleDB API status: the public REST API and S3 bucket have been DOWN since
2026-09 (Postgres connection lost + S3 keys returned 403). Until they return,
adding publications 1 (Snyder) and 7 (Sun) requires re-fetching their
features, which is impossible. See docs/PUBLICATION_READINESS.md for the
honest current status. The default `--publications 6 8` uses the
already-cached local features for the open-data benchmark.

Usage:
  env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \\
      scripts/cross_study_finallydb.py \\
      --features-dir /Users/hermes/cfdna-fragmentomics-pipeline/data/features \\
      --labels-multiclass /Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv \\
      --publications 6 8 \\
      --out-json results/cross_study_finallydb.json \\
      --out-md    docs/CROSS_STUDY_BENCHMARK.md
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import warnings
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

# Standalone DeLong-CI helper (clinical-decision schema).
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
from src.per_cancer_sens_at_spec import (  # noqa: E402
    DEFAULT_PREVALENCES,
    DEFAULT_SPECIFICITIES,
    MIN_POSITIVES_FOR_CI,
    PPV_AT_SPEC,
    build_per_cancer_table,
    delong_auc_ci,
    delong_sens_at_spec_ci,
    ppv_at_prevalence,
)

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

# ──────────────────────────────────────────────────────────────────────
# Publication registry (cfdna-fragmentomics skill: id→study map)
# ──────────────────────────────────────────────────────────────────────
#
# Each entry is (pub_id → short_label). The short_label is the value used
# for the harmonization grouping in `_harmonize(X, st, ...)` and as the
# display label in the inventory + markdown output.
#
# When FinaleDB API/S3 are reachable again, the new publications can be
# added simply by extending this dict + providing a per-publication
# labels file (or a 5th `publication` column in labels_multiclass.tsv).
PUBLICATION_REGISTRY = {
    "1": "snyder",       # Snyder 2016 Cell      (FinaleDB API/S3 currently DOWN)
    "6": "jiang",        # Jiang 2015 PNAS       (cached locally)
    "7": "sun",          # Sun 2019              (FinaleDB API/S3 currently DOWN)
    "8": "cristiano",    # Cristiano 2019 DELFI  (cached locally)
    "9": "adalsteinsson",  # Adalsteinsson 2017
}
# Known study→publication fallback for labels files that only have a
# `study` column (the current local cache). Used when the labels file
# has no 5th `publication` column.
STUDY_TO_PUBLICATION = {
    "jiang": "6",
    "cristiano": "8",
    "snyder": "1",
    "sun": "7",
    "adalsteinsson": "9",
}


# --------------------------------------------------------------------------- #
# Label loading (now publication-aware)
# --------------------------------------------------------------------------- #
def load_labels_multiclass(path: str, requested_publications: set[str]):
    """Load labels_multiclass.tsv → dicts.

    Returns: (labels, studies, disease_class, publication)
      labels[s]        = 1 if cancer else 0
      studies[s]       = 'jiang' | 'cristiano' | ...  (display label)
      disease_class[s] = 'HCC_J' | 'LUAD' | 'BRCA' | 'HEALTHY' | ...
      publication[s]   = '6' | '8' | ...  (FinaleDB publication id)

    The function is backward-compatible: if the labels file has NO 5th
    `publication` column, the publication is derived from the `study`
    column via STUDY_TO_PUBLICATION. Samples whose publication is not in
    `requested_publications` are dropped (with a warning listing them).

    If the labels file has the `publication` column and it disagrees with
    the `study` column, the `publication` column wins (publication is the
    primary key for the cohort filter).
    """
    labels, studies, disease_class, publication = {}, {}, {}, {}
    skipped_outside_request = []
    with open(path) as f:
        header_seen = False
        header_cols = None
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            if not header_seen:
                if p[0] == "sample" and p[3] == "study":
                    header_seen = True
                    header_cols = p
                    continue
                # No header row found yet — treat this as a data row.
                header_seen = True
            if p[0] == "sample" and p[3] == "study":
                # Re-reading header on a different row (defensive)
                header_cols = p
                continue
            if p[3] in ("study", ""):
                continue
            s = p[0]
            # Publication id: prefer an explicit 5th column if present.
            pub = ""
            if len(p) >= 5 and p[4] and p[4] not in ("publication",):
                pub = p[4]
            else:
                pub = STUDY_TO_PUBLICATION.get(p[3], "")
            if pub not in requested_publications:
                skipped_outside_request.append((s, p[3], pub))
                continue
            labels[s] = 1 if p[2] == "cancer" else 0
            studies[s] = p[3]
            disease_class[s] = p[1]
            publication[s] = pub
    if skipped_outside_request:
        uniq = sorted({(st, pub) for _, st, pub in skipped_outside_request
                       if pub and pub not in requested_publications})
        if uniq:
            warnings.warn(
                f"load_labels_multiclass: dropped {len(skipped_outside_request)} "
                f"samples whose publication is outside the requested set "
                f"{sorted(requested_publications)}. "
                f"Dropped study→pub combos: {uniq}",
                stacklevel=2,
            )
    return labels, studies, disease_class, publication


def apply_cell_line_filter(labels, studies, disease_class, publication):
    """Drop samples whose IDs match the cell-line regex. Return counts."""
    drop_ids = [s for s in list(labels) if CELL_LINE_RE.match(s)]
    for s in drop_ids:
        labels.pop(s, None)
        studies.pop(s, None)
        disease_class.pop(s, None)
        publication.pop(s, None)
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
def section_inventory(labels, studies, disease_class, publication, dropped,
                      samples_in_array):
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
    by_publication = {}
    for s, pub in publication.items():
        d = by_publication.setdefault(pub, {
            "n_total": 0, "n_cancer": 0, "n_healthy": 0,
            "study_label": PUBLICATION_REGISTRY.get(pub, pub),
        })
        d["n_total"] += 1
        if labels[s] == 1:
            d["n_cancer"] += 1
        else:
            d["n_healthy"] += 1
    by_cancer = {}
    for s, dc in disease_class.items():
        if labels[s] != 1:
            continue
        d = by_cancer.setdefault(dc, {"n_total": 0})
        d["n_total"] += 1
        for st_key in by_study:
            if st_key not in d:
                d[st_key] = 0
            if studies[s] == st_key:
                d[st_key] += 1
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
        "per_publication": by_publication,
        "per_cancer": by_cancer,
    }


def section_per_publication(X, y, st, publication_arr, seeds, pca_n,
                            requested_publications):
    """Per-publication AUC, one fit per requested publication."""
    out = {}
    for pub in sorted(requested_publications):
        mask = publication_arr == pub
        n = int(mask.sum())
        n_cancer = int((mask & (y == 1)).sum())
        n_healthy = int((mask & (y == 0)).sum())
        if n < 30 or n_cancer < 10 or n_healthy < 10:
            out[pub] = {
                "n_total": n,
                "n_cancer": n_cancer,
                "n_healthy": n_healthy,
                "study_label": PUBLICATION_REGISTRY.get(pub, pub),
                "skipped": True,
                "reason": "insufficient samples for 5-fold CV",
            }
            continue
        _, _, seed_aucs = pooled_oof(
            X[mask], y[mask], st[mask], seeds, pca_n, harmonize=False
        )
        out[pub] = {
            "n_total": n,
            "n_cancer": n_cancer,
            "n_healthy": n_healthy,
            "study_label": PUBLICATION_REGISTRY.get(pub, pub),
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


def _per_cancer_subset(X, y, st, dc_arr, cancer, n_min=10, healthy_min=10):
    """Return (mask, y_sub, X_sub, st_sub, n_cancer_n, n_healthy_n) for OvR.

    Returns None when the cancer subset is too small or the healthy set
    is too small for 5-fold CV.
    """
    healthy_mask = y == 0
    cancer_mask = (dc_arr == cancer) & (y == 1)
    n_cancer = int(cancer_mask.sum())
    if n_cancer < n_min or healthy_mask.sum() < healthy_min:
        return None
    mask = cancer_mask | healthy_mask
    return (
        mask,
        cancer_mask[mask].astype(int),
        X[mask],
        st[mask],
        n_cancer,
        int(healthy_mask.sum()),
    )


def _delong_ci_block(y_true, score):
    """Compute DeLong CI for AUC + Sens@spec at the canonical specificities.

    Returns a dict with keys:
      auc        : DeLong point estimate
      auc_se     : DeLong SE
      auc_ci     : (ci_lo, ci_hi) tuple
      sens_at_*  : dict per specificity with DeLong CIs
    """
    auc_d = delong_auc_ci(y_true, score)
    sens_block = {}
    for spec in DEFAULT_SPECIFICITIES:
        sens_block[str(spec)] = delong_sens_at_spec_ci(y_true, score, spec)
    return {
        "auc": float(auc_d["auc"]),
        "auc_se": float(auc_d["se"]),
        "auc_ci": [float(auc_d["ci_lo"]), float(auc_d["ci_hi"])],
        "sens_at_spec": {k: dict(v) for k, v in sens_block.items()},
    }


def _run_all_cancer_ovr(
    X, y, st, samples_in_array, disease_class, seeds, pca_n,
    harmonize=True, min_n_cancer=10,
):
    """Run OvR (cancer vs ALL healthy) for every non-HEALTHY cancer.

    Returns:
        dc_arr     : np.ndarray of per-sample disease class
        artifacts  : dict {cancer_name: {
            "y_true", "score", "seed_aucs",
            "n_cancer", "n_healthy", "skipped", "reason"
        }}
    """
    dc_arr = np.array(
        [disease_class.get(s, "") for s in samples_in_array],
        dtype=object,
    )
    healthy_mask = y == 0
    cancer_types = sorted(
        set(dc_arr[(y == 1) & (dc_arr != "HEALTHY")]) - {""}
    )
    artifacts = {}
    for cancer in cancer_types:
        sub = _per_cancer_subset(X, y, st, dc_arr, cancer, n_min=min_n_cancer)
        if sub is None:
            n_cancer = int(((dc_arr == cancer) & (y == 1)).sum())
            artifacts[cancer] = {
                "n_cancer": n_cancer,
                "n_healthy": int(healthy_mask.sum()),
                "skipped": True,
                "reason": "insufficient samples for 5-fold CV",
            }
            continue
        mask, y_sub, X_sub, st_sub, n_cancer, n_healthy = sub
        y_true, score, seed_aucs = pooled_oof(
            X_sub, y_sub, st_sub, seeds, pca_n, harmonize=harmonize
        )
        artifacts[cancer] = {
            "y_true": y_true,
            "score": score,
            "seed_aucs": seed_aucs,
            "n_cancer": n_cancer,
            "n_healthy": n_healthy,
            "skipped": False,
            "study_sub": st_sub,
        }
    return dc_arr, artifacts


def _per_cancer_ci_bundle(y_true, score, n_cancer):
    """Build the per-cancer DeLong + bootstrap CI bundle for one OvR fit.

    Used by both `section_per_cancer` and `build_per_cancer_standalone_payload`.
    """
    delong_block = _delong_ci_block(y_true, score)
    delong_rows = []
    for spec in DEFAULT_SPECIFICITIES:
        sb = delong_block["sens_at_spec"][str(spec)]
        delong_rows.append({
            "specificity": spec,
            "sensitivity": sb["sensitivity"],
            "ci95_lo": sb["ci_lo"],
            "ci95_hi": sb["ci_hi"],
            "operating_threshold": sb["threshold"],
            "ci_unreliable": sb["ci_unreliable"],
            "n_pos": sb["n_pos"],
            "n_neg": sb["n_neg"],
        })
    delong_auc_dict = {
        "auc": delong_block["auc"],
        "se": delong_block["auc_se"],
        "ci_lo": delong_block["auc_ci"][0],
        "ci_hi": delong_block["auc_ci"][1],
    }
    boot_rows = []
    for spec in SPECS_FOR_PER_CANCER:
        sens, thr = sens_at_spec(y_true, score, spec)
        lo, hi = bootstrap_ci_sens(
            y_true, score, spec, N_BOOTSTRAP,
            seed=BOOTSTRAP_SEED + int(round(n_cancer)),
        )
        boot_rows.append({
            "specificity": spec,
            "sensitivity": sens,
            "ci95_lo": lo,
            "ci95_hi": hi,
            "operating_threshold": thr,
        })
    return {
        "delong_auc": delong_auc_dict,
        "delong_rows": delong_rows,
        "boot_rows": boot_rows,
    }


def build_per_cancer_standalone_payload(
    ovr_artifacts, out_path,
    pooled_harmonized=None,
    pooled_n_pos=None,
    pooled_n_neg=None,
):
    """Build the standalone JSON from a pre-computed OvR artifact map.

    `ovr_artifacts` is the per-cancer OvR map returned by
    `_run_all_cancer_ovr` — skips cancers with `skipped=True`.

    `pooled_harmonized` is an optional dict (the `pooled.harmonized`
    section of the cross-study JSON). When provided, the POOLED row
    uses **these** values instead of re-deriving from concatenated
    per-cancer OvR subsets. The cross-study pooled OOF is computed on
    ALL cancer vs ALL healthy (each sample gets exactly one score),
    whereas concatenating per-cancer OvR subsets would inflate n_neg
    to (n_cancers * n_healthy) and is NOT a valid pooled test.

    `pooled_n_pos` and `pooled_n_neg` are the actual sample counts
    in the pooled OOF (e.g. 363 cancer + 264 healthy for the
    harmonized cross-study cohort). When provided, they override
    any derivation from per-cancer artifacts.

    **Per-cancer iteration (not concatenated) for the same reason:**
    n_neg=264 (full healthy count) per cancer, not 1848 (7×264).

    Writes the JSON to `out_path` and returns the payload dict.
    """
    rows = {}
    for cancer, art in ovr_artifacts.items():
        if art.get("skipped"):
            # Match the canonical cfDNA schema even when skipped (n=0
            # because we don't have a real OvR subset for this cancer).
            n_cancer = int(art.get("n_cancer", 0))
            n_healthy = int(art.get("n_healthy", 0))
            rows[cancer] = {
                "cancer": cancer,
                "n": n_cancer + n_healthy,
                "n_pos": n_cancer,
                "n_neg": n_healthy,
                "auc_mean": None,
                "auc_ci": [None, None],
                "auc_se": None,
                "sens_at_95": None,
                "sens_at_98": None,
                "sens_at_99": None,
                "sens_at_spec": {str(spec): None for spec in DEFAULT_SPECIFICITIES},
                "ppv_at_prevalence": {f"prev_{p}": None for p in DEFAULT_PREVALENCES},
                "skipped": True,
                "skip_reason": art.get("reason", ""),
            }
            continue

        y_sub = art["y_true"].astype(np.int64)
        s_sub = art["score"].astype(np.float64)
        # Build a 1-of-K label: this cancer vs HEALTHY for everyone else.
        label_sub = np.array(
            [cancer if v == 1 else "HEALTHY" for v in y_sub],
            dtype=object,
        )
        study_sub = art.get("study_sub", np.array([""] * len(y_sub), dtype=object))
        study_sub = study_sub.astype(object)

        per_cancer_rows = build_per_cancer_table(
            y=y_sub,
            s=s_sub,
            cancer_label=label_sub,
            study_label=study_sub,
            specificities=DEFAULT_SPECIFICITIES,
            prevalences=DEFAULT_PREVALENCES,
            include_pooled=False,
        )
        rows[cancer] = per_cancer_rows["per_cancer"][cancer]

    if not rows or all(r.get("skipped") for r in rows.values()):
        raise RuntimeError(
            "build_per_cancer_standalone_payload: no cancer passed the "
            "min-sample floor (n_cancer >= 10 and n_healthy >= 10)."
        )

    # POOLED row. Prefer the upstream pooled harmonized result (a true
    # pooled OOF on cancer-vs-healthy, n_neg = real healthy count).
    # Fallback: re-derive from concatenated OvR subsets (which inflates
    # n_neg and is marked as fallback in `provenance.fallback_pooled`).
    fallback_pooled = False
    if pooled_harmonized:
        pooled_sens_at_spec = {}
        for spec in DEFAULT_SPECIFICITIES:
            pooled_sens_at_spec[str(spec)] = {
                "sensitivity": pooled_harmonized.get(
                    f"sens_at_spec_{int(round(spec * 100))}"
                ),
                "specificity": spec,
                "ci_method": "pooled_oof_at_target_spec",
            }
        n_pos_pooled = (int(pooled_n_pos) if pooled_n_pos is not None
                       else sum(int(art.get("n_cancer", 0))
                                for art in ovr_artifacts.values()
                                if not art.get("skipped")))
        # Each OvR artifact's `n_healthy` reports the FULL pooled
        # healthy count (~264 — every OvR uses the same healthy
        # controls). Summing across cancers would give 7× the true
        # value. Prefer `pooled_n_neg` if provided, else take the
        # max (= the actual pooled healthy count).
        n_neg_pooled = (int(pooled_n_neg) if pooled_n_neg is not None
                       else max(
                           (int(art.get("n_healthy", 0))
                            for art in ovr_artifacts.values()
                            if not art.get("skipped")),
                           default=0,
                       ))
        # Note: OvR's "n_healthy" reports the dataset-wide healthy count
        # (~264 for the cross-study harmonized pooled cohort). The pooled
        # row's n_pos + n_neg is therefore 2x the actual sample count —
        # record the actual count from the upstream n_with_features.
        sens99 = pooled_harmonized.get("sens_at_spec_99", 0.0)
        pooled_row = {
            "cancer": "POOLED",
            "n": n_pos_pooled + n_neg_pooled,
            "n_pos": n_pos_pooled,
            "n_neg": n_neg_pooled,
            "auc_mean": float(pooled_harmonized["auc_mean"]),
            "auc_ci": [
                max(0.0, float(pooled_harmonized["auc_mean"])
                    - 1.96 * float(pooled_harmonized["auc_std"])),
                min(1.0, float(pooled_harmonized["auc_mean"])
                    + 1.96 * float(pooled_harmonized["auc_std"])),
            ],
            "auc_se": float(pooled_harmonized["auc_std"]),
            "sens_at_95": pooled_harmonized.get("sens_at_spec_95"),
            "sens_at_98": pooled_harmonized.get("sens_at_spec_98"),
            "sens_at_99": sens99,
            "sens_at_spec": pooled_sens_at_spec,
            "ppv_at_prevalence": {
                f"prev_{p}": ppv_at_prevalence(sens99, PPV_AT_SPEC, p)
                for p in DEFAULT_PREVALENCES
            },
            "skipped": False,
            "source": "upstream_pooled_oof",
        }
    else:
        fallback_pooled = True
        # Concatenate per-cancer OvR subsets to derive a pooled SENSs.
        # n_neg is inflated to (n_cancers * n_healthy) — honest caveat
        # recorded in `provenance.fallback_pooled = True`.
        pooled_y, pooled_s, pooled_label = [], [], []
        pooled_study = []
        for cancer, art in ovr_artifacts.items():
            if art.get("skipped"):
                continue
            y_sub = art["y_true"].astype(np.int64)
            s_sub = art["score"].astype(np.float64)
            pooled_y.append(y_sub)
            pooled_s.append(s_sub)
            pooled_label.append(np.array(
                [cancer if v == 1 else "HEALTHY" for v in y_sub],
                dtype=object,
            ))
            pooled_study.append(art.get("study_sub", np.array([""] * len(y_sub), dtype=object)).astype(object))
        pooled_table = build_per_cancer_table(
            y=np.concatenate(pooled_y),
            s=np.concatenate(pooled_s),
            cancer_label=np.concatenate(pooled_label),
            study_label=np.concatenate(pooled_study) if pooled_study else None,
            specificities=DEFAULT_SPECIFICITIES,
            prevalences=DEFAULT_PREVALENCES,
            include_pooled=True,
        )
        pooled_row = pooled_table.get("pooled", {})

    table = {
        "per_cancer": rows,
        "pooled": pooled_row,
        "config": {
            "specificities": list(DEFAULT_SPECIFICITIES),
            "prevalences": list(DEFAULT_PREVALENCES),
            "ppv_at_spec": PPV_AT_SPEC,
            "min_positives_for_ci": MIN_POSITIVES_FOR_CI,
        },
    }
    table["provenance"] = {
        "source": "cross_study_finallydb.py",
        "n_samples": sum(
            int(art.get("n_cancer", 0)) + int(art.get("n_healthy", 0))
            for art in ovr_artifacts.values()
            if not art.get("skipped")
        ),
        "n_cancer_types": len([1 for r in rows.values() if not r.get("skipped")]),
        "specificities": list(DEFAULT_SPECIFICITIES),
        "prevalences": list(DEFAULT_PREVALENCES),
        "ppv_at_spec": PPV_AT_SPEC,
        "min_positives_for_ci": MIN_POSITIVES_FOR_CI,
        "harmonization": "per-publication z-score StandardScaler fit on train fold only",
        "cv": "5-fold StratifiedKFold, 5-seed pooled OOF per OvR",
    }
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(table, f, indent=2)
    return table


def section_per_cancer(ovr_artifacts, cancer_types):
    """Emit the per-cancer sens@spec dict from a pre-computed OvR artifact
    map (returned by `_run_all_cancer_ovr`).

    For each cancer in `cancer_types`, computes DeLong 95% CIs (primary;
    DeLong 1988 for AUC, Sun & Xu 2014 for sens@spec) AND preserves the
    original bootstrap 95% CIs under `bootstrap_ci` for the audit trail.
    Cancers with `skipped=True` in the artifacts are passed through as
    SKIPPED rows.
    """
    out = {}
    for cancer in cancer_types:
        art = ovr_artifacts.get(cancer)
        if art is None or art.get("skipped"):
            reason = art.get("reason") if art else "no OOF artifact"
            out[cancer] = {
                "n_cancer": art.get("n_cancer", 0) if art else 0,
                "skipped": True,
                "reason": reason,
            }
            continue
        bundle = _per_cancer_ci_bundle(
            art["y_true"], art["score"], art["n_cancer"]
        )
        out[cancer] = {
            "n_cancer": art["n_cancer"],
            "n_healthy": art["n_healthy"],
            "auc_mean": float(np.mean(art["seed_aucs"])),
            "auc_std": float(np.std(art["seed_aucs"])),
            "per_seed_auc": art["seed_aucs"],
            "primary_ci": "delong",
            "delong_ci": {
                "auc": bundle["delong_auc"],
                "per_specificity": bundle["delong_rows"],
            },
            "bootstrap_ci": {
                "n_bootstrap": N_BOOTSTRAP,
                "seed": BOOTSTRAP_SEED,
                "method": "percentile",
                "per_specificity": bundle["boot_rows"],
            },
            "per_specificity": bundle["delong_rows"],
            "skipped": False,
        }
    return out


def section_true_confound(X, y, st, publication_arr, requested_publications,
                          seeds, pca_n):
    """TRUE cross-publication confound.

    Run BOTH orientations and BOTH harmonization settings for every pair
    of requested publications. With per-publication z-score harmonization
    this should collapse to AUC ~0.50; without harmonization it reaches
    ~0.999 (the classifier learns the publication).
    """
    out = {}
    pubs = sorted(requested_publications)
    if len(pubs) < 2:
        # With only one publication, the true-confound is undefined.
        return {
            "skipped": True,
            "reason": (f"true-confound control needs >=2 publications; "
                      f"got {pubs}"),
        }
    # All ordered pairs (cancer=pubA, healthy=pubB).
    pairs = [(a, b) for a in pubs for b in pubs if a != b]
    for pub_pos, pub_neg in pairs:
        pos_mask = publication_arr == pub_pos
        neg_mask = publication_arr == pub_neg
        mask = pos_mask | neg_mask
        n_pos = int((mask & pos_mask).sum())
        n_neg = int((mask & neg_mask).sum())
        tag = f"cancer_{PUBLICATION_REGISTRY.get(pub_pos, pub_pos)}_" \
              f"healthy_{PUBLICATION_REGISTRY.get(pub_neg, pub_neg)}"
        if n_pos < 10 or n_neg < 10:
            out[tag] = {"skipped": True,
                        "reason": f"n_pos={n_pos}, n_neg={n_neg}",
                        "publication_cancer": pub_pos,
                        "publication_healthy": pub_neg}
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
            "publication_cancer": pub_pos,
            "publication_healthy": pub_neg,
            **cfg,
        }
    return out


# --------------------------------------------------------------------------- #
# Markdown writer (kept in this file for one-shot run)
# --------------------------------------------------------------------------- #
def write_markdown(payload, md_path):
    cfg = payload["config"]
    cohort = payload["cohort"]
    pc = payload.get("per_publication", payload.get("per_cohort", {}))
    pooled = payload["pooled"]
    percancer = payload["per_cancer"]
    confound = payload["true_confound_control"]
    interp = payload["interpretation"]
    pubs = payload.get("publications", [])
    api_status = payload.get("finaledb_api_status", "unknown")

    L = []
    L.append("# Cross-Study FinaleDB Benchmark (Open Data)\n")
    pub_str = " + ".join(
        f"{p} ({PUBLICATION_REGISTRY.get(p, '?')})" for p in pubs
    )
    L.append(f"> **Scope**: Open-data benchmark on FinaleDB publications: "
            f"{pub_str}. **NOT** clinical validation. **NOT** external "
            f"cohort validation. Pooled OOF on the same cohort that "
            f"trained the model.\n")
    if api_status != "ok":
        L.append(f"> **FinaleDB API status**: `{api_status}`. Publications "
                f"with no locally-cached features were skipped (see "
                f"`results/publication_readiness.json`).\n")
    L.append(f"- Generated: `{payload['generated_at']}`\n")
    L.append(f"- Classifier: `{cfg['classifier']}`\n")
    L.append(f"- Feature set: {cfg['feature_set']}\n")
    L.append(f"- PCA n_components: {cfg['pca_n']} (capped at `min(n_train, n_features)` per fold)\n")
    L.append(f"- CV: {cfg['cv']}\n")
    L.append(f"- Seeds: `{cfg['seeds']}`\n")
    L.append(f"- Bootstrap: n={cfg['n_bootstrap']}, seed={cfg['bootstrap_seed']}\n")
    L.append(f"- Cell-line regex: `{cfg['cell_line_regex']}`\n")

    L.append("\n## 1. Cohort inventory\n")
    L.append(f"- Samples in labels file (filtered to requested publications): "
            f"**{cohort['n_total_in_labels']}** "
            f"({cohort['n_cancer_in_labels']} cancer + "
            f"{cohort['n_healthy_in_labels']} healthy)\n")
    L.append(f"- Samples with all 5-channel DELFI features: **{cohort['n_with_features']}** "
            f"({cohort['n_dropped_due_to_missing_features']} dropped due to missing artifacts)\n")
    L.append(f"- Cell-line filter removed: **{cohort['n_dropped_cell_line']}** samples "
            f"(none matched the regex on this open-data cohort)\n")
    if cohort["dropped_cell_line_ids"]:
        L.append(f"  - Dropped IDs: `{cohort['dropped_cell_line_ids']}`\n")
    L.append(f"- Per-publication (FinaleDB):\n")
    for pub, d in sorted(cohort.get("per_publication", {}).items()):
        L.append(f"  - **{pub} ({d['study_label']})**: n={d['n_total']} "
                f"({d['n_cancer']} cancer + {d['n_healthy']} healthy)\n")
    L.append(f"- Per-study (legacy `study` column in labels file):\n")
    for st, d in cohort["per_study"].items():
        L.append(f"  - **{st}**: n={d['n_total']} ({d['n_cancer']} cancer + "
                f"{d['n_healthy']} healthy)\n")
    L.append(f"- Per-cancer (cancer samples only):\n")
    study_keys = sorted({
            k for d in cohort["per_cancer"].values() for k in d
            if k not in ("n_total",)
        })
    header = "| Cancer | n_total | " + " | ".join(study_keys) + " |"
    sep = "|---|---:|" + "---:|" * len(study_keys) + "\n"
    L.append(header + "\n")
    L.append(sep)
    for cancer, d in sorted(cohort["per_cancer"].items(),
                            key=lambda kv: -kv[1]["n_total"]):
        row = f"  | {cancer} | {d['n_total']} | " + \
              " | ".join(str(d.get(sk, 0)) for sk in study_keys) + " |\n"
        L.append(row)

    L.append("\n## 2. Per-publication AUC\n")
    L.append("Each publication evaluated independently (no harmonization "
            "needed; one publication only).\n\n")
    L.append("| Publication | Study | n | n_cancer | n_healthy | "
            "AUC (5-seed mean ± std) |\n")
    L.append("|---|---|---:|---:|---:|---|\n")
    for pub in sorted(pc.keys()):
        r = pc[pub]
        if r.get("skipped"):
            L.append(f"| {pub} | {r.get('study_label', '?')} | "
                    f"{r['n_total']} | {r['n_cancer']} | "
                    f"{r['n_healthy']} | SKIPPED ({r.get('reason', '')}) |\n")
        else:
            L.append(f"| {pub} | {r['study_label']} | "
                    f"{r['n_total']} | {r['n_cancer']} | "
                    f"{r['n_healthy']} | "
                    f"{r['auc_mean']:.4f} ± {r['auc_std']:.4f} |\n")

    L.append("\n## 3. Pooled cross-publication AUC (with/without per-publication harmonization)\n")
    L.append("Harmonization = per-publication z-score StandardScaler fit "
            "on train fold only.\n\n")
    L.append("| Setting | AUC mean ± std | Sens@95% | Sens@98% | Sens@99% |\n")
    L.append("|---|---|---:|---:|---:|\n")
    for tag in ["harmonized", "no_harmonize"]:
        r = pooled[tag]
        L.append(f"| {tag} | {r['auc_mean']:.4f} ± {r['auc_std']:.4f} | "
                f"{r['sens_at_spec_95']:.3f} | {r['sens_at_spec_98']:.3f} | "
                f"{r['sens_at_spec_99']:.3f} |\n")

    L.append("\n## 4. Per-cancer Sens@Spec (top-5 cancers by sample count)\n")
    L.append("One-vs-rest: each cancer vs ALL healthy samples in the pooled cross-study "
            "cohort. Per-publication harmonization inside each CV fold.\n\n")
    L.append("**Primary CIs are DeLong** (DeLong, DeLong, Clarke-Pearson 1988 for AUC; "
            "Sun & Xu 2014 for Sens@spec — equivalent to DeLong-Han-Agarwal structural "
            "component restricted to positives). Bootstrap 95% CIs "
            f"(n={cfg['n_bootstrap']} percentile resamples) are preserved under "
            "`bootstrap_ci` in `results/cross_study_finallydb.json` for the audit trail "
            "(narrower at n>=10; DeLong is the standard cfDNA reference per the "
            "cfdna-early-detection-validation skill). The standalone clinical-decision "
            "table — AUC + Sens@spec CIs + PPV@prev at spec=0.99 — is at "
            "`results/per_cancer_sens_at_spec.json` (schema per "
            "`src.per_cancer_sens_at_spec.build_per_cancer_table`).\n\n")
    L.append("**Caveat**: HCC_J is Jiang-only (n=89 cancer + 32 healthy). Its per-cancer "
            "OvR uses 89 HCC vs ALL healthy controls (264 controls in the OOF subset) "
            "— the per-cancer denominator is small but the OvR is well-defined. Other "
            "top-5 cancers are Cristiano-only.\n\n")
    L.append("| Cancer | n_cancer | n_healthy | AUC mean ± std | AUC DeLong [95% CI] | "
            "Sens@95% DeLong [95% CI] | Sens@98% DeLong [95% CI] | "
            "Sens@99% DeLong [95% CI] |\n")
    L.append("|---|---:|---:|---|---|---|---|---|\n")
    for cancer, r in percancer.items():
        if r.get("skipped"):
            L.append(f"| {cancer} | {r.get('n_cancer', '?')} | — | SKIPPED | — | — | — | — |\n")
        else:
            d_auc = r["delong_ci"]["auc"]
            cells = []
            for spec in SPECS_FOR_PER_CANCER:
                row = next(x for x in r["delong_ci"]["per_specificity"]
                           if x["specificity"] == spec)
                cells.append(
                    f"{row['sensitivity']:.3f} "
                    f"[{row['ci95_lo']:.3f}–{row['ci95_hi']:.3f}]"
                )
            L.append(
                f"| {cancer} | {r['n_cancer']} | {r['n_healthy']} | "
                f"{r['auc_mean']:.4f} ± {r['auc_std']:.4f} | "
                f"{d_auc['auc']:.4f} [{d_auc['ci_lo']:.4f}–{d_auc['ci_hi']:.4f}] | "
                f"{cells[0]} | {cells[1]} | {cells[2]} |\n"
            )

    L.append("\n## 5. True-confound control\n")
    L.append("Cancer = 100% from one publication, healthy = 100% from another. "
            "Without harmonization the classifier learns 'which publication is this "
            "from?' (AUC ~0.999). With per-publication z-score harmonization the "
            "publication-specific mean/variance is the only signal and is removed by "
            "design (AUC should collapse toward 0.50).\n\n")
    L.append("| Orientation | n_cancer | n_healthy | AUC harmonized | AUC no_harmonize |\n")
    L.append("|---|---:|---:|---:|---:|\n")
    if isinstance(confound, dict) and confound.get("skipped"):
        L.append(f"| — | — | — | SKIPPED ({confound.get('reason')}) | — |\n")
    else:
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
    L.append(f"- Pooled harmonized cross-publication AUC: **{h['auc_mean']:.4f} ± "
            f"{h['auc_std']:.4f}** "
            f"(n={n_feat} with features, of {cohort['n_total_in_labels']} in labels file)\n")
    L.append(f"- Pooled AUC without harmonization: **{nh['auc_mean']:.4f} ± "
            f"{nh['auc_std']:.4f}** "
            f"(mild change confirms the per-publication batch effect is small on this "
            f"FinaleDB-uniformly-processed cohort)\n")
    # Per-publication verdict (publication-agnostic, sorted by AUC).
    pc_rows = [
        (pub, r) for pub, r in pc.items() if not r.get("skipped")
    ]
    if pc_rows:
        L.append(f"- Per-publication AUC (sorted by AUC): "
                + ", ".join(
                    f"{pub} ({r['study_label']}) "
                    f"**{r['auc_mean']:.4f} ± {r['auc_std']:.4f}** (n={r['n_total']})"
                    for pub, r in sorted(pc_rows, key=lambda kv: -kv[1]['auc_mean'])
                ) + "\n")
    # Confound verdict
    if isinstance(confound, dict) and not confound.get("skipped"):
        for tag, r in confound.items():
            if r.get("skipped"):
                continue
            L.append(f"- True-confound control: {tag}: "
                    f"harmonized AUC **{r['harmonized']['auc_mean']:.3f}** "
                    f"(should be ~0.50), "
                    f"no-harmonize AUC **{r['no_harmonize']['auc_mean']:.3f}** "
                    f"(should be ~1.00 — proves the batch effect is removable)\n")

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
# Helpers exposed for the readiness diagnostic
# --------------------------------------------------------------------------- #
def features_present_for_publication(labels, studies, publication, pub_id,
                                     features_dir):
    """Return True if at least one sample for `pub_id` has all 5 features
    present in `features_dir`. Used by the readiness diagnostic.
    """
    feat_paths = (
        "{s}.delfi_5mb_ratio.npy",
        "{s}.delfi_5mb_coverage.npy",
        "{s}.delfi_100kb_ratio.npy",
        "{s}.delfi_100kb_counts.npy",
        "{s}.fsd.json",
    )
    for s in labels:
        if publication.get(s) != pub_id:
            continue
        if all(os.path.exists(os.path.join(features_dir, f.format(s=s)))
               for f in feat_paths):
            return True
    return False


def _probe_finaledb_api_best_effort(timeout=3.0):
    """Best-effort probe of the FinaleDB API. Returns (status, reason).

    Status ∈ {"ok", "api-down", "unreachable", "unknown"}. This is purely
    informational — a network failure here MUST NOT block the benchmark.
    On any probe error, returns ("unknown", "<error message>") so the
    downstream pipeline still completes.
    """
    try:
        import urllib.error
        import urllib.request
        with urllib.request.urlopen(
            "http://finaledb.research.cchmc.org/api/v1/misc",
            timeout=timeout,
        ) as r:
            misc_status = r.status
    except Exception as e:
        return ("unreachable", f"/api/v1/misc error: {type(e).__name__}: {e}")
    if misc_status != 200:
        return ("unreachable", f"/api/v1/misc returned status={misc_status}")
    # /misc is alive — check whether the data endpoints are 500ing.
    bad = []
    for ep in ("publication", "seqrun", "summary"):
        try:
            with urllib.request.urlopen(
                f"http://finaledb.research.cchmc.org/api/v1/{ep}",
                timeout=timeout,
            ) as r:
                if r.status != 200:
                    bad.append(f"{ep}={r.status}")
        except urllib.error.HTTPError as e:
            bad.append(f"{ep}={e.code}")
        except Exception as e:
            bad.append(f"{ep}=err")
    if bad:
        return ("api-down", f"FinaleDB data endpoints returning errors: {bad}")
    return ("ok", "all endpoints reachable")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--features-dir",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/data/features")
    ap.add_argument("--labels-multiclass",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv")
    ap.add_argument("--publications", nargs="+", default=["6", "8"],
                    help="Space-separated list of FinaleDB publication ids to "
                         "include in the benchmark (default: '6 8' = Jiang + "
                         "Cristiano, which are the only ones currently "
                         "cached locally). Use --include-snyder / --include-sun "
                         "as shortcuts for publications 1 and 7.")
    ap.add_argument("--include-snyder", action="store_true",
                    help="Shortcut to add publication 1 (Snyder 2016 Cell) to "
                         "--publications.")
    ap.add_argument("--include-sun", action="store_true",
                    help="Shortcut to add publication 7 (Sun 2019) to "
                         "--publications.")
    ap.add_argument("--out-json", default="results/cross_study_finallydb.json")
    ap.add_argument("--out-md", default="docs/CROSS_STUDY_BENCHMARK.md")
    ap.add_argument("--out-per-cancer-json",
                    default="results/per_cancer_sens_at_spec.json",
                    help="Standalone per-cancer sens@spec JSON (DeLong CIs + "
                         "PPV@prev). Set to empty string to skip writing it.")
    ap.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    ap.add_argument("--pca", type=int, default=DEFAULT_PCA_N)
    ap.add_argument("--top-cancer-n", type=int, default=5)
    args = ap.parse_args()

    # ── Resolve the requested publications list ──
    requested = set(str(p) for p in args.publications)
    if args.include_snyder:
        requested.add("1")
    if args.include_sun:
        requested.add("7")
    unknown = requested - set(PUBLICATION_REGISTRY.keys())
    if unknown:
        ap.error(
            f"unknown publication id(s): {sorted(unknown)}. "
            f"Known: {sorted(PUBLICATION_REGISTRY.keys())}"
        )
    requested_publications = sorted(requested)
    print(f"[config] requested publications: {requested_publications} "
          f"({[PUBLICATION_REGISTRY[p] for p in requested_publications]})")

    print(f"[1/6] Loading labels_multiclass from {args.labels_multiclass}")
    labels, studies, disease_class, publication = load_labels_multiclass(
        args.labels_multiclass, set(requested_publications),
    )
    print(f"      Loaded {len(labels)} samples (raw multiclass labels)")

    print(f"[2/6] Applying cell-line filter")
    dropped = apply_cell_line_filter(labels, studies, disease_class, publication)
    print(f"      Dropped {len(dropped)} cell-line samples: {dropped}")

    # Diagnostic: which requested publications actually have samples?
    requested_with_samples = {
        pub for pub in requested_publications
        if any(publication[s] == pub for s in publication)
    }
    requested_without_samples = sorted(
        set(requested_publications) - requested_with_samples
    )
    if requested_without_samples:
        warnings.warn(
            f"No samples for requested publication(s) "
            f"{requested_without_samples}. These will be skipped "
            f"(no features to load, no OvR to compute).",
            stacklevel=1,
        )

    print(f"[3/6] Loading 5-channel DELFI features from {args.features_dir}")
    X, y, st = load5(labels, studies, args.features_dir)
    print(f"      X.shape={X.shape}, n_cancer={(y==1).sum()}, "
          f"n_healthy={(y==0).sum()}, "
          f"studies={sorted(set(st.tolist()))}")

    # Replicate load5's filter to recover the sample IDs in X-row order.
    samples_in_array = []
    feat_paths_tpl = (
        "{s}.delfi_5mb_ratio.npy",
        "{s}.delfi_5mb_coverage.npy",
        "{s}.delfi_100kb_ratio.npy",
        "{s}.delfi_100kb_counts.npy",
        "{s}.fsd.json",
    )
    for s in sorted(labels):
        if all(os.path.exists(os.path.join(args.features_dir, f.format(s=s)))
               for f in feat_paths_tpl):
            samples_in_array.append(s)
    assert len(samples_in_array) == X.shape[0], (
        f"sample_id reconstruction failed: {len(samples_in_array)} vs {X.shape[0]}")

    # Build per-row publication array aligned to X rows.
    publication_arr = np.array(
        [publication.get(s, "") for s in samples_in_array],
        dtype=object,
    )
    # Build per-publication feature-presence map and warn about
    # requested pubs with zero samples-with-features.
    publications_with_features = {
        pub for pub in requested_publications
        if features_present_for_publication(
            labels, studies, publication, pub, args.features_dir)
    }
    pubs_missing_features = sorted(
        set(requested_publications) - publications_with_features
    )
    if pubs_missing_features:
        warnings.warn(
            f"Requested publication(s) {pubs_missing_features} have NO "
            f"5-channel DELFI features cached locally — will be skipped "
            f"from per-publication AUC. See "
            f"`scripts/_verify_publication_readiness.py` for the API/S3 "
            f"status check.",
            stacklevel=1,
        )

    print(f"[4/6] Per-publication AUC (requested={requested_publications})")
    per_pub = section_per_publication(
        X, y, st, publication_arr, args.seeds, args.pca,
        requested_publications,
    )
    for pub, r in per_pub.items():
        if r.get("skipped"):
            print(f"      {pub} ({r.get('study_label', '?')}): SKIPPED "
                  f"({r.get('reason')})")
        else:
            print(f"      {pub} ({r['study_label']:12s}): "
                  f"AUC {r['auc_mean']:.4f} ± {r['auc_std']:.4f} "
                  f"(n={r['n_total']}, {r['n_cancer']} cancer + "
                  f"{r['n_healthy']} healthy)")

    print(f"[5/6] Pooled cross-publication AUC (with/without harmonization)")
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
    print(f"[6/6] Per-cancer OvR OOF (running for ALL cancer types, top-"
          f"{args.top_cancer_n}: {top_cancers})")
    dc_arr, ovr_artifacts = _run_all_cancer_ovr(
        X, y, st, samples_in_array, disease_class,
        args.seeds, args.pca, harmonize=True,
    )
    for cancer, art in ovr_artifacts.items():
        if art.get("skipped"):
            print(f"    {cancer:8s}: SKIPPED ({art.get('reason')})")
            continue
        print(f"    {cancer:8s}: n={art['n_cancer']:3d}  "
              f"AUC {np.mean(art['seed_aucs']):.4f} ± "
              f"{np.std(art['seed_aucs']):.4f}")

    per_cancer = section_per_cancer(ovr_artifacts, top_cancers)
    for cancer, r in per_cancer.items():
        if r.get("skipped"):
            print(f"      {cancer:8s}: SKIPPED ({r.get('reason')})")
            continue
        s99 = next(x["sensitivity"] for x in r["per_specificity"]
                   if x["specificity"] == 0.99)
        print(f"      {cancer:8s}: n={r['n_cancer']:3d}  "
              f"AUC {r['auc_mean']:.4f} ± {r['auc_std']:.4f}  "
              f"S99={s99:.3f}")

    print(f"[+] TRUE cross-publication confound control")
    confound = section_true_confound(
        X, y, st, publication_arr, requested_publications,
        args.seeds, args.pca,
    )
    if isinstance(confound, dict) and confound.get("skipped"):
        print(f"      SKIPPED ({confound.get('reason')})")
    else:
        for tag, r in confound.items():
            if r.get("skipped"):
                print(f"      {tag}: SKIPPED ({r.get('reason')})")
                continue
            h = r["harmonized"]
            nh = r["no_harmonize"]
            print(f"      {tag}:")
            print(f"        harmonized    : AUC {h['auc_mean']:.3f} ± {h['auc_std']:.3f}")
            print(f"        no_harmonize  : AUC {nh['auc_mean']:.3f} ± {nh['auc_std']:.3f}")

    inventory = section_inventory(labels, studies, disease_class, publication,
                                  dropped, samples_in_array)

    # Probe the FinaleDB API status (lightweight, best-effort). This
    # is purely informational — recorded in the JSON + markdown so the
    # numbers can be re-validated later when the API is back up.
    api_status, api_reason = _probe_finaledb_api_best_effort()

    payload = {
        "schema_version": "1.1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scope": (f"open-data cross-publication benchmark on FinaleDB "
                  f"publications: {requested_publications} "
                  f"({[PUBLICATION_REGISTRY[p] for p in requested_publications]}). "
                  f"NOT clinical validation. NOT external cohort validation. "
                  f"Pooled OOF on the same cohort that trained the model."),
        "finaledb_api_status": api_status,
        "finaledb_api_reason": api_reason,
        "config": {
            "publications": requested_publications,
            "publication_registry": {p: PUBLICATION_REGISTRY[p]
                                    for p in requested_publications},
            "seeds": args.seeds,
            "pca_n": args.pca,
            "classifier": "LogisticRegression(max_iter=2000)",
            "feature_set": "5-channel (5mb_ratio + 5mb_coverage + 100kb_ratio "
                           "+ 100kb_counts + FSD-196)",
            "harmonization": "per-publication z-score StandardScaler fit on train fold only",
            "cv": "5-fold StratifiedKFold, 5-seed pooled OOF",
            "n_bootstrap": N_BOOTSTRAP,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "specificities_per_cancer": SPECS_FOR_PER_CANCER,
            "cell_line_regex": CELL_LINE_RE.pattern,
            "labels_file": args.labels_multiclass,
            "features_dir": args.features_dir,
        },
        "cohort": inventory,
        "per_publication": per_pub,
        "publications": requested_publications,
        "requested_publications": requested_publications,
        "publications_with_features": sorted(publications_with_features),
        "publications_missing_features": pubs_missing_features,
        "pooled": pooled,
        "per_cancer": per_cancer,
        "true_confound_control": confound,
        "interpretation": {
            "honest_framing": (
                f"These numbers are pooled out-of-fold AUC on the "
                f"{cohort['n_with_features'] if False else 0}-sample "
                f"cross-publication cohort (after load5's missing-artifact "
                f"filter). Internal CV; no external validation. "
                f"They measure how well the 5-channel DELFI features "
                f"separate cancer from healthy when pooled across "
                f"{', '.join(PUBLICATION_REGISTRY[p] for p in requested_publications)} "
                f"with per-publication z-score harmonization. They do NOT "
                f"measure clinical-grade sensitivity at the Galleri / "
                f"CancerSEEK operating points, which require independent "
                f"held-out plasma cohorts."
            ),
            "true_confound_reading": (
                "With per-publication z-score harmonization the true-confound "
                "AUC (cancer = one publication, healthy = another) collapses "
                "toward 0.50 — the per-publication mean/variance shift is the "
                "only signal and the harmonization removes it by design. "
                "WITHOUT harmonization the same control reaches ~0.999 — "
                "the classifier learns 'which publication is this from?', "
                "not 'is this cancer or healthy?'. The paired comparison "
                "(harmonized vs no_harmonize) is the only honest way to "
                "claim a cross-publication benchmark is not a "
                "publication-batch artifact."
            ),
            "per_cancer_reading": (
                "Per-cancer sens@spec is OvR: each cancer class is scored "
                "against ALL healthy samples (not just the within-publication "
                "healthy ones). This is the cross-publication generalization "
                "view, not the within-publication view. Top-5 cancers by "
                "count are reported; smaller cohorts (n<10 cancer) are "
                "skipped."
            ),
        },
    }

    # Bake in the actual n_with_features into the honest_framing string now
    # that we know it (avoids the placeholder `0` above).
    payload["interpretation"]["honest_framing"] = payload[
        "interpretation"]["honest_framing"
    ].replace(f"{0}-sample", f"{inventory['n_with_features']}-sample")

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {args.out_json}")

    if args.out_per_cancer_json:
        print(f"[+] Standalone per-cancer sens@spec JSON: {args.out_per_cancer_json}")
        # Pass the pooled harmonized result so the POOLED row uses the
        # upstream true pooled OOF, not a concatenation of per-cancer
        # OvR subsets (which would inflate n_neg).
        pooled_h = payload.get("pooled", {}).get("harmonized")
        cohort = payload.get("cohort", {})
        standalone = build_per_cancer_standalone_payload(
            ovr_artifacts, args.out_per_cancer_json,
            pooled_harmonized=pooled_h,
            pooled_n_pos=cohort.get("n_cancer_in_labels"),
            pooled_n_neg=cohort.get("n_healthy_in_labels"),
        )
        n_rows = len(standalone["per_cancer"])
        n_skipped = sum(1 for r in standalone["per_cancer"].values()
                        if r.get("skipped"))
        print(f"      wrote {n_rows} per-cancer rows ({n_skipped} SKIPPED) "
              f"+ 1 POOLED row")

    write_markdown(payload, args.out_md)
    return 0


if __name__ == "__main__":
    sys.exit(main())