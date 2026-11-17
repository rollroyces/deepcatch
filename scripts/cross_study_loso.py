#!/usr/bin/env python3
"""Leave-one-study-out (LOSO) cross-publication benchmark for DeepCatch.

WHY THIS EXISTS
---------------
The pooled cross-study benchmark (``scripts/cross_study_finallydb.py``)
never holds a publication out. Its CV splits the *merged* cohort, so a
fold can train on samples from publication B and test on samples from
publication A. The pooled AUC therefore measures "separate cancer from
healthy when both classes are present in both publications", NOT
"generalise to a publication the model has never seen".

This script measures the thing the pooled number does not: train on one
publication, test on the other.

PROTOCOL (mirrors the pooled pipeline)
--------------------------------------
1. Load the same 5-channel DELFI features (``load5``) and the same
   labels file as the headline benchmark.
2. Split by publication: train = publication A, test = publication B.
3. Harmonization (per-publication z-score), two modes:
   - ``self_zscore`` (primary): each publication is z-scored with its
     own statistics. The held-out publication is normalized by its own
     data (unsupervised, no labels) -- the "new site normalizes itself"
     deployment assumption, and the closest analogue to the pooled
     benchmark's per-publication z-score.
   - ``train_global`` (sensitivity): one scaler fit on the training
     publication is applied to the test publication; no test-set
     statistics at all.
4. PCA (n_components = min(pca_n, n_train)) fit on train only, then
   LogisticRegression(max_iter=2000) fit on train only, then score the
   test publication.
5. Reported per direction: pooled test AUC + DeLong 95% CI, sens@spec
   grid, per-cancer one-vs-rest AUC on the test publication, and a
   PCA-component sweep (the top-k components of one fitted PCA are
   nested, so the sweep is exact rather than a re-fit).

CONTROLS
--------
``within_study_cross_type_control``: leave one cancer type out of the
test publication, train on the remaining types + healthy, test on the
held-out type. This separates "transfers to an unseen STUDY" from
"detects an unseen CANCER TYPE" -- the Jiang -> Cristiano direction
alone cannot, because Jiang is HCC-only while Cristiano contains no HCC.

Also reported: in-domain (within-publication) 5-seed x 5-fold CV
references (seeded, so the artifact reproduces exactly), and a
within-study leave-one-cancer-type-out control.

Implementation note
-------------------
The sibling-pipeline import (``honest_benchmark``) and this repo's
``cross_study_finallydb`` helpers are imported lazily, inside
``_external()``. ``cross_study_finallydb`` imports the sibling checkout
at module import time, so importing it eagerly would make this module
unimportable without that checkout -- and
``test/test_cross_study_loso.py`` would then fail at collection on CI
rather than at run time.

Outputs:
  results/cross_study_loso.json
  docs/CROSS_STUDY_LOSO.md
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path

# The repository's own path layout, resolved from this file: the
# sibling pipeline, the feature cache and the labels are looked up
# next to the repository, then in the home directory, with the
# environment variables overriding both.
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "scripts"))
from _paths import FEATURES_DIR, LABELS_TSV, SIBLING_PIPELINE  # noqa: E402


import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))

# Sister pipeline (cfdna-fragmentomics-pipeline): 5-channel loader +
# per-study harmonizer. Same imports as the headline benchmark --
# but deferred, because that checkout is absent from a fresh clone and
# from CI (see _external below).
from src.per_cancer_sens_at_spec import (  # noqa: E402
    DEFAULT_SPECIFICITIES,
    SCREENING_SPECIFICITIES,
    delong_auc_ci,
)


def _external():
    """Import the sibling pipeline and this repo's cross-study helpers.

    ``scripts/cross_study_finallydb.py`` imports the sibling
    ``cfdna-fragmentomics-pipeline`` checkout at module import time, and
    that checkout is not present in a fresh clone or on a CI runner.
    Importing it eagerly would make this module unimportable there, and
    the unit tests in ``test/test_cross_study_loso.py`` would then fail
    at collection rather than at run time. Both imports are therefore
    deferred to the call sites that actually need the data.

    Returns (load5, helpers) where helpers is the
    ``cross_study_finallydb`` module.
    """
    import importlib

    pipeline = os.environ.get(
        "DEEPCATCH_SISTER_PIPELINE",
        str(SIBLING_PIPELINE),
    )
    for path in (pipeline, os.path.join(pipeline, "scripts")):
        if path not in sys.path:
            sys.path.insert(0, path)
    load5 = importlib.import_module("honest_benchmark").load5
    helpers = importlib.import_module("cross_study_finallydb")
    return load5, helpers

FEATURE_SUFFIXES = (
    "delfi_5mb_ratio.npy",
    "delfi_5mb_coverage.npy",
    "delfi_100kb_ratio.npy",
    "delfi_100kb_counts.npy",
    "fsd.json",
)


# --------------------------------------------------------------------------- #
# Cohort
# --------------------------------------------------------------------------- #
def feature_order(labels: dict, feat_dir: str) -> list[str]:
    """Sample ids in the exact row order ``load5`` produces.

    ``load5`` iterates ``sorted(labels)`` and skips samples missing any
    of the 5 channels, so the row order is reproducible without calling
    it twice.
    """
    out = []
    for s in sorted(labels):
        if all(os.path.exists(os.path.join(feat_dir, f"{s}.{suffix}"))
               for suffix in FEATURE_SUFFIXES):
            out.append(s)
    return out


def load_cohort(features_dir: str, labels_path: str, publications: set[str]):
    load5, helpers = _external()
    labels, studies, disease_class, publication = helpers.load_labels_multiclass(
        labels_path, publications
    )
    dropped = helpers.apply_cell_line_filter(
        labels, studies, disease_class, publication
    )
    X, y, st = load5(labels, studies, feat_dir=features_dir)
    order = feature_order(labels, features_dir)
    if len(order) != X.shape[0]:
        raise RuntimeError(
            f"row alignment failed: order={len(order)} vs X={X.shape[0]}"
        )
    disease = np.array([disease_class.get(s, "UNKNOWN") for s in order])
    return {
        "X": X,
        "y": np.asarray(y, dtype=int),
        "study": np.asarray(st),
        "disease": disease,
        "samples": order,
        "dropped_cell_line": dropped,
    }


# --------------------------------------------------------------------------- #
# Model helpers
# --------------------------------------------------------------------------- #
def harmonize(Xtr: np.ndarray, Xte: np.ndarray, mode: str):
    """Return (Xtr, Xte) under the requested harmonization mode."""
    if mode == "self_zscore":
        # Each publication normalized by its own statistics (unsupervised).
        sc_tr = StandardScaler().fit(Xtr)
        sc_te = StandardScaler().fit(Xte)
        return sc_tr.transform(Xtr), sc_te.transform(Xte)
    if mode == "train_global":
        sc = StandardScaler().fit(Xtr)
        return sc.transform(Xtr), sc.transform(Xte)
    raise ValueError(f"unknown harmonization mode: {mode}")


def fit_pca_lr(Xtr, ytr, pca_n: int, seed: int = 0):
    """Fit PCA + LogisticRegression on train only. Returns (pca, model)."""
    n_components = int(min(pca_n, Xtr.shape[0], Xtr.shape[1]))
    pca = PCA(n_components=n_components, random_state=seed).fit(Xtr)
    model = LogisticRegression(max_iter=2000).fit(pca.transform(Xtr), ytr)
    return pca, model


def score_with(pca, model, X, n_components: int | None = None) -> np.ndarray:
    """Score X with a fitted PCA/LR, optionally truncating to n_components."""
    Z = pca.transform(X)
    if n_components is not None and n_components < Z.shape[1]:
        Z = Z[:, :n_components]
    return model.predict_proba(Z)[:, 1]


def component_sweep(pca, Xtr, ytr, Xte, yte, ks) -> dict:
    """Test AUC as a function of the number of PCA components.

    PCA components are nested, so one PCA fit serves every k; only the
    logistic regression is refit (cheap) per k.
    """
    Ztr = pca.transform(Xtr)
    Zte = pca.transform(Xte)
    out = {}
    for k in ks:
        if k > Ztr.shape[1]:
            continue
        model = LogisticRegression(max_iter=2000).fit(Ztr[:, :k], ytr)
        s = model.predict_proba(Zte[:, :k])[:, 1]
        out[str(k)] = round(float(delong_auc_ci(yte, s)["auc"]), 4)
    return out


def sens_at_spec(y_true, y_score, spec: float) -> float:
    fpr, tpr, _ = roc_curve(y_true, y_score)
    idx = np.where(fpr <= (1.0 - spec))[0]
    if len(idx) == 0:
        return 0.0
    return float(tpr[idx[-1]])


def spec_grid(y_true, y_score) -> dict:
    grid = {}
    for spec in tuple(DEFAULT_SPECIFICITIES) + tuple(SCREENING_SPECIFICITIES):
        grid[f"{spec:g}"] = round(sens_at_spec(y_true, y_score, spec), 4)
    return grid


def auc_with_ci(y_true, y_score) -> dict:
    d = delong_auc_ci(np.asarray(y_true), np.asarray(y_score))
    return {
        "auc": round(float(d["auc"]), 4),
        "ci_lo": round(float(d["ci_lo"]), 4),
        "ci_hi": round(float(d["ci_hi"]), 4),
        "n_pos": int(d["n_pos"]),
        "n_neg": int(d["n_neg"]),
    }


def per_cancer_ovr(y_test, score_test, disease_test) -> dict:
    """One-vs-rest AUC per cancer type on the test publication."""
    out = {}
    healthy_idx = np.where(y_test == 0)[0]
    for cancer in sorted({d for d in disease_test if d != "HEALTHY"}):
        idx = np.where(disease_test == cancer)[0]
        if len(idx) < 5:
            out[cancer] = {"skipped": True, "reason": f"n={len(idx)}"}
            continue
        yy = np.concatenate([np.ones(len(idx)), np.zeros(len(healthy_idx))])
        ss = np.concatenate([score_test[idx], score_test[healthy_idx]])
        out[cancer] = {"n_pos": int(len(idx)), **auc_with_ci(yy, ss)}
    return out


# --------------------------------------------------------------------------- #
# Experiments
# --------------------------------------------------------------------------- #
def loso_direction(cohort, train_study: str, test_study: str, pca_n: int,
                   mode: str, sweep_ks=()) -> dict:
    X, y, st = cohort["X"], cohort["y"], cohort["study"]
    tr = st == train_study
    te = st == test_study
    if tr.sum() == 0 or te.sum() == 0:
        return {"skipped": True,
                "reason": f"missing study(s) {train_study}/{test_study}"}
    Xtr, ytr = X[tr].copy(), y[tr]
    Xte, yte, dte = X[te].copy(), y[te], cohort["disease"][te]
    Xtr, Xte = harmonize(Xtr, Xte, mode)
    pca, model = fit_pca_lr(Xtr, ytr, pca_n)
    score = score_with(pca, model, Xte)

    res = {
        "train_study": train_study,
        "test_study": test_study,
        "harmonization": mode,
        "n_train": int(tr.sum()),
        "n_test": int(te.sum()),
        "train_cancer": int((ytr == 1).sum()),
        "train_healthy": int((ytr == 0).sum()),
        "test_cancer": int((yte == 1).sum()),
        "test_healthy": int((yte == 0).sum()),
        "train_disease_types": sorted({d for d in cohort["disease"][tr]
                                       if d != "HEALTHY"}),
        "pca_n_components": int(pca.components_.shape[0]),
        "pooled_test_auc": auc_with_ci(yte, score),
        "sens_at_spec": spec_grid(yte, score),
        "per_cancer_ovr_on_test": per_cancer_ovr(yte, score, dte),
    }
    if sweep_ks:
        res["pca_component_sweep_auc"] = component_sweep(
            pca, Xtr, ytr, Xte, yte, sweep_ks
        )
    return res


def in_domain_cv(X, y, seeds, pca_n: int):
    """5-seed x 5-fold OOF inside a single publication.

    Mirrors ``scripts/cross_study_finallydb.py:pooled_oof`` -- same
    folds, same per-study z-score harmonization fit on the train fold
    only -- but with a SEEDED PCA. The headline benchmark's PCA was
    unseeded until v2.3.0, which made the 506-sample publication's
    in-domain AUC drift by ~0.003 between runs; both are seeded now, so
    this artifact and the headline benchmark agree to within that
    former drift.

    Returns (auc_mean, auc_std, per_seed_aucs).
    """
    seed_aucs = []
    for i, sd in enumerate(seeds):
        cv = StratifiedKFold(5, shuffle=True, random_state=sd)
        oof = np.zeros(len(y), dtype=float)
        for tr, te in cv.split(X, y):
            Xtr, Xte = harmonize(X[tr].copy(), X[te].copy(), "train_global")
            n_components = int(min(pca_n, Xtr.shape[0], Xtr.shape[1]))
            pca = PCA(n_components=n_components,
                      random_state=1000 + i).fit(Xtr)
            model = LogisticRegression(max_iter=2000).fit(pca.transform(Xtr),
                                                           y[tr])
            oof[te] = model.predict_proba(pca.transform(Xte))[:, 1]
        seed_aucs.append(float(roc_auc_score(y, oof)))
    return (float(np.mean(seed_aucs)), float(np.std(seed_aucs)), seed_aucs)


def in_domain_reference(cohort, study: str, seeds, pca_n: int) -> dict:
    """Within-publication CV reference for one study."""
    st = cohort["study"]
    m = st == study
    auc_mean, auc_std, seed_aucs = in_domain_cv(
        cohort["X"][m], cohort["y"][m], seeds, pca_n
    )
    return {
        "study": study,
        "n": int(m.sum()),
        "n_cancer": int((cohort["y"][m] == 1).sum()),
        "n_healthy": int((cohort["y"][m] == 0).sum()),
        "cv_auc_mean": round(auc_mean, 4),
        "cv_auc_std": round(auc_std, 4),
        "per_seed_auc": [round(float(a), 4) for a in seed_aucs],
    }


def within_study_cross_type_control(cohort, study: str, pca_n: int) -> dict:
    """Leave-one-cancer-type-out inside one publication.

    Train on (healthy + every cancer type except T), test on T vs the
    same healthy controls. Isolates 'unseen cancer type' transfer from
    'unseen study' transfer.
    """
    st, y, dis = cohort["study"], cohort["y"], cohort["disease"]
    m = st == study
    Xs, ys, ds = cohort["X"][m], y[m], dis[m]
    healthy = np.where(ys == 0)[0]
    if len(healthy) < 10:
        return {"skipped": True, "reason": "too few healthy controls"}
    types = sorted({d for d in ds if d != "HEALTHY"})
    out = {}
    for T in types:
        t_idx = np.where(ds == T)[0]
        if len(t_idx) < 5:
            out[T] = {"skipped": True, "reason": f"n={len(t_idx)}"}
            continue
        train_idx = np.concatenate([healthy, np.where(
            (ds != "HEALTHY") & (ds != T))[0]])
        test_pos = t_idx
        if len(np.unique(ys[train_idx])) < 2:
            out[T] = {"skipped": True,
                      "reason": "training set is single-class (this study "
                                "has fewer than two cancer types)"}
            continue
        Xtr, ytr = Xs[train_idx].copy(), ys[train_idx]
        Xte = np.concatenate([Xs[test_pos], Xs[healthy]]).copy()
        yte = np.concatenate([np.ones(len(test_pos)),
                              np.zeros(len(healthy))])
        sc = StandardScaler().fit(Xtr)
        Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
        pca, model = fit_pca_lr(Xtr, ytr, pca_n)
        score = score_with(pca, model, Xte)
        out[T] = {"n_pos": int(len(test_pos)),
                  "n_train": int(len(train_idx)), **auc_with_ci(yte, score)}
    return out


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_markdown(payload: dict, md_path: str) -> None:
    L = []
    L.append("# Leave-One-Study-Out (LOSO) Cross-Publication Benchmark\n")
    L.append("> **Scope**: open-data cross-publication transfer on FinaleDB "
             "publications 6 (jiang) + 8 (cristiano). **NOT** clinical "
             "validation. **NOT** external cohort validation. The test "
             "publication is fully held out of training.\n")
    L.append(f"- Generated: `{payload['generated_at']}`\n")
    L.append(f"- Features: `{payload['config']['feature_set']}`\n")
    L.append(f"- PCA n_components: `{payload['config']['pca_n']}` "
             f"(capped at min(pca_n, n_train) per direction)\n")
    L.append(f"- Classifier: `{payload['config']['classifier']}`\n")

    L.append("\n## 1. Headline: train on one publication, test on the other\n")
    L.append("| Direction | Harmonization | n_train (ca/hea) | n_test (ca/hea) "
             "| Test AUC [DeLong 95% CI] | sens@95% | sens@99% |")
    L.append("|---|---|---|---|---|---|---|\n")
    for d in payload["loso"]:
        if d.get("skipped"):
            continue
        p = d["pooled_test_auc"]
        L.append(f"| {d['train_study']} -> {d['test_study']} | "
                 f"{d['harmonization']} | {d['n_train']} "
                 f"({d['train_cancer']}/{d['train_healthy']}) | "
                 f"{d['n_test']} ({d['test_cancer']}/{d['test_healthy']}) | "
                 f"**{p['auc']}** [{p['ci_lo']}-{p['ci_hi']}] | "
                 f"{d['sens_at_spec'].get('0.95')} | "
                 f"{d['sens_at_spec'].get('0.99')} |")

    L.append("\n## 2. In-domain references (within-publication CV)\n")
    L.append("| Publication | n (ca/hea) | 5-seed x 5-fold OOF AUC |")
    L.append("|---|---|---|\n")
    for s in payload["in_domain"]:
        L.append(f"| {s['study']} | {s['n']} ({s['n_cancer']}/{s['n_healthy']}) "
                 f"| {s['cv_auc_mean']} +/- {s['cv_auc_std']} |")

    L.append("\n## 3. Per-cancer one-vs-rest AUC on the held-out publication\n")
    for d in payload["loso"]:
        if d.get("skipped") or d["harmonization"] != "self_zscore":
            continue
        L.append(f"\n**{d['train_study']} -> {d['test_study']}** "
                 f"(train types: {', '.join(d['train_disease_types'])})\n")
        L.append("| Cancer | n | AUC [DeLong 95% CI] |")
        L.append("|---|---:|---|\n")
        for cancer, row in d["per_cancer_ovr_on_test"].items():
            if row.get("skipped"):
                L.append(f"| {cancer} | - | skipped ({row['reason']}) |")
            else:
                L.append(f"| {cancer} | {row['n_pos']} | "
                         f"{row['auc']} [{row['ci_lo']}-{row['ci_hi']}] |")

    L.append("\n## 4. PCA-component sweep (primary direction)\n")
    L.append("| n_components | Test AUC |")
    L.append("|---:|---:|\n")
    for d in payload["loso"]:
        sweep = d.get("pca_component_sweep_auc")
        if not sweep:
            continue
        L.append(f"\n**{d['train_study']} -> {d['test_study']}** "
                 f"({d['harmonization']})\n")
        L.append("| n_components | Test AUC |")
        L.append("|---:|---:|\n")
        for k, v in sweep.items():
            L.append(f"| {k} | {v} |")

    ctrl = payload.get("within_study_cross_type_control", {})
    if ctrl:
        L.append("\n## 5. Within-study cross-cancer-type control\n")
        L.append("Leave one cancer type out of a single publication; train on "
                 "the remaining types + healthy; test on the held-out type. "
                 "This separates 'unseen study' transfer from 'unseen cancer "
                 "type' transfer.\n")
        L.append("| Study | Cancer | n | AUC [DeLong 95% CI] |")
        L.append("|---|---|---:|---|\n")
        for study, rows in ctrl.items():
            if not isinstance(rows, dict):
                continue
            for cancer, row in rows.items():
                if not isinstance(row, dict) or row.get("skipped"):
                    continue
                L.append(f"| {study} | {cancer} | {row['n_pos']} | "
                         f"{row['auc']} [{row['ci_lo']}-{row['ci_hi']}] |")

    L.append("\n## Interpretation\n")
    for line in payload["interpretation"]:
        L.append(f"- {line}\n")

    repro = payload.get("reproducibility", {})
    if repro:
        L.append("\n## Reproducibility\n")
        for key, text in repro.items():
            L.append(f"- **{key}**: {text}\n")

    # Every element is one line; joining with "\n" (after stripping any
    # trailing newline) keeps the markdown tables well-formed regardless
    # of whether an individual append already carried a newline.
    with open(md_path, "w") as f:
        f.write("\n".join(line.rstrip("\n") for line in L))
        f.write("\n")


def main():
    _, helpers = _external()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features-dir",
                    default=str(FEATURES_DIR))
    ap.add_argument("--labels-multiclass",
                    default=str(LABELS_TSV))
    ap.add_argument("--publications", nargs="+", default=["6", "8"])
    ap.add_argument("--pca", type=int, default=helpers.DEFAULT_PCA_N)
    ap.add_argument("--seeds", type=int, nargs="+",
                    default=list(helpers.DEFAULT_SEEDS))
    ap.add_argument("--out-json", default="results/cross_study_loso.json")
    ap.add_argument("--out-md", default="docs/CROSS_STUDY_LOSO.md")
    args = ap.parse_args()

    pubs = set(args.publications)
    cohort = load_cohort(args.features_dir, args.labels_multiclass, pubs)
    studies = sorted({s for s in cohort["study"]})

    loso = []
    for tr_s in studies:
        for te_s in studies:
            if tr_s == te_s:
                continue
            for mode in ("self_zscore", "train_global"):
                sweep = (5, 10, 20, 50, 100, 200) if (
                    mode == "self_zscore" and tr_s == "jiang" and te_s == "cristiano"
                ) else ()
                loso.append(loso_direction(
                    cohort, tr_s, te_s, args.pca, mode, sweep
                ))

    in_domain = [in_domain_reference(cohort, s, args.seeds, args.pca)
                 for s in studies]
    cross_type = {
        s: within_study_cross_type_control(cohort, s, args.pca)
        for s in studies
    }

    primary = next(
        (d for d in loso if d["train_study"] == "jiang"
         and d["test_study"] == "cristiano"
         and d["harmonization"] == "self_zscore"), None)
    reverse = next(
        (d for d in loso if d["train_study"] == "cristiano"
         and d["test_study"] == "jiang"
         and d["harmonization"] == "self_zscore"), None)

    interpretation = []
    if primary:
        p = primary["pooled_test_auc"]
        interpretation.append(
            f"Train-Jiang / test-Cristiano pooled test AUC = {p['auc']} "
            f"(DeLong 95% CI {p['ci_lo']}-{p['ci_hi']}), n_test={primary['n_test']}."
        )
        interpretation.append(
            "Jiang is HCC-only (train types: "
            f"{', '.join(primary['train_disease_types'])}) while Cristiano "
            "contains no HCC, so this direction conflates two transfers: "
            "unseen STUDY and unseen CANCER TYPE. Read it together with the "
            "within-study cross-type control before attributing the result to "
            "either."
        )
        interpretation.append(
            "The pooled cross-study benchmark in "
            "scripts/cross_study_finallydb.py never holds a publication out, "
            "so it cannot answer this question; these numbers are the "
            "transfer measurement the pooled number is often read as."
        )
    if reverse:
        p = reverse["pooled_test_auc"]
        interpretation.append(
            f"Reverse direction (train-Cristiano / test-Jiang) pooled test AUC "
            f"= {p['auc']} (DeLong 95% CI {p['ci_lo']}-{p['ci_hi']})."
        )
    interpretation.append(
        "self_zscore lets the held-out publication normalize itself using its "
        "own feature statistics (unsupervised, no labels). train_global "
        "applies the training publication's scaler to the test publication "
        "with no test-set statistics; it is reported as a sensitivity check "
        "because the two publications were pre-processed differently."
    )

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "config": {
            "feature_set": "5-channel (5mb_ratio + 5mb_coverage + 100kb_ratio "
                           "+ 100kb_counts + FSD-196)",
            "pca_n": args.pca,
            "classifier": "LogisticRegression(max_iter=2000)",
            "seeds": args.seeds,
            "features_dir": args.features_dir,
            "labels_multiclass": args.labels_multiclass,
            "publications": sorted(pubs),
        },
        "reproducibility": {
            "transfer_numbers": "deterministic: the LOSO directions use a "
                                "seeded PCA (random_state=0) and a fixed "
                                "train/test split by publication, so the "
                                "transfer AUCs reproduce exactly",
            "in_domain_reference": "computed by scripts/cross_study_loso.py:"
                                   "in_domain_cv, which mirrors "
                                   "cross_study_finallydb.py:pooled_oof (same "
                                   "folds, same per-study z-score "
                                   "harmonization fit on the train fold) but "
                                   "with a SEEDED PCA. The headline "
                                   "benchmark's PCA was unseeded before v2.3.0, so its own "
                                   "per-publication AUC drifts ~0.003 "
                                   "run-to-run on the 506-sample publication; "
                                   "these references can therefore differ "
                                   "from results/cross_study_finallydb.json's "
                                   "per_cohort values by that much",
        },
        "cohort": {
            "n_samples_with_features": int(cohort["X"].shape[0]),
            "n_features": int(cohort["X"].shape[1]),
            "n_cancer": int((cohort["y"] == 1).sum()),
            "n_healthy": int((cohort["y"] == 0).sum()),
            "per_study": {
                s: {
                    "n": int((cohort["study"] == s).sum()),
                    "n_cancer": int(((cohort["study"] == s)
                                     & (cohort["y"] == 1)).sum()),
                    "n_healthy": int(((cohort["study"] == s)
                                      & (cohort["y"] == 0)).sum()),
                    "disease_types": sorted({d for d in cohort["disease"][
                        cohort["study"] == s] if d != "HEALTHY"}),
                }
                for s in studies
            },
        },
        "loso": loso,
        "in_domain": in_domain,
        "within_study_cross_type_control": cross_type,
        "interpretation": interpretation,
    }

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"[LOSO] wrote {args.out_json}")

    md_dir = os.path.dirname(args.out_md)
    if md_dir:
        os.makedirs(md_dir, exist_ok=True)
    write_markdown(payload, args.out_md)
    print(f"[LOSO] wrote {args.out_md}")

    for d in loso:
        if d.get("skipped"):
            continue
        p = d["pooled_test_auc"]
        print(f"  {d['train_study']:>10} -> {d['test_study']:<10} "
              f"[{d['harmonization']:<13}] AUC={p['auc']:.4f} "
              f"[{p['ci_lo']}, {p['ci_hi']}]")
    return payload


if __name__ == "__main__":
    main()
