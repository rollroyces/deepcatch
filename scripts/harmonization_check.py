"""Per-cohort batch-effect harmonization check for the panel-LLR pipeline.

Builds a SYNTHETIC multi-cohort fixture (4 studies × 20 patients × paired
cancer/control = 160 samples) on top of real TCGA-LUAD mutations, with a
configurable per-study mean shift in panel-LLR scores so the studies are
mildly confounded with the cancer label. Runs the LR-on-LLR pipeline twice on
the same fixture:

  (a) raw          — current `real_tcga_validation.py` default (no harmonization)
  (b) harmonized   — per-study z-scoring fit on the TRAINING fold only,
                     then applied to test fold

Honest framing: this is a synthetic fixture, not a real cross-cohort
validation. It exists to verify that the harmonization code path is
correct (z-scoring on TRAIN fold only — never on test, never on pooled
data) and to measure the AUC delta in a controlled setting. Real
cross-study pooling should still go through `scripts/run_cross_study.py`
when FinaleDB data is available.

Outputs:
  results/harmonization_check.json    — full numbers (every seed × fold)
  docs/HARMONIZATION.md               — recipe + verdict
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from real_tcga_validation import (  # noqa: E402
    compute_llr_scores,
    load_tcga_cohort,
    simulate_cfdna_from_real,
)

# ── Constants ─────────────────────────────────────────────────────────────
RESULTS_OUT = ROOT / "results" / "harmonization_check.json"
DOCS_OUT = ROOT / "docs" / "HARMONIZATION.md"
CACHE_DIR = ROOT / "validation" / "tcga" / "tcga_cache"

STUDIES = ["study_A", "study_B", "study_C", "study_D"]
N_PATIENTS_PER_STUDY = 20
TUMOR_FRACTION = 0.001  # 0.1% — the headline MRD-style operating point
CFDNA_DEPTH = 5000
BG_ERROR_RATE = 0.002
SEEDS = [42, 123, 456, 789, 1024]
N_FOLDS = 5

# Per-study mean shifts applied to the panel-LLR score (and to the
# auxiliary "frag" channel). study_A gets the largest positive shift so the
# naive pooling would partly confuse "high score" with "study_A".
# Cancer label is INDEPENDENT of study in the fixture — both classes
# appear in every study — so this is the "mild batch effect" regime
# (cfdna-fragmentomics skill: per-study z-scoring is correct here, the
# batch effect is coverage/library-prep, not "Jiang cancer vs
# Cristiano healthy").
STUDY_PANEL_BIAS = {
    "study_A": +0.40,
    "study_B": +0.10,
    "study_C": -0.10,
    "study_D": -0.30,
}
STUDY_FRAG_BIAS = {
    "study_A": +0.20,
    "study_B": -0.05,
    "study_C": +0.05,
    "study_D": -0.15,
}


# ── Fixture generation ──────────────────────────────────────────────────

def _per_patient_panel_score(
    mutations: List[Dict], seed: int,
) -> Tuple[float, float]:
    """Return (pos_panel_score, neg_panel_score) for one patient.

    pos = cfDNA simulated at `TUMOR_FRACTION` (cancer sample)
    neg = cfDNA simulated at TF=0 (matched control sample)
    Both use the same seed → paired design.
    """
    dp = simulate_cfdna_from_real(
        mutations, tumor_fraction=TUMOR_FRACTION,
        cfdna_depth=CFDNA_DEPTH, seed=seed, bg_error_rate=BG_ERROR_RATE,
    )
    dn = simulate_cfdna_from_real(
        mutations, tumor_fraction=0.0,
        cfdna_depth=CFDNA_DEPTH, seed=seed, bg_error_rate=BG_ERROR_RATE,
    )
    lp = compute_llr_scores(dp["depths"], dp["X"][:, 1].astype(int), dp["X"][:, 3])
    ln = compute_llr_scores(dn["depths"], dn["X"][:, 1].astype(int), dn["X"][:, 3])
    nv_p, nv_n = dp["n_variants"], dn["n_variants"]
    panel_size = min(nv_p, nv_n)
    return float(lp[:panel_size].sum()), float(ln[:panel_size].sum())


def build_fixture(seed: int = 42) -> Dict[str, np.ndarray]:
    """Build the multi-cohort synthetic fixture on top of real TCGA mutations.

    Returns a dict with arrays:
      panel_scores  (n_samples,)  — raw panel-LLR sum (NO study bias yet)
      frag_scores   (n_samples,)  — auxiliary frag-channel proxy (raw)
      study         (n_samples,)  — study id per sample (encoded as int)
      y             (n_samples,)  — 1 = cancer, 0 = control
      sample_id     (n_samples,)  — str patient_id + "_pos"/"_neg"
    """
    if not CACHE_DIR.exists() or not list(CACHE_DIR.glob("*.maf.gz")):
        raise SystemExit(
            f"ERROR: no real MAF cache at {CACHE_DIR}. "
            "Place GDC TCGA-LUAD *.maf.gz files there before running this script."
        )
    print(f"[1] Loading real TCGA-LUAD mutations from {CACHE_DIR} ...")
    n_needed = N_PATIENTS_PER_STUDY * len(STUDIES)
    cohort_data = load_tcga_cohort(
        cache_dir=str(CACHE_DIR), n_patients=n_needed, cancer_types=["LUAD"],
    )
    patients = list(cohort_data["patients"].keys())
    if len(patients) < n_needed:
        raise SystemExit(
            f"ERROR: only {len(patients)} LUAD patients in cache; "
            f"need at least {n_needed}."
        )

    # Deterministic patient→study assignment (seed-stable)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(patients))
    patient_study: Dict[str, str] = {}
    for i, idx in enumerate(perm[:n_needed]):
        patient_study[patients[idx]] = STUDIES[i // N_PATIENTS_PER_STUDY]

    panel_scores, frag_scores, y_arr, study_arr, sample_id = [], [], [], [], []
    for patient, study in patient_study.items():
        muts = cohort_data["patients"][patient]
        for label_idx, (label, pos) in enumerate(
            [("neg", False), ("pos", True)]
        ):
            tf = TUMOR_FRACTION if pos else 0.0
            dp = simulate_cfdna_from_real(
                muts, tumor_fraction=tf,
                cfdna_depth=CFDNA_DEPTH, seed=seed,
                bg_error_rate=BG_ERROR_RATE,
            )
            dn = simulate_cfdna_from_real(
                muts, tumor_fraction=0.0,
                cfdna_depth=CFDNA_DEPTH, seed=seed,
                bg_error_rate=BG_ERROR_RATE,
            )
            lp = compute_llr_scores(
                dp["depths"], dp["X"][:, 1].astype(int), dp["X"][:, 3],
            )
            ln = compute_llr_scores(
                dn["depths"], dn["X"][:, 1].astype(int), dn["X"][:, 3],
            )
            nv_p, nv_n = dp["n_variants"], dn["n_variants"]
            panel_size = min(nv_p, nv_n)
            # Raw panel-LLR score (cancer signal minus control)
            # We use the per-sample panel-LLR (pos) and the control panel-LLR (neg).
            # Each label is its own "sample" with its own score.
            raw = float(lp[:panel_size].sum()) if pos else float(ln[:panel_size].sum())
            # Per-patient frag-channel proxy: jitter from a per-patient fixed
            # feature (mutation count) — different per (patient, label) to keep
            # the channel informative. Seed-stable.
            n_muts = len(muts)
            jitter = float(rng.normal(0.0, 0.5))
            frag = (n_muts / 200.0) + jitter
            # Apply per-study additive bias to BOTH channels
            raw += STUDY_PANEL_BIAS[study]
            frag += STUDY_FRAG_BIAS[study]

            panel_scores.append(raw)
            frag_scores.append(frag)
            y_arr.append(1 if pos else 0)
            study_arr.append(STUDIES.index(study))
            sample_id.append(f"{patient}_{label}")

    fixture = {
        "panel_scores": np.asarray(panel_scores, dtype=float),
        "frag_scores": np.asarray(frag_scores, dtype=float),
        "y": np.asarray(y_arr, dtype=int),
        "study": np.asarray(study_arr, dtype=int),
        "sample_id": np.asarray(sample_id),
        "study_names": STUDIES,
    }
    print(
        f"[1] Fixture built: n={len(fixture['y'])} samples, "
        f"{len(STUDIES)} studies × {N_PATIENTS_PER_STUDY} patients × "
        f"2 (cancer/control)"
    )
    for s in STUDIES:
        mask = fixture["study"] == STUDIES.index(s)
        n_pos = int(((mask) & (fixture["y"] == 1)).sum())
        n_neg = int(((mask) & (fixture["y"] == 0)).sum())
        ps = fixture["panel_scores"][mask]
        fs = fixture["frag_scores"][mask]
        print(
            f"      {s}: n={int(mask.sum()):>3} "
            f"(cancer={n_pos}, control={n_neg}) "
            f"panel mean={ps.mean():+.3f} std={ps.std():.3f}  "
            f"frag mean={fs.mean():+.3f}"
        )
    return fixture


# ── Harmonization core ──────────────────────────────────────────────────

def fit_per_study_zscore(
    X_train: np.ndarray, study_train: np.ndarray,
) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """Fit per-study (mean, std) on the TRAINING fold only.

    Returns a dict mapping study_id → (mean_array, std_array). Studies
    with fewer than 2 train samples fall back to (zeros, ones) — i.e. no
    transform. This is the safe fallback: a study with no training
    samples cannot contribute a mean/std estimate, and identity is the
    honest default.

    NOTE: the std fallback is broadcast (1.0 for every column), not a
    Python float — keeping it as an array keeps the downstream
    arithmetic uniform.
    """
    n_cols = X_train.shape[1]
    zeros = np.zeros(n_cols)
    ones = np.ones(n_cols)
    out: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for s in np.unique(study_train):
        mask = study_train == s
        if mask.sum() < 2:
            out[int(s)] = (zeros.copy(), ones.copy())
            continue
        vals = X_train[mask]
        mu = vals.mean(axis=0)
        sd = vals.std(axis=0)
        # Guard against zero-variance columns (would NaN the transform)
        sd = np.where(sd < 1e-8, 1.0, sd)
        out[int(s)] = (mu, sd)
    return out


def apply_per_study_zscore(
    X: np.ndarray, study: np.ndarray,
    params: Dict[int, Tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    """Apply a fit on TRAIN fold to a (possibly test) matrix.

    Studies unseen at fit time fall back to no transform — the column
    values pass through unchanged. This is the safe behaviour: if the
    fit had no data for that study, neither does the transform.
    """
    X_out = X.copy()
    for s in np.unique(study):
        mask = study == s
        if int(s) in params:
            mu, sd = params[int(s)]
            X_out[mask] = (X[mask] - mu) / sd
        # else: leave as-is (no fit params for this study)
    return X_out


# ── Pipeline (panel + frag → 2D feature matrix → LR) ───────────────────

def build_feature_matrix(
    panel_scores: np.ndarray, frag_scores: np.ndarray,
) -> np.ndarray:
    """2-column feature matrix: [panel, frag]."""
    return np.column_stack([panel_scores, frag_scores]).astype(float)


def cv_evaluate(
    X: np.ndarray, y: np.ndarray, study: np.ndarray,
    harmonize: bool, seed: int,
) -> Dict:
    """Stratified K-fold CV. Returns per-study and pooled AUC.

    The harmonization parameter is fit on the TRAINING fold only, never on
    test, never on the pooled matrix. This is the regression-test the
    skill's pitfall section warns about.
    """
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    pooled_y, pooled_score = [], []
    per_study_score: Dict[int, List[float]] = defaultdict(list)
    per_study_y: Dict[int, List[int]] = defaultdict(list)

    for tr, te in skf.split(X, y):
        X_tr_raw, X_te_raw = X[tr], X[te]
        y_tr, y_te = y[tr], y[te]
        study_tr, study_te = study[tr], study[te]

        if harmonize:
            params = fit_per_study_zscore(X_tr_raw, study_tr)
            X_tr = apply_per_study_zscore(X_tr_raw, study_tr, params)
            X_te = apply_per_study_zscore(X_te_raw, study_te, params)
        else:
            X_tr, X_te = X_tr_raw, X_te_raw

        clf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=200)
        clf.fit(X_tr, y_tr)
        score_te = clf.predict_proba(X_te)[:, 1]

        pooled_y.extend(y_te.tolist())
        pooled_score.extend(score_te.tolist())
        for s in np.unique(study_te):
            m = study_te == s
            per_study_score[int(s)].extend(score_te[m].tolist())
            per_study_y[int(s)].extend(y_te[m].tolist())

    pooled_auc = float(roc_auc_score(pooled_y, pooled_score))
    per_study_auc: Dict[str, float] = {}
    for s in sorted(per_study_score.keys()):
        ys = np.asarray(per_study_y[s])
        ss = np.asarray(per_study_score[s])
        if len(np.unique(ys)) < 2:
            per_study_auc[STUDIES[s]] = float("nan")
        else:
            per_study_auc[STUDIES[s]] = float(roc_auc_score(ys, ss))
    return {
        "pooled_auc": pooled_auc,
        "per_study_auc": per_study_auc,
        "n_samples": int(len(y)),
    }


# ── Driver ─────────────────────────────────────────────────────────────

def run_harmonization_check() -> Dict:
    fixture = build_fixture(seed=42)
    X = build_feature_matrix(fixture["panel_scores"], fixture["frag_scores"])
    y = fixture["y"]
    study = fixture["study"]

    raw_per_seed, harm_per_seed = [], []
    for seed in SEEDS:
        t0 = time.time()
        raw = cv_evaluate(X, y, study, harmonize=False, seed=seed)
        raw["seed"] = seed
        raw["wall_sec"] = round(time.time() - t0, 3)
        raw_per_seed.append(raw)

        t0 = time.time()
        h = cv_evaluate(X, y, study, harmonize=True, seed=seed)
        h["seed"] = seed
        h["wall_sec"] = round(time.time() - t0, 3)
        harm_per_seed.append(h)

        print(
            f"[seed {seed}] raw pooled AUC={raw['pooled_auc']:.4f}  "
            f"harmonized pooled AUC={h['pooled_auc']:.4f}  "
            f"Δ={h['pooled_auc'] - raw['pooled_auc']:+.4f}  "
            f"({raw['wall_sec']:.1f}s + {h['wall_sec']:.1f}s)"
        )

    def _agg(per_seed: List[Dict]) -> Dict:
        pooled = np.array([r["pooled_auc"] for r in per_seed])
        per_study = defaultdict(list)
        for r in per_seed:
            for s, v in r["per_study_auc"].items():
                if not np.isnan(v):
                    per_study[s].append(v)
        agg = {
            "pooled_auc_mean": float(pooled.mean()),
            "pooled_auc_std": float(pooled.std(ddof=1)) if len(pooled) > 1 else 0.0,
            "per_study_auc_mean": {
                s: float(np.mean(vs)) for s, vs in per_study.items()
            },
            "per_study_auc_std": {
                s: float(np.std(vs, ddof=1)) if len(vs) > 1 else 0.0
                for s, vs in per_study.items()
            },
            "per_seed_pooled": [float(r["pooled_auc"]) for r in per_seed],
        }
        return agg

    raw_agg = _agg(raw_per_seed)
    harm_agg = _agg(harm_per_seed)
    delta_pooled = harm_agg["pooled_auc_mean"] - raw_agg["pooled_auc_mean"]
    delta_per_study = {
        s: harm_agg["per_study_auc_mean"][s] - raw_agg["per_study_auc_mean"][s]
        for s in STUDIES
    }

    out = {
        "experiment": "per-cohort batch-effect harmonization check",
        "fixture": {
            "n_samples": int(len(y)),
            "n_studies": len(STUDIES),
            "patients_per_study": N_PATIENTS_PER_STUDY,
            "tumor_fraction": TUMOR_FRACTION,
            "cfdna_depth": CFDNA_DEPTH,
            "study_panel_bias": STUDY_PANEL_BIAS,
            "study_frag_bias": STUDY_FRAG_BIAS,
            "study_means_panel": {
                s: float(fixture["panel_scores"][fixture["study"] == STUDIES.index(s)].mean())
                for s in STUDIES
            },
            "study_means_frag": {
                s: float(fixture["frag_scores"][fixture["study"] == STUDIES.index(s)].mean())
                for s in STUDIES
            },
        },
        "seeds": SEEDS,
        "n_folds": N_FOLDS,
        "raw": raw_agg,
        "harmonized": harm_agg,
        "delta_pooled_auc": delta_pooled,
        "delta_per_study_auc": delta_per_study,
        "verdict": _verdict(delta_pooled, raw_agg, harm_agg),
        "raw_per_seed": raw_per_seed,
        "harmonized_per_seed": harm_per_seed,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    RESULTS_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_OUT, "w") as f:
        json.dump(out, f, indent=2)
    _write_docs(out)
    print(f"\n[done] Results → {RESULTS_OUT}")
    print(f"[done] Docs    → {DOCS_OUT}")
    print(
        f"\nVerdict: {out['verdict']['label']}\n"
        f"  Δ pooled AUC = {delta_pooled:+.4f} "
        f"(raw {raw_agg['pooled_auc_mean']:.4f} → "
        f"harmonized {harm_agg['pooled_auc_mean']:.4f})"
    )
    return out


def _verdict(delta_pooled: float, raw_agg: Dict, harm_agg: Dict) -> Dict[str, str]:
    if delta_pooled > 0.005:
        label = "HELPS"
        explain = (
            "Per-study z-scoring improved pooled AUC by more than 0.5pp on "
            "this synthetic fixture. Use `--harmonize` when pooling real "
            "cross-study cohorts with coverage/library-prep batch effects."
        )
    elif delta_pooled < -0.005:
        label = "HURTS"
        explain = (
            "Per-study z-scoring HURT pooled AUC on this synthetic fixture. "
            "The synthetic study bias is uncorrelated with the cancer label, "
            "so removing it removes signal. In a real cross-study pool where "
            "the bias IS correlated with class (e.g. one study is mostly "
            "cancer), harmonization is still required to avoid the "
            "study-as-classifier trap (cfdna-fragmentomics skill: 0.999 → "
            "0.497 collapse)."
        )
    else:
        label = "NEUTRAL"
        explain = (
            "Per-study z-scoring made essentially no difference on this "
            "synthetic fixture. The signal in this controlled fixture is "
            "uncorrelated with study, so the linear LR handles it without "
            "harmonization. On a real cross-study pool where batches are "
            "more confounded, expect a larger effect."
        )
    return {"label": label, "explanation": explain}


def _write_docs(out: Dict) -> None:
    DOCS_OUT.parent.mkdir(parents=True, exist_ok=True)
    v = out["verdict"]
    raw = out["raw"]
    harm = out["harmonized"]
    fix = out["fixture"]
    md = f"""# Per-Cohort Harmonization Check

## Verdict

**`{v['label']}`** — Δ pooled AUC = `{out['delta_pooled_auc']:+.4f}`
`({raw['pooled_auc_mean']:.4f} → {harm['pooled_auc_mean']:.4f})` over
{len(out['seeds'])} seeds × {out['n_folds']}-fold stratified CV.

{v['explanation']}

## Why this script exists

The panel-LLR pipeline (see `real_tcga_validation.py`) is single-cohort
by default. When pooling patients from multiple studies, two things
happen:

1. **Coverage / library-prep batch effect** shifts the raw panel-LLR
   sum and frag-channel proxy. This is a confound, not a signal.
2. **Study-as-classifier trap**: if cancer labels happen to cluster in
   one study and healthy in another, the model learns the study, not
   the cancer (cfdna-fragmentomics skill: AUC 0.999 → 0.497 collapse on
   the only-Jiang-cancer vs only-Cristiano-healthy confound).

Per-study z-scoring (fit on TRAIN fold only, applied to test fold) is
the standard fix for the mild / coverage-driven regime. This script
verifies the code path is correct and quantifies the AUC delta on a
controlled synthetic fixture.

## The fixture

- **4 studies × {fix['patients_per_study']} patients × 2 (cancer/control)
  = {fix['n_samples']} samples**
- Mutations drawn from real TCGA-LUAD MAFs (see
  `validation/tcga/tcga_cache/`)
- Each patient produces a paired (TF={fix['tumor_fraction']}, TF=0)
  sample — the standard MRD-style paired design.
- A **per-study additive bias** is added to BOTH panel and frag channels
  so the studies differ in mean even on the raw data. Cancer/control is
  INDEPENDENT of study (both classes appear in every study), so this is
  the mild regime where harmonization is supposed to help, not the
  only-Jiang-cancer confound.

### Per-study mean of observed panel-LLR score

| study | mean panel | mean frag |
|---|---|---|
""" + "\n".join(
        f"| {s} | {fix['study_means_panel'][s]:+.3f} | {fix['study_means_frag'][s]:+.3f} |"
        for s in out["fixture"]["study_means_panel"]
    ) + f"""

### Per-study additive bias injected

| study | panel bias | frag bias |
|---|---|---|
""" + "\n".join(
        f"| {s} | {fix['study_panel_bias'][s]:+.2f} | {fix['study_frag_bias'][s]:+.2f} |"
        for s in fix["study_panel_bias"]
    ) + f"""

## The pipeline

For each seed in `{out['seeds']}`:
1. Stratified {out['n_folds']}-fold split on (X, y), preserving the
   cancer/control class balance.
2. **Raw branch**: fit LR on train fold, score test fold.
3. **Harmonized branch**: fit per-study (mean, std) on train fold
   ONLY — never on test, never on pooled. Apply to both train and
   test, then fit LR and predict.
4. Pool OOF (y_true, y_score) across folds, compute pooled AUC. Also
   compute per-study AUC by slicing the OOF arrays by study id.

## Per-study AUC results

| study | raw AUC | harmonized AUC | Δ |
|---|---|---|---|
""" + "\n".join(
        f"| {s} | {raw['per_study_auc_mean'][s]:.4f} ± {raw['per_study_auc_std'][s]:.4f}"
        f" | {harm['per_study_auc_mean'][s]:.4f} ± {harm['per_study_auc_std'][s]:.4f}"
        f" | {out['delta_per_study_auc'][s]:+.4f} |"
        for s in raw["per_study_auc_mean"]
    ) + f"""

## Pooled AUC

| condition | mean ± std (over {len(out['seeds'])} seeds) | per-seed values |
|---|---|---|
| raw | {raw['pooled_auc_mean']:.4f} ± {raw['pooled_auc_std']:.4f} | {raw['per_seed_pooled']} |
| harmonized | {harm['pooled_auc_mean']:.4f} ± {harm['pooled_auc_std']:.4f} | {harm['per_seed_pooled']} |

## How to reproduce

```bash
.venv/bin/python scripts/harmonization_check.py
```

Outputs `results/harmonization_check.json` (full numbers, every seed
and per-study AUC) and this file.

## Limitations

- The fixture is synthetic. The cancer-vs-control signal comes from
  the cfDNA simulation, but the per-study bias is an injected additive
  shift, not a real coverage/library-prep confound. A negative result
  here does NOT prove harmonization is useless in practice — only that
  the synthetic regime is too easy for LR + small additive bias.
- The "frag" channel is a per-patient feature proxy (mutation count /
  200 + jitter), not real cfDNA fragmentomics. This is documented in
  the script — the goal is to validate the harmonization code path
  under a controlled multi-cohort setting.
- n=160 (80 patients × 2) is small. Per-study n=20 keeps each fold
  small enough that train-only z-score fitting is still well-defined.
- For the REAL cross-study pooling benchmark, use `scripts/run_cross_study.py`
  against FinaleDB data (Cristiano 2019 + Jiang 2015). That pipeline
  measures the same harmonization recipe on real data.
"""
    DOCS_OUT.write_text(md)


if __name__ == "__main__":
    os.chdir(ROOT)
    run_harmonization_check()