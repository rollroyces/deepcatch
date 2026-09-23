"""Per-cohort batch-effect harmonization check for the panel-LLR pipeline.

Builds a SYNTHETIC multi-cohort fixture (4 studies × 20 patients × paired
cancer/control = 160 samples) purely in numpy — no real mutation files, no
MAF cache, no network. A per-study mean shift in panel-LLR scores is
injected so the studies are mildly confounded with the cancer label,
mimicking a coverage / library-prep batch effect. Runs the LR-on-LLR
pipeline twice on the same fixture:

  (a) raw          — current `real_tcga_validation.py` default (no harmonization)
  (b) harmonized   — per-study z-scoring fit on the TRAINING fold only,
                     then applied to test fold

Honest framing: this is a synthetic fixture, not a real cross-cohort
validation. It exists to verify that the harmonization code path is
correct (z-scoring on TRAIN fold only — never on test, never on pooled
data) and to measure the AUC delta in a controlled setting. Real
cross-study pooling should still go through `scripts/run_cross_study.py`
when FinaleDB data is available.

The fixture is intentionally cheap to build (~10 ms with numpy only) so
the end-to-end tests in `test/test_harmonization_check.py` can run on a
developer laptop in seconds without needing the 80 MB TCGA-LUAD MAF
cache. CI behavior is unchanged: the unit tests for the z-score
contract never required a cache.

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

# ── Constants ─────────────────────────────────────────────────────────────
RESULTS_OUT = ROOT / "results" / "harmonization_check.json"
DOCS_OUT = ROOT / "docs" / "HARMONIZATION.md"

STUDIES = ["study_A", "study_B", "study_C", "study_D"]
N_PATIENTS_PER_STUDY = 20
# Conceptual reference: paired-design tumor_fraction = 0.001 (0.1%), the
# MRD-style operating point. Kept for documentation only — build_fixture()
# no longer simulates cfDNA so this value is no longer read at runtime.
TUMOR_FRACTION = 0.001
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

# Cancer signal magnitude on the panel-LLR sum (raw units). Set so that the
# pooled cancer-vs-control AUC sits in the ~0.99 regime once the LR has the
# features — the bias is intentionally a small fraction of the signal so
# the verdict remains NEUTRAL (per-study z-scoring barely moves the pooled
# metric when the bias is uncorrelated with the cancer label).
PANEL_CANCER_MEAN = 280.0
PANEL_NOISE_STD = 80.0
FRAG_CANCER_MEAN = 0.50
FRAG_NOISE_STD = 0.10
FRAG_JITTER_STD = 0.05


# ── Fixture generation ──────────────────────────────────────────────────

def build_fixture(seed: int = 42) -> Dict[str, np.ndarray]:
    """Build the multi-cohort synthetic fixture in pure numpy.

    Returns a dict with arrays:
      panel_scores  (n_samples,)  — synthetic panel-LLR sum per sample
      frag_scores   (n_samples,)  — synthetic frag-channel proxy per sample
      study         (n_samples,)  — study id per sample (encoded as int)
      y             (n_samples,)  — 1 = cancer, 0 = control
      sample_id     (n_samples,)  — str "<study>_<patient_idx>_<pos|neg>"
      study_names   list[str]     — STUDIES (preserves the public contract)

    Design (deterministic given `seed`):

      For each study, N_PATIENTS_PER_STUDY patients × 2 labels (cancer +
      matched control at TF=0):
        panel_raw  = N(PANEL_CANCER_MEAN, PANEL_NOISE_STD) if cancer
                   else N(0, PANEL_NOISE_STD)                  (paired control)
        frag_raw   = N(FRAG_CANCER_MEAN, FRAG_NOISE_STD) if cancer
                   else N(0, FRAG_NOISE_STD)                  (paired control)
      Per-study additive bias shifts the mean (the coverage/library-prep
      batch effect). Cancer label is INDEPENDENT of study — both classes
      appear in every study — so this is the "mild batch effect" regime
      where the LR handles the bias on its own and the verdict is NEUTRAL.

    No MAF cache, no `load_tcga_cohort`, no network. End-to-end cost is
    dominated by `StratifiedKFold` over 160 rows, which is sub-second.
    """
    rng = np.random.default_rng(seed)
    n_per_study = N_PATIENTS_PER_STUDY

    panel_scores: List[float] = []
    frag_scores: List[float] = []
    y_arr: List[int] = []
    study_arr: List[int] = []
    sample_id: List[str] = []

    for s_idx, study in enumerate(STUDIES):
        # Per-study batch-effect bias (the coverage / library-prep confound).
        bias_panel = STUDY_PANEL_BIAS[study]
        bias_frag = STUDY_FRAG_BIAS[study]

        # Draw N_PATIENTS_PER_STUDY patients × 2 arms (cancer, control) per study.
        # cancer_raw:     N(PANEL_CANCER_MEAN, PANEL_NOISE_STD)
        # control_raw:    N(0, PANEL_NOISE_STD)               ← paired at TF=0
        cancer_panel = rng.normal(PANEL_CANCER_MEAN, PANEL_NOISE_STD,
                                  size=n_per_study)
        control_panel = rng.normal(0.0, PANEL_NOISE_STD, size=n_per_study)
        cancer_frag = rng.normal(FRAG_CANCER_MEAN, FRAG_NOISE_STD,
                                 size=n_per_study)
        control_frag = rng.normal(0.0, FRAG_NOISE_STD, size=n_per_study)

        # Tiny within-study jitter so LR has some per-sample separability
        # beyond the bias shift (mirrors the per-arm jitter used in
        # scripts/foundation_real_smoke.py).
        frag_jitter = rng.normal(0.0, FRAG_JITTER_STD, size=2 * n_per_study)

        for i in range(n_per_study):
            # cancer sample
            panel_scores.append(cancer_panel[i] + bias_panel)
            frag_scores.append(cancer_frag[i] + bias_frag + frag_jitter[2 * i])
            y_arr.append(1)
            study_arr.append(s_idx)
            sample_id.append(f"{study}_{i:03d}_pos")
            # matched control at TF=0
            panel_scores.append(control_panel[i] + bias_panel)
            frag_scores.append(control_frag[i] + bias_frag + frag_jitter[2 * i + 1])
            y_arr.append(0)
            study_arr.append(s_idx)
            sample_id.append(f"{study}_{i:03d}_neg")

    fixture = {
        "panel_scores": np.asarray(panel_scores, dtype=float),
        "frag_scores": np.asarray(frag_scores, dtype=float),
        "y": np.asarray(y_arr, dtype=int),
        "study": np.asarray(study_arr, dtype=int),
        "sample_id": np.asarray(sample_id),
        "study_names": STUDIES,
    }
    print(
        f"[1] Fixture built (pure numpy, no MAF cache): n={len(fixture['y'])} "
        f"samples, {len(STUDIES)} studies × {N_PATIENTS_PER_STUDY} patients "
        f"× 2 (cancer/control)"
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
            f"panel mean={ps.mean():+.2f} std={ps.std():.2f}  "
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
            # Synthetic fixture — record the params actually used by build_fixture.
            "panel_cancer_mean": PANEL_CANCER_MEAN,
            "panel_noise_std": PANEL_NOISE_STD,
            "frag_cancer_mean": FRAG_CANCER_MEAN,
            "frag_noise_std": FRAG_NOISE_STD,
            "frag_jitter_std": FRAG_JITTER_STD,
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
- **Pure numpy synthesis** (no TCGA MAF cache, no network). Cancer panel-LLR
  is drawn from `N(panel_cancer_mean={fix['panel_cancer_mean']:.2f},
  panel_noise_std={fix['panel_noise_std']:.2f})`, matched controls from
  `N(0, {fix['panel_noise_std']:.2f})` — the standard MRD-style paired
  design (cancer at TF={TUMOR_FRACTION}, control at TF=0).
- The frag channel uses the same paired design with `frag_cancer_mean=
  {fix['frag_cancer_mean']:.3f}` and a per-sample jitter of
  ±{fix['frag_jitter_std']:.3f}.
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

- The fixture is fully synthetic — pure numpy draws, no real mutations or
  cfDNA simulation. The per-study bias is an injected additive shift, not
  a real coverage/library-prep confound. A NEUTRAL verdict here does NOT
  prove harmonization is useless in practice — only that the synthetic
  regime is too easy for LR + small additive bias that is uncorrelated
  with the cancer label.
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