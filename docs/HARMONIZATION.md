# Per-Cohort Harmonization Check

## Verdict

**`NEUTRAL`** — Δ pooled AUC = `+0.0032`
`(0.9968 → 1.0000)` over
5 seeds × 5-fold stratified CV.

Per-study z-scoring made essentially no difference on this synthetic fixture. The signal in this controlled fixture is uncorrelated with study, so the linear LR handles it without harmonization. On a real cross-study pool where batches are more confounded, expect a larger effect.

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

- **4 studies × 20 patients × 2 (cancer/control)
  = 160 samples**
- **Pure numpy synthesis** (no TCGA MAF cache, no network). Cancer panel-LLR
  is drawn from `N(panel_cancer_mean=280.00,
  panel_noise_std=80.00)`, matched controls from
  `N(0, 80.00)` — the standard MRD-style paired
  design (cancer at TF=0.001, control at TF=0).
- The frag channel uses the same paired design with `frag_cancer_mean=
  0.500` and a per-sample jitter of
  ±0.050.
- A **per-study additive bias** is added to BOTH panel and frag channels
  so the studies differ in mean even on the raw data. Cancer/control is
  INDEPENDENT of study (both classes appear in every study), so this is
  the mild regime where harmonization is supposed to help, not the
  only-Jiang-cancer confound.

### Per-study mean of observed panel-LLR score

| study | mean panel | mean frag |
|---|---|---|
| study_A | +143.461 | +0.440 |
| study_B | +134.995 | +0.203 |
| study_C | +146.424 | +0.308 |
| study_D | +142.201 | +0.086 |

### Per-study additive bias injected

| study | panel bias | frag bias |
|---|---|---|
| study_A | +0.40 | +0.20 |
| study_B | +0.10 | -0.05 |
| study_C | -0.10 | +0.05 |
| study_D | -0.30 | -0.15 |

## The pipeline

For each seed in `[42, 123, 456, 789, 1024]`:
1. Stratified 5-fold split on (X, y), preserving the
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
| study_A | 0.9965 ± 0.0014 | 1.0000 ± 0.0000 | +0.0035 |
| study_B | 1.0000 ± 0.0000 | 1.0000 ± 0.0000 | +0.0000 |
| study_C | 0.9920 ± 0.0048 | 1.0000 ± 0.0000 | +0.0080 |
| study_D | 1.0000 ± 0.0000 | 1.0000 ± 0.0000 | +0.0000 |

## Pooled AUC

| condition | mean ± std (over 5 seeds) | per-seed values |
|---|---|---|
| raw | 0.9968 ± 0.0012 | [0.995, 0.99796875, 0.99734375, 0.99734375, 0.99625] |
| harmonized | 1.0000 ± 0.0000 | [1.0, 1.0, 1.0, 1.0, 1.0] |

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
