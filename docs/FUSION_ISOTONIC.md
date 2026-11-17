# Isotonic-Calibrated LR Fusion — Honest Report

**Date:** 2026-09-19
**Branch:** `main` @ 86a57e2
**Cohort:** 627 cross-study samples (363 cancer / 264 healthy),
pan-cancer FinaleDB + DELFI-derived fragmentomics features
(`features-dir=cfdna-fragmentomics-pipeline/data/features`,
`pca_n=200`, `target_auc=0.92` for the synthetic mutation channel).
**Seeds:** 10.

## What we tested

`src/fragmentomics/fusion_ablation.py` already had four fusion
strategies:

| Strategy | Method |
|---|---|
| `tumor_naive` | LR on PCA(features) — single channel baseline |
| `mutation_only` | Synthetic calibrated LLR (AUC ≈ 0.92) — single channel |
| `naive_average` | `(tn_score + mut_score) / 2` |
| `lr_fusion` | LR on `[tn_score, mut_score]`, per-fold trained |
| `lr_fusion_isotonic` | **NEW**: LR on `[tn_score, mut_score]`, then per-fold isotonic calibration fit on the training fold only |

The new method's per-fold pipeline:

1. Fit LR on `(tn_tr, mut_tr) → y_tr`.
2. Score the **training** fold in-sample → `lr_tr_scores`.
3. Fit `sklearn.isotonic.IsotonicRegression(y_min=0, y_max=1,
   out_of_bounds="clip")` on `(lr_tr_scores, y_tr)`.
4. Score the **test** fold → `lr_te_scores`.
5. Apply the isotonic mapping to `lr_te_scores` → calibrated scores.

Step (3) uses only training-fold data, so there is no test-side
leakage. Step (2) is in-sample *on purpose*: the isotonic step is
fitting the calibration curve, not the LR weights.

## Results (10 seeds, 5-fold CV, pooled OOF)

Source: `results/fusion_ablation_isotonic.json`.

| Strategy | AUC (mean ± std) | Sens@95 | Sens@99 |
|---|---|---|---|
| Tumor-naive only | 0.9740 ± 0.0013 | 0.880 | 0.742 |
| Mutation-only (calibrated AUC 0.92) | 0.9020 ± 0.0000 | 0.595 | 0.267 |
| Naive average | **0.9884 ± 0.0010** | **0.927** | **0.841** |
| LR fusion (uncalibrated) | **0.9884 ± 0.0009** | **0.932** | 0.833 |
| **LR fusion (isotonic)** | **0.9849 ± 0.0014** | 0.928 | 0.646 |

### DeLong (seed 0, vs tumor_naive)

| Strategy | Δ AUC | z | p (two-sided) | significant @ 0.05 |
|---|---|---|---|---|
| Mutation-only | −0.0709 | −5.30 | 1.2e−07 | yes (worse) |
| Naive average | +0.0145 | 3.63 | 0.00028 | yes |
| LR fusion | +0.0144 | 3.53 | 0.00042 | yes |
| LR fusion (isotonic) | +0.0103 | 2.09 | 0.036 | marginal |

### Per-seed AUC for LR vs LR-isotonic

Every single seed shows an AUC regression for isotonic, ranging from
−0.0014 (seed 1) to −0.0051 (seed 6):

```
seed 0:  lr=0.9872   iso=0.9831   delta=-0.0041
seed 1:  lr=0.9882   iso=0.9868   delta=-0.0014
seed 2:  lr=0.9889   iso=0.9863   delta=-0.0026
seed 3:  lr=0.9882   iso=0.9846   delta=-0.0036
seed 4:  lr=0.9882   iso=0.9840   delta=-0.0043
seed 5:  lr=0.9873   iso=0.9850   delta=-0.0023
seed 6:  lr=0.9878   iso=0.9827   delta=-0.0051
seed 7:  lr=0.9898   iso=0.9850   delta=-0.0048
seed 8:  lr=0.9883   iso=0.9843   delta=-0.0040
seed 9:  lr=0.9901   iso=0.9869   delta=-0.0032
```

## Honest bottom line

**Isotonic calibration did NOT help, and in fact it hurt.** Specifically:

1. **AUC regressed by ~0.35 percentage points** (0.9884 → 0.9849),
   uniformly across all 10 seeds.
2. **Sens@99 regressed by ~19 percentage points** (0.833 → 0.646),
   and in two seeds (6 and 7) isotonic produced **Sens@99 = 0.000** —
   the step function collapsed the test-fold scores onto a single
   value, eliminating high-specificity detection entirely.
3. **DeLong significance against tumor-naive dropped from p=0.00042
   (uncalibrated LR) to p=0.036 (isotonic)** — isotonic made the
   fusion's gain over the single-channel baseline barely detectable.

### Why did this happen?

LR fusion on two already-[0,1]-bounded scores (tumor-naive
`predict_proba` and the sigmoid-mapped mutation LLR) produces scores
that are *already well-calibrated* in the sense that matters here:
they live on a comparable probability scale, and the empirical
sensitivity at fixed specificity is dominated by the ranking, not
the absolute scores. Isotonic regression on such inputs tends to
introduce **flat regions** (where the in-sample LR score doesn't
change much but the empirical class fraction does) and these flat
regions can collapse an entire test fold onto one calibrated value
when the test fold has narrow score variance — exactly what we saw
in seeds 6 and 7.

In other words: the LR fusion's output is already a monotone
function of the underlying class-posterior signal, so re-fitting a
monotone mapping on a single training fold removes information
rather than adding it.

### Recommendation

**Keep the naive average (or the uncalibrated LR fusion) as the
default.** The README's existing recommendation — that equal-weight
fusion of two well-calibrated scores is already optimal in this
regime — is **vindicated** by this experiment, not overturned.

The `lr_fusion_isotonic` strategy is **kept in the code and in the
JSON output** for completeness and for future experimentation (e.g.
if DeepCatch ever ships an uncalibrated black-box mutation score, a
post-hoc calibrator may become useful). But it should not be the
default.

### Reproducing this experiment

```bash
cd deepcatch
env -u PYTHONPATH .venv/bin/python -m src.fragmentomics.fusion_ablation \
    --features-dir ../cfdna-fragmentomics-pipeline/data/features \
    --labels ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
    --seeds 10 --pca-n 200 \
    --out results/fusion_ablation_isotonic.json
```

Total wall-time on a 2024 Mac mini: ~6 minutes for 10 seeds × 5 folds
× 5 inner-fold + 5 strategies (per seed).

### Code & test changes

- `src/fragmentomics/fusion_ablation.py`:
  - new function `_fusion_lr_isotonic_score(...)` (line 178+)
  - `_evaluate_seed` now appends `lr_fusion_isotonic` per fold and
    runs DeLong against it
  - the `strategies` list (and module docstring) now lists all 5
    strategies

- `test/test_fusion_ablation.py`:
  - imports the new function
  - updated `test_evaluate_seed_returns_all_strategies` to require
    all 5 strategies and 4 DeLong comparisons
  - new `test_isotonic_returns_same_shape_as_lr_fusion`
  - new `test_isotonic_is_monotone_in_raw_lr_scores`
  - new `test_cli_with_seeds1_emits_all_four_fusion_methods`

- `README.md` / `README.zh-CN.md`: test badge bumped 67/67 → 70/70,
  fusion section updated to mention `lr_fusion_isotonic` row.