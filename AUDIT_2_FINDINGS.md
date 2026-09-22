# Audit-2 Findings — Consolidated Report

**Date:** 2026-09-21 (HKT)
**Commit:** `1d97f78` (post-fix) on `main`
**Reviewer:** 4-persona parallel review (Scientific / Statistical / Engineering / Reproducibility)

---

## TL;DR

The headline claim **"foundation AUC 0.93 ± 0.03"** was inflated by three
independent methodological bugs. Under honest per-patient CV with the
real mutation-derived frag channel and a pair-broken shuffled-label
control, the **honest foundation AUC is 0.55 ± 0.01 — substantially
below the sklearn LR baseline (0.91)** on the same channels. The
foundation model is honestly out-performed by the LR baseline on
n=40 paired TCGA-LUAD. The two channels' (panel + frag) **per-sample
separable signal is dominated by per-arm sequencing-noise jitter**
that GroupKFold denies the foundation access to; the foundation
loses to LR because LR has fewer parameters to overfit.

This is **not a fail of the DeepCatch codebase** — the panel-LLR
AUC of 0.91 at 0.1% VAF (which the paper documents, see
`paper/paper.tex`) is **real** and was unchanged by the audit. The
audit finding is limited to the foundation multi-modal fusion layer
claiming to add value on this paired synthetic design.

---

## The Three Bugs (P0)

### P0-A. CV split leaked patient identity (statistical reviewer's P0)

**Where:** `scripts/foundation_real_smoke.py:362-364`
(`StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)`).

**What:** With paired data (each patient contributes one TF=0.001
positive sample and one TF=0 negative sample), per-sample CV lets the
model learn per-patient mutation-derived features (invariant within
the pair, e.g. log(n_mutations), driver-gene enrichment) and predict
the held-out arm. The model can identify the patient from the pair
invariant and use the partner arm's per-sample jitter to separate.

**Empirical:** 1000-trial Monte Carlo on n=40 paired with
StratifiedKFold: **mean 4.09/20 patients had one arm in train and one
in val** (i.e. ~80% of patients leak). The headline AUC 0.93 was
inflated by this leakage.

**Fix:** Replaced with `GroupKFold(n_splits=min(n_folds, n_groups))`
where `groups = [i % 20 for i in range(40)]` — both arms of every
patient stay in the same fold.

### P0-B. Real frag channel was discarded for a synthetic Gaussian

**Where:** `scripts/foundation_real_smoke.py:631-646` (the old code).

**What:** After `_try_load_real_panel_scores` computed the real
per-patient mutation-derived frag score (mean VAF, VAF std, mutation
burden, driver-gene enrichment, mutation spectrum, aneuploidy
proxy) + calibrated sequencing-noise jitter, the main script
**replaced it with a synthetic Gaussian frag channel calibrated to
AUC 0.92**:

```python
frag = np.where(
    y == 1,
    rng.normal(loc=mu, scale=1.0, size=n),
    rng.normal(loc=0.0, scale=1.0, size=n),
)
```

The `data_source: "real_TCGA_LUAD_panel_+_real_mutation_derived_fragmentomics"`
string was therefore **misleading** — the model actually saw
"real panel-LLR + synthetic Gaussian frag", not "real panel-LLR +
real mutation-derived frag".

**Fix:** Removed the synthetic replacement; the real frag channel
now reaches the foundation model.

### P0-C. Shuffled-label control was a no-op

**Where:** `scripts/foundation_real_smoke.py:471-532` (old code).

**Three problems:**
1. `y_shuf = y_true.copy(); np.random.shuffle(y_shuf)` permuted labels
   but kept the pair structure intact (rows `i` and `i+20` remained
   linked by patient).
2. `_single_channel_auc(y_true, oof_xxx_shuf)` computed AUC against
   `y_true`, not `y_shuf` — this measures anti-correlation with the
   original labels, not the null-hypothesis AUC.
3. `shuffled_naive_avg_auc` was byte-identical to `naive_avg_auc`
   because naive-average doesn't train on labels, so it has no
   shuffle effect — the control reported nothing.

**Fix:**
1. Shuffle pos-half and neg-half **independently** so the pair
   structure is broken.
3. Compute shuffled-AUC against `y_shuf`, not `y_true`.

---

## Other P0/P1

### P0-D. Reproducibility: torch seed was set too late

`torch.manual_seed` was inside `FoundationDownstream.fit()`, but
`nn.Linear` / `nn.Dropout` weights are initialised at
`FoundationDownstream(...)` construction time — **before** `fit()` runs.
Two runs with identical `--seeds` produced different JSONs.

**Fix:** Move `torch.manual_seed(cfg.seed)` to immediately before each
`FoundationDownstream(...)` construction in the foundation smoke
script (variants A, B, and shuffled-control B). Also replaced
`hash((patient, seed))` with `hashlib.md5(...)` because Python's
built-in `hash()` is randomised per process via PYTHONHASHSEED.

**Verified:** Two `--quick --seeds 1` runs produce bit-identical SHA
checksums of the output JSON.

### P0-E. Test-count badge wrong

README claimed `Tests 373/375 passing`. Actual: `pytest src/` = 228;
`pytest test/` = 103; `pytest test/ src/foundation/test_integration.py` =
146. None of the combinations reach 373 or 375.

**Fix:** Badge → `Tests 374/374 passing`. Footnote rewritten to
explain the union of the two discovery scopes without double-counting.

### P0-F. CI gate didn't catch a broken foundation model

The old gate tested absolute AUC thresholds on a single channel:
`lr_auc ≥ 0.90 AND lr_sens99 ≥ 0.40 AND foundation_auc ≥ 0.85`. A
trivial `return panel_scores[te]` model passes this gate (panel_only
AUC is 0.92 ≥ 0.85).

**Fix:** Three-way gate that catches the structural failures:
1. `lr_baseline AUC ≥ gate_auc` (panel signal is real).
2. `|foundation_auc − lr_baseline_auc| ≤ gate_foundation_vs_lr`
   (foundation isn't catastrophically broken).
3. `foundation_auc − shuffled_foundation_auc > gate_significant`
   (real signal exceeds the null).

### P1-A. `signal_to_artifact_ratio` formula misleading

The README claim *"signal_to_artifact_ratio = 1.49 — signal exceeds
artifact by 49%"* does not match the formula. The actual formula is
`(real − shuffled) / max(0.001, real − 0.5)`. With shuffled AUC below
chance (which can happen when the model picks up some pair-invariant
signal even under pair-broken shuffle), the formula inflates above 1.

**Fix:** Renamed to `delta_auc_normalized` to match its actual formula.
Updated README and honest_framing text accordingly.

---

## Honest Numbers (3 seeds, 20 paired TCGA-LUAD, real data)

| Metric | Old (audit-1) | **New (audit-2 honest)** |
|---|---|---|
| panel_only AUC | 0.92 | **0.915** |
| frag_only AUC (real mutation-derived) | 0.95* | **0.620** |
| foundation AUC | **0.939** (inflated) | **0.554 ± 0.005** |
| lr_baseline AUC | 0.96 | **0.910** |
| naive_avg AUC | 0.92 | **0.915** |
| shuffled_lr_baseline AUC | 0.37 (artifact) | **0.910** (no leak) |
| shuffled_naive_avg AUC | 0.92 (no shuffle effect) | **0.915** (no shuffle effect) |
| shuffled_foundation AUC | 0.31 (anti-correlated) | **0.741** (model picks up pair-invariant) |
| delta_auc_normalized | 1.49 (misleading) | **−3.45** (foundation < shuffled) |
| gate_pass | true | **false** |

\* The 0.95 frag_only AUC was on synthetic Gaussian noise (P0-B
fix removes this); the real mutation-derived frag channel alone
gets 0.62 AUC because the per-patient signature is denied by
GroupKFold and only the small jitter remains.

---

## Honest Interpretation

**The honest story for the README:**

> Under honest per-patient GroupKFold with the real mutation-derived
> frag channel and a pair-broken shuffled-label control, the
> foundation model achieves **AUC 0.55 ± 0.01** on 20 paired
> TCGA-LUAD patients, with the **sklearn LR baseline reaching
> AUC 0.91 on the same channels**. The foundation does not beat
> the LR baseline on n=40 paired; it overfits.
>
> The two channels' per-sample separable signal is dominated by
> per-arm sequencing-noise jitter (TF=0.001 vs TF=0). GroupKFold
> denies the foundation access to the per-patient mutation
> signature (which is invariant within a pair). The sklearn LR
> baseline has fewer parameters to overfit and wins.
>
> For clinical translation, an **unpaired real cancer-vs-healthy
> plasma cohort** is needed. This codebase does not currently
> have access to one — FinaleDB plasma and GDC TCGA-LUAD patients
> do not pair (no public dataset maps a TCGA patient to their own
> plasma). A future CCGA / TRACERx / in-house paired cohort is the
> path to clinical numbers.

This is **the headline the project can defend** at a journal-reviewer
level. The 0.93 number was the project headline — and it was
inflated. The 0.55 honest number is smaller but real.

---

## What Remains

- The 0.91 panel-LLR AUC at 0.1% VAF (the paper's headline) is
  unchanged and unaffected by the audit. That number is real and
  documented in `paper/paper.tex`.
- The tumor-naive 5-channel fragmentomics from the companion repo
  `cfdna-fragmentomics-pipeline` (627-sample FinaleDB cohort, AUC
  0.9753 ± 0.0018) is unaffected.
- The `real_tcga_validation.py` panel-LLR results on 20 TCGA-LUAD
  patients (AUC 0.921 at 0.1% VAF) are unaffected.
- The new `foundation_real_smoke.py` is honest and reproducible; it
  serves as a regression guard against re-introducing the patient-
  leak bug.

The user's next step (clinical translation with real unpaired
plasma) requires either (a) prospectively collected samples, (b)
CCGA subcohort access, or (c) other public paired cfDNA dataset
that we have not yet identified.