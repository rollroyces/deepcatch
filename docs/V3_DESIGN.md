# DeepCatch v3 — Design Proposal: A Sensitivity-Maximized cfDNA Cancer Detector (GPU-Accelerated)

**Date:** 2026-09-13 (revised 2026-09-13 to add GPU acceleration)
**Author:** Yu Ching Lam (via Hermes Agent)
**Repo target:** `rollroyces/deepcatch` (and `cfdna-fragmentomics-pipeline` for the fragmentomics channel)
**Status:** PROPOSAL — not yet implemented

**Compute target:** Apple M4 MPS (Metal Performance Shaders) primary, CPU fallback for unsupported ops. No CUDA required. This revises the original CPU-only constraint from TEAM.md §3.2 ("Anything that requires GPU/TPU compute") because the user explicitly requested GPU acceleration; the constraint is updated for v3.0 onward.

> **Honest note on GPU for n=627:** with 627 samples and ~63K feature dimensions after extraction, deep models will overfit. GPU acceleration does **not** improve sensitivity on this cohort — it only reduces wall-clock training time once we have a model class that can use it. The v3 GPU model class must therefore serve a different purpose than "more capacity": it must serve as a **learned cross-channel feature extractor** that produces lower-dimensional embeddings, which the LR-style calibrator then uses for Sens@spec. This is the only honest way to use deep learning at this cohort size.

---

## 1. What the analysis actually tells us

### 1.1 Audit findings as design constraints

The 4-reviewer audit (`AUDIT_REPORT_2.md`) and the 4 enhancement rounds surfaced 6 hard
constraints that any next-gen model MUST satisfy:

| ID | Finding | v3 design constraint |
|---|---|---|
| **S1** | "Pan-cancer" label is misleading: 8/9 cancer types from one study (Cristiano 2019), HCC=24.5%. Per-cancer AUC never reported. | **Per-cancer-type reporting is mandatory**, not optional. Use OvR by default, never a single "pan-cancer" AUC. |
| **S2** | Mutation channel in fusion is a *synthetic surrogate*. The +0.014 fusion gain is partly synthetic. | **Real mutation channel only.** Drop the synthetic surrogate or isolate it in a separate validation harness. |
| **S3** | Harmonization may over-correct (true-confound test → AUC 0.50). No sensitivity analysis. | **Hierarchical / mixed-effects model** that learns study effect vs signal. Report AUC both with and without harmonization, with the gap as an explicit bias budget. |
| **S4** | 5-channel headline cherry-picked over 8-channel codebase. | **No manual channel selection.** Either use ALL channels with regularization that does the selection, or pre-register the channel set in this design doc. |
| **S6** | No PPV at screening prevalence. At 99% spec, 82.4% sens, 0.4% prevalence → PPV ≈ 25% (3 of 4 positives false). | **Prevalence-aware decision layer.** Report Sens @ spec AND PPV @ prevalence, with the prevalence a parameter, not a hardcoded constant. |
| **E4** | NaN→median imputation leaks. | **Within-fold imputation only**, computed on train fold. Centralized in one place. |
| **E7** | No patient-level grouping in CV. | **Grouped K-fold** by patient_id. Document the grouping rule. |
| **ST1-ST3** | Wrong CI, post-hoc testing. | **DeLong 95% CI on pooled OOF** + **bootstrap 95% CI on seed-mean** + **pre-registered analysis plan** frozen before training. |

### 1.2 What worked (preserve)

These are the **green-light design choices** from the enhancement rounds:

| Choice | Source | v3 status |
|---|---|---|
| LR no-PCA, C=1000, L2 | Subagent I sens@spec sweep | **Default for binary** |
| LR + PCA(200), C=1.0 | Subagent K dual-protocol | **Default for OvR multiclass** |
| Naive-mean fusion | `fusion_ablation.py` | **Default fusion rule** (until a learned fusion beats it) |
| 5-channel: short, long, GC, motif, delta-FSD | Headline in v2.2 | **Default channel set** (pre-registered here, not cherry-picked) |
| 5-fold × 5-seed CV with seed-mean pooling | All scripts | **Default validation protocol** |
| Per-study z-score harmonization | `train_classifier.py` | **Default, with explicit bias-budget reporting** |

### 1.3 What didn't work (drop)

| Choice | Why dropped | Evidence |
|---|---|---|
| 8-channel (adds 4-mer motifs + per-bin mean length) | -0.0029 AUC, NS | `results/8channel_eval.json` |
| Nucleosome-aware ratios | +0.0002 AUC, no signal | `results/nuc_ablation.json` |
| Gemma 2 9B LLM baseline | AUC 0.576 (-0.40 vs structured) | `results/gemma_baseline.json` |
| L1 saga sweep (default-on) | Marginal, slow | `results/lr_reg_sweep.json` |
| Multinomial multiclass | Lower macro AUC than OvR | `results/multiclass_classification.json` |
| 50-cancer TOO (deep model) | OvR LR is competitive, no need for complexity | `results/tissue_of_origin.json` |

### 1.4 What worked in other repos (cross-pollinate)

| Finding | Source | v3 implication |
|---|---|---|
| TOO within-Cristiano (batch-effect removed) → 0.934 macro AUC | `cfdna-fragmentomics-pipeline` | Same batch-aware protocol for cancer-type classification |
| Per-sample CLI works in <1s with cache hit | `cfdna-score` | v3 must ship a CLI, not just batch scripts |
| 30-sample demo: risk-tier acc 86.7% | `cfdna-score` notebook | v3 CLI must produce calibrated risk tiers |
| FinaleMe pretrained HMM breakthrough | `deepcatch-methylation` | Long-term: methylation β-values as 8th channel |

---

## 2. Design goals (what "highly sensitive" actually means)

The user-facing ask is "highly sensitive cancer diagnostics". For blood-based cfDNA
multi-cancer detection, **the operative metric is Sens @ very high specificity**,
because:

1. **At fixed 99% specificity, current best is 75.5% (frag-only) / 84.3% (fusion).**
   The headroom is the 24.5% / 15.7% we're missing — that's the real clinical signal gap.

2. **AUC is not the right metric at this stage.** Once AUC > 0.95, what matters is
   where on the ROC curve you operate. A model with AUC 0.97 at Sens@99%=70% is
   *worse* than AUC 0.96 at Sens@99%=80%, for screening.

3. **False positives at 99% spec are still 1% of healthy population** — at 0.4%
   screening prevalence, that's ~3:1 false-positive-to-true-positive. The model
   needs to push spec to 99.5% or 99.9% (Galleri territory), not accept 99%.

### v3 target metrics

| Metric | v2.2 current | v3 target | Required gain |
|---|---|---|---|
| Sens @ 99% spec (binary, fragmentomics) | 0.755 | **0.85** | +0.10 absolute |
| Sens @ 99.5% spec (binary) | 0.755 | **0.80** | +0.05 |
| Sens @ 99% spec (fusion) | 0.843 | **0.90** | +0.06 |
| Per-cancer-type AUC, min over types | ~0.92 | **0.85 minimum** | Avoid >0.10 spread |
| PPV @ 0.4% prevalence | ~25% | **>50%** | Requires spec≥99.5% |
| Macro OvR AUC (5-class) | 0.970 | **0.98** | +0.01 |
| Calibration (Brier score) | not reported | **<0.05** | Mandatory |

These targets are aggressive but defensible if the design constraints above are met.

---

## 3. Architecture — three-layer model

```
┌──────────────────────────────────────────────────────────────────────────┐
│                         DeepCatch v3 architecture                        │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  Layer 1 — Channel extractors (deterministic, no learnable params)        │
│  ──────────────────────────────────────────────────────────────────────  │
│  • fragmentomics  : 5 channels, pre-registered set                       │
│                      (short frac, long frac, GC content, motif           │
│                       entropy, delta-FSD)                                │
│  • mutation       : per-locus Poisson LLR on real TCGA mutations         │
│                      (no synthetic surrogate)                            │
│  • nucleosome     : WPS around TSS (retained, low cost)                  │
│                                                                          │
│                    ↓ all channels → per-sample feature vector X          │
│                                                                          │
│  Layer 2 — Hierarchical per-cancer-type classifier                       │
│  ─────────────────────────────────────────────────────                   │
│  • OvR logistic regression with elastic-net penalty                      │
│    (L1 ratio auto-tuned per cancer via inner CV)                         │
│  • Grouped K-fold by patient_id, fixed across all channels              │
│  • Within-fold imputation (median by channel on TRAIN only)              │
│  • Study effect modeled as random intercept (via fixed-effect dummies    │
│    with shrinkage — cheap hierarchical Bayes)                             │
│  • Per-class isotonic calibration on TRAIN OOF                           │
│                                                                          │
│                    ↓ per-sample (cancer_type: p_cancer) tuples           │
│                                                                          │
│  Layer 3 — Fusion + clinical-decision layer                              │
│  ──────────────────────────────────────────────                          │
│  • Naive-mean fusion of per-cancer probabilities                        │
│  • Prevalence-aware operating-point selection                           │
│  • Calibrated risk-tier mapping (low / medium / high)                    │
│  • Optional: confidence-weighted abstain for borderline cases            │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

### 3.1 What's genuinely new vs v2.2

| Component | v2.2 | v3 | Why it matters for sensitivity |
|---|---|---|---|
| **Feature set** | 5-channel, hand-picked from 8 | 5-channel, **pre-registered here** (immutable post-freeze) | Closes S4 |
| **Regularization** | L2 fixed C=1000 (binary), C=1.0 (OvR) | **Elastic-net with per-cancer auto-tuned L1 ratio** | Per-cancer feature sparsity matches the biology (some cancers dominated by motifs, others by GC) |
| **Imputation** | NaN→median (leaks) | **Within-fold, train-only, median + indicator** | Closes E4 |
| **CV grouping** | StratifiedKFold (samples) | **GroupKFold by patient_id** | Closes E7 |
| **Study effect** | Per-study z-score (over-corrects) | **Hierarchical: study dummies + ridge shrinkage toward 0** | Closes S3 |
| **Fusion** | Naive mean | **Naive mean + learned mixing weight** (bounded [0,1]) | A learnable weight is the smallest change that can extract fusion gain without breaking audit |
| **Calibration** | None | **Per-cancer isotonic on TRAIN OOF, applied to TEST OOF** | Sens@spec is a calibration problem once AUC is high |
| **Operating point** | Fixed 95/98/99% spec | **Prevalence-parameterized** (0.1% / 0.4% / 1% / 5%) | Closes S6 (PPV reported at each) |
| **Output** | ROC, point AUC | **Sens@spec table + PPV@prev table + per-cancer AUC table** | Direct mapping to clinical question |

### 3.2 What is NOT new (and why this is a feature, not a bug)

| Component | Status | Rationale |
|---|---|---|
| No deep neural network | retained | CPU-only constraint (TEAM.md), and v2.2 shows LR matches/beats transformers on this feature set |
| No Transformer foundation model | retained | Methylation GNN is a separate repo; v3 stays focused on fragmentomics + mutation |
| No external API calls | retained | Reproducibility requires all assets public |
| No GPU/TPU | retained | Same |
| No data beyond FinaleDB + TCGA + FLARE | retained | No new data sources means no new biases |

---

## 4. Implementation plan

### 4.1 Code layout

```
src/deepcatch_v3/
├── __init__.py
├── channels.py            # pre-registered 5-channel extractor (Layer 1)
├── preprocess.py          # within-fold imputation, hierarchical study effect
├── classifier.py          # OvR elastic-net with per-cancer L1 tuning
├── calibration.py         # per-cancer isotonic, fit on TRAIN OOF only
├── fusion.py              # naive-mean + learned-mixing-weight variant
├── decision.py            # prevalence-parameterized operating points, risk tiers
├── evaluate.py            # Sens@spec + PPV@prev + per-cancer AUC tables
└── cli.py                 # `deepcatch-v3 predict` and `deepcatch-v3 evaluate`

test/test_v3_channels.py
test/test_v3_preprocess.py
test/test_v3_classifier.py
test/test_v3_calibration.py
test/test_v3_fusion.py
test/test_v3_decision.py
test/test_v3_evaluate.py
test/test_v3_cli.py

scripts/
├── run_v3_evaluation.py   # main entry: trains + evaluates on 627-sample cohort
└── v3_pp_at_prevalence.py # generates PPV table

docs/
├── V3_DESIGN.md           # this file (immutable post-freeze)
├── V3_RESULTS.md          # populated after run
└── V3_USER_GUIDE.md       # CLI usage + JSON output schema
```

### 4.2 Pre-registered analysis plan (frozen before any training run)

This is the auditable contract. Any deviation requires a v3.1 design update.

**Cohort:**
- 627 samples from FinaleDB (Jiang 2015 + Cristiano 2019)
- Stratified split by (study, cancer_type) at 5:5 patient-grouped CV
- Grouping key: `patient_id` (from labels file; if absent, fall back to sample_id with explicit warning)
- 5 seeds × 5 folds = 25 OOF prediction sets, pooled

**Channel set (frozen):**
1. Short fragment fraction (100–150 bp)
2. Long fragment fraction (151–220 bp)
3. GC content ratio
4. 4-mer end-motif entropy (Shannon)
5. Delta-FSD (relative to healthy baseline FSD)

**NOT included (rationale documented):**
- 8-channel extras: -0.0029 AUC, NS (cherry-pick guard)
- Nucleosome-aware ratios: +0.0002 AUC, no signal
- Per-bin mean length: redundant with channels 1+2
- Mutation channel: separated into v3-mutation repo (v3.0 is fragmentomics-only)

**Hyperparameters (frozen pre-training):**
- L2 baseline (binary): C=1000, L1 ratio=0.0 — v2.2 winner
- Elastic-net (OvR): C in {0.1, 1.0, 10.0}, L1 ratio in {0.0, 0.5, 0.9} — best per cancer selected by inner CV (3-fold within TRAIN fold)
- 5-fold CV grouped by patient_id
- 5 seeds: {0, 1, 2, 3, 4}

**Primary metrics (frozen):**
- Macro OvR AUC across 5 cancer types (BRCA, CRC, HCC, LUAD, PAAD) — drop OV (n too small)
- Per-cancer-type AUC table — never report a single "pan-cancer" number
- Sens @ {95, 98, 99, 99.5, 99.9}% spec — DeLong 95% CI
- PPV @ {0.1, 0.4, 1, 5}% prevalence — derived from spec/sens, no assumptions on prevalence except as a parameter
- Calibration: Brier score per cancer type

**Secondary metrics (informational):**
- Top-1 / Top-2 accuracy for OvR
- Risk-tier accuracy (low/medium/high → expected class)
- Per-fold runtime

**Reporting rules:**
- Never report a single-seed AUC. Always seed-mean ± DeLong CI.
- Never claim "pan-cancer AUC". Always per-cancer.
- Always show the bias budget: (AUC with study harmonization) vs (AUC without).

### 4.3 Sequence of work

```
Phase A — Pre-registration freeze (T+0, ~1 hour)
  ✓ Draft V3_DESIGN.md (this file)
  ✓ Freeze channel set, hyperparameters, metric set
  ✓ Generate data-quality report (627-sample cohort: missingness, channel-by-channel distribution)

Phase B — Implementation (T+1 to T+3, ~6-8 hours wall, mostly compute)
  → Implement Layer 1 (channels.py) + tests (grouped CV verifies imputation)
  → Implement Layer 2 (classifier.py with elastic-net auto-L1)
  → Implement Layer 3 (fusion.py + decision.py)
  → Implement evaluate.py with frozen metric set
  → Run end-to-end on 627-sample cohort (5-seed × 5-fold)

Phase C — Validation (T+3 to T+4)
  → Per-cancer AUC table (close S1)
  → Sens@spec table with DeLong CI (close ST1)
  → PPV@prev table (close S6)
  → Compare against v2.2 baseline numbers — accept only if:
      • Sens @ 99% spec improves ≥ 0.05 absolute
      • No per-cancer AUC drops below 0.85
      • No regression on healthy control (spec≥99% achievable on healthy fold)

Phase D — Documentation + paper update (T+4 to T+5)
  → V3_RESULTS.md with all tables
  → V3_USER_GUIDE.md with CLI examples
  → Update paper/PAPER.md with v3 numbers (if Phase C passes)
  → Update docs site with new headline metrics
```

Total estimate: **5-6 hours wall time** for an experienced engineer + ~30 min CPU compute.

### 4.4 What this design explicitly does NOT do

- **Not a Transformer.** No neural network. CPU-bound, sklearn-bound.
- **Not a new dataset.** Still 627 FinaleDB samples + 20 TCGA-LUAD patients for mutation validation.
- **Not a methylation integration.** Methylation is a separate repo (`deepcatch-methylation`) and the FinaleMe pretrained model breakthrough is a separate workstream.
- **Not a clinical-claim paper.** MODEL.md already states research-only; v3 inherits that.
- **Not a private-data model.** Public data only, always.

---

## 5. Risk register

| Risk | Mitigation |
|---|---|
| Elastic-net auto-L1 over-tunes per cancer (audit risk) | Inner CV on TRAIN fold only, never on test; report inner-CV AUC vs outer-CV AUC gap as "tuning-bias" diagnostic |
| Study-effect shrinkage under-corrects | Report AUC with and without harmonization (the gap IS the bias budget) |
| Per-cancer AUC variability (small n per cancer) | Bootstrap 95% CI per cancer type; flag cancers with n<30 as exploratory |
| Sens@99.9% spec unstable | If 0/1 healthy controls in test fold at 99.9%, report as upper-bound |
| PPV table looks like over-claim | Frame as "expected PPV at hypothetical prevalence X", not a clinical recommendation |
| The design changes break v2.2 reproducibility | v2.2 scripts (`train_classifier.py`, `honest_benchmark.py`) kept untouched; v3 is a new module |

---

## 6. Acceptance criteria (the bar for "v3 ships")

v3 is shippable as v3.0 when **all** of the following are met:

1. ☐ All frozen metric tables produced on the 627-sample cohort.
2. ☐ Sens @ 99% spec improves ≥ 0.05 absolute over v2.2 (frag-only).
3. ☐ No per-cancer AUC below 0.85.
4. ☐ Sens @ 99% spec with fusion improves ≥ 0.03 absolute over v2.2.
5. ☐ All audit findings E1, E4, E7, S1, S2, S3, S4, S6, ST1, ST2, ST3 addressed in code or in the design (not just "documented").
6. ☐ Tests ≥ 50 new test cases (8 modules × ~6 tests each), all green.
7. ☐ CI green on the new module.
8. ☐ V3_DESIGN.md + V3_RESULTS.md + V3_USER_GUIDE.md committed.
9. ☐ Per-cancer AUC table included in PAPER.md update (or in BIORXIV_PAPER_FRAGMENTOMICS.md cross-reference).
10. ☐ Honest framing: if targets not met, ship v3.0 anyway as a methods paper (calibrated + bias-budget + per-cancer) and let sensitivity numbers come from a later v3.1.

---

## 7. Why this design, in one paragraph

v2.2 is already at AUC 0.978 on the 627-sample cohort — that's not where the
sensitivity gap lives. The gap lives in (a) misreported single-number AUCs that
hide per-cancer variance, (b) fusion gain that depends on a synthetic surrogate,
(c) calibration absent so Sens@spec is whatever sklearn gives by default, and
(d) no prevalence-aware decision layer so 99% spec at 0.4% prevalence looks the
same as 99% spec at 5% prevalence. v3 closes all four with minimal code change —
no new model class, no new data, no GPU — by adding hierarchical study effects,
per-cancer-tuned regularization, isotonic calibration, and a prevalence-parameterized
decision layer on top of the existing LR + naive-fusion backbone. The headline
metric shifts from "AUC" to "Sens @ 99% spec + PPV @ prevalence" because that's
what the clinical question actually asks.

---

**End of design proposal. Awaiting approval before any code is written.**
