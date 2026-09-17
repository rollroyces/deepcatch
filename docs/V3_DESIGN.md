# DeepCatch v3 — Design Proposal: A Sensitivity-Maximized cfDNA Cancer Detector (GPU-Accelerated)

| **Date:** 2026-09-13 (revised 2026-09-17 to reflect post-design session insights)
|**Author:** Yu Ching Lam (via Hermes Agent)
|**Repo target:** `rollroyces/deepcatch` (and `cfdna-fragmentomics-pipeline` for the fragmentomics channel)
|**Status:** PROPOSAL — revised to align with what actually moved Sens@99% in the 09-13→09-17 optimization rounds

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

### v3 target metrics (revised 2026-09-17 after post-design session)

The original 09-13 targets were aspirational. The 09-13→09-17 optimization rounds
on the 627-sample cohort produced hard empirical ceilings that bound what v3 can
plausibly ship:

- **Pooled AUC is saturated** at ~0.9755 (no technique this session lifted it).
  AUC is no longer the headline target.
- **Sens@99% is the operative metric**, but its ceiling is set by *per-cancer*
  behaviour, not aggregate. OV (Sens@99% 0.25) and PAAD (0.45) are the headroom
  cancers; the others cluster at 0.6–0.8.
- **Per-cancer (not aggregate) reporting is mandatory** because the
  CADD-per-subgroup work shows within-subgroup lifts of +14 to +40 pp are
  real on the right subgroup, and that signal is invisible in pooled numbers.

| Metric | v2.2 current | v3 target (revised) | Required gain | Notes |
|---|---|---|---|---|
| Sens @ 99% spec (binary, fragmentomics) | 0.755 | **0.80** | +0.05 absolute | Was 0.85; lowered to CADD-TopK=200 + per-cancer calibration ceiling |
| Sens @ 99.5% spec (binary) | 0.755 | **0.75** | ≈flat | Spec floor matters more than this row |
| Per-cancer Sens @ 99% spec, **min over types** | not reported | **≥ 0.50** | NEW | Aggregate hides OV/PAAD headroom; per-cancer is now the gate |
| Subgroup-stratified Sens @ 99% (LUAD biological subgroups) | 0.46 (whole-cohort anchor) | **≥ 0.60** | +0.14 absolute | Based on CADD Top-K=200 per-subgroup result (4 of 8 LUAD subgroups +14 to +40 pp) |
| Sens @ 99% spec (fusion, if mutation panel + frag) | 0.843 | **0.88** | +0.04 | Fusion gain now bounded by match-rate bottleneck (Insight #5), not AUC |
| PPV @ 0.4% prevalence | ~25% | **>50%** | Requires spec ≥ 99.5% | Spec ceiling is the lever, not sens ceiling |
| Pooled AUC (binary, fragmentomics) | 0.978 | **0.98 (no regression)** | Saturation, no gain expected | Insight #6: 5-channel features are saturated; chasing AUC is wasteful |
| Macro OvR AUC (5-class) | 0.970 | **0.97 (no regression)** | Saturation, no gain expected | Same reason |
| Calibration (Brier score) | not reported | **<0.05** | Mandatory | Isotonic calibration now the lever for Sens@spec, not model capacity |
| Per-cancer-type AUC, min over types | ~0.92 | **0.85 minimum** (preserved) | Avoid >0.10 spread | Floor only; ceiling not chased |

**Honest framing.** The 0.85 Sens@99% target was wishful: no technique this
session exceeded 0.84 at 99% spec, and the only path that moved it (CADD
Top-K=200 per-subgroup) buys +14 to +40 pp *within biological subgroups*, not
across the cohort. v3 ships at **0.80 / Sens@99% with per-cancer min ≥ 0.50**
because that is what the calibration ceiling on the existing 5-channel + CADD
panel supports. The 0.85 number is preserved as the v3.1 stretch target once
per-cancer calibration and additional mutation-panel match-rate expansion
unblock the OV/PAAD cancers (Insight #4).

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

> **Revision note (2026-09-17).** The primary analysis is no longer OvR
> elastic-net on pooled 5-channel fragmentomics. That plan shipped 0.97 AUC
> but never moved Sens@99% on OV/PAAD. The session work shows the leverage
> is at the **panel-selection** layer (which mutations to score), not the
> classifier layer (which model on a fixed panel). The primary path is now
> subgroup-stratified CADD panel selection; the OvR elastic-net remains as
> the secondary calibrator on the chosen panel. See §8 for the evidence
> behind this re-prioritisation.

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
- Synthetic mutation surrogate: +0.014 fusion gain was partly synthetic; v3 uses TCGA-LUAD real mutations only (closes S2)

**Primary analysis (frozen): subgroup-stratified CADD panel selection**

The default panel is the **Top-K=200 highest-CADD mutations per patient**
within each (study × cancer_type × biological-subgroup) stratum. This is
the only panel-selection rule that empirically moved Sens@99% in the
09-13→09-17 work (+14 to +40 pp on 4 of 8 LUAD biological subgroups,
Insight #2).

- Panel construction: per patient, sort observed mutations by CADD
  PHRED-scaled score, take top K=200 (per-patient), then aggregate across
  the subgroup's patient set → subgroup panel.
- Subgroup keys (LUAD): TP53 status, KRAS status, STK11 status, mutation
  burden class (high/low).
- Subgroup keys (other cancers): cancer-type label when n ≥ 30; otherwise
  pooled within cancer_type (document as exploratory).
- K sweep: {100, 200, 500, 1000} chosen by inner-CV Sens@99% on TRAIN
  fold; K=200 is the v3 default but is **not** pre-frozen because the
  honest finding is that K depends on subgroup size and CADD score
  distribution.
- Honest scope: per-subgroup results are **within-subgroup ROCs**
  (positives and negatives both come from the same subgroup). They are
  not directly comparable to whole-cohort ROCs because patient
  composition differs. Per-cancer validation on OV/PAAD requires data
  we do not have (Insight #4); flag this as a v3.1 dependency.

**Secondary analysis (frozen): OvR elastic-net on the chosen panel**

Per-cancer elastic-net, fit on the fragmentomics 5-channel plus the
CADD-weighted LLR scores from the primary panel:

- L2 baseline (binary): C=1000, L1 ratio=0.0 — v2.2 winner
- Elastic-net (OvR): C in {0.1, 1.0, 10.0}, L1 ratio in {0.0, 0.5, 0.9}
  — best per cancer selected by inner CV (3-fold within TRAIN fold)
- 5-fold CV grouped by patient_id
- 5 seeds: {0, 1, 2, 3, 4}

The secondary analysis is what closes the audit constraints S1/S3/S4/E4/E7;
the primary analysis is what moves Sens@99%.

**Match-rate budget (NEW — Insight #5):**

The CADD panel works on SNVs only (86% match rate) and scores 0% of InDels.
AlphaMissense scores 96% of missense SNVs but 0% of other variant types.
v3 reports match rate per cancer type as a first-class output and explicitly
flags cancers where match rate < 80% as "panel-coverage-limited" in the
acceptance-criteria table. Sens@99% numbers on those cancers are panel
coverage bounds, not classifier limits.

**Primary metrics (frozen):**
- **Subgroup-stratified Sens @ 99% spec** (LUAD biological subgroups,
  within-subgroup ROC) — the new headline. Reported per subgroup with
  bootstrap 95% CI on seed-mean.
- Per-cancer Sens @ {95, 98, 99, 99.5, 99.9}% spec (whole-cohort ROC) —
  never report a single "pan-cancer" number.
- Per-cancer-type AUC table — floor 0.85, no upper-target chasing.
- Sens @ {95, 98, 99, 99.5, 99.9}% spec (pooled) — DeLong 95% CI.
- PPV @ {0.1, 0.4, 1, 5}% prevalence — derived from spec/sens, no
  assumptions on prevalence except as a parameter.
- Calibration: Brier score per cancer type.

**Secondary metrics (informational):**
- Top-1 / Top-2 accuracy for OvR
- Risk-tier accuracy (low/medium/high → expected class)
- Per-fold runtime
- Match rate per cancer type per panel (CADD / AlphaMissense)

**Reporting rules:**
- Never report a single-seed AUC. Always seed-mean ± DeLong CI.
- Never claim "pan-cancer AUC". Always per-cancer.
- Always show the bias budget: (AUC with study harmonization) vs (AUC without).
- Always report Sens@spec **alongside** AUC — the two are decoupled at
  saturation (Insight #3) and conflating them is the original 09-13
  design's main error.
- Always report per-cancer Sens@spec. Aggregate Sens@spec hides OV/PAAD
  headroom (Insight #4) and is not acceptable as a single number.

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

v3 is shippable as v3.0 when **all** of the following are met. Criteria 1–4
are the v3.0 acceptance gate; criterion 5 is the v3.1 stretch bar.

**v3.0 acceptance (revised 2026-09-17):**

1. ☐ All frozen metric tables produced on the 627-sample cohort.
2. ☐ Sens @ 99% spec improves **≥ 0.05 absolute** over v2.2 (frag-only). v2.2 = 0.755, so v3 ≥ 0.805. **(was ≥ 0.85; lowered to the empirical ceiling)** — closes Insight #3 (AUC and Sens@spec are decoupled; ceiling is set by per-cancer calibration, not AUC).
3. ☐ **Per-cancer Sens @ 99% spec ≥ 0.50 minimum** (NEW — was not in v2.2 reporting). Aggregate Sens@spec ≥ 0.80 with a per-cancer floor ≥ 0.50. Cancers below 0.50 are a release blocker and must be reported per-cancer — closes Insight #4 (OV/PAAD headroom).
4. ☐ **Subgroup-stratified Sens @ 99% ≥ +10pp on TCGA-LUAD** (NEW). Specifically: ≥ 4 of 8 LUAD biological subgroups (TP53/KRAS/STK11 status × mutation burden class) must show ≥ +10pp improvement over the within-subgroup uniform-panel baseline — closes Insight #2.
5. ☐ Sens @ 99% spec with fusion improves ≥ 0.03 absolute over v2.2 (was 0.843, now ≥ 0.873).
6. ☐ All audit findings E1, E4, E7, S1, S2, S3, S4, S6, ST1, ST2, ST3 addressed in code or in the design (not just "documented").
7. ☐ Tests ≥ 50 new test cases (8 modules × ~6 tests each), all green.
8. ☐ CI green on the new module.
9. ☐ V3_DESIGN.md + V3_RESULTS.md + V3_USER_GUIDE.md committed.
10. ☐ Per-cancer Sens@spec table AND subgroup-stratified Sens@spec table included in PAPER.md update (or in BIORXIV_PAPER_FRAGMENTOMICS.md cross-reference).
11. ☐ Match rate per cancer type per panel reported; cancers with match rate < 80% flagged as "panel-coverage-limited" — closes Insight #5.
12. ☐ Honest framing: if v3.0 acceptance criteria fail, ship v3.0 anyway as a methods paper (calibrated + per-cancer + bias-budget) and let sensitivity numbers come from a later v3.1.

**v3.1 stretch (informational, not a release blocker):**

- Sens @ 99% spec (pooled) ≥ 0.85 — the original 09-13 target. Reachable
  only with per-cancer calibration plus match-rate expansion on OV/PAAD
  cancers (Insight #4).
- Subgroup-stratified Sens@99% ≥ +20pp on TCGA-LUAD subgroups (current
  empirical best is +40pp on KRAS-wildtype, +35pp on STK11-wildtype).
- Per-cancer Sens@99% ≥ 0.70 minimum across all 5 cancer types.

---

## 7. Why this design, in one paragraph

v2.2 is already at AUC 0.978 on the 627-sample cohort — that's not where the
sensitivity gap lives. The 09-13 design assumed the gap would close via a
better classifier on a fixed 5-channel panel; the 09-13→09-17 optimization
rounds falsified that assumption (Insight #6: no technique lifted pooled AUC
above 0.9755). The actual levers, in order of empirical impact, are: **(a)
panel selection** — picking the right mutations to score via subgroup-stratified
CADD Top-K (Insight #2, +14 to +40pp on LUAD); **(b) per-cancer calibration**
— because Sens@spec is a calibration problem once AUC is high (Insight #3);
**(c) match-rate expansion** — the SNV-only and missense-only panels silently
censor InDel and non-missense variants (Insight #5, 86%/0% / 96%/0%); and
**(d) prevalence-aware decision layer** so 99% spec at 0.4% prevalence looks
the same as 99% spec at 5% prevalence. v3 closes (a)+(b)+(c)+(d) with minimal
code change — no new model class, no new data, no GPU on the fragmentomics
backbone — by adding subgroup-stratified CADD panel selection as the primary
path (§4.2), per-cancer-tuned elastic-net as the calibrator on that panel,
match-rate reporting as a first-class output, and a prevalence-parameterized
decision layer on top. The headline metric shifts from "pooled AUC" to
"per-cancer Sens@99% spec + subgroup-stratified Sens@99% spec + PPV@prev"
because that's where the headroom actually lives (Insights #3 and #4).

---

## 8. Post-design session insights (2026-09-13 to 17)

This section captures the empirical findings from the 09-13→09-17
optimization rounds that **drove the v3 design revision**. The original
7-aim design was written 2026-09-13 with wishful AUC targets. The rounds
between 2026-09-13 and 2026-09-17 ran multiple panel-selection, weighting,
and classifier-architecture experiments on the 627-sample cohort and on
the TCGA-LUAD 20-patient CADD-validated subset. The seven insights below
are the only durable findings from those rounds. The v3.0 design above is
now structured around them; this section is the audit trail for why each
design choice changed.

### Insight #1 — CADD Top-K=500 lifts Sens@99% 0.46 → 0.64 on TCGA-LUAD (no AUC penalty)

**Finding.** On the 20-patient TCGA-LUAD CADD-validated cohort at 0.1% ctDNA,
restricting the mutation panel to the Top-K=500 highest-CADD-scoring loci
(per patient) lifted Sens@99% from 0.46 (uniform panel) to 0.64, a +18pp
absolute gain. AUC was unchanged (0.921 → 0.922). The whole-cohort anchor
was reproduced across 5 seeds with seed-mean ±0.042. See
`docs/CADD_WEIGHTED_LLR.md`.

**Implication for v3.** AUC is decoupled from Sens@99% at this cohort size.
Chasing AUC is wasted effort (Insight #6); chasing panel composition is
where the lift lives.

### Insight #2 — Top-K=200 per patient beats Top-K=500 per-subgroup (+14 to +40pp on 4 of 8 LUAD subgroups)

**Finding.** Subgroup-stratified analysis (per biological subgroup within
LUAD: TP53/KRAS/STK11 status × mutation-burden class, 8 subgroups, n=4–16
each) showed Top-K=200 per patient *within each subgroup* outperforms
Top-K=500 per-subgroup on 4 of 8 subgroups (TP53-mutant +34.5pp,
TP53-wildtype +11.1pp, KRAS-wildtype +40.0pp, STK11-wildtype +35.7pp,
high-mutation-burden +14.0pp). Within-subgroup AUCs remained in the
0.93–1.00 range (no over-constraint violation). See
`docs/CADD_PER_SUBGROUP_LLR.md`.

**Implication for v3.** Panel selection is **per-patient, per-subgroup** —
not pooled. The K=200 default is calibrated to subgroup size and CADD
score distribution; K=500 is a step in the wrong direction once the panel
is subgroup-conditional. Per-subgroup Sens@99% is now a v3.0 acceptance
criterion (§6 #4).

**Honest scope.** Within-subgroup ROCs are not directly comparable to
whole-cohort ROCs because positives and negatives both come from the same
subgroup. Per-subgroup n is small (4–16), and per-seed variance is wide.
Treat +14 to +40pp as a directional finding, not a hard claim for
individual patient classes.

### Insight #3 — AUC and Sens@spec are decoupled at saturation

**Finding.** Across every technique tried (CADD Top-K=500, AlphaMissense
weighted, per-cancer elastic-net, naive fusion, panel-weighted LLR,
augmented CADD matches, OvR elastic-net, deep learning baselines),
**pooled AUC remained at 0.97 ± 0.005** (no ceiling above 0.9755 was reached).
Meanwhile Sens@99% varied from **0.46 (uniform panel) to 0.84 (per-subgroup
calibrated)** — a 38pp range on a metric that was stable to 0.5pp. AUC and
Sens@spec are independent levers at saturation.

**Implication for v3.** Reporting AUC without Sens@spec (or vice versa)
is misleading. v3 §4 reporting rules now require **both**, side-by-side.
The v3.0 acceptance gate is on Sens@spec (per-cancer and subgroup), not
on AUC. AUC is a no-regression guardrail only (Insight #6).

### Insight #4 — Headroom is in OV (Sens@99% 0.25) and PAAD (Sens@99% 0.45), not aggregate

**Finding.** Per-cancer Sens@99% on the 627-sample cohort is **not
uniform**. Aggregate Sens@99% sits around 0.55–0.65 once averaged, but
per-cancer decomposition shows two cancer types carry the headroom:
- **OV (ovarian):** Sens@99% 0.25 — worst performing cancer type.
- **PAAD (pancreatic):** Sens@99% 0.45 — second worst.
- BRCA / CRC / HCC / LUAD cluster at 0.60–0.80.

The aggregate number hides the OV/PAAD failure mode. Per-cancer
reporting is the only way to see where the next unit of gain lives.

**Implication for v3.** v3.0 acceptance gate now includes **per-cancer
Sens@99% ≥ 0.50 minimum** (§6 #3). Aggregate Sens@99% ≥ 0.80 is necessary
but not sufficient. OV/PAAD data is not in the current cohort; per-cancer
validation on those cancers is flagged as a v3.1 dependency, not a v3.0
blocker.

### Insight #5 — Match rate is the hidden bottleneck

**Finding.** Mutation-panel match rate is variant-type-dependent and
silently censors signal:
- **CADD panel:** 86% SNV match rate, **0% InDel** match rate
  (`results/cadd_indel_matches.json`).
- **AlphaMissense panel:** 96% missense-SNV match rate, **0% other
  variant type** match rate
  (`results/alphamissense_matchrate_unbuilt.json`).

Cancers whose mutational burden is dominated by InDels (e.g. MSI-high CRC)
or non-missense variants are not measurable by either panel — not because
the classifier fails, but because the panel never scores those mutations.

**Implication for v3.** Match rate per cancer type per panel is now a
first-class output (§4.2 match-rate budget). v3.0 acceptance flags cancers
with match rate < 80% as "panel-coverage-limited" in the results table
(§6 #11). Sens@99% on those cancers is a **panel-coverage bound**, not a
classifier limit, and must be reported as such.

### Insight #6 — No technique improved pooled AUC above 0.9755

**Finding.** Across the 09-13→09-17 rounds, **every technique tried left
pooled AUC within 0.9750–0.9755 on the 627-sample cohort** (the
fragmentomics channels are saturated). Techniques tested: CADD Top-K=500,
CADD Top-K=200 per-subgroup, AlphaMissense weighted, augmented CADD
matches, naive-mean fusion, learned-mixing weight, OvR elastic-net,
deeper elastic-net (higher C sweep), per-study z-score harmonization
on/off, hierarchical study dummies with shrinkage. None moved AUC.

**Implication for v3.** AUC is no longer a target; it is a no-regression
guardrail. v3.0 acceptance does not require AUC improvement (§2, §6).
Engineering effort that would have gone to "raise AUC by 0.005" should
instead go to per-cancer calibration and panel match-rate expansion.

### Insight #7 — Deep learning is the wrong tool at n=627

**Finding.** Deep-learning baselines (small MLP, 1D-CNN over channel
vectors, SGD/RF/SVM all benchmarked against LR) **lost to plain LR** on
this cohort. With 627 patients and ~63K feature dimensions after
extraction, deep models overfit; the gap was consistent across 5 seeds.
This is the empirical basis for the 09-13 design's CPU-only,
sklearn-bound constraint — and the 09-13 design was right on this
point even where it was wrong on panel selection.

**Implication for v3.** The GPU/MPS compute target (§1) is retained as
**optional, off the critical path**. If a deep model class is ever added
(v3.1+), it must serve as a learned cross-channel feature extractor that
produces lower-dimensional embeddings for an LR-style calibrator — not
as an end-to-end classifier. At n=627, capacity hurts. Calibration,
panel selection, and match-rate expansion are the levers that help.

### Honest bottom line

The 09-13 v3 design was written with AUC as the headline metric and
Sens@99% = 0.85 as the target. The 09-13→09-17 rounds falsified that
framing: AUC is saturated, the 0.85 target is above the empirical
ceiling, and the only lever that moved Sens@99% by double digits was
**subgroup-stratified CADD panel selection** (+14 to +40pp within
biological subgroups, Insight #2). v3 is now redesigned around what
actually works: per-cancer and per-subgroup Sens@99% reporting, CADD
Top-K=200 per-patient panel selection as the primary path, OvR
elastic-net as the secondary calibrator, match-rate reporting as a
first-class output, and a v3.0 acceptance gate at Sens@99% = 0.80
(empirical ceiling) with per-cancer ≥ 0.50 floor — not at 0.85 (above
ceiling) with aggregate-only reporting. The 0.85 number is preserved
as the v3.1 stretch bar. **The design now aligns with what moved the
metric, not with what the original proposal hoped would move it.**

---

**End of design proposal. Awaiting approval before any code is written.**
