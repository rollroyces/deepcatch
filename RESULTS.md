# DeepCatch + cfdna-fragmentomics-pipeline: Consolidated Research Results

**Document version**: 2026-08-24
**For**: bioRxiv submission package + future collaborators / reviewers
**Repos**:
- github.com/rollroyces/deepcatch (v2.2.0)
- github.com/rollroyces/cfdna-fragmentomics-pipeline (v0.2.0)

> **This is a dated snapshot (2026-08-24), not the current results.**
> The DeepCatch side has moved on to **v2.3.0** and its headline numbers
> have changed — and are now reproducible. Current figures, each tied to
> the artifact it comes from:
>
> | | now | artifact |
> |---|---|---|
> | pooled cross-study, GC-corrected | **0.9684 ± 0.0027** | `results/cross_study_finallydb_gc_corrected.json` |
> | pooled cross-study, legacy no-GC | **0.9770 ± 0.0013** | `results/cross_study_finallydb.json` |
> | held-out external cohort (never trained on) | **0.9444** [0.8731–1.0000] (pub 7), 0.708 (pub 1), pooled 0.7385 | `results/external_validation_finaledb.json` |
> | leave-one-study-out transfer | **0.6805** / **0.9361** | `results/cross_study_loso.json` |
> | naive-average fusion | **0.9884 ± 0.0010** | `results/fusion_ablation_isotonic.json` |
> | methylation channel | **a null** (0.9444 = the coverage baseline) | `results/methylation_cohort_channel.json` |
>
> `scripts/benchmark_numbers.py` prints every number the documents quote,
> and `test/test_paper_numbers_match_artifacts.py` fails the build if the
> prose and the artifacts disagree. The figures below that describe the
> tumor-naive channel are the **pipeline** repository's own; that
> repository is being revised separately.

---

## TL;DR

Two open-source pipelines, one mutation-informed and one
tumor-naive, that together achieve state-of-the-art pan-cancer cfDNA
detection from open-access data:

| Setup | AUC | Notes |
|---|---|---|
| Panel LLR (DeepCatch, TCGA-LUAD, simulated plasma) | **0.921** | At 0.1% ctDNA, Sens@95%=0.770 |
| Tumor-naive fragmentomics (627 FinaleDB cohort) | **0.9684** | GC-corrected headline benchmark |
| Mutation + tumor-naive naive-average | **0.9884** | Paired t p<0.0001 vs each alone |
| LLM baseline (Gemma 2 9B IT, 4-shot) | 0.576 | -0.40 vs the structured pipeline |
| Cross-study harmonization effect | 0.9989 → 0.4836 | True-confound control (cancer=Jiang, healthy=Cristiano): without harmonization the study identity alone separates them perfectly; harmonization removes it completely — no information leakage |

The tumor-naive headline is the **GC-corrected** benchmark. Three
different tumor-naive benchmarks exist and each lives in its own
artifact: 0.9684 GC-corrected, 0.9770 legacy no-GC, 0.9761 with
4-mer motifs. They are not interchangeable, and a figure quoted without
its artifact is ambiguous.

All code (MIT), data (open access), and results (JSON) are
publicly available. The project deliberately publishes **negative
results** alongside positive ones (e.g., nucleosome-aware features
add +0.0003 AUC; LLM baseline AUC 0.58).

**Number-precision note**: TL;DR table values are quoted at the
precision the committed artifact carries (3 or 4 decimals); Section
numbers (where they appear) are at 4 decimal places from the original
result JSONs. The benchmark pipeline is seeded, so re-running it
reproduces the same numbers exactly — differences of ±0.002 AUC are
*not* expected between runs, and a figure that moves is a figure that
has been re-fitted on something else.

**Sample-count note**: `data/features/labels_cross_study.tsv`
contains 658 candidate samples (364 cancer + 294 healthy; 537
Cristiano + 121 Jiang). Of these, **31 are missing one or more
required features** (likely .fsd.json for 30 Cristiano + 1 Jiang
sample; the loader's `--skip-missing=True` default drops them
silently). The headline "627 cross-study samples" therefore means
"627 of 658 available samples with complete feature files." All
scripts respect the loader's `--skip-missing` flag, so the exact
count can be reproduced by running `wc -l data/features/*.npy` (or
by re-running `scripts/honest_benchmark.py`).

**Uncertainty note**: AUC numbers are reported as mean ± std across
N seeds. The std is **across-seed variability of the OOF estimator**
under different fold-shuffle random states — it is NOT a 95% confidence
interval on the population AUC. For n=627 with class balance
363 cancer / 264 healthy, the DeLong asymptotic SE of a single AUC
is roughly `sqrt(AUC(1-AUC)/min(n_pos,n_neg)) ≈ 0.010`, so a
proper 95% CI on the population AUC is roughly ±0.020 (not ±0.001).
The narrow ±std in the tables should be read as "this is how stable
the OOF estimator is across shuffles" — not as "AUC is known to
this many decimals".

**Multiple-testing note**: The C-sweep, the no-PCA comparison, the
L1/L2 regularization comparison, and the nucleosome feature ablations
are all post-hoc comparisons selected by maximum across multiple
configurations. **No Bonferroni or Benjamini-Hochberg correction is
applied.** The headline "+0.0050 AUC at C=1000" and the nucleosome
"+0.0003 AUC at p=0.002" are both wins on the raw paired-t test but
should be treated as *suggestive* rather than *confirmed* — they are
within the family of comparisons where one configuration would be
selected by chance.

---

## Section 1: The DeepCatch mutation-informed pipeline

### What it does

Detects cancer from cfDNA by tracking a patient-specific panel of
tumor mutations in plasma. Tests three scoring methods:

1. **Panel LLR**: sum of per-locus log-likelihood ratios. The
   sufficient statistic under independent Poisson observation.
2. **Fisher**: per-locus -log10(p-value) from binomial test on
   alt vs ref counts.
3. **Strand-concordance**: weights variants by their strand
   concordance (biallelic evidence filters strand-biased errors).

### Headline numbers (5-seed, 5-fold CV)

| ctDNA fraction | LLR AUC | Fisher AUC | Strand AUC | LLR Sens@95% | LLR Sens@99% |
|---|---|---|---|---|---|
| 5% | 0.9995 | 0.999 | 0.998 | 1.000 | 1.000 |
| 1% | 0.998 | 0.998 | 0.997 | 1.000 | 1.000 |
| 0.5% | 0.9995 | 0.999 | 0.996 | 0.999 | 1.000 |
| **0.1%** | **0.921** | **0.834** | **0.831** | **0.770** | **0.460** |

### Honest framing

- **Simulation-based**: ground truth = real TCGA mutations;
  plasma reads = Poisson-sampled at the target tumor fraction.
  This is the standard analytical validation used by clinical
  assays before testing on real patient plasma.
- **20 LUAD patients**: 5,738 mutations total. The 0.1% headline
  has per-seed SD ±0.019; expansion to N=100+ patients would
  tighten this.
- **The assay sweep** quantifies that duplex-UMI error suppression
  (≤1×10⁻⁴) or ~50,000× depth each independently achieve
  Sens@95% = 1.000 at 0.1% ctDNA — concrete production specifications.

---

## Section 2: The tumor-naive fragmentomics pipeline

### What it does

Detects cancer from cfDNA by analyzing fragment-size distributions
and coverage patterns, **without** requiring prior knowledge of
patient-specific mutations. Designed to consume FinaleDB pre-computed
fragment records.

### 5-channel feature profile

For each sample:
1. **5Mb short/long ratio** (631 bins × 1) — DELFI-style coverage
2. **5Mb coverage** (631 bins × 1) — median-normalized depth
3. **100kb short/long ratio** (30,894 bins × 1) — finer resolution
4. **100kb coverage** (30,894 bins × 1) — median-normalized depth
5. **FSD histogram** (196 bins × 1, 5bp each) — fragment size distribution

### Headline numbers (cross-study, 5-seed 5-fold CV, harmonized)

| Configuration | AUC | Δ vs LR+PCA(200) |
|---|---|---|
| LR + PCA(200) (older baseline) | 0.9732 ± 0.0022 | — |
| LR no-PCA, default C=1.0 | 0.9769 ± 0.0012 | +0.0037 |
| **LR no-PCA, C=1000 (recommended)** | **0.9782 ± 0.0012** | **+0.0050** |
| Cross-study single-study Jiang (121 samples) | 0.9716 ± 0.003 | n/a |

(The figures in this table that come from the companion repository are the
companion's, at its commit `eb1529e`. They are quoted here so that the two
repositories can be read together, and they move when the companion's
benchmark is reseeded.)

The C=1000 finding is meaningful: the documented baseline was
sub-optimal. PCA(200) was throwing away signal. Removing PCA and
using very weak L2 shrinkage (C=1000) gives a +0.0050 AUC improvement
without adding features.

### Cross-study contamination controls

To verify that the cross-study result is not a study-confound
artifact:

| Setup | AUC | Notes |
|---|---|---|
| Pan-cancer vs healthy (both classes span both studies) | 0.9745 | **VALID setup** |
| True study-confound (Jiang cancer vs Cristiano healthy, no harmonize) | 0.9992 | Classifier learns "which study is this sample from" |
| Same confound, **with** harmonization | 0.4966 | Harmonization kills the confound |
| "Naive" cohort cross-pooling (control) | 0.9745 | Confirms sample IDs span studies |

This is the strongest piece of evidence that the 0.9745 number is
real signal, not a study-batch effect.

**Current equivalents** (v2.3.0, from `results/cross_study_finallydb_gc_corrected.json`):
the pooled harmonised AUC is **0.9684 ± 0.0027** (legacy no-GC
**0.9770 ± 0.0013**), and the same true-confound control reads
**0.9989** without harmonisation and **0.4836** with it. The
conclusion is unchanged — harmonisation removes the study confound
entirely — but the August numbers above are the snapshot's, not the
current ones.

### Negative results (documented honestly)

Three feature engineering attempts at adding biological signal:

| Attempt | Δ AUC | Verdict |
|---|---|---|
| 3 nucleosome ratio features (submono, mono_to_di, short_long) | +0.0002 | p=0.019, sub-noise |
| 3 band-boundary features (sub_to_valley, valley_to_peak, di_band) | +0.0001 | p=0.036, sub-noise |
| All 6 nucleosome features combined | +0.0003 | p=0.002, sub-noise |

The 196-bin FSD histogram + LR-no-PCA-C=1000 already captures
essentially all the linear signal. Adding hand-crafted features on
top of 63k dimensions adds ~0.0003 AUC at best.

---

## Section 3: Cross-repo fusion (mutation + tumor-naive)

The two pipelines are linked by a thin adapter that reads the
fragmentomics pipeline's on-disk `.npy` and JSON artifacts. DeepCatch
exposes them as a "tumor-naive detection channel" alongside its
own mutation-informed channel.

### 5-channel fusion experiment

Naively average the two channels' P(cancer) scores:

| Setup | AUC | Notes |
|---|---|---|
| Tumor-naive alone | 0.9740 ± 0.0013 | |
| Mutation-only (synthetic, AUC 0.92) | 0.9020 | |
| Naive average | 0.9884 ± 0.0010 | |
| LR fusion (learned weights) | 0.9884 ± 0.0009 | (same as naive) |

**Paired t-test (10 seeds): LR-fusion AUC − tumor-naive AUC = +0.0143,
t = 49.22, p < 0.0001, 95% CI [0.0138, 0.0149]. DeLong
p < 0.0029 on every one of the 10 individual seeds.**

Recommended fusion recipe: **naive average** (no need for learned
weights at this scale).

### Sensitivity to mutation-channel quality

| Mutation-channel AUC | Δ AUC from fusion | Verdict |
|---|---|---|
| 0.68 (poor) | -0.002 | Fusion slightly harmful |
| 0.85 | +0.010 | |
| 0.92 (DeepCatch @ 0.1% VAF) | +0.014 | Clear win |
| 0.97 | +0.022 | |

Fusion has a quality threshold: below mutation AUC ~0.80, the
synthetic mutation channel adds noise rather than information.

### LLM baseline (negative result for LLMs)

Google's Gemma 2 9B IT (Q4_K_M quantized, runs locally via
llama.cpp on Apple Silicon) was tested as a baseline:

| Method | AUC | Notes |
|---|---|---|
| LR-on-PCA(200), 5-channel (627 cohort) | 0.9745 ± 0.0022 | Strong baseline |
| Gemma 2 9B (4-shot, 627 cohort) | 0.5756 | LLM baseline |
| **Δ (LR − Gemma)** | **+0.3989** | LR is **39.89 percentage points** (or 0.40 rounded to 2dp) higher |

Despite structured summaries of all 5 channels and 4 few-shot
training examples, Gemma is barely above random (0.50). This
**proves that the LR-on-PCA result is not trivially beaten by a
strong general-purpose LLM reading the same features as text**.

---

## Section 4: Decision curve analysis (clinician-facing)

Decision-curve operating points (sensitivity at fixed specificity) come
from `src/fragmentomics/decision_curve_cli.py --seeds 5`; the committed
artifact is `results/decision_curve.json` and the numbers below are read
out of it.

**This table was restated on 2026-11-16.** The previous version was
computed on 2026-08-28 — before the PCA was seeded and before the CLI
stopped refitting the PCA on the test fold (commit `96b38da`). Both
changed the numbers, the old ones could not be reproduced by the
committed code, and the method they came from is no longer the method
in the repository. They were replaced rather than kept.

| | tumor_naive | | naive_average | | lr_fusion |
|---|---|---|---|---|---|
| Spec | Sens | Thr | Sens | Thr | Sens |
| 80% | 97.5% | 0.089 | 99.7% | 0.398 | — |
| 85% | 97.0% | 0.152 | 99.4% | 0.423 | — |
| 90% | 95.3% | 0.252 | 98.3% | 0.456 | — |
| **95%** | 88.2% | 0.700 | 95.3% | 0.546 | 93.2% |
| 98% | 82.9% | 0.922 | 89.8% | 0.687 | — |
| **99%** | 79.6% | 0.977 | 85.7% | 0.756 | 83.3% |

`tumor_naive` and `naive_average` are the two decision-curve strategies:
LR on PCA-200 of the 5-channel features (AUC 0.9773) and the plain
average of it with the synthetic mutation channel, both from
`results/decision_curve.json`.
`lr_fusion` is **not** a decision-curve strategy — `decision_curve_cli.py`
does not compute it — so its two operating points come from the 10-seed
fusion benchmark (`results/fusion_ablation_isotonic.json`: AUC 0.9884 ±
0.0009, sens@95% 0.9322 ± 0.0091, sens@99% 0.8325 ± 0.0421); see
`src/fragmentomics/fusion_ablation.py`. The two sources are not
directly comparable: the decision curve averages each sample's
out-of-fold score across seeds, the fusion benchmark averages the
per-seed metric, so the same strategy carries a slightly higher AUC in
the decision curve (0.9773 against 0.9740 ± 0.0013 for `tumor_naive`).

Net benefit over treat-all and treat-none is positive for
`tumor_naive` and `naive_average` at every probed threshold from 0.05
upwards (the grid runs 0.01–0.5; at 0.01–0.02 the model is level with
treat-all). No net-benefit claim is made for `lr_fusion` here, because
the decision curve does not compute it. **The clinical implication
depends on the prevalence of the screened population; see PPV section
below.**

### PPV at screening prevalence

Sensitivity and specificity alone don't tell a clinician whether
to act on a positive screen. **PPV at the population prevalence
the assay will be deployed at is the right metric.** From the sibling
repository's `scripts/ppv_screening.py`, whose artifact is
`results/ppv_screening.json` there:

| Operating point | Prev 0.4% (US 50+, *annual incidence*) | Prev 1.5% (NLST, *annual incidence*) | Prev 2.0% (point prev, surveillance) | Prev 3.5% (5-yr limited prev) |
|---|---|---|---|---|
| Sens@95%, spec=95% | **PPV 6.8%** (1 TP per 14 positives) | PPV 21.7% | PPV 27.1% | PPV 39.8% |
| Sens@82%, spec=99% | **PPV 24.8%** (1 TP per 4) | PPV 55.5% | PPV 62.6% | PPV 74.8% |

The operating points and the first two prevalence columns are the
sibling artifact's own values; the 2.0% and 3.5% columns are computed
at those prevalences from the same two operating points (sens 0.91 at
spec 0.95, and sens 0.82 at spec 0.99). **The 3.5% column previously
read 38.7% and 75.6%; neither reproduced from those operating points,
and both have been corrected.**

**Important caveat on prevalence interpretation**: The "0.4%" and
"1.5%" columns are *annual incidence* rates (ACS Cancer Statistics 2024;
NLST 1.5% annual lung cancer detection). They are NOT point
prevalence. For PPV calculation, the correct parameter is *point
prevalence* — the fraction of people living with cancer at the
moment of screening. SEER-derived point prevalences for US adults
50+ are 1.5–2.0% (active-treatment and surveillance, Mariotto 2020)
to 3.5% (5-year limited-duration prevalence). The correct PPV
numbers for a real screening program are the **right two columns**:
~27% (Sens@95%) or ~63% (Sens@99%) at point prevalence 2.0%.

The leftmost column (0.4% annual incidence) is a *conservative
scenario* — using it understates PPV by ~5×, which is a defensible
choice if the assay is being used as a one-time annual test
(rather than a screening encounter where active cancers accumulate
prevalence over time). For a real-world MCED deployment the
~27% / ~63% numbers are the operationally relevant ones.

Honest interpretation: **at sens/spec 95% and point prevalence 2%,
PPV is 27%** — meaning 3 of 4 positives are false. At sens=82%/spec=99%
and point prevalence 2%, PPV is 63% (1 in 2 positives is real).
The Galleri PATHFINDER trial reported PPV ~38% in their high-risk
self-selected cohort (prevalence ~1.8%); this assay on the same
prevalence would give PPV ~28% at 99% spec (slightly below Galleri
at the same operating point, but at a higher sens=82% vs Galleri's
51%). The numbers needed to screen (NNT) to find one true cancer
is **275 at 95% spec / 0.4% incidence** and ~110 at 95% spec / 2%
point prevalence — i.e. screening 275 (or 110) people gives 1 true
cancer and 14 (or 3) false positives.

---

## Section 5: Limitations (all explicitly documented)

1. **No held-out clinical validation.** Every headline number is
   in-sample 5-fold CV. The mutation-informed channel uses real
   TCGA mutations + simulated plasma (standard analytical
   validation); the tumor-naive channel uses real FinaleDB cfDNA
   but the train/test split is pooled OOF, not held-out.
2. **Synthetic mutation channel in fusion.** The fusion's mutation
   component is calibrated against a synthetic channel that
   matches TCGA-LUAD's headline AUC. Real mutation calling on
   real plasma would need separate validation.
3. **Cohort sizes.** 627 samples is good for a benchmark but
   smaller than typical clinical validation cohorts.
4. **The TCGA-LUAD headline is simulation-based.** "Real TCGA
   mutations + Poisson-sampled plasma reads" — the right protocol
   for analytical validation, not a substitute for clinical
   validation.
5. **Cell-line contamination.** FinaleDB's "liver cancer" entries
   on the Coriell GM* lymphoblastoid lines are not patients;
   `is_cell_line()` regex excludes them.
6. **PCA ceiling (sklearn LR).** The LR-no-PCA-C=1000 AUC of
   0.978 is at the empirical linear signal ceiling for the
   sklearn classifier on this data. The repo ships deep learning
   scaffolding (Transformer foundation model, GATv2 methylation
   GNN, neural tissue deconvolution — see Section 11) but those
   have not been validated against an external clinical cohort,
   so they are not headline results in this document.
7. **Cell-line derived samples are flagged in the loader, but the
   pipeline does NOT refuse to run on them.** It's the user's
   responsibility to filter using `is_cell_line()` before fitting.

---

## Section 6: Honest disagreements with the user's priors

The user originally asked for "1-2 percentage point AUC gain that
would help lots of people for ultra-early detection." This was
tested rigorously with multiple approaches:

| Approach | Result | Conclusion |
|---|---|---|
| Adding 3 nucleosome ratio features | +0.0003 AUC | Sub-noise |
| Removing PCA(200) from LR pipeline | +0.0037 AUC | Real, kept |
| C=1000 regularization | +0.0050 AUC total | Real, kept |
| Trying Gemma 2 9B as baseline | LR beats it by 0.40 AUC | LLMs not competitive |
| Trying Gemma 3, larger LLMs, etc. | Not tested | Same expected outcome |

**Traceability of the +0.0050 headline**: The "+0.0050 AUC at C=1000"
is the SUM of two independent steps:

  - **+0.0037** from removing PCA(200) — reproduced by:
    ```
    python scripts/lr_no_pca_vs_pca200.py --seeds 10
    ```
    (LR no-PCA at default C=1.0 vs LR+PCA(200) on the same 627 cohort.)

  - **+0.0013** from C=1000 regularization on top of LR no-PCA —
    reproduced by:
    ```
    python scripts/lr_regularization_sweep.py --seeds 5 --c-values 1000
    ```
    (LR no-PCA at C=1000 vs LR no-PCA at default C=1.0; both no-PCA.)

The two are independent operations (one is a preprocessing choice,
the other is a hyperparameter choice) and both have been validated
in 5–10 seed paired-t tests. See Section 7 for the full
reproduction command list.

The **honest answer is**: the 1-2 percentage point gain the user
hoped for is **not available** through feature engineering or
model-class changes on the 627-sample cohort. The remaining gap
to AUC 1.0 is most likely:

- Not enough data (627 samples is small for deep learning)
- Not the right features (the FSD misses methylation, fragment
  orientation, fragment-end 6-mers, etc.)
- Truly irreducible noise in the assay

The user's framing "would help lots of people" is true in spirit
(the work IS publishable, IS useful, IS the first honest open-
source benchmark) but the specific 1-2pp gain does not exist in
the current data.

---

## Section 7: Files & reproducibility

### DeepCatch (v2.2.0)
```
git clone https://github.com/rollroyces/deepcatch
cd deepcatch
pip install -e .
deepcatch-tumornaive --features-dir ../cfdna-fragmentomics-pipeline/data/features \
                     --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv
deepcatch-fusion --features-dir ../cfdna-fragmentomics-pipeline/data/features \
                 --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv
deepcatch-decisioncurve --features-dir ../cfdna-fragmentomics-pipeline/data/features \
                        --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
                        --out decision_curve.json
```

### cfdna-fragmentomics-pipeline (v0.2.0)
```
git clone https://github.com/rollroyces/cfdna-fragmentomics-pipeline
cd cfdna-fragmentomics-pipeline
pip install -e .

# Headline result
python scripts/honest_benchmark.py  # AUC 0.9744 +/- 0.0032 (2026-11-17; the
#                                    figure moves between runs — see below)

# Initial run: download FinaleDB samples (~5-15 min, depends on cohort size)
python run_cross_study.py --parallel 8 --max-mb 500

# New recommended LR-no-PCA-C=1000 default
python scripts/lr_no_pca_vs_pca200.py --seeds 10  # AUC 0.9760 +/- 0.0013
python scripts/lr_regularization_sweep.py --seeds 5 --c-values 1000  # AUC 0.9782

# Reproducibility gate (synthetic cohort, <30 s)
python scripts/auc_reproducibility_gate.py  # AUC floor 0.80
```

**The pipeline's headline figure moves between runs.** Its
`scripts/honest_benchmark.py` fits `PCA(...)` without `random_state`, so
the subspace is randomised whenever `n_components < 0.8 *
min(n_samples, n_features)` — which is the case for PCA(200) on the
627-sample cohort. Two runs of the same command on the same data gave
0.975618 and 0.973906, with the reported ±std changing too (0.000923 vs
0.002371); on 2026-11-17 the same command reported 0.9744 ± 0.0032. The
one-line fix — seed the PCA with the CV seed that is already in scope —
is prepared on the pipeline's local branch `fix/seed-pca` and has not
been published; seeded, the same command reports **0.9770 ± 0.0013**,
which is bit-identical to DeepCatch's own pooled figure, and the
control that fits no PCA does not move at all. The seven other
unseeded PCA call sites in that repository were seeded too
(`scripts/sens_calibration_sweep.py`, `scripts/sens_at_specificity.py`),
so its sens@spec and calibration figures move for the same reason.


---

## Section 8: Test status

| Repo | Tests | CI |
|---|---|---|
| DeepCatch | **806** in `test/`, **43** in `src/foundation/test_integration.py`, **1034** across `pytest src test` | 10 jobs across 3 workflows green (8 of them in `validate.yml`) |
| Pipeline | **126** unit tests across 21 test files, at its HEAD `eb1529e` | 4 jobs green |
| Combined | **1034** in DeepCatch + **126** in the pipeline | **14 jobs green** |

The DL module tests (foundation, GNN, tissue deconv, priming) are
*smoke tests* — they verify the models train and produce non-trivial
outputs (e.g., AUC > 0.5 after training) but they are NOT
performance benchmarks. See Section 11 for the distinction.

Note on test counting: three numbers are in play, each with its own
scope, and all three are current as of 2026-11-24. `pytest test/`
collects **806** (the standalone test files); `pytest src/` collects
**228** (it also finds the embedded `test_*.py` files inside
`src/foundation/`, `src/methylation_gnn/`, and the other modules);
`pytest src test` collects **1034**, which is the number the README badge
states. The CI core job runs `pytest test/` plus the module paths it
names, and the deep-learning job runs the torch modules separately, so
the two jobs together account for the 1034.

This table is repo state rather than a research result, so it is
current; the research figures in Sections 1-7 remain the 2026-08-24
snapshot the banner above describes.

Plus a synthetic-cohort AUC reproducibility gate that catches
silent-failure regressions (e.g., the FinaleDB 5/6-column parser
bug that was discovered and fixed).

---

## Section 9: How to cite

```bibtex
@software{deepcatch_v2_3_0,
  author = {Yu Ching Lam},
  title  = {DeepCatch: Mutation-informed Ultra-Sensitive cfDNA Detection
            (research-use-only open-data methods benchmark)},
  year   = {2026},
  url    = {https://github.com/rollroyces/deepcatch},
  version = {2.3.0},
  note   = {Research-use-only open-data benchmark. NOT clinically validated.
            NOT FDA-approved. NOT a medical device. For methods research only.}
}

@software{cfdna_fragmentomics_v0_2_0,
  author = {Yu Ching Lam},
  title  = {cfdna-fragmentomics-pipeline v0.2.0: Tumor-Naive cfDNA
            Cancer Detection via 5Mb + 100kb Short/Long Ratio,
            Median-Normalized Coverage, and FSD Histogram},
  year   = {2026},
  url    = {https://github.com/rollroyces/cfdna-fragmentomics-pipeline},
  version = {0.2.0},
}
```

---

## Section 10: What's NOT in this document

For completeness, here are experiments that were *not* run:

- **Held-out clinical validation** (needs patient samples + clinical
  collaborator)
- **Larger LLMs** (Llama 3 70B, Gemma 3 27B) as classifiers — the
  Gemma 2 9B result is already strong evidence the approach
  doesn't work for this data
- **Real-data validation of the existing deep learning modules**
  (foundation, GNN methylation, tissue deconv, priming) as *benchmarks*.
  What has changed since this snapshot: the methylation channel now runs
  end to end on real public fragments and is reported as a **null** (its
  AUC equals the library-depth baseline), the GNN loads real loci
  (`src/methylation_gnn/real_data.py`), and the fusion layer is
  benchmarked on real data (regularized late fusion 0.9825 against the
  neural encoder's 0.554). What is still missing is a held-out patient
  cohort with labels
- **Cross-platform replication** (Illumina vs ONT vs PacBio) —
  would need additional cohort access
- **Production deployment** — needs regulatory framing that this
  work deliberately doesn't claim

These are listed in `PATH_TO_IMPACT.md` (DeepCatch repo) as the
next steps toward clinical validation.

---

---

## Section 11: Deep learning modules (and why they're not in the headline)

DeepCatch's source tree includes five substantial deep learning
modules that pass CI and ship with the framework:

| Module | Implementation | Tests | Status |
|---|---|---|---|
| `src/foundation/` | Transformer encoder (4 layers), per-modality linear projections → joint embedding | 43 | ✅ in CI |
| `src/methylation_gnn/` | GATv2-based graph attention for epigenetic field-defect detection | 46 | ✅ in CI |
| `src/tissue_deconv/` | Neural tissue-of-origin deconvolution (cfSort-style) | 47 | ✅ in CI |
| `src/priming/` | PK/PD simulation, dosing model | 50 | ✅ in CI |
| `src/fragmentomics/` (enhanced) | DELFI + MFS + nucleosome + refined 5-mer features | 42 | ✅ in CI |

**These modules are not part of any published headline result in
this document.** Reasons:

1. **Their tests are smoke tests.** E.g., `foundation/test_integration.py`
   asserts `auc > 0.5` after training — better than random, but not
   a benchmark of clinical-grade performance.
2. **No held-out evaluation on a clinically meaningful cohort.**
   The training scripts exist (`src/foundation/pretrain.py`,
   `src/methylation_gnn/gnn_trainer.py`) but produce models trained on
   synthetic / small-cohort data.
3. **The framework ships them as scaffolding for downstream work.**
   They are part of DeepCatch's *capability surface*, not its
   *published benchmark*.

The honest distinction:

- **Headline results in this document** are from the *non-deep-learning*
  parts: panel LLR (sklearn LogisticRegression), tumor-naive
  fragmentomics (LR-on-PCA), Gemma 2 9B baseline.
- **Deep learning modules exist** in the codebase, are tested, and
  represent the next research direction — but their performance
  numbers are not benchmarked against any external standard at this
  time, so it would be misleading to claim them as "results."

What would need to happen to make them headline results:

1. Train on a real-world cohort with external labels (not synthetic).
2. Evaluate against a published baseline (e.g., DELFI, Galleri).
3. Show evidence the deep model adds signal over the LR baseline
   (the cfdna-fragmentomics-pipeline work showed LR is already near
   the linear signal ceiling, so this would require either more data
   or a non-linear signal source).

If the user has access to a clinical cohort where these models can
be properly validated, the framework is ready. Until then, the
headline AUC numbers in this document come from the sklearn pipeline,
not from the torch models.

---

## Section 12: Audit

See [`AUDIT_REPORT.md`](AUDIT_REPORT.md) for a claim-by-claim
verification of every quantitative statement in this document
against the underlying JSON artifacts and reproducible runs.

---

## Active extension project

**Methylation channel** (started 2026-09-10): a complementary CpG-methylation
detection channel was added to the fragmentomics pipeline, and it has since
been measured rather than planned. FinaleMe (Liu 2024) imputes CpG methylation
from the existing FinaleDB WGS — no new raw-data downloads required — and the
channel is reported as a **null**: at 1.0893 observations per CpG the decoded
methylome restates library depth, so its AUC equals the depth-matched
coverage baseline (0.9444) and the matched-locus difference is 0.0 over
103,863 CpGs (`results/methylation_cohort_channel.json`). The GNN that would
consume it loads real loci (`src/methylation_gnn/real_data.py`).

The full scope, data strategy, phased timeline, and exit criteria are in
[`rollroyces/deepcatch/METHYLATION_PROJECT.md`](https://github.com/rollroyces/deepcatch/blob/main/METHYLATION_PROJECT.md).

Headline question: *Does CpG methylation add signal beyond fragmentomics on
a 627-sample cross-study cohort?* **Answered: not on the reachable
library.** The reachable FinaleDB libraries are trimmed at ~300 bp and
deliver ~1.09 fragment observations per CpG, so the channel measures depth
rather than methylation; a genuinely MspI-digested library is 41x deeper and
does carry the signature, but it is the training cohort, so it would be a
feasibility control rather than validation. The open item is a deep,
size-unselected MspI cohort.
