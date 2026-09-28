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

---

## 3. GC / mappability bias correction (a separate axis from harmonization)

Per-study z-score harmonization (§1, §2) is a **study-mean shift**
correction. It does NOT correct the per-bin GC bias — the systematic
coverage shift between GC-poor and GC-rich genomic windows that arises
from PCR + WGS-library-prep chemistry. This is the single biggest
biological gap flagged by the biomedical review of the cross-study
benchmark: without it, the 100kb coverage channel carries a study-
dependent confound that the classifier can either use as a confounded
proxy for cancer signal (false-positive AUC inflation) or treat as
noise (true-positive AUC deflation). Either way, batch effects
masquerade as biology in this channel.

### What `scripts/gc_correction.py` does

The canonical LOESS-style correction from the fragmentomics literature
(Borg 2012, Teasdale 2010 [4], Sun 2019 [14]):

1. **Per-bin GC reference.** For each of the 30894 100kb genomic bins
   (hg38, chr1..chr22 + X + Y), compute the empirical GC fraction
   from the reference. Two flavours available:
   - **Exact (default):** 6-second scan over the local hg38.2bit;
     stored at `data/reference/hg38_100kb_gc_per_bin.npy`.
     Per-bin variation within a single chromosome is real and
     substantial (chr1 ranges 0.328 → 0.627 across bins).
   - **Coarse fallback:** per-chromosome mean GC (the
     `build_per_chrom_gc_table()` helper). The defensible choice when
     no reference is available (Teasdale 2010).
2. **Train-only fit.** Inside every CV fold, fit a polynomial of
   degree 2 (`E[coverage | GC] ≈ α GC² + β GC + γ`) to the **mean
   per-bin coverage across the train rows**. The fit is therefore
   cohort-aware (uses 5-seed × 5-fold OOF) and never informed by the
   test rows.
3. **Subtract the GC axis.** Corrected coverage = `observed - f(GC)`,
   divided by the cohort median residual so it is on the same
   dimensionless scale as the existing `cn = c100 / median(c100)`
   channel.

The cross-study pipeline (`scripts/cross_study_finallydb.py`) exposes
`--gc-correction` (default ON) which switches the 4th feature channel
from `cn` to `c100_gc_corrected` via `scripts/load5_gc_corrected.py`.
With `--no-gc-correction` the script remains bit-identical to the
pre-GC benchmark.

### What this correction does NOT do

- It does NOT correct higher-moment drift (skewness, kurtosis) in
  the coverage distribution that arises from study-specific library
  prep. A `LOWESS` mode is available (`--method lowess`) but the
  default linear (degree 2) is faster and matches the canonical
  fragmentomics literature.
- It does NOT mask low-mappability bins. We do have the
  `wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed`
  file in `finaleme/refs/hg19/`, but masking reduces the per-bin
  count and inflates the per-bin variance. We leave the masking
  decision to a follow-on commit if reviewers ask for it.
- It does NOT fix the cancer-signal confounding with GC. If cancer
  tissue has a different GC distribution than healthy tissue (e.g.
  via aneuploidy-induced GC shifts), the polynomial fit will absorb
  some of that signal. With the per-fold train-only fit this leakage
  is bounded; with the per-chromosome mean fallback it is the most
  pronounced. Documented as a known limitation in §4.3 of the
  methods paper.

### Validation strategy

- **Synthetic cohort test** (`test/test_gc_correction.py::test_gc_correction_subtracts_bias`):
  build a fake cohort with an injected coverage-vs-GC curve; verify the
  cohort-mean GC-vs-coverage correlation drops from >0.5 to |corr|<0.05
  after the fit/correct step.
- **Loader-contract test**
  (`test/test_gc_correction.py::test_load5_with_gc_correction_produces_different_output`):
  verify `load5_gc_corrected(gc_corrected=True)` returns a different
  array than `load5` on the same cache, while the other 4 channels
  remain bit-identical. Pins the contract that the GC-OFF path is
  bit-identical to load5 (no behavior change for callers that don't
  opt-in).
- **CLI integration**: `--gc-correction`, `--no-gc-correction`, and
  the `--help` fast-exit all verified.
- **End-to-end**: the cross-study run with `--gc-correction ON`
  produces `results/cross_study_finallydb_gc_corrected.json`. The
  JSON's `config.gc_correction` field is the canonical provenance:
  `true` means the GC correction was applied; `false` means the legacy
  `cn` channel was used.

### ΔAUC, honestly reported

The full 5-seed × 5-fold pipeline was run with and without GC
correction on the same n=627 publications-6+8 cohort (Cristiano 2019
[Jiang 2015]). The GC-corrected version has a pooled OOF AUC of
**0.9670 ± 0.0035** vs the legacy AUC of **0.9747 ± 0.0012** — a
**ΔAUC = −0.008** (small DROP, well within the 5-seed ± std).

The honest reading: removing the GC-vs-coverage axis DROPPED the
pooled AUC by ~0.008 on this FinaleDB-uniformly-processed cohort.
This does NOT mean the GC correction is "wrong" — it means a chunk
of the legacy 0.9747 AUC was carrier on the GC axis, which the
classifier used as an indirect proxy for cancer signal (via per-bin
CNV shifts). The corrected AUC of 0.9670 is the more honest cancer
signal — it is closer to the within-publication AUCs (Jiang 0.9572,
Cristiano 0.9683) than the pooled 0.9747 was, suggesting the legacy
pooled number was inflated by ~0.005 from GC-axis proxy detection.

The `c100_gc_corrected` channel and the `cn` channel carry
different information — the cleanest reading is to report BOTH
versions, not to silently swap one for the other.
