# User Guide — DeepCatch + cfdna-fragmentomics-pipeline

**For**: anyone who wants to use these tools — bioinformaticians running
benchmarks, clinicians evaluating the assay, students learning fragmentomics,
developers extending the code.

**Two repos, one workflow**:
- **`cfdna-fragmentomics-pipeline`** (the *companion* repo) — data download,
  feature extraction, and the headline tumor-naive fragmentomics
  classifier. **This is where most users will start.**
- **`DeepCatch`** (this repo) — the panel-LLR mutation-informed detector, plus
  cross-repo fusion of mutation + tumor-naive channels. **You will need
  this if you want to do fusion or use the panel LLR demo.**

This document covers: install, data, common workflows, the CLI
reference, output interpretation, and troubleshooting. The
**scientific audit trail** is in [`RESULTS.md`](RESULTS.md) and
[`AUDIT_REPORT_2.md`](AUDIT_REPORT_2.md); the **bioRxiv submission
package** is in [`paper/BIORXIV_SUBMISSION.md`](paper/BIORXIV_SUBMISSION.md);
the **model card** is in [`MODEL.md`](MODEL.md).

> **Recommended workflow**: clone and run the **pipeline** repo first
> (it produces the features). Then clone this repo (`DeepCatch`) to
> run the fusion and the panel LLR demo. The two repos share the
> same data directory.

---

## 30-second TL;DR

If you have 5 minutes, run this:

```bash
# 1. Get the code
git clone https://github.com/rollroyces/cfdna-fragmentomics-pipeline
cd cfdna-fragmentomics-pipeline
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e .

# 2. Smoke test (no data needed; ~30 seconds)
python scripts/auc_reproducibility_gate.py
# Expect: "PASS: AUC 0.9981 >= floor 0.80"

# 3. Get data (this is the time-consuming step: 1-3 hours + 100 GB)
python run_cross_study.py --parallel 8 --max-mb 500

# 4. Run the headline benchmark
python scripts/honest_benchmark.py
# Expect: "Cross-study AUC 0.974 ± 0.002 (5-seed pooled OOF)"
```

If anything fails, see the [Troubleshooting](#troubleshooting) section.

---

## 1. Installation

### 1.1 Python version

- **Required**: Python 3.10 or newer (3.11 recommended; CI uses 3.11).
- **Tested**: 3.11 (CI) and 3.14 (local development venv).
- macOS: Xcode CLI tools (`xcode-select --install`).
- Linux: standard build essentials (`apt install build-essential` on Debian/Ubuntu).
- Windows: not officially supported; use WSL2.

### 1.2 Virtual environment (recommended)

```bash
python3.11 -m venv .venv
source .venv/bin/activate   # macOS / Linux
# .venv\Scripts\activate   # Windows (use WSL2)
pip install --upgrade pip
```

### 1.3 Install the pipeline (FIRST)

```bash
git clone https://github.com/rollroyces/cfdna-fragmentomics-pipeline
cd cfdna-fragmentomics-pipeline
pip install -e .
```

This installs the package and 4 console scripts: `cfdna-fetch`,
`cfdna-classify`, `cfdna-fsd`, `cfdna-delfi`. See [§4 CLI reference](#4-cli-reference).

### 1.4 Install DeepCatch (after the pipeline)

```bash
cd ..  # back to wherever you keep your repos
git clone https://github.com/rollroyces/deepcatch
cd deepcatch
pip install -e .
```

This installs 3 more console scripts: `deepcatch-tumornaive`,
`deepcatch-fusion`, `deepcatch-decisioncurve`. The DeepCatch
`USAGE.md` explains how to use them with the pipeline's feature
directory.

### 1.5 Verify installation

```bash
# 57 + 28 = 85 unit tests across both repos
cd cfdna-fragmentomics-pipeline && env -u PYTHONPATH ./.venv/bin/python -m pytest test/
cd ../deepcatch && env -u PYTHONPATH ./.venv/bin/python -m pytest test/

# Smoke test for the headline pipeline
cd ../cfdna-fragmentomics-pipeline
env -u PYTHONPATH ./.venv/bin/python scripts/auc_reproducibility_gate.py
# Expect: "[auc_gate] synthetic cohort: 80 samples, AUC=0.9981"
# Expect: "[auc_gate] PASS: AUC 0.9981 >= floor 0.8000"
```

---

## 2. Data

### 2.1 What data do I need?

**For the headline 627-sample cross-study benchmark**, you need:
1. **Jiang 2015 (publication 6 in FinaleDB)**: 121 samples (89 HCC + 32 healthy) from low-pass WGS.
2. **Cristiano 2019 (publication 8 in FinaleDB)**: 537 samples (8 cancer types + 274 healthy) from deep WGS.

**For the panel LLR demo (DeepCatch)**: 20 TCGA-LUAD patients with somatic mutations (downloaded from GDC automatically).

**For custom use**: any `.frag.tsv` files (chrom, start, end, mapq, strand) — yours, a colleague's, or FinaleDB's.

### 2.2 Important data-availability warning

**The 627-sample feature set is NOT in the repository** (~600 MB, gitignored).
A fresh clone requires you to re-download from FinaleDB:

```bash
python run_cross_study.py --parallel 8 --max-mb 500
```

This takes **1-3 hours and ~100 GB bandwidth** (the 627 raw `.frag.tsv.bgz`
files from FinaleDB sum to ~100 GB; feature extraction compresses them
to ~600 MB of `.npy` and `.fsd.json` in `data/features/`).

### 2.3 Stage your own data (no FinaleDB)

If you have your own `.frag.tsv` files:

```bash
# 1. Create labels file (tab-separated, no header)
#    format: <sample_id>\t<label: cancer|healthy>\t<study: any_tag>
#    e.g.:
#      patient1   cancer  study_a
#      patient2   healthy study_a
echo -e "patient1\tcancer\tstudy_a\npatient2\thealthy\tstudy_a" > data/features/labels.tsv

# 2. Stage one .frag.tsv per sample (5 columns: chrom, start, end, mapq, strand)
ls data/raw/*.frag.tsv
# patient1.frag.tsv  patient2.frag.tsv

# 3. Extract features
cfdna-fsd --features-dir data/features --samples patient1 patient2
cfdna-delfi --features-dir data/features --samples patient1 patient2

# 4. Classify
cfdna-classify --features-dir data/features --labels data/features/labels.tsv
```

---

## 3. Common workflows

### 3.1 "I just want the headline numbers, fast" (30 seconds, no data)

Run the synthetic AUC gate:

```bash
python scripts/auc_reproducibility_gate.py
```

Output:
```
[auc_gate] synthetic cohort: 80 samples, AUC=0.9981
[auc_gate] PASS: AUC 0.9981 >= floor 0.8000
```

This proves the pipeline is installed and working. AUC 0.9981 is on a
synthetic cohort; real-data AUC is documented in §3.4 below.

### 3.2 "I want to reproduce the headline benchmark" (3-5 minutes after data is staged)

After running `run_cross_study.py` (§2.2), run:

```bash
python scripts/honest_benchmark.py
```

Output (abbreviated):
```
=== A: SINGLE-STUDY (Jiang 2015), 5-channel (PCA n=80) ===
  Cohort: 89 cancer + 32 healthy = 121 total
5-channel (PCA n=80, harmonized)         N=121  AUC 0.9716±0.0032
=== C: CROSS-STUDY (pan-cancer), 5-channel (PCA n=200) ===
  Cohort: 363 cancer + 264 healthy = 627 total
5-channel (PCA n=200, harmonized)        N=627  AUC 0.9749±0.0022
=== E: CROSS-STUDY 8-channel (98-subset) (PCA n=200) ===
  Cohort: 50 cancer + 48 healthy = 98 total
8-channel (PCA n=200, harmonized)        N= 98  AUC 0.8761±0.0090
```

**Interpreting output**:
- **Section A**: per-study single-study benchmark. 0.9716 is the
  per-study baseline before any cross-study pooling.
- **Section C**: cross-study pan-cancer headline. **0.9749 ± 0.0022**
  is the documented tumor-naive baseline AUC.
- **Section D**: naive cross-study as a negative control (HCC vs all
  healthy, cancer class only). AUC should be similar to Section C.
- **Section E**: 8-channel (with motifs + mean-length) on the 98-sample
  subset where those features exist. Lower AUC because the subset is
  much smaller.

### 3.3 "I want to know which LR config to use" (~2 minutes)

```bash
# C-sweep on LR no-PCA (fast, with --skip-l1 to avoid the slow L1 saga sweep)
python scripts/lr_regularization_sweep.py --seeds 5 --c-values 1000 --skip-l1
```

Output:
```
  AUC 0.9782 ± 0.0012  (82.8s)
============================================================
penalty         C        AUC   std
------------------------------------------------------------
l2       1000.000     0.9782 0.0012
============================================================
BEST: l2 C=1000.0  AUC=0.9782
```

**Why C=1000?** With 60k features and 627 samples, the default
sklearn `C=1.0` over-shrinks. The optimum sits in a wide plateau
(300-1500) at AUC ≈ 0.978.

### 3.4 "I want the cross-repo fusion AUC 0.989" (DeepCatch)

```bash
cd ../deepcatch
python src/fragmentomics/fusion_ablation.py \
    --features-dir ../cfdna-fragmentomics-pipeline/data/features \
    --labels ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
    --seeds 10 --pca-n 200 \
    --out /tmp/fusion.json
```

Output (abbreviated):
```
[fusion] loaded X.shape=(627, 63246)
[fusion] class balance: 363 cancer / 264 healthy
[fusion] mutation-only sanity: AUC=0.902 Sens@95=0.595 (target ~0.92 / 0.77)
  tumor_naive    : mean=0.974  std=0.0014
  naive_average  : mean=0.9885 std=0.0009
  lr_fusion      : mean=0.9891 std=0.0009
```

**Interpreting fusion output**:
- `tumor_naive` (~0.974): 5-channel fragmentomics alone.
- `naive_average` (~0.989): equal-weighted average of tumor-naive and
  synthetic mutation.
- `lr_fusion` (~0.989): LR-learned fusion weights.
- `lr_fusion > naive_average` indicates the LR is learning useful
  fusion weights; if they're equal, use naive average (simpler).
- The +0.0143 gain over tumor_naive alone is **partly synthetic**:
  the mutation channel is calibrated to a target AUC 0.92, not real
  variant calling on real plasma.

### 3.5 "I want to compute PPV at screening prevalence" (1 second, no data)

```bash
python scripts/ppv_screening.py --out results/ppv_screening.json
```

Output:
```
... Prev 0.4% (US 50+, annual incidence) | PPV 6.8% ...
... Prev 2.0% (point prev, surveillance) | PPV 27.1% ...
... Prev 3.5% (5-yr limited prev)       | PPV 38.7% ...
```

**Interpreting PPV**:
- **0.4% (annual incidence)** is conservative; understates PPV.
- **2.0-3.5% (point prevalence)** is realistic for US adults 50+.
- For a screening program, use the point-prevalence numbers.
- See [`RESULTS.md` §4.1](RESULTS.md) for the full table and discussion.

---

## 4. CLI reference

### 4.1 Pipeline (`cfdna-fragmentomics-pipeline`)

| Script | Console script | Purpose | Typical use |
|---|---|---|---|
| `run_cross_study.py` | — | Full pipeline: FinaleDB → features → labels | First-time setup |
| `scripts/fetch_finaledb.py` | `cfdna-fetch` | Download `.frag.tsv` from FinaleDB | Add more samples |
| `scripts/extract_fsd.py` | `cfdna-fsd` | Compute fragment size distribution | Re-extract after a fix |
| `scripts/extract_delfi.py` | `cfdna-delfi` | Compute 5Mb + 100kb DELFI ratio/coverage | Re-extract after a fix |
| `scripts/extract_motifs.py` | — | Compute 4-mer end-motif frequencies | Add motif features |
| `scripts/train_classifier.py` | `cfdna-classify` | Generic LR / RF / LightGBM classifier | Custom experiments |
| `scripts/honest_benchmark.py` | — | Full 5-section benchmark (A-E) | Re-run headline |
| `scripts/lr_no_pca_vs_pca200.py` | — | LR no-PCA vs LR+PCA(200) | Reproduce +0.0037 gain |
| `scripts/lr_regularization_sweep.py` | — | LR C-sweep on no-PCA features | Reproduce +0.0013 gain |
| `scripts/model_ablation.py` | — | Compare LR / RF / LightGBM (no --help; runs main) | Model-class ablation |
| `scripts/nuc_ablation.py` | — | 5-channel vs 5+nuc vs 5+band vs 5+all | Nucleosome-feature ablation |
| `scripts/eval_8channel.py` | — | 5-channel vs 8-channel on the 98-subset | Re-verify the 8-channel test |
| `scripts/gemma_baseline.py` | — | Gemma 2 9B vs LR baseline (needs `gemma-2-9b-it-Q4_K_M.gguf`) | Reproduce the LLM baseline |
| `scripts/auc_reproducibility_gate.py` | — | Synthetic-cohort AUC floor (CI gate) | Smoke test in CI |
| `scripts/ppv_screening.py` | — | PPV at 5 screening prevalences × 4 operating points | Clinical framing |
| `scripts/gc_correction.py` | — | GC-bias correction on coverage vectors | Optional preprocessing |
| `scripts/build_gc_reference.py` | — | Build the GC reference track from hg38.2bit | First-time setup |

**All scripts accept `--help`** (except `model_ablation.py`, which
runs its main block on import — a known bug; use `python -c "import sys;
sys.path.insert(0, 'scripts'); import model_ablation; model_ablation.main()"`).

### 4.2 DeepCatch

| Script | Console script | Purpose |
|---|---|---|
| `real_tcga_validation.py` | — | Run the headline panel-LLR benchmark on TCGA-LUAD (downloads from GDC) |
| `scripts/run_jiang_pipeline.py` | — | Run the Jiang 4-mer analysis (needs `data/deepcatch_data.xlsx` from Prof. Jiang's lab) |
| `scripts/smoke_tumor_naive_integration.py` | — | CI smoke test for the tumor-naive adapter + fusion |
| `scripts/adapter_auc_gate.py` | — | Synthetic-cohort AUC floor for the adapter |
| `src/fragmentomics/train_tumor_naive.py` | `deepcatch-tumornaive` | Run the tumor-naive adapter on a features dir |
| `src/fragmentomics/fusion_ablation.py` | `deepcatch-fusion` | Run the cross-repo fusion ablation |
| `src/fragmentomics/decision_curve_cli.py` | `deepcatch-decisioncurve` | Compute decision-curve operating points |

### 4.3 Common flags

Most scripts accept:
- `--features-dir DIR` — directory of `{sample}.delfi_*.npy` + `{sample}.fsd.json` (default: `data/features`)
- `--labels FILE` — TSV file with `<sample_id>\t<label>\t<study>` (default: `data/features/labels_cross_study.tsv`)
- `--seeds N` — number of CV random states (default: 5; the headline used 10)
- `--out FILE` — output JSON path

---

## 5. Interpreting the outputs

### 5.1 AUC conventions

- **AUC** = area under the ROC curve, ∈ [0, 1]; 0.5 = random, 1.0 = perfect.
- All AUC numbers in this project use **5-fold stratified cross-validation, pooled OOF predictions** (each test prediction made on data never seen during training).
- `mean ± std` over seeds: the std is across-seed variability of the OOF estimator, **NOT** a 95% confidence interval. The proper DeLong CI on n=627 is roughly ±0.020.

### 5.2 Sens / Spec conventions

- **Sens@95%** = sensitivity (true positive rate) at specificity ≥ 95%. The threshold is chosen as the largest FPR ≤ 0.05, not by optimizing on test data.
- **Sens@99%** = same at specificity ≥ 99%. With n=264 healthy, FPR ≤ 0.01 means ≤ 2.6 false positives, often quantized to 0 or 1.

### 5.3 PPV / NPV conventions

- **PPV** = P(cancer | test positive). Use **point prevalence** (not annual incidence) for a real screening program.
- US 50+ point prevalence: ~1.5-2.0% (active treatment/surveillance), ~3.5% (5-year limited-duration).
- See [RESULTS.md §4.1](RESULTS.md) for the full PPV table.

### 5.4 When results look wrong

- **AUC ~ 0.5** for all configs: the labels file is wrong (label column is reversed, or "cancer" and "healthy" are swapped). Check with `head -3 data/features/labels_cross_study.tsv`.
- **AUC ~ 0.5 for one cohort, ~1.0 for another**: severe batch effect with no overlap. Harmonization is the answer; check `_harmonize` is being called (it's automatic in `train_classifier.evaluate_cv`).
- **NAs / NaNs in the output**: a feature column is all-NaN. The current code drops columns with std=0; columns with mixed-NaN/valid need the per-fold median imputation (already in `train_classifier.evaluate_cv` as of v0.2.0).
- **Runtime > 10 min for `honest_benchmark.py`**: probably the deep-WGS samples in Cristiano 2019. Re-run with `--pca-n 80` for a faster config.
- **`--help` triggers a 5-minute benchmark**: this is the original `honest_benchmark.py` bug. The current version (post round-3 fix) runs `--help` in 2.8s.

---

## 6. Common use cases

### 6.1 "I want to add a new cancer type to the cohort"

1. Edit `data/features/labels_cross_study.tsv` to add the new samples:
   ```
   <new_sample_id>   cancer   <study>
   <new_sample_id>   healthy  <study>
   ```
2. Drop the corresponding `.frag.tsv` files into `data/raw/`.
3. Run `python run_cross_study.py --parallel 8 --max-mb 500` (it will only process the new samples).
4. Re-run the benchmark: `python scripts/honest_benchmark.py`.

### 6.2 "I want to add a new feature"

1. Find the right hook:
   - For FSD: edit `scripts/extract_fsd.py` (function `fsd_summarize` or `fsd_json`).
   - For DELFI: edit `scripts/extract_delfi.py` (`compute_delfi`).
   - For motifs: edit `scripts/extract_motifs.py` (`compute_motifs`).
2. Add the feature to `_build_X` in the script's caller (e.g., `train_classifier._build_full_profile`).
3. Re-extract features and re-run the benchmark.

### 6.3 "I want to compare my assay to this one" (head-to-head)

1. Train on this data using the documented config (`lr_no_pca_vs_pca200.py` / `lr_regularization_sweep.py`).
2. Train on your data using the **same feature engineering pipeline** (call the `_load` functions in `train_classifier.py`).
3. Report both AUCs with 95% DeLong CIs; the comparison is fair only if feature engineering is identical.

### 6.4 "I want to add DeepCatch's mutation channel to my pipeline"

1. Install DeepCatch: `pip install -e ../deepcatch`.
2. Generate the synthetic mutation score (calibrated to AUC 0.92) using `src/fragmentomics/fusion_ablation.py:_simulate_mutation_scores`.
3. Add fusion logic following the `naive_average` or `lr_fusion` patterns in `fusion_ablation.py`.
4. **Honest caveat**: the synthetic mutation channel is not real variant calling. For real fusion experiments, replace it with an actual Signatera-like or CAPP-Seq-like pipeline.

---

## 7. Troubleshooting

### 7.1 "ModuleNotFoundError: No module named 'fragmentomics'"

The repo's internal modules aren't on the Python path. Either:
- `pip install -e .` (recommended), or
- Set `PYTHONPATH=src` for the deepcatch repo or `PYTHONPATH=scripts` for the pipeline.

### 7.2 "All AUCs are ~0.5"

The labels file is wrong. Check:
```bash
head -5 data/features/labels.tsv
head -5 data/features/labels_cross_study.tsv
```
Verify `<sample_id>\t<label: cancer|healthy>\t<study>` (any study tag).

### 7.3 "FileNotFoundError: data/features/labels_cross_study.tsv"

You need to download data first (§2.2). Or use the synthetic AUC gate
to verify the install without data (§3.1).

### 7.4 "Runtime too long"

- `lr_regularization_sweep.py` defaults to a 35-min L1 saga sweep.
  Add `--skip-l1` to skip it (run takes ~80s instead).
- `honest_benchmark.py` runs 5 sections; each is ~30s. The bottleneck
  is the cross-study section (627 samples × 5 folds × 5 seeds).
  Reduce `--seeds` to 3 for a faster smoke test.
- `gemma_baseline.py` requires a 5.5 GB GGUF model file at
  `~/models/gemma-2-9b-it-Q4_K_M.gguf`. The script also takes
  15-20 min for 627 samples; pass `--limit 40` for a faster smoke test.

### 7.5 "LR + PCA(200) gave 0.9732 but a re-run gave 0.9750"

This is normal LR convergence drift (±0.002). The DeLong asymptotic SE
on a single AUC at n=627 is ~0.010; the proper 95% CI is ~±0.020.
Don't panic if you get a slightly different number on re-run.

### 7.6 "My added feature doesn't change AUC"

The 5-channel baseline is near-optimal for the linear signal in the
196-bin FSD. Adding a few features to 60k dimensions is below the
LR's effective signal-to-noise. Try (a) increasing the number of new
features to ≥50, or (b) using a non-linear model (RF with
`--max-features=sqrt`, LightGBM).

### 7.7 "Galleri PATHFINDER is better (PPV 38%); why is your PPV lower?"

Galleri PATHFINDER used a methylation panel (Liu 2020). The cfDNA
fragmentomics-only approach here is fundamentally different and has
**lower signal** for early-stage cancer. The PPV numbers here are
honest; the Galleri comparison is appropriate only for sensitivity
at fixed specificity, not for PPV directly.

---

## 8. Glossary

- **FSD**: Fragment Size Distribution. The 196-bin histogram of
  cfDNA fragment lengths (5bp bins, 20-1000bp range). Cancer-derived
  cfDNA tends to be shorter (~145bp) and shows a different short-fragment
  enrichment profile.
- **DELFI** (DEep Learning for the FrIgmentome): a fragmentomics
  method that bins the genome into 5Mb or 100kb windows and
  computes the short/long ratio per bin.
- **OOF** (Out-Of-Fold): predictions made on data held out during
  training. Standard practice for honest cross-validation.
- **5-seed pooled OOF**: 5 cross-validation runs, each with a
  different random seed for the train/test split, then all OOF
  predictions are pooled. This reduces split-induced variance.
- **Harmonization**: removing study-specific batch effects so that
  features are comparable across studies. This project uses
  per-study z-score StandardScaler, fit on the train fold only.
- **DeLong test** (DeLong, DeLong, Clarke-Pearson 1988): a non-
  parametric test for the difference between two correlated AUCs.
  Returns a z-statistic and p-value.
- **PPV** (Positive Predictive Value): P(cancer | test positive).
  Critical for screening: a 99%-specificity test has PPV ~25% at
  population prevalence (~25% at 0.4% prev), meaning 3 in 4 positives
  are false alarms.
- **cfDNA** (cell-free DNA): DNA fragments in the bloodstream, ~167bp
  on average (one nucleosome + linker). Tumor-derived cfDNA is the
  basis of "liquid biopsy" cancer detection.
- **FinaleDB**: a public cfDNA WGS database hosted by CCHMC
  (Cincinnati Children's Hospital Medical Center). >2,500 datasets
  with uniform preprocessing. Source: https://pubmed.ncbi.nlm.nih.gov/33258919/
- **HEAD/tumor-naive** (this project): the fragmentomics-only
  classifier that doesn't require patient-specific mutation panels.
  Contrast with the "mutation-informed" panel-LLR detector in
  DeepCatch.
- **MAPQ** (Mapping Quality): the Phred-scaled probability that a
  read is mapped to the correct position. cfDNA analysis typically
  filters reads with MAPQ <30 or <60 depending on the application.

---

## 9. Where to go next

- **For audit trail**: [`RESULTS.md`](RESULTS.md) and [`AUDIT_REPORT_2.md`](AUDIT_REPORT_2.md)
- **For bioRxiv submission**: [`paper/BIORXIV_SUBMISSION.md`](paper/BIORXIV_SUBMISSION.md)
- **For DeepCatch's mutation-informed channel**: see the [DeepCatch repo's USAGE.md](../DeepCatch/USAGE.md)
- **For a real-world MCED framing**: this work is honest methods, not
  clinical validation. For clinical-impact claims, see the limitations
  in [RESULTS.md §5](RESULTS.md) and the round-4 journal-reviewer
  analysis in `~/JOURNAL_REVIEW_REJECTION_ANALYSIS.md`.

---

*Last updated: 2026-08-28* (round-4 audit consolidation)
*See git log for the precise commit.*
