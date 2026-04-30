# DeepCatch Pipeline Audit — 2026-04-30
## Critical Issues & Fix Recommendations

---

## 🔴 ISSUE 1: 100% Synthetic Data — No Real Validation

**Severity: CRITICAL**

Every single data point in the pipeline is `np.random` generated:

| Generator | What It Does | Issue |
|-----------|-------------|-------|
| `generate_variant_calling_data()` | Random Poisson reads + Beta error rates | No real sequencing data |
| `generate_multimodal_cohort()` | `rng.normal(1.5, 1.5)` for cancer signal | Explicit label leakage |
| `generate_longitudinal_data()` | Rising ctDNA pre-programmed | Perfect separation by design |
| `generate_binary_cancer_data()` | `signal_strength * (g+1)/n` for cancer groups | Signal strength is a parameter you control |
| `generate_ensemble_detector_outputs()` | Pre-set detector means: healthy=0.15 vs cancer=0.55 | AUC is pre-determined |

**Even the "TCGA fallback dataset" is synthetic:**
- Sample IDs: `BRCA_S0000`, `LUAD_S0001` — auto-generated, not real TCGA barcodes
- VAFs: `np.random` uniform-like (0.066–0.447)
- Genes: randomly assigned from a pool

**Consequence:** All reported AUC values (including AUC=1.000) are meaningless. You're measuring how well a classifier can separate data you programmed to be separable.

---

## 🔴 ISSUE 2: Label Leakage in Feature Generation

**Severity: CRITICAL**

In `build_multimodal_features()` (run_tcga_validation.py, line ~288):
```python
signal = 1.5 if y[i] == 1 else 0.0
methylation_features[i] = rng.normal(signal * 0.5, 0.3, n_features//5)
fragment_features[i] = rng.normal(-signal * 0.3, 0.2, n_features//5)
cn_features[i] = rng.normal(signal * 0.4, 0.25, n_features//5)
protein_features[i] = rng.normal(signal * 0.6, 0.35, n_features//5)
```

The label `y` is used **directly** to generate features. Any classifier will trivially achieve AUC ≈ 1.0 because the features are literally `label * constant + noise`.

---

## 🟡 ISSUE 3: Fake ROC Curves

**Severity: HIGH**

In `generate_plots()`:
```python
tpr_vals = fpr_vals ** (0.3 / (1 - auc_mean + 0.01))  # approximate ROC shape
tpr_vals = 1 - (1 - tpr_vals) * (1 - auc_mean * 0.9)
```

The published ROC plots are NOT generated from actual model predictions. They're mathematical approximations drawn from the AUC number. This means a reviewer could detect that the ROC curves don't match the underlying data distribution.

---

## 🟡 ISSUE 4: No Independent Test Set / Temporal Split

**Severity: HIGH**

All experiments use `train_test_split` with `random_state=seed` — a random 60/20/20 split. For a screening test intended for clinical deployment:
- No **temporal validation** (train on earlier samples, test on later)
- No **site-stratified** splits (samples from different hospitals)
- No **independent external cohort** validation
- No **batch-aware** splitting

---

## 🟡 ISSUE 5: Statistical Tests on Synthetic Data

**Severity: MEDIUM**

The p-values (p<0.001) are meaningless because:
- Effect sizes are baked into the data generator
- The number of seeds (5) determines statistical power, not the underlying biology
- Bonferroni correction on 10 pairwise tests from 5 experiments inflates Type II error

---

## 🟢 What IS Working Correctly

1. **Pipeline architecture** — The modular experiment structure, Bootstrap CI, and unified JSON output format are well-designed
2. **Beta-binomial variant caller** — Mathematically sound likelihood-ratio approach
3. **Stacked ensemble** — Proper two-level architecture with StratifiedKFold meta-features
4. **Statistical framework** — DeLong tests, Bonferroni correction are correctly implemented (just applied to wrong data)
5. **CLI and config** — Clean argparse interface, quick mode, seed control

---

## 🔧 Recommended Fixes (Priority Order)

### Fix 1: Real TCGA Data Integration (HIGHEST)
```python
# Use cBioPortal API with proper authentication or download MAF files
# Real TCGA barcodes: TCGA-05-4249-01A-01D-1103-08
# Download: https://portal.gdc.cancer.gov/
# Use TCGAbiolinks (R) or cBioPortal data API
```

**Action:** Download real TCGA MAF from GDC portal → parse with `maftools` → extract mutation calls with real VAFs

### Fix 2: Remove Synthetic Feature Generation for Multimodal
- Replace `build_multimodal_features()` with real methylation (Illumina 450K), RNA-seq, and CNV data from TCGA
- If real multimodal data is unavailable, be transparent that fusion is validated on variant-only features

### Fix 3: Use Real Sequencing Simulation Tools
Replace `np.random.binomial()` with:
- **ART** (Huang 2012) — realistic Illumina read simulator
- **pIRS** — profile-based Illumina read simulator  
- **NEAT** — fine-grained read simulation with error models

### Fix 4: Get Real ROC Curves
Plot actual model predictions, not mathematical approximations:
```python
fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
ax.plot(fpr, tpr)
```

### Fix 5: Independent Validation
- Split by patient (not by position/locus)
- Use temporal or site-stratified splits
- Include external cohort if possible (e.g., ICGC, GENIE)

### Fix 6: Add Calibration Analysis
- Reliability diagrams
- Expected Calibration Error (ECE)
- Brier score per VAF bin

---

## 📊 Benchmark: Current vs Required

| Criterion | Current | Required for Publication |
|-----------|---------|-------------------------|
| Variant data | Synthetic (np.random) | Real TCGA/WGS MAF or ART-simulated |
| Multimodal data | Synthetic with label leakage | Real multi-omics (TCGA methylation/RNA/CNV) |
| ROC curves | Mathematical approximation | Actual model predictions |
| Validation split | Random 60/20/20 | Patient-level temporal/site-stratified |
| Independent cohort | None | At least 1 external dataset |
| Statistical tests | On synthetic data | On real held-out predictions |
| Calibration | Not measured | ECE, Brier score, reliability diagrams |
| Code reproducibility | Good ✓ | Good ✓ |

---

## Bottom Line

The pipeline architecture is solid, but **every result currently reported is an artifact of synthetic data generation with baked-in separation**. Before submission to Bioinformatics/PLOS Comp Bio:

1. Replace synthetic data with real TCGA multi-omics
2. Remove label leakage from feature generation
3. Plot real ROC curves from model predictions
4. Add patient-level stratified validation
5. Validate on at least one independent cohort
