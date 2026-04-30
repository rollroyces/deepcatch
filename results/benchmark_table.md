# DeepCatch vs Published Liquid Biopsy Assays: Benchmark Comparison

**Generated:** 2026-04-30 05:29:05

## ⚠ Important Caveat

**Our DeepCatch models are evaluated via SIMULATION STUDIES**, not prospective clinical trials. Results from simulation studies reflect modeled performance under controlled, idealized conditions and are NOT directly comparable to results from clinical validation studies on real patient cohorts. Clinical studies face real-world challenges (sample heterogeneity, collection variability, pre-analytical factors, diverse patient populations) that simulations cannot fully capture.

## Summary Table

| Study / Model | Type | N | Cancer Types | Detection Method | AUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|---|
| DeepCatch: variant_calling | 💻 SIMULATION STUDY ⚠ | 2000 | Simulated (TCGA-parameterized) | Computational: variant_calling pipeline | 0.758 [0.677–0.839] | 0.550 [0.386–0.714] | 0.779 [0.627–0.932] |
| DeepCatch: multimodal_fusion | 💻 SIMULATION STUDY ⚠ | 300 | Simulated (TCGA-parameterized) | Computational: multimodal_fusion pipeline | 0.978 [0.965–0.991] | 0.911 [0.845–0.977] | 0.905 [0.835–0.975] |
| DeepCatch: longitudinal_cet | 💻 SIMULATION STUDY ⚠ | 100 | Simulated (TCGA-parameterized) | Computational: longitudinal_cet pipeline | 1.000 [1.000–1.000] | 0.900 [0.725–1.075] | 1.000 [1.000–1.000] |
| DeepCatch: temporal_transformer | 💻 SIMULATION STUDY ⚠ | 200 | Simulated (TCGA-parameterized) | Computational: temporal_transformer pipeline | 0.999 [0.998–1.001] | 0.920 [0.834–1.006] | 1.000 [1.000–1.000] |
| DeepCatch: ensemble_integration | 💻 SIMULATION STUDY ⚠ | 400 | Simulated (TCGA-parameterized) | Computational: ensemble_integration pipeline | 1.000 [1.000–1.000] | 0.950 [0.862–1.038] | 1.000 [1.000–1.000] |
| DeepCatch: TCGA variant_caller | 💻 SIMULATION STUDY (TCGA-parameterized) ⚠ | 87 | LUAD | TCGA-parameterized: variant_caller with downs... | 0.806 [0.778–0.835] | 0.863 [0.812–0.913] | 0.578 [0.538–0.618] |
| DeepCatch: TCGA multimodal_fusion | 💻 SIMULATION STUDY (TCGA-parameterized) ⚠ | 87 | LUAD | TCGA-parameterized: multimodal_fusion with do... | 1.000 [1.000–1.000] | 0.920 [0.780–1.060] | 1.000 [1.000–1.000] |
| CancerSEEK (Cohen 2018) | 🔬 CLINICAL STUDY | 1005 | 8 types (ovary, liver, stomach, pancreas... | ctDNA mutations (61 amplicons) + 8 protein bi... | N/R | 0.700 | 0.990 |
| GRAIL/Galleri (CCGA, Liu 2020) | 🔬 CLINICAL STUDY | 6689 | 50+ cancer types | Methylation pattern analysis (targeted bisulf... | N/R | 0.518 | 0.995 |
| DELFI (Cristiano 2019) | 🔬 CLINICAL STUDY | 236 | 7 types | Fragment size + coverage patterns (low-covera... | 0.940 | 0.730 | 0.980 |
| LUNAR (Guardant Health) | 🔬 CLINICAL STUDY | 10000 | Colorectal cancer (CRC) screening | ctDNA mutations + methylation (targeted panel... | N/R | 0.830 | 0.970 |
| DETECT-A (Lennon 2020) | 🔬 CLINICAL STUDY | 9911 | 26 cancer types | ctDNA (CancerSEEK) + PET-CT confirmation | N/R | 0.271 (1st blood) | 0.989 |
| PanSeer (Chen 2020) | 🔬 CLINICAL STUDY | 605 | 5 types (stomach, esophagus, colorectal,... | Methylation markers (595 genomic regions) | N/R | 0.880 | 0.960 |
| cfMeDIP-seq (Shen 2018) | 🔬 CLINICAL STUDY | 358 | 7 types | Cell-free methylated DNA immunoprecipitation ... | 0.950 | 0.800 | 0.960 |
| ctDNA (Phallen 2017) | 🔬 CLINICAL STUDY | 238 | Colorectal, breast, lung, ovarian | 58-gene panel with error-correction | N/R | 0.620 | 1.000 |
| CAPP-Seq (Newman 2014/2016) | 🔬 CLINICAL STUDY | 100 | NSCLC (lung) | Deep sequencing with molecular barcoding | N/R | 0.850 (NSCLC) | 0.960 |
| UroSEEK (Springer 2018) | 🔬 CLINICAL STUDY | 570 | Bladder cancer | Urine DNA testing (TERT promoter + 10 genes) | N/R | 0.830 | 0.930 |

## Detailed Study Information

### DeepCatch: variant_calling

- **Type:** SIMULATION STUDY ⚠
- **Assay:** ML Model: Variant Calling
- **Cohort:** 2000 patients
- **Cancer Types:** Simulated (TCGA-parameterized)
- **Method:** Computational: variant_calling pipeline
- **Reported Metrics:**
  - Auc: 0.758 [0.677–0.839]
  - Sensitivity: 0.550 [0.386–0.714]
  - Specificity: 0.779 [0.627–0.932]
  - F1: 0.068 [0.026–0.111]
  - Auprc: 0.468 [0.300–0.636]
  - Accuracy: 0.777 [0.627–0.927]
- **DOI/Ref:** This study (in preparation)
- **Notes:** SIMULATION STUDY — NOT directly comparable to clinical studies. Results reflect modeled performance under idealized conditions.

### DeepCatch: multimodal_fusion

- **Type:** SIMULATION STUDY ⚠
- **Assay:** ML Model: Multimodal Fusion
- **Cohort:** 300 patients
- **Cancer Types:** Simulated (TCGA-parameterized)
- **Method:** Computational: multimodal_fusion pipeline
- **Reported Metrics:**
  - Auc: 0.978 [0.965–0.991]
  - Sensitivity: 0.911 [0.845–0.977]
  - Specificity: 0.905 [0.835–0.975]
  - F1: 0.858 [0.813–0.903]
  - Auprc: 0.957 [0.931–0.982]
  - Accuracy: 0.907 [0.869–0.944]
- **DOI/Ref:** This study (in preparation)
- **Notes:** SIMULATION STUDY — NOT directly comparable to clinical studies. Results reflect modeled performance under idealized conditions.

### DeepCatch: longitudinal_cet

- **Type:** SIMULATION STUDY ⚠
- **Assay:** ML Model: Longitudinal Cet
- **Cohort:** 100 patients
- **Cancer Types:** Simulated (TCGA-parameterized)
- **Method:** Computational: longitudinal_cet pipeline
- **Reported Metrics:**
  - Auc: 1.000 [1.000–1.000]
  - Sensitivity: 0.900 [0.725–1.075]
  - Specificity: 1.000 [1.000–1.000]
  - F1: 0.933 [0.816–1.050]
  - Auprc: 1.000 [1.000–1.000]
  - Accuracy: 0.970 [0.917–1.023]
- **DOI/Ref:** This study (in preparation)
- **Notes:** SIMULATION STUDY — NOT directly comparable to clinical studies. Results reflect modeled performance under idealized conditions.

### DeepCatch: temporal_transformer

- **Type:** SIMULATION STUDY ⚠
- **Assay:** ML Model: Temporal Transformer
- **Cohort:** 200 patients
- **Cancer Types:** Simulated (TCGA-parameterized)
- **Method:** Computational: temporal_transformer pipeline
- **Reported Metrics:**
  - Auc: 0.999 [0.998–1.001]
  - Sensitivity: 0.920 [0.834–1.006]
  - Specificity: 1.000 [1.000–1.000]
  - F1: 0.956 [0.908–1.003]
  - Auprc: 0.998 [0.995–1.001]
  - Accuracy: 0.980 [0.959–1.001]
- **DOI/Ref:** This study (in preparation)
- **Notes:** SIMULATION STUDY — NOT directly comparable to clinical studies. Results reflect modeled performance under idealized conditions.

### DeepCatch: ensemble_integration

- **Type:** SIMULATION STUDY ⚠
- **Assay:** ML Model: Ensemble Integration
- **Cohort:** 400 patients
- **Cancer Types:** Simulated (TCGA-parameterized)
- **Method:** Computational: ensemble_integration pipeline
- **Reported Metrics:**
  - Auc: 1.000 [1.000–1.000]
  - Sensitivity: 0.950 [0.862–1.038]
  - Specificity: 1.000 [1.000–1.000]
  - F1: 0.971 [0.921–1.022]
  - Auprc: 1.000 [1.000–1.000]
  - Accuracy: 0.990 [0.972–1.008]
- **DOI/Ref:** This study (in preparation)
- **Notes:** SIMULATION STUDY — NOT directly comparable to clinical studies. Results reflect modeled performance under idealized conditions.

### DeepCatch: TCGA variant_caller

- **Type:** SIMULATION STUDY (TCGA-parameterized) ⚠
- **Assay:** TCGA-based Validation: Variant Caller
- **Cohort:** 87 patients
- **Cancer Types:** LUAD
- **Method:** TCGA-parameterized: variant_caller with downsampled VAFs
- **Reported Metrics:**
  - Auc: 0.806 [0.778–0.835]
  - Sensitivity: 0.863 [0.812–0.913]
  - Specificity: 0.578 [0.538–0.618]
  - F1: 0.169 [0.156–0.182]
  - Auprc: 0.236 [0.186–0.287]
  - Accuracy: 0.591 [0.554–0.629]
- **DOI/Ref:** This study (in preparation)
- **Notes:** TCGA real data, downsampled to cfDNA VAF levels. Uses real mutation data with simulated cfDNA characteristics. NOT directly comparable to prospective clinical studies.

### DeepCatch: TCGA multimodal_fusion

- **Type:** SIMULATION STUDY (TCGA-parameterized) ⚠
- **Assay:** TCGA-based Validation: Multimodal Fusion
- **Cohort:** 87 patients
- **Cancer Types:** LUAD
- **Method:** TCGA-parameterized: multimodal_fusion with downsampled VAFs
- **Reported Metrics:**
  - Auc: 1.000 [1.000–1.000]
  - Sensitivity: 0.920 [0.780–1.060]
  - Specificity: 1.000 [1.000–1.000]
  - F1: 0.950 [0.862–1.038]
  - Auprc: 1.000 [1.000–1.000]
  - Accuracy: 0.978 [0.939–1.017]
- **DOI/Ref:** This study (in preparation)
- **Notes:** TCGA real data, downsampled to cfDNA VAF levels. Uses real mutation data with simulated cfDNA characteristics. NOT directly comparable to prospective clinical studies.

### CancerSEEK (Cohen 2018)

- **Type:** CLINICAL STUDY
- **Assay:** Multi-analyte blood test (protein + ctDNA)
- **Cohort:** 1005 patients
- **Cancer Types:** 8 types (ovary, liver, stomach, pancreas, esophagus, colorectal, lung, breast)
- **Method:** ctDNA mutations (61 amplicons) + 8 protein biomarkers
- **Reported Metrics:**
  - Sensitivity: 0.700
  - Specificity: 0.990
- **DOI/Ref:** 10.1126/science.aar3247
- **Notes:** Landmark multi-cancer early detection study. Sensitivity varies strongly by cancer type and stage.

### GRAIL/Galleri (CCGA, Liu 2020)

- **Type:** CLINICAL STUDY
- **Assay:** Targeted methylation (bisulfite sequencing)
- **Cohort:** 6689 patients
- **Cancer Types:** 50+ cancer types
- **Method:** Methylation pattern analysis (targeted bisulfite-seq)
- **Reported Metrics:**
  - Sensitivity: 0.518
  - Specificity: 0.995
  - Ppv: 0.440
- **DOI/Ref:** 10.1016/j.annonc.2020.02.011
- **Notes:** Sub-study of Circulating Cell-free Genome Atlas (CCGA). Tissue of origin accuracy: 93%.

### DELFI (Cristiano 2019)

- **Type:** CLINICAL STUDY
- **Assay:** Genome-wide fragmentation profiling
- **Cohort:** 236 patients
- **Cancer Types:** 7 types
- **Method:** Fragment size + coverage patterns (low-coverage WGS)
- **Reported Metrics:**
  - Auc: 0.940
  - Sensitivity: 0.730
  - Specificity: 0.980
- **DOI/Ref:** 10.1038/s41586-019-1272-6
- **Notes:** Fragmentomics approach; AUC 0.94 for cancer detection. Requires only low-coverage WGS.

### LUNAR (Guardant Health)

- **Type:** CLINICAL STUDY
- **Assay:** ctDNA + epigenomic targeted panel
- **Cohort:** 10000 patients
- **Cancer Types:** Colorectal cancer (CRC) screening
- **Method:** ctDNA mutations + methylation (targeted panel)
- **Reported Metrics:**
  - Sensitivity: 0.830
  - Specificity: 0.970
- **DOI/Ref:** NCT04136002
- **Notes:** ECLIPSE study. CRC screening in average-risk population.

### DETECT-A (Lennon 2020)

- **Type:** CLINICAL STUDY
- **Assay:** Multi-analyte: ctDNA + PET-CT
- **Cohort:** 9911 patients
- **Cancer Types:** 26 cancer types
- **Method:** ctDNA (CancerSEEK) + PET-CT confirmation
- **Reported Metrics:**
  - Sensitivity: 0.271 (1st blood)
  - Specificity: 0.989
  - Ppv: 0.195
- **DOI/Ref:** 10.1126/science.abb9601
- **Notes:** Prospective interventional study in women 65-75. First blood test alone: 27.1% sensitivity.

### PanSeer (Chen 2020)

- **Type:** CLINICAL STUDY
- **Assay:** DNA methylation (targeted bisulfite-seq)
- **Cohort:** 605 patients
- **Cancer Types:** 5 types (stomach, esophagus, colorectal, lung, liver)
- **Method:** Methylation markers (595 genomic regions)
- **Reported Metrics:**
  - Sensitivity: 0.880
  - Specificity: 0.960
- **DOI/Ref:** 10.1038/s41467-020-17316-z
- **Notes:** Detected cancer up to 4 years before conventional diagnosis. Pre-diagnosis plasma samples.

### cfMeDIP-seq (Shen 2018)

- **Type:** CLINICAL STUDY
- **Assay:** Immunoprecipitation-based methylation
- **Cohort:** 358 patients
- **Cancer Types:** 7 types
- **Method:** Cell-free methylated DNA immunoprecipitation + sequencing
- **Reported Metrics:**
  - Auc: 0.950
  - Sensitivity: 0.800
  - Specificity: 0.960
- **DOI/Ref:** 10.1038/s41586-018-0708-8
- **Notes:** Enrichment-based methylation profiling. Low input requirement.

### ctDNA (Phallen 2017)

- **Type:** CLINICAL STUDY
- **Assay:** Targeted error correction sequencing (TEC-Seq)
- **Cohort:** 238 patients
- **Cancer Types:** Colorectal, breast, lung, ovarian
- **Method:** 58-gene panel with error-correction
- **Reported Metrics:**
  - Sensitivity: 0.620
  - Specificity: 1.000
- **DOI/Ref:** 10.1126/scitranslmed.aan2415
- **Notes:** Early demonstration of targeted error-correction for ctDNA detection.

### CAPP-Seq (Newman 2014/2016)

- **Type:** CLINICAL STUDY
- **Assay:** CAncer Personalized Profiling by deep Sequencing
- **Cohort:** 100 patients
- **Cancer Types:** NSCLC (lung)
- **Method:** Deep sequencing with molecular barcoding
- **Reported Metrics:**
  - Sensitivity: 0.850 (NSCLC)
  - Specificity: 0.960
- **DOI/Ref:** 10.1038/nm.3519
- **Notes:** iDES-enhanced version improved sensitivity to 0.01% VAF. Requires tumor-informed panel.

### UroSEEK (Springer 2018)

- **Type:** CLINICAL STUDY
- **Assay:** Urine-based multi-gene assay
- **Cohort:** 570 patients
- **Cancer Types:** Bladder cancer
- **Method:** Urine DNA testing (TERT promoter + 10 genes)
- **Reported Metrics:**
  - Sensitivity: 0.830
  - Specificity: 0.930
- **DOI/Ref:** 10.7554/eLife.32143
- **Notes:** Non-invasive urine-based detection. Combined with cytology improves sensitivity.

---

### Interpretation Guide

1. **AUC values from simulation studies** serve as upper-bound estimates of model capability under controlled conditions.
2. **Sensitivity at fixed specificity** (e.g., 99%) is the most clinically relevant metric for screening.
3. **Direct comparison between simulation and clinical studies is NOT statistically valid.** Simulations do not capture pre-analytical variability, population heterogeneity, or clinical confounders.
4. **Future work:** Prospective clinical validation on real cfDNA samples is required before any clinical claims.
