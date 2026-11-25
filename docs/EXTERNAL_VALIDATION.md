# External Validation on Held-Out FinaleDB Publications

> Train: FinaleDB 6 + 8 (the cross-study cohort). Test: publications 1, 7, 9 -- never trained on. Public data, no IRB required.

- Generated: `2026-11-25T12:55:28.567738+00:00`

- Train n=627 (cancer 363) | External n=322 (cancer 300, non-cancer 22)

- Classifier: PCA(n_components=200) + LogisticRegression(max_iter=2000), fitted on the training cohort only

- Training-derived threshold (99% specificity on the training OOF): 0.976602 (train specificity at that threshold: 0.9924)


## Pooled external AUC

**0.7097** (DeLong 95% CI 0.604-0.8154), n_pos=300, n_neg=22


**Pooled over the publications carrying both classes**: **0.7385** (DeLong 95% CI 0.6281-0.8488), n_pos=61, n_neg=22


| Sens@spec | 0.95 | 0.98 | 0.99 | 0.995 | 0.999 |
|---|---:|---:|---:|---:|---:|
| Sensitivity | 0.3967 | 0.1267 | 0.1267 | 0.1267 | 0.1267 |

## Per-publication AUC

| Publication | n | AUC [DeLong 95% CI] |
|---|---:|---|
| 1 | 52 | 0.708 [0.5148, 0.9012] |
| 7 | 31 | 0.9444 [0.8731, 1.0] |
| 9 | 239 (case-only) | sensitivity at the training threshold = 0.6569 |

## Per-disease AUC (vs the external non-cancer controls)

| Disease | n | AUC [DeLong 95% CI] |
|---|---:|---|
| Bladder cancer | - | skipped (n=2) |
| Breast cancer | 125 | 0.6938 [0.5823, 0.8054] |
| Colorectal cancer | 20 | 0.9182 [0.8295, 1.0] |
| Esophageal cancer | - | skipped (n=2) |
| Head and neck cancer | - | skipped (n=2) |
| Healthy | - | skipped (n=1) |
| Inflammatory bowel disease | - | skipped (n=4) |
| Kidney cancer | - | skipped (n=2) |
| Liver cancer | - | skipped (n=2) |
| Liver transplant | 13 | 0.5542 [0.3579, 0.7504] |
| Lung cancer | 7 | 0.7403 [0.4389, 1.0] |
| Ovarian cancer | - | skipped (n=2) |
| Pancreatic cancer | 8 | 0.6023 [0.3327, 0.8718] |
| Prostate cancer | 122 | 0.7127 [0.59, 0.8355] |
| Skin cancer | - | skipped (n=2) |
| Systemic lupus erythematosus | - | skipped (n=4) |
| Testicular cancer | - | skipped (n=2) |
| Uterine cancer | - | skipped (n=2) |

## Is the separation the sequencing depth?


| Publication | AUC (depth alone) | AUC (model) | AUC (model, depth removed) |
|---|---:|---:|---:|
| 1 | 0.5271 | 0.708 | 0.6718 |
| 7 | 0.5342 | 0.9444 | 0.7436 |
| 9 | - | - | case-only publication: an AUC against depth is not computable, and neither is the depth-removed reading |

Does sequencing depth separate the classes on its own, and does the model's separation survive it? The score is regressed on the depth and what is left is scored again.


## What stage is this cohort at?


**Late stage, where the public records say so.** FinaleDB's sample
record carries no stage field, so the stage had to be read from the
records that describe each publication (`results/external_cohort_stage_annotation.json`, read 2026-11-23):


| Publication | Stage | Where the statement comes from |
|---|---|---|
| 1 — Snyder 2016 | **Stage IV (all cancer samples in the cohort)** | GEO GSE71378, Sample_extract_protocol_ch1 |
| 7 — Sun 2019 | not resolved here | the paper's supplementary tables (Supplemental Table S1, PMC6396422) and the EGA submission of the same cohort |
| 9 — Adalsteinsson 2017 | **metastatic tumors (the cohort the paper describes)** | the paper's title and abstract (PubMed 29109393) |

The external AUCs (0.9444 on publication 7, 0.708 on publication 1, 0.7385 pooled over the publications carrying both classes) measure LATE-STAGE discrimination between cancer plasma and non-cancer plasma. They are a cross-study transfer result and NOT an early-detection result, and they must not be read as one: the held-out cohort contains no case that the public records place at Stage I.


## Honest framing

The test publication was never used for training, feature extraction statistics, or hyper-parameter selection. This is an external validation on public data, not a clinical validation: the non-cancer arm is small and comes from different publications than the cancer arm in places, so the pooled AUC mixes a study shift with a cohort-composition shift. The cancer arm is advanced-stage where the public records characterise it (see the section above), so these numbers say nothing about early detection. Read the per-publication and per-disease rows before quoting the pooled number.

