# External Validation on Held-Out FinaleDB Publications

> Train: FinaleDB 6 + 8 (the cross-study cohort). Test: publications 1, 7, 9 -- never trained on. Public data, no IRB required.

- Generated: `2026-11-15T17:51:13.438394+00:00`

- Train n=627 (cancer 363) | External n=77 (cancer 63, non-cancer 14)

- Classifier: PCA(n_components=200) + LogisticRegression(max_iter=2000), fitted on the training cohort only


## Pooled external AUC

**0.6304** (DeLong 95% CI 0.514-0.7468), n_pos=63, n_neg=14


| Sens@spec | 0.95 | 0.98 | 0.99 | 0.995 | 0.999 |
|---|---:|---:|---:|---:|---:|
| Sensitivity | 0.5556 | 0.5556 | 0.5556 | 0.5556 | 0.5556 |

## Per-publication AUC

| Publication | n | AUC [DeLong 95% CI] |
|---|---:|---|
| 1 | 9 | 1.0 [1.0, 1.0] |
| 7 | 30 | 0.9864 [0.9564, 1.0] |
| 9 | 38 (case-only) | sensitivity at the training threshold = 0.4211 |

## Per-disease AUC (vs the external non-cancer controls)

| Disease | n | AUC [DeLong 95% CI] |
|---|---:|---|
| Bladder cancer | - | skipped (n=1) |
| Breast cancer | 20 | 0.4179 [0.2086, 0.6271] |
| Colorectal cancer | 18 | 0.9405 [0.8439, 1.0] |
| Esophageal cancer | - | skipped (n=1) |
| Head and neck cancer | - | skipped (n=1) |
| Inflammatory bowel disease | - | skipped (n=1) |
| Kidney cancer | - | skipped (n=1) |
| Liver transplant | 13 | 0.5302 [0.305, 0.7554] |
| Lung cancer | - | skipped (n=1) |
| Ovarian cancer | - | skipped (n=1) |
| Prostate cancer | 19 | 0.5075 [0.2856, 0.7294] |

## Honest framing

The test publication was never used for training, feature extraction statistics, or hyper-parameter selection. This is an external validation on public data, not a clinical validation: the non-cancer arm is small and comes from different publications than the cancer arm in places, so the pooled AUC mixes a study shift with a cohort-composition shift. Read the per-publication and per-disease rows before quoting the pooled number.

