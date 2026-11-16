# External Validation on Held-Out FinaleDB Publications

> Train: FinaleDB 6 + 8 (the cross-study cohort). Test: publications 1, 7, 9 -- never trained on. Public data, no IRB required.

- Generated: `2026-11-16T03:19:50.694947+00:00`

- Train n=627 (cancer 363) | External n=278 (cancer 264, non-cancer 14)

- Classifier: PCA(n_components=200) + LogisticRegression(max_iter=2000), fitted on the training cohort only

- Training-derived threshold (99% specificity on the training OOF): 0.976602 (train specificity at that threshold: 0.9924)


## Pooled external AUC

**0.4177** (DeLong 95% CI 0.3562-0.4793), n_pos=264, n_neg=14


| Sens@spec | 0.95 | 0.98 | 0.99 | 0.995 | 0.999 |
|---|---:|---:|---:|---:|---:|
| Sensitivity | 0.3674 | 0.3674 | 0.3674 | 0.3674 | 0.3674 |

## Per-publication AUC

| Publication | n | AUC [DeLong 95% CI] |
|---|---:|---|
| 1 | 9 | 0.875 [0.63, 1.0] |
| 7 | 31 | 0.9915 [0.9712, 1.0] |
| 9 | 238 (case-only) | sensitivity at the training threshold = 0.4454 |

## Per-disease AUC (vs the external non-cancer controls)

| Disease | n | AUC [DeLong 95% CI] |
|---|---:|---|
| Bladder cancer | - | skipped (n=1) |
| Breast cancer | 120 | 0.4119 [0.3235, 0.5003] |
| Colorectal cancer | 19 | 0.9474 [0.8572, 1.0] |
| Esophageal cancer | - | skipped (n=1) |
| Head and neck cancer | - | skipped (n=1) |
| Inflammatory bowel disease | - | skipped (n=1) |
| Kidney cancer | - | skipped (n=1) |
| Liver transplant | 13 | 0.5302 [0.305, 0.7554] |
| Lung cancer | - | skipped (n=1) |
| Ovarian cancer | - | skipped (n=1) |
| Prostate cancer | 119 | 0.3193 [0.2362, 0.4024] |

## Honest framing

The test publication was never used for training, feature extraction statistics, or hyper-parameter selection. This is an external validation on public data, not a clinical validation: the non-cancer arm is small and comes from different publications than the cancer arm in places, so the pooled AUC mixes a study shift with a cohort-composition shift. Read the per-publication and per-disease rows before quoting the pooled number.

