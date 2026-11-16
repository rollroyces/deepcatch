# Leave-One-Study-Out (LOSO) Cross-Publication Benchmark
> **Scope**: open-data cross-publication transfer on FinaleDB publications 6 (jiang) + 8 (cristiano). **NOT** clinical validation. **NOT** external cohort validation. The test publication is fully held out of training.
- Generated: `2026-11-16T03:19:15.331592+00:00`
- Features: `5-channel (5mb_ratio + 5mb_coverage + 100kb_ratio + 100kb_counts + FSD-196)`
- PCA n_components: `200` (capped at min(pca_n, n_train) per direction)
- Classifier: `LogisticRegression(max_iter=2000)`

## 1. Headline: train on one publication, test on the other
| Direction | Harmonization | n_train (ca/hea) | n_test (ca/hea) | Test AUC [DeLong 95% CI] | sens@95% | sens@99% |
|---|---|---|---|---|---|---|
| cristiano -> jiang | self_zscore | 506 (274/232) | 121 (89/32) | **0.9361** [0.8909-0.9813] | 0.7753 | 0.5056 |
| cristiano -> jiang | train_global | 506 (274/232) | 121 (89/32) | **0.928** [0.8806-0.9754] | 0.7416 | 0.5169 |
| jiang -> cristiano | self_zscore | 121 (89/32) | 506 (274/232) | **0.6805** [0.6344-0.7267] | 0.2555 | 0.1825 |
| jiang -> cristiano | train_global | 121 (89/32) | 506 (274/232) | **0.6653** [0.6185-0.7121] | 0.2044 | 0.1533 |

## 2. In-domain references (within-publication CV)
| Publication | n (ca/hea) | 5-seed x 5-fold OOF AUC |
|---|---|---|
| cristiano | 506 (274/232) | 0.969 +/- 0.0021 |
| jiang | 121 (89/32) | 0.9789 +/- 0.0027 |

## 3. Per-cancer one-vs-rest AUC on the held-out publication

**cristiano -> jiang** (train types: BRCA, CRC, LUAD, OTHER_C, OV, PAAD)
| Cancer | n | AUC [DeLong 95% CI] |
|---|---:|---|
| HCC_J | 89 | 0.9361 [0.8909-0.9813] |

**jiang -> cristiano** (train types: HCC_J)
| Cancer | n | AUC [DeLong 95% CI] |
|---|---:|---|
| BRCA | 53 | 0.6051 [0.5219-0.6883] |
| CRC | 27 | 0.6078 [0.4659-0.7496] |
| LUAD | 79 | 0.7301 [0.6642-0.7961] |
| OTHER_C | 27 | 0.4393 [0.2991-0.5796] |
| OV | 28 | 0.6421 [0.5152-0.769] |
| PAAD | 60 | 0.8412 [0.7859-0.8965] |

## 4. PCA-component sweep (primary direction)
| n_components | Test AUC |
|---:|---:|

**jiang -> cristiano** (self_zscore)
| n_components | Test AUC |
|---:|---:|
| 5 | 0.6617 |
| 10 | 0.7074 |
| 20 | 0.7223 |
| 50 | 0.667 |
| 100 | 0.6442 |

## 5. Within-study cross-cancer-type control
Leave one cancer type out of a single publication; train on the remaining types + healthy; test on the held-out type. This separates 'unseen study' transfer from 'unseen cancer type' transfer.
| Study | Cancer | n | AUC [DeLong 95% CI] |
|---|---|---:|---|
| cristiano | BRCA | 53 | 0.9862 [0.9688-1.0] |
| cristiano | CRC | 27 | 0.8922 [0.8053-0.9791] |
| cristiano | LUAD | 79 | 0.9869 [0.9699-1.0] |
| cristiano | OTHER_C | 27 | 0.9563 [0.9096-1.0] |
| cristiano | OV | 28 | 1.0 [1.0-1.0] |
| cristiano | PAAD | 60 | 0.857 [0.784-0.9299] |

## Interpretation
- Train-Jiang / test-Cristiano pooled test AUC = 0.6805 (DeLong 95% CI 0.6344-0.7267), n_test=506.
- Jiang is HCC-only (train types: HCC_J) while Cristiano contains no HCC, so this direction conflates two transfers: unseen STUDY and unseen CANCER TYPE. Read it together with the within-study cross-type control before attributing the result to either.
- The pooled cross-study benchmark in scripts/cross_study_finallydb.py never holds a publication out, so it cannot answer this question; these numbers are the transfer measurement the pooled number is often read as.
- Reverse direction (train-Cristiano / test-Jiang) pooled test AUC = 0.9361 (DeLong 95% CI 0.8909-0.9813).
- self_zscore lets the held-out publication normalize itself using its own feature statistics (unsupervised, no labels). train_global applies the training publication's scaler to the test publication with no test-set statistics; it is reported as a sensitivity check because the two publications were pre-processed differently.

## Reproducibility
- **transfer_numbers**: deterministic: the LOSO directions use a seeded PCA (random_state=0) and a fixed train/test split by publication, so the transfer AUCs reproduce exactly
- **in_domain_reference**: computed by scripts/cross_study_loso.py:in_domain_cv, which mirrors cross_study_finallydb.py:pooled_oof (same folds, same per-study z-score harmonization fit on the train fold) but with a SEEDED PCA. The headline benchmark's PCA was unseeded before v2.3.0, so its own per-publication AUC drifts ~0.003 run-to-run on the 506-sample publication; these references can therefore differ from results/cross_study_finallydb.json's per_cohort values by that much
