# CANCER-SEEK -- Data Science & AI Module

**LEGACY SUPERSEDED DESIGN DOC** -- preserved here for archival context.

The CANCER-SEEK Data Science / ML / RL design document is **historically
informative** but is not the current operating spec for `deepcatch`. The
current operating spec lives in:

- **[MODEL.md](../MODEL.md)** -- up-to-date model card for `deepcatch`.
- **[USAGE.md](../USAGE.md)** -- runbook for the repo.
- **[CANCER_SEEK_ARCHITECTURE.md](../CANCER_SEEK_ARCHITECTURE.md)** -- system
  architecture for the ultra-sensitive screening path.
- **[DATA_SCIENCE_ML.md](../../legacy/DATA_SCIENCE_ML.md)** -- this file (legacy).
- **[METHYLATION_PROJECT.md](../../METHYLATION_PROJECT.md)** -- methylation
  sub-project plan.

**Ported from `rollroyces/cancer-screening` DATA_SCIENCE_ML.md on 2026-10-16
without content edits** so the historical reference is byte-for-byte faithful.

---

## 1. Overview -- Why AI is Essential Here

```
The Problem at 0.001% VAF:
  - Raw signal is indistinguishable from noise by eye
  - 1 mutation among 100,000 molecules
  - Traditional thresholds (1%, 5% VAF) completely miss Stage I

Why AI Wins:
  - Pattern recognition across 100s of genomic positions simultaneously
  - Learns subtle co-mutation signatures invisible to statistics
  - Combines multi-modal signals (DNA + methylation + protein) optimally
  - Adapts to new cancer subtypes without manual rule-writing
  - Reinforcement learning optimizes WHO to screen and WHEN
```

## 2. Machine Learning Subsystem

### 2.1 Feature Engineering Pipeline

GENOMIC_FEATURES:
- variant_allele_frequency_{gene}_{position}
- total_mutations_in_{cancer_type}_pathway
- mutational_burden_score
- driver_mutation_present_{gene}
- variant_phylogic_af_class
- cnv_log2_ratio_{chromosome_segment}
- microsatellite_instability_score

METHYLATION_FEATURES:
- beta_value_{cpg_island}_{gene_promoter}
- methylation_similarity_to_{cancer_type}_reference
- epigenetic_age_acceleration
- tissue_of_origin_methylation_score
- differentially_methylated_region_count
- methylation_entropy_score
- hypomethylation_block_score

PROTEIN_FEATURES:
- ca125, cea, ca19-9, psa, afp, cyfra21-1, nse, proGRP
- protein_ratio_to_age_adjusted_reference
- anomaly_score_vs_known_biomarkers_panel

CLINICAL_FEATURES:
- age, sex, bmi, smoking_status, family_history_cancer_count
- prior_malignancy, chronic_inflammation_markers

TEMPORAL_FEATURES:
- delta_VAF_vs_prior_screen
- velocity_of_mutation_appearance
- time_since_last_negative_screen

### 2.2 Classical ML -- First-Line Triage

Models:
- **XGBoost** (primary) -- handles mixed types, SHAP, robust to imbalance.
- **Random Forest** -- baseline, feature importance.
- **Logistic Regression** (calibrated) -- interpretable, calibrated probs.
- **LightGBM** -- fast, large n.

XGBoost params (representative; full set in `src/ultrasensitive/xgboost_classifier.py`):
- `objective: binary:logistic`
- `eval_metric: aucpr`
- `scale_pos_weight: 50`
- `max_depth: 6, learning_rate: 0.05`
- `subsample: 0.8, colsample_bytree: 0.8`
- `min_child_weight: 3`
- `reg_alpha: 0.1, reg_lambda: 1.0`
- `tree_method: gpu_hist` (with `cpu_hist` fallback)

### 2.3 Training Strategy

Dataset (representative):
- TCGA -- 20,000 tumors
- ICGC -- 25,000 tumors
- Published ctDNA studies -- ~10,000 samples
- Healthy controls -- 20,000
- Synthetic spike-ins -- 5,000
- Total: ~80,000 samples

Class imbalance:
- SMOTE oversampling.
- Focal loss `FL(p) = -alpha * (1-p)^gamma * log(p)` with gamma=2.
- Class weights inversely proportional to frequency.
- Positive:Negative = 1:10 for training.

CV:
- 5-fold stratified CV.
- Cancer type held-out validation (no leakage).
- Age-matched healthy controls per cancer type.

## 3. Deep Learning Subsystem

### 3.1 Deep Methylation Encoder (DME)
- Transformer-based, 6 layers / 8 heads.
- Pre-training: Masked Methylation Prediction (BERT-style).
- Fine-tuning: Cancer classification + tissue of origin + tumor fraction regression.

(Ported to `src/ultrasensitive/methylation_encoder.py`; real-data validation
deferred to Phase 4 of `METHYLATION_PROJECT.md`.)

### 3.2 ctDNA Signal Detector (CSD-Net)
- Multi-head attention over UMI-tagged reads.
- Per-position mutation probability + sample-level cancer probability + driver classification.

### 3.3 Multi-Modal Fusion Network (FusionNet)
- Cross-modal attention between (genomic, methylation, protein, clinical) modalities.
- Modality-importance gating.
- Multi-task heads: cancer, tissue, tumor fraction.

(Ported to `src/ultrasensitive/fusion_net.py`; future routing through
`src/multimodal_fusion/advanced_fusion.py`.)

### 3.4 Variational Autoencoder (beta-VAE)
- Unsupervised, trained on healthy samples.
- Anomaly score = beta * KL + reconstruction error.
- Catches rare / novel cancer subtypes missing from the training set.

(Ported to `src/ultrasensitive/vae_anomaly.py`.)

## 4. Reinforcement Learning Subsystem

### 4.1 Adaptive Screening Policy (ASP-Agent)

```
State (s):  Patient profile -- age, sex, risk factors, last screen results,
            mutations detected, methylation pattern, family history.
Action (a): a1 = No screening (wait N months),
            a2 = Standard cfDNA,
            a3 = Full deep screening (cfDNA + methyl + imaging),
            a4 = Immediate clinical referral.
Reward (r): + early cancer detected (QALY),
            + true negative confirmed,
            - missed cancer (QALY loss),
            - false positive (cost).
```

SAC (Soft Actor-Critic) is the agent choice -- off-policy, maximum entropy,
handles continuous action space.

(Ported to `src/screening_policy/rl_agent.py`. The demo Q-function is linear;
upgrading to a neural-network Q head is Phase 5 of
`src/screening_policy/RL_5_UPGRADE.md`.)

### 4.2 Treatment Response Agent
- Post-diagnosis monitoring and treatment selection.
- State: ctDNA levels over time, mutation profile, methylation response.
- Action: Imaging schedule, biopsy timing, treatment escalation/de-escalation.
- Reward: Progression-free survival, QoL, treatment cost.

(Future phase; outside the scope of the v2.3.0 port.)

## 5. Model Interpretability

### 5.1 SHAP Values
Every prediction ships with per-feature SHAP attribution. (Implemented in
`src/ultrasensitive/xgboost_classifier.py::XGBoostCancerClassifier._generate_explanations`.)

### 5.2 Attention Visualization
- Which CpG sites / genomic positions are driving the prediction?
- Browser track export for methylation model.

### 5.3 Uncertainty Quantification
- MC Dropout over 100 forward passes.
- mean +/- std -- high std => flag for expert review.

## 6. Evaluation Metrics

```
Primary:
  - AUC-ROC >= 0.98
  - AUC-PR   >= 0.85 (imbalanced)
  - Sens at 99% spec >= 75% (Stage I)
  - Sens at 99.5% spec >= 60% (Stage I)
  - LOD confirmed at 0.001% VAF

Clinical:
  - NPV >= 99.9%
  - PPV >= 50% (low prevalence)
  - Detection >= 6 months earlier than clinical dx
  - NNS < 500 for Stage I
```

## 7. Technical Implementation Stack

Data Engineering: Airflow, PostgreSQL+Redis, Hadoop/Spark, MLflow.
Model Training: PyTorch, XGBoost/LightGBM, PyTorch Lightning, A100 GPU, W&B.
Inference: ONNX Runtime, TensorRT, FastAPI, Docker+Kubernetes.
Interpretability: SHAP, Captum, HiPlot.

## 8. Data Sources for Training

| Dataset          | Samples  | Key use                          |
|------------------|----------|----------------------------------|
| TCGA             | 20,000   | Pre-train methylation + genomic  |
| ICGC             | 25,000   | External validation              |
| GDC Data Portal  | Ongoing  | Continuous updates               |
| GRAIL CCGA       | 15,000   | ctDNA-specific training          |
| UK Biobank       | 50,000   | Healthy controls + longitudinal  |
| In-house spike-ins | 5,000   | LOD validation                   |
| **Total**        | **~80,000** |                              |

## 9. Deployment Architecture

```
Clinical API (FastAPI + ONNX)
        |
        v
[Triage XGBoost] <1s
        |
        v
[FusionNet (DL)]
        |
        v
[Report Generator (PDF)]
        |
        v
[Result DB]
```

## 10. Priority Implementation Order

```
Step 1 -- Foundation (Month 1-3): XGBoost baseline, SHAP, first AUC baseline.
Step 2 -- Deep Learning (Month 3-6): DME pretrain+finetune, CSD-Net, FusionNet.
Step 3 -- Ensemble + Calibration (Month 6-8): Isotonic regression calibration.
Step 4 -- RL Screening Policy (Month 8-12): SAC agent training + A/B test.
Step 5 -- Production + Validation (Month 12+): ONNX export + clinical API + prospective validation.
```

---

*AI doesn't replace the oncologist -- it gives them a precision instrument to
find the signal in the noise.*

*Ported from `rollroyces/cancer-screening` DATA_SCIENCE_ML.md on 2026-10-16
without editorial change. The current operating spec is `deepcatch/MODEL.md`.*
