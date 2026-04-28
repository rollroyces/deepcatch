# Code Issues Found in DeepCatch Project

## Critical Issues (Must Fix Before Publication)

### C1: Data Leakage in Contrastive Learner Evaluation
**File:** `agent1-variant-calling/evaluate.py`, lines 104-106
**Problem:** The contrastive learner is trained on the entire dataset including test data:
```python
contrastive_caller.fit(data, labels, epochs=50, batch_size=64)
```
**Impact:** All contrastive learner sensitivity/specificity numbers are uninterpretable and likely massively inflated. The model is evaluated on data it was trained on.
**Fix:** Split data into train/val/test BEFORE training. Train only on train set. Evaluate only on test set.

### C4: No Cross-Validation for CET 100% Sensitivity Claim
**File:** `agent3-longitudinal/run_final.py`, line 84
**Problem:** CET threshold (6.0) is calibrated on the same cohort used for reporting. Single random seed (42). No bootstrap CIs. No sensitivity analysis across seeds.
**Impact:** The headline result (100% sensitivity, 99.95% specificity) is unvalidated. This is the single most important finding in the paper.
**Fix:** Implement k-fold cross-validation. Calibrate threshold on a held-out calibration set. Report results across 10+ random seeds with 95% CIs.

### C3: CET Threshold Optimized on Test Data
**File:** `agent3-longitudinal/run_final.py`, line 84
**Problem:** `threshold=6.0` was chosen to maximize performance on the same cohort. This is circular.
**Impact:** Guarantees inflated performance. The threshold would perform worse on unseen data.
**Fix:** Use a proper train/calibration/test split. Tune threshold only on calibration set.

### C10: Circular PoN Validation (MAML Circular Training)
**File:** `agent1-variant-calling/bayesian_caller.py`, lines 70-108 and `agent6-ensemble/ensemble_core.py`
**Problem 1:** Bayes caller PoN is generated from the same simulation framework being validated.
**Problem 2:** MAML meta-learner trains and tests on same data — no held-out meta-test tasks.
**Impact:** Both components produce artificially good results. MAML "99% 1-shot accuracy" is completely invalid.
**Fix:** Use independent data sources for PoN. Split cancer subtypes for MAML meta-train/meta-test.

## Major Issues (Should Fix Before Publication)

### C2: MAML Tested on Training Data
**File:** `agent6-ensemble/ensemble_core.py`, MAMLMetaLearner class
**Problem:** `X_maml = detector_outputs.reshape(1, n, self.n_detectors)` — single task used for both meta-training and testing.
**Impact:** The "99%+ balanced accuracy even with 1-shot adaptation" claim is based on circular evaluation.
**Fix:** Hold out cancer subtypes for meta-testing.

### C5: GNN Fusion Test Set Too Small
**File:** `agent2-multimodal-fusion/evaluate.py`, line 56
**Problem:** 600 patients split 420/90/90. Test set of only 90 patients (~45 cancer).
**Impact:** 95% CI on AUC is approximately [0.58, 0.80]. The "11.9% improvement" is within noise.
**Fix:** Increase cohort to ≥3000 patients with ≥500 in test set. Or use bootstrap within larger test set.

### C6: Temporal Transformer STABLE Accuracy = 0%
**File:** `agent3-longitudinal/temporal_transformer.py`, `agent3-longitudinal/results/final_results.json`
**Problem:** Per-class accuracy: STABLE=0%, RISING=100%, FALLING=100%. Model classifies ALL healthy patients as FALLING.
**Impact:** The "100% specificity" claim is for the wrong task. Model cannot distinguish healthy from benign.
**Fix:** Redesign classification task. Report clinical sensitivity/specificity, not trajectory-class accuracy.

### C7: Deliberately Degraded Feature in Synthetic Data
**File:** `agent2-multimodal-fusion/synthetic_data.py`, lines 159-160
**Problem:** Feature 14 ("latent-correlated") was explicitly replaced with noise:
```python
latent_correlated = self.rng.normal(0, 0.3)  # Now just noise, not directly correlated
```
**Impact:** One of 16 variant features provides no signal. Biases results downward but represents questionable experimental design.
**Fix:** Either remove the feature entirely or restore its correlation. Document the decision.

### C12: No Confidence Intervals on Main Results
**File:** `agent3-longitudinal/improved_methods.py`, lines 369-400
**Problem:** `benchmark_method()` returns point estimates only. No bootstrap, no CI computation.
**Impact:** 100% sensitivity reported without acknowledging n=1000 → CI [99.6%, 100%]. Misleading presentation.
**Fix:** Add bootstrap CIs to all reported metrics.

### C14: Dependent Measurements in Single-Timepoint Baseline
**File:** `agent3-longitudinal/run_final.py`, lines 52-63
**Problem:** Single-timepoint baseline pools ALL measurements across ALL timepoints, treating them as independent. Each patient contributes 8 measurements.
**Impact:** Inflates effective sample size. CIs are artificially narrow. The 64.3% TPR is likely slightly inflated.
**Fix:** Use only one measurement per patient (e.g., first, last, or random).

## Moderate Issues

### C8: Random GNN Edges, Not Biological
**File:** `agent2-multimodal-fusion/models/gnn_fusion.py`, lines 88-160
**Problem:** `rng.choice(self.n_cpg, n_connections, replace=False)` — edges are RANDOM, not based on genomic coordinates or pathway databases.
**Impact:** The paper's claim that edges "encode biological relationships" is misleading. The GNN is learning from random connectivity.
**Fix:** Use real genomic coordinates for proximity edges. Use pathway databases (KEGG, Reactome) for cross-modality edges.

### C9: Arbitrary CET Bonus Weights
**File:** `agent3-longitudinal/improved_methods.py`, lines 135-157
**Problem:** 
```python
streak_bonus = 0.5 * min(n, 5)  # Arbitrary
trend_bonus = max(0, slope) * 3.0  # Arbitrary
```
No calibration, no ablation, no justification for these weights.
**Impact:** These are essentially free parameters tuned to make CET look good. Could be p-hacking.
**Fix:** Either calibrate on held-out data or remove bonuses and use pure SPRT.

### C13: Unstable Sensitivity-at-Specificity Calculation
**File:** `agent2-multimodal-fusion/train.py`, lines 166-170
**Problem:** `np.argmin(np.abs(specificity - 0.99))` picks the closest point on the ROC curve without interpolation.
**Impact:** With a small test set, this can be unstable. One patient can flip the sensitivity-at-99%-specificity from 0% to 10%.
**Fix:** Use interpolation between ROC points. Or report full ROC curves rather than single operating points.

### C15: Single Random Seed Throughout
**Files:** Multiple (all experiments use seed=42)
**Problem:** Every experiment uses the same random seed. No demonstration of reproducibility.
**Impact:** All results could be due to a lucky seed.
**Fix:** Report mean ± SD across ≥5 random seeds.

### C16: Random Graph Structure Sensitivity
**File:** `agent2-multimodal-fusion/models/gnn_fusion.py`
**Problem:** Graph structure is created once and cached. No analysis of how different random graph initializations affect results.
**Impact:** The GNN's performance may depend on the specific random graph edges.

## Minor Issues

### C11: Separate Patient RNG in Simulation
**File:** `agent3-longitudinal/simulation.py`, benign spike section
**Problem:** `spike_rng = np.random.RandomState(trajectory.patient_id * 1000)` — deterministic given patient ID.
**Impact:** Minor. Documented for reproducibility. Not a bug.

### C17: Hardcoded Hyperparameters
**Files:** Multiple
**Problem:** CET threshold=6.0, BOCD hazard=0.06, shared_signal α=0.3 — chosen without grid search.
**Impact:** Optimal values may differ for different cohorts or conditions.

---

## Summary

| Severity | Count | Must Fix Before Publication? |
|----------|-------|------------------------------|
| CRITICAL | 4 | **YES** — C1, C3, C4, C10 |
| MAJOR | 7 | **YES** — C2, C5, C6, C7, C12, C14 |
| MODERATE | 4 | **Recommended** — C8, C9, C13, C15 |
| MINOR | 2 | **Optional** — C11, C17 |
