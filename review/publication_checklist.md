# DeepCatch Publication Readiness Checklist

## CRITICAL — Must Complete Before Any Submission

- [ ] **Fix data leakage in contrastive learner** (C1): Re-run with proper train/test split. All current contrastive numbers are invalid.
- [ ] **Cross-validate CET results** (C4): Run with ≥5 random seeds. Report mean ± 95% CI. Calibrate threshold on held-out data.
- [ ] **Fix CET threshold optimization** (C3): Use proper train/calibration/test split.
- [ ] **Fix MAML circular evaluation** (C10): Hold out cancer subtypes for meta-testing. Remove the 99% 1-shot claim until validated.
- [ ] **Add confidence intervals** to ALL reported metrics: sensitivity, specificity, AUC, detection time.
- [ ] **Address the 0/30 @ 0.01 VAF anomaly**: Why is contrastive performance worse at HIGHER VAF? This is physically impossible and suggests model instability.

## MAJOR — Should Complete Before Submission

- [ ] **Increase fusion test set**: n ≥ 300 for reliable AUC estimation. Current n=90 gives CIs wider than the claimed effect.
- [ ] **Fix temporal transformer reporting**: Report CLINICAL sensitivity/specificity, not trajectory-class accuracy. Document that STABLE=0% accuracy.
- [ ] **Document the degraded feature** (C7): Explain why feature 14 was deliberately replaced with noise.
- [ ] **Fix dependent measurements** in single-timepoint baseline (C14): Use one measurement per patient.
- [ ] **Add sensitivity analysis** for key biological parameters: ctDNA shedding rate, doubling time, sequencing depth, CHIP prevalence.
- [ ] **Tone down all claims**: Add qualifying language. "100% sensitivity" → "perfect detection in a single simulation run." "1.55× improvement" → more precise definition.
- [ ] **Model CHIP**: At minimum, analyze the expected impact of CHIP on specificity with analytical bounds.

## MODERATE — Improve Before Submission

- [ ] Replace random GNN edges with real biological connections (or stop claiming biological structure)
- [ ] Calibrate CET bonus weights on held-out data or remove them
- [ ] Use interpolation for sensitivity-at-specificity calculations
- [ ] Run ALL experiments with multiple random seeds and report variance
- [ ] Analyze GNN graph structure sensitivity

## CLAIMS TO REMOVE OR DOWNGRADE

- [ ] **DELETE OR COMPLETELY REWRITE**: "5.3× Grail Stage I" (C5) — compares simulation to clinical trial
- [ ] **DELETE OR COMPLETELY REWRITE**: "MAML 99%+ 1-shot accuracy" (C10) — circular evaluation
- [ ] **DELETE OR COMPLETELY REWRITE**: "17-34× earlier detection" (C6) — unvalidated extrapolation
- [ ] **QUALIFY HEAVILY**: "100% sensitivity, 99.95% specificity" (C4) — add CIs, cross-validation, qualifying language
- [ ] **QUALIFY**: "92% sensitivity, 98.5% specificity" (C8) — note this is projected, not measured
- [ ] **QUALIFY**: "3mm³ tumor detection" (C7) — depends on unvalidated VAF↔volume model
- [ ] **QUALIFY**: "6-18 month lead time" (C11) — model-based estimate, not empirically validated
- [ ] **QUALIFY**: "40% cost reduction" (C9) — simplified cost model

## PAPER REVISIONS NEEDED

- [ ] **Abstract**: Add "in simulation" qualifiers. Report ranges instead of point estimates.
- [ ] **Results Section 2.1**: Fix contrastive learner numbers with proper train/test split. Add CIs.
- [ ] **Results Section 2.2**: Note that GNN improvement may not be statistically significant at current sample size.
- [ ] **Results Section 2.3**: Add cross-validation results. Report bootstrap CIs. Calibrate threshold properly.
- [ ] **Results Section 2.4**: Remove MAML 99% claim. Add sensitivity analysis for ρ.
- [ ] **Discussion**: Strengthen limitations section. Add CHIP discussion. Be more honest about simulation-only validation.
- [ ] **Add new section**: "Robustness Analysis" showing performance degradation under parameter perturbations.
- [ ] **Methods**: Document all random seeds used. Add cross-validation methodology.

## ADDITIONAL EXPERIMENTS NEEDED

- [ ] CET with 5+ random seeds and proper train/calibration/test split
- [ ] Contrastive learner with proper train/test split across 10 random seeds
- [ ] GNN fusion with ≥500 test patients across 5 random seeds
- [ ] CET performance with CHIP modeled (even analytically)
- [ ] Sensitivity analysis sweeping: ctDNA shedding (0.3-3×), doubling time (100-500d), depth (10K-100K×), blood volume (5-10mL)
- [ ] Ablation study: CET with/without streak bonus, trend bonus
- [ ] Calibration analysis: how well do CET scores calibrate to actual cancer probability?

## TARGET VENUE REASSESSMENT

After fixes:

| Journal | Readiness After Fixes |
|---------|----------------------|
| Nature Medicine | Still not ready — requires real patient data or at minimum external validation cohort |
| Lancet Oncology | Still not ready — same requirements |
| Cancer Discovery | Borderline — could be submitted if all CRITICAL issues fixed |
| Clinical Cancer Research | Borderline — methods paper with simulation might work |
| **Bioinformatics** | **Ready** — methods journal; simulation + proper validation acceptable |
| **PLOS Computational Biology** | **Ready** — would require fixing all CRITICAL + most MAJOR issues |
| **BMC Bioinformatics** | **Ready** — lower bar for methods papers |
| **medRxiv (preprint)** | **Ready now** — as preprint to solicit feedback; just fix data leakage first |

## ESTIMATED WORK REQUIRED

| Phase | Tasks | Estimated Time |
|-------|-------|---------------|
| Critical fixes | Train/test splits, cross-validation, remove invalid claims | 2-3 weeks |
| Major improvements | Larger test sets, CIs, sensitivity analysis, CHIP modeling | 3-4 weeks |
| Paper revisions | Rewrite claims, add qualifiers, strengthen limitations | 1-2 weeks |
| New experiments | Parameter sweeps, ablation studies, calibration analysis | 2-3 weeks |
| **TOTAL** | | **8-12 weeks** |

## FINAL CHECKLIST ITEM

- [ ] **Have all claims been independently verified by someone not involved in the original development?** ← THIS REVIEW
