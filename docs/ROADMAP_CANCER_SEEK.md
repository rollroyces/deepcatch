# CANCER-SEEK v1.0 -- Project Roadmap (port to DeepCatch umbrella)

**Sourced from** `rollroyces/cancer-screening` PROJECT_PLAN.md (2026-04-27)
and migrated to `rollroyces/deepcatch` as v2.3.0 on 2026-10-16.

## Vision
> Detect cancer at Stage I -- before symptoms, before spread -- when intervention is most effective.
> Our detection floor: **0.001% variant allele frequency** (1 in 100,000 molecules).

## Mission Goals

1. **Sensitivity**: Detect >=80% of Stage I tumors at 0.001% VAF.
2. **Specificity**: <=0.5% false positive rate.
3. **Speed**: Whole pipeline from blood -> report in <72 hours.
4. **Scalability**: Support 10,000+ screens/year.
5. **Accessibility**: Make early screening affordable and accessible globally.

## Workstreams

### Workstream 1 -- Sample Processing & Molecular Biology
- [ ] Blood collection protocol (10-20 mL Streck tube).
- [ ] Plasma cfDNA extraction (QIAamp Circulating Nucleic Acid Kit).
- [ ] UMI library preparation (NEBnext / Twist Duplex Sequencing kit).
- [ ] Methylation conversion (TAPS or bisulfite).
- [ ] QC metrics: DV200 >30%, concentration >=1.7 ng/uL.

### Workstream 2 -- Sequencing & Detection
- [x] **UMI error-correction pipeline** -- ported to `src/ultrasensitive/umi.py`.
- [ ] Variant calling (MuTect2 + VarDict for low-VAF); scaffold in `src/ultrasensitive/variant_calling.py`.
- [ ] Methylation analysis (Bismark + methylKit); reference in `src/ultrasensitive/methylation.py`.
- [ ] CNV profiling (ichorCNA); scaffold in `src/ultrasensitive/cnv_analyzer.py`.
- [ ] ddPCR assay design for hotspot mutations.
- [x] **Limit of Detection (LOD) validation** -- ported to `src/ultrasensitive/spike_in_simulator.py`.

### Workstream 3 -- AI/ML Classifier
- [x] **XGBoost baseline** with SHAP interpretability -- ported to `src/ultrasensitive/xgboost_classifier.py`.
- [x] **Deep Methylation Encoder** (Transformer) -- ported to `src/ultrasensitive/methylation_encoder.py`. (Real-data validation deferred to Phase 4 of `METHYLATION_PROJECT.md`.)
- [x] **Multi-Modal FusionNet** (cross-attention, 4-modality) -- ported to `src/ultrasensitive/fusion_net.py`. (Deprecation: route through `src/multimodal_fusion/` post-v2.4.)
- [x] **beta-VAE anomaly detector** -- ported to `src/ultrasensitive/vae_anomaly.py`.
- [x] Training data generator (synthetic) -- ported.

### Workstream 4 -- Clinical Reporting
- [ ] Secure data portal (HIPAA/GDPR compliant).
- [ ] PDF report generation.
- [ ] Physician dashboard with trend tracking.
- [ ] Alert system for high-risk results.
- [ ] Longitudinal tracking (rescreening intervals).

### Workstream 5 -- Validation & Regulatory
- [ ] Analytical validation (LOD, precision, reproducibility).
- [ ] Retrospective clinical study (>=2,000 samples).
- [ ] IRB/Ethics submission.
- [ ] LDT registration or FDA 510(k) pathway.
- [ ] Multi-site prospective trial.

### Workstream 6 -- Adaptive Screening Policy (NEW after merge)
- [x] **Screening environment simulator** -- ported to `src/screening_policy/rl_agent.py`.
- [x] **SAC-like RL agent** (linear Q-function demo) -- ported.
- [ ] Phase 5: upgrade linear Q-function to a neural-network Q head (`src/foundation/FoundationDownstream` extension).
- [ ] Phase 6: A/B-test vs. fixed-interval screening on synthetic patient cohorts.

## Phase 1 Milestones (Months 1-6)

```
Month 1-2: Platform Setup
  - [x] Sequencing pipeline validated (UMI scaffold in src/ultrasensitive/).
  - [x] UMI library protocol optimized (algorithm ported).
  - [ ] ddPCR assays for top 20 hotspot mutations.
  - [x] Spike-in LOD experiments -- simulator shipped.

Month 3-4: Bioinformatic Pipeline
  - [x] Variant calling scaffold (CHIP filter included).
  - [x] Methylation reference implementation (TAPS / bisulfite-ready loaders).
  - [x] CNV detection scaffold (CBS segmentation).
  - [ ] Production-readiness gate: real-data validation per workstream 5.

Month 5-6: ML Classifier v1
  - [x] XGBoost baseline + SHAP.
  - [x] Deep Methylation Encoder scaffold.
  - [x] Multi-Modal FusionNet reference implementation.
  - [x] beta-VAE anomaly detector.

Month 7-8: Alpha Testing
  - [ ] 500 clinical samples (blinded).
  - [ ] Performance assessment: sensitivity, specificity.
  - [ ] Iterate on panel design based on findings.

Month 9-12: Beta + Validation
  - [ ] 2,000-sample retrospective study.
  - [ ] First prospective screening cohort enrolled.
  - [ ] Regulatory consultation complete.
```

## Estimated Budget (Year 1)

| Category | Estimated Cost |
|----------|----------------|
| Sequencing equipment (NovaSeq, ddPCR) | $300,000 - $500,000 |
| Lab setup & consumables | $150,000 |
| Bioinformatics infrastructure | $80,000 |
| ML development & cloud compute | $60,000 |
| Clinical validation studies | $200,000 |
| Regulatory & legal | $100,000 |
| Personnel (4 FTE) | $500,000 |
| **Total** | **$1.39M - $1.59M** |

## Success Metrics (KPIs)

| Metric | Target |
|--------|--------|
| Stage I sensitivity | >=80% |
| Overall specificity | >=99.5% |
| Limit of Detection (LOD) | 0.001% VAF confirmed |
| Report turnaround | <72 hours |
| User satisfaction (physicians) | >=4.5/5.0 |
| Cost per screen | <$500 (target: <$200) |

## Collaboration Opportunities

- **Academic partners**: Hospital oncology departments for sample access.
- **Research institutions**: TCGA, ICGC for genomic data.
- **Tech partners**: Google Health, DeepMind for ML infrastructure.
- **Pharma**: Companion diagnostics for clinical trials.
- **Public health**: Government screening programs (low-resource settings).

## Next Action Items

1. **Define scope**: Which cancer types are priority? (Lung? Colon? Pan-cancer?)
2. **Confirm resources**: Budget, lab space, team size.
3. **Select first target mutation panel**: Top 20-50 hotspots for proof-of-concept.
4. **Start with spike-in validation**: Confirm 0.001% LOD in your lab.
5. **Reach out to me** with questions, refinements, or when you're ready for Phase 2 design.

---

*This is a living document. Update as scope and resources become clearer.*
