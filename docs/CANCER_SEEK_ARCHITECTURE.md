# Ultra-Sensitive Stage 1 Cancer Screening System
## Technical Architecture & Design Document

**Version:** 1.0
**Source:** Ported from `rollroyces/cancer-screening` ARCHITECTURE.md on 2026-10-16.
**Goal:** Detect cancer signals at 0.001% allele frequency (1 in 100,000 cells).
**Designation:** CANCER-SEEK v1.0.

> Modules referenced below were ported from `cancer-screening/src/core/*` into
> `deepcatch/src/ultrasensitive/`. See `ROADMAP_CANCER_SEEK.md` for the integration
> roadmap and `ultra-early-readiness` for binding-gap coverage.

---

## 1. Problem Statement

**Target:** Stage I cancer detection -- tumors typically <2cm, <10^8 cells.
**Challenge:** In early stage, tumor DNA fraction in blood is extremely low:
- Stage I: ~0.01% -- 0.1% (1-100 molecules per 100,000).
- Our goal: **0.001% sensitivity** (ultra-early detection floor).

This requires **digital detection** -- counting individual DNA molecules rather than bulk signals.

## 2. Detection Modalities (Multi-Analyte Panel)

### 2.1 ctDNA -- Circulating Tumor DNA (Primary Signal)
- Tumor-derived DNA fragments shed into bloodstream.
- Carries **somatic mutations** (SNVs, indels) unique to the tumor.
- Half-life: ~15-120 minutes (real-time snapshot of tumor status).

### 2.2 Methylation Signatures (Epigenetic)
- Abnormal methylation patterns are among the earliest molecular changes.
- **Cancer-specific methylation signatures** detectable before mutations.
- Can also identify **tissue of origin** (lung, colon, breast, etc.).

### 2.3 Copy Number Variation (CNV) Profiling
- Chromosomal amplifications/deletions in ctDNA.
- Detectable at low tumor fractions with deep sequencing.

### 2.4 microRNA (miRNA) Panels
- Stable in plasma, cancer-specific expression profiles.
- Complementary signal to DNA-based detection.

### 2.5 Protein Biomarkers (Multi-Omics)
- e.g., CA-125, CEA, PSA in combination.
- Lower sensitivity alone, but high specificity in panel.

## 3. Detection Technologies

### 3.1 Molecular Barcode Error Suppression (UMB & UMI)

```
Problem: NGS error rate ~0.1% -- masks true 0.001% signals.
Solution: Unique Molecular Identifiers (UMIs)
  Each cfDNA molecule tagged with random 12-bp barcode before PCR.
  Duplex sequencing: read both strands, require consensus.
  True mutations survive; PCR/NGS errors are filtered out.
  Effective error rate: 10^-8 to 10^-9 (enables 0.001% detection).
```

### 3.2 Droplet Digital PCR (ddPCR) -- Targeted Verification
- Partition each reaction into ~20,000 nanoliter droplets.
- Absolute quantification without standard curves.
- **Sensitivity: 0.01-0.1%** for known mutations.
- Use as secondary confirmation for positive ctDNA signals.

### 3.3 Next-Generation Sequencing (NGS) -- Untargeted Screen
- **Shallow whole-genome sequencing** (sWGS): CNV detection.
- **Targeted deep sequencing**: 500+ cancer genes panel.
- **Methylation sequencing**: TAPS (TET-assisted pyridine borane sequencing).
- **Depth requirement**: 100,000x coverage for 0.001% sensitivity.

### 3.4 Machine Learning Classifier

```
Input features:
  - Mutation allele frequencies (UMI-corrected).
  - Methylation beta-values (tumor vs. normal).
  - CNV log2 ratios.
  - miRNA expression levels.
  - Protein biomarker concentrations.

Model: Ensemble classifier
  - XGBoost + Neural Network (methylation encoder).
  - Trained on: 50,000+ cancer + 50,000 healthy samples.
  - Output: Probability score + cancer type prediction.
```

## 4. System Architecture

```
+----+   +-------+   +----------------+   +-----------------+
| 1  |-->|   2   |-->|       3        |-->|        4        |
|Intake| |Library| |   Sequencing    |  | AI Analytics    |
|10-20mL| Prep    |   (Illumina)      |  |  (PyTorch)      |
|Blood | |Streck | |  NovaSeq 100kx   |  |                 |
+-----+ +-------+   +----------------+   +-----------------+
                                                 |
                                                 v
                                          +----------------+
                                          |  Clinical      |
                                          |  Report (PDF)  |
                                          +----------------+
```

### 4.1 Sample Intake
- Blood draw (10-20 mL) -> Plasma isolation -> cfDNA extract.
- Expected yield: ~5-30 ng cfDNA / mL plasma.

### 4.2 Molecular Processing
- **ctDNA UMI Library Prep** -- NEBnext / Twist Duplex Sequencing.
- **Methylation Conversion** -- TAPS / bisulfite.
- **miRNA/Protein Extraction** -- plasma aliquot.

### 4.3 Sequencing
- Illumina NovaSeq (NGS) target depth 100,000x.
- Bio-Rad QX200 (ddPCR) for targeted verification.

### 4.4 AI Analytics Pipeline
1. Raw reads -> Quality control (FastQC).
2. Adapter trimming -> Barcode demultiplexing.
3. UMI deduplication -> Error correction (`src/ultrasensitive/umi.py`).
4. Alignment (BWA-MEM2) -> Variant calling (`src/ultrasensitive/variant_calling.py`).
5. Methylation calling (Bismark/Nanopolish; structure in `src/ultrasensitive/methylation.py`).
6. CNV profiling (`src/ultrasensitive/cnv_analyzer.py`).
7. ML inference: Cancer/No-Cancer + Tissue of Origin
   (`src/ultrasensitive/{xgboost_classifier,fusion_net,methylation_encoder}.py`).
8. Report generation (PDF).

### 4.5 Clinical Report
- Cancer Probability Score (0-100%).
- Detected mutation list with VAF.
- Tissue of origin prediction.
- Recommended follow-up (imaging, biopsy).
- Historical trend (if repeat screening).

## 5. Sensitivity & Specificity Targets

| Stage     | Tumor Size | ctDNA Fraction | Target Sensitivity |
|-----------|-----------|-----------------|-------------------|
| Stage I   | <2cm      | 0.01-0.1%       | **>=80%**          |
| Stage II  | 2-5cm     | 0.1-1%          | >=90%              |
| Stage III | >5cm      | 1-10%           | >=95%              |
| Stage IV  | Metastatic| 10-50%          | >=99%              |

**Overall Specificity target: >=99.5%** (false positive rate <0.5%).
**Critical threshold: 0.001% VAF detection floor**.

## 6. Technical Specifications for 0.001% Sensitivity

```
Detection Floor Calculation:
  - Input: 10 mL blood -> ~5,000 genome equivalents (GE) / mL.
  - Total: ~50,000 GE in 10 mL plasma.
  - 0.001% = 0.5 mutant molecules per 50,000 GE.
  - To detect 0.5 molecules reliably:
      Need 50,000x coverage = 50,000 x 30 bp read depth, OR
      Consensus from multiple UMI tags per molecule.

UMI Strategy:
  - 12-bp random UMI: 16.7M unique tags.
  - Expected duplicates per molecule: ~3-5x.
  - Consensus building: 3/3 or 2/3 strand agreement.
  - Final error rate: <10^-8.
```

## 7. Key Implementation Challenges & Mitigations

| Challenge | Mitigation |
|-----------|------------|
| Low ctDNA input in Stage I | Large volume plasma (15-20 mL), ultra-low input protocols |
| Germline variant confusion | Matched white blood cell sequencing (filter germline); `variant_calling.LowVAFVariantFilter.filter_chip_variants` |
| Chip/fiber contamination | Dedicated clean rooms, gowning protocols |
| Variants of uncertain significance | ACMG classification, literature linking |
| Tumor heterogeneity | Multi-gene panel (500+ genes), methylation triangulation |
| False positives from clonal hematopoiesis (CHIP) | Simultaneous WBC sequencing, CHIP flagging pipeline |

## 8. Cancer Gene Panel Design (Targeted Sequencing)

**Core 500-gene panel focusing on:**
- **Driver mutations**: TP53, KRAS, EGFR, BRAF, PIK3CA, ALK, ROS1, RET.
- **DNA repair genes**: BRCA1/2, ATM, PALB2, MLH1, MSH2.
- **Methylation drivers**: CDKN2A, RASSF1A, SHOX2, SEPT9.
- **Tissue-specific markers**:
  - Lung: EGFR, KRAS, ALK, ROS1, TP53.
  - Breast: PIK3CA, ESR1, HER2 (ERBB2).
  - Colorectal: KRAS, NRAS, BRAF, MSI genes.
  - Liver: TP53, TERT promoter, CTNNB1.
  - Pancreas: KRAS, TP53, CDKN2A, SMAD4.

## 9. Clinical Validation Roadmap

```
Phase 1 (In-silico):
  - Simulate ctDNA signals at 0.001-0.1% VAF
    (implementation: ``src/ultrasensitive/spike_in_simulator.py``).
  - Benchmark UMI error correction algorithms.

Phase 2 (Analytical):
  - Spike-in experiments: mix tumor DNA into healthy plasma.
  - Validate limit of detection (LOD) = 0.001% VAF.

Phase 3 (Retrospective):
  - Test on stored samples from Stage I patients + matched controls.
  - Expected n >= 2,000 (1,000 cancer / 1,000 healthy).

Phase 4 (Prospective):
  - Multi-site clinical trial.
  - Annual screening cohort.
```

## 10. Next Steps & Recommendations

1. **Define cancer type priority** -- lung, colorectal, breast, or pan-cancer first?
2. **Access to sequencing infrastructure?** -- Illumina, Oxford Nanopore, or hybrid approach.
3. **Sample source** -- clinical samples, cell lines (HCC827, A549), or synthetic spike-ins?
4. **Budget & hardware** -- ddPCR only (feasible ~$50K), full NGS lab (~$500K+).
5. **Regulatory pathway** -- LDT (Laboratory Developed Test) vs. FDA clearance.

---

## Appendix: Module Roadmap Inside DeepCatch

| Architecture component                | Module                                            | Port status |
|---------------------------------------|---------------------------------------------------|-------------|
| UMI error correction                  | `src/ultrasensitive/umi.py`                       | ported |
| Low-VAF variant calling + CHIP filter | `src/ultrasensitive/variant_calling.py`           | ported |
| Methylation beta-value processor      | `src/ultrasensitive/methylation.py`               | ported |
| CNV CBS segmentation                  | `src/ultrasensitive/cnv_analyzer.py`              | ported (bug fix: `__Logger` -> `__name__`) |
| XGBoost baseline + SHAP               | `src/ultrasensitive/xgboost_classifier.py`        | ported |
| FusionNet (4-modality g+m+p+c)        | `src/ultrasensitive/fusion_net.py`                | ported as reference; route through `src/multimodal_fusion/` post-v2.4 |
| Deep Methylation Encoder Transformer  | `src/ultrasensitive/methylation_encoder.py`       | ported (architecture-only; real-data validation is Phase 4 of `METHYLATION_PROJECT.md`) |
| beta-VAE anomaly detector             | `src/ultrasensitive/vae_anomaly.py`               | ported |
| Adaptive Screening RL agent           | `src/screening_policy/rl_agent.py`                | ported (with `screening_policy/__init__.py`) |
| Spike-in LOD validation               | `src/ultrasensitive/spike_in_simulator.py`        | ported (with `out_csv=...` for committed `results/`) |
| cpG FinaleMe methylation              | per `METHYLATION_PROJECT.md` Phase 0/1            | deferred until collaborator (Phase 4 blocked) |

---

*This document is the technical design framework. Clinical implementation requires
regulatory compliance (FDA/NMPA), IRB approval, and professional laboratory
certification. Sourced from `rollroyces/cancer-screening` and adapted to the
`deepcatch/src/ultrasensitive/` location after the 2026-10-16 consolidation.*
