# DeepCatch: An Open-Source Panel-Based Ultra-Sensitive MRD Detection Benchmark with Real Tumor Mutations

**Authors:** Royce et al.

**Repository:** github.com/rollroyces/deepcatch (v2.2.0)

## Abstract

**Background.**  Molecular residual disease (MRD) detection from cell-free DNA
requires tracking a patient's tumor mutations in plasma at ctDNA fractions as
low as 0.1%.  Commercial assays (Signatera, CAPP-Seq) demonstrate this
clinically, but their validation data is not open-access.  Benchmarking
computational methods for ultra-early detection requires real mutations as
ground truth and an honest error model — both of which are rarely available
together.

**Results.**  We present DeepCatch, an open-source panel-based MRD detection
pipeline that uses real TCGA tumor mutations (20 LUAD patients, 5,738
mutations from GDC) as ground truth and simulates plasma cfDNA with a
context-aware sequencing error model (CpG sites: 10×, homopolymers: 5×,
clean baseline).  Three scoring methods are benchmarked at five ctDNA
fractions (10% to 0.1%) with fixed-specificity metrics (no threshold
optimization on test data).  At 0.1% ctDNA (ultra-early regime), panel-level
LLR aggregation achieves AUC 0.922 (Sens@95% = 0.600).  The CAPP-Seq
standard Fisher-method panel score achieves AUC 0.849, and a strand-
concordance-weighted variant reaches AUC 0.836.  An assay parameter sweep
(error rate × depth) shows that duplex-UMI consensus (error ≤ 1e-4) or
sequencing depth ≥ 50,000× each independently achieve Sens@95% = 1.000 at
0.1% ctDNA.

Fragmentomics features were cross-validated against two independent real
cfDNA datasets: FLARE/GSE317007 (ONT sequencing, 12 HNSCC samples) and the
published Jiang 4-mer HCC signature (129 plasma samples).  The CpG-island-
aware methylation module was validated against GSE185307 (ONT methylation,
13 LUAD+healthy), recovering cancer-specific signal at CpG islands
(AUC 0.738 vs 0.357 genome-wide).

**Conclusions.**  DeepCatch provides a reproducible, honest panel-benchmark
for ultra-sensitive MRD detection using open-access data.  All code (228/228
tests), all data (TCGA GDC, GEO), and all results are freely available.
The assay sweep identifies actionable production specifications: duplex-UMI
error suppression and 50,000× depth are sufficient for AUC 1.000 at 0.1%
ctDNA.

**Keywords:** cfDNA, cell-free DNA, MRD, molecular residual disease, liquid
biopsy, early cancer detection, panel-based detection, duplex sequencing,
fragmentomics, methylation, TCGA, open science.
