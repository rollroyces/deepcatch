<!-- GENERATED FILE - do not hand-edit.
     Source: paper/paper.tex (the manuscript the PDF is compiled from).
     Regenerate: python3 scripts/generate_preprint_abstract.py
     Guarded by: test/test_generated_docs_match_their_generators.py
     -->

# An Open-Source Benchmark for Panel-Based Ultra-Sensitive Detection of Circulating Tumor DNA from Real Tumor Mutations

**Authors:** Yu Ching Lam — Independent Researcher (ORCID 0009-0008-9113-769X)

**Repository:** https://github.com/rollroyces/deepcatch (v2.3.0)

## Abstract

Background. Molecular residual disease (MRD) detection from cell-free DNA (cfDNA) requires tracking a patient's tumor mutations in plasma at circulating tumor DNA (ctDNA) fractions as low as 0.1%. Commercial assays such as Signatera and CAPP-Seq demonstrate this clinically, but their validation data and pipelines are not openly accessible.

Results. We present DeepCatch, an open-source panel-based MRD detection pipeline benchmarked against real tumor mutations. Using 20 real TCGA-LUAD patients (5,738 somatic mutations with real read counts, GDC open access) as ground truth and a context-aware error model (CpG sites ~10×, homopolymers ~5× error), we simulate plasma cfDNA at five tumor fractions (10%–0.1%) and evaluate three scoring methods at fixed specificity, no test-set threshold optimisation. At 0.1% ctDNA, panel-level log-likelihood-ratio (LLR) aggregation achieves AUC 0.921 (Sens@95% = 0.770, Sens@99% = 0.460), with paired cancer-vs-control win rate 1.000; the CAPP-Seq-standard Fisher-method panel score achieves AUC 0.834 and a strand-concordance-weighted variant 0.831. An assay sweep shows duplex-UMI consensus (error ≤ 10⁻⁴) or depth ≥ 50,000× each independently achieve Sens@95% = 1.000 at 0.1% ctDNA. Fragmentomics features were cross-validated against independent real cfDNA (FLARE/GSE317007), reproducing the published CG-depletion / AT-enrichment signature. A tumor-naive fragmentomics channel (627 cross-study FinaleDB samples, Jiang 2015 + Cristiano 2019) naively averaged with a mutation channel calibrated to AUC 0.902 reaches AUC 0.9884 (10 seeds: +0.0144 over the tumor-naive channel alone, t=54.50, positive on every seed); the mutation channel there is a synthetic proxy, not the real panel LLR: a complementarity result, not a performance claim. Because a pooled cross-study AUC is measured inside the cohort that trained it, we also report transfer on never-trained publications: AUC 0.9444 (95% CI 0.8731–1.0000) on one held-out whole-genome study (n=31) and 0.708 (0.5148–0.9012, n=52) on a second, pooling to 0.7385 (0.6281–0.8488, n=83) over held-out studies with both classes. Transfer between whole-genome studies is real but far from uniform; an exome-targeted study only partly transfers. Holding one study out of training costs most of the pooled AUC single-type (0.6805). The pooled cross-study figure is a within-cohort separation number, not a generalisation claim; we release the held-out cohort script so it is checkable.

Conclusions. DeepCatch provides a reproducible, honest benchmark for ultra-sensitive MRD detection using open-access data. The assay sweep yields concrete production specifications: duplex-UMI error suppression and ~50,000× depth suffice for AUC 1.000 at 0.1% ctDNA. Clinical validation on real plasma cfDNA remains the essential next step.

Keywords: cfDNA; ctDNA; MRD; liquid biopsy; early cancer detection; panel-based detection; duplex sequencing; open science

<!-- 2940 characters against bioRxiv's 3000-character limit;
     this is the text to paste into the submission form. -->
