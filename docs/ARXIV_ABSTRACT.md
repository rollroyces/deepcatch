<!-- GENERATED FILE - do not hand-edit.
     Source: paper/paper.tex (the manuscript the PDF is compiled from).
     Regenerate: python3 scripts/generate_preprint_abstract.py --variant arxiv
     Guarded by: test/test_generated_docs_match_their_generators.py
     -->

# An Open-Source Benchmark for Panel-Based Ultra-Sensitive Detection of Circulating Tumor DNA from Real Tumor Mutations

**Authors:** Yu Ching Lam — The Chinese University of Hong Kong (ORCID 0009-0008-9113-769X)

**Repository:** https://github.com/rollroyces/deepcatch (v2.3.0)

## Abstract for arXiv's metadata field

Paste this into the form's Abstract field, exactly as it stands:
the two spaces after each line break are what keep the paragraph
breaks from being stripped by arXiv's own text processing.

```
Background. Molecular residual disease (MRD) detection from cell-free DNA (cfDNA) requires tracking a patient's tumor mutations in plasma at circulating tumor DNA (ctDNA) fractions as low as 0.1%. Commercial assays such as Signatera and CAPP-Seq demonstrate this clinically, but their validation data are not open.
  Results. We present DeepCatch, an open-source panel-based MRD detection pipeline benchmarked against real tumor mutations. Using 20 real TCGA-LUAD patients (5,738 somatic mutations with real read counts, GDC open access) as ground truth and a context-aware error model, we simulate plasma cfDNA at five tumor fractions (10%-0.1%) and score three methods at fixed specificity, without test-set threshold optimisation. At 0.1% ctDNA, panel-level log-likelihood-ratio aggregation achieves AUC 0.921 (Sens@95% = 0.770), and an assay sweep shows duplex-UMI consensus (error <= 10^-4) or depth >= 50,000x each reach Sens@95% = 1.000. A tumor-naive fragmentomics channel (627 cross-study FinaleDB samples) averaged with a mutation channel reaches AUC 0.9884, but its mutation channel is a synthetic proxy: a complementarity result, not a performance claim. On never-trained publications, transfer reaches AUC 0.9444 (95% CI 0.8731-1.0000) on one held-out whole-genome study and 0.708 on a second, pooling to 0.7385; transfer is real but not uniform, an exome-targeted study only partly transfers, and holding a study out of training costs most of the pooled single-type AUC (0.6805). The pooled figure is a within-cohort separation number, not a generalisation claim.
  Conclusions. DeepCatch provides a reproducible, honest benchmark for ultra-sensitive MRD detection on open-access data and turns the assay sweep into production specifications: duplex-UMI error suppression and ~50,000x depth suffice for AUC 1.000 at 0.1% ctDNA. Clinical validation on real plasma cfDNA remains the essential next step.
```

<!-- 1915 characters against arXiv's 1920-character limit on the metadata field.
     The paper keeps the manuscript's own abstract: arXiv's limit is on the form field,
     not on the manuscript. -->
