---
title: 'DeepCatch: an open-source, artifact-first benchmark for cell-free DNA fragmentomics'
tags:
  - Python
  - bioinformatics
  - reproducibility
  - fragmentomics
  - cancer-detection
authors:
  - name: Yu Ching Lam (Independent Researcher)
    orcid: 0009-0008-9113-769X
    affiliation: 1
affiliations:
  - name: DeepCatch Project
    index: 1
date: 16 November 2026
bibliography: paper/METHODS_PAPER_CITATIONS.bib
---

# Summary

DeepCatch is a Python package and benchmark suite for cell-free DNA
(cfDNA) fragmentomics research on public data. It turns the harmonised
public fragment archives (FinaleDB, GDC TCGA, GEO) into a single
reproducible pipeline: per-sample fragment features, a tumor-naive
cancer-versus-healthy classifier, a tumor-informed mutation panel
detector, a methylation channel built on FinaleMe, and the validation
controls (shuffled-label null, true-confound control, leave-one-study-out
transfer) that decide whether a reported AUC means anything.

The package exists because the published fragmentomics methods are hard to
replicate. Their training data and feature pipelines are largely
proprietary, and a researcher who wants to check whether a claim transfers
usually has to rebuild the whole thing from a methods section. DeepCatch
keeps the artifact layer first: every headline number is a committed JSON
file that the code in the same repository regenerates, and two regression
tests enforce that the pipeline is deterministic and that the narrative
documents cannot drift away from those files.

The measured results are deliberately mixed, and the negative ones are
reported with the same prominence as the positive ones. A tumor-naive
5-channel fragmentomics classifier reaches 0.9684 ± 0.0027 pooled AUC on
627 open-data samples and 0.9444 [0.8731–1.0000] on a publication that
contributes no training data — but the same classifier collapses to 0.6805
when a single-type training set is held out, a second held-out publication
reaches only 0.708, and a held-out exome study transfers only partly. The methylation channel, run end-to-end on real
public fragment files, is a **null**: its AUC (0.9444) is exactly the AUC
of the library-depth baseline, and at the 103,863 loci where both arms have
coverage the median difference between cancer and control is 0.0 (each
locus is tested on the samples that cover it; permutation p = 1.0). Both findings come with the mechanism and the
control that make them interpretable.

# Statement of need

Fragmentomics-based multi-cancer detection has produced headline-grade
AUCs on discovery cohorts (DELFI reported ≈0.94 for HCC versus healthy on
232 plasma samples), but the underlying training data, feature pipelines
and pretrained models are not released. Independent replication is rare,
and negative results about the *architecture* — which ideas help and which
do not at a given sample size — are rarely published at all, because
publishing them is not rewarded and the data sits behind access
agreements.

DeepCatch targets that gap from the direction an independent researcher
can actually work: no IRB, no institutional data-use agreement, no
proprietary cohort. Every dataset it uses is public and re-fetchable by
the scripts that ship with it. Its contribution is not a new classifier;
it is the infrastructure and the honest measurements.

# Features

The package is organised as a library (`src/`), a set of runnable
experiments (`scripts/`, `real_tcga_validation.py`), a published
artifact layer (`results/`) and the documentation that ties them together.

**Fragmentomics features.** A self-contained extractor for the 5-channel
DELFI profile (5 Mb ratio and coverage, 100 kb ratio and counts) plus the
196-dimensional fragment-size distribution, accepting both FinaleDB
fragment schemas. It is bit-identical to the companion pipeline's
extractor, so results are comparable across the two repositories.

**Tumor-informed panel detection.** A per-locus log-likelihood-ratio
aggregation over a patient's own somatic mutations, benchmarked on real
TCGA-LUAD somatic mutations with real read counts (5,738 mutations across
20 patients, GDC open access) and a context-aware sequencing-error model.
At 0.1% ctDNA the panel score reaches AUC 0.921 at fixed specificity, and
an assay sweep turns that into production specifications: duplex-UMI
error suppression (≤10⁻⁴) or ≥50,000× depth each independently reach
Sens@95% = 1.000.

**Tumor-naive cross-study benchmark.** A 5-seed × 5-fold pooled
out-of-fold evaluation across FinaleDB publications, with per-study
z-score harmonisation fitted on the training fold only, a per-bin GC and
mappability bias correction, and a true-confound control (cancer from one
study, healthy from another) that separates cancer signal from study
identity.

**External validation on never-trained publications.** Publications that
contribute no training data are fetched from the public FinaleDB API and
S3 bucket and scored by a model trained only on publications 6 and 8.
Results are reported per publication and per disease, never as a single
pooled number, and a case-only publication is reported as sensitivity at
a training-fixed threshold rather than as an AUC.

**Leave-one-study-out transfer.** Training on one publication and testing
on the other, in both directions, with a within-study
leave-one-cancer-type-out control to separate "unseen study" from "unseen
cancer type".

**Methylation channel.** FinaleMe is executed end-to-end on public
fragment files — fragments to a coordinate-sorted BAM, per-CpG feature
matrices, a WGBS prior, one HMM trained on the cohort and each sample
decoded with it — so that `src/methylation_gnn/` consumes real β values
rather than synthetic ones. The channel is then measured rather than
assumed: a model-free estimator (the share of a CpG's observations that
sit at a fragment end) is compared against the decoder output, every
score is paired with its library-depth baseline, and the per-CpG
difference is tested at loci both arms cover.

**Validation controls.** Shuffled-label and true-confound controls,
leave-one-out and leave-one-study-out cross-validation, DeLong
confidence intervals, sensitivity at fixed specificity, and a
prevalence-floor projection that states what a screening PPV would be at
0.1% prevalence.

# Quality and reproducibility

Reproducibility is enforced rather than asserted:

- the PCA inside every cross-validation fold is seeded with the
  fold-split seed, and an AST-based regression test fails the build if
  any `PCA(...)` call in the package is fitted without a `random_state`;
- a second regression test fits the benchmark twice and requires the two
  runs to agree exactly, and separately requires that a different seed
  still gives a different answer, so the seed cannot be quietly pinning
  everything;
- every test that runs a pipeline writes its outputs to a temporary
  directory, so a test run cannot silently overwrite a published
  artifact;
- `scripts/benchmark_numbers.py` prints every figure quoted in the
  manuscript directly from the committed JSON, so the prose and the
  artifacts cannot drift apart.

The suite is 978 tests (`pytest src test -q`), split across two CI jobs
so that the core fragmentomics path stays installable without a deep
learning stack.

# Documentation

`README.md` covers installation, the quick start and the artifact map;
`MODEL.md` is a model card stating scope and limitations; `USAGE.md` is
the command reference; `docs/` holds the per-experiment write-ups
(cross-study benchmark, leave-one-study-out, external validation,
methylation pipeline, harmonization, the two readiness reports);
`paper/METHODS_PAPER.md` is the long-form manuscript and
`paper/README.md` the compilation and submission guide.

# Code availability

MIT-licensed, at <https://github.com/rollroyces/deepcatch>. Install with
`pip install -e .` (Python ≥3.10; the CPU fragmentomics path needs only
NumPy, SciPy, scikit-learn, pandas and statsmodels). Docker images and
the full CI configuration are included. The experiment scripts refetch
their public inputs, and the artifacts they produce are committed to the
repository so that the numbers can be checked without running anything.

# Limitations

The benchmark measures separation on public, pre-processed fragment
archives, not on clinical plasma: the external cohort's control arm is
liver-transplant patients rather than healthy donors, and the panel
experiment uses simulated plasma reads over real tumour mutations. The
methylation channel does not carry signal on the publicly reachable
libraries, because they are size-selected at ~300 bp and deliver about
one fragment observation per CpG; a genuinely MspI-digested library is 41×
deeper and does carry the signature, but the reachable one is the
training cohort, so it is a feasibility control rather than a validation.
Stage labels are absent from the open-data cohort, so no stage-I claim is
made.

# Acknowledgements

The author thanks the public data providers whose archives make
replication possible: the FinaleDB and FinaleMe teams (Zenodo
10.5281/zenodo.14013719 and 10.5281/zenodo.6914806), the GDC for TCGA
open-access somatic mutations, and the GEO depositors of GSE317007 and
GSE185307.

# References

- Cristiano S, Leal A, Phallen J, et al. Genome-wide cell-free DNA
  fragmentation in patients with cancer. *Nature*. 2019;570(7761):385–389.
  doi:[10.1038/s41586-019-1272-6](https://doi.org/10.1038/s41586-019-1272-6)
- Jiang P, Chan KCA, Lo YMD. Lengthening and shortening of plasma DNA in
  hepatocellular carcinoma patients. *PNAS*. 2015;112(11):E1317–E1325.
  doi:[10.1073/pnas.1500076112](https://doi.org/10.1073/pnas.1500076112)
- Jiang P, Sun K, Peng W, et al. Plasma DNA end-motif profiling as a
  fragmentomic marker in cancer, pregnancy, and transplantation.
  *Cancer Discovery*. 2020;10(5):664–673.
  doi:[10.1158/2159-8290.CD-19-0622](https://doi.org/10.1158/2159-8290.CD-19-0622)
- Newman AM, Bratman SV, To AT, et al. An ultrasensitive method for
  quantitating circulating tumor DNA with broad patient coverage.
  *Nature Medicine*. 2014;20(4):548–554.
  doi:[10.1038/nm.3519](https://doi.org/10.1038/nm.3519)
- Newman AM, Lovejoy AF, Ahrendt S, et al. Integrated digital error
  suppression for sensitive detection of circulating tumor DNA.
  *Nature Biotechnology*. 2016;34(3):313–320.
  doi:[10.1038/nbt.3520](https://doi.org/10.1038/nbt.3520)
- Liu MC, Oxnard GR, Klein EA, et al. Sensitive and specific multi-cancer
  detection and localization using methylation signatures in cell-free
  DNA. *Annals of Oncology*. 2020;31(6):745–759.
  doi:[10.1016/j.annonc.2020.02.011](https://doi.org/10.1016/j.annonc.2020.02.011)
- Liu X, Liu Z, Zhang Y, et al. FinaleMe: an efficient and accurate method
  for quantifying methylation patterns in cell-free DNA. *bioRxiv*.
  2018. doi:[10.1101/478222](https://doi.org/10.1101/478222)
- Straub J, Mathios D, Foda Z, et al. Finaletoolkit: an end-to-end
  infrastructure for the analysis of cell-free DNA fragmentomics data.
  *bioRxiv*. 2023. The identifier previously recorded for this preprint
  does not resolve, so none is claimed here.
- DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two
  or more correlated ROC curves: a nonparametric approach. *Biometrics*.
  1988;44(3):837–845.
  doi:[10.2307/2531595](https://doi.org/10.2307/2531595)
