# Phase 0 Plan -- FinaleMe smoke validation

**Goal:** Verify methylation features produce a meaningful cancer-vs-healthy signal
(AUC > 0.65) before scaling to the 627-cohort Phase 1.

**Source:** Ported from `rollroyces/deepcatch-methylation` PHASE0_PLAN.md (SHA 69f3ba32)
on 2026-10-16, alongside the v0.11 cfRNA coupling module ported to
`src/methylation/cfrna_coupling.py` and `src/methylation/methylation_reader.py`.

## Pivot from original plan

The original Phase 0 plan was to run **FinaleMe** on 50 raw FinaleDB WGS BAMs.
This requires:
- Downloading ~50 BAM files from FinaleDB (multi-GB).
- Java 21 + Maven + FinaleMe JAR build (~1 hour setup).
- Reference downloads (~5-10 GB).
- BAM -> fragment BED preprocessing.

That's heavy for a smoke test. **Phase 0 is being pivoted** to use **GSE122126**
(Moss 2018, *Nature Communications*), a public Illumina 450K cfDNA methylation
array dataset that:

| Property        | Value                                                                                |
|-----------------|-------------------------------------------------------------------------------------|
| Cancer / healthy| 11 cancer / 59 cfDNA total (with cell-type reference atlas)                         |
| Assay           | Illumina 450K + EPIC                                                                |
| Access          | Open via GEO (no DAC)                                                               |
| Size            | ~1 GB                                                                                |
| Already processed| Yes (SeriesMatrix + IDATs)                                                          |
| Cancer types    | Colorectal (4), lung (4), breast (3), CUP (4)                                       |

Using already-processed methylation arrays validates the **feature pipeline +
LR baseline + 5-fold CV protocol** without requiring FinaleMe yet. FinaleMe becomes
the **Phase 1 component** for the 627-cohort scale-up, where we need imputation
because FinaleDB has no methylation calls.

## Exit criterion (unchanged)

**Phase 0 succeeds if methylation features + LR baseline achieve AUC > 0.65** on
the GSE122126 cancer-vs-healthy classification task. The expected AUC for
methylation-only classifiers on this dataset is **0.85-0.95** (per Moss 2018
paper); we are looking for the lower bound 0.65 to confirm the pipeline works.

## Tasks

1. **Download GSE122126** from GEO (~1 GB) into `data/raw/GSE122126/`.
2. **Parse SeriesMatrix** to extract per-CpG beta-values + sample labels.
3. **Filter probes:** drop low-variance probes + cross-reactive probes + non-CpG probes.
4. **Cross-study split:** hold out 30% of samples, train on 70%, evaluate.
5. **LR baseline:** L2-LR with C in {0.001, 0.01, 0.1, 1.0, 10.0} sweep; 5-fold pooled CV.
6. **Report:** per-seed AUC, mean +/- std, 95% DeLong CI, paired t-test vs fragmentomics AUC.
7. **Negative-result protocol:** if AUC < 0.65, document the failure mode and revisit the strategy before Phase 1.

## Deliverables

| File                                       | Description                                                 |
|--------------------------------------------|-------------------------------------------------------------|
| `scripts/fetch_gse122126.py`               | Download + parse SeriesMatrix from GEO                      |
| `src/methylation/methylation_baseline.py`  | Methylation -> feature matrix -> LR baseline                 |
| `test/test_methylation/test_baseline.py`   | Unit tests for the baseline pipeline                          |
| `results/phase0_smoke_validation.json`     | Per-seed AUCs + summary statistics                           |
| `docs/ultrasensitive/PHASE0_RESULTS.md`    | Honest write-up of what worked, what didn't                  |

## Timeline

- Download: 30 min (1 GB over typical connection)
- Parse + filter: 1 hour (manual scripting)
- Baseline + 5-fold CV: 2-4 hours
- Write-up: 1-2 hours

**Total: 1-2 days wall-clock** for Phase 0.

## What this validates vs what it doesn't

### Validated
- Methylation feature extraction pipeline (450K array -> feature matrix).
- LR baseline + 5-fold CV protocol on methylation data.
- Cross-study split methodology.
- Per-seed AUC reporting with DeLong CI.

### Not validated
- FinaleMe performance on FinaleDB WGS (deferred to Phase 1).
- Tissue-of-origin classifier (Phase 3).
- Methylation GNN training (Phase 4, blocked on collaborator).
- Methylation + fragmentomics fusion (Phase 2).

## Honest read

This pivot is **better than the original plan** for Phase 0 because:

1. GSE122126 methylation arrays are **higher quality** than FinaleMe-imputed
   methylation from WGS. If methylation features don't work on the arrays, they
   certainly won't work on the imputation.
2. The download is **30x smaller** (1 GB vs 30+ GB for FinaleDB BAMs).
3. The pipeline is **simpler** (no BAM preprocessing, no Java install).
4. We get a **known-good benchmark** (Moss 2018 reports AUC 0.85-0.95 for
   cancer-vs-healthy on this exact dataset) so we can sanity-check our pipeline.

Phase 1 will then add FinaleMe as the imputation method for FinaleDB WGS samples.
If FinaleMe is comparable to array methylation in accuracy (auROC ~0.85-0.91 per
FinaleMe paper), the 627-cohort baseline should be in the AUC 0.80-0.90 range.
If it's much worse (<0.70), we know the imputation is the bottleneck, not the
features.

## Sister-repo context

The `cfRNA coupling` (`src/methylation/cfrna_coupling.py`) and
`methylation reader` (`src/methylation/methylation_reader.py`) modules ported
alongside this plan work together with `mrnavax_rare` -- they are orthogonal to
the cancer-vs-healthy head-to-head and operate on observed cfRNA persistence in
plasma samples. They do NOT consume FinaleDB or FinaleMe data; their contract is
the structural `CfRNAPersistenceInput` typed-dict shape.
