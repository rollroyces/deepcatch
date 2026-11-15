---
name: External cohort validation (collaboration)
about: Offer a cfDNA cohort for external validation, or ask for the analysis pipeline
title: "[collaboration] "
labels: clinical-data
assignees: ''
---

## What you have

Thanks for considering this. This template exists so that a cohort can
be discussed without any patient data leaving your institution.

**Please do not attach patient data, identifiers, or read files to an
issue.** Nothing in the workflow needs them at this stage.

| Field | Your answer |
|---|---|
| Cohort name / internal study ID | |
| Sample type | plasma / serum / urine / other |
| Sequencing | WGS / targeted panel / shallow / other |
| Genome build | hg19 / hg38 |
| Approximate samples | cancer: __  non-cancer: __  serial draws: __ |
| Labels available | cancer/non-cancer · stage · treatment · tissue of origin · matched tumour |
| Data access | public repository (accession) / EGA or dbGaP (accession) / institutional DUA |
| Can it be shared with a third party? | yes / no / under review |

## What we would do with it

The analysis side already exists and is open source:

- the 5-channel DELFI + fragment-size feature extractor
  (`src/fragmentomics/finale_features.py`, and the same schema the
  public FinaleDB cohorts use),
- the per-study harmonisation and the statistical protocol (DeLong
  CIs, sens@spec grid, PPV at screening prevalence, leave-one-cohort-out),
- the external-validation harness
  (`scripts/external_validation_finaledb.py`),
- and the leave-one-study-out benchmark (`scripts/cross_study_loso.py`).

The result would be a held-out cohort the model has never seen, scored
with the same code that produced the public-cohort numbers, with every
figure traceable to a committed JSON artifact.

## What is useful at each scope

1. **Unpaired plasma, cancer + non-cancer** — enough to run the same
   external-validation protocol as the public cohort.
2. **Serial draws** (same patient over time) — the longitudinal
   Stage-2 tracker cannot be evaluated on anything else.
3. **Stage I vs late-stage labels** — the stage breakdown is implemented
   and currently waiting for data.
4. **Matched tumour tissue** — enables tumour-specific fragmentomics.

## Questions

- Is the data usable under the terms you are able to offer?
- Who should be on the paper, and in which order?

Questions about the analysis are welcome in
[Discussions](https://github.com/rollroyces/deepcatch/discussions) instead —
this form is for cohort offers.

<!--
Co-authoring model: see TEAM.md §2 for the contribution and governance
model, including what is and is not accepted.
-->
