# EGA / dbGaP — Data Access Request Templates (Public cfDNA Datasets)

**Audience:** Royce (the user)
**Why this exists:** The repo's `PRODUCTION_ROADMAP.md §4` and `NEXT_STEPS.md §3` both flag
real plasma cfDNA WGS access as **the binding constraint** for clinical-grade ultra-early
cancer signaling validation. Each of these datasets needs an access application to a
Data Access Committee (DAC), and turnaround is **2–6 months**. **Start the paperwork
today** so the clock is already running when the rest of the readiness work matures.

The four datasets below are the realistic targets for DeepCatch's tumor-naive /
multi-modal validation work. The dataset descriptions are reproduced from the
literature references already cited in `docs/CROSS_PLATFORM_BENCHMARK.md` and
`results/real_tcga_validation.json`. Edit the bracketed `[...]` placeholders for
your specific institution + contact details before sending.

---

## Quick reference: targets

| Dataset | Accession | Samples | Has controls? | Use case | Effort |
|---|---|---|---|---|---|
| **Cristiano 2019 (DELFI)** | EGA: `EGAS00001003828` / dbGaP: `phs0034536` | 545 (215 cancer, 330 healthy) | ✅ Yes | Tumor-naive fragmentomics on plasma WGS, multi-cancer. Direct match for DeepCatch's 5-channel pipeline. | 2–6 months |
| **Newman 2016 (CAPP-Seq NSCLC MRD)** | EGA (CAPP-Seq) | 40 patients, serial draws | N/A (MRD) | Tumor-informed MRD benchmark on DeepCatch's panel-LLR layer. | 3–6 months |
| **Abbosh 2017 (TRACERx lung MRD)** | EGA: `EGAD00001002469` | 24 patients, 96 serial plasma samples | N/A (MRD) | Higher-resolution MRD benchmark than CAPP-Seq; phase-matched tumor + plasma. | 3–6 months |
| **Mouliere 2018 (Sci Transl Med)** | EGA (per-study) | cfDNA size distributions + matched tumor | ✅ Yes | Fragment-size-only profile (CfDNA fragmentomics orthogonal). | 3–6 months |

**Not recommended for this round** (require collaboration / sponsorship):
- GRAIL CCGA1 / CCGA3 — by collaboration/agreement only; not an open-access
  application. Defer until after the open-data Tier 1 results are published.
- TCGA — does NOT contain plasma cfDNA (tissue only). It defines panels, not
  plasma truth. Already in use for the panel-LLR detector (5,738 LUAD
  mutations, no access required).

---

## Template A — dbGaP Data Access Request (DAR)

**Where:** https://dbgap.ncbi.nlm.nih.gov/aa/dar
**Format:** NIH's online form; the prose below is for the "Research Use Statement"
field (≤ 5,000 characters).

```
[PI NAME], [INSTITUTION], requests controlled access to the cfDNA whole-genome
sequencing (WGS) data from the Cristiano et al. 2019 (Nature) DELFI cohort
(dbGaP phs0034536 / EGA EGAS00001003828). The dataset contains 545 plasma
WGS samples (215 cancer, 330 healthy controls) with paired fragment-size and
coverage profiles.

This work supports the open-source DeepCatch project
(github.com/rollroyces/deepcatch), a multi-modal fragmentomics + methylation
framework for ultra-early cancer detection. DeepCatch's 5-channel pooled
cross-study benchmark on FinaleDB publications 6 + 8 (n = 627) reports
AUC 0.967 ± 0.004 with shuffled-label null 0.466 — a clinically meaningful
signal that requires independent validation on the DELFI cohort. The proposed
analysis will (a) re-run DeepCatch's tumor-naive 5-channel pipeline on the
DELFI cohort using only the open-source code, (b) report sens@spec at the
clinical-decision grid (95% / 98% / 99%) and the screening-grade grid
(99.5% / 99.9%), and (c) audit for batch effects via per-study z-score
harmonization. No clinical claims will be made without matched-WBC CHIP
subtraction.

Research use statement: cancer-detection methods research, not clinical
decision support. Outputs will be released as open-source code + open-data
benchmarks. The principal investigator has prior experience with TCGA-LUAD
data access (phs000178.v11.p8) under dbGaP project #XXXXX.

[IRB / ethics approval: not required — dataset is de-identified per dbGaP
policy.]
```

---

## Template B — EGA Data Access Agreement

**Where:** https://ega-archive.org/datasets/[ACCESSION]
**Format:** EGA's standard Data Access Agreement form; submit per study.

```
[APPLICANT NAME]        [INSTITUTION]        [EMAIL]        [DATE]

To: EGA Data Access Committee — [STUDY PI / CONTACT, see EGA metadata]

Re: Data Access Agreement for [STUDY TITLE] ([EGA ACCESSION])

I, [APPLICANT NAME], request access to the controlled-tier human genomic
data in the above EGA study under the terms of the EGA Data Access
Agreement. The data will be used solely for the research purpose stated
below; no re-identification, redistribution, or attempted contact with
study participants is permitted.

Research purpose: independent validation of the open-source DeepCatch
fragmentomics pipeline (github.com/rollroyces/deepcatch) on real plasma
whole-genome sequencing data from the [STUDY] cohort. The work will
produce open-source benchmarks, not clinical decision support.

Specifically:
  1. Re-run DeepCatch's 5-channel tumor-naive pipeline on the [STUDY]
     cohort using only open-source code, on data already extracted to
     /Users/hermes/cfdna-fragmentomics-pipeline/data/features/.
  2. Report sens@spec at 95%/98%/99%/99.5%/99.9% specificity and
     per-cancer-type DeLong 95% CIs.
  3. Audit for batch effects via per-study z-score harmonization and
     the existing shuffled-label null control pipeline.
  4. Publish results as an open-access preprint + open-source benchmark,
     citing the original [STUDY] paper.

This project does not constitute human-subjects research on the
applicant's side (de-identified per EGA policy). The applicant agrees to:
  (a) use the data only for the stated research purpose,
  (b) not attempt to re-identify any participant,
  (c) destroy the data after the access period ends,
  (d) acknowledge the original [STUDY] publication in any derived work.

Signature:        ___________________________
Date:            ____ / ____ / ________
Institution:     ___________________________
```

---

## Pre-flight checklist before sending

Before submitting any of the above, make sure these are in order:

- [ ] **NIH eRA Commons account** active (for dbGaP) — register at https://commons.era.nih.gov/
- [ ] **GA4GH Passport** or **EGA CEGA account** active (for EGA) — register at https://ega-archive.org/
- [ ] **Local IRB** aware of the project (even if dataset is de-identified, an internal IRB notification is required by most institutions)
- [ ] **Data-security plan** signed by institution IT: where will the controlled-tier data live, who has access, how is access logged
- [ ] **One-line commitment** that no clinical decision will be informed by these datasets (this is the most common DAC objection — be explicit)
- [ ] **Citing the original paper** correctly (DACs reject applications that don't)
- [ ] **Acknowledgment plan** spelled out in the application (this is non-negotiable for most DACs)

---

## What to expect on the timeline

| Stage | Typical duration | Action |
|---|---|---|
| Submission → Initial DAC review | 2–4 weeks | EGA: usually auto-receipt; dbGaP: 1–2 weeks for completeness check |
| DAC review → Approval/Reject | 4–12 weeks | Yes — most DACs reject on the first pass. Be prepared to revise. |
| Approval → Data access granted | 1–4 weeks | Credentials + S3 / Aspera credentials issued |
| **Total** | **2–6 months** | Often longer on the first try |

---

## Pro tip from NEXT_STEPS.md §3

> **Apply for ≥ 2 datasets in parallel** — turnaround is non-deterministic and
> the worst-case for one dataset is often not the worst-case for the other.
> Apply for Cristiano 2019 (DELFI) + Newman 2016 (CAPP-Seq) on the same day;
> if both come through, you'll have two independent cohorts for validation.
> If only one comes through, you'll still have something to publish on.

---

*This document is the apply-today counterpart to `NEXT_STEPS.md §3` and
`docs/PRODUCTION_ROADMAP.md §4 Tier 2`. The templates above are reusable
across all open-access cfDNA WGS studies; edit the bracketed `[...]`
placeholders before sending. The goal is to start the data-access clock
today, not to wait until the rest of the readiness work matures.*