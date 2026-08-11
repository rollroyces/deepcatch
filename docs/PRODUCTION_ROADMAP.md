# DeepCatch — Path to Production: Design & Validation Plan

**Status:** Research-stage simulation benchmark (honest). This document is the
concrete, staged plan to reach clinically meaningful deployment. It is written
to be executable, falsifiable, and honest about cost/time/risk at every step.

---

## 1. Where we stand (verified, 2026-08-10)

Real TCGA-LUAD mutations (20 patients, 5,738 mutations, GDC open access) as
ground truth; plasma reads simulated. **Panel-based (MRD-style) detection** is
the architecture that works:

| ctDNA fraction | Panel AUC | Sens @ 95% spec | Sens @ 99% spec |
|---|---|---|---|
| 1% | 1.000 | 1.000 | 1.000 |
| 0.5% | 1.000 | 1.000 | 1.000 |
| **0.1%** | **0.935** | **0.770** | **0.490** |

Per-position detection saturates at AUC 0.64 / sens 0.18 @ 0.1% ctDNA — it is
information-limited and should be retired as the headline metric.

Assay sweep (0.1% ctDNA, panel detection): duplex-UMI error suppression (~1e-4)
or ~50k× depth each reach **sens@95% = 1.000**. These are the production
specifications, not hopes.

**The single binding constraint from here on is real plasma cfDNA data.**
Everything below is organized around acquiring it and validating on it.

---

## 2. Product strategy: MRD first, MCED later

The clinical value of cfDNA detection splits into two products with very
different difficulty:

| | **MRD (tumor-informed)** | **MCED (screening)** |
|---|---|---|
| What | Track a patient's known mutations after treatment; detect recurrence | Detect any cancer in an asymptomatic person |
| Panel | Patient's own tumor mutations (what DeepCatch already simulates) | Fixed multi-cancer panel, methylation + fragmentomics |
| Ground truth | Surgical pathology, imaging follow-up | Longitudinal outcomes (years) |
| Typical ctDNA at detection | 0.01–1% | 0.001–0.1% |
| Specificity needed | 95–99% (adjuvant decisions) | 99.5%+ (population screening) |
| Evidence path | Retrospective cohort → prospective | NHS-Galleri-scale trials (140k+ patients) |
| Cost to validate | ~$1–5M, 2–3 years | $100M+, 5–10 years |

**Recommendation: DeepCatch should go to production as an MRD platform first.**
It is (a) scientifically what the panel detector already is, (b) the fastest
route to regulatory approval and clinical use, (c) the revenue that funds the
MCED program later. MCED is the long-horizon mission; MRD is the on-ramp.

---

## 3. Assay design (derived from the sweep)

| Parameter | Specification | Rationale |
|---|---|---|
| Panel | 100–500 loci per patient (tumor-informed; median LUAD ≈ 199 in our cohort) | Panel aggregation is what reaches 0.935 AUC @ 0.1% |
| Error suppression | **Duplex-UMI consensus** (target ≤1e-4) | Sweep: 1e-4 → sens@95% 1.000 @ 0.1% ctDNA |
| Depth | 20,000–50,000× on panel | Sweep: 50k× alone → 1.000 even at 2e-3 error |
| Input | 2 × 10 mL blood in cfDNA-stabilizing tubes (Streck/cfDNA BCT), ≤ 48 h processing | cfDNA yield & fragment integrity |
| Extraction | Column-based cfDNA kit + QC (yield, 167-bp peak, contamination) | Garbage in, garbage out |
| Sequencing | Illumina NovaSeq X / NextSeq 2k, PE150, 50M–100M reads/sample | Cost ~$500–1,500/sample |
| Matched WBC | Buffy coat / WBC DNA for **CHIP subtraction** (mandatory) | CHIP alone costs 5–10% specificity in >60 y |
| Multiplexing | Duplex barcodes (unique dual indexes) | Index hopping control |

The bench work is standard molecular biology; the differentiator is the
detector software (DeepCatch) and the validation rigor.

---

## 4. Data acquisition (the binding constraint — start now)

Tier 1 — **Already in hand: Jiang lab (CUHK) 129-sample 4-mer dataset.**
Finish the nested-CV re-estimate (data file needed at `data/deepcatch_data.xlsx`),
then hold out a strict external split. This is the fastest real-plasma result.
Coordinate with Prof. Jiang on publication/attribution terms.

Tier 2 — **Public WGS cfDNA datasets** (each needs an access application;
start the paperwork immediately — EGA/dbGaP review takes 2–6 months):
- Cristiano et al. 2019, *Nature* (DELFI; 545 samples incl. 215 cancer) — EGA: EGAS00001003828
- Mouliere et al. 2018, *Sci Transl Med* (fragment size profiles) — EGA
- Newman et al. 2016, *Nat Med* (CAPP-Seq NSCLC MRD, 40 patients) — EGA
- Abbosh et al. 2017 (TRACERx lung MRD, 24 patients, 96 serial samples) — EGA: EGAD00001002469
- GRAIL CCGA1/CCGA3 — by collaboration/agreement only
- TCGA does **not** contain plasma cfDNA (tissue only) — it defines panels,
  not plasma truth.

Tier 3 — **Own clinical cohort (the real product evidence):**
- Partner with 1–2 oncology centers (or Jiang lab) for a retrospective MRD
  cohort: resected NSCLC/CRC patients, serial post-op plasma draws, recurrence
  outcomes from imaging follow-up. 100–300 patients is a publishable,
  approvable start.
- IRB + data-use agreements + sample SOPs (2–4 months to first draws).

---

## 5. Validation ladder (each rung gates the next)

1. **Analytical validation** (CLIA-style, ~6 months once assay runs):
   - LoD dilution series: healthy plasma spiked with tumor DNA at 0.1%, 0.05%,
     0.01%, 0.005% ctDNA × 3 replicates × 3 runs → report LoD at 95% detection
   - Precision: intra-/inter-run, inter-operator, inter-site
   - Reproducibility: ≥95% concordance on replicate draws
   - Contamination/carryover, index-hopping audits, CHIP subtraction efficacy
2. **Clinical validity (MRD)** (~12 months):
   - Retrospective cohort: sensitivity/specificity of DeepCatch panel score vs
     recurrence at 3/6/12/24 months; lead-time vs imaging
   - Compare against published Signatera/CAPP-Seq numbers on the SAME cohort
     if possible (head-to-head is the gold-standard evidence)
   - Longitudinal model: the redesigned hierarchical-Bayes tracker (see §6)
3. **Clinical utility** (~24 months): prospective interventional study — MRD
   status guides adjuvant therapy escalation/de-escalation (cf. DYNAMIC,
   CIRCULATE trials). This is where "helps the world" is actually demonstrated.
4. **Regulatory**: CLIA/CAP lab accreditation first; then FDA De Novo or
   Breakthrough Device (or EU IVDR Class C/D) using the analytical + clinical
   validity package. Plan 12–24 months and $2–10M for the regulatory program.

---

## 6. Software productionization

1. **Pipeline**: wrap the detector in Nextflow/Snakemake + containers;
   pinned references (hg38, panel BED, PoN); mandatory QC gates; versioned
   outputs; one-command reproducibility (`make validate`).
2. **Calibration + monitoring**: per-run calibration (reliability diagram, ECE),
   batch-effect detection (control samples per run), drift monitoring in
   production; alerting when accuracy drifts.
3. **Longitudinal Stage 2 redesign (open work, don't ship the current CET)**:
   the honest baseline is AUC 0.49 — replace the single-patient VAF SPRT with a
   hierarchical-Bayes model across panel loci (cf. Setty et al. 2022), model
   CHIP trajectories as a distinct state, and require improvement over the
   current honest baseline before any claim.
4. **Tumor-agnostic path for MCED (later)**: fixed multi-cancer panel +
   methylation/fragmentomics features; cancer-type-first architecture
   (HCC-class assays are the easiest entry, per the Jiang results).
5. **Deployment**: HIPAA/GDPR-compliant processing; clinical report generation
   (already exists in `src/clinical`); audit logging; model versioning.

---

## 7. 12-month execution plan

| Months | Milestone | Exit criterion |
|---|---|---|
| 1–3 | Finalize Jiang nested-CV; submit EGA/dbGaP access applications (Cristiano, CAPP-Seq, TRACERx); draft assay SOP + IRB | Access approvals or rejection letters; SOP v1 |
| 4–6 | Run detector on 1–2 real WGS cfDNA datasets (fragmentomics + panel modules); build duplex-UMI demo pipeline; LoD dilution study (spiked plasma) | Real-plasma ROC ≥ simulation within stated gap; LoD table |
| 7–9 | Retrospective MRD cohort (100–300 patients); hierarchical-Bayes longitudinal redesign; CHIP filtering with matched WBC | MRD sens/spec with CIs; longitudinal AUC > 0.49 honest baseline |
| 10–12 | Analytical validation package; CLIA lab partnership; paper(s) + pre-print; MCED scoping doc | Validation report; submission-ready manuscript |

## 8. Risks (honest)

- **Simulation→reality gap is the #1 risk.** Every simulation number above will
  degrade on real plasma (the audit's degradation estimates were hand-waved —
  measure, don't assume). The design mitigates this by validating at each rung.
- **Non-shedders**: 10–30% of early cancers shed no detectable ctDNA — MRD
  sensitivity has a biological ceiling; report it, don't hide it.
- **CHIP**: without matched WBC, specificity collapses in older populations.
- **Tumor heterogeneity**: panel mutations can be lost after treatment —
  include subclonal/truncal prioritization in panel design.
- **Data access delays**: EGA/dbGaP routinely take 3–6 months; start now.
- **Regulatory cost/time**: real; the MRD-first strategy minimizes it.
- **Do not overclaim**: every public claim must trace to a computation in the
  repo (the project's existing honesty culture is a feature — keep it).

---

*Authored by Hermes Agent, 2026-08-10. All performance numbers trace to
`results/real_tcga_validation.json` (20 LUAD patients, 5 seeds) and
`real_tcga_validation.py` on branch `p0-fixes`.*
