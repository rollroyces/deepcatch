# DeepCatch — Next Steps: From Research Simulation to Clinical Validation

**Status:** All P0 code fixes committed on `p0-fixes` branch (228/228 tests green).
This document is the prioritized, executable action plan.

---

## 0. Push the branch (TODAY — 5 minutes)

```bash
cd /Users/hermes/deepcatch
git push origin p0-fixes
```

If GitHub auth isn't set up:
```bash
gh auth login
# or: export GH_TOKEN=<your PAT with repo scope>
```

Then open a PR on github.com/rollroyces/deepcatch from `p0-fixes` → `main`.
The CI will run (core + DL + guard jobs) and confirm 228/228 green.

---

## 1. Re-run Jiang nested-CV (THIS WEEK — 1 hour once data is on disk)

**The single highest-value real-data validation you can do RIGHT NOW.**

What you need: Prof. Jiang's Table S1 xlsx (129 samples × 256 4-mer frequencies,
38 controls, 91 cancer across 6 types). This file was in the repo but removed
for privacy (commit `8c812c0`).

Steps:
```bash
# 1. Put the file at any of these locations
cp /path/to/Table_S1.xlsx ~/deepcatch/data/deepcatch_data.xlsx
# OR set the env var
export DEEPCATCH_DATA_DIR=/path/to/directory/containing/deepcatch_data.xlsx

# 2. Run the nested-CV Jiang pipeline (feature selection inside folds now!)
cd ~/deepcatch
source .venv/bin/activate
python run_jiang_analysis.py \
    -i data/deepcatch_data.xlsx \
    -o results/jiang_nested_cv/ \
    --top-k 50 --seed 42 --lr-C 10.0

# 3. Compare with the old AUC (0.9845). Expect 0.96-0.98 with nested CV.
#    The old number had feature selection leakage (MWU top-50 on full data).
#    The new number is HONEST.
```

Expected output: `results/jiang_nested_cv/summary_report.md` with per-cancer-type
AUCs, nested-CV CIs, and selection stability analysis.

If you don't have the xlsx handy, contact the Jiang lab — they already know
about DeepCatch (see `results/prof_jiang_4mer_analysis/summary_for_professor_jiang.md`).

---

## 2. Validate the FLARE dataset longitudinal tracking (2-4 weeks)

We downloaded **GSE317007** (FLARE pipeline, 6 HNSCC patients × 2 time points,
ONT sequencing, 256 4-mer end-motif frequencies). The data is at
`/tmp/GSE317007_motifs.txt.gz`. DeepCatch's fragmentomics module already
validates against it (CG-depletion, AT-enrichment pattern confirmed).

**What's needed for full validation:**
- Patient-to-timepoint mapping (which QA IDs are the same patient at C1D1 vs C5D1)
- Clinical response labels (RECIST: CR/PR/SD/PD)

**Action:** Email the FLARE authors (contact below) requesting the clinical
metadata. Template at `data/FLARE_CONTACT_TEMPLATE.md`.

Once you have the labels:
```bash
# Build a labels file (response 0/1 per sample)
cat > data/flare_labels.csv << EOF
sample,timepoint,response
QA08,baseline,1
QA14,C5D1,1
...
EOF

# Run fragmentomics monitoring analysis
python -c "
import sys; sys.path.insert(0,'.')
from scripts.run_jiang_pipeline import *
# Load FLARE data + labels → compute pre-post feature shifts
# → validate longitudinal tracking
"
```

---

## 3. Acquire real plasma cfDNA WGS data (2-6 months — start APPLICATIONS NOW)

These datasets have cancer + control labels AND raw cfDNA sequencing reads.
Each requires a Data Access Committee application (2-6 month turnaround).

| Dataset | Samples | Has Controls? | Access | Application |
|---|---|---|---|---|
| **Cristiano 2019 (DELFI)** | 545 (215 cancer, 330 healthy) | ✅ Yes | EGA: EGAS00001003828 | dbGaP: phs0034536 |
| **CAPP-Seq NSCLC MRD** (Newman 2016) | 40 patients, serial draws | N/A (MRD) | EGA | Contact authors |
| **PanSeer / Taizhou** | 123,115 enrolled | ✅ Yes | EGA DAC | Requires DAC approval |
| **Snyder et al. 2016** | Healthy cfDNA nucleosome maps | N/A | GEO: GSE71378 | Open access |
| **FLARE / GSE317007** | 6 pts × 2 time points | ❌ No | GEO | Open access ✅ Already downloaded |
| **GSE185307** (cfDNA meth) | 24 samples, cancer + controls? | Possibly | GEO | Check → download if accessible |

**Templates needed:**
- dbGaP Data Access Request (DAR) — use NIH's online system
- EGA Data Access Agreement — per-study, contact DAC
- Letter to dataset authors — see `data/CONTACT_TEMPLATES.md`

---

## 4. IRB + own clinical cohort (3-6 months to first samples)

The data that makes DeepCatch a product rather than a research project.

**Minimum viable clinical study:**
- Retrospective MRD cohort
- 100-300 resected NSCLC/CRC/HCC patients
- Serial post-op plasma draws (q3 months, 2 years follow-up)
- Recurrence outcomes from imaging (CT/MRI, RECIST)
- Matched WBC for CHIP subtraction

**Partners to approach:**
- Prof. Jiang's lab at CUHK (already collaborating) — HCC patients
- IRCSS Istituto Nazionale Tumori, Milan (FLARE authors) — HNSCC
- Any thoracic/GI oncology center with a biobank

**Template:** `data/IRB_STUDY_OUTLINE.md` (draft a 1-page study concept)

---

## 5. Assay development (parallel to clinical, 3-6 months)

The sweep tells you the production specs. Build it:

| Spec | Target | Current benchmark |
|---|---|---|
| Panel size | 100-500 loci (tumor-informed, per-patient) | 199 median (LUAD) ✅ |
| Error rate | ≤ 1e-4 (duplex UMI consensus) | Sim at 1e-4 → AUC 1.0 @ 0.1% |
| Depth | 20k-50k× on panel | Sim at 50k× → AUC 1.0 @ 0.1% |
| Matched WBC | Buffy coat WGS/WES (mandatory for CHIP) | Not yet modeled |
| Blood volume | 2 × 10 mL in Streck/cfDNA BCT tubes | — |
| Processing | ≤ 48h from draw to plasma isolation | — |

**Vendor contacts:** Twist Bioscience (custom panels), IDT (xGen cfDNA), Qiagen
(QIAamp cfDNA), New England Biolabs (NEBNext duplex UMI).

---

## 6. Longitudinal model redesign (Stage 2, 2-4 months)

The current CET honest baseline is AUC 0.49 — the longitudinal model needs
a redesign before it can claim a clinical benefit. The 2026 literature points to:

**Censored-Poisson Bayesian Latent-Growth Change-Point Detector**
("Seeing Below the Limit of Detection", 2026-06-10)
- Models measurements BELOW the assay LoD (flickering detects/non-detects)
- Jointly estimates tumor growth trajectory across serial draws
- Bayesian framework with change-point detection for emerging subclones
- Directly applicable to DeepCatch's Stage 2 — replaces the SPRT-based CET

**Implementation plan:**
1. Port the censored-Poisson model into `src/longitudinal/`
2. Calibrate with real serial cfDNA data (FLARE or own cohort)
3. Benchmark against the current honest baseline (AUC 0.49) — REQUIRE improvement
4. Do NOT make any performance claim until the benchmark is beaten

---

## 7. Papers & preprints (when you have at least one real-data result)

**Minimum publishable unit:**
"Panel-based ultra-sensitive MRD detection from cfDNA: a spike-in benchmark
and real-plasma fragmentomics validation"

Content:
- TCGA spike-in benchmark (panel AUC 0.92 @ 0.1% ctDNA, full sweep)
- Jiang 4-mer real plasma (HCC AUC 0.98 nested CV)
- FLARE fragmentomics cross-validation (independent dataset)
- Assay sweep → production design guidance (duplex UMI + 50k× depth)

**Target journals:** Bioinformatics, PLOS Computational Biology, BMC Genomics,
or JCO Clinical Cancer Informatics (if you add clinical validation).

---

## 8. Weekly operating rhythm

| Day | Action |
|---|---|
| Monday | Push branch, check CI, review open issues |
| Tuesday | Code: one improvement from the P1/P2 list |
| Wednesday | Data: contact one dataset author / apply for one access |
| Thursday | Analysis: re-run pipeline with latest data, update README |
| Friday | Review: verify all claims trace to computation; update PRODUCTION_ROADMAP |

---

## Quick reference: files created in this session

| File | Purpose |
|---|---|
| `p0-fixes` branch (5 commits) | All code fixes, panel detection, Fisher/Strand scoring, context-aware sim, CHIP wiring, nested CV |
| `review/agent_review_2026-08-10.md` | Full external review + fix log |
| `docs/PRODUCTION_ROADMAP.md` | Production validation plan with literature references |
| `results/real_tcga_validation.json` | Panel AUC 0.9215 verified (20 LUAD patients) |
| `/tmp/GSE317007_motifs.txt.gz` | Real cfDNA fragmentomics data (FLARE, downloaded & validated) |

---

*This document is meant to be checked off. Update it as you complete items.*
