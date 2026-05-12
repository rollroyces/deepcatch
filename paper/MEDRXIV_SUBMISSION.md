# DeepCatch — medRxiv Submission Guide

**Date:** 2026-05-12
**Version:** v2.0 (Preprint)
**Target:** medRxiv (https://www.medrxiv.org/)

---

## About medRxiv

medRxiv (pronounced "med-archive") is a free online archive for **unpublished preprints** in the medical, clinical, and related health sciences. Operated by openRxiv (non-profit), founded by Cold Spring Harbor Laboratory, Yale University, and BMJ.

### Key Facts
- **No peer review** — screened for scientific content, ethics, and health risk only
- **No formatting requirements** — any PDF is accepted
- **Free to submit** — no fees
- **Multiple license options** — CC BY, CC BY-NC, CC BY-ND, CC BY-NC-ND, CC0, or no reuse
- **DOI assigned immediately** upon posting
- **Screening time** — typically 2–5 business days (slightly longer than bioRxiv due to health-risk check)
- **Versioned** — revisions retain same DOI

### What medRxiv Accepts ✅
- **Original research articles** ← DeepCatch fits here
- Systematic reviews and meta-analyses
- Data articles
- **Methodological research/investigations** ← Also fits here
- Clinical research design protocols

### What medRxiv Declines ❌
- Narrative reviews and hypotheses
- Case reports
- Editorial/commentary/opinion articles
- Laboratory protocols not part of a research paper

---

## ⚠️ Critical: Simulation-Only Study — medRxiv Considerations

DeepCatch is a **fully simulation-based study** with no real patient data. This is important for medRxiv submission:

### Potential Screening Flags (Be Prepared)
1. **"Clinical" claims without clinical data** — Screening staff may flag claims that could be misconstrued as clinically validated. Your paper already handles this well (explicit caveats on every table/claim).

2. **Health risk screening** — medRxiv checks for material that "might pose a health risk." Since your paper clearly states "simulation only, not for clinical use," this shouldn't be an issue.

3. **Category fit** — Methodological research is explicitly accepted. Choose "Bioinformatics" or "Computational Biology" as primary category.

### What to Emphasize in Submission
- ✅ Paper explicitly states "simulation-based" in abstract
- ✅ Every comparison table notes "simulation" vs "clinical"
- ✅ No patient data, no PHI, no IRB needed
- ✅ Open-source, fully reproducible
- ✅ Methodological contribution (performance-weighted fusion, CET framework)

---

## Pre-Submission Checklist

### ✅ Paper Modifications Done (2026-05-12)

| Item | Status | Notes |
|------|--------|-------|
| Target updated to medRxiv | ✅ | In `.tex` header |
| Date updated to May 12, 2026 | ✅ | In `\date{}` |
| bioRxiv references changed to medRxiv | ✅ | In body text |
| Abstract finalized | ✅ | 248 words, ≤ 300 limit |
| Competing interests declared | ✅ | In `.tex` |
| Funding declared | ✅ | "No specific funding" |
| Data availability stated | ✅ | GitHub, COSMIC, TCGA |
| Author contributions listed | ✅ | In `.tex` |
| Ethics statement | ✅ | "No real patient data, no IRB required" |
| Clinical trial registration | N/A | Simulation study |

### ⚠️ Before You Submit — PDF Required

**You must compile the LaTeX to PDF on your local machine.** The sandbox doesn't have LaTeX.

```bash
cd deepcatch/paper

# 3-pass compilation for resolved cross-references
pdflatex deepcatch_final.tex
bibtex deepcatch_final
pdflatex deepcatch_final.tex
pdflatex deepcatch_final.tex

# Check for errors
grep -c "Error" deepcatch_final.log    # Must be 0
grep "LaTeX Warning" deepcatch_final.log  # Should be minimal

# Output
ls -lh deepcatch_final.pdf
```

---

## Submission Walkthrough (Step by Step)

### Step 1: Register at medRxiv

1. Go to **https://submit.medrxiv.org/**
2. Click **"Create Account"** (or sign in with ORCID)
3. Fill in your profile:
   - **Full name:** Royce Lam
   - **Email:** roycelam@umich.edu
   - **Institution:** Independent Researcher
   - **ORCID:** Add if you have one (recommended)
4. Verify your email

---

### Step 2: Start New Submission

Click **"Submit a New Manuscript"**

---

### Step 3: Upload Files

**Required:**
- `deepcatch_final.pdf` — the compiled PDF (single file, all figures/tables inline)

**Optional but recommended:**
- `deepcatch_references.bib` — bibliography source
- `supplementary.pdf` — supplementary materials (compile from `supplementary.tex`)

**Note:** LaTeX `.tex` source upload is optional for medRxiv (unlike arXiv). Upload PDF only is the standard path.

---

### Step 4: Enter Manuscript Metadata

#### Title
```
DeepCatch: Performance-Weighted Multi-Modal Fusion with Cumulative Evidence Tracking for Pan-Cancer Detection from Cell-Free DNA
```

#### Authors
Add in order:
1. **Royce Lam** — Independent Researcher — roycelam@umich.edu

(If there are additional contributors, add them here with affiliations)

#### Abstract
Copy the full abstract from the PDF (or use `abstract_final.txt`). It should auto-extract from the PDF.

#### Subject Categories
Select **primary + secondary** (medRxiv allows multiple):

**Primary:**
- **Oncology** or **Cancer Biology**

**Secondary (pick 1-2):**
- **Bioinformatics**
- **Health Informatics**
- **Diagnostics**

**Alternative approach:** If there's a "Methodology" / "Computational Methods" subcategory, select that.

#### Keywords (copy these exactly)
```
liquid biopsy, circulating tumor DNA, early cancer detection, multi-modal fusion, 
longitudinal screening, fragmentomics, cell-free DNA, simulation study, 
cumulative evidence tracking, pan-cancer
```

---

### Step 5: Ethics & Declarations

#### IRB / Ethics Approval
```
Not applicable. This is a simulation study using only publicly available data 
(COSMIC v99, TCGA PanCancer Atlas). No human subjects, patient data, or 
protected health information were used. No institutional review board 
approval was required.
```

#### Clinical Trial Registration
```
Not applicable. This is a computational simulation study, not a clinical trial.
```

#### Competing Interests
```
The authors declare no competing interests.
```

#### Funding Statement
```
This work received no specific funding.
```

#### Data Availability Statement
```
Source code: https://github.com/rollroyces/deepcatch (MIT license)
Public data: COSMIC v99 (https://cancer.sanger.ac.uk/cosmic), 
             TCGA PanCancer Atlas (https://gdc.cancer.gov)
All simulation parameters and synthetic data generators are included in the repository.
No real patient data were used. Full reproducibility via RUN_ALL.sh.
```

---

### Step 6: License Selection

medRxiv license options:

| License | Recommend? | Why |
|---------|------------|-----|
| **CC BY 4.0** | ✅ **Recommended** | Maximum openness, required by many funders, allows text mining |
| CC BY-NC 4.0 | ⚠️ Alternative | Restricts commercial use |
| CC BY-ND 4.0 | ❌ | No derivatives — limits reuse |
| CC BY-NC-ND 4.0 | ❌ | Most restrictive |
| CC0 | ⚠️ | Public domain — no attribution required |
| No reuse | ❌ | Defeats purpose of preprint |

**Recommendation:** **CC BY 4.0** — standard for open science, required by many journals (PLOS, eLife, etc.), and compatible with your MIT code license.

---

### Step 7: Review & Confirm

Before clicking submit, double-check:
- [ ] Author list is complete and correct
- [ ] Author order matches the manuscript
- [ ] All authors have consented to submission
- [ ] No patient-identifiable information in the PDF
- [ ] Abstract matches the PDF exactly
- [ ] Keywords are accurate
- [ ] "Simulation study" caveat is prominent in abstract
- [ ] License selection is correct

Click **Submit**.

---

### Step 8: Post-Submission

1. **Confirmation email** — received immediately with temporary manuscript ID
2. **Screening** — 2–5 business days
3. **Possible outcomes:**
   - **Approved** — DOI assigned, publicly posted
   - **Revision requested** — minor changes needed (e.g., clarify ethics statement)
   - **Declined** — if screening criteria not met (rare for well-prepared submissions)
4. **Once posted:**
   - DOI is active immediately
   - Indexed in Google Scholar, Europe PMC
   - Share the DOI link

---

## medRxiv Screening — What They Check

Based on medRxiv's published criteria, screening includes:

1. **Scientific content** — Is it research? ✅ (Yes)
2. **Health risk** — Could it cause harm if misinterpreted? ⚠️ (Mitigated by prominent "simulation only" disclaimers)
3. **Ethical compliance** — IRB stated (N/A for simulation) ✅
4. **Competing interests** — Declared ✅
5. **Plagiarism check** — Original work ✅
6. **Scope** — Within health sciences? ✅ (Cancer detection methodology)

**Risk assessment:** LOW rejection risk. The paper is well-prepared with honest limitations and clear "simulation-only" labeling throughout.

---

## After Posting — What to Do

1. **Update GitHub README** — add medRxiv DOI badge
   ```markdown
   [![medRxiv](https://img.shields.io/badge/medRxiv-10.1101%2FXXXXXX-blue)](https://www.medrxiv.org/content/10.1101/XXXXXX)
   ```

2. **Update the paper's \date{}** — replace placeholder with actual DOI

3. **Share on social media** — medRxiv provides altmetrics, download counts

4. **Enable/disable commenting** — medRxiv allows public comments (you can disable)

5. **Monitor downloads & altmetrics** — via the medRxiv article page

---

## Journal Compatibility

Most journals accept medRxiv preprints. Confirm before formal submission:

| Journal | Accepts medRxiv? | Notes |
|---------|------------------|-------|
| Bioinformatics (Oxford) | ✅ Yes | Preprint must be acknowledged |
| PLOS Computational Biology | ✅ Yes | Encourages preprint posting |
| Nature Cancer | ✅ Yes | No restriction |
| JAMA Oncology | ✅ Yes | Check current policy |
| Cancer Discovery | ✅ Yes | Verify policy |
| Lancet Oncology | ✅ Yes | No restriction |

---

## Quick Reference: Submission URL & Contact

- **Submit:** https://submit.medrxiv.org/
- **Support:** medrxiv@medrxiv.org
- **Your email:** roycelam@umich.edu
- **Your GitHub:** https://github.com/rollroyces/deepcatch

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| PDF too large | Compress figures; medRxiv limit is generous (~40 MB) |
| Abstract extraction failed | Manually paste abstract in the form |
| Category not found | Pick closest match; screening staff can re-categorize |
| "Already posted elsewhere" rejection | Ensure no prior preprint on bioRxiv/arXiv — if previously on bioRxiv, withdraw first |
| Screening delayed > 5 days | Email medrxiv@medrxiv.org with manuscript ID |
| Need to update after posting | Use "Submit a Revision" in Author Area |

---

*Created: 2026-05-12 — Ironman 🦾*
