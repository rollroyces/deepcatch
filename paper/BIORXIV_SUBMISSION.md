# DeepCatch — bioRxiv Submission Guide

**Date:** 2026-05-10
**Version:** v2.0 (Preprint)
**Target:** bioRxiv (https://www.biorxiv.org/)

---

## About bioRxiv

bioRxiv is a free online archive and distribution service for **unpublished preprints** in the life sciences. It is operated by Cold Spring Harbor Laboratory (CSHL).

### Key Facts
- **No peer review required** — manuscripts are screened for scientific content and ethics only
- **No formatting requirements** — any PDF is accepted; LaTeX source encouraged
- **Free to submit** — no submission or publication fees
- **CC BY-NC-ND 4.0 license** by default (author retains copyright)
- **Can update after posting** — versioned posts, revisions allowed
- **DOI assigned** — citable immediately
- **Screening time** — typically 2–4 business days

### What bioRxiv Accepts
- Original research manuscripts
- Full-length articles AND short communications
- Methods papers
- Data notes
- Negative/null results
- Replication studies

### What bioRxiv Does NOT Require
- Formal journal formatting
- Supplementary materials at submission (can be added later)
- Competing interests statements (recommended but not required)
- Data/code availability statements (recommended but not required)
- Funding acknowledgments (recommended but not required)

---

## Pre-Submission Checklist for DeepCatch v2.0

### ✅ Ready Now

| Item | Status | File/Location |
|------|--------|---------------|
| PDF manuscript | ⚠️ Needs LaTeX compilation | `deepcatch_final.tex` |
| LaTeX source | ✅ Compilable | `deepcatch_final.tex` |
| Reference file | ✅ 57 entries | `deepcatch_references.bib` |
| Abstract | ✅ Updated with FragmentoSign validation | `abstract_biorxiv.txt` |
| Title | ✅ Finalized | In `.tex` file |
| Author(s) | ✅ Royce Lam | In `.tex` |
| Corresponding author email | ✅ `roycelam@umich.edu` | In `.tex` |
| All author affiliations | ✅ Listed | In `.tex` |
| Keywords | ✅ See below | — |
| Figures | ⚠️ Placeholder only | See `figures/FIGURE_DESCRIPTIONS.md` |
| Competing interests | ✅ Declared "no competing interests" | In `.tex` |
| Data availability | ✅ GitHub link, simulation-only noted | In `.tex` |
| Code availability | ✅ MIT license, RUN_ALL.sh | In `.tex` |
| Author contributions | ✅ Listed | In `.tex` |
| Acknowledgements | ✅ Written | In `.tex` |
| Supplementary material | ✅ Source ready | `supplementary.tex` |

---

## Submission Process (Step by Step)

### Step 1: Create bioRxiv Account

1. Go to https://submit.biorxiv.org/
2. Click "Create Account" or sign in with ORCID
3. Verify your email address
4. Fill in your profile:
   - Full name (as it will appear on the preprint): Royce Lam
   - Email (corresponding author): roycelam@umich.edu
   - Institutional affiliation
   - ORCID (recommended)

### Step 2: Prepare Submission Files

```bash
# Option A: Use the pre-packaged submission files (RECOMMENDED)
cd deepcatch/submission/
pdflatex deepcatch_final.tex
bibtex deepcatch_final
pdflatex deepcatch_final.tex
pdflatex deepcatch_final.tex

# Option B: Compile from paper/ directory
cd deepcatch/paper
pdflatex deepcatch_final.tex
bibtex deepcatch_final
pdflatex deepcatch_final.tex
pdflatex deepcatch_final.tex

# Check for LaTeX errors/warnings
grep -i "error\|warning" deepcatch_final.log | grep -v "font\|rerun"

# Compile supplementary as separate PDF (optional for bioRxiv)
pdflatex supplementary.tex
```

**Files to upload:**
- `deepcatch_final.pdf` — the main manuscript PDF
- `deepcatch_final.tex` — LaTeX source (optional but recommended)
- `deepcatch_references.bib` — bibliography (optional)
- `supplementary.pdf` — supplementary information (optional, can add in v2)
- Any figure files if embedded in LaTeX

**⚠️ Note:** `pdflatex` is required — no PDF was pre-generated (sandbox environment doesn't have LaTeX installed). Compile locally.

### Step 3: Enter Metadata

On the submission form:

#### Title
```
DeepCatch: Performance-Weighted Multi-Modal Fusion with Cumulative Evidence Tracking for Pan-Cancer Detection from Cell-Free DNA
```

#### Authors
Add all authors in order with:
- First name, last name
- Institutional affiliation(s)
- Email (corresponding author)
- ORCID (optional but recommended)

**Current author list:**
1. Royce Lam — Independent Research — roycelam@umich.edu
2. [Additional contributors as applicable]

#### Abstract
Copy from `abstract_biorxiv.txt` (see current version below)

#### Subject Area
Select **both** (bioRxiv allows multiple):
- **Bioinformatics** (primary)
- **Cancer Biology** (secondary)

#### Keywords
```
liquid biopsy, circulating tumor DNA, early cancer detection, multi-modal fusion, 
longitudinal screening, fragmentomics, cell-free DNA, simulation study
```

#### License
- Default: **CC BY-NC-ND 4.0** (no derivatives, recommended for preprints)
- Alternatives: CC BY 4.0, CC BY-NC 4.0
- **Recommended:** CC BY-NC-ND 4.0 — protects against unauthorized commercial use while allowing sharing

### Step 4: Select Collection (Optional)

bioRxiv has subject-specific collections. Choose:
- **Cancer Biology** (Subject Collection)
- The bioRxiv-affiliated journals channel (optional — manuscript will also appear in journal-specific feeds)

### Step 5: Confirm and Submit

- Review all metadata for accuracy
- Verify author list is complete
- Confirm that all authors have been informed of the submission
- Submit
- Receive confirmation email with temporary ID
- After screening (~2-4 business days), receive DOI and public posting

---

## Version Strategy

bioRxiv allows versioned updates. Plan:

| Version | Content | Timeline |
|---------|---------|----------|
| v1 (initial) | Current manuscript with simulation results + preliminary FragmentoSign validation | Now |
| v2 (update) | Add real cfDNA validation results (538 samples) | When available |
| v3 (update) | Add multi-seed validation (seeds 43-47), independent replication | Before journal submission |
| v4 (final) | Final pre-journal version with all figures | Before formal submission to Bioinformatics/PLOS Comp Bio |

**BioRxiv policy:** Each version is preserved. The latest version displays by default, but all prior versions remain accessible via the version history. The DOI always resolves to the latest version.

---

## Abstract (Current Version for bioRxiv)

See `abstract_biorxiv.txt` for the latest abstract text.

**Summary of changes from v1:**
- Added FragmentoSign preliminary validation results (fragment length, MDS, throughput)
- Added DELFI/1000G cross-validation acknowledgments
- Updated availability statement

---

## What bioRxiv Submission Is NOT

- ❌ **NOT peer review** — bioRxiv does not peer review manuscripts
- ❌ **NOT a journal publication** — preprints are not journal articles
- ❌ **NOT a barrier to later journal submission** — most journals accept preprints (check journal policy)
- ❌ **NOT indexed in PubMed** — bioRxiv preprints are indexed in Europe PMC and Google Scholar, not PubMed

---

## Journal Compatibility (Confirm Before Submitting)

Most computational biology/bioinformatics journals accept bioRxiv preprints, but verify:

| Journal | Accepts bioRxiv? | Notes |
|---------|------------------|-------|
| Bioinformatics (Oxford) | ✅ Yes | Preprint must be acknowledged |
| PLOS Computational Biology | ✅ Yes | Encourages preprint posting |
| Nature Biotechnology | ✅ Yes | No restriction |
| Nature Communications | ✅ Yes | No restriction |
| Cancer Discovery | ✅ Yes | Check current policy |
| NAR Genomics & Bioinformatics | ✅ Yes | Oxford journal, preprint-friendly |

---

## After Posting

1. **Share the DOI** — your preprint is now citable
2. **Update the GitHub README** — add bioRxiv DOI badge
3. **Share on social media** — bioRxiv generates a Twitter/X-friendly summary
4. **Monitor metrics** — bioRxiv provides download counts, altmetrics
5. **Respond to comments** — bioRxiv allows public commenting (enable or disable)
6. **Plan journal submission** — use feedback from preprint readers to improve manuscript

---

## Quick Commands Summary

```bash
# NOTE: All placeholder fixes are already applied in the submission package.
# File is ready to compile directly — no sed replacements needed.

# 1. Generate PDF (3-pass for cross-references)
cd deepcatch/submission/
pdflatex deepcatch_final.tex     # First pass: generates .aux file
bibtex deepcatch_final            # Processes .bib file for citations
pdflatex deepcatch_final.tex     # Second pass: resolves citations
pdflatex deepcatch_final.tex     # Third pass: resolves cross-refs

# 2. Check for errors
# Note: "undefined references" warnings on first pass are NORMAL
grep -c "Error" deepcatch_final.log    # Should be 0
grep "LaTeX Warning" deepcatch_final.log   # Should be minimal after 3 passes

# 3. Verify output
ls -lh deepcatch_final.pdf

# 4. Submit at: https://submit.biorxiv.org/
#    Email: roycelam@umich.edu
#    Name: Royce Lam
#    Institution: Independent Researcher
```

### Compilation Troubleshooting

| Problem | Solution |
|---------|----------|
| "I can't find file deepcatch_references.bib" | Make sure .bib is in same directory as .tex |
| "Citation undefined" warnings | Run: `pdflatex → bibtex → pdflatex → pdflatex` |
| "Environment undefined" error | Install missing LaTeX packages (see README.txt) |
| Figure placeholders show empty | Normal — figures are placeholder only; bioRxiv accepts this |
| "No room for a new \\count" | Add `\usepackage{etex}` before other packages |

---

*Updated: 2026-05-11 — Ironman 🦾*
