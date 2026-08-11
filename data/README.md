# Data Acquisition Guide

This directory holds real data files needed to run DeepCatch validation.
Files are NOT committed to the repo (privacy/licensing).  Place each
dataset here and set `DEEPCATCH_DATA_DIR` or follow the instructions below.

---

## 1. Jiang 4-mer end-motif frequencies (CUHK, Table S1)

**File to place:** `deepcatch_data.xlsx`
**Source:** Prof. Jiang's lab, CUHK — 129 plasma DNA samples × 256 4-mer motifs
**Status:** Already analyzed (results in `results/prof_jiang_4mer_analysis/`)
**Pipeline:** `run_jiang_analysis.py -i data/deepcatch_data.xlsx`
            `scripts/run_jiang_pipeline.py` (full pipeline verison)
**Note:** This file was removed from git tracking for patient privacy
(commit 8c812c0). Contact Prof. Jiang to re-provision.

---

## 2. Real TCGA mutation MAFs (GDC open access)

**Files to place:** `validation/tcga/tcga_cache/*.maf.gz`
**Source:** GDC API — downloaded automatically by `real_tcga_validation.py`
when run with `--n-patients 20 --cancer-types LUAD`
**Alternative:** Place GDC MAF files manually, or run the downloader:
```bash
python real_tcga_validation.py --n-patients 20 --cancer-types LUAD
```
The pipeline will:
1. Look for `*.maf.gz` in `validation/tcga/tcga_cache/`
2. If none found, download from GDC open-access API (30 aliquot MAFs)
3. Save normalized copies as `gdc_TCGA-LUAD_*.maf.gz`
4. Subsequent runs are offline (cached MAF files are re-used)
**Licensing:** TCGA open-access data — cite TCGA publication guidelines.
30 aliquot MAFs are committed in `p0-fixes` branch for reproducibility.

---

## 3. FLARE fragmentomics (GSE317007)

**File:** Already downloaded to `/tmp/GSE317007_motifs.txt.gz`
**Source:** GEO GSE317007 — FLARE pipeline, 6 HNSCC patients × 2 time points
**Status:** Fragmentomics features validated; clinical metadata (response labels)
needed for full longitudinal validation
**To copy permanently:**
```bash
cp /tmp/GSE317007_motifs.txt.gz data/GSE317007_flare_motifs.txt.gz
```
**Contact for clinical metadata:**
Loris De Cecco (loris.dececco@istitutotumori.mi.it)
IRCSS Istituto Nazionale Tumori, Milan, Italy
See `FLARE_CONTACT_TEMPLATE.md` for email draft.

---

## 4. Public cfDNA datasets for future validation

| Dataset | Access point | Timeline | Action |
|---|---|---|---|
| Cristiano 2019 (DELFI) | dbGaP phs0034536 | 2-6 months | File DAR at dbGaP |
| CAPP-Seq NSCLC MRD | EGA / contact Newman lab | Variable | Email authors |
| PanSeer (Taizhou) | EGA DAC | 3-6 months | Requires DAC approval |
| GSE185307 (cfDNA meth) | GEO | Immediate | Download + check labels |
| TRACERx (lung MRD) | EGA: EGAD00001002469 | 3-6 months | File EGA application |

---

## 5. How to add a new dataset

1. Place the file in `data/` (or set `DEEPCATCH_DATA_DIR` to its location)
2. Add a section to this README documenting:
   - What the file is (format, columns, sample count)
   - Where it came from (GEO/GDC accession, lab, contact)
   - What license/terms govern its use
   - Which DeepCatch script loads it
3. Update `NEXT_STEPS.md` with the validation task it enables

Never commit raw patient data without explicit consent and a data-sharing
agreement in place. The repo's `.gitignore` already excludes `data/`.
