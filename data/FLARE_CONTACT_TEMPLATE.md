# FLARE Dataset — Contact Template

**Recipient:** Loris De Cecco (loris.dececco@istitutotumori.mi.it)
**Institution:** IRCSS Istituto Nazionale Tumori, Milan, Italy
**Dataset:** GSE317007 — FLARE: Long-Read Nanopore Fragmentomics for Liquid Biopsy

---

## Email Draft

> Subject: DeepCatch validation on FLARE dataset (GSE317007) — clinical metadata request

> Dear Dr. De Cecco,

> I'm a researcher working on DeepCatch (github.com/rollroyces/deepcatch), an
> open-source computational framework for multi-cancer early detection from
> cfDNA. We develop panel-based MRD detection and fragmentomics analysis tools,
> including 4-mer end-motif profiling similar to the Jiang et al. (2020) approach.

> I recently downloaded your FLARE supplementary data (GSE317007 —
> `GSE317007_normalized_motifs_matrix.txt.gz`) and validated our fragmentomics
> feature extraction pipeline against it. The 256 4-mer motif profiles from
> your 12 samples (CG-depletion and AT-enrichment patterns) are consistent with
> published pan-cancer cfDNA fragmentation signatures, confirming our pipeline
> works on independent real data.

> To extend this to a clinically meaningful validation, I would be grateful if
> you could share the patient-to-timepoint mapping for the 12 samples. Based on
> the GEO description, the study design appears to be:

> - 6 patients with recurrent/metastatic HNSCC on Nivolumab
> - 2 time points each: baseline (C1D1) and on-treatment (C5D1)

> If you can share which QA IDs correspond to which patient and time point
> (e.g., QA08 = Patient 1, C1D1; QA14 = Patient 1, C5D1), and any clinical
> response labels (RECIST: CR/PR/SD/PD or PFS), we could validate our
> longitudinal fragmentomics tracker — which specifically tests whether
> 4-mer profiles shift between baseline and treatment in a direction that
> correlates with clinical response.

> This would be properly cited in any resulting publication. I'm happy to
> share our analysis results with you before publication.

> Thank you for making the data open access — it's a valuable resource for
> the cfDNA fragmentomics community.

> Best regards,
> Royce
> DeepCatch Project
> github.com/rollroyces/deepcatch

---

## What to do with the response

If they provide the metadata, save it as `data/flare_metadata.csv`:
```csv
sample,patient,timepoint,response
QA08,1,baseline,PR
QA09,2,baseline,SD
...
QA14,1,C5D1,PR
QA19,6,C5D1,PD
```

Then run the longitudinal validation (add to `run_jiang_analysis.py` or a
dedicated script). The key analysis: pre-post change in MDS, GC bias, CG
ratio, and motif diversity — stratified by clinical response.
