# CADD Top-K=200 validation on FLARE / GSE317007 — honest no-data report

**Date:** 2026-09-17
**Branch / commit:** `main` @ `fd1c3e3`
**Author:** DeepCatch subagent (delegated task)
**Bottom line:** **Not validatable as shipped.** FLARE/GSE317007 does not contain per-sample somatic variant calls (VCF/MAF). CADD Top-K=200 panel-LLR requires per-sample truth sets and per-allele coverage that this cohort does not provide.

---

## 1. What was checked

Searched every plausible location for FLARE-related data:

| Path | Result |
|---|---|
| `/Users/hermes/deepcatch/validation/flare/` | does not exist |
| `/Users/hermes/deepcatch/data/flare/` | does not exist |
| `/Users/hermes/deepcatch/data/` | contains `FLARE_CONTACT_TEMPLATE.md` and `JIANG_DATA_REQUEST.md` only — no FLARE matrices on disk |
| `/Users/hermes/deepcatch/data/GSE317007_*` | not present |
| `/tmp/GSE317007*` | not present (per `data/README.md` the motif matrix once lived at `/tmp/GSE317007_motifs.txt.gz`; that path is now empty) |
| `find . -iname "*flare*"` | only the metadata template, no data files |
| `find . -iname "*.vcf*" -o -iname "*.maf*"` | only TCGA-LUAD MAFs (the training cohort); zero FLARE/HNSCC variant files anywhere |

The GEO record for `GSE317007` was inspected directly. The deposit ships exactly one processed file — a 12×256 normalized 5′-end motif matrix — plus an SRA pointer (BioProject `PRJNA1405652`) to the raw ONT reads. There is no per-sample somatic SNV file (no VCF, no MAF, no JSON of variant calls). The associated paper (Frontiers in Genetics 2026, PMID 41908150) describes the FLARE pipeline as **fragmentomics only**: ichorCNA CNV profiles (window-level, not per-locus), cfTools/CancerDetector methylation-based tumor fraction, and 5′-end motif profiles. SNV calling is not a stated output.

## 2. Why CADD Top-K=200 cannot run on this cohort

`scripts/cadd_per_subgroup_llr.py` (the per-subgroup CADD Top-K=200 module) needs three things per sample:

1. **Truth set** — a list of `(chrom, pos, ref, alt)` somatic mutations labeled by patient. We have this for the 20 TCGA-LUAD patients (5,738 mutations from `validation/tcga/tcga_cache/`); we do **not** have it for any of the 12 FLARE samples.
2. **Per-allele coverage** at each mutation site — typically produced by simulating cfDNA reads from the truth set, or by pulling real coverage from a BAM. We have neither FLARE BAMs nor FLARE FASTQs locally.
3. **A binary cancer-vs-control label** to compute AUC and Sens@spec. All 12 FLARE samples are HNSCC (6 patients × baseline + on-treatment). There are zero healthy controls in GSE317007.

Without any of those, there is nothing to compute a CADD panel-LLR over.

## 3. Alternatives considered and rejected

- **Use CNV segments from FLARE as "mutation truth."** Rejected: CADD scores per-locus SNVs. CNV segments are window-level copy-number states with no CADD PHRED equivalent, and panel-LLR over CNV windows is a model DeepCatch does not implement.
- **Use ichorCNA tumor fraction per sample as a surrogate.** Rejected: there are no healthy controls in GSE317007, so AUC and Sens@spec cannot be computed.
- **Download SRA raw reads and run a Nanopore somatic caller (ClairS / Longshot) ourselves.** Technically possible but out of scope for a single subagent session: SRA fetch of 12 ONT cfDNA runs is multi-hour, multi-100-GB work; tumor-only somatic calling on cfDNA without matched germline is unreliable; and we would still need a CADD-scored truth set to evaluate the Top-K=200 lift. Honest no-go here. Logged as a future-engineering option in the JSON output.
- **Treat HNSCC as "similar to LUAD" and reuse the 5,738 TCGA-LUAD mutations against FLARE coverage.** Rejected: no FLARE coverage exists locally; biologically, a LUAD truth set on HNSCC cfDNA would be a category error (driver overlap is small).

## 4. Does the +X pp lift from TCGA-LUAD transfer to ONT-sequenced HNSCC samples?

**Cannot be determined from GSE317007 as published.** The per-subgroup Top-K=200 lift was measured on simulated cfDNA from TCGA-LUAD ground truth (real mutations, simulated reads, real CADD scores). Transferring that lift to ONT-sequenced HNSCC cfDNA would require:

1. A per-sample HNSCC truth set of somatic SNV positions (none shipped).
2. Per-allele coverage at those positions in real ONT cfDNA BAMs (no BAMs locally).
3. HNSCC-matched healthy controls for specificity (no controls in GSE317007).

None of those three are satisfiable from the data on disk or from the GEO supplement.

## 5. What IS already validated against FLARE

DeepCatch's **4-mer end-motif fragmentomics** feature extraction (a different channel from CADD panel-LLR) was cross-validated on GSE317007 and reproduces the published CG-depletion / AT-enrichment signature. See `paper/PAPER.md` and `data/FLARE_CONTACT_TEMPLATE.md`. That is a fragmentomics result, not a CADD panel result, and it does not speak to the Top-K=200 question.

## 6. No metrics fabricated

This report contains **no AUC, no Sens@spec, no per-cancer lift**. All numbers in `results/cadd_flare_validation.json` are metadata (file paths, sample counts, accession IDs, paper reference). Computing fake numbers would violate the explicit honesty constraint on this task and on the DeepCatch project.

## 7. Recommended next steps (out of scope here)

A future, longer-horizon task could perform the actual cross-platform validation:

1. `prefetch` all 12 SRA runs from `PRJNA1405652`, align to hg38 with `minimap2 -ax map-ont`.
2. Run a Nanopore somatic SNV caller per tumor sample (ClairS-TO or Longshot; no matched germline, so expect low precision and treat results as an upper bound).
3. Use the union of called SNVs as the per-sample "truth," score with CADD, re-run `cadd_per_subgroup_llr.py` logic on those 12 samples × ~dozens-of-SNVs.
4. Compare to the TCGA-LUAD baseline; report whether the +X pp lift survives the platform + cancer-type change.

Realistically this is a multi-day engineering task, not a single subagent scope.

## 8. Files written this session

- `/Users/hermes/deepcatch/results/cadd_flare_validation.json` — structured honest-no-data record (paths checked, GEO inspection, alternatives evaluated, no fabricated metrics).
- `/Users/hermes/deepcatch/docs/CADD_FLARE_VALIDATION.md` — this file.

No source code, no test files, no CI-affecting paths were modified. Existing 51/51 test pass count on `main` is unchanged.