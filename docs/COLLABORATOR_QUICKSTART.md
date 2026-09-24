# Collaborator Quickstart — 3 steps to contribute your cohort

> **One-page recipe.** Full spec at [`docs/COLLABORATOR_DATA_INTERFACE.md`](COLLABORATOR_DATA_INTERFACE.md). Errors? See the troubleshooting table in §8 of the full spec.

---

## Step 1 — make the directory

```bash
mkdir -p my_cohort/features
touch my_cohort/labels.tsv my_cohort/manifest.json
```

Result:

```
my_cohort/
├── features/        ← your 5 .npy files per sample land here
├── labels.tsv       ← sample IDs + disease class + cancer/healthy label
└── manifest.json    ← cohort metadata (lab, IRB, sample counts)
```

---

## Step 2 — copy your files

For each sample in your cohort, drop 5 `.npy` files into `my_cohort/features/`:

```
my_cohort/features/
├── H342.delfi_5mb_ratio.npy       (631 bins)
├── H342.delfi_5mb_coverage.npy    (631 bins)
├── H342.delfi_100kb_ratio.npy     (30894 bins)
├── H342.delfi_100kb_counts.npy    (30894 bins)
├── H342.fsd_histogram.npy         (196 bins, 5bp, normalized to sum=1)
├── H201.<same 5 channels>
└── ...
```

Then fill in `labels.tsv` (tab-separated, header row required):

```
sample_id	disease_class	label	study_id	publication_id	tissue_of_origin
H342	HCC	cancer	my_lab_2026		liver
H201	HCC	cancer	my_lab_2026		liver
H001	HEALTHY	healthy	my_lab_2026		
H002	HEALTHY	healthy	my_lab_2026		
```

And `manifest.json` (see §5 of the full spec for all fields):

```json
{
  "schema_version": "1.0",
  "cohort_name": "my_lab_2026",
  "contributing_lab": "Smith Lab",
  "contact_email": "lab-pi@example.org",
  "irb_number": "IRB-2026-001",
  "sample_count_cancer": 2,
  "sample_count_healthy": 2,
  "channels_present": ["delfi_5mb_ratio", "delfi_5mb_coverage",
                       "delfi_100kb_ratio", "delfi_100kb_counts",
                       "fsd_histogram"],
  "vector_lengths": {"delfi_5mb_ratio": 631, "delfi_5mb_coverage": 631,
                      "delfi_100kb_ratio": 30894, "delfi_100kb_counts": 30894,
                      "fsd_histogram": 196}
}
```

---

## Step 3 — run the adapter

```bash
# Sanity check (no files written)
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/adapter_local_cohort.py \
    --cohort-root my_cohort \
    --cohort-name my_lab_2026 \
    --dry-run

# Convert + write normalized cohort
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/adapter_local_cohort.py \
    --cohort-root my_cohort \
    --cohort-name my_lab_2026

# Optional: merge with existing open-data labels + run cross-study benchmark
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/adapter_local_cohort.py \
    --cohort-root my_cohort \
    --cohort-name my_lab_2026 \
    --merge-with /Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv \
    --cross-study
```

The output lives at `results/local_cohort/my_lab_2026/`:

```
results/local_cohort/my_lab_2026/
├── features/                  (symlinks to your original .npy files)
├── labels.tsv                 (canonical 6-column format)
└── adapter_manifest.json      (provenance: n_in, n_out, n_dropped, ...)
```

The adapter never mutates your original cohort. Cell-line samples (GM*, HeLa, HepG2, K562, …) are filtered automatically.

---

## Need help?

- **Spec:** [`docs/COLLABORATOR_DATA_INTERFACE.md`](COLLABORATOR_DATA_INTERFACE.md) — every field, every error
- **Reference cohort:** `data/synthetic_collaborator_cohort/` — exact layout you can copy
- **Reference implementation:** [cfdna-fragmentomics-pipeline](https://github.com/rollroyces/cfdna-fragmentomics-pipeline) — produces the 5 channels from raw `.bam`/`.frag.tsv.bgz`
- **Issue:** open one at [github.com/rollroyces/deepcatch/issues](https://github.com/rollroyces/deepcatch/issues) with title prefix `[cohort]`

---

**Honest scope:** this is a research benchmark, not a clinical validation. The pooled AUC is internal pooled-OOF on the same cohort that trained the model; for clinical-grade claims you need an independent held-out cohort at a clinical site with IRB approval. Per-study harmonization is applied at the benchmark level.