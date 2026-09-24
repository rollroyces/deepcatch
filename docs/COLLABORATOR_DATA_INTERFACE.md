# Collaborator Data Interface — Plasma Cohort Contribution Spec

> **TL;DR.** Drop your files into the directory layout below, run one Python command, and your cohort joins the cross-study cfDNA fragmentomics benchmark. No code changes required on your side; no institutional integration needed beyond signing the standard data-use agreement.

This document specifies the **directory layout** and **manifest schema** any clinical collaborator can produce to contribute their plasma cfDNA cohort to the DeepCatch cross-study benchmark. The conversion is a single command (`scripts/adapter_local_cohort.py`) and the resulting cohort is automatically compatible with the existing `scripts/cross_study_finallydb.py` benchmark.

The full stack is open-source (MIT), the adapter is documented, and every cohort we accept is acknowledged in `RESULTS.md` with citation + IRB reference.

For a 3-step recipe, see [`docs/COLLABORATOR_QUICKSTART.md`](COLLABORATOR_QUICKSTART.md). For help or to flag a layout we haven't accounted for, open an issue or email the contact below.

---

## 1. What we ask for

| Item | Required? | Notes |
|---|---|---|
| `.npy` files for the 5 standard channels | **Yes** | One file per sample × channel, see §3 |
| `labels.tsv` | **Yes** | Sample IDs + disease class + binary label, see §4 |
| `manifest.json` | **Yes** | Cohort metadata, see §5 |
| IRB / ethics approval number | **Yes** | In `manifest.json` |
| Contact email + contributing lab | **Yes** | In `manifest.json` (we'll send the benchmark report) |
| Per-sample raw `.bam` / `.frag.tsv.bgz` | **No** | Only the pre-extracted channels are required; raw data stays at your site |
| Per-patient clinical metadata | **No** | Optional — only if your IRB allows it |

We do NOT need: raw sequencing reads, alignment files, or anything beyond the 5 pre-extracted channels. If you can produce the 5 channels but not the raw reads, you're still eligible to contribute.

---

## 2. Directory layout

```
<cohort_root>/                          ← your choice of path
├── features/                           ← 5 .npy files per sample
│   ├── <sample_id>.delfi_5mb_ratio.npy
│   ├── <sample_id>.delfi_5mb_coverage.npy
│   ├── <sample_id>.delfi_100kb_ratio.npy
│   ├── <sample_id>.delfi_100kb_counts.npy
│   └── <sample_id>.fsd_histogram.npy
├── labels.tsv                          ← sample IDs + disease class + label
└── manifest.json                       ← cohort metadata
```

**`<sample_id>`** is a string you choose (e.g. `H342`, `PT-001-2026`, `CRC_2026_42`). It must be unique within your cohort and must not contain characters that confuse filenames (no `/`, `..`, spaces, or shell-special chars). We strongly recommend alphanumeric + underscore + hyphen.

The **`<cohort_root>`** directory name is arbitrary — name it after your lab, study, or IRB protocol. The adapter only cares about what's inside it.

---

## 3. Feature files (the 5 channels)

Each sample needs **5 `.npy` files**, one per channel. All five must have `dtype=float32` and the documented vector length.

| Channel filename suffix | Shape (vector length) | What it is |
|---|---|---|
| `.delfi_5mb_ratio.npy` | `(631,)` | DELFI short/long ratio per 5 Mb bin across the genome |
| `.delfi_5mb_coverage.npy` | `(631,)` | 5 Mb coverage, **median-normalized** per sample to median = 1.0 |
| `.delfi_100kb_ratio.npy` | `(30894,)` | DELFI short/long ratio per 100 kb bin |
| `.delfi_100kb_counts.npy` | `(30894,)` | 100 kb coverage, raw counts (NOT normalized — the adapter handles this) |
| `.fsd_histogram.npy` | `(196,)` | 5 bp-resolution fragment-length histogram, normalized to sum = 1.0 (range 20–1000 bp) |

### How to produce these from raw data

If you're starting from aligned `.bam` files or FinaleDB-style `.frag.tsv.bgz` files, the [`cfdna-fragmentomics-pipeline`](https://github.com/rollroyces/cfdna-fragmentomics-pipeline) reference implementation produces all five channels. The pipeline's `--out-dir` argument should point at a directory whose contents are the `<sample_id>.<channel>.npy` files you need — these are the same artifact shape we expect.

If your pipeline produces a single 5-channel `.npy` per sample, split it along the documented vector lengths (the per-channel vectors concatenate in this order: 5mb_ratio, 5mb_coverage, 100kb_ratio, 100kb_counts, fsd_histogram). The 5-channel concatenation is exactly what the pipeline's `load5()` loader produces.

### Why 5 channels and not 8?

The full pipeline also produces `5mb_meanlen`, `100kb_meanlen`, and a 256-dim 4-mer motif histogram per sample. Our honest ablation (5-seed paired t-test, 627 cross-study samples) showed these add **< +0.001 AUC** vs the 5-channel baseline — well below the +0.005 AUC threshold that earns a feature its keep. So we only require the 5 channels, but the adapter will accept the 8-channel cohort if you have it (TBD; v1 requires exactly 5).

---

## 4. `labels.tsv`

Tab-separated, one row per sample. The header row is required (it tells the adapter how to parse the columns).

### Columns

| Position | Name | Required? | Example values | Notes |
|---|---|---|---|---|
| 1 | `sample_id` | **Yes** | `H342`, `CRC_42` | Must match a `features/<sample_id>.*.npy` prefix |
| 2 | `disease_class` | **Yes** | `HCC`, `BRCA`, `LUAD`, `CRC`, `OV`, `PAAD`, `HEALTHY` | Free-text per-cancer label — used for per-cancer OvR sens@spec reporting |
| 3 | `label` | **Yes** | `cancer` or `healthy` | Binary label for cancer-vs-healthy AUC |
| 4 | `study_id` | Optional | `jiang_2015`, `cristiano_2019`, `my_lab_2026` | Short identifier for your cohort; becomes the harmonization grouping key |
| 5 | `publication_id` | Optional | `6` (FinaleDB pub id), or leave blank | Only needed if you're contributing samples from a previously-published cohort |
| 6 | `tissue_of_origin` | Optional | `breast`, `lung`, `colon` | Free-text; used for tissue-of-origin subgroup reporting |

A minimal valid `labels.tsv`:

```
sample_id	disease_class	label	study_id	publication_id	tissue_of_origin
H342	HCC	cancer	my_lab_2026		liver
H201	HCC	cancer	my_lab_2026		liver
CRC_001	CRC	cancer	my_lab_2026		colon
CRC_002	CRC	cancer	my_lab_2026		colon
H001	HEALTHY	healthy	my_lab_2026		
H002	HEALTHY	healthy	my_lab_2026		
```

### Critical rules

1. **`sample_id` in `labels.tsv` MUST match a `features/<sample_id>.*.npy` prefix exactly.** If the prefix doesn't exist, the adapter rejects the cohort with a clear error listing every missing file.
2. **`label` is case-insensitive but must be one of `cancer` or `healthy`.** Anything else triggers a schema error.
3. **`disease_class` is free text** — we use it to group samples for per-cancer OvR analyses. Choose any consistent labels (e.g. `HCC`, `CRC`, `BRCA`, `HEALTHY`); avoid mixing `CRC` and `Colorectal` in the same cohort — pick one and stick with it.
4. **Cell-line IDs are filtered out automatically.** Samples whose IDs match the cell-line regex (e.g. `GM1100`, `HeLa`, `HepG2`, `K562`) are dropped before the benchmark runs. If you have cell-line controls in your cohort, leave them in — the adapter will remove them.

---

## 5. `manifest.json`

A small JSON file recording the cohort metadata. Required fields:

```json
{
  "schema_version": "1.0",
  "cohort_name": "my_lab_2026_hcc",
  "contributing_lab": "Smith Lab, Anytown University",
  "contact_email": "lab-pi@example.org",
  "irb_number": "IRB-2026-001",
  "sample_count_cancer": 45,
  "sample_count_healthy": 30,
  "channels_present": [
    "delfi_5mb_ratio", "delfi_5mb_coverage",
    "delfi_100kb_ratio", "delfi_100kb_counts",
    "fsd_histogram"
  ],
  "vector_lengths": {
    "delfi_5mb_ratio": 631,
    "delfi_5mb_coverage": 631,
    "delfi_100kb_ratio": 30894,
    "delfi_100kb_counts": 30894,
    "fsd_histogram": 196
  },
  "generation_date": "2026-09-15",
  "scope": "Plasma cfDNA WGS, HCC vs healthy, IRB-approved prospective collection.",
  "citation": "Smith et al., 2026, J. Hepatol. (in prep).",
  "license": "Data use agreement; redistribution prohibited without permission."
}
```

| Field | Required? | Notes |
|---|---|---|
| `schema_version` | **Yes** | Always `"1.0"` for this revision |
| `cohort_name` | **Yes** | Short identifier; becomes the output subdirectory name |
| `contributing_lab` | **Yes** | Human-readable lab/center name |
| `contact_email` | **Yes** | Where we send the benchmark report |
| `irb_number` | **Yes** | Or "N/A (deidentified data, IRB exemption granted)" |
| `sample_count_cancer` | **Yes** | Must match the labels.tsv count |
| `sample_count_healthy` | **Yes** | Must match the labels.tsv count |
| `channels_present` | **Yes** | Subset of the 5 standard channels you produced |
| `vector_lengths` | **Yes** | Per-channel vector length; cross-checked against the .npy shapes |
| `generation_date` | Optional | ISO date; we use it for citation ordering |
| `scope` | Optional | 1–2 sentence description of what the cohort is |
| `citation` | Optional | How to cite the cohort (paper, accession, or "synthetic — no citation") |
| `license` | Optional | Defaults to "research use only" if blank |

---

## 6. Worked example — the synthetic reference cohort

A canonical reference cohort ships at `data/synthetic_collaborator_cohort/`. It is regenerated by:

```bash
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/_synthetic_collaborator_cohort.py \
    --out data/synthetic_collaborator_cohort \
    --cohort-name synthetic_reference_2026
```

Layout after generation:

```
data/synthetic_collaborator_cohort/
├── features/
│   ├── SYN_C001.delfi_5mb_ratio.npy
│   ├── SYN_C001.delfi_5mb_coverage.npy
│   ├── ...                                                (155 .npy files = 31 samples × 5 channels)
│   ├── GM1100.delfi_5mb_ratio.npy                          ← cell-line sample (will be dropped)
│   └── GM1100.fsd_histogram.npy
├── labels.tsv                                              (31 rows = 20 cancer + 10 healthy + 1 cell-line)
└── manifest.json
```

The synthetic cohort's values are NOT biologically valid signal — they are deterministic random arrays. A 5-fold CV on the synthetic cohort produces an honest AUC ~0.65 (above 0.5 chance, far below clinical-grade 0.97+). The cohort is for **layout testing** only — verifying the directory structure before you go to the trouble of producing real fragmentomics features.

---

## 7. Running the adapter

After producing `<cohort_root>/`, run the adapter from the deepcatch repo root:

```bash
# 1. Sanity check + JSON report (does not write anything)
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/adapter_local_cohort.py \
    --cohort-root /path/to/<cohort_root> \
    --cohort-name my_lab_2026 \
    --dry-run

# 2. Convert + write normalized cohort to results/local_cohort/my_lab_2026/
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/adapter_local_cohort.py \
    --cohort-root /path/to/<cohort_root> \
    --cohort-name my_lab_2026

# 3. (Optional) Merge with existing labels + run cross-study benchmark
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/adapter_local_cohort.py \
    --cohort-root /path/to/<cohort_root> \
    --cohort-name my_lab_2026 \
    --merge-with /Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv \
    --cross-study
```

After step 2, your normalized cohort lives at:

```
results/local_cohort/my_lab_2026/
├── features/                           (symlinks to your original .npy files)
├── labels.tsv                          (canonical 6-column format)
└── adapter_manifest.json               (provenance: n_in, n_out, n_dropped, ...)
```

The adapter does NOT mutate your original cohort — it symlinks (or copies) the `.npy` files and writes a fresh `labels.tsv` + `adapter_manifest.json`. Re-running is safe and idempotent.

---

## 8. Common errors and what they mean

| Error message | Cause | Fix |
|---|---|---|
| `cohort_root does not exist or is not a directory` | Path typo or wrong cwd | Check `--cohort-root` is an absolute path; the adapter requires `<cohort_root>/{features/,labels.tsv,manifest.json}` |
| `missing features/ subdirectory` | You forgot the `features/` dir | Create it; put all per-sample .npy files inside |
| `missing labels.tsv` | You forgot `labels.tsv` | Create it with the 6-column header from §4 |
| `missing manifest.json` | You forgot `manifest.json` | See §5 for the schema |
| `<sample_id>.fsd_histogram.npy (shape (50,), expected (196,))` | Wrong FSD vector length | Your FSD should be 196 bins (5bp, 20–1000bp), not summary stats |
| `manifest.json channels_present includes unknown channels` | Channel name typo | Use exactly: `delfi_5mb_ratio`, `delfi_5mb_coverage`, `delfi_100kb_ratio`, `delfi_100kb_counts`, `fsd_histogram` |
| `schema validation failed ... <N> problem(s)` | Vector lengths don't match across samples | All samples MUST have the same vector length per channel; check that you didn't trim any channels for outlier samples |
| `schema error: labels.tsv has no data rows` | Empty TSV | Add at least one sample row |
| `--merge-with: file not found` | Wrong path to existing labels | Use an absolute path to a `.tsv` (e.g. the cfdna-fragmentomics-pipeline `labels_multiclass.tsv`) |

If you hit an error not listed here, copy the error message into an issue and we'll update this table.

---

## 9. Honest scope

What this interface DOES:

- Let you contribute a cohort of N≥10 plasma samples without writing any Python.
- Run the cross-study benchmark on your cohort merged with the open-data cohorts (FinaleDB publications 6 + 8, 627 samples).
- Get back a JSON report with per-cancer sens@spec, pooled AUC, and a true-confound control.

What this interface DOES NOT:

- Constitute clinical validation. The pooled AUC is pooled-OOF on a research cohort, not a held-out clinical validation. For clinical-grade claims, you need an independent held-out cohort at a clinical site with IRB approval.
- Replace your own quality control. Garbage in → garbage out. Run your own per-sample QC before submitting.
- Provide a hosting service. We do not store your raw data. Only the pre-extracted 5-channel vectors go through the adapter; your `.bam` / `.frag.tsv.bgz` files stay at your site.

The schema validation is rigorous but **not magic** — it catches structural problems (wrong shapes, missing files, typo'd channels) but cannot catch biological batch effects (depth drift, batch prep variation, different fragmentomics pipelines). Per-study harmonization is applied at the benchmark level, but we strongly recommend you tell us about known batch effects in `manifest.json`'s `scope` field so we can flag them.

---

## 10. Contact

- **Open an issue:** [github.com/rollroyces/deepcatch/issues](https://github.com/rollroyces/deepcatch/issues) — title prefix `[cohort]`
- **Email:** see `CITATION.cff` at the repo root for the maintainer email
- **Reference implementation:** [cfdna-fragmentomics-pipeline](https://github.com/rollroyces/cfdna-fragmentomics-pipeline) — for producing the 5 channels from raw data
- **Cross-study benchmark:** [`scripts/cross_study_finallydb.py`](../scripts/cross_study_finallydb.py) — what runs on the merged cohort
- **Existing benchmarks:** [`docs/CROSS_STUDY_BENCHMARK.md`](CROSS_STUDY_BENCHMARK.md) — what the baseline numbers look like

---

## Appendix: Schema versions

| Version | Date | Notes |
|---|---|---|
| `1.0` | 2026-09 | Initial 5-channel layout. `delfi_5mb_ratio`, `delfi_5mb_coverage`, `delfi_100kb_ratio`, `delfi_100kb_counts`, `fsd_histogram`. |

Backwards compatibility: the adapter honors `schema_version` — if we ship a `2.0` schema with extra channels, the `1.0` interface will continue to work and the adapter will translate.