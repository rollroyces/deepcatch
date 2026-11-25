# Cross-Platform Readiness (FinaleDB fragmentomics + FinaleMe methylation)

Status: scaffolded. The tool suite is installed and the probe passes;
what is missing is the decode input and its output, so the pipeline
still cannot print a cross-platform AUC. See
`docs/CROSS_PLATFORM_FINALEME.md` for the prerequisite recipe and
`results/cross_platform_readiness_probe.json` for the measured state.

## What the cross-platform pipeline is

Three scripts orchestrate a methylation channel ON TOP OF the
existing FinaleDB fragmentomics baseline (627-sample cohort, AUC
0.9782 ± 0.001 with LR-no-PCA C=1000):

```bash
scripts/_cross_platform_readiness_probe.py     # Always-safe diagnostic (exit 0)
scripts/finaleme_pipeline.py                  # FinaleMe Step 1+3 orchestrator (JAR / models probe)
scripts/finaleme_to_deepcatch_bridge.py        # FinaleMe per-sample β-values → DeepCatch JSON
scripts/cross_platform_finaledb_validate.py    # Cross-platform orchestrator
```

The orchestrator refuses to compute a `cross_platform_auc` when
either channel is missing — the JSON `data_source` field reads
one of `both`, `methylation_only`, `fragmentomics_only`, or `no_data`.

## Current state (FINALEME)

Measured by `scripts/_cross_platform_readiness_probe.py`, committed as
`results/cross_platform_readiness_probe.json`. The table this document
carried until 2026-11-25 said the JAR, the models and the reference
files were MISSING; they are not — they were installed on 2026-11-16 and
the table predated them. What is actually missing is the decoded
output, and the reason is the input, not the tool.

| Component                                       | Status (this machine)         | Why |
|------------------------------------------------|-------------------------------|-----|
| Java 21 runtime                                | present                       | `java -version` exits 0     |
| FinaleMe JAR (hybrid v0.58.1)                   | **present** (12.4 MB)        | `~/.hermes/.local/finaleme/FinaleMe-0.58.1-jar-with-dependencies.jar` |
| Pretrained HMM models                           | **present**                   | `~/.hermes/.local/finaleme/models` — FinaleMe's example cancer/healthy WGS HMMs; the per-chromosome model from the chr22 test run is alongside them |
| Reference files (~1.4 GB: hg19.2bit + prior bw + CpG motif + mappability + chrom sizes) | **present** | all five resolve in `~/.hermes/.local/finaleme/refs/hg19/` |
| FinaleMe output directory                       | **MISSING**                   | the directory exists but holds no decoded output; only a chr22 test run's `BH01.chr22.deep.decoded.txt.gz` was ever produced |
| FinaleDB features cache                          | **PRESENT**                    | 630 samples × 5 channels in `../cfdna-fragmentomics-pipeline/data/features/` |

**The one thing that is actually missing is the decode input.** The
training cohort's raw fragment files are not on this machine
(`../cfdna-fragmentomics-pipeline/data/raw/` is empty), and the decoded
per-CpG β-values that `results/methylation_cohort_channel.json` was
measured on (31 samples) are not either — they were an input to that
measurement, produced by an earlier run, exactly as the SAUNA tarball is
an input to its cohort. So the methylation channel's numbers are
reproducible given the decoded files, and the decoded files are not
versioned; the probe artifact records which of them exist right now.

Running the decode is: `scripts/finaleme_pipeline.py features --input
<sample.bam> --output-dir ~/.hermes/.local/finaleme/output`, with the
BAM built from fragments by `scripts/finaleme_bam_from_frag.py`. The
probe's "Next steps" prints the same command.

`data_source = fragmentomics_only` today. Cross-platform AUC is
**NOT** printed — the validator emits `cross_platform_auc: null`
with a refusal reason pointing at the missing FinaleMe suite.

## What the pipeline will emit when runnable

When the FinaleMe suite is installed, the JSON at
`results/cross_platform_readiness.json` carries:

- `methylation_channel_auc` (single-channel AUC from β-imputation),
- `fragmentomics_channel_auc` (single-channel AUC from the 5-channel
  fragmentomics baseline — identical to the published 0.9782 / 0.9732
  within noise),
- `fusion_strategies` (naive_average / logit_average / lr_fusion
  matching the methodology demo's
  `cross_platform_methylation_fusion.json`),
- `cross_platform_auc` (the best of the three fusion strategies,
  mean over 5 random seeds with 5-fold StratifiedKFold pooled OOF).

Honest upper-bound: methylation adds **≤ +0.005 AUC** over
fragmentomics alone on low-pass cfDNA (cfDNA biology is the
ceiling, not feature engineering). The methylation channel is
more useful for tissue-of-origin classification than
cancer-vs-healthy.

## What we deliberately DID NOT do

- **Did not** synthesize per-sample β-values from any published
  AUC. The bridge rejects `BetaValues.tsv` files whose β-values
  fall outside [0, 1] or whose sample IDs do not overlap with the
  FinaleDB features cache.
- **Did not** install the FinaleMe JAR or pretrained models.
  That's a ~3-hour one-time setup documented in
  `docs/CROSS_PLATFORM_FINALEME.md §Prerequisites`. The pipeline
  is honest about being scaffolded.
  **Superseded 2026-11-25:** the JAR, the models and the reference
  files ARE installed on this machine (see the measured table above);
  what is missing is the decoded output, because the training cohort's
  raw fragments are not here. The probe is the source of truth and it
  is committed as `results/cross_platform_readiness_probe.json`.
- **Did not** modify the methodology demo's
  `results/cross_platform_methylation_fusion.json` (DISJOINT
  cohorts, synthesized scores — see its `honest_caveats`).

## Path forward

Operator checklist:

1. Run `_cross_platform_readiness_probe.py` — verifies the current
   `data_source` verdict. ✅ done: it passes, and its output is
   committed as `results/cross_platform_readiness_probe.json`.
2. ~~Install FinaleMe prerequisites per
   `docs/CROSS_PLATFORM_FINALEME.md §Prerequisites`~~ — ✅ done on
   this machine (2026-11-16): the JAR, the models and the reference
   files are all present.
3. **Fetch the decode input**, which is the step that is actually open:
   the training cohort's raw fragment files are not on this machine
   (`../cfdna-fragmentomics-pipeline/data/raw/` is empty), and the
   decoded per-sample `*.decoded.bed.gz` / `BetaValues.tsv` files have
   to come from them. `scripts/finaleme_bam_from_frag.py` turns a
   fragment file into the BAM the decoder takes. ~32 sec per sample
   for chr22; full genome ~5-10 minutes per sample.
4. Run `cross_platform_finaledb_validate.py` — emits
   `results/cross_platform_readiness.json` with `cross_platform_auc`
   (and ALL three fusion strategies, not just the headline).
5. Cite the JSON (path + filename) in any paper or PR that
   uses the cross-platform number, per the
   pre-flight-result-cross-check pattern.

## Companion scripts for the install path

```bash
# From the sibling repo, the FinaleMe install was validated end-to-end:
cd deepcatch-methylation
ls scripts/download_hg19_parallel.sh
ls scripts/download_tcga_methylation.py
ls scripts/parse_finaleme_output.py
```

These cover the reference-file download (~90 sec parallel), the
TCGA methylation-array baseline (Phase 0 of the methylation
project plan), and the FinaleMe output parser — all validated on
this M4 Mac mini. The deepcatch-side bridge reuses that pre-flight
pattern and adds the per-sample validation that the missingness
state requires.

## Diagnostics

The probe is always safe — re-run it freely:

```bash
env -u PYTHONPATH ./.venv/bin/python scripts/_cross_platform_readiness_probe.py
```

It exits 0, prints the checklist, and is the canonical way to
check the operator-readable state of the inputs.
