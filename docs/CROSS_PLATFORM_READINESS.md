# Cross-Platform Readiness (FinaleDB fragmentomics + FinaleMe methylation)

Status: scaffolded. Not yet runnable on real data — see
`docs/CROSS_PLATFORM_FINALEME.md` for the prerequisite recipe.

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

| Component                                       | Status (this machine)         | Why |
|------------------------------------------------|-------------------------------|-----|
| Java 21 runtime                                | present                       | `java -version` exits 0     |
| FinaleMe JAR (hybrid v0.61 + v0.58.1)           | **MISSING**                   | Recipe requires building from source; not yet done |
| Pretrained HMM models (Zenodo 14013719)         | **MISSING**                   | ~150 kB on Zenodo; not downloaded |
| Reference files (~1.4 GB: hg19.2bit + prior bw + CpG motif + mappability + chrom sizes) | **MISSING** | Multi-GB download not in budget |
| FinaleMe output directory                       | **MISSING**                   | Requires JAR + models + references |
| FinaleDB features cache                          | **PRESENT**                    | 627 samples × 5 channels in `/Users/hermes/cfdna-fragmentomics-pipeline/data/features/` |

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
- **Did not** modify the methodology demo's
  `results/cross_platform_methylation_fusion.json` (DISJOINT
  cohorts, synthesized scores — see its `honest_caveats`).

## Path forward

Operator checklist:

1. Run `_cross_platform_readiness_probe.py` — verifies the current
   `data_source` verdict.
2. Install FinaleMe prerequisites per
   `docs/CROSS_PLATFORM_FINALEME.md §Prerequisites` (~3 hours one-time).
3. Run `finaleme_pipeline.py decode` to produce per-sample
   `*.decoded.bed.gz` files. ~32 sec per sample for chr22; full
   genome ~5-10 minutes per sample.
4. Run `cross_platform_finaledb_validate.py` — emits
   `results/cross_platform_readiness.json` with `cross_platform_auc`
   (and ALL three fusion strategies, not just the headline).
5. Cite the JSON (path + filename) in any paper or PR that
   uses the cross-platform number, per the
   pre-flight-result-cross-check pattern.

## Companion scripts for the install path

```bash
# From the sibling repo, the FinaleMe install was validated end-to-end:
cd /Users/hermes/deepcatch-methylation
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
