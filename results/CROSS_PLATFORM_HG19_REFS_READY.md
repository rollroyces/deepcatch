# hg19 reference install — FinaleMe cross-platform pipeline

`generated_at: 2026-09-27T15:05:00Z` (corrected: 2026-09-27T15:15:00Z)

This change installs the hg19 reference files that `scripts/finaleme_pipeline.py decode`
and `scripts/cross_platform_finaledb_validate.py` need for the methylation channel.

## Files installed (~624 MB total)

Under `~/.hermes/.local/finaleme/refs/hg19/`:

| File | Size | Source |
|---|---|---|
| `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` | 147 MB | Zenodo 19392525 |
| `CpG_index.hg19.bed.gz` | 124 MB | Zenodo 19392525 |
| `CpG_index.hg19.bed.gz.csi` | 1 MB | Zenodo 19392525 |
| `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` | 309 MB | Zenodo 19392525 |
| `wgEncodeDuke...wgEncodeDacMapabilityConsensusExcludable.hg19.bed` | 86 KB | Zenodo 19392525 |
| `hg19.chrom.sizes` | 2 KB | UCSC goldenPath |

All sizes match Zenodo HEAD `Content-Length`. `hg19.2bit` was NOT downloaded:
the existing FinaleDB 5-channel features cache does not need raw reads, and the
probe explicitly notes it is only used for chrom.sizes extraction (which UCSC
already covers).

## What the methylation channel still needs

The methylation channel requires **per-sample FinaleMe decoded output** in
`~/.hermes/.local/finaleme/output/` (e.g. `BetaValues.tsv`). To produce this
requires running FinaleMe on the fragment BEDs of each FinaleDB sample.
Estimated runtime: 5+ days for the full cohort (630 samples).

Until those per-sample decodes exist, `cross_platform_finaledb_validate.py`
correctly refuses to print a `cross_platform_auc` even though all references
are installed — the methylation channel has nothing to read.

## Validate output (corrected)

```
data_source: fragmentomics_only
methylation_available: False
fragmentomics_available: True
n_overlap_samples: 0
cross_platform_auc: <refused — fragmentomics_only>
```

This is the **honest** state: methylation references are installed but no
decoded output exists yet, so the bridge has no real inputs and the validate
script refuses to fabricate a cross-platform number.

## Probe verdict (real)

```
[OK]    Java runtime            java version "21.0.12.1" 2026-08-18 LTS
[OK]    FinaleMe JAR            /Users/hermes/.hermes/.local/finaleme/FinaleMe-0.58.1-jar-with-dependencies.jar (12.4 MB)
[OK]    Pretrained HMM models   /Users/hermes/.hermes/.local/finaleme/models
[OK]    Reference files         all 5 reference files present
[MISS]  FinaleMe output         no BetaValues.tsv / *.decoded.bed.gz (decoded runs need real BAMs)
[OK]    FinaleDB features       630 samples × 5 channels (3564 total .npy files)
```

## Why an earlier version of this doc said "data_source: both"

The original `d9fc205` commit briefly produced a JSON with `data_source: 'both'`
+ `methylation_channel_auc: 1.0` driven by a 20-sample `synth_BetaValues.tsv`
in the output directory. This was a synthetic-fixture leak: the bridge found
that file by pattern-matching and treated it as real FinaleMe output. The
corrective patch (commit `<this commit>`) hardens `find_finaleme_per_sample`
to skip files prefixed with `synth_` and adds a regression test. The
synthetic file is deleted; `data_source` is now correctly `fragmentomics_only`
with a refusal reason.

The hard limit on cross-platform AUC remains: **real FinaleMe decoding
requires per-sample fragment BEDs that the FinaleDB features cache does not
expose.** Closing this gap is a 5-day batch decode job, not a code change.