# hg19 reference install — FinaleMe cross-platform pipeline

`generated_at: 2026-09-27T15:05:00Z`

This change installs the hg19 reference files that `scripts/finaleme_pipeline.py decode`
and `scripts/cross_platform_finaledb_validate.py` need for the methylation channel.

## Files installed (~624 MB total)

Under `~/.hermes/.local/finaleme/refs/hg19/`, symlinked to `~/.hermes/.local/finaleme/`
where the probe + pipeline expect them:

| File | Size | Source |
|---|---|---|
| `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` | 155,137,808 | Zenodo 19392525 |
| `CpG_index.hg19.bed.gz` | 130,111,624 | Zenodo 19392525 |
| `CpG_index.hg19.bed.gz.csi` | 1,635,723 | Zenodo 19392525 |
| `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` | 324,295,195 | Zenodo 19392525 |
| `wgEncodeDuke...wgEncodeDacMapabilityConsensusExcludable.hg19.bed` | 85,795 | Zenodo 19392525 |
| `hg19.chrom.sizes` | 1,971 (93 chr) | UCSC goldenPath |

All sizes match Zenodo HEAD `Content-Length`. `hg19.2bit` was NOT downloaded:
the existing FinaleDB 5-channel features cache does not need raw reads, and the
probe explicitly notes it is only used for chrom.sizes extraction (which UCSC
already covers).

The methylation channel still needs **per-sample FinaleMe decoded output** in
`~/.hermes/.local/finaleme/output/` (e.g. `BetaValues.tsv` from
`scripts/finaleme_pipeline.py decode` on each fragment BED). For this run,
the validate runs in `--synthetic-fixture` mode against a 20-sample
synth_BetaValues.tsv covering 10 cancer + 10 healthy FinaleDB sample IDs
that already have a 5-channel features cache. The synthetic β-values carry
a 7pp cancer-vs-healthy gap (healthy μ=0.65, cancer μ=0.58, σ=0.10) which
matches the published cfDNA-WGBS hypomethylation direction.

## Validate output

```
data_source: both
methylation_available: True
fragmentomics_available: True
n_overlap_samples: 20
methylation_channel_auc: 1.0
is_synthetic_fixture: True
cross_platform_auc: None
```

- `data_source='both'` confirms BOTH channels are wired through the probe +
  bridge + validate pipeline. Previously it was `fragmentomics_only`
  (commit 506e14a) because the methylation channel could not find any
  decoded output.
- `methylation_channel_auc=1.0` is the synthetic-fixture consequence —
  low within-class std + a 7pp gap gives perfect separation. Real
  FinaleMe-decoded β-values will regress this toward the published
  cfDNA WL ceiling of ~0.85–0.91 auROC.
- `cross_platform_auc=null` is a pre-existing loader bug:
  `scripts/_finaledb_feature_loader.py` requires a `fsd_histogram.npy`
  channel that does not exist in the on-disk features dir (which has
  `delfi_*` + `motifs.npy` + `gc_corrected.npy` instead). This is
  orthogonal to the hg19 reference install and is left for a follow-up.

## Probe verdict

```
[OK]    Java runtime            java version "21.0.12.1" 2026-08-18 LTS
[OK]    FinaleMe JAR            /Users/hermes/.hermes/.local/finaleme/FinaleMe-0.58.1-jar-with-dependencies.jar (12.4 MB)
[OK]    Pretrained HMM models   /Users/hermes/.hermes/.local/finaleme/models
[OK]    Reference files         all 5 reference files present
[MISS]  FinaleMe output         no BetaValues.tsv / *.decoded.bed.gz (decoded runs need real BAMs)
[OK]    FinaleDB features       630 samples × 5 channels (3564 total .npy files)
```
