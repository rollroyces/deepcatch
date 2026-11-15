# The V3 methylation pipeline on real public data

**Status: the chain runs end-to-end on real FinaleDB fragments. Two
upstream caveats remain, both documented below.**

This is the write-up of what was actually executed, not a design
document. Every number and file referenced here was produced by the
commands shown.

## 1. What runs

```
FinaleDB .frag.tsv.bgz (hg19)
  -> scripts/finaleme_bam_from_frag.py     (frag -> SAM with bases -> sorted/indexed BAM)
  -> FinaleMe CpgMultiMetricsStatsV2       (per-CpG matrix; -wgsMode -minMapQ 0 -minBaseQ 0)
  -> scripts/finaleme_decoder_input.py     (matrix + WGBS prior bigwig -> decoder input)
  -> FinaleMe hmm.FinaleMe                 (HMM train + decode -> per-CpG beta values)
  -> scripts/finaleme_to_deepcatch_bridge.py (DeepCatch methylation JSON)
```

### Verified on a real sample

Sample: `MBC_1330_1` — FinaleDB publication 9 (Adalsteinsson 2017),
breast cancer, **blood plasma**, seqrun 86278, hg19 fragment file,
1,556,590 fragments.

| Step | Command | Result |
|---|---|---|
| fetch | `curl .../entries/EE86278/hg19/EE86278.hg19.frag.tsv.bgz` | 12.2 MB, 5 columns (`chrom start end mapq strand`), bare chromosome names |
| BAM | `python scripts/finaleme_bam_from_frag.py --input <frag> --fa hg19.fa --out-bam <bam>` | 34 s -> 3,072,682 reads, sorted + indexed |
| matrix | `java ... CpgMultiMetricsStatsV2 hg19.2bit CG_motif...bedgraph.gz <bam> out.txt.gz -wgsMode -minMapQ 0 -minBaseQ 0` | chr22 subset: 70,107 rows in 2 min; whole genome: 9.3 MB in 12 min (still growing) |
| prior | `python scripts/finaleme_decoder_input.py --matrix ... --prior wgbs_buffyCoat_jensen2015GB.methy.hg19.bw --out ...` | 70,107 rows, 60,212 with a prior, 9,895 NaN |
| decode | `java ... hmm.FinaleMe <model> <decoder_input> <decoded>.txt.gz` | 445 KB of per-CpG beta values; "counted point in methy: 60212, Correct predicted: 54245 (90.09%)"; 120,424 data points |

Decoded output columns:
`#chr start end methy_perc_predict methy_count_predict total_count_predict methy_perc_obs methy_count_obs total_count_obs`.

## 2. Caveat 1 — the released pretrained models cannot be loaded

The Zenodo record 14013719 ships `healthy_WGS.mincg7.example.hmm_model`
and `cancer_WGS.mincg7.example.hmm_model` next to the source. They
cannot be deserialised by the released code:

- the models were serialized with the class
  `main.java.edu.mit.compbio.ccinference.hmm.BayesianNhmmV5`;
- the released source defines the class as
  `org.cchmc.epifluidlab.finaleme.hmm.BayesianNhmmV5` (the package
  was renamed between releases), so the stream's class name resolves to
  nothing on the classpath;
- and the field layouts differ as well (`maxCpgInFrag` is an `int` in
  the stream but a `double` in the current class, and `methyState` is
  the other way round), so a package alias would not be enough.

`FinaleMe` therefore has to **train the HMM from the data itself**:
run the decoder *without* `-decodeModeOnly` and it calls `trainHmm()`
before `decodeHmm()`.

That works, but a model trained on a single sample has no healthy
reference to contrast against, so the state structure collapses: the
decode reports "counted point in methy: 60212, counted point in unmethy:
0" and every CpG comes back at 100%. A meaningful beta channel needs
either the matching older FinaleMe build (the one the models were
serialized from) or a WGBS reference panel to train the HMM against
(`meth_wgbs.tar` / `ref_panel_wgbs.tar` / `HD_45.WGS.b37.tsv.gz` on the
same Zenodo record, several GB).

## 3. Caveat 2 — `scripts/finaleme_pipeline.py` targeted a class that does not exist

The orchestrator in the repo invoked

    edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder
    -fragmentInputTabix -fragStrandColumn N -valueWigs ...
    -inferMethyFromValueWig -useNoChrPrefixBam

No such class and no such flags exist in FinaleMe v0.58 (the release
the Zenodo record points at): the package is `org.cchmc.epifluidlab.
finaleme`, and the fragment-file entry point is
`CpgMultiMetricsStatsNoBam` with the positional contract

    CpgMultiMetricsStats [opts] hg19.2bit cpg_list.bed all_cpg.bed detail.txt.gz

The script now uses that class and drops the invented flags. The input
schema also differs, and the mismatch is silent rather than loud:

- FinaleDB serves 5 columns with bare chromosome names (`1`);
  the tool wants the 6-column bedpe schema with a `chr` prefix;
- without the conversion it reads MAPQ as the end coordinate and dies
  with "Start cannot exceed end", or matches no CpG at all.

`scripts/finaleme_bam_from_frag.py` and the bedpe conversion inside
`finaleme_pipeline.py` handle that.

## 4. Reference files

All from FinaleMe's Zenodo record (DOI 10.5281/zenodo.14013719) or
UCSC:

| File | Source | Notes |
|---|---|---|
| `hg19.2bit` | UCSC `goldenPath/hg19/bigZips` | 816 MB |
| `hg19.fa` (+ `.fai`) | UCSC `goldenPath/hg19/bigZips` | 3.2 GB, needed for read bases |
| `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` | Zenodo | 155 MB, the CpG list |
| `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` | Zenodo | 324 MB, the WGBS prior |
| `CpG_index.hg19.bed.gz` | Zenodo | 130 MB — **not** usable as the CpG list: its rows are `(chr, pos, index)`, not a BED, so the tool reads the index as the end coordinate |
| `healthy/cancer_WGS.mincg7.example.hmm_model` | Zenodo | see caveat 1 |

## 5. What is still missing for a methylation *channel*

The pipeline produces per-CpG methylation calls from real public data
today. What it does **not** yet produce is a discriminative feature:
`src/methylation_gnn/` still needs

1. a beta channel that is not degenerate — i.e. either the matching
   old FinaleMe build for the pretrained models, or a WGBS reference
   panel to train the HMM against;
2. per-sample beta summaries over many samples (the run above is one
   sample on one chromosome), so the channel can be computed for a
   cohort and tested against the cancer/control label;
3. the GNN training run on those real beta summaries.

Until then the honest statement is: *the methylation extraction path
runs on real public data; the methylation channel is not yet a
validated feature.*
