# The V3 methylation pipeline on real public data

**Status: the chain runs end-to-end on real FinaleDB fragments, and the
resulting channel has been measured. The measurement is a null: the
channel does not separate cancer from control beyond library complexity.
Sections 5 and 6 give the evidence; sections 2 and 3 give the two
upstream caveats that remain.**

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

## 5. The channel was measured, and it is a null

The pipeline was run over a **31-sample cohort** — FinaleDB publication 7
(Sun 2019), blood plasma, 17 colorectal cancer and 14 liver-transplant
controls, hg19 fragments, chr22 CpGs. The HMM was trained once on the
cohort and each sample decoded with that model.

```
python scripts/methylation_cohort_channel.py \
    --matrices-dir  ../methylation_cohort/matrices \
    --labels        ../methylation_cohort/cohort_labels.tsv \
    --work-dir      ../methylation_cohort/final_work \
    --out-json      results/methylation_cohort_channel.json
```

| Quantity | LOO AUC |
|---|---|
| decoded β summary (4 features) | **0.9444** |
| number of covered CpGs alone | **0.9444** |
| observations per CpG alone | **0.9444** |
| `beta_mean` alone, no model | 0.9530 |
| β summary, 1,550 CpGs covered by *both* arms | 0.8120 |
| β summary, depth-matched to 11,003 CpGs | 0.9444 (coverage-only 0.8034) |

Three different quantities — the decoded methylome, the number of CpGs
the library reached, and how many times each CpG was observed — give the
**identical** leave-one-out AUC of 0.9444. The channel's discrimination
is library complexity; the methylome adds nothing to it.

The matched-locus test is the one that cannot be explained away. Every
CpG is compared between the arms on the samples that actually cover it,
so composition cannot manufacture the difference:

```
n_cpgs_tested                                     103863
median_beta_difference_cancer_minus_control       0.0
fraction_of_cpgs_higher_in_cancer                 0.0245
permutation_p_value                               1.0
null_median_sd                                    0.0
```

The median cancer-minus-control difference is **exactly zero** over
103,863 CpGs, and the permutation null is also exactly zero, so the
observed value is unremarkable (p = 1.0). At the loci both arms cover,
the two groups carry the same methylation calls.

### 5.1 Why: the libraries do not carry the signal FinaleMe reads

FinaleMe predicts methylation from the *fragmentation pattern*: an
unmethylated CpG is cut out of the library, so it sits at a fragment
end; a methylated one is not cut, so it stays in the interior
(`FinaleMe.java` reads columns 4, 7, 9, 10 and 11 of the feature file —
`FragLen`, `Norm_Frag_cov`, `Offset_frag`, `Dist_frag_end` and the
prior — and never reads `methy_stat`).

Two things in the FinaleDB libraries destroy that information, and both
are measurable without running the decoder:

```
python scripts/methylation_library_check.py \
    --frag-files <finaledb.frag.tsv.bgz> <crag.frag.bed.gz> \
    --cpg-bedgraph CG_motif.hg19.chr22.bedgraph.gz \
    --out-json results/methylation_library_check.json
```

| | FinaleDB pub 7 (Sun 2019) | CRAG BH01 (deep MspI library) |
|---|---|---|
| chr22 fragments | 330,975 | 13,449,528 (**41x**) |
| median / mean length | 126 / 112 bp | 175 / 404 bp |
| longest fragment | **296 bp** | 30.6 Mb |
| fraction >= 300 bp | **0.000%** | 12.9% |
| fragment ends landing on a CpG | 1.11% / 1.03% (**0.68x / 0.63x chance**) | 2.34% / 1.32% (1.43x / 0.81x chance) |

1. **Size selection.** The FinaleDB library is trimmed at ~300 bp — not
   one fragment in 330,975 is longer. Its fragment ends are truncation
   points, not digestion sites, and they land on CpGs at *below* the
   chance rate. There is no methylation signature left to read.
2. **Depth.** A per-CpG state is only identified where a CpG is observed
   repeatedly. Across the cohort the libraries deliver **1.09
   observations per CpG** (range 1.00–1.22): the controls reach
   54,922–64,897 CpGs and the cancers 11,003–17,894, each CpG seen once.
   At one observation per CpG the decoder is extrapolating, not
   measuring, and what it returns is a restatement of how deep each
   library went.

The third measurement is the model-free one. A CpG that lies exactly at
a fragment end was cut out of the library, so the share of a CpG's
observations that sit at a fragment end estimates its unmethylated
probability without an HMM, a prior or any training.
`scripts/methylation_channel_diagnosis.py` computes it and pairs every
score with its coverage baseline.

## 6. What is still missing for a methylation *channel*

The extraction path runs on real public data and the channel has been
measured on 31 real samples. What remains is a library that can carry
the signal:

1. a deep, size-unselected MspI library — the CRAG/Jiang chr22 BAMs on
   Zenodo 6914806 are 41x deeper per chromosome than the FinaleDB
   libraries and keep their digestion boundaries; those samples are
   publication 6, i.e. the training cohort, so they would be a
   feasibility test and not an independent validation;
2. an *independent* MspI cohort for the actual validation, because no
   FinaleDB publication in the reachable corpus is both deep and
   untrimmed;
3. the GNN training run on the resulting real beta summaries.

Until then the honest statement is: *the methylation extraction path
runs on real public data; on the publicly available FinaleDB libraries
the channel measures library complexity and nothing else, and it is not
a validated feature.*
