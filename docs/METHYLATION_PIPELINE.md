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
  -> src/methylation_gnn/real_data.py      (real cohort -> the GNN's graph)
```

### The GNN branch now runs on real data

`src/methylation_gnn/` used to be architecture-only. Its graph builder
filled the node set with randomly generated regions when none were loaded
(`_ensure_regions`), and its data module filled the beta matrix with
random values, so the "methylation" the GNN consumed had no relationship
to any real sample.

`src/methylation_gnn/real_data.py` is the seam. It reads the same
decoded output the channel was measured on and returns the structure
`RegulatoryGraphBuilder.build_graph` expects, with the real CpG
coordinates as the graph's nodes:

```python
from src.methylation_gnn import (
    RegulatoryGraphBuilder, load_finaleme_cohort, load_real_regions,
    build_real_cohort_graphs,
)

cohort = load_finaleme_cohort("<final_work>", "<cohort_labels.tsv>")
builder = RegulatoryGraphBuilder(n_nodes=len(cohort["regions"]))
load_real_regions(builder, cohort)          # real loci, no synthetic fallback
graphs = build_real_cohort_graphs(builder, cohort)
```

On the published cohort this loads the **31 real samples** over the
**1,550 most widely covered loci** — the same locus pool the channel
measured on, because the loader applies the same selection rule (and
records which rule it used in `locus_selection`).

Two things are deliberately *not* faked:

* `cpg_density`, `gc_content` and `cpg_obs_exp` are absent from the
  FinaleMe output, so they are absent from the returned dict and the
  graph builder's documented zero-fill applies. The synthetic path used
  to fill them with random numbers, which would make a real-data run look
  better calibrated than it is.
* A run pointed at real loci can no longer fall back to synthetic
  regions: `require_real_regions` turns the old silent fallback into an
  error.

`test/test_methylation_gnn_real_data.py` covers the loader, the
no-fallback rule and the end-to-end path, and runs the real cohort when
the decoded output is present.

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
| β summary, 1,550 most widely covered loci | 0.803–0.829 across repeats — not reproducible (§5.4) |
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

The median difference and the permutation p-value are identical on every
run. `fraction_of_cpgs_higher_in_cancer` moves only in its third digit
(0.0245–0.0248 over six runs) for the reason given in §5.4.

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

### 5.2 The control: the same measurement on a library that carries the signal

A null is only interpretable if the measurement is shown to be able to
produce an effect, so the identical estimator was run on a genuine
MspI-digested library — CRAG BH01 (Zenodo 6914806), a healthy donor,
chr22, 13.4 M fragments — using the same pipeline and the same
5 Mb window of chr22:

| | FinaleDB pub 7 (control) | FinaleDB pub 7 (cancer) | **CRAG BH01 (deep MspI)** |
|---|---|---|---|
| observations per CpG | 1.00 | 1.20 | **21.40** (median 21, max 31) |
| CpGs with >= 2 observations | 26 | 392 | 158 (of the window so far) |
| ... of those, with a boundary fraction of exactly 0 | **100%** | **97.2%** | **88.6%** |
| ... with a graded, non-zero boundary fraction | 0% | 2.8% | **11.4%** (0.042–0.062) |

At one observation per CpG the estimate is either "interior" or
nothing: 97–100% of the CpGs that are seen at all carry a boundary
fraction of exactly zero, which is why the decoder returns hard calls
and why the per-CpG comparison between the arms has a median difference
of zero. At 21 observations per CpG the same estimator returns graded
values for one CpG in nine, and the decoder has something to work
with. The pipeline is therefore not broken; the FinaleDB libraries are
too shallow and too trimmed for the feature they are being asked to
carry.

The CRAG samples are publication 6 — the DeepCatch training cohort — so
this is a feasibility control, not an independent validation. No
independent, deep, untrimmed MspI cohort is reachable from the public
corpora surveyed in `docs/PUBLIC_DATA_SOURCES.md`.

### 5.3 Depth is the cause, on one library

The cleanest control uses one library and changes only the depth. The
same CRAG matrix, scored at every depth it happens to reach:

| observations per CpG | 1 | 2 | 3 | 5 | 8 | 10 | 15 | 21 |
|---|---|---|---|---|---|---|---|---|
| CpGs whose methylation estimate is graded | 0.8% | 1.6% | 3.1% | 3.7% | 5.7% | 8.1% | 11.1% | 11.5% |

"Graded" means the CpG's boundary fraction is not a constant zero, i.e.
the estimate carries information about whether that CpG was cut out of
the library. The resolution rises monotonically with depth, and at the
one-observation depth the FinaleDB libraries sit at, it is 0.8% — the
degenerate regime the channel is stuck in. Reproduce with:

```
python scripts/methylation_channel_diagnosis.py \
    --matrices-dir <window matrix dir> \
    --out-json results/methylation_channel_diagnosis.json
```

and read `depth.resolution_by_depth`.

### 5.4 One row of the §5 table is not reproducible — FinaleMe's HMM is not deterministic

The table in §5 is reproducible except for one row: the β summary
restricted to the 1,550 most widely covered loci (they are not all
covered by both arms, which is why the per-CpG test above is run on the
103,863 loci where both arms have coverage). Over five runs of the same
cohort that row came back as **0.8034 / 0.8120 / 0.812 / 0.8162 /
0.8291**.

The cause is upstream, in FinaleMe. On byte-identical decoder inputs
(the per-sample matrices and the WGBS prior), the **trained HMM model
differs between runs** — the decoded calls then flip for 16–20 of ~51,000
CpGs per sample, and the shared-locus AUC moves with them. This was
isolated by re-scoring the *same* decoded files under two different
`PYTHONHASHSEED` values: with the shared loci read in sorted order both
seeds return the identical value, so the scoring code is deterministic
and the HMM is not.

Everything else in §5 was identical on every run: the decoded-β AUC
(0.9444), the covered-CpG baseline (0.9444), the observations-per-CpG
baseline (0.9444), the depth-matched pair (0.9444 / 0.8034), and the
per-CpG median difference of exactly 0.0 over 103,863 CpGs with
permutation p = 1.0. The share of CpGs higher in cancer moves only in
its third digit (0.0245–0.0248). The null rests on the median
difference and the p-value, not on the shared-locus row — and the
shared-locus row does not carry the argument in either direction, since
every value in its range is far below the 0.9444 headline.

The artifact records this in its `reproducibility_note` field. A
secondary order-dependence in the scoring code (the shared loci were
read straight out of a `set`, so the float summation behind `np.mean`
depended on `PYTHONHASHSEED`) was fixed at the same time;
`test/test_methylation_hashseed_stability.py` pins it.

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
