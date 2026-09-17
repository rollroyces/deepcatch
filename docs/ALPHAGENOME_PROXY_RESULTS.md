# AlphaGenome Atlas — Proxy Run Results (NOT REAL ALPHAGENOME)

> **Bottom line first: this is a deterministic-proxy run, not a real
> AlphaGenome run. The numbers below cannot tell us whether real
> AlphaGenome AVI scores would beat the existing CADD Top-K=200
> baseline. They tell us only that the weighted-aggregator code path
> runs end-to-end on this machine with a crude per-variant-class
> stand-in.**

## Why this run exists

The AlphaGenome Atlas API was unreachable from this environment at the
time of the run (no `ALPHAGENOME_API_KEY` was set, and no Tabix cache
for `avi_snvs_*.tsv.gz` was present). To honour the [AlphaGenome Terms
of Use](../alphagenome_weights.py) (no training of other ML models on
AVI outputs; non-commercial use only) and to keep the experiment
reproducible, the runner (`alphagenome_panel_run.py`) falls back to a
**deterministic proxy** that is clearly flagged in the output JSON.

## What the proxy actually does

The proxy is **not** a fabricated AVI score. It is a fixed table
mapping TCGA `Variant_Classification` labels to Phred-like anchor
values. The values are an order-of-magnitude ranking that *mirrors the
biological ranking AVI is trained to recover*, not AVI itself:

| Variant class          | Proxy value |
|------------------------|-------------|
| Nonsense_Mutation      | 30.0        |
| Frame_Shift_Del        | 29.5        |
| Frame_Shift_Ins        | 29.0        |
| Splice_Site            | 27.0        |
| Nonstop_Mutation       | 26.0        |
| Translation_Start_Site | 25.0        |
| Missense_Mutation      | 18.0        |
| In_Frame_Del           | 16.0        |
| In_Frame_Ins           | 15.5        |
| Splice_Region          | 14.0        |
| 5'UTR                  | 10.0        |
| 3'UTR                  | 8.0         |
| Silent                 | 3.0         |
| Intron                 | 2.0         |
| IGR / RNA              | 1.0         |
| _default (unknown)_    | 5.0         |

These values are then min-max normalised to [0.05, 1.0] across the
cohort, mirroring CADD-style normalisation. Every entry the proxy
emits is tagged `source="proxy"` in `alphagenome_weights.AVIWeight`,
and the run-level metadata sets `avi_primary_source: "proxy"`.

## What the proxy is NOT

- **Not** an AlphaGenome AVI score for any specific variant.
- **Not** correlated with the actual AlphaGenome model's per-position
  outputs. The proxy discards all genomic-position information; it
  depends only on coarse variant class.
- **Not** a meaningful comparison to the existing
  `cadd_vs_alphamissense_topk_20patient.json` results, which used
  real AlphaMissense scores. A proxy-vs-real comparison would be
  apples-to-oranges.

## Honest numbers (proxy run, 20 LUAD patients, 5 seeds, 5 TFs)

| Method                | AUC @ 10% | AUC @ 5% | AUC @ 1% | AUC @ 0.5% | **AUC @ 0.1%** | Δ vs uniform @ 0.1% |
|-----------------------|-----------|----------|----------|------------|----------------|---------------------|
| panel_llr_uniform     | 1.0000    | 1.0000   | 1.0000   | 1.0000     | **0.9890 ± 0.0045** | — |
| panel_llr_avi         | 1.0000    | 1.0000   | 1.0000   | 1.0000     | **0.9895 ± 0.0037** | **+0.0005** |
| panel_llr_topk_500    | 1.0000    | 1.0000   | 1.0000   | 1.0000     | **0.9925 ± 0.0050** | **+0.0035** |
| panel_llr_topk_1000   | 1.0000    | 1.0000   | 1.0000   | 1.0000     | **0.9895 ± 0.0037** | **+0.0005** |
| panel_llr_topk_2000   | 1.0000    | 1.0000   | 1.0000   | 1.0000     | **0.9895 ± 0.0037** | **+0.0005** |

`avi_primary_source: "proxy"`, `avi_n_real: 0`, `avi_n_proxy: 113833`,
`pipeline_type: REAL_MUTATIONS_+_SIMULATED_PLASMA_READS`, wall time
≈ 444 s. Full per-seed AUCs are in
[`results/alphagenome_weighted_llr.json`](../results/alphagenome_weighted_llr.json).

### What this proxy run can — and cannot — tell us

- **Can**: confirm the weighted-aggregator code path runs end-to-end
  without regression below the 0.92 floor at 0.1% ctDNA. (All five
  methods are at ≥ 0.989.)
- **Can**: confirm the proxy yields a stable, deterministic ranking
  (LoF > missense > silent) that, with the very coarse Top-K=500
  selection, gives a tiny lift over uniform at 0.1%.
- **Cannot**: tell us whether **real AlphaGenome AVI** would beat the
  existing CADD Top-K=200 baseline. The existing CADD Top-K=200 result
  was +0.057 AUC over uniform on a separate cohort; that comparison
  used real CADD scores, not a proxy. We do not have any real AVI
  number to compare against here.

### Why the proxy Δ is small

This is expected and honest. The proxy depends only on coarse
`Variant_Classification`; it cannot distinguish e.g. a damaging
missense in a hotspot from a benign passenger at the same gene, which
is exactly the discrimination a real AlphaGenome model is built for.
Real AVI is expected (per the AlphaGenome preprint) to add a
non-trivial Δ over both uniform and CADD; this proxy will not.

## Recommendation

1. **Re-run with a real AlphaGenome API key when available.** Set
   `ALPHAGENOME_API_KEY` (free, non-commercial tier at
   <https://deepmind.google.com/science/alphagenome>) and rerun:
   ```bash
   env -u PYTHONPATH .venv/bin/python alphagenome_panel_run.py \
       --n-patients 50 --seeds 5 --api-key "$ALPHAGENOME_API_KEY"
   ```
   The runner will detect the key and switch to `avi_primary_source:
   "atlas_api"` automatically.
2. **Do not cite the proxy numbers above as evidence about real
   AlphaGenome performance.** They are presented here only to prove
   the runner works end-to-end and to give the next runner a sanity
   baseline.
3. **Compare against CADD Top-K=200 only on a real-AlphaGenome run.**
   The honest question is "does real AlphaGenome beat real CADD at
   Top-K=200?", and that experiment is still open.

## Files

| File | Purpose |
|------|---------|
| `alphagenome_weights.py` | AVI loader with 3-tier fallback (API → Tabix → proxy). ToS-compliance comments at top. |
| `alphagenome_panel_run.py` | 5-strategy panel-LLR comparison runner. |
| `test/test_alphagenome_weights.py` | 6 proxy-path tests + 1 skipped real-API test. |
| `results/alphagenome_weighted_llr.json` | Proxy-run artefact, `avi_primary_source: "proxy"`. |
| `results/alphagenome_run.log` | stdout from the proxy run. |
