# AlphaGenome Atlas — Proxy Run Results (NOT REAL ALPHAGENOME)

> **Bottom line first: this is a deterministic-proxy run, not a real
> AlphaGenome run. The numbers below cannot tell us whether real
> AlphaGenome AVI scores would beat the existing CADD Top-K=200
> baseline. They tell us only that the weighted-aggregator code path
> runs end-to-end on this machine with a crude per-variant-class
> stand-in.**

> **Three-fix patch (this revision):** indels recovered (992 / 124,841
> now kept instead of dropped), `panel_llr_topk_by_avi_{K}` methods
> added, and tests updated. See "Three-fix patch" section below.

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

## Three-fix patch (2026-09-19, this revision)

The original proxy run had three honest issues, fixed in this
revision. Full design is in [`alphagenome_fix_proposal.md`](../alphagenome_fix_proposal.md).

### Fix 1 — keep indels (was: silently dropped)

`alphagenome_weights.load_tcga_mutations_with_ref_alt` used to drop
any variant with `len(ref) != 1 or len(alt) != 1` because the Atlas
SNV table only covers SNVs. That silently discarded 992 of 124,841
mutations (~0.79%, mostly `Frame_Shift_Del` / `Frame_Shift_Ins` /
`In_Frame_Del` / `In_Frame_Ins`).

After the patch, indels are kept. They receive deterministic proxy
weights via `_PROXY_BY_CLASS` (Frame_Shift_Del=29.5, etc.) and are
clearly tagged `source="proxy"` in the AVIWeight container — same as
every other proxy-weighted variant class. The proxy is honest because
it is clearly flagged; silent deletion is not.

Verification:

```text
Total: 124841
SNV:   123849
Indel:    992
```

### Fix 2 — top-K-by-AVI selector (`panel_llr_topk_by_avi_{K}`)

The original `panel_llr_topk_{K}` method selects the top-K positions
by `weight * LLR` ranking — a joint criterion that conflates AVI
prior with LLR magnitude. Fix 2 adds a new selector that picks the
top-K positions by **raw AVI alone**, then sums their weighted LLR:

```python
def _panel_score_with_weights_by_avi(per_pos_llr, avi_norm, weights, top_k):
    idx = np.argsort(avi_norm)[::-1][:top_k]
    return float((weights[idx] * per_pos_llr[idx]).sum())
```

This is the "panel prioritisation by prior" framing — closer to the
published CADD Top-K=20 framing, where the panel is *first* chosen
on biological importance and *then* scored. Both `topk_values` (500,
1000, 2000) are exposed.

### Fix 3 — tests

Added four tests to `test/test_alphagenome_weights.py`:

1. `test_load_tcga_mutations_keeps_indels` — indel count > 0 in loaded
   mutations, < 5% of total, and all indels have a non-empty
   `variant_class`.
2. `test_indel_keys_have_source_proxy_when_weighted` — indels in the
   proxy path come out as `source="proxy"` with positive scores (not
   silently dropped or `source="missing"`).
3. `test_panel_score_with_weights_by_avi_is_deterministic` —
   `_panel_score_with_weights_by_avi` returns identical output for
   identical input, and raises `ValueError` on length mismatch.
4. `test_run_comparison_includes_topk_by_avi_methods` — the new
   methods family appears in `alphagenome_panel_run.run_comparison`
   alongside the existing `panel_llr_topk_{K}` family.

The pre-existing `test_results_json_schema_and_proxy_flag` test was
extended to expect all 8 method names in the JSON output.

## Honest numbers (proxy run, 20 LUAD patients, 5 seeds, 5 TFs)

After the three-fix patch, the methods compared are 8 (was 5):

| Method                     | AUC @ 10% | AUC @ 5% | AUC @ 1% | AUC @ 0.5% | **AUC @ 0.1%**     | Δ vs uniform @ 0.1% |
|----------------------------|-----------|----------|----------|------------|--------------------|---------------------|
| panel_llr_uniform          | 1.0000    | 1.0000   | 1.0000   | 1.0000     | 0.9890 ± 0.0045    | —                   |
| panel_llr_avi              | 1.0000    | 1.0000   | 1.0000   | 1.0000     | 0.9895 ± 0.0037    | +0.0005             |
| panel_llr_topk_500         | 1.0000    | 1.0000   | 1.0000   | 1.0000     | 0.9930 ± 0.0054    | +0.0040             |
| panel_llr_topk_1000        | 1.0000    | 1.0000   | 1.0000   | 1.0000     | 0.9895 ± 0.0037    | +0.0005             |
| panel_llr_topk_2000        | 1.0000    | 1.0000   | 1.0000   | 1.0000     | 0.9895 ± 0.0037    | +0.0005             |
| **panel_llr_topk_by_avi_500**  | 1.0000 | 1.0000 | 1.0000 | 1.0000     | **1.0000 ± 0.0000**| **+0.0110**     |
| **panel_llr_topk_by_avi_1000** | 1.0000 | 1.0000 | 1.0000 | 1.0000     | **1.0000 ± 0.0000**| **+0.0110**     |
| panel_llr_topk_by_avi_2000 | 1.0000    | 1.0000   | 1.0000   | 1.0000     | 0.9895 ± 0.0037    | +0.0005             |

`avi_primary_source: "proxy"`, `avi_n_real: 0`, `avi_n_proxy: 114735`,
`n_mutations_with_ref_alt: 124841` (was 123849; +992 indels now
included), `pipeline_type: REAL_MUTATIONS_+_SIMULATED_PLASMA_READS`,
wall time ≈ 430 s. Full per-seed AUCs are in
[`results/alphagenome_weighted_llr.json`](../results/alphagenome_weighted_llr.json).

### `avi_n_proxy` is not 992 — what gives?

The fix proposal suggested `avi_n_proxy` would become exactly 992 after
Fix 1. That was a misreading: `avi_n_proxy` is the **total** count of
variants weighted by the proxy path (because no Atlas API or Tabix cache
is available), not just the indels recovered. In this run:

- `avi_n_proxy: 114735` — every variant the runner saw, because no
  Atlas API key is set and no Tabix cache exists.
- 992 of those 114,735 are the indels previously dropped.
- `avi_n_real: 0` — clearly flagged: not a single weight came from a
  real AlphaGenome source. This run is not evidence about real
  AlphaGenome performance.

### What this proxy run can — and cannot — tell us

- **Can**: confirm the weighted-aggregator code path runs end-to-end
  without regression below the 0.92 floor at 0.1% ctDNA. (All eight
  methods are at ≥ 0.989.)
- **Can**: confirm the proxy yields a stable, deterministic ranking
  (LoF > missense > silent) that, with the very coarse Top-K=500
  selection, gives a tiny lift over uniform at 0.1%.
- **Can**: report the comparison of the two top-K selectors. Top-K by
  raw AVI (`panel_llr_topk_by_avi_{500,1000}`) reaches AUC = 1.0000 ±
  0.0000 at 0.1% ctDNA; top-K by `w*LLR` (`panel_llr_topk_{500}`)
  reaches only 0.9930 ± 0.0054. The two methods are NOT equivalent;
  decoupling prior from likelihood matters even on this proxy path.
- **Cannot**: tell us whether **real AlphaGenome AVI** would beat the
  existing CADD Top-K=200 baseline. The existing CADD Top-K=200 result
  was +0.057 AUC over uniform on a separate cohort; that comparison
  used real CADD scores, not a proxy. We do not have any real AVI
  number to compare against here.

### Bottom-line comparison: top-K-by-w\*LLR vs top-K-by-AVI

The two selectors are different ranking rules. For each `K` we report
Δ over `panel_llr_uniform` at 0.1% ctDNA:

| K    | Top-K by `w*LLR` (Δ AUC) | Top-K by raw AVI (Δ AUC) | Winner |
|------|--------------------------|--------------------------|--------|
| 500  | +0.0040 (0.9930 ± 0.0054) | **+0.0110** (1.0000 ± 0.0000) | **raw AVI** |
| 1000 | +0.0005 (0.9895 ± 0.0037) | **+0.0110** (1.0000 ± 0.0000) | **raw AVI** |
| 2000 | +0.0005 (0.9895 ± 0.0037) | +0.0005 (0.9895 ± 0.0037)    | tie    |

Honest reading: at K=500 and K=1000, top-K-by-raw-AVI wins clearly
(+0.0110 vs +0.0040 and +0.0005). At K=2000 both methods tie (the
panel is large enough that AVI ordering and `w*LLR` ordering produce
the same effective subset). This is a real signal — the two methods
are not equivalent on this cohort — but the result is on the proxy
path and may not survive a swap to real AVI.

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
| `alphagenome_panel_run.py` | 8-strategy panel-LLR comparison runner. |
| `test/test_alphagenome_weights.py` | 10 proxy-path tests + 1 skipped real-API test (was 6 + 1). |
| `results/alphagenome_weighted_llr.json` | Proxy-run artefact, `avi_primary_source: "proxy"`. |
| `results/alphagenome_run.log` | stdout from the proxy run. |
