# AlphaGenome code fixes — design proposal

**Date:** 2026-09-17
**Code under review:** `alphagenome_weights.py`, `alphagenome_panel_run.py`

## Issues identified

### Issue 1: Indels dropped (992 / 124,841 = 0.79%)

**Current code** (`alphagenome_weights.py:344`):
```python
if not ref or not alt or len(ref) != 1 or len(alt) != 1:
    # AVI Atlas SNV table covers SNVs only; skip indels here
    continue
```

**Impact:** 992 mutations (Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, etc.) silently dropped. ~2.9% of patients on average lose 1 indel.

**Why:** The Atlas SNV table is SNV-only. Indels have a separate AVIScore scoring track (different genomic coordinate space). The code conservatively drops them rather than mishandle.

**Fix options:**
- **A: Keep dropping** with explicit count + log message ("dropped 992 indels; no AVI source available for indels")
- **B: Use proxy weights for indels** (already in `_PROXY_BY_CLASS`: Frame_Shift_Del=29.5, In_Frame_Del=16.0, etc.). Indels contribute to LLR aggregation as proxy-weighted even when other mutations are real AVI-weighted.
- **C: Fetch AVI for indels** via separate Atlas endpoint (`requested_scorers=["AVIScore"]` requires Interval objects, not Variant — needs different code path).

**Recommendation:** **Option B.** Indels get proxy weights (same class anchors as today). The proxy is honest because it's clearly tagged `source="proxy"`. Real AVI is reserved for SNVs where it actually helps. Combined run will say "real=113833, proxy=992" instead of dropping indels entirely.

### Issue 2: Top-K uses largest LLR contributions, not highest-AVI

**Current code** (`alphagenome_panel_run.py:67`):
```python
if top_k is not None and top_k < len(contrib):
    idx = np.argsort(contrib)[::-1][:top_k]
    return float(contrib[idx].sum())
```

**Impact:** Top-K selects based on `weight * LLR` ranking. A high-LLR mutation with low weight still ranks high. A high-AVI mutation with low LLR might rank low (because LLR is the dominant factor at most positions).

**Why:** This is what the CADD Top-K=200 finding (commit `ec16e0d`) showed works best on the 20-patient cohort. The intuition: weighting + selection both look at "does the model think this position is informative?".

**But it conflates two things:** CADD/AVI priority (priors) vs LLR magnitude (likelihood evidence). Selecting top-K by `w * LLR` is a joint criterion.

**Fix options:**
- **A: Keep as-is** (current behavior). It's what gave us the +0.057 lift.
- **B: Add `panel_llr_topk_by_avi_{K}`** that selects top-K by raw AVI score alone, then sums their weighted LLR. This is "panel prioritization by prior" — closer to the published CADD Top-K=20 framing.
- **C: Add `panel_llr_topk_by_llr_{K}`** that selects top-K by raw LLR (uniform selection on likelihood), then weights.

**Recommendation:** **Add Option B as a new method.** The honest comparison is then:
- `panel_llr_uniform` (baseline)
- `panel_llr_avi` (weighted, full panel)
- `panel_llr_topk_{K}` (top-K by `w*LLR`, current)
- `panel_llr_topk_by_avi_{K}` (top-K by raw AVI, then weight × LLR) ← **NEW**

If Option B (AVI-only Top-K) wins, that's the published CADD Top-K=20 finding re-discovered with AlphaGenome. If current (w*LLR) wins, that's a new finding worth reporting.

### Issue 3: No fallback to top-K-by-AVI

This is the same as Issue 2. Adding Option B addresses both.

## Recommended plan

### Patch 1: `alphagenome_weights.py` — keep indels, give them proxy weights

**File:** `alphagenome_weights.py`, function `load_tcga_mutations_with_ref_alt`

**Before** (lines 344-346):
```python
if not ref or not alt or len(ref) != 1 or len(alt) != 1:
    # AVI Atlas SNV table covers SNVs only; skip indels here
    continue
```

**After:**
```python
if not ref or not alt:
    continue  # malformed allele string
# Note: indels (len != 1) are KEPT — they will receive proxy weights via
# `_PROXY_BY_CLASS` since the Atlas SNV table is SNV-only. Real AVI is
# only available for SNVs; indels still contribute via deterministic proxy.
```

### Patch 2: `alphagenome_panel_run.py` — add `panel_llr_topk_by_avi_K` method

**File:** `alphagenome_panel_run.py`, function `_panel_score_with_weights`

Add a new helper:
```python
def _panel_score_with_weights_by_avi(
    per_pos_llr: np.ndarray,
    avi_norm: np.ndarray,
    weights: np.ndarray,
    top_k: int,
) -> float:
    """Select top-K by raw AVI priority, then sum weighted LLR."""
    if len(avi_norm) != len(per_pos_llr):
        raise ValueError("avi_norm length mismatch")
    idx = np.argsort(avi_norm)[::-1][:top_k]
    return float((weights[idx] * per_pos_llr[idx]).sum())
```

Then add the method to the methods list:
```python
methods = ["panel_llr_uniform", "panel_llr_avi"]
methods += [f"panel_llr_topk_{k}" for k in topk_values]      # top-K by w*LLR (current)
methods += [f"panel_llr_topk_by_avi_{k}" for k in topk_values]  # top-K by AVI (NEW)
```

And the scoring branch:
```python
elif method.startswith("panel_llr_topk_by_avi_"):
    k = int(method.split("_")[-1])
    sp = _panel_score_with_weights_by_avi(lp, patient_avi[p], w, top_k=k)
    sn = _panel_score_with_weights_by_avi(ln, patient_avi[p], w, top_k=k)
```

(Need to also track `patient_avi[p]` — the raw AVI norm vector — separately from `patient_weights[p]`.)

### Patch 3: Tests for new behavior

Add to `test/test_alphagenome_weights.py`:
- Test that indels are now kept (count of loaded mutations > old count)
- Test that indel keys have `source="proxy"` in weighted output
- Test that `_panel_score_with_weights_by_avi` returns deterministic output
- Test that the new method appears in the methods list

### Effort estimate

- Patch 1: 5 lines change + 1 docstring update + 1 test
- Patch 2: 15 lines change (helper + dispatch) + 1 test
- Patch 3: 4 new tests

Total: ~25 LOC + ~80 LOC tests. Should land in 1 subagent dispatch.

## Bottom line

3 honest issues, 3 simple fixes. No API key required (proxies handle indels; new method works on whatever weights are present). After fix:
- Indel count: 0 → 992 (~0.79%)
- Methods: 5 → 8 (uniform + AVI + 3 top-K-by-w*LLR + 3 top-K-by-AVI)
- All 3 fixes preserve the "honest proxy, clearly flagged" principle
