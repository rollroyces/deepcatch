# Synthetic-Bypass Bug — Foundation Pretrainer (real-data audit)

**Bug ID:** `pretrain.synthetic_bypass.v1`
**Discovered:** 2026-09-23 (PR-5 audit)
**Status:** FIXED in commit (uncommitted at time of audit)

## Summary

`src/foundation/pretrain.py` exposed a public `FoundationPretrainer`
class whose three phase methods (`pretrain_phase1_mmp`,
`pretrain_phase2_contrastive`, `pretrain_phase3_joint`) each called
`self.data_generator.generate_dataset(...)` unconditionally and
ignored any real modality dict the caller might have supplied.

The bug invalidated the claim made in `docs/PRETRAINING.md` that
the existing checkpoint
(`checkpoints/foundation_pretrained_finaledb.pt`) was "real-data-derived":
the encoder had **never** seen a real FinaleDB cfDNA sample — only the
synthetic generator's hash-deterministic output.

## Root cause

```python
# pre-fix src/foundation/pretrain.py (line ~170)
def pretrain_phase1_mmp(self, n_samples=5000, ...):
    ...
    modalities_np, _ = self.data_generator.generate_dataset(  # ← synthetic!
        n_samples=n_samples, prefix="phase1",
    )
    modalities = self._modalities_to_tensors(modalities_np)
    ...
```

The same pattern appeared in phase 2 (line ~272) and phase 3
(line ~363). The constructor instantiated `MultiModalDataGenerator`
unconditionally and stored it as `self.data_generator`; there was
no API path for callers to inject a real cohort.

The unit tests in `src/foundation/test_integration.py` (tests
20-24, 30) passed because they used the synthetic path, masking
the bug.

## Evidence at time of discovery

`scripts/pretrain_real_finaledb.py` assembled a real 6-modality
dict from the FinaleDB cohort and passed it to
`FoundationPretrainer(config=...)` — but `FoundationPretrainer`'s
constructor had no `modalities` parameter. The dict was
silently discarded. The phases then sampled from the synthetic
generator (hash-deterministic, label-free, not derived from any
real cfDNA fragments).

The end-to-end verification at the bottom of
`pretrain_real_finaledb.py` *did* produce finite outputs, but
those outputs reflected synthetic-data weights, not real-data
weights.

## Fix

### 1. Constructor accepts a real cohort

```python
FoundationPretrainer(
    config=PRODUCTION_CONFIG,
    device="cpu",
    modalities=real_modalities,         # NEW
    use_real_modalities=True,           # NEW (default True)
)
```

`use_real_modalities=True` makes phases draw mini-batches from
the real cohort instead of calling `self.data_generator`.
If no cohort is supplied the flag is silently degraded to `False`
(a warning is logged) so legacy callers don't crash.

### 2. Per-phase overrides

Each phase method gained two optional keyword arguments:

```python
pretrainer.pretrain_phase1_mmp(
    n_epochs=10,
    batch_size=32,
    modalities=cohort,                # per-phase override
    use_real_modalities=True,         # per-phase override (None = use default)
)
```

The per-phase arguments override the constructor-level defaults,
enabling mixed regimes (real-data phase 1, synthetic phase 2, …)
and forward-compatible call sites.

### 3. Internal helper

A new private method `_resolve_modalities(n_samples, prefix,
override_modalities, override_use_real)` returns the resolved
modality tensor dict and a `used_real: bool` flag, replacing the
three independent `self.data_generator.generate_dataset(...)` call
sites.

## Regression test

`test/test_pretrain_bug_fix.py` (11 tests, all passing) pins the
fix:

- `test_phase1_with_real_modalities_skips_synthetic_generator`
- `test_phase2_with_real_modalities_skips_synthetic_generator`
- `test_phase3_with_real_modalities_skips_synthetic_generator`
- `test_encoder_input_is_real_modalities_when_flag_is_true`
- `test_synthetic_path_still_calls_generator_when_flag_is_false`
- `test_synthetic_fallback_when_flag_true_but_no_modalities`
- `test_constructor_modalities_used_by_default_in_phase1`
- `test_per_phase_modalities_override_constructor`
- `test_per_phase_flag_false_overrides_constructor_true`
- `test_full_pretrain_pipeline_with_real_modalities`
- `test_encoder_forward_on_real_modalities_is_finite`

The three "skips_synthetic_generator" tests monkey-patch
`self.data_generator.generate_dataset` with a call counter and
assert it stays at zero across all three phases when
`use_real_modalities=True` is set. Pre-fix, the counter would
have shown 3 calls (one per phase).

The "encoder_input_is_real_modalities_when_flag_is_true" test
spies on `encoder.forward`, captures the first sample of each
modality the encoder receives, and asserts equality against the
real cohort. Pre-fix, the encoder would have seen synthetic
hash-derived values that share no overlap with the real cohort.

## Impact

### Before fix

- `scripts/pretrain_real_finaledb.py` saved checkpoints with
  synthetic-data weights labelled as "real-data-derived".
- `FoundationDownstream(pretrained=True, ...)` downstream of those
  checkpoints loaded real-data-shaped weights that had been trained
  on hash-deterministic synthetic cfDNA. Downstream cancer-detection
  results were uninformative at best and misleading at worst.

### After fix

- New checkpoints trained via the fixed pretrainer with the
  FinaleDB cohort (`checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt`)
  are real-data-trained end to end.
- Synthetic path is preserved (controlled by the flag) for
  CI / unit tests / scenarios without real data.
- Pretraining logs now distinguish "real" vs "synthetic" sources.

## Files touched

- `src/foundation/pretrain.py` — bug fix (the core change).
- `test/test_pretrain_bug_fix.py` — new regression test file.
- `docs/PRETRAIN_BUG.md` — this audit document.
- `docs/PRETRAINING.md` — CRITICAL UPDATE section + PRODUCTION_CONFIG
  pretraining documentation.
- `scripts/pretrain_real_finaledb.py` — not modified (it already
  built the cohort correctly; the bug was in the pretrainer it
  called). Re-running it via the new
  `scripts/pretrain_production_finaledb.py` driver produces the
  real-data-trained checkpoint.

## Lessons learned

1. **Don't let a single utility silently own a side effect.**
   `self.data_generator.generate_dataset()` was called with no
   guard, no override, no log line. The unit tests passed because
   they used the same synthetic path. The integration test was
   missing — there was no test that confirmed real data flowed
   through to the encoder.
2. **Test the data path, not just the loss shape.** A test that
   asserts `loss.item() > 0` after pretraining is necessary but
   not sufficient. The bug-fix test asserts what the encoder
   *saw*, not what it computed.
3. **Public API surface matters.** The constructor had no
   `modalities` parameter, so callers had no way to inject real
   data — even if they'd noticed the bug, they couldn't fix it
   from outside.
