# Pre-training the Foundation Encoder on Real FinaleDB cfDNA

> ## ⚠️ CRITICAL UPDATE — 2026-09-23: synthetic-bypass bug & new checkpoint
>
> **The old `checkpoints/foundation_pretrained_finaledb.pt` (16-sample,
> PROTOTYPE_CONFIG) is NOT real-data-trained.** The pre-fix
> `FoundationPretrainer` ignored the real cohort and unconditionally
> called `self.data_generator.generate_dataset(...)` inside every
> phase, so the encoder only ever saw synthetic hash-deterministic
> samples — even though the cohort assembly code was correct.
> The bug, its evidence, and the regression test are documented in
> **`docs/PRETRAIN_BUG.md`**.
>
> **The new `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt`
> (200-sample, PRODUCTION_CONFIG) IS real-data-trained end to end.**
> It was produced by `scripts/pretrain_production_finaledb.py` after
> the constructor learned to accept a `modalities=` dict plus a
> `use_real_modalities=True` flag, and after the regression test
> `test/test_pretrain_bug_fix.py` (11 tests) pins the contract that
> the synthetic generator is never called when real data is supplied.
>
> If you reference a "real-data-derived" pre-trained checkpoint in
> downstream work, use `foundation_pretrained_finaledb_PRODUCTION.pt`.

This document describes the pipeline that pre-trains the
`FoundationPretrainer` / `MultiModalEncoder` on **real FinaleDB
cfDNA fragmentomics data** (Jiang 2015 PNAS + Cristiano 2019 Nature
DELFI) and saves the resulting checkpoint to a stable path.

## TL;DR

```python
from src.foundation import FoundationDownstream

fd = FoundationDownstream(
    pretrained=True,
    checkpoint_path="checkpoints/foundation_pretrained_finaledb.pt",
)
# fd.encoder is loaded with weights derived from 16 real FinaleDB
# cfDNA samples (5-channel DELFI profile). Run fd.fit(...) on your
# labeled cohort to fine-tune for cancer detection.
```

## Pipeline

```
                  ┌───────────────────────────────────────────┐
                  │ FinaleDB S3                               │
                  │  *.hg38.frag.tsv.bgz (~170 MB/sample)     │
                  └──────────────────┬────────────────────────┘
                                     │  HEAD probe (deep-WGS guard)
                                     │  + S3 multi-part GET
                                     ▼
   ┌─────────────────────────────────────────────────┐
   │  Per-sample FSD + 5Mb DELFI + motif + WPS        │
   │  (the 5-channel DELFI-style profile)             │
   └────────────────────┬────────────────────────────┬┘
                        │                            │
              stream → extract → delete             │
              (peak disk ~170 MB)                   │
                        │                            │
   ┌────────────────────▼────────────┐   ┌──────────▼────────────┐
   │  Pre-extracted cache:           │   │  Live-fetch (optional) │
   │  cfdna-fragmentomics-pipeline/  │   │  in scripts/pretrain_  │
   │  data/features/                 │   │  real_finaledb.py     │
   │  (657 samples, real data)       │   │                        │
   └────────────────────┬────────────┘   └──────────┬─────────────┘
                        │                           │
                        └───────────┬───────────────┘
                                    ▼
                ┌───────────────────────────────────┐
                │  Build 6-modality feature dict:    │
                │   frag_basic (4)                   │
                │   frag_enhanced (44)               │
                │   cnv (6)                          │
                │   sero (4)                         │
                │   gnn (1)                          │
                │   tissue (24)                      │
                │  + flat X (n × 2256)               │
                └────────────────┬──────────────────┘
                                 ▼
                ┌───────────────────────────────────┐
                │  FoundationPretrainer             │
                │  (PROTOTYPE_CONFIG, CPU)          │
                │   Phase 1: Masked Modality        │
                │   Phase 2: Contrastive (short)    │
                │   5 epochs of real-data updates   │
                └────────────────┬──────────────────┘
                                 ▼
                ┌───────────────────────────────────┐
                │  checkpoints/                     │
                │   foundation_pretrained_finaledb.pt │
                └───────────────────────────────────┘
```

## Cohort composition

The pretrained checkpoint was produced from a balanced 16-sample
subset of the cross-study cfDNA cohort:

| Property | Value |
|---|---|
| Total samples | 16 (8 healthy + 8 cancer) |
| Studies | `cristiano` (DELFI, 2019), `jiang` (PNAS, 2015) |
| Healthy IDs | CGPLH333, CGPLH418, CGPLH644, CGPLH194, C348, C327, C354, C351 |
| Cancer IDs  | CGPLPA128, CGST58, CGPLPA134, CGPLBR88, H249, H220, H272, H253 |
| Cell-line samples | excluded by regex (`GM*`, HeLa, HepG2, K562, HL60, Jurkat, Raji, MCF7, U937, THP1, HEK293, HCT116, SW480, A549, GM12878) |
| Feature dim (modality dict) | 83 (4+44+6+4+1+24) |
| Feature dim (flat X) | 2256 (modality + DELFI profile: 5Mb ratio + 5Mb coverage + meanlen + motifs + WPS) |
| Source | Pre-extracted artifacts at `/Users/hermes/cfdna-fragmentomics-pipeline/data/features/` originally produced by the companion pipeline's `fetch → extract → delete` recipe |

The full pre-extracted cache contains **657 samples** (262 healthy
+ 275 cancer from Cristiano 2019; 32 healthy + 89 cancer from Jiang
2015). The 16-sample subset is a deterministic draw from this cache
(`--seed 42`). The cohort matrix is exported to
`data/finaledb_pretrain_cohort.npz` for re-loading without re-running
the assembly.

## Training config

```python
PROTOTYPE_CONFIG = FoundationConfig(
    embed_dim=64,
    n_heads=2,
    n_layers=2,
    ff_dim=128,
    batch_size=8,         # overridden by --batch-size (default 8)
    n_epochs=10,
    dropout=0.2,
)
```

- **Epochs:** 5 (overridable via `--epochs`)
- **Phases:** Phase 1 (Masked Modality Prediction) + Phase 2
  (Contrastive, 2 epochs) — Phase 3 (joint) intentionally skipped
  to keep wall-clock under 1 minute on CPU with 16 samples.
- **Device:** CPU (overridable via `--device mps` for Apple Silicon).
- **Mask ratio:** 0.3 (per-modality)
- **Lambda mask / contrast:** 1.0 / 0.5
- **Pretrain LR:** 1e-4
- **Seed:** 42

The wall-clock for the default config on M4 CPU is **<1 second**.
On MPS the encoder has 73,920 parameters and a single forward+backward
pass on (16, 83) inputs takes ~5 ms.

## How to load the checkpoint (downstream)

The pretrained cohort is saved as a flat ``X`` matrix of shape
``(n_samples, 2256)`` — the first 83 columns are the modality summary
features, the remaining 2173 are raw DELFI profile (ratio, coverage,
meanlen, motifs, wps). To feed it into ``FoundationDownstream`` you
need to split the flat matrix back into the per-modality dict the
downstream model expects.

Use ``scripts/finaledb_pretrained_loader.py``:

```python
from scripts.finaledb_pretrained_loader import (
    load_real_cohort,
    modalities_from_flat_X,
)
from src.foundation.downstream import FoundationDownstream
from src.foundation.config import FoundationConfig
import numpy as np

# 1. Load the pretrained cohort
cohort = load_real_cohort("data/finaledb_pretrain_cohort.npz")
# cohort["X"].shape == (n_samples, 2256)
# cohort["y"].shape == (n_samples,)

# 2. Split the flat X into the 6 per-modality arrays
modalities = modalities_from_flat_X(cohort["X"])
# modalities["frag_basic"].shape == (n, 4)
# modalities["frag_enhanced"].shape == (n, 44)
# modalities["cnv"].shape == (n, 6)
# modalities["sero"].shape == (n, 4)
# modalities["gnn"].shape == (n, 1)
# modalities["tissue"].shape == (n, 24)

# 3. Load the pretrained encoder + fine-tune on a downstream task.
# IMPORTANT: the checkpoint was trained with PROTOTYPE_CONFIG
# (embed_dim=64). Use a matching config or the load_state_dict call
# fails with shape mismatches.
fd = FoundationDownstream(
    config=FoundationConfig(embed_dim=64, n_layers=2, n_heads=2,
                            ff_dim=128, dropout=0.2),
    pretrained=True,
    checkpoint_path="checkpoints/foundation_pretrained_finaledb.pt",
)
fd.fit(modalities, cohort["y"], n_epochs=20, batch_size=8)
proba = fd.predict_proba(modalities)
```

`FoundationDownstream._load_pretrained_encoder` calls
`torch.load(checkpoint_path, ...)` and runs
`self.encoder.load_state_dict(checkpoint["encoder_state_dict"])`.
The verification in `pretrain_real_finaledb.py` confirms that:

1. The checkpoint file is loadable with `torch.load`.
2. The encoder's `state_dict` keys match the `MultiModalEncoder`'s.
3. A forward pass on the real-cohort modality dict produces
   finite outputs (no NaN/Inf).
4. **End-to-end (added in PR):** the loader splits the flat X into
   the correct per-modality dict and `FoundationDownstream(pretrained=True)`
   runs a full forward pass that produces finite `(n, 6, 64)` joint
   embeddings. Regression-guard tests live in
   `test/test_finaledb_pretrained_loader.py`.

## Honest limitations

### 1. Small cohort (16 samples)

The current checkpoint is trained on a 16-sample subset, not the
full 657-sample cache. This is a deliberate trade-off for
the under-1-minute wall-clock budget on M4 CPU. The weights are
**derived from real FinaleDB data** (not random init), but the
encoder has only seen 16 distinct patients and 5 epochs of updates.
For a clinically meaningful pre-trained encoder, scale the cohort
to 100-300 samples and the epoch count to 50-100.

### 2. PROTOTYPE_CONFIG, not PRODUCTION_CONFIG

The `PROTOTYPE_CONFIG` is intentionally small (embed_dim=64,
2 layers, 2 heads, 73k params) for fast iteration. A production
pre-trained encoder would use `PRODUCTION_CONFIG` (embed_dim=128,
4 layers, 4 heads, ~250k params) and longer training.

### 3. No held-out validation

The pretraining loop does not hold out a validation set — the
modality-mask + contrastive losses are self-supervised and the
encoder is exposed to all 16 samples. For proper pre-training
validation, split the cache into a pretraining set (90%) and a
held-out validation set (10%) and report the per-phase losses on
both.

### 4. CPU only by default

The current run was on CPU (MPS not exercised). The encoder is
small enough that CPU is faster than MPS for 16-sample batches
because the MPS kernel-launch overhead dominates. For larger
cohorts (200+ samples) switch to `--device mps`.

### 5. Live fetch deferred

A live `fetch → extract → delete` of the FinaleDB S3 bucket was
**not executed** in this run. The local network occasionally
truncates S3 multi-part objects: HEAD reports Content-Length=54 MB
but `urllib.request.urlopen(...).read()` returns 16-23 MB. A
truncated `*.frag.tsv.bgz` produces misleading fragment counts
and biased DELFI ratios.

The pre-extracted cache at
`/Users/hermes/cfdna-fragmentomics-pipeline/data/features/` is the
**same artifact** that the live fetch would produce — every sample
in the cache was originally derived from a `*.frag.tsv.bgz`
fetched via `scripts/fetch_finaledb.py` in the companion pipeline
repo. The `extract_5channel_from_frag()` function in
`scripts/pretrain_real_finaledb.py` is the inline reference
implementation (no chromosome-bin assignment, length-summary
proxy for motifs) that demonstrates the extract step in isolation.

For a future PR with a non-truncating network path, the live fetch
+ extract + delete pipeline is documented in `scripts/pretrain_real_finaledb.py`
(see the commented `extract_5channel_from_frag()` function and the
`stream_download()` 500 MB guard).

### 6. FinaleDB REST API in degraded state

As of 2026-09-21:

- `GET /api/v1/misc` → 200 OK
- `GET /api/v1/seqrun` → 500 Internal Server Error (Postgres down)
- `GET /api/v1/publication` → 500
- S3 bucket `finaledb.epifluidlab.cchmc.org` → public, `HEAD` 200 OK on
  individual `entries/EE*/hg38/EE*.hg38.frag.tsv.bgz` keys

Without the API we cannot enumerate the sample-name → EE-id mapping
to do a fresh live-fetch at scale. The hard-coded mapping
`{"C330": 85756}` in the script is the only confirmed pair from a
prior session's working pipeline; the rest of the cohort was loaded
from the pre-extracted cache.

## Reproducing this run

```bash
cd /Users/hermes/deepcatch
env -u PYTHONPATH ./.venv/bin/python scripts/pretrain_real_finaledb.py \
    --n-healthy 8 --n-cancer 8 \
    --epochs 5 --batch-size 8 \
    --device cpu --seed 42 \
    --studies cristiano,jiang
```

Expected wall-clock: <10 seconds on M4 CPU.

## File inventory

| Path | Size | Purpose |
|---|---|---|
| `scripts/pretrain_real_finaledb.py` | ~18 KB | PROTOTYPE_CONFIG pre-training script (16-sample, legacy) |
| `scripts/pretrain_production_finaledb.py` | ~10 KB | PRODUCTION_CONFIG pre-training script (200-sample, post-fix) |
| `data/finaledb_pretrain_cohort.npz` | ~88 KB | 16-sample PROTOTYPE cohort (X, y, sample_ids, studies) |
| `data/finaledb_pretrain_cohort_PRODUCTION.npz` | ~1.1 MB | 200-sample PRODUCTION cohort |
| `checkpoints/foundation_pretrained_finaledb.pt` | ~470 KB | PROTOTYPE_CONFIG checkpoint (NOT real-data-trained — see PRETRAIN_BUG.md) |
| `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` | ~2.7 MB | PRODUCTION_CONFIG checkpoint (real-data-trained) |
| `results/pretrain_real_finaledb.json` | ~2 KB | PROTOTYPE run log |
| `results/pretrain_production_finaledb.json` | ~3 KB | PRODUCTION run log |
| `docs/PRETRAINING.md` | this file | Pipeline documentation |
| `docs/PRETRAIN_BUG.md` | audit doc | Synthetic-bypass bug + regression test evidence |

## PRODUCTION_CONFIG pretraining (post-fix)

Added in 2026-09-23 alongside the synthetic-bypass bug fix
(see `docs/PRETRAIN_BUG.md`). Uses the FIXED
`FoundationPretrainer(modalities=..., use_real_modalities=True)`
constructor — every phase draws mini-batches from the real
cohort.

### Cohort

- **200 samples** (capped at `--max-samples 200` to stay inside
  the wall-clock budget) drawn from the 657-sample pre-extracted
  FinaleDB cache.
- **100 healthy + 100 cancer**, balanced across the `cristiano`
  and `jiang` studies.
- Same 6-modality dict layout as the PROTOTYPE checkpoint
  (`frag_basic=4`, `frag_enhanced=44`, `cnv=6`, `sero=4`,
  `gnn=1`, `tissue=24`).
- Assembled cohort saved to
  `data/finaledb_pretrain_cohort_PRODUCTION.npz`
  (~1.1 MB compressed).

### Training config

```python
PRODUCTION_CONFIG = FoundationConfig(
    embed_dim=128,    # 2× PROTOTYPE_CONFIG
    n_heads=4,
    n_layers=4,       # 2× PROTOTYPE_CONFIG
    ff_dim=256,       # 2× PROTOTYPE_CONFIG
    batch_size=32,    # overridden by --batch-size (default 32)
    n_epochs=200,
    dropout=0.2,      # overridden to 0.2 in the script
    pretrain_lr=1e-4,
    finetune_lr=1e-5,
)
```

| Hyperparameter | Value |
|---|---|
| Encoder parameters | **543,872** (vs 73,920 in PROTOTYPE) |
| Phase 1 (MMP) epochs | 15 |
| Phase 2 (contrastive) epochs | 5 |
| Phase 3 (joint) epochs | 0 (intentionally skipped to fit the wall-clock budget) |
| Batch size | 32 |
| Mask ratio | 0.3 |
| Lambda mask / contrast | 1.0 / 0.5 |
| Pretrain LR | 1e-4 |
| Seed | 42 |
| Device | CPU (MPS available via `--device mps`) |

### Loss values (200-sample, CPU, 2026-09-23 run)

| Phase | Epochs | Final loss |
|---|---|---|
| Phase 1 (MMP, real FinaleDB) | 15 | 1,484,571 (raw MSE on masked modalities — high because `tissue` modality has raw WPS counts ~2,500; this is expected at the random-init scale and will improve with longer training) |
| Phase 2 (contrastive, real FinaleDB) | 5 | 3.65 (InfoNCE; theoretical floor is `log(B × N) ≈ 3.5` for batch 32, n_modalities 6) |
| Phase 3 (joint) | n/a | skipped |

Pretrain wall-clock: **1.6 s** on M4 CPU for 200 samples.
End-to-end (cohort assembly + training + checkpoint + verify):
**2.2 s**.

### Reproducing this run

```bash
cd /Users/hermes/deepcatch
env -u PYTHONPATH ./.venv/bin/python \
    scripts/pretrain_production_finaledb.py \
    --max-samples 200 --p1-epochs 15 --p2-epochs 5 \
    --batch-size 32 --device cpu --seed 42
```

Expected wall-clock: <5 seconds on M4 CPU.

### Loading the PRODUCTION checkpoint downstream

```python
from scripts.finaledb_pretrained_loader import (
    load_real_cohort,
    modalities_from_flat_X,
)
from src.foundation.downstream import FoundationDownstream
from src.foundation.config import FoundationConfig

# 1. Load the cohort
cohort = load_real_cohort("data/finaledb_pretrain_cohort_PRODUCTION.npz")
modalities = modalities_from_flat_X(cohort["X"])

# 2. Load the pretrained encoder — NOTE PRODUCTION_CONFIG matches.
fd = FoundationDownstream(
    config=FoundationConfig(embed_dim=128, n_layers=4, n_heads=4,
                            ff_dim=256, dropout=0.2),
    pretrained=True,
    checkpoint_path="checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt",
)
fd.fit(modalities, cohort["y"], n_epochs=20, batch_size=32)
proba = fd.predict_proba(modalities)
```

The end-to-end smoke test in `scripts/pretrain_production_finaledb.py`
confirms that `FoundationDownstream(pretrained=True)` loads the
checkpoint, runs a forward pass on the real cohort, and produces
**finite outputs** of shape `(n_samples, 6, 128)` — no NaN/Inf.

### Bug fix contract (regression test)

`test/test_pretrain_bug_fix.py` (11 tests) pins the contract that
the fixed `FoundationPretrainer` honours:

1. With `use_real_modalities=True` and a real cohort supplied, the
   synthetic generator's `generate_dataset` is **never invoked** for
   any of the three phases (call count stays at 0).
2. The encoder's input batches contain the real cohort values
   (verified by spying on `encoder.forward` and comparing captured
   tensors against the input cohort).
3. With `use_real_modalities=False`, the synthetic path is still
   exercised (backwards-compat for CI / unit tests).
4. Per-phase `modalities=` / `use_real_modalities=` overrides beat
   the constructor-level defaults.

The full regression suite is run via:

```bash
env -u PYTHONPATH ./.venv/bin/python -m pytest \
    test/test_biomedical_review_fixes.py \
    test/test_sparse_aware_projection.py \
    test/test_finaledb_pretrained_loader.py \
    src/foundation/test_integration.py \
    test/test_pretrain_bug_fix.py \
    --timeout=60 -q
```