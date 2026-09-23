#!/usr/bin/env python3
"""PRODUCTION_CONFIG pretraining on real FinaleDB cfDNA cohort.

Driver for the real-data-trained checkpoint produced after the
synthetic-bypass bug fix (see ``docs/PRETRAIN_BUG.md``).

Differences vs ``scripts/pretrain_real_finaledb.py``:

- Uses ``FoundationConfig(embed_dim=128, n_layers=4, n_heads=4,
  ff_dim=256, dropout=0.2)`` — the ``PRODUCTION_CONFIG`` preset
  from ``src.foundation.config`` (≈250k encoder parameters vs
  73k for PROTOTYPE_CONFIG).
- Caps the cohort at ``--max-samples`` (default 200) so the run
  finishes inside the wall-clock budget.
- Calls the FIXED ``FoundationPretrainer(modalities=..., use_real_modalities=True)``
  constructor — the encoder sees the real FinaleDB cfDNA samples,
  not synthetic hash-deterministic noise.
- Phase 1 (MMP) for ``--p1-epochs`` epochs (default 15) +
  Phase 2 (contrastive) for ``--p2-epochs`` epochs (default 5).
  Phase 3 (joint) is intentionally skipped to stay under the
  wall-clock budget — the MMP + contrastive objectives capture
  the bulk of the cross-modal signal.
- Saves the checkpoint to
  ``checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt``
  so it does not overwrite the 16-sample PROTOTYPE checkpoint.

Run::

    env -u PYTHONPATH ./.venv/bin/python scripts/pretrain_production_finaledb.py \\
        --max-samples 200 --p1-epochs 15 --p2-epochs 5 \\
        --batch-size 32 --device cpu --seed 42

Default wall-clock target: <30 minutes on M4 CPU for 200 samples.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Reuse the cohort builder from pretrain_real_finaledb.py.
from pretrain_real_finaledb import (  # noqa: E402
    PIPELINE_FEAT_DIR,
    load_labels,
    pick_cohort,
    build_modality_dict,
)

COHORT_OUT = REPO_ROOT / "data" / "finaledb_pretrain_cohort_PRODUCTION.npz"
CHECKPOINT_OUT = (
    REPO_ROOT / "checkpoints" / "foundation_pretrained_finaledb_PRODUCTION.pt"
)
LOG_OUT = REPO_ROOT / "results" / "pretrain_production_finaledb.json"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "PRODUCTION_CONFIG pretraining on the real FinaleDB cfDNA "
            "cohort (post synthetic-bypass fix)."
        )
    )
    ap.add_argument("--max-samples", type=int, default=200,
                    help="Hard cap on cohort size (default 200).")
    ap.add_argument("--n-healthy", type=int, default=0,
                    help="0 = auto (half of --max-samples).")
    ap.add_argument("--n-cancer", type=int, default=0,
                    help="0 = auto (half of --max-samples).")
    ap.add_argument("--studies", default="cristiano,jiang")
    ap.add_argument("--p1-epochs", type=int, default=15,
                    help="Phase 1 (MMP) epochs.")
    ap.add_argument("--p2-epochs", type=int, default=5,
                    help="Phase 2 (contrastive) epochs.")
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--device", default="cpu",
                    choices=["cpu", "mps", "cuda"])
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    # Resolve cohort size.
    if args.n_healthy == 0 and args.n_cancer == 0:
        n_healthy = args.max_samples // 2
        n_cancer = args.max_samples - n_healthy
    else:
        n_healthy = args.n_healthy or args.max_samples // 2
        n_cancer = args.n_cancer or (args.max_samples - n_healthy)
    n_healthy = min(n_healthy, args.max_samples // 2)
    n_cancer = min(n_cancer, args.max_samples - n_healthy)

    print("=" * 70)
    print("PRODUCTION_CONFIG Pretraining on Real FinaleDB cfDNA")
    print("  (post synthetic-bypass fix — see docs/PRETRAIN_BUG.md)")
    print("=" * 70)
    t0 = time.time()
    print(f"max_samples={args.max_samples}  n_healthy={n_healthy}  "
          f"n_cancer={n_cancer}  p1_epochs={args.p1_epochs}  "
          f"p2_epochs={args.p2_epochs}  batch_size={args.batch_size}  "
          f"device={args.device}")

    studies = tuple(s.strip() for s in args.studies.split(",") if s.strip())

    # ── Load + assemble cohort ─────────────────────────────────────
    labels = load_labels()
    healthy_ids, cancer_ids = pick_cohort(
        labels, n_healthy, n_cancer, studies, args.seed,
    )
    sample_ids = healthy_ids + cancer_ids
    y = np.array(
        [0] * len(healthy_ids) + [1] * len(cancer_ids),
        dtype=np.int64,
    )
    print(f"\nCohort: {len(sample_ids)} samples "
          f"({len(healthy_ids)} healthy + {len(cancer_ids)} cancer)  "
          f"studies={studies}")
    print(f"  healthy: {healthy_ids[:5]}{'…' if len(healthy_ids) > 5 else ''}")
    print(f"  cancer:  {cancer_ids[:5]}{'…' if len(cancer_ids) > 5 else ''}")

    modalities, X, ordered_ids = build_modality_dict(sample_ids)
    print(f"\nModality dict:")
    for k, v in modalities.items():
        print(f"  {k}: shape={v.shape}  mean={v.mean():.3f}  "
              f"std={v.std():.3f}")
    print(f"Flat X.shape = {X.shape}")

    # ── Save cohort .npz ───────────────────────────────────────────
    COHORT_OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        COHORT_OUT,
        X=X,
        y=y,
        sample_ids=np.array(ordered_ids),
        studies=np.array([labels[s]["study"] for s in ordered_ids]),
        feature_dim=X.shape[1],
        n_healthy=int((y == 0).sum()),
        n_cancer=int((y == 1).sum()),
    )
    print(f"\nWrote cohort → {COHORT_OUT}  "
          f"({COHORT_OUT.stat().st_size/1024:.1f} KB)")

    # ── Pretrain on REAL data ──────────────────────────────────────
    from src.foundation import (
        FoundationPretrainer,
        FoundationDownstream,
        PRODUCTION_CONFIG,
    )
    import torch

    config = PRODUCTION_CONFIG
    config.batch_size = args.batch_size
    config.device = args.device
    # Override dropout to match the brief (0.2).
    config.dropout = 0.2
    # Embed dim 128, n_layers 4, n_heads 4, ff_dim 256 are already set
    # in PRODUCTION_CONFIG.

    # CRITICAL: pass the real modalities to the pretrainer and let
    # the post-fix constructor use them. This is the whole point of
    # the bug-fix — without these two arguments the pretrainer would
    # silently fall back to synthetic data.
    pretrainer = FoundationPretrainer(
        config=config,
        device=args.device,
        verbose=True,
        modalities=modalities,
        use_real_modalities=True,
    )
    print(f"\nEncoder params: {pretrainer.encoder.num_params:,}")
    print(f"  modality_dims: {list(pretrainer.encoder.modality_dims.keys())}")

    train_time = 0.0
    losses: Dict[str, List[float]] = {"phase1": [], "phase2": [], "phase3": []}

    t_train = time.time()
    print("\n[Phase 1] MMP on real FinaleDB cohort")
    losses["phase1"] = pretrainer.pretrain_phase1_mmp(
        n_samples=len(sample_ids),
        n_epochs=args.p1_epochs,
        batch_size=args.batch_size,
    )
    print(f"  phase1 final loss: {losses['phase1'][-1]:.6f}")

    print("\n[Phase 2] Contrastive on real FinaleDB cohort")
    losses["phase2"] = pretrainer.pretrain_phase2_contrastive(
        n_samples=len(sample_ids),
        n_epochs=args.p2_epochs,
        batch_size=args.batch_size,
    )
    print(f"  phase2 final loss: {losses['phase2'][-1]:.6f}")

    train_time = time.time() - t_train
    print(f"\nPretrain wall-clock: {train_time:.1f}s")

    # ── Save checkpoint ────────────────────────────────────────────
    CHECKPOINT_OUT.parent.mkdir(parents=True, exist_ok=True)
    pretrainer.save_checkpoint(str(CHECKPOINT_OUT))
    ckpt_size_mb = CHECKPOINT_OUT.stat().st_size / (1024 * 1024)
    print(f"\nCheckpoint → {CHECKPOINT_OUT}  ({ckpt_size_mb:.2f} MB)")

    # ── Verify FoundationDownstream loads it ───────────────────────
    print(f"\n[Verify] FoundationDownstream(pretrained=True, "
          f"checkpoint_path='{CHECKPOINT_OUT}')")
    fd = FoundationDownstream(
        config=config,
        pretrained=True,
        checkpoint_path=str(CHECKPOINT_OUT),
        device=args.device,
    )
    fd._build_classifier(n_classes=2)
    fd.encoder.eval()
    with torch.no_grad():
        mod_t = {
            k: torch.from_numpy(v).to(args.device)
            for k, v in modalities.items()
        }
        joint = fd.encoder(mod_t)
    out_finite = bool(torch.isfinite(joint).all().item())
    print(f"  forward output: shape={tuple(joint.shape)}  finite={out_finite}")
    if not out_finite:
        raise RuntimeError(
            "Forward pass produced NaN/Inf — checkpoint is broken"
        )

    # ── Log ────────────────────────────────────────────────────────
    LOG_OUT.parent.mkdir(parents=True, exist_ok=True)
    log = {
        "checkpoint_path": str(CHECKPOINT_OUT),
        "checkpoint_size_mb": round(ckpt_size_mb, 3),
        "cohort_npz": str(COHORT_OUT),
        "n_samples": int(X.shape[0]),
        "n_healthy": int((y == 0).sum()),
        "n_cancer": int((y == 1).sum()),
        "studies": list(studies),
        "feature_dim": int(X.shape[1]),
        "config": config.to_dict(),
        "encoder_num_params": int(pretrainer.encoder.num_params),
        "p1_epochs": args.p1_epochs,
        "p2_epochs": args.p2_epochs,
        "batch_size": args.batch_size,
        "device": args.device,
        "seed": args.seed,
        "pretrain_losses": {
            k: [round(x, 6) for x in v] for k, v in losses.items()
        },
        "train_time_s": round(train_time, 1),
        "wall_clock_total_s": round(time.time() - t0, 1),
        "forward_pass_finite": out_finite,
        "data_source": "real FinaleDB cfDNA cohort (post-fix)",
        "bug_fix": "see docs/PRETRAIN_BUG.md",
        "sample_ids": sample_ids,
    }
    with open(LOG_OUT, "w") as fh:
        json.dump(log, fh, indent=2)
    print(f"\nLog → {LOG_OUT}")
    print(f"Total wall-clock: {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
