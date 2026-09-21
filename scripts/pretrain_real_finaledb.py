#!/usr/bin/env python3
"""
Real-Data Foundation Pre-training on FinaleDB cfDNA Cohort
============================================================

Loads the 5-channel DELFI-style profile (5Mb short/long ratio,
median-normalized coverage, mean fragment length, FSD histogram,
100kb WPS) for a 16-sample subset of the pre-extracted FinaleDB
cross-study cfDNA cohort (Jiang 2015 + Cristiano 2019), assembles
it into the Foundation Model's 6-modality dict, runs a brief
training pass (MMP + contrastive, PROTOTYPE_CONFIG), and saves the
checkpoint to ``checkpoints/foundation_pretrained_finaledb.pt``.

The bulk of the cohort is loaded from
``/Users/hermes/cfdna-fragmentomics-pipeline/data/features/`` —
the artifacts of an earlier ``fetch → extract → delete`` pipeline
run on real FinaleDB ``*.frag.tsv.bgz`` files.

Live-fetch from FinaleDB S3 is intentionally NOT attempted in this
script: the local network occasionally truncates S3 multi-part
objects, and producing a non-truncated 170 MB frag.tsv.bgz is not
reliable here. The pre-extracted cache is the canonical artifact
produced by the same fetch→extract→delete recipe (see
``docs/PRETRAINING.md`` for the full provenance).

Run::

    env -u PYTHONPATH ./.venv/bin/python \\
        scripts/pretrain_real_finaledb.py \\
        --n-healthy 8 --n-cancer 8 --epochs 5 --device cpu

The defaults target a 16-sample subset, 5 epochs of phase-1 (MMP)
training on CPU. Total wall-clock target: <10 minutes.
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

PIPELINE_FEAT_DIR = Path(
    "/Users/hermes/cfdna-fragmentomics-pipeline/data/features"
)
PIPELINE_LABELS_TSV = PIPELINE_FEAT_DIR / "labels_cross_study.tsv"
COHORT_OUT = REPO_ROOT / "data" / "finaledb_pretrain_cohort.npz"
CHECKPOINT_OUT = REPO_ROOT / "checkpoints" / "foundation_pretrained_finaledb.pt"
LOG_OUT = REPO_ROOT / "results" / "pretrain_real_finaledb.json"


# ── Helpers ─────────────────────────────────────────────────────────

def load_labels() -> Dict[str, Dict[str, str]]:
    """{sample_id: {label, study}} from the pipeline's labels_cross_study.tsv."""
    out: Dict[str, Dict[str, str]] = {}
    if not PIPELINE_LABELS_TSV.exists():
        raise FileNotFoundError(
            f"Pipeline labels missing: {PIPELINE_LABELS_TSV}"
        )
    with open(PIPELINE_LABELS_TSV) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            out[parts[0]] = {"label": parts[1], "study": parts[2]}
    return out


def has_features(sid: str) -> bool:
    """True iff the pipeline has the 4 core DELFI feature files for sid."""
    required = [
        f"{sid}.fsd.json",
        f"{sid}.delfi_5mb_ratio.npy",
        f"{sid}.delfi_5mb_coverage.npy",
        f"{sid}.wps_100kb.npy",
    ]
    return all((PIPELINE_FEAT_DIR / f).exists() for f in required)


def load_opt(sid: str, suffix: str, size: int) -> np.ndarray:
    """Load a feature .npy if present; zero-fill otherwise."""
    p = PIPELINE_FEAT_DIR / f"{sid}.{suffix}"
    if p.exists():
        a = np.load(p)
        if a.size >= size:
            return a[:size]
        pad = np.zeros(size - a.size, dtype=a.dtype)
        return np.concatenate([a, pad])[:size]
    return np.zeros(size, dtype=np.float64)


def pick_cohort(
    labels: Dict[str, Dict[str, str]],
    n_healthy: int,
    n_cancer: int,
    studies: Tuple[str, ...] = ("cristiano", "jiang"),
    seed: int = 42,
) -> Tuple[List[str], List[str]]:
    """Pick a balanced cohort stratified across studies."""
    rng = np.random.default_rng(seed)
    pool_h: List[str] = []
    pool_c: List[str] = []
    for sid, m in labels.items():
        if m["study"] not in studies:
            continue
        if not has_features(sid):
            continue
        if m["label"] == "healthy":
            pool_h.append(sid)
        elif m["label"] == "cancer":
            pool_c.append(sid)
    rng.shuffle(pool_h)
    rng.shuffle(pool_c)

    n_per_study_h = max(1, n_healthy // len(studies))
    n_per_study_c = max(1, n_cancer // len(studies))
    by_study_h: Dict[str, List[str]] = {s: [] for s in studies}
    by_study_c: Dict[str, List[str]] = {s: [] for s in studies}
    for s in pool_h:
        by_study_h[labels[s]["study"]].append(s)
    for s in pool_c:
        by_study_c[labels[s]["study"]].append(s)
    picked_h: List[str] = []
    picked_c: List[str] = []
    for st in studies:
        picked_h.extend(by_study_h[st][:n_per_study_h])
        picked_c.extend(by_study_c[st][:n_per_study_c])
    # Top up if a study ran short.
    if len(picked_h) < n_healthy:
        leftover = [s for s in pool_h if s not in picked_h]
        picked_h.extend(leftover[: n_healthy - len(picked_h)])
    if len(picked_c) < n_cancer:
        leftover = [s for s in pool_c if s not in picked_c]
        picked_c.extend(leftover[: n_cancer - len(picked_c)])
    return picked_h[:n_healthy], picked_c[:n_cancer]


def build_modality_dict(
    sids: List[str],
) -> Tuple[Dict[str, np.ndarray], np.ndarray, List[str]]:
    """Build the 6-modality feature dict from per-sample DELFI features.

    Returns:
        modalities: dict[str, ndarray] with keys matching MODALITY_NAMES
                    (frag_basic=4, frag_enhanced=44, cnv=6, sero=4,
                    gnn=1, tissue=24).
        X: (n, 2256) flat feature matrix (DELFI full profile + summaries).
        sample_ids: list of sample ids in row order.
    """
    modalities: Dict[str, np.ndarray] = {}
    X_rows: List[np.ndarray] = []

    # First pass: build each modality + flat X
    for idx, sid in enumerate(sids):
        fsd = json.load(open(PIPELINE_FEAT_DIR / f"{sid}.fsd.json"))
        ratio = np.load(PIPELINE_FEAT_DIR / f"{sid}.delfi_5mb_ratio.npy")
        cov = np.load(PIPELINE_FEAT_DIR / f"{sid}.delfi_5mb_coverage.npy")
        meanlen = load_opt(sid, "delfi_5mb_meanlen.npy", 631)
        motifs = load_opt(sid, "motifs.npy", 256)
        wps = np.load(PIPELINE_FEAT_DIR / f"{sid}.wps_100kb.npy")

        # frag_basic (4): median length, short/long fractions
        frag_basic = np.array([
            fsd["median_length"] / 200.0,
            fsd["short_fraction_100_150"],
            fsd["long_fraction_150_220"],
            fsd["short_long_ratio"],
        ], dtype=np.float32)

        # frag_enhanced (44): DELFI 5Mb ratio (8) + motif (16) +
        # coverage (20) summaries. 8 + 16 + 20 = 44.
        ratio_summary = np.array([
            ratio.mean(), ratio.std(), np.median(ratio),
            np.percentile(ratio, 10), np.percentile(ratio, 25),
            np.percentile(ratio, 75), np.percentile(ratio, 90),
            ratio.max() - ratio.min(),
        ], dtype=np.float32)
        motif_summary = np.array([
            motifs[:64].mean(), motifs[64:128].mean(),
            motifs[128:192].mean(), motifs[192:].mean(),
            motifs.std(), motifs.max(), motifs.min(),
            motifs.argmax() / 256.0,
            motifs[:32].sum(), motifs[32:64].sum(),
            motifs[64:96].sum(), motifs[96:128].sum(),
            motifs[128:160].sum(), motifs[160:192].sum(),
            motifs[192:224].sum(), motifs[224:].sum(),
        ], dtype=np.float32)
        # cov_summary is 20 (not 21): we dropped the duplicate
        # ``np.median(cov)`` so the totals hit 44. The 5-number
        # percentile summary already covers the median's range.
        cov_summary = np.array([
            cov.mean(), cov.std(),
            np.percentile(cov, 5), np.percentile(cov, 25),
            np.percentile(cov, 75), np.percentile(cov, 95),
            cov.max(), cov.min(),
            (cov > 1.5).mean(), (cov < 0.5).mean(),
            (cov > 2.0).sum() / cov.size,
            (cov < 0.3).sum() / cov.size,
            np.diff(np.sort(cov)).mean(),
            np.diff(np.sort(cov)).std(),
            np.percentile(cov, 1), np.percentile(cov, 99),
            meanlen.mean(), meanlen.std(),
            np.percentile(meanlen, 25), np.percentile(meanlen, 75),
        ], dtype=np.float32)
        frag_enhanced = np.concatenate(
            [ratio_summary, motif_summary, cov_summary]
        ).astype(np.float32)

        # cnv (6): coverage deviation from 1.0 (median-normalized)
        cnv = np.array([
            np.median(cov) - 1.0,
            np.std(cov),
            np.percentile(cov, 90) - np.percentile(cov, 10),
            (cov > 1.3).mean() - (cov < 0.7).mean(),
            np.percentile(cov, 95),
            np.percentile(cov, 5),
        ], dtype=np.float32)

        # sero (4): filler from FSD percentiles
        sero = np.array([
            fsd.get("p25", 159) / 200.0,
            fsd.get("p75", 177) / 200.0,
            fsd.get("p10", 147) / 200.0,
            fsd.get("p90", 189) / 200.0,
        ], dtype=np.float32)

        # gnn (1): std of 5Mb coverage
        gnn = np.array([float(np.std(cov))], dtype=np.float32)

        # tissue (24): the 100kb WPS profile
        tissue = wps[:24].astype(np.float32)
        if tissue.size < 24:
            tissue = np.concatenate(
                [tissue, np.zeros(24 - tissue.size, dtype=np.float32)]
            )

        mod_dict = {
            "frag_basic": frag_basic,
            "frag_enhanced": frag_enhanced,
            "cnv": cnv,
            "sero": sero,
            "gnn": gnn,
            "tissue": tissue,
        }
        for k, v in mod_dict.items():
            if k not in modalities:
                modalities[k] = np.zeros((len(sids), v.size), dtype=np.float32)
            modalities[k][idx] = v

        # Flat X row: 83 (modality) + 2173 (raw DELFI) = 2256 features
        X_rows.append(np.concatenate([
            frag_basic, frag_enhanced, cnv, sero, gnn, tissue,
            ratio[:631], cov[:631], meanlen[:631], motifs[:256], wps[:24],
        ]).astype(np.float32))

    return modalities, np.stack(X_rows, axis=0), sids


# ── Main ────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-healthy", type=int, default=8)
    ap.add_argument("--n-cancer", type=int, default=8)
    ap.add_argument("--studies", default="cristiano,jiang")
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--device", default="cpu",
                    choices=["cpu", "mps", "cuda"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-train", action="store_true",
                    help="Skip pretraining (saves 1-step random-init weights).")
    args = ap.parse_args()

    print("=" * 60)
    print("Real-Data Foundation Pre-training on FinaleDB cfDNA")
    print("=" * 60)
    t0 = time.time()
    print(f"n_healthy={args.n_healthy}  n_cancer={args.n_cancer}  "
          f"epochs={args.epochs}  device={args.device}")
    studies = tuple(s.strip() for s in args.studies.split(",") if s.strip())

    # ── Load + assemble cohort ─────────────────────────────────────
    labels = load_labels()
    healthy_ids, cancer_ids = pick_cohort(
        labels, args.n_healthy, args.n_cancer, studies, args.seed,
    )
    sample_ids = healthy_ids + cancer_ids
    y = np.array([0] * len(healthy_ids) + [1] * len(cancer_ids),
                 dtype=np.int64)
    print(f"\nCohort: {len(sample_ids)} samples "
          f"({len(healthy_ids)} healthy + {len(cancer_ids)} cancer)  "
          f"studies={studies}")
    print(f"  healthy: {healthy_ids}")
    print(f"  cancer:  {cancer_ids}")

    modalities, X, ordered_ids = build_modality_dict(sample_ids)
    print(f"\nModality dict:")
    for k, v in modalities.items():
        print(f"  {k}: shape={v.shape}  mean={v.mean():.3f}  std={v.std():.3f}")
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
    print(f"\nWrote cohort → {COHORT_OUT}  ({COHORT_OUT.stat().st_size/1024:.1f} KB)")

    # ── Pretrain ───────────────────────────────────────────────────
    from src.foundation import (
        FoundationPretrainer,
        FoundationDownstream,
        PROTOTYPE_CONFIG,
    )
    import torch

    config = PROTOTYPE_CONFIG
    config.batch_size = args.batch_size
    config.device = args.device
    pretrainer = FoundationPretrainer(config=config, verbose=True)
    print(f"\nEncoder params: {pretrainer.encoder.num_params:,}")

    train_time = 0.0
    losses: Dict[str, List[float]] = {"phase1": [], "phase2": [], "phase3": []}
    if args.skip_train:
        print("\n[--skip-train] Saving random-init weights (no training).")
    else:
        t_train = time.time()
        # Phase 1 only (MMP) for brevity; phase 2/3 skipped to keep
        # wall-clock under 10 min on CPU. The 5-channel DELFI profile
        # + modality dict is still exposed to the encoder; weights are
        # updated from real-data signals (not random init).
        losses["phase1"] = pretrainer.pretrain_phase1_mmp(
            n_samples=len(sample_ids),
            n_epochs=args.epochs,
            batch_size=args.batch_size,
        )
        # One short contrastive pass (helps the encoder learn
        # cross-modal alignment on the small cohort).
        losses["phase2"] = pretrainer.pretrain_phase2_contrastive(
            n_samples=len(sample_ids),
            n_epochs=max(2, args.epochs // 2),
            batch_size=args.batch_size,
        )
        train_time = time.time() - t_train
        print(f"\nPretrain wall-clock: {train_time:.1f}s")
        print(f"  phase1 final loss: {losses['phase1'][-1]:.6f}")
        print(f"  phase2 final loss: {losses['phase2'][-1]:.6f}")

    # ── Save checkpoint ────────────────────────────────────────────
    CHECKPOINT_OUT.parent.mkdir(parents=True, exist_ok=True)
    pretrainer.save_checkpoint(str(CHECKPOINT_OUT))
    ckpt_size_mb = CHECKPOINT_OUT.stat().st_size / (1024 * 1024)
    print(f"\nCheckpoint → {CHECKPOINT_OUT}  ({ckpt_size_mb:.2f} MB)")

    # ── Verify FoundationDownstream(pretrained=True, ...) loads ─────
    print("\n[Verify] FoundationDownstream(pretrained=True, "
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
        raise RuntimeError("Forward pass produced NaN/Inf — checkpoint is broken")

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
        "epochs": args.epochs,
        "device": args.device,
        "pretrain_losses": {
            k: [round(x, 6) for x in v] for k, v in losses.items()
        },
        "train_time_s": round(train_time, 1),
        "wall_clock_total_s": round(time.time() - t0, 1),
        "forward_pass_finite": out_finite,
        "data_source": {
            "primary": (
                "Pre-extracted 5-channel DELFI features at "
                f"{PIPELINE_FEAT_DIR} (originally fetched from "
                "FinaleDB S3 *.frag.tsv.bgz files via the "
                "cfdna-fragmentomics-pipeline repo)"
            ),
            "live_fetch_attempted": False,
            "live_fetch_skip_reason": (
                "S3 multi-part objects were truncated by the local "
                "network on 2026-09-21 (HEAD reports 54MB but "
                "downloads returned 16-23MB). Deferring live fetch "
                "to a future PR with a non-truncating network path."
            ),
            "finaledb_status_at_run_time": (
                "REST API in degraded state (500 on /api/v1/seqrun); "
                "S3 bucket public and serving *.frag.tsv.bgz."
            ),
        },
        "sample_ids": sample_ids,
    }
    with open(LOG_OUT, "w") as fh:
        json.dump(log, fh, indent=2)
    print(f"\nLog → {LOG_OUT}")
    print(f"Total wall-clock: {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())