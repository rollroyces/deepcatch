#!/usr/bin/env python3
"""
sens_at_spec ablation: paired t-test of per-seed FoundationDownstream
AUCs under loss="ce" vs loss="sens_at_spec" on the 20-patient TCGA-LUAD
panel at TF=0.1%.

Reads the per-seed arrays from the two smoke JSONs (CE and sens_at_spec)
and emits the head-to-head comparison with the paired t-test, mean diff,
std diff, and a CI on the difference.

Output: results/sens_at_spec_ablation.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))


def _paired_t(a: np.ndarray, b: np.ndarray) -> Dict[str, float]:
    """Two-sided paired t-test. Returns mean_diff, std_diff, t, df, p, ci95."""
    from scipy import stats

    d = b - a
    n = len(d)
    mean_d = float(np.mean(d))
    std_d = float(np.std(d, ddof=1))
    se_d = std_d / math.sqrt(n) if n > 0 else float("nan")
    t_stat = mean_d / se_d if se_d > 0 else float("nan")
    df = n - 1
    p = float(2 * (1 - stats.t.cdf(abs(t_stat), df=df))) if n > 1 else float("nan")
    t_crit = float(stats.t.ppf(0.975, df=df)) if df > 0 else float("nan")
    ci_lo = mean_d - t_crit * se_d
    ci_hi = mean_d + t_crit * se_d
    return {
        "n": int(n),
        "mean_diff": mean_d,
        "std_diff": std_d,
        "se_diff": float(se_d),
        "t_stat": float(t_stat),
        "df": int(df),
        "p_value": float(p),
        "ci95_lo": float(ci_lo),
        "ci95_hi": float(ci_hi),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--ce-json",
        default=str(_ROOT / "results" / "sens_at_spec_ce.json"),
    )
    ap.add_argument(
        "--sens-json",
        default=str(_ROOT / "results" / "sens_at_spec_sens.json"),
    )
    ap.add_argument(
        "--out",
        default=str(_ROOT / "results" / "sens_at_spec_ablation.json"),
    )
    args = ap.parse_args()

    with open(args.ce_json) as f:
        ce = json.load(f)
    with open(args.sens_json) as f:
        sens = json.load(f)

    # Per-seed arrays must come from the same 5 seeds, run on the same
    # 20 patients, with the same shuffle (StratifiedKFold(random_state=seed)).
    # The synthetic frag channel is also seeded by `seed + 9999`, so the
    # two runs share the per-seed panel_scores and frag_scores — the
    # ONLY difference is the loss passed to FoundationDownstream.
    ce_seed = list(range(ce["seeds"]))
    sens_seed = list(range(sens["seeds"]))
    assert ce_seed == sens_seed, "CE and sens_at_spec runs must use the same seed set"

    n_seeds = ce["seeds"]
    out: Dict[str, object] = {
        "n_seeds": n_seeds,
        "n_patients": ce["n_samples"] // 2,  # 40 samples = 20 cancer + 20 control
        "n_samples_total": ce["n_samples"],
        "n_cancer": ce["n_cancer"],
        "n_healthy": ce["n_healthy"],
        "tumor_fraction": 0.001,
        "alpha_pos_sens": sens["alpha_pos"],
        "ce_loss_name": "ce",
        "sens_loss_name": "sens_at_spec",
        "data_source": ce["data_source"],
        "ce_run_path": str(args.ce_json),
        "sens_run_path": str(args.sens_json),
    }

    # Primary: foundation AUC (the main metric the loss is designed to
    # improve at the high-specificity operating point).
    ce_fnd = np.asarray(ce["foundation_aucs"], dtype=np.float64)
    sens_fnd = np.asarray(sens["foundation_aucs"], dtype=np.float64)
    out["foundation_aucs_ce"] = ce_fnd.tolist()
    out["foundation_aucs_sens_at_spec"] = sens_fnd.tolist()
    out["foundation_auc_mean_ce"] = float(np.mean(ce_fnd))
    out["foundation_auc_mean_sens_at_spec"] = float(np.mean(sens_fnd))
    out["foundation_auc_std_ce"] = float(np.std(ce_fnd, ddof=1))
    out["foundation_auc_std_sens_at_spec"] = float(np.std(sens_fnd, ddof=1))
    out["foundation_auc_paired_t"] = _paired_t(ce_fnd, sens_fnd)

    # Secondary: sens@99 (the actual metric the focal loss is designed
    # to improve). This is the headline MRD operating point.
    ce_s99 = np.asarray(ce["foundation_sens_at_99"], dtype=np.float64)
    sens_s99 = np.asarray(sens["foundation_sens_at_99"], dtype=np.float64)
    out["foundation_sens_at_99_ce"] = ce_s99.tolist()
    out["foundation_sens_at_99_sens_at_spec"] = sens_s99.tolist()
    out["foundation_sens_at_99_mean_ce"] = float(np.mean(ce_s99))
    out["foundation_sens_at_99_mean_sens_at_spec"] = float(np.mean(sens_s99))
    out["foundation_sens_at_99_paired_t"] = _paired_t(ce_s99, sens_s99)

    # Tertiary: sens@95
    ce_s95 = np.asarray(ce["foundation_sens_at_95"], dtype=np.float64)
    sens_s95 = np.asarray(sens["foundation_sens_at_95"], dtype=np.float64)
    out["foundation_sens_at_95_ce"] = ce_s95.tolist()
    out["foundation_sens_at_95_sens_at_spec"] = sens_s95.tolist()
    out["foundation_sens_at_95_mean_ce"] = float(np.mean(ce_s95))
    out["foundation_sens_at_95_mean_sens_at_spec"] = float(np.mean(sens_s95))
    out["foundation_sens_at_95_paired_t"] = _paired_t(ce_s95, sens_s95)

    # Sanity check: the LR baseline (sklearn LR on [panel, frag]) does
    # NOT take the loss flag — should be identical across runs since
    # panel_scores + frag_scores are seeded the same way. If they
    # differ, something else is contaminating.
    ce_lr = np.asarray(ce["lr_baseline_aucs"], dtype=np.float64)
    sens_lr = np.asarray(sens["lr_baseline_aucs"], dtype=np.float64)
    out["lr_baseline_aucs_ce"] = ce_lr.tolist()
    out["lr_baseline_aucs_sens_at_spec"] = sens_lr.tolist()
    out["lr_baseline_aucs_match"] = bool(np.allclose(ce_lr, sens_lr, atol=1e-9))

    # Honest classification of the lift.
    p = out["foundation_auc_paired_t"]["p_value"]
    diff = out["foundation_auc_paired_t"]["mean_diff"]
    if abs(diff) < 0.005:
        verdict_auc = "null"
    elif diff > 0:
        verdict_auc = "positive"
    else:
        verdict_auc = "negative"

    s99_diff = out["foundation_sens_at_99_paired_t"]["mean_diff"]
    s99_p = out["foundation_sens_at_99_paired_t"]["p_value"]
    if abs(s99_diff) < 0.05:
        verdict_s99 = "null"
    elif s99_diff > 0:
        verdict_s99 = "positive"
    else:
        verdict_s99 = "negative"

    out["verdict_auc"] = verdict_auc
    out["verdict_sens_at_99"] = verdict_s99
    out["verdict_summary"] = {
        "auc": (
            f"{verdict_auc} (Δ {diff:+.4f} AUC, p={p:.3f})"
        ),
        "sens_at_99": (
            f"{verdict_s99} (Δ {s99_diff:+.3f} sens@99, p={s99_p:.3f})"
        ),
        "smoke_gate_pass_ce": ce["gate_pass"],
        "smoke_gate_pass_sens_at_spec": sens["gate_pass"],
        "signal_to_artifact_ce": ce.get("signal_to_artifact_ratio"),
        "signal_to_artifact_sens_at_spec": sens.get("signal_to_artifact_ratio"),
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[ablation] wrote {args.out}")
    print(json.dumps(out["verdict_summary"], indent=2))
    print(f"\nfoundation AUC: CE={np.mean(ce_fnd):.4f} ± {np.std(ce_fnd, ddof=1):.4f}, "
          f"sens_at_spec={np.mean(sens_fnd):.4f} ± {np.std(sens_fnd, ddof=1):.4f}")
    print(f"  paired Δ = {diff:+.4f}, 95% CI [{out['foundation_auc_paired_t']['ci95_lo']:+.4f}, "
          f"{out['foundation_auc_paired_t']['ci95_hi']:+.4f}], "
          f"t={out['foundation_auc_paired_t']['t_stat']:.3f}, p={p:.4f}")
    print(f"foundation sens@99: CE={np.mean(ce_s99):.3f} ± {np.std(ce_s99, ddof=1):.3f}, "
          f"sens_at_spec={np.mean(sens_s99):.3f} ± {np.std(sens_s99, ddof=1):.3f}")
    print(f"  paired Δ = {s99_diff:+.3f}, 95% CI [{out['foundation_sens_at_99_paired_t']['ci95_lo']:+.3f}, "
          f"{out['foundation_sens_at_99_paired_t']['ci95_hi']:+.3f}], "
          f"t={out['foundation_sens_at_99_paired_t']['t_stat']:.3f}, p={s99_p:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
