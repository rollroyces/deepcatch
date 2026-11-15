#!/usr/bin/env python3
"""Print the benchmark numbers that the manuscript and docs quote.

One place to read them from, so the prose and the artifacts cannot drift
apart: the script reads the committed ``results/*.json`` and prints every
figure the documents quote, with the artifact it came from.
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(name: str):
    path = os.path.join(ROOT, "results", name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def line(label, value, source):
    print(f"  {label:<46} {value}   [{source}]")


def main() -> int:
    gc = load("cross_study_finallydb_gc_corrected.json")
    legacy = load("cross_study_finallydb.json")
    pc = load("per_cancer_sens_at_spec.json")

    for tag, art, src in (("GC-corrected (headline)", gc,
                           "cross_study_finallydb_gc_corrected.json"),
                          ("legacy (no GC correction)", legacy,
                           "cross_study_finallydb.json")):
        if not art:
            continue
        print(f"== {tag}")
        pooled = art["pooled"]["harmonized"]
        line("pooled harmonized AUC",
             f"{pooled['auc_mean']:.4f} +- {pooled['auc_std']:.4f}", src)
        line("pooled no_harmonize AUC",
             f"{art['pooled']['no_harmonize']['auc_mean']:.4f}", src)
        line("sens@95 / @98 / @99",
             " / ".join(f"{pooled[f'sens_at_spec_{s}']:.3f}"
                        for s in (95, 98, 99)), src)
        for pub, row in (art.get("per_publication")
                         or art.get("per_cohort", {})).items():
            name = row.get("study", pub)
            line(f"per-publication {name} (pub {pub})",
                 f"{row['auc_mean']:.4f} +- {row['auc_std']:.4f}", src)
        tc = (art.get("true_confound_control", {})
              .get("cancer_jiang_healthy_cristiano", {}))
        if tc and isinstance(tc.get("harmonized"), dict):
            line("true-confound harm / no-harm",
                 f"{tc['harmonized']['auc_mean']:.4f} / "
                 f"{tc['no_harmonize']['auc_mean']:.4f}", src)
        print()

    if pc:
        print("== per-cancer standalone")
        line("pooled AUC (DeLong CI)",
             f"{pc['pooled']['auc_mean']:.4f} "
             f"[{pc['pooled']['auc_ci'][0]:.4f}, {pc['pooled']['auc_ci'][1]:.4f}]",
             "per_cancer_sens_at_spec.json")
        for cancer, row in sorted(pc.get("per_cancer", {}).items()):
            line(f"{cancer} AUC",
                 f"{row['auc_mean']:.4f} "
                 f"[{row['auc_ci'][0]:.4f}, {row['auc_ci'][1]:.4f}]",
                 "per_cancer_sens_at_spec.json")
            s = row.get("sens_at_spec", {})
            if s:
                line(f"{cancer} sens@95/98/99",
                     " / ".join(f"{s[k].get('sensitivity', float('nan')):.3f}"
                                for k in ("0.95", "0.98", "0.99") if k in s),
                     "per_cancer_sens_at_spec.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
