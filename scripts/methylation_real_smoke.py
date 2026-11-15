#!/usr/bin/env python3
"""Run ``src/methylation_gnn`` on real FinaleMe-derived methylation data.

Until now the methylation layer was architecture-only: every number it
produced came from synthetic beta values. This script closes that gap
by taking the per-CpG methylation calls that FinaleMe decodes from real
public fragment data and pushing them through the actual
``RegulatoryGraphBuilder`` -> ``MethylationGNN`` path.

Input
-----
A FinaleMe decoder output (``*.decoded.txt.gz``), whose columns are

    #chr start end methy_perc_predict methy_count_predict
    total_count_predict methy_perc_obs methy_count_obs total_count_obs

Output
------
``results/methylation_real_smoke.json``:
  * the real beta summary (n CpGs, mean/median, methylation level);
  * the graph the builder produced from those real betas;
  * the shapes a forward pass through the real GNN produced.

No synthetic beta values are generated anywhere in this file.

Usage:
    python scripts/methylation_real_smoke.py \
        --decoded /path/to/sample.decoded.txt.gz \
        --label 1 --sample-name MBC_1330_1 \
        --out-json results/methylation_real_smoke.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)


def parse_decoded(path: str):
    """Parse a FinaleMe decoder output into (positions, betas)."""
    positions = []
    betas = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                start = int(parts[1])
                pct = float(parts[3])
            except ValueError:
                continue
            positions.append((parts[0], start))
            betas.append(pct / 100.0)
    return positions, np.asarray(betas, dtype=float)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--decoded", required=True,
                    help="FinaleMe decoder output (.decoded.txt.gz)")
    ap.add_argument("--sample-name", default="sample")
    ap.add_argument("--label", type=int, default=0)
    ap.add_argument("--disease", default="unknown")
    ap.add_argument("--max-nodes", type=int, default=8000,
                    help="cap on the number of CpG nodes (the graph is "
                         "quadratic in the co-methylation edge computation)")
    ap.add_argument("--out-json", default="results/methylation_real_smoke.json")
    args = ap.parse_args()

    positions, betas = parse_decoded(args.decoded)
    if betas.size == 0:
        print(f"no beta values parsed from {args.decoded}", file=sys.stderr)
        return 1

    summary = {
        "n_cpgs": int(betas.size),
        "beta_mean": round(float(betas.mean()), 4),
        "beta_median": round(float(np.median(betas)), 4),
        "beta_std": round(float(betas.std()), 4),
        "fraction_methylated_below_0.5": round(float((betas < 0.5).mean()), 4),
        "n_chromosomes": len({c for c, _ in positions}),
    }
    print(f"[1/3] parsed {summary['n_cpgs']} real CpG beta values "
          f"(mean {summary['beta_mean']})")

    print("[2/3] building the regulatory graph from the real betas ...", flush=True)
    import torch
    from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder
    from src.methylation_gnn.gnn_model import MethylationGNN

    idx = np.arange(betas.size)
    if idx.size > args.max_nodes:
        idx = np.linspace(0, idx.size - 1, args.max_nodes).astype(int)
    beta_subset = betas[idx]

    builder = RegulatoryGraphBuilder(
        n_nodes=int(beta_subset.size), edge_k=20, reference_genome="hg19")
    graph = builder.build_graph(
        sample_name=args.sample_name,
        methylation_data={"beta_values": beta_subset},
        label=args.label,
        cancer_type=args.disease,
    )

    n_edges = int(graph.edge_index.shape[1]) if hasattr(graph, "edge_index") else 0
    graph_summary = {
        "n_nodes": int(graph.x.shape[0]) if hasattr(graph, "x") else None,
        "node_feature_dim": int(graph.x.shape[1]) if hasattr(graph, "x") else None,
        "n_edges": n_edges,
        "edge_types": int(graph.edge_type.shape[0])
        if hasattr(graph, "edge_type") and graph.edge_type is not None else None,
    }

    print("[3/3] forward pass through the real GNN ...", flush=True)
    model = MethylationGNN(n_node_features=graph_summary["node_feature_dim"] or 20)
    with torch.no_grad():
        out = model(graph.x, graph.edge_index,
                    getattr(graph, "edge_type", None))
    forward = {k: list(v.shape) for k, v in out.items() if hasattr(v, "shape")}

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {
            "decoded_file": os.path.abspath(args.decoded),
            "origin": "FinaleMe decode of real FinaleDB fragment data "
                      "(see docs/METHYLATION_PIPELINE.md)",
            "synthetic_beta_used": False,
        },
        "sample": {
            "name": args.sample_name,
            "label": args.label,
            "disease": args.disease,
        },
        "beta_summary": summary,
        "graph": graph_summary,
        "graph_structure_is_synthetic": True,
        "forward_pass_output_shapes": forward,
        "honest_framing": (
            "The methylation VALUES are real: they are decoded by FinaleMe "
            "from real public FinaleDB fragment data, and no synthetic beta "
            "is generated anywhere in this file. The graph TOPOLOGY is not: "
            "no reference region set or Hi-C contact map is loaded, so "
            "RegulatoryGraphBuilder falls back to synthetic region "
            "annotations and random edges -- the builder logs that as "
            "'No reference regions loaded' and 'No edges produced'. This is "
            "therefore a demonstration that the module consumes real "
            "methylation data end to end, not a performance result: the "
            "GNN is untrained, and a channel that discriminates cancer from "
            "non-cancer needs a cohort-trained HMM, a trained GNN, and a "
            "reference region set."),
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"wrote {args.out_json}")
    print(f"  beta summary : {summary}")
    print(f"  graph        : {graph_summary}")
    print(f"  forward      : {forward}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
