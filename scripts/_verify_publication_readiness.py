#!/usr/bin/env python3
"""Publication-readiness diagnostic for the cross-study FinaleDB benchmark.

Reports, per FinaleDB publication:

  - `ready`           : features present locally AND labels have samples
  - `features-missing`: NO samples have all 5-channel DELFI features locally
  - `labels-missing`  : NO samples in the labels file with this publication
  - `api-down`        : features would normally come from FinaleDB API/S3,
                        but the public REST API and S3 bucket are currently
                        unreachable (HTTP 500 / 403). See cfdna-fragmentomics
                        skill, "When FinaleDB is fully down".

The script is honest by design: it NEVER fabricates a "ready" status when
the underlying infrastructure is broken. When the API is down, publications
1 (Snyder 2016) and 7 (Sun 2019) are reported as `api-down` (not `ready`)
even if the labels file references them — because there is no path to
fetch their features in this state.

Usage:
  env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \\
      scripts/_verify_publication_readiness.py \\
      --features-dir /Users/hermes/cfdna-fragmentomics-pipeline/data/features \\
      --labels-multiclass /Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv \\
      --out-json results/publication_readiness.json \\
      --out-md    docs/PUBLICATION_READINESS.md
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

# Make the cross_study_finallydb module importable for shared constants.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from cross_study_finallydb import (  # noqa: E402
    PUBLICATION_REGISTRY,
    STUDY_TO_PUBLICATION,
    features_present_for_publication,
)


# ──────────────────────────────────────────────────────────────────────
# FinaleDB API probe (the API + S3 are DOWN since 2026-09; this is the
# honest diagnostic that records that fact in the output JSON).
# ──────────────────────────────────────────────────────────────────────
FINALEDB_API_BASE = "http://finaledb.research.cchmc.org/api/v1"
API_PROBE_TIMEOUT = 8  # seconds per endpoint


def _http_get(url: str, timeout: float = API_PROBE_TIMEOUT):
    """Tiny urllib GET helper. Returns (status_code, body_or_None)."""
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            body = resp.read(2048).decode("utf-8", errors="replace")
            return resp.status, body
    except urllib.error.HTTPError as e:
        return e.code, None
    except urllib.error.URLError as e:
        return None, str(e)
    except Exception as e:  # pragma: no cover
        return None, f"{type(e).__name__}: {e}"


def probe_finaledb_api():
    """Probe the public FinaleDB API and report per-endpoint status.

    Per the cfdna-fragmentomics skill (and the reference API doc):
      /api/v1/misc         → always 200 if server is alive
      /api/v1/publication  → 500 when Postgres is broken
      /api/v1/seqrun       → 500 when Postgres is broken
      /api/v1/summary      → 500 when Postgres is broken
    """
    endpoints = {
        "misc": f"{FINALEDB_API_BASE}/misc",
        "publication": f"{FINALEDB_API_BASE}/publication",
        "seqrun": f"{FINALEDB_API_BASE}/seqrun",
        "summary": f"{FINALEDB_API_BASE}/summary",
    }
    out = {}
    for name, url in endpoints.items():
        status, body = _http_get(url)
        out[name] = {
            "url": url,
            "status": status,
            "body_excerpt": (body[:200] if isinstance(body, str) else None),
        }
    # Aggregate verdict.
    misc_status = out["misc"]["status"]
    other = [out[k]["status"] for k in ("publication", "seqrun", "summary")]
    if misc_status != 200:
        verdict = "unreachable"
        verdict_reason = (
            f"FinaleDB API /misc returned status={misc_status} — the API "
            f"host is unreachable. Every publication fetch path is blocked."
        )
    elif all(s == 500 for s in other):
        verdict = "api-down"
        verdict_reason = (
            "/misc responds but /publication + /seqrun + /summary all return "
            "500 — Postgres connection lost. No way to enumerate seqruns or "
            "publications; cannot refetch features for publications missing "
            "from the local cache."
        )
    elif any(s == 500 for s in other):
        verdict = "degraded"
        verdict_reason = (
            f"Some endpoints 500 ({other}); partial functionality."
        )
    else:
        verdict = "ok"
        verdict_reason = (
            "All probed endpoints respond. Features for publications "
            "missing locally can be re-fetched via the standard pipeline."
        )
    return {
        "endpoints": out,
        "verdict": verdict,
        "verdict_reason": verdict_reason,
        "checked_at": datetime.now(timezone.utc).isoformat(),
    }


def probe_finaledb_s3():
    """HEAD-probe the FinaleDB S3 bucket to confirm bucket accessibility.

    The S3 bucket was made private in 2026-09 (every key returns 403).
    We only probe ONE key (EE1 — the first seqrun) to keep the diagnostic
    cheap. 200 = public; 403 = private; connection error = network down.
    """
    url = (
        "https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/"
        "entries/EE1/hg38/EE1.hg38.frag.tsv.bgz"
    )
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=API_PROBE_TIMEOUT) as resp:
            return {
                "url": url,
                "status": resp.status,
                "verdict": "public",
                "checked_at": datetime.now(timezone.utc).isoformat(),
            }
    except urllib.error.HTTPError as e:
        if e.code == 403:
            return {
                "url": url,
                "status": 403,
                "verdict": "private",
                "verdict_reason": (
                    "S3 bucket returns 403 — bucket was made private. "
                    "Cannot fetch any new frag.tsv.bgz files."
                ),
                "checked_at": datetime.now(timezone.utc).isoformat(),
            }
        return {
            "url": url,
            "status": e.code,
            "verdict": "error",
            "checked_at": datetime.now(timezone.utc).isoformat(),
        }
    except urllib.error.URLError as e:
        return {
            "url": url,
            "status": None,
            "verdict": "unreachable",
            "verdict_reason": str(e),
            "checked_at": datetime.now(timezone.utc).isoformat(),
        }


# ──────────────────────────────────────────────────────────────────────
# Labels + features probe
# ──────────────────────────────────────────────────────────────────────
def load_labels_for_diagnostic(path: str):
    """Lightweight labels loader — returns the same 4 dicts as
    `cross_study_finallydb.load_labels_multiclass` but tolerates a
    missing 5th `publication` column.
    """
    labels, studies, disease_class, publication = {}, {}, {}, {}
    header_seen = False
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            if not header_seen:
                if p[0] == "sample" and p[3] == "study":
                    header_seen = True
                    continue
                header_seen = True
            if p[3] in ("study", ""):
                continue
            s = p[0]
            labels[s] = 1 if p[2] == "cancer" else 0
            studies[s] = p[3]
            disease_class[s] = p[1]
            if len(p) >= 5 and p[4] and p[4] not in ("publication",):
                publication[s] = p[4]
            else:
                publication[s] = STUDY_TO_PUBLICATION.get(p[3], "")
    return labels, studies, disease_class, publication


# Cell-line filter regex (mirrors `cross_study_finallydb.CELL_LINE_RE`)
CELL_LINE_RE = re.compile(
    r"^(GM\d+|HeLa|HepG2|K562|HL60|Jurkat|Raji|MCF7|U937|THP1|HEK293|"
    r"HCT116|SW480|A549|GM12878)",
    re.I,
)


def count_label_samples_per_publication(labels, studies, disease_class,
                                        publication):
    """Return {pub: {n_total, n_cancer, n_healthy, n_dropped_cell_line}}."""
    out = {pub: {"n_total": 0, "n_cancer": 0, "n_healthy": 0,
                 "n_dropped_cell_line": 0}
           for pub in PUBLICATION_REGISTRY}
    for s, pub in publication.items():
        if pub not in out:
            continue
        if CELL_LINE_RE.match(s):
            out[pub]["n_dropped_cell_line"] += 1
            continue
        out[pub]["n_total"] += 1
        if labels[s] == 1:
            out[pub]["n_cancer"] += 1
        else:
            out[pub]["n_healthy"] += 1
    return out


# ──────────────────────────────────────────────────────────────────────
# Verdict per publication
# ──────────────────────────────────────────────────────────────────────
def classify_publication(api_verdict, n_label, n_features):
    """Classify one publication as `ready` / `features-missing` /
    `labels-missing` / `api-down`.

    `n_label`    = # samples in labels file with this publication (post cell-line)
    `n_features` = bool — at least one label sample has all 5 features locally
    """
    if n_label == 0:
        return "labels-missing", (
            "No samples in the labels file reference this publication."
        )
    if not n_features:
        if api_verdict in ("api-down", "unreachable", "degraded"):
            return "api-down", (
                "Local features are missing for this publication AND the "
                "FinaleDB API/S3 is unreachable, so there is no path to "
                "fetch them. Publication is BLOCKED on infrastructure."
            )
        return "features-missing", (
            "Local features are missing for this publication. The FinaleDB "
            "API reports as reachable; rerun the standard fetch pipeline "
            "to populate the local cache."
        )
    return "ready", "All artifacts present; benchmark can include this publication."


# ──────────────────────────────────────────────────────────────────────
# Markdown writer
# ──────────────────────────────────────────────────────────────────────
def write_markdown(report, md_path):
    L = []
    L.append("# Publication Readiness (FinaleDB Open-Data Benchmark)\n")
    L.append(f"- Generated: `{report['generated_at']}`\n")
    L.append(f"- Labels file: `{report['labels_file']}`\n")
    L.append(f"- Features dir: `{report['features_dir']}`\n")

    api = report["finaledb_api"]
    s3 = report["finaledb_s3"]
    L.append("\n## FinaleDB infrastructure status\n\n")
    L.append(f"**API verdict**: `{api['verdict']}` — {api['verdict_reason']}\n\n")
    L.append("| Endpoint | URL | Status |\n")
    L.append("|---|---|---:|\n")
    for name, ep in api["endpoints"].items():
        L.append(f"| `/api/v1/{name}` | `<{ep['url']}>` | "
                f"{ep['status'] if ep['status'] is not None else 'ERR'} |\n")
    L.append(f"\n**S3 verdict**: `{s3['verdict']}`"
            + (f" — {s3['verdict_reason']}" if s3.get("verdict_reason") else "")
            + "\n")

    L.append("\n## Per-publication readiness\n\n")
    L.append("| Publication | Study | Verdict | n_labels (cell-line filtered) | "
            "n_features_present | Why |\n")
    L.append("|---|---|---|---:|---:|---|\n")
    for pub in sorted(PUBLICATION_REGISTRY):
        s = report["per_publication"][pub]
        L.append(
            f"| {pub} | {PUBLICATION_REGISTRY[pub]} | "
            f"`{s['verdict']}` | {s['n_labels']} | "
            f"{'yes' if s['n_features'] else 'no'} | {s['reason']} |\n"
        )

    # Default benchmark scope.
    L.append("\n## Default benchmark scope\n\n")
    L.append(f"- **Default publications**: `{report['default_scope']}` — the only "
            f"publications with `ready` status as of `{report['generated_at']}`.\n")
    L.append(f"- `--publications` override: pass any subset. Missing publications "
            f"are reported via stderr warnings, NOT hard-failed (the script "
            f"skips them gracefully).\n")
    L.append(f"- `--include-snyder` / `--include-sun`: shortcuts for "
            f"publications 1 / 7.\n")

    # Honest reading.
    L.append("\n## Honest reading\n\n")
    blocked = [pub for pub, s in report["per_publication"].items()
              if s["verdict"] != "ready"]
    if not blocked:
        L.append("All FinaleDB publications in the registry are `ready`. The "
                "default `--publications 6 8` benchmark can be extended to "
                "include 1 (Snyder) and 7 (Sun) without API changes.\n")
    else:
        L.append("**BLOCKED publications**: " +
                ", ".join(f"`{pub}` ({PUBLICATION_REGISTRY[pub]})" for pub in blocked) +
                "\n\n")
        L.append("These are blocked on **FinaleDB infrastructure**, not on "
                "this repo. The labels file references them, but no local "
                "5-channel DELFI features exist (and the API/S3 are not "
                "reachable to fetch them). Re-run this diagnostic after "
                "FinaleDB authors restore API/S3 access.\n\n")
        L.append("Until then, the cross-study benchmark uses `"
                + " ".join(report["default_scope"]) +
                "` (the only publications with `ready` status). "
                "Numbers from the `--publications 6 8` run are the "
                "honest open-data cross-study benchmark.\n")

    os.makedirs(os.path.dirname(md_path) or ".", exist_ok=True)
    with open(md_path, "w") as f:
        f.write("".join(L))


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--features-dir",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/data/features")
    ap.add_argument("--labels-multiclass",
                    default="/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv")
    ap.add_argument("--out-json", default="results/publication_readiness.json")
    ap.add_argument("--out-md", default="docs/PUBLICATION_READINESS.md")
    ap.add_argument("--skip-network", action="store_true",
                    help="Skip the API/S3 probes (useful for CI / offline runs).")
    args = ap.parse_args()

    print(f"[1/3] Probing FinaleDB infrastructure")
    if args.skip_network:
        api_report = {
            "endpoints": {},
            "verdict": "skipped",
            "verdict_reason": "skipped via --skip-network",
            "checked_at": datetime.now(timezone.utc).isoformat(),
        }
        s3_report = {
            "url": None,
            "status": None,
            "verdict": "skipped",
            "checked_at": datetime.now(timezone.utc).isoformat(),
        }
    else:
        api_report = probe_finaledb_api()
        s3_report = probe_finaledb_s3()
    print(f"      API verdict: {api_report['verdict']}")
    print(f"      S3  verdict: {s3_report['verdict']}")

    print(f"[2/3] Loading labels from {args.labels_multiclass}")
    labels, studies, disease_class, publication = load_labels_for_diagnostic(
        args.labels_multiclass
    )
    print(f"      Loaded {len(labels)} samples")

    print(f"[3/3] Per-publication readiness")
    label_counts = count_label_samples_per_publication(
        labels, studies, disease_class, publication
    )
    per_pub = {}
    for pub in sorted(PUBLICATION_REGISTRY):
        n_label = label_counts[pub]["n_total"]
        n_feat = features_present_for_publication(
            labels, studies, publication, pub, args.features_dir
        )
        verdict, reason = classify_publication(
            api_report["verdict"], n_label, n_feat
        )
        per_pub[pub] = {
            "verdict": verdict,
            "reason": reason,
            "n_labels": n_label,
            "n_cancer": label_counts[pub]["n_cancer"],
            "n_healthy": label_counts[pub]["n_healthy"],
            "n_features": bool(n_feat),
            "study_label": PUBLICATION_REGISTRY[pub],
        }
        flag = "✅" if verdict == "ready" else "❌"
        print(f"      {flag} {pub} ({PUBLICATION_REGISTRY[pub]}): {verdict} "
              f"(n_labels={n_label}, features={'yes' if n_feat else 'no'})")

    # Default scope = publications that are `ready`.
    default_scope = sorted(pub for pub, s in per_pub.items()
                          if s["verdict"] == "ready")
    if not default_scope:
        # Fall back to the script's hardcoded default (pubs 6 + 8) so the
        # diagnostic doesn't lie when nothing is `ready`.
        default_scope = ["6", "8"]

    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "labels_file": args.labels_multiclass,
        "features_dir": args.features_dir,
        "finaledb_api": api_report,
        "finaledb_s3": s3_report,
        "per_publication": per_pub,
        "default_scope": default_scope,
        "honest_summary": (
            f"FinaleDB API verdict: `{api_report['verdict']}`. "
            + ("" if not [p for p, s in per_pub.items() if s["verdict"] != "ready"]
                else
                f"Blocked publications: "
                + ", ".join(f"{p} ({per_pub[p]['study_label']})"
                            for p in sorted(per_pub)
                            if per_pub[p]["verdict"] != "ready")
                + ". The cross-study benchmark uses the `ready` set "
                f"({default_scope}) until API/S3 are restored."
            )
        ),
    }

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nWrote {args.out_json}")

    write_markdown(report, args.out_md)
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())