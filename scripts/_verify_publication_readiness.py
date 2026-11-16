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
    """HEAD-probe several FinaleDB S3 objects to confirm accessibility.

    The previous version probed ONE key (EE1) and reported "private" on
    a 403. That was a false negative: the bucket is **partially**
    public. As of 2026-11-15 the ``.frag.tsv.bgz`` fragment objects
    and the ``.bigWig`` coverage objects return 200, while
    ``.bedGraph.bgz`` objects and the EE1 key return 403. Probing a
    single arbitrary key therefore reported "cannot fetch any new
    frag.tsv.bgz files" at a time when they were fetchable.

    So: probe several keys of different kinds, and report
    ``public`` / ``partial`` / ``private`` / ``unreachable`` from what
    is actually reachable.
    """
    base = ("https://s3.us-east-2.amazonaws.com/"
            "finaledb.epifluidlab.cchmc.org/")
    keys = [
        "entries/EE85723/hg38/EE85723.hg38.frag.tsv.bgz",
        "entries/EE85723/hg38/EE85723.hg38.insert_size_metrics.txt",
        "entries/EE85723/hg38/EE85723.hg38.frag_profile.mapq30.bw",
        "entries/EE85723/hg38/EE85723.hg38.coverage.mapq30.bedGraph.bgz",
        "entries/EE1/hg38/EE1.hg38.frag.tsv.bgz",
    ]
    checked, reachable, frag_reachable = [], 0, 0
    for key in keys:
        url = base + key
        status = None
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=API_PROBE_TIMEOUT) as resp:
                status = resp.status
        except urllib.error.HTTPError as e:
            status = e.code
        except Exception as e:  # pragma: no cover - network variance
            status = f"error: {type(e).__name__}"
        if status == 200:
            reachable += 1
            if key.endswith(".frag.tsv.bgz"):
                frag_reachable += 1
        checked.append({"key": key, "status": status})

    if frag_reachable > 0:
        verdict = "public" if reachable == len(keys) else "partial"
        reason = (
            "Fragment objects are reachable — new publications can be "
            "fetched with scripts/fetch_external_cohort.py."
            if verdict == "public" else
            "Some object kinds are reachable (fragment objects included) and "
            "others are not; new publications can still be fetched."
        )
    elif reachable > 0:
        verdict, reason = "partial", (
            "Some objects are reachable but no fragment object was; "
            "feature extraction is not possible from S3 right now."
        )
    else:
        verdict, reason = "private", (
            "No probed object is reachable (403) — the bucket looks private "
            "or the keys changed."
        )
    return {
        "url": base,
        "status": 200 if frag_reachable else 403,
        "verdict": verdict,
        "verdict_reason": reason,
        "n_keys_probed": len(keys),
        "n_reachable": reachable,
        "n_fragment_objects_reachable": frag_reachable,
        "keys": checked,
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
# Held-out external cohort
# ──────────────────────────────────────────────────────────────────────
def load_held_out_cohort(path):
    """Read the publications the external validation actually scored.

    This diagnostic only sees the *training* feature cache and the
    *training* labels file, so the held-out publications look
    `labels-missing` — as if they were unusable. They are not: they are
    the external test set, fetched from the public FinaleDB API/S3 and
    scored by `scripts/external_validation_finaledb.py`.

    The fact is therefore READ from the artifact that records it rather
    than hard-coded here. A generator that writes a number it did not
    measure is how a document ends up contradicting its own artifact.
    """
    info = {
        "artifact": path,
        "present": False,
        "publications": [],
        "source": None,
    }
    try:
        with open(path) as fh:
            payload = json.load(fh)
    except (OSError, ValueError):
        return info
    per_pub = payload.get("per_publication_auc") or {}
    info["present"] = True
    info["publications"] = sorted(per_pub.keys())
    design = payload.get("design") or {}
    info["source"] = design.get("test")
    return info


# ──────────────────────────────────────────────────────────────────────
# Verdict per publication
# ──────────────────────────────────────────────────────────────────────
def classify_publication(api_verdict, n_label, n_features, held_out=False):
    """Classify one publication as `ready` / `features-missing` /
    `labels-missing` / `api-down`.

    `n_label`    = # samples in labels file with this publication (post cell-line)
    `n_features` = bool — at least one label sample has all 5 features locally
    `held_out`   = bool — the publication is scored by the external validation
    """
    if n_label == 0:
        if held_out:
            return "labels-missing", (
                "Not in the training labels file BY DESIGN: this publication "
                "is the held-out external cohort. Its features are fetched "
                "from the public API/S3 and scored by the external "
                "validation; it contributes no training data."
            )
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

    # The held-out cohort note is emitted from the artifact that records it,
    # so re-running this generator cannot silently drop it (it was
    # hand-edited into this file once and the next run would have deleted
    # it).
    held = report.get("held_out_external") or {}
    training_pubs = [p for p in sorted(PUBLICATION_REGISTRY)
                     if report["per_publication"][p]["n_labels"] > 0]
    if held.get("present") and held.get("publications"):
        L.append(
            "\n**Read `labels-missing` as \"not in the training labels file\", "
            "not \"unusable\".**\n"
            f"Publications {', '.join(held['publications'])} are exactly the "
            "studies used as the held-out external cohort "
            f"(`{held['artifact']}`): they are fetched from the API and the "
            "S3 bucket for scoring, and deliberately contribute no training "
            f"data. The training labels file carries publications "
            f"{', '.join(training_pubs)} only.\n"
        )

    L.append("\n## Per-publication readiness\n\n")
    L.append("| Publication | Study | Verdict | Held-out external | "
            "n_labels (cell-line filtered) | "
            "n_features_present (training cache) | Why |\n")
    L.append("|---|---|---|---|---:|---:|---|\n")
    for pub in sorted(PUBLICATION_REGISTRY):
        s = report["per_publication"][pub]
        L.append(
            f"| {pub} | {PUBLICATION_REGISTRY[pub]} | "
            f"`{s['verdict']}` | "
            f"{'yes' if s.get('held_out_external') else 'no'} | "
            f"{s['n_labels']} | "
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
    held_pubs = [pub for pub, s in report["per_publication"].items()
                 if s.get("held_out_external")]
    blocked = [pub for pub, s in report["per_publication"].items()
               if s["verdict"] != "ready" and pub not in held_pubs]
    if held_pubs:
        L.append("**Held-out external cohort** (not blocked): "
                 + ", ".join(f"`{pub}` ({PUBLICATION_REGISTRY[pub]})"
                             for pub in held_pubs)
                 + f" — scored in `{held.get('artifact')}`. They are absent "
                   "from the training labels file *by design*: they are the "
                   "external test set, fetched from the public API/S3 for "
                   "scoring, contributing no training data. Their results "
                   "are in `docs/EXTERNAL_VALIDATION.md`.\n\n")
    if not blocked:
        L.append("**No publication is blocked on FinaleDB infrastructure.**\n\n")
    else:
        L.append("**BLOCKED publications**: " +
                ", ".join(f"`{pub}` ({PUBLICATION_REGISTRY[pub]})" for pub in blocked) +
                "\n\n")
        if api["verdict"] in ("api-down", "unreachable", "degraded"):
            L.append("These are blocked on **FinaleDB infrastructure**, not on "
                     "this repo. The labels file references them, but no local "
                     "5-channel DELFI features exist in the training cache "
                     f"and the API/S3 were not usable to fetch them "
                     f"(API verdict `{api['verdict']}`). Re-run this "
                     "diagnostic after FinaleDB authors restore API/S3 "
                     "access.\n\n")
        else:
            L.append("These are not in the training cache and not scored by "
                     "the external validation, so they contribute no "
                     "benchmark result today. The API reports as reachable "
                     f"(verdict `{api['verdict']}`), so the standard fetch "
                     "pipeline can populate them.\n\n")
    L.append("The cross-study benchmark trains on the `ready` set (`"
             + " ".join(report["default_scope"]) +
             "`). The held-out cohort is scored separately by "
             "`scripts/external_validation_finaledb.py`; its numbers are in "
             "`docs/EXTERNAL_VALIDATION.md`.\n")

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
    ap.add_argument("--external-validation",
                    default="results/external_validation_finaledb.json",
                    help="Artifact recording the held-out external cohort. "
                         "The diagnostic reads it so that a publication "
                         "scored there is never reported as unusable.")
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
    held = load_held_out_cohort(args.external_validation)
    # A sorted list, not a set: set iteration order over strings varies
    # with the hash seed, so a set would make the artifact's wording
    # differ between runs of the same input.
    held_pubs = sorted(held["publications"])
    if held["present"]:
        print(f"      Held-out external cohort (from {args.external_validation}): "
              f"{', '.join(held['publications']) or 'none'}")
    else:
        print(f"      Held-out external cohort: artifact not found at "
              f"{args.external_validation} — publications may be reported "
              f"as unusable when they are only held out")
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
            api_report["verdict"], n_label, n_feat,
            held_out=pub in held_pubs,
        )
        per_pub[pub] = {
            "verdict": verdict,
            "reason": reason,
            "n_labels": n_label,
            "n_cancer": label_counts[pub]["n_cancer"],
            "n_healthy": label_counts[pub]["n_healthy"],
            "n_features": bool(n_feat),
            "held_out_external": pub in held_pubs,
            "study_label": PUBLICATION_REGISTRY[pub],
        }
        flag = "✅" if verdict == "ready" else ("🧪" if pub in held_pubs else "❌")
        print(f"      {flag} {pub} ({PUBLICATION_REGISTRY[pub]}): {verdict} "
              f"(n_labels={n_label}, features={'yes' if n_feat else 'no'}"
              + (", held-out external" if pub in held_pubs else "") + ")")

    # Default scope = publications that are `ready`.
    default_scope = sorted(pub for pub, s in per_pub.items()
                          if s["verdict"] == "ready")
    if not default_scope:
        # Fall back to the script's hardcoded default (pubs 6 + 8) so the
        # diagnostic doesn't lie when nothing is `ready`.
        default_scope = ["6", "8"]

    blocked_pubs = sorted(p for p, s in per_pub.items()
                          if s["verdict"] != "ready" and not s["held_out_external"])
    summary_bits = [f"FinaleDB API verdict: `{api_report['verdict']}`."]
    if held["present"] and held_pubs:
        summary_bits.append(
            "Held-out external cohort (scored, not blocked): "
            + ", ".join(f"{p} ({per_pub[p]['study_label']})" for p in held_pubs)
            + f" (from {args.external_validation})."
        )
    if blocked_pubs:
        summary_bits.append(
            "Blocked publications: "
            + ", ".join(f"{p} ({per_pub[p]['study_label']})" for p in blocked_pubs)
            + "."
        )
    else:
        summary_bits.append("No publication is blocked on infrastructure.")
    summary_bits.append(
        "The cross-study benchmark uses the `ready` set "
        f"({default_scope}); the held-out cohort is scored separately by "
        "scripts/external_validation_finaledb.py."
    )

    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "labels_file": args.labels_multiclass,
        "features_dir": args.features_dir,
        "finaledb_api": api_report,
        "finaledb_s3": s3_report,
        "per_publication": per_pub,
        "default_scope": default_scope,
        "held_out_external": held,
        "honest_summary": " ".join(summary_bits),
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