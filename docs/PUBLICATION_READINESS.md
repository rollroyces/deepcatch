# Publication Readiness (FinaleDB Open-Data Benchmark)
- Generated: `2026-09-23T15:33:49.019966+00:00`
- Labels file: `/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv`
- Features dir: `/Users/hermes/cfdna-fragmentomics-pipeline/data/features`

## FinaleDB infrastructure status

**API verdict**: `api-down` — /misc responds but /publication + /seqrun + /summary all return 500 — Postgres connection lost. No way to enumerate seqruns or publications; cannot refetch features for publications missing from the local cache.

| Endpoint | URL | Status |
|---|---|---:|
| `/api/v1/misc` | `<http://finaledb.research.cchmc.org/api/v1/misc>` | 200 |
| `/api/v1/publication` | `<http://finaledb.research.cchmc.org/api/v1/publication>` | 500 |
| `/api/v1/seqrun` | `<http://finaledb.research.cchmc.org/api/v1/seqrun>` | 500 |
| `/api/v1/summary` | `<http://finaledb.research.cchmc.org/api/v1/summary>` | 500 |

**S3 verdict**: `private` — S3 bucket returns 403 — bucket was made private. Cannot fetch any new frag.tsv.bgz files.

## Per-publication readiness

| Publication | Study | Verdict | n_labels (cell-line filtered) | n_features_present | Why |
|---|---|---|---:|---:|---|
| 1 | snyder | `labels-missing` | 0 | no | No samples in the labels file reference this publication. |
| 6 | jiang | `ready` | 121 | yes | All artifacts present; benchmark can include this publication. |
| 7 | sun | `labels-missing` | 0 | no | No samples in the labels file reference this publication. |
| 8 | cristiano | `ready` | 537 | yes | All artifacts present; benchmark can include this publication. |
| 9 | adalsteinsson | `labels-missing` | 0 | no | No samples in the labels file reference this publication. |

## Default benchmark scope

- **Default publications**: `['6', '8']` — the only publications with `ready` status as of `2026-09-23T15:33:49.019966+00:00`.
- `--publications` override: pass any subset. Missing publications are reported via stderr warnings, NOT hard-failed (the script skips them gracefully).
- `--include-snyder` / `--include-sun`: shortcuts for publications 1 / 7.

## Honest reading

**BLOCKED publications**: `1` (snyder), `7` (sun), `9` (adalsteinsson)

These are blocked on **FinaleDB infrastructure**, not on this repo. The labels file references them, but no local 5-channel DELFI features exist (and the API/S3 are not reachable to fetch them). Re-run this diagnostic after FinaleDB authors restore API/S3 access.

Until then, the cross-study benchmark uses `6 8` (the only publications with `ready` status). Numbers from the `--publications 6 8` run are the honest open-data cross-study benchmark.
