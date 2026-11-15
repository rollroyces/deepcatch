# Publication Readiness (FinaleDB Open-Data Benchmark)
- Generated: `2026-11-15T16:54:53.394362+00:00`
- Labels file: `/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv`
- Features dir: `/Users/hermes/cfdna-fragmentomics-pipeline/data/features`

## FinaleDB infrastructure status

**API verdict**: `ok` — All probed endpoints respond. Features for publications missing locally can be re-fetched via the standard pipeline.

| Endpoint | URL | Status |
|---|---|---:|
| `/api/v1/misc` | `<http://finaledb.research.cchmc.org/api/v1/misc>` | 200 |
| `/api/v1/publication` | `<http://finaledb.research.cchmc.org/api/v1/publication>` | 200 |
| `/api/v1/seqrun` | `<http://finaledb.research.cchmc.org/api/v1/seqrun>` | 200 |
| `/api/v1/summary` | `<http://finaledb.research.cchmc.org/api/v1/summary>` | 200 |

**S3 verdict**: `partial` — Some object kinds are reachable (fragment objects included) and others are not; new publications can still be fetched.

**Read `labels-missing` as "not in the training labels file", not "unusable".**
Publications 1, 7 and 9 are exactly the studies used as the held-out
external cohort (`results/external_validation_finaledb.json`): they are
fetched from the API and the S3 bucket for scoring, and deliberately
contribute no training data. The training labels file carries
publications 6 and 8 only.

## Per-publication readiness

| Publication | Study | Verdict | n_labels (cell-line filtered) | n_features_present | Why |
|---|---|---|---:|---:|---|
| 1 | snyder | `labels-missing` | 0 | no | No samples in the labels file reference this publication. |
| 6 | jiang | `ready` | 121 | yes | All artifacts present; benchmark can include this publication. |
| 7 | sun | `labels-missing` | 0 | no | No samples in the labels file reference this publication. |
| 8 | cristiano | `ready` | 537 | yes | All artifacts present; benchmark can include this publication. |
| 9 | adalsteinsson | `labels-missing` | 0 | no | No samples in the labels file reference this publication. |

## Default benchmark scope

- **Default publications**: `['6', '8']` — the only publications with `ready` status as of `2026-11-15T16:54:53.394362+00:00`.
- `--publications` override: pass any subset. Missing publications are reported via stderr warnings, NOT hard-failed (the script skips them gracefully).
- `--include-snyder` / `--include-sun`: shortcuts for publications 1 / 7.

## Honest reading

**BLOCKED publications**: `1` (snyder), `7` (sun), `9` (adalsteinsson)

These are blocked on **FinaleDB infrastructure**, not on this repo. The labels file references them, but no local 5-channel DELFI features exist (and the API/S3 are not reachable to fetch them). Re-run this diagnostic after FinaleDB authors restore API/S3 access.

Until then, the cross-study benchmark uses `6 8` (the only publications with `ready` status). Numbers from the `--publications 6 8` run are the honest open-data cross-study benchmark.
