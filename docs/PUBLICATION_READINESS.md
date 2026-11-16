# Publication Readiness (FinaleDB Open-Data Benchmark)
- Generated: `2026-11-16T14:40:29.983122+00:00`
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
Publications 1, 7, 9 are exactly the studies used as the held-out external cohort (`results/external_validation_finaledb.json`): they are fetched from the API and the S3 bucket for scoring, and deliberately contribute no training data. The training labels file carries publications 6, 8 only.

## Per-publication readiness

| Publication | Study | Verdict | Held-out external | n_labels (cell-line filtered) | n_features_present (training cache) | Why |
|---|---|---|---|---:|---:|---|
| 1 | snyder | `labels-missing` | yes | 0 | no | Not in the training labels file BY DESIGN: this publication is the held-out external cohort. Its features are fetched from the public API/S3 and scored by the external validation; it contributes no training data. |
| 6 | jiang | `ready` | no | 121 | yes | All artifacts present; benchmark can include this publication. |
| 7 | sun | `labels-missing` | yes | 0 | no | Not in the training labels file BY DESIGN: this publication is the held-out external cohort. Its features are fetched from the public API/S3 and scored by the external validation; it contributes no training data. |
| 8 | cristiano | `ready` | no | 537 | yes | All artifacts present; benchmark can include this publication. |
| 9 | adalsteinsson | `labels-missing` | yes | 0 | no | Not in the training labels file BY DESIGN: this publication is the held-out external cohort. Its features are fetched from the public API/S3 and scored by the external validation; it contributes no training data. |

## Default benchmark scope

- **Default publications**: `['6', '8']` — the only publications with `ready` status as of `2026-11-16T14:40:29.983122+00:00`.
- `--publications` override: pass any subset. Missing publications are reported via stderr warnings, NOT hard-failed (the script skips them gracefully).
- `--include-snyder` / `--include-sun`: shortcuts for publications 1 / 7.

## Honest reading

**Held-out external cohort** (not blocked): `1` (snyder), `7` (sun), `9` (adalsteinsson) — scored in `results/external_validation_finaledb.json`. They are absent from the training labels file *by design*: they are the external test set, fetched from the public API/S3 for scoring, contributing no training data. Their results are in `docs/EXTERNAL_VALIDATION.md`.

**No publication is blocked on FinaleDB infrastructure.**

The cross-study benchmark trains on the `ready` set (`6 8`). The held-out cohort is scored separately by `scripts/external_validation_finaledb.py`; its numbers are in `docs/EXTERNAL_VALIDATION.md`.
