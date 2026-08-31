# DeepCatch Paper — Compilation Guide

This directory contains the LaTeX source for the DeepCatch manuscript
(DeepCatch v2.2.0, August 2026).

## Files

| File | Description | Submission status |
|------|-------------|-------------------|
| `paper.tex` | **Primary manuscript (submitted to bioRxiv)** | ✓ Submission PDF: `biorxiv_submission_v2.2.0.pdf` |
| `biorxiv_submission_v2.2.0.pdf` | Compiled submission PDF (201 KB, 5 pages) | ✓ Ready for bioRxiv |
| `BIORXIV_SUBMISSION.md` | Submission-form metadata + author info | ✓ Required fields present |
| `deepcatch_final.tex` | OLD DRAFT (different title, different scope, by "Yu Ching Lam") | ✗ Deprecated; do not use |
| `references_final.bib` | Bibliography (50 entries; 12 keys cited in `paper.tex`) | ✓ |
| `supplementary.tex` | Supplementary materials (16 sections) | Self-contained; not auto-compiled by CI |
| `figures/` | (empty — see notes below) | |

## How to Compile

### Primary (`paper.tex` → `biorxiv_submission_v2.2.0.pdf`)

The CI workflow `compile-paper.yml` automatically compiles `paper.tex` to
PDF on every push. To compile locally:

```bash
cd paper/
pdflatex paper.tex        # produces paper.pdf
```

The output `paper.pdf` is byte-identical to `biorxiv_submission_v2.2.0.pdf`
when compiled with pdfTeX 1.40+ and TeX Live 2023+. The CI artifact is the
authoritative version.

### Old draft (`deepcatch_final.tex`)

`deepcatch_final.tex` is a **deprecated** earlier draft (v1.0, April 2026).
It still compiles to `deepcatch_final.pdf` via the CI but is **not the
submitted manuscript** and has 4 unresolved `natbib` citations. It is
retained for reference only.

## Dependencies

All standard LaTeX packages — no exotic dependencies. The submission
PDF compiles cleanly under pdfTeX 1.40.25 (TeX Live 2023) without any
`??` warnings or missing references.

## Submission Checklist (bioRxiv)

Before submitting:

- [x] PDF compiles cleanly
- [x] All 12 `\cite{}` keys resolve to `references_final.bib`
- [x] Abstract ≤ 3000 characters
- [x] Submission metadata in `BIORXIV_SUBMISSION.md`
- [ ] **Register ORCID** (5 min, free) and update `BIORXIV_SUBMISSION.md`
- [ ] **Replace `[your email]` placeholder** with the corresponding-author email
- [ ] Optionally add 2 missing figures (see Notes below)

## Notes on figures

The submitted PDF embeds **no figures** at this time. PAPER.md and
`paper.tex` describe Figure 1 (panel detection performance across ctDNA
fractions) and Figure 2 (ultra-early assay sweep) by caption only. The
PNG assets for these are committed at:

- `results/real_tcga_performance.png` (66 KB) → Figure 1
- `results/benchmark_comparison.json` (data for the assay-sweep heatmap) → Figure 2

To add the figures, edit `paper.tex` to include:
```latex
\usepackage{graphicx}
...
\begin{figure}[h]
  \includegraphics[width=0.9\textwidth]{../results/real_tcga_performance.png}
  \caption{...}
  \label{fig:1}
\end{figure}
```

This was left for a future revision because bioRxiv accepts
figure-less preprints (rare but allowed) and the underlying PNGs
require a 30-second matplotlib regeneration step.

## Target Journal Recommendation

Based on the round-4 reviewer analysis (Q1-Q7 in
`/Users/hermes/ROUND4_STATISTICAL_AUDIT.md`), the realistic target
journals for this work, in order:

1. **Bioinformatics** (Oxford) — methods-focused, accepts
   single-author, less concerned with clinical validation
2. **NAR Genomics & Bioinformatics** — same
3. **PLOS Computational Biology** — accepts negative results,
   methods-focused
4. **Briefings in Bioinformatics** — review-with-methods format

For Nature Medicine / Cancer Discovery (clinical impact journals),
the work would need: per-cancer-type AUC table, held-out clinical
cohort, and re-running the fusion with real (not synthetic) mutation
channel. Estimated 6-12 months additional work.

## Current Status

- **Manuscript (v2.2.0)**: ✓ Compiled to PDF, ready for bioRxiv
  submission (after ORCID + email fill-in).
- **Companion repo** (cfdna-fragmentomics-pipeline): at HEAD; both
  repos have synchronized `RESULTS.md`, `AUDIT_REPORT_2.md`.
- **Audit trail**: 4 rounds of multi-reviewer audits completed.
  See `AUDIT_REPORT_2.md` and the round-4 reports in `~/`:
  - `JOURNAL_REVIEW_REJECTION_ANALYSIS.md` (round 4 journal reviewer)
  - `ROUND4_STATISTICAL_AUDIT.md` (round 4 statistical reviewer)
  - `cfdna-fragmentomics-pipeline/AUDIT_REPORT_4.md` (round 4 engineering)

*Last updated: 2026-08-28*
