# DeepCatch Paper — Compilation and Submission Guide

This directory holds the DeepCatch manuscripts and the artefacts needed to
submit them. **Status as of v2.3.0 (2026-11-16): the preprint is compiled
and current; it has NOT been uploaded** — the remaining steps need the
author's arXiv account and the Zenodo DOI.

## Which manuscript is which

| File | What it is | State |
|------|------------|-------|
| `paper.tex` | **Primary manuscript** (preprint). Panel-based MRD detection, the tumor-naive fragmentomics channel, the held-out external validation, the leave-one-study-out transfer measurement and the methylation null. | ✓ compiles; `biorxiv_submission_v2.3.0.pdf` is the submission artefact |
| `biorxiv_submission_v2.3.0.pdf` | Compiled submission PDF (9 pages, 3 figures) | ✓ current — regenerated from `paper.tex` after the v2.3.0 results landed |
| `METHODS_PAPER.md` | The long-form methods manuscript (the PLOS Computational Biology / Bioinformatics Advances article; JOSS gets the root `paper.md`). Same evidence, more detail, per-cancer tables and a reproducibility paragraph. | ✓ numbers regenerated from the seeded benchmark |
| `METHODS_PAPER.tex` | The methods manuscript as LaTeX — **generated** from `METHODS_PAPER.md` by `scripts/markdown_to_latex_paper.py`, not transcribed. Do not hand-edit it: the next generation erases the edit. | ✓ compiles under pdfTeX and tectonic; CI compiles it and uploads the PDF. **Not yet a PLOS submission file**: it uses the `article` class, and PLOS Comp Bio requires its own template (see `docs/JOURNAL_SUBMISSION.md`) |
| `BIORXIV_SUBMISSION.md` | Submission-form metadata: what the preprint claims, what it must not claim, and what still needs the author | ✓ (gitignored — it holds personal details) |
| `METHODS_PAPER_STATUS.md` | Section-by-section status meta-doc for the methods manuscript | ✓ updated with the v2.3.0 numbers |
| `PAPER.md` | The earlier panel-only manuscript, kept for reference | superseded by `paper.tex` |
| `deepcatch_final.tex` | Deprecated v1.0 draft (April 2026), different title and scope | ✗ do not use |
| `supplementary.tex` | The **deprecated v1.0 draft's** supplementary (the multi-modal product framing: CTC/miRNA detectors, contrastive variant calling, MAML, longitudinal CET monitoring). It describes figures S1–S10 and tables S1–S5 that the current manuscript does not contain, and it ships no image files. | ✗ do not upload with `paper.tex` |
| `references_final.bib` | BibTeX file from the Bioinformatics-targeted version of the manuscript. The current manuscript carries its references inline, so this file is not its bibliography and nothing in the repository builds from it | reference only |

## What the preprint now reports

Everything below is backed by a committed artifact under `results/`:

| Claim | Number | Artifact |
|---|---|---|
| Panel LLR at 0.1% ctDNA | AUC 0.921 (Sens@95% 0.770) | `results/real_tcga_validation.json` |
| Tumor-naive fragmentomics, pooled, GC-corrected | 0.9684 ± 0.0027 | `results/cross_study_finallydb_gc_corrected.json` |
| Held-out external cohorts (never trained on) | 0.9444 [0.8731–1.0000] (pub 7), 0.708 (pub 1), pooled 0.7385 | `results/external_validation_finaledb.json` |
| Leave-one-study-out transfer | 0.6805 / 0.9361 | `results/cross_study_loso.json` |
| Late fusion vs the neural encoder | 0.9825 ± 0.0094 vs 0.554 | `results/fusion_regularized_smoke.json` |
| Methylation channel | **null** (0.9444 = the coverage baseline) | `results/methylation_cohort_channel.json` |

## How to compile

CI compiles `paper.tex` on every push (`.github/workflows/compile-paper.yml`).
Locally, with tectonic:

```bash
cd paper/
tectonic -X compile paper.tex     # -> paper.pdf
```

or with pdfTeX:

```bash
cd paper/
pdflatex paper.tex && pdflatex paper.tex
```

The CI artefact is the authoritative PDF. The two figure files the
manuscript embeds are regenerated from the committed artifacts:

```bash
python scripts/plot_figures.py    # -> docs/figures/fig1..fig4*.png
```

### The methods manuscript's LaTeX

`paper/METHODS_PAPER.tex` is generated from `paper/METHODS_PAPER.md`:

```bash
python scripts/markdown_to_latex_paper.py            # writes the .tex
python scripts/markdown_to_latex_paper.py --check    # checks it instead
```

`--check` does two things: it compares the committed `.tex` with a fresh
generation, and it compares the two documents' texts word for word. The
second check is the one that catches a converter bug — a dropped sentence
or a changed number — which the first cannot see, because a converter
that loses a sentence produces the same `.tex` every time.

Both checks run in CI, and the LaTeX is compiled there too:

```bash
cd paper/
pdflatex METHODS_PAPER.tex && pdflatex METHODS_PAPER.tex
```

The markdown is the source; the LaTeX is the submission artefact. Edit the
markdown and regenerate. The compiled PDF is a CI artefact rather than a
committed file — `paper/*.pdf` is gitignored — so this manuscript's page
count is not stated in this file: a number about a PDF that the
repository does not keep is a number that cannot be checked, and the
submission PDF's page count is checked against the PDF that *is* kept
(`test/test_submission_description_matches_the_pdf.py`).

## Submission checklist (the preprint)

- [x] PDF compiles cleanly under tectonic and pdfTeX
- [x] All three figures embedded and rendered from the committed artifacts
- [x] All `\cite{}` keys in `paper.tex` resolve against the manuscript's own inline bibliography
- [x] Abstract within both limits: 2,940 characters for bioRxiv's 3,000, and a 1,915-character abridged text for arXiv's 1,920-character metadata field (`docs/ARXIV_ABSTRACT.md`, generated from the manuscript)
- [x] Submission metadata in `BIORXIV_SUBMISSION.md` §3a (the arXiv form)
- [x] ORCID registered (`0009-0008-9113-769X`) and present in `CITATION.cff` — and in all four author blocks (`test/test_author_affiliation_is_consistent.py`)
- [ ] **Upload the source archive to arXiv** — arXiv takes TeX source, not a PDF: the `deepcatch-arxiv-source` artifact of the latest "Compile Paper PDF" run (it holds `deepcatch_v2_3_0.tar.gz`, one `.tex` plus the three flattened figures, and it compiles standalone to the same 9 pages as the committed PDF). Everything is ready; the last steps are the account (username, password, captcha, then the email verification), the "Link my arXiv account with ORCID" click on arXiv's account page, and the "Submit Article" click — arXiv asks the author to submit, so those are only the author's to make
- [x] **AI-use disclosure inserted** into `paper.tex` § Acknowledgements and
      `METHODS_PAPER.md` §6 (Block B in `BIORXIV_SUBMISSION.md` §7) -- it names
      the providers and the author's confirmation of review, and
      `test/test_author_affiliation_is_consistent.py` fails if any of the four
      manuscripts loses it
- [ ] **Add the corresponding-author email** to the submission form (kept out of the repo on purpose)
- [x] **Zenodo's GitHub integration enabled and `v2.3.0` tagged** — the release is archived as **10.5281/zenodo.24453044** (see `ZENODO_SETUP.md`)

## Figures

| Figure | Source | Regenerated by |
|---|---|---|
| Figure 1 — panel ROC across ctDNA fractions | `results/real_tcga_performance.png` | `real_tcga_validation.py` |
| Figure 2 — pooled cross-study ROC + sens@spec | `docs/figures/fig1_pooled_roc.png`, `docs/figures/fig4_sens_operating_points.png` | `scripts/plot_figures.py` |

## Target journals

The venue decision and its per-venue readiness checklist live in
[`docs/JOURNAL_SUBMISSION.md`](../docs/JOURNAL_SUBMISSION.md). In short:
the **preprint first** (arXiv, to establish priority), then the
**methods manuscript** to a journal that values software-engineering
quality — JOSS for the software, PLOS Computational Biology or
Bioinformatics Advances for the science.

## A note on the numbers quoted from the companion repository

The abstract's tumor-naive and fusion numbers are the ones this
repository produces: **0.9740 ± 0.0013** for the channel and **0.9884 ±
0.0010** for the naive average with a calibrated synthetic mutation
channel, both from `results/fusion_ablation_isotonic.json`
(`src/fragmentomics/fusion_ablation.py`, seeded PCA, 10 seeds).

The companion repository's `scripts/honest_benchmark.py` fits its PCA
without `random_state`, so its own headline moved on every run. Measured
on 2026-11-17, on an isolated worktree at that repository's HEAD (its
uncommitted work untouched): the unseeded run gave **0.9744 ± 0.0032**
and raised an lbfgs "failed to converge" warning, while the seeded run
gives **0.9770 ± 0.0013** with no warning — the same number this
repository reports for its legacy pooled figure, which is the check that
the two compute the same thing once both are seeded. Seven more PCA
calls in that repository were seeded on the same visit, and its
sens@spec and calibration figures now carry the same pooled OOF AUC.
The changes are prepared on that repository's local branch
`fix/seed-pca` (not pushed).
