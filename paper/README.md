# DeepCatch Paper — Compilation Guide

This directory contains the LaTeX source for the DeepCatch manuscript.

## Files

| File | Description |
|------|-------------|
| `deepcatch_final.tex` | Primary manuscript (Bioinformatics target) |
| `references_final.bib` | Bibliography in BibTeX format |
| `supplementary.tex` | Supplementary materials (extended methods, additional figures) |
| `figures/` | Figure files (PNG, PDF, or TikZ sources)

## How to Compile

### Option 1: pdflatex + bibtex (traditional)

```bash
cd paper/
pdflatex deepcatch_final.tex
bibtex deepcatch_final
pdflatex deepcatch_final.tex
pdflatex deepcatch_final.tex
```

### Option 2: latexmk (recommended)

```bash
cd paper/
latexmk -pdf deepcatch_final.tex
```

### Option 3: Overleaf

Upload the entire `paper/` directory to Overleaf. Set the main document to `deepcatch_final.tex`.

### Option 4: GitHub CI

Push changes to `paper/deepcatch_final.tex` or `paper/supplementary.tex` — the CI workflow will automatically compile and upload the PDF as an artifact.

## Dependencies

- `natbib` for citation formatting
- `graphicx` for figures
- `booktabs` for tables
- `hyperref` for hyperlinks
- `lineno` for line numbers (submission requirement)

All are standard LaTeX packages available in TeX Live / MiKTeX.

## Submission Checklist

Before submitting:

- [ ] Remove `linenumbers` for camera-ready version
- [ ] Fill in all author names and affiliations
- [ ] Add corresponding author email
- [ ] Ensure all figures are in vector format (PDF preferred)
- [ ] Check reference formatting matches journal requirements
- [ ] Add data availability statement
- [ ] Add code availability statement (point to this GitHub repo)
- [ ] Add competing interests declaration
- [ ] Add author contributions section
- [ ] Add acknowledgments (funding, collaborators)

## Target Journals (Ranked)

1. **Nature Medicine** — IF 82.9 — Clinical translation focus
2. **Cancer Discovery** — IF 29.7 — High-impact cancer research
3. **Nature Communications** — IF 14.7 — Open access, good for methods
4. **Clinical Cancer Research** — IF 11.5 — Translational focus
5. **JCO Precision Oncology** — IF 5.3 — Precision medicine focus

## Current Status

- **Manuscript**: Draft complete, pending wet-lab validation results
- **Figures**: Descriptions written, figures to be generated
- **References**: 21+ papers reviewed and cited
- **Supplementary**: Structure in place, content to be expanded

*Last updated: 2026-04-28*
