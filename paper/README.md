# DeepCatch Paper — Compilation Guide

This directory contains the LaTeX source for the DeepCatch manuscript.

## Files

| File | Description |
|------|-------------|
| `main.tex` | Primary manuscript (Nature Medicine / Lancet Oncology / Cancer Discovery target) |
| `references.bib` | Bibliography in BibTeX format |
| `supplementary.tex` | Supplementary materials (extended methods, additional figures) |
| `abstract.txt` | Plain-text abstract for submission systems |
| `cover_letter.txt` | Template cover letter for journal submission |
| `figures/` | Figure files (PNG, PDF, or TikZ sources) |

## How to Compile

### Option 1: pdflatex + bibtex (traditional)

```bash
cd paper/
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

### Option 2: latexmk (recommended)

```bash
cd paper/
latexmk -pdf main.tex
```

### Option 3: Overleaf

Upload the entire `paper/` directory to Overleaf. Set the main document to `main.tex`.

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
