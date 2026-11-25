# DeepCatch documentation

> 📘 **Looking for the project landing page?**
> → [https://rollroyces.github.io/deepcatch/](https://rollroyces.github.io/deepcatch/)
>
> That page is the public-facing entry point. The files below are the source
> documents — Markdown source for the paper, design notes, and reviewer pack.

This directory holds the canonical docs for the DeepCatch project:

## Public-facing docs site (deployed to GitHub Pages)

- **[site/index.html](site/index.html)** — landing page (rendered at
  [https://rollroyces.github.io/deepcatch/](https://rollroyces.github.io/deepcatch/))
- **[site/assets/style.css](site/assets/style.css)** — site stylesheet

## Manuscript + preprints (sources)

- **[paper/paper.tex →](../paper/paper.tex)** — LaTeX source (the current manuscript)
- **[paper/biorxiv_submission_v2.3.0.pdf →](../paper/biorxiv_submission_v2.3.0.pdf)** — pre-compiled PDF (the submission artefact; rebuilt from `paper/paper.tex` after the v2.3.0 results landed)
- **[PAPER.md →](../paper/PAPER.md)** — the earlier panel-only manuscript, kept for reference (superseded by `paper/paper.tex`; it predates the external-validation, leave-one-study-out and methylation-channel results)
- **[PREPRINT_ABSTRACT.md](PREPRINT_ABSTRACT.md)** — the paste-ready abstract for the bioRxiv form, generated from `paper/paper.tex`
- **[ARXIV_ABSTRACT.md](ARXIV_ABSTRACT.md)** — the paste-ready abstract for arXiv's metadata field (1,913 characters, against arXiv's 1,920), generated from the same manuscript; the paper itself keeps its full abstract
- **[PATH_TO_IMPACT.md](PATH_TO_IMPACT.md)** — companion-piece narrative

## Design + planning

- **[V3_DESIGN.md](V3_DESIGN.md)** — next-generation model design (GPU-accelerated, MPS)
- **[PORTFOLIO_ROADMAP.md](PORTFOLIO_ROADMAP.md)** — unified roadmap covering all 3 repos (supersedes the older NEXT_STEPS.md / PRODUCTION_ROADMAP.md)
- **[PRODUCTION_ROADMAP.md](PRODUCTION_ROADMAP.md)** — productionization plan (historical — superseded by PORTFOLIO_ROADMAP.md §2 Phase D)

## Long-form guides

- **[USER_GUIDE.md](USER_GUIDE.md)** — full user-facing guide (50 KB)
- **[../USAGE.md](../USAGE.md)** — repo-root CLI + API usage reference

---

**Author:** Yu Ching Lam · ORCID [0009-0008-9113-769X](https://orcid.org/0009-0008-9113-769X) ·
Independent Researcher, Hong Kong.

DeepCatch is research-stage software, not for clinical diagnosis. See
[../MODEL.md](../MODEL.md) for the model card.
