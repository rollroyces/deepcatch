# Contributing to DeepCatch

Thank you for considering it. This repository is a single-author
open-source project, so the bar for a contribution is mostly **the
numbers have to stay honest** — everything else is ordinary practice.

## The one rule that is unusual here

**Every number in a document must be traceable to a committed artifact
under `results/`, and the artifact must be the one the code actually
produced.**

The repository enforces this with a drift guard
(`scripts/benchmark_numbers.py` plus `test/test_paper_numbers_match_artifacts.py`).
If you change a benchmark, regenerate the artifact, update the documents
that quote it, and let the guard tell you which ones you missed.

Two traps this project has already fallen into, so you do not have to:

- **A generator that writes a number it did not measure.** If you add a
  figure to a generated document, read it from the artifact rather than
  hard-coding it, and make the generator record the *effective*
  configuration it used — not the flags it was given.
- **A document that quotes a value from a different run than the one it
  names.** If a table's two sides come from different pipelines, that is
  a finding, not a comparison. Keep the rows on one artifact.

## Getting set up

```bash
git clone https://github.com/rollroyces/deepcatch
cd deepcatch
python3 -m venv .venv
./.venv/bin/pip install -e ".[dev,api]"   # [dev] adds pytest, [api] adds torch + the DL suite
```

The fragmentomics benchmarks read a **feature cache** extracted from
public FinaleDB data (per-sample `.npy` files) and a **labels TSV**.
Both live outside the repository; point the scripts at them:

```bash
export DEEPCATCH_FEATURES_DIR=/path/to/features
export DEEPCATCH_LABELS=/path/to/labels_multiclass.tsv
```

The sibling repository `cfdna-fragmentomics-pipeline` supplies the
5-channel loader and the harmonisation helper, so that the training and
external cohorts are assembled by exactly the same code path:

```bash
export DEEPCATCH_SISTER_PIPELINE=/path/to/cfdna-fragmentomics-pipeline
```

`paper/REPRODUCE.md` lists every step and every input the full
reproduction needs.

## Running the tests

```bash
# The CI's core job
python3 -m pytest src/fragmentomics src/preprocessing src/clinical \
    src/multimodal_fusion test/ -q --no-header -m "not slow"

# The CI's deep-learning job
python3 -m pytest src/foundation src/tissue_deconv src/priming \
    src/methylation_gnn test/test_foundation_smoke.py -q --no-header
```

The suite count is quoted in a dozen documents, and
`test/test_test_count_consistency.py` fails if a document and the real
collection disagree. **Adding or removing a test means updating that
count everywhere** — the guard names the documents.

## Contributing a dataset

Public, IRB-free data only: no patient-identifiable material belongs in
this repository, and the data-access templates in `data/` describe how
to request anything controlled.

To bring a cohort into the benchmark:

1. Fetch the features into the feature cache directory (the loader
   expects `{sample}.{channel}.npy` plus `{sample}.fsd.json`).
2. Write a labels TSV with the canonical columns
   (`sample`, `disease_class`, `label`, `study`, optionally `publication`).
   The label vocabulary is closed — `cancer`/`tumor`/… or
   `healthy`/`control`/… — and an unrecognised label is an error, never
   a control. See `docs/COLLABORATOR_DATA_INTERFACE.md`.
3. Run the cross-study benchmark with the cohort's publications and
   check the confound controls: the shuffled-label null must sit near
   chance and the true-confound control must collapse once harmonised.
4. Quote the new artifact's numbers in the documents that discuss it,
   and add it to the drift guard's quoted-number list.

## Style

- Python: follow the surrounding code; the project targets Python 3.10+.
- Keep new modules importable without the deep-learning extras if they
  live under `src/fragmentomics`, `src/preprocessing`, `src/clinical` or
  `src/multimodal_fusion` — the core CI job installs only the light
  dependencies.
- Write the honest result. A null is a result; a negative result that
  is documented beats a positive result that is not reproducible.

## Review and release

- Changes to a published number need the artifact regenerated in the
  same change, and the reproduction driver green.
- Tags and releases are cut by the maintainer after the reproduction
  check passes; a PR that changes a headline number should say which
  artifact it changed.
- Zenodo's GitHub integration is enabled for this repository, so a
  published release is archived automatically and mints a DOI. The DOI
  is declared in `CITATION.cff`, in `MODEL.md`'s BibTeX entry and in
  both manuscripts' availability sections, and
  `test/test_release_doi_is_consistent.py` checks that they agree
  (see `ZENODO_SETUP.md`).
- The version a document states is named by its **release tag**, never
  by a commit hash: a tag cannot drift, and a commit hash in a version
  pointer quietly stops being the commit the document describes.
- See `TEAM.md` for the governance model and the roles that are still
  open, and `docs/COLLABORATOR_DATA_INTERFACE.md` for what a clinical
  cohort needs to look like.
