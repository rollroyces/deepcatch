# Zenodo DOI Setup

Zenodo gives every GitHub release a permanent DOI for citation. This
file explains the setup so a reviewer (or future me) can finish it.

**Current release to tag: `v2.3.0`** (the tag commands below are
version-specific — an older copy of this file pointed at `v2.2.0`,
which would have put the DOI on the wrong version).

## What's already in place

- **`CITATION.cff`** at repo root (citation metadata in CFF format,
  which GitHub natively understands for the "Cite this repository"
  button). Both `deepcatch` and `cfdna-fragmentomics-pipeline` have one.
- **README** has citation blocks.
- **MODEL.md** has a `@software` BibTeX entry ready for the DOI.
- **paper.md** and `paper/PAPER.md` have `@misc` BibTeX entries.
- The submission artefact is `paper/biorxiv_submission_v2.3.0.pdf`
  (compiled from the current `paper/paper.tex`).

## What needs to happen (one-time, ~5 minutes)

### 1. Link your GitHub account to Zenodo

- Go to https://zenodo.org/account/settings/github/
- Find `rollroyces` in the list (under "GitHub repositories")
- For each repo (DeepCatch, cfdna-fragmentomics-pipeline), toggle it ON
- Zenodo will automatically archive future releases

### 2. Make a release tag

For `rollroyces/deepcatch` (**v2.3.0** is the current version — tag
the tip of `main`, not an older commit):

```bash
# From the local repo
cd /path/to/deepcatch
git checkout main && git pull
git tag -a v2.3.0 -m "DeepCatch v2.3.0: held-out external validation (0.9444), leave-one-study-out transfer, regularized late fusion (0.9825), methylation channel measured as a null, reproducible benchmark"
git push origin v2.3.0
```

For `rollroyces/cfdna-fragmentomics-pipeline` (check the tip before
tagging — that repository has its own release history):

```bash
cd /path/to/cfdna-fragmentomics-pipeline
git tag -a v0.2.0 -m "v0.2.0: cross-study AUC 0.9770, AUC gate, HEAD-size guard, 18 unit tests"
git push origin v0.2.0
```

> **Seed that repository's PCA before you tag it.** Its
> `scripts/honest_benchmark.py` fits `PCA(...)` without `random_state`,
> so its headline moves on every run: two runs of the same command on the
> same data gave **0.975618** and **0.973906**, with a reported ±std that
> also changed (0.000923 vs 0.002371). The fix is one line — give the PCA
> the CV seed that is already in scope:
>
> ```python
> pca = PCA(n_components=min(pca_n, max_pca),
>           random_state=sd).fit(Xtr)
> ```
>
> Measured with that one line, the same command returns **0.977018 ±
> 0.001316** on every run — identical to DeepCatch's own legacy pooled
> number (0.977018 ± 0.001316), which confirms the two repositories run
> the same computation once both are seeded. The sibling's documents
> (`BENCHMARK.md`, `docs/MULTICANCER_RESULTS.md`) still quote the
> unseeded 0.9745 and need updating at the same time; a DOI minted on a
> tag whose headline number is known to drift is exactly the failure
> this file warns about.

### 3. Create the GitHub Release

- Go to https://github.com/rollroyces/deepcatch/releases/new
- Choose the tag `v2.3.0`, title `DeepCatch v2.3.0`, description with
  the headline numbers (paste from MODEL.md), attach
  `paper/biorxiv_submission_v2.3.0.pdf` as the release artifact
- Click "Publish release"

Repeat for the pipeline repo.

### 4. Zenodo archives the release

- Within ~5 minutes of publishing, Zenodo will create a DOI for the
  release. It will appear at:
  https://zenodo.org/record/<id>
- Zenodo will also POST a comment on the GitHub release with the DOI.

### 5. Wire the DOI into the repo

Replace `10.5281/zenodo.to-be-assigned-on-first-release` in
`CITATION.cff` and any other placeholder DOIs with the real ones
Zenodo assigned. Commit. Done.

## What the DOI gives you

- **Permanence**: Zenodo DOIs don't break if GitHub does.
- **Citation tracking**: Zenodo tracks how many times the DOI is
  resolved.
- **Versioning**: each GitHub release tag gets its own DOI, so
  reviewers can cite the exact version they tested.

## What the DOI does NOT give you

- **Peer review.** A Zenodo DOI is not the same as a journal
  acceptance. For peer review, follow up with bioRxiv → journal
  submission.
- **Indexing in PubMed.** Zenodo DOIs are not PubMed-indexed. bioRxiv
  DOIs are also not PubMed-indexed. Only journal-acceptance DOIs
  end up in PubMed.

## Alternative: just use bioRxiv's DOI

Once the bioRxiv submission is accepted (~24-48 hour screening),
bioRxiv will assign a DOI. That DOI is just as good as Zenodo's for the
purposes of being cited, and it's more discoverable to the
bioinformatics community. The choice is yours; both work.

## Honest recommendation

- **Submit to bioRxiv first** — the PDF is
  `paper/biorxiv_submission_v2.3.0.pdf` and the submission metadata is
  in `paper/BIORXIV_SUBMISSION.md` (see also
  `docs/JOURNAL_SUBMISSION.md` for the venue map).
- **Once the bioRxiv DOI is in hand, use it as the primary citation
  in the README** rather than a Zenodo DOI. bioRxiv is more visible
  in the bioinformatics community.
- **Tag the GitHub release with `v2.3.0`** (and, if you want the
  preprint code pinned separately, `v2.3.0-bioRxiv-2026-XX-XX`) so
  anyone reproducing the result uses the exact same code as the
  pre-print.
- **Zenodo is a backup**, not the primary. Skip it if you don't want
  the extra 5 minutes — but if you do tag, tag **v2.3.0**, because a
  DOI minted on an older tag would misrepresent which code produced the
  reported numbers.
