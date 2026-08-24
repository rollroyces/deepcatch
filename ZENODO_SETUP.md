# Zenodo DOI Setup

Zenodo gives every GitHub release a permanent DOI for citation. This
file explains the setup so a reviewer (or future me) can finish it.

## What's already in place

- **`CITATION.cff`** at repo root (citation metadata in CFF format,
  which GitHub natively understands for the "Cite this repository"
  button). Both `deepcatch` and `cfdna-fragmentomics-pipeline` have one.
- **README** has citation blocks.
- **MODEL.md** has a `@software` BibTeX entry ready for the DOI.
- **paper/PAPER.md** has a `@misc` BibTeX entry.

## What needs to happen (one-time, ~5 minutes)

### 1. Link your GitHub account to Zenodo

- Go to https://zenodo.org/account/settings/github/
- Find `rollroyces` in the list (under "GitHub repositories")
- For each repo (DeepCatch, cfdna-fragmentomics-pipeline), toggle it ON
- Zenodo will automatically archive future releases

### 2. Make a release tag

For `rollroyces/deepcatch` (v2.2.0 is the current version):

```bash
# From the local repo
cd /Users/hermes/deepcatch
git tag -a v2.2.0 b59ebc4 -m "DeepCatch v2.2.0: tumor-naive + fusion + DeLong + decision curve + AUC gate + MODEL.md"
git push origin v2.2.0
```

For `rollroyces/cfdna-fragmentomics-pipeline` (current tip is
`058c1e1`):

```bash
cd /Users/hermes/cfdna-fragmentomics-pipeline
git tag -a v0.2.0 058c1e1 -m "v0.2.0: cross-study AUC 0.9745, AUC gate, HEAD-size guard, 18 unit tests"
git push origin v0.2.0
```

### 3. Create the GitHub Release

- Go to https://github.com/rollroyces/deepcatch/releases/new
- Choose the tag, title `DeepCatch v2.2.0`, description with the
  headline numbers (paste from MODEL.md), attach
  `paper/biorxiv_submission_v2.2.0.pdf` as the release artifact
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
bioRxiv will assign a DOI. That DOI is just as good as Zenodo's for
the purposes of being cited, and it's more discoverable to the
bioinformatics community. The choice is yours; both work.

## Honest recommendation

- **Submit to bioRxiv first** (already prepared, see
  `paper/BIORXIV_SUBMISSION.md`).
- **Once the bioRxiv DOI is in hand, use it as the primary citation
  in the README** rather than a Zenodo DOI. bioRxiv is more visible
  in the bioinformatics community.
- **Tag the GitHub release with `v2.2.0-bioRxiv-2026-XX-XX`** so
  anyone reproducing the result uses the exact same code as the
  pre-print.
- **Zenodo is a backup**, not the primary. Skip it if you don't want
  the extra 5 minutes.

