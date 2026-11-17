# Where to publish this, and what is still missing

**Purpose.** Decide the submission sequence and track exactly what each
venue needs. The short version: **preprint first, then the methods
manuscript**, and the two candidates for the latter are JOSS (if the
software is the contribution) and PLOS Computational Biology /
Bioinformatics Advances (if the measurement is).

Everything in the "ready" column below is verifiable in this repository;
the "missing" column is the honest gap.

---

## 1. The submission sequence

```
        bioRxiv preprint  ──►  JOSS  (software paper, DOI via Zenodo)
              │                    └─ OR ─►  PLOS Computational Biology
              │                              / Bioinformatics Advances
              └─►  (data + DOI) available for whichever journal follows
```

The preprint is not optional: bioRxiv assigns the DOI and the priority,
it needs no peer review to accept, and it is what makes the *journal*
submission credible for an author without an institutional affiliation.
Everything in the manuscript is already inside it, so the preprint costs
one upload.

---

## 2. Venue comparison

| | **bioRxiv** (preprint) | **JOSS** | **PLOS Computational Biology** | **Bioinformatics Advances** | **GigaScience** |
|---|---|---|---|---|---|
| Reviews | none (moderated) | software + docs | science | science | science + data |
| Fits this work? | yes — priority + DOI | **yes** — the contribution *is* the reproducible pipeline | yes — the honest external validation and the nulls | yes — methods emphasis | partly — wants more biological interpretation |
| Length | no limit | ~2,500 words | no limit | no limit | no limit |
| Data/code required | code URL | **public repo, installable, versioned** | data availability statement | data availability | data + workflow |
| Needs clinical cohort | no | no | no (but reviewers will ask) | no | likely |
| Time to decision | days | 1–3 months | 2–4 months | 1–3 months | 2–4 months |
| DOI | yes (auto) | **yes (Zenodo, needs the GitHub integration)** | yes | yes | yes |

**Recommendation.** Submit the preprint to bioRxiv now, then target
**JOSS** with the software framing (the reproducibility guarantees,
the artifact-as-tests design and the honest-null discipline are the
reviewable contribution), and keep **PLOS Computational Biology** as the
parallel science submission if a second journal is wanted.

Rationale for JOSS specifically: it is the only one of these venues that
*reviews the software*, which is exactly the axis the work is strongest
on. It also requires exactly the things this repository now has — a
versioned public repo, an installable package, a test suite, and a
DOI — and it does not ask for prospective clinical data that no
independent researcher can produce.

---

## 3. Readiness, venue by venue

### bioRxiv — ready except for the upload

| Requirement | State |
|---|---|
| Compiled PDF | ✅ `paper/biorxiv_submission_v2.3.0.pdf` (8 pages, 2 figures embedded) |
| Abstract within limits | ✅ |
| Submission metadata | ✅ `paper/BIORXIV_SUBMISSION.md` (gitignored: it carries personal details) |
| ORCID | ✅ `0009-0008-9113-769X`, also in `CITATION.cff` |
| Code URL in the manuscript | ✅ |
| **Affiliation consistent with the ORCID record** | ⚠️ **decide before uploading** — the manuscript's author block and `CITATION.cff` declare "Independent Researcher (no institutional affiliation)", but the public ORCID record for `0009-0008-9113-769X` (checked 2026-11-16) carries **one public employment: University of Michigan** (Ann Arbor, MI, US; ROR `00jmfr291`) with **no start date, no end date, no department and no role title** — so ORCID renders it as a *current* employment — and no education, works or other affiliations of any kind. bioRxiv asks for an affiliation, and a reviewer who clicks the ORCID sees the record. Two coherent resolutions, both a few minutes' work: **(a)** the ORCID entry is stale or a placeholder → end-date or delete it and keep "Independent Researcher"; **(b)** it is current → declare the affiliation in `paper/paper.tex` (author block, line 13) and `CITATION.cff`, which makes the paper *stronger*, not weaker, for a reviewer who distrusts an unaffiliated author. Which of the two is true is the author's fact, not the repository's — but the two must not contradict each other on the day of upload. |
| **Upload** | ⬜ author action — needs the bioRxiv account |
| **Corresponding-author email** | ⬜ author action — deliberately not in the repo |

### JOSS — one structural gap

| Requirement | State |
|---|---|
| Public, versioned repository | ✅ MIT-licensed, tagged releases |
| Installable package | ✅ `pyproject.toml` |
| Test suite | ✅ 835 tests (`pytest src test -q`), split across two CI jobs |
| Documentation | ✅ `README.md`, `MODEL.md`, `USAGE.md`, `docs/` |
| Citation metadata | ✅ `CITATION.cff` |
| Contribution guidelines | ✅ `CONTRIBUTING.md` (setup, the drift-guard rules, how to bring a cohort in) |
| Author list | ✅ `AUTHORS` |
| Code of conduct | ✅ `CODE_OF_CONDUCT.md` |
| **DOI** | ⬜ requires Zenodo's GitHub integration + the `v2.3.0` tag |
| **Paper file discoverable from the repo root** | ✅ `paper.md` in the repo root, in JOSS format (summary, statement of need, features, quality/reproducibility, documentation, code availability, limitations, acknowledgements) |
| **Short description** (the JOSS submission form's one-paragraph field) | ✅ the Summary section of `paper.md` is the one-paragraph description |

### PLOS Computational Biology / Bioinformatics Advances

| Requirement | State |
|---|---|
| Methods manuscript | ✅ `paper/METHODS_PAPER.md` (all sections filled, numbers regenerated) |
| Data availability statement | ✅ §5 |
| Code availability | ✅ MIT, public, versioned |
| **Acknowledgements section** | ⬜ **neither manuscript has one.** The data providers are named in the data-availability statements and `TEAM.md` credits the reviewers, but PLOS and Bioinformatics Advances both expect an acknowledgements section — whom to name, and how to disclose that the multi-persona audits were AI-assisted, is the author's call |
| Figure files | ✅ `docs/figures/fig1..fig4*.png`, rendered from the committed artifacts |
| **Statistical reviewer concerns** | ⚠️ expect a request for the per-cancer denominators and the single-cohort HCC_J reading — both are already stated in §3.3 and §4.2 |
| **DOI** | ⬜ same Zenodo dependency |

---

## 4. The one thing that gates every journal

**Zenodo's GitHub integration.** Every venue above wants a DOI for the
software, and Zenodo only mints one automatically if the repository is
linked to the account. The steps are in `ZENODO_SETUP.md`:

1. Enable the repository at <https://zenodo.org/account/settings/github/>
2. Tag the release: `git tag -a v2.3.0 -m "v2.3.0" && git push origin v2.3.0`
3. Zenodo archives the release and issues the DOI; paste it into
   `CITATION.cff`, `MODEL.md` and `paper/paper.tex`.

Without it, the preprint can still go out (bioRxiv issues its own DOI),
but the journal submissions would cite a commit hash instead of a
permanent identifier.

---

## 5. What a reviewer will press on, and the answer

| Likely objection | Where the answer is |
|---|---|
| "The pooled AUC is measured inside the cohort that trained it." | §3.2.1 and the external validation in §3.2.2: 0.9444 on a held-out whole-genome publication, 0.708 on a second, 0.6805 in the single-type direction, with the reason spelled out. |
| "Is the methylation channel real?" | §3.2.4: it is a **null**, with the mechanism (size selection at 300 bp, 1.09 observations per CpG) and a control (depth ablation on a genuine MspI library). |
| "Are the numbers reproducible?" | §2.5: the PCA is seeded, and `test/test_benchmark_reproducibility.py` asserts that two runs agree exactly. |
| "Do the artifacts match the text?" | `scripts/benchmark_numbers.py` prints every quoted figure straight from the committed JSON. |
| "No clinical cohort." | §4.2 and §4.3: this is an open-data methods benchmark, not a clinical validation, and it says so. |
