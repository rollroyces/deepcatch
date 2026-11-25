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
| DOI | yes (auto) | **yes (Zenodo: 10.5281/zenodo.24453044)** | yes | yes | yes |

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
| Compiled PDF | ✅ `paper/biorxiv_submission_v2.3.0.pdf` (9 pages, 3 figures embedded). **Upload this file, not an older one:** the submission that bioRxiv rejected on 2026-11-25 carried the v2.0-era title "DeepCatch: Performance-Weighted Multi-Modal Fusion with Cumulative Evidence Tracking for Pan-Cancer Detection from Cell-Free DNA", so it was not this manuscript. The current title is "An Open-Source Benchmark for Panel-Based Ultra-Sensitive Detection of Circulating Tumor DNA from Real Tumor Mutations" |
| Abstract within limits | ✅ |
| Submission metadata | ✅ `paper/BIORXIV_SUBMISSION.md` (gitignored: it carries personal details) |
| ORCID | ✅ `0009-0008-9113-769X`, also in `CITATION.cff` |
| Code URL in the manuscript | ✅ |
| **Affiliation consistent with the ORCID record** | ⚠️ **decide before uploading** — the manuscript's author block and `CITATION.cff` declare "Independent Researcher (no institutional affiliation)", but the public ORCID record for `0009-0008-9113-769X` (re-checked 2026-11-24 directly against the ORCID public API, `https://pub.orcid.org/v3.0/0009-0008-9113-769X/record`) carries **one public employment: University of Michigan** (Ann Arbor, MI, US; ROR `00jmfr291`) with **no start date, no end date, no department and no role title** — so ORCID renders it as a *current* employment — and **zero** education groups, works or other affiliations. The record is the author's own published statement of affiliation, and a reviewer who clicks the ORCID sees exactly this. bioRxiv asks for an affiliation, and a reviewer who clicks the ORCID sees the record. Two coherent resolutions, both a few minutes' work, and **neither of them is a default — which of the two is true is the author's fact, not the repository's**: **(a) the employment is current** → declare "University of Michigan, Ann Arbor, Michigan, USA" in `paper/paper.tex` (author block, line 13), in `CITATION.cff` and in every other surface that declares an affiliation on the same day (`docs/SCHOLAR_VISIBILITY.md` §4 lists them all), which makes the paper *stronger*, not weaker, for a reviewer who distrusts an unaffiliated author; **(b) the record is stale** → end-date or delete the ORCID employment first and keep "Independent Researcher". Either way the two must not contradict each other on the day of upload. `paper/BIORXIV_SUBMISSION.md` §7 records the same two options under the same letters, and the local checklist there is what the form is filled from. |
| **Affiliation is a registration requirement, not a formality** | ⚠️ **bioRxiv's own submission guide says so**: "Authors wishing to deposit manuscripts must first register on the site and supply an affiliation so there is an entity providing oversight we can contact in the event concerns are raised about research misconduct" (https://www.biorxiv.org/submit-a-manuscript, checked 2026-11-17). openRxiv's 2026-10-20 statement repeats it: "All bioRxiv and medRxiv authors must therefore be verifiable human beings and affiliated with an organization to which we can turn if there are allegations of fraud or misconduct." "Independent Researcher" is not an entity, so the row above is a real gate on the upload, not a cosmetic inconsistency. **How hard it is in practice — answered, and the earlier answer was wrong:** counted from bioRxiv's own API (`http://api.biorxiv.org/details/biorxiv/2026-01-01/2026-01-31/0/json`, all 5,480 preprints posted in January 2026), **115 of them, 2.1%, are single-authored**, and the earlier version of this row read that as evidence that the rule is not enforced. That inference was wrong, because it measured the wrong population: a single-authored preprint can perfectly well carry an institutional affiliation (one PhD student is a single author), so "single-authored" and "affiliation-less" are not the same thing. **A submission with no institutional affiliation was rejected on 2026-11-25**, and bioRxiv's own words were: *"your manuscript is inappropriate for bioRxiv because bioRxiv requires authors to add an institutional affiliation. As a channel for rapid communication of science that has not undergone formal peer review, we require submissions to be associated with an organization that can provide oversight of research activities."* (manuscript ID BIORXIV/2026/725086; nothing was posted, no DOI was minted, and the record does not appear in bioRxiv's API). The rule is therefore a hard gate, and the open question is no longer *whether* to declare an affiliation but **which one**, and whether the one declared is true |
| **Disclosure of AI use** | ⚠️ **expected, and cheap to satisfy.** bioRxiv's FAQ: "Authors are responsible for the content they submit, including material produced by generative AI and related technologies. As with other computational methodologies, use of these tools should be detailed in manuscripts, and the output must not plagiarize other work or breach copyright." openRxiv's 2026-10-20 statement adds that an author who uses AI "must shape, review, stand behind and have intellectual input into the interpretation of the work they submit", and that a submission whose AI-generated material shows no critical review "may receive further review". The multi-persona audits this project ran are documented in the repository itself (`docs/`, `AUDIT_2_FINDINGS.md`, `CHANGELOG.md`), so a one-sentence disclosure in the acknowledgements is the consistent choice. Ready-to-paste wording: the local (gitignored) checklist `paper/BIORXIV_SUBMISSION.md` §7 |
| **Upload** | ⬜ author action — but it is the *account* that is missing, not the upload. The steps, in order, verified against the live site on 2026-11-24: **1)** open `submit.biorxiv.org/cgi/registration` and enter the corresponding author's email — bioRxiv checks whether an account exists for that address, so this one field is the whole of step 1; **2)** if no account exists, create one (name, password, affiliation); if one exists, sign in with it — **Bench>Press shares ONE registration across ~70 journals** (medRxiv, the CSH journals, *Genetics*, *Genome Research*, *The Plant Cell* and ~65 more), so an address already registered with any of them is already a bioRxiv account: do not repeat the registration, just sign in with the same email and password and complete the bioRxiv profile page; **3)** in the Personal Info area, associate the ORCID iD `0009-0008-9113-769X` with the account — **until that is done, ORCID sign-in cannot work at all**: bioRxiv's callback answers *"No account is associated with orcid, please login and associate orcid via the edit personal information page"*, which is what the `orcidcode=0003` in the returned URL records, so the ORCID side can succeed while the bioRxiv side still refuses the author; **4)** the author area opens and the manuscript can be submitted. **The root cause is in the ORCID record itself**: it exposes **no email address** (ORCID keeps email addresses private by default — read from the public API on 2026-11-24), so "Create an account with ORCiD" has no address to register with, and the iD has to be attached to an account made with an email address instead. The affiliation given in step 2 is the affiliation the preprint must declare, so it is the same decision as the row above and the repository has to be changed to match it on the same day |
| **arXiv fallback, source bundle prepared** | ✅ **built and verified.** arXiv takes a LaTeX source archive rather than a PDF, and this manuscript reaches outside the repository for three figures (`results/real_tcga_performance.png`, `docs/figures/fig1_pooled_roc.png`, `docs/figures/fig4_sens_operating_points.png`), so a bundle made by copying `paper/paper.tex` alone would compile with three missing figures. `python scripts/make_arxiv_bundle.py --verify` flattens them into `dist/deepcatch_v2_3_0.tar.gz`, compiles it standalone and compares the result with the committed PDF (9 pages, the same). It is the fallback that needs no institutional affiliation — which is exactly the field bioRxiv's registration requires. The paper workflow now uploads the bundle as the `deepcatch-arxiv-source` artifact, so the archive that is actually submitted is the one built from the committed manuscript rather than a local copy that can be stale |
| **Direct transfer to the journal (B2J)** | ✅ **available for the target venues, and it saves re-typing the paper.** bioRxiv's submission guide: "bioRxiv can save authors time in submitting papers to journals or peer review services by transmitting their manuscript files and metadata directly from the bioRxiv site." The partner list includes **PLOS Computational Biology, PLOS Genetics, PLOS Biology, PLOS Medicine, PLOS One, Bioinformatics, Bioinformatics Advances, GigaScience, eLife, PNAS, Nature and Science**. So the preprint is not just a priority record: it is also the submission form, already filled in, for every venue on the list — which makes "preprint first, then PLOS Comp Bio / Bioinformatics Advances / GigaScience" a one-click handoff rather than a re-entry |
| **A posted preprint is permanent** | ℹ️ "After an article is deposited in bioRxiv, it is citable and therefore cannot be removed." Revisions keep the same DOI and must carry their supplemental material again, so the PDF uploaded is the version that stays — which is the reason the uploaded PDF is the committed, CI-built `paper/biorxiv_submission_v2.3.0.pdf` and not a locally rebuilt one |
| **Corresponding-author email** | ⬜ author action — deliberately not in the repo |
| **Discoverability (ORCID, Google Scholar)** | ✅ the metadata is consistent and guarded; the three account-side steps are in [`docs/SCHOLAR_VISIBILITY.md`](SCHOLAR_VISIBILITY.md) (post the preprint, add the two works to the ORCID record, claim one Scholar profile) |

### JOSS — one structural gap

| Requirement | State |
|---|---|
| Public, versioned repository | ✅ MIT-licensed, tagged releases |
| Installable package | ✅ `pyproject.toml` |
| Test suite | ✅ 1037 tests (`pytest src test -q`), split across two CI jobs |
| Documentation | ✅ `README.md`, `MODEL.md`, `USAGE.md`, `docs/` |
| Citation metadata | ✅ `CITATION.cff` |
| Contribution guidelines | ✅ `CONTRIBUTING.md` (setup, the drift-guard rules, how to bring a cohort in) |
| Author list | ✅ `AUTHORS` |
| Code of conduct | ✅ `CODE_OF_CONDUCT.md` |
| **DOI** | ✅ **10.5281/zenodo.24453044** (v2.3.0, archived and published; declared in `CITATION.cff`, `MODEL.md` and both manuscripts) |
| **Paper file discoverable from the repo root** | ✅ `paper.md` in the repo root, in the **current** JOSS format: summary, statement of need, **state of the field**, **software design**, features, quality/reproducibility, documentation, code availability, **research impact statement**, limitations, **AI usage disclosure**, acknowledgements. The four bold sections are new: JOSS changed its submission format in late 2025 and its own review bot now checks for each of them, so a paper in the old format is rejected in pre-review. `test/test_author_affiliation_is_consistent.py` requires all twelve headings, and checks that the AI disclosure is a real disclosure (it must name a model provider and confirm the human review), because JOSS's AI usage policy treats an incomplete one as a research-integrity matter |
| **Short description** (the JOSS submission form's one-paragraph field) | ✅ the Summary section of `paper.md` is the one-paragraph description |
| **Submission** | ⬜ author action — the JOSS form, and it is **not** blocked on the affiliation decision: JOSS asks for no institutional affiliation, so this submission can go ahead while the preprint's own affiliation is still open. The flow is one ORCID consent click (ORCID asks "Authorize access" for the client named **Open Journals**, which asks only for the ORCID iD) and then the form: the repository URL, the archive DOI `10.5281/zenodo.24453044` (which already points at the `v2.3.0` tag) and the Summary section as the description; the paper file is `paper.md` in the repository root, and the review takes 1–3 months |
| **Pre-review gate 1 — public development history** | ✅ **met.** JOSS requires the repository to have been public for **more than six months** with development spanning that period: the GitHub repository was created on **2026-04-28**, so the six-month threshold was cleared on 2026-10-28, and the commit history spans April, May, August, September, October and November 2026 (358 commits at the time of writing). JOSS's own check is an automated commit-distribution analysis, so the gap between June and July is visible in it and is worth being able to explain |
| **Pre-review gate 2 — demonstrated research impact** | ⚠️ **partly met, and it moves with the preprint.** JOSS asks for evidence that the software is being used for research, "at minimum by the developers themselves, and ideally by others", and says aspirational statements about future use are not sufficient. What exists today: the developers' own use (the package's artifacts are the evidence base of the methods manuscript), a working integration with the companion pipeline (the extractor is bit-identical, so the numbers describe the same features of the same samples), and public artifacts with DOIs on both ends. What does not: any documented external research use. The paper's research impact statement says exactly that rather than implying otherwise, and the preprint is named there as in preparation. JOSS also asks whether related publications exist — a preprint server is explicitly **not** a previous publication, so the bioRxiv preprint can be posted before, during or after the JOSS review |
| **Pre-review gate 3 — open source practices** | ✅ **met.** For a single-author project JOSS wants *multiple* indicators present at submission time: a meaningful public commit history over time ✅, tagged releases and a changelog ✅ (`v2.3.0`, `CHANGELOG.md`), tests and CI ✅ (1037 tests, two CI jobs), clear documentation ✅ (`README.md`, `MODEL.md`, `USAGE.md`, `docs/`), a `CONTRIBUTING.md` ✅, and stated support or governance expectations ✅ (`CONTRIBUTING.md` § "Support and governance"). JOSS states that a solo project which is otherwise clearly open will not be rejected for lacking a PR workflow |
| **Pre-review gate 4 — iterative development over time** | ✅ **met.** The history shows continued iteration rather than one burst: 31 commits in April, 30 in May, 39 in August, 61 in September, 8 in October and 189 in November 2026, with a tagged release in between (v2.2.0 on 2026-08-21, v2.3.0 on 2026-11-23) |
| **Scope — "must not focus on new research results"** | ✅ **addressed.** JOSS requires the paper not to focus on new research results accomplished with the software. The Summary now says explicitly that the measurements are what the package produces on its own benchmark, reported to show what the software does and how it behaves under its controls, and that the research findings themselves are discussed in the separate methods manuscript |
| **Community engagement — the weakest signal, and it is not a gate** | ⚠️ **no external use yet, and JOSS says so plainly.** The collaborative-effort criterion is scored "Good / OK / Not acceptable", and for a single-author project with no evidence of community engagement the answer is "not acceptable" — although the same paragraph says that for single-author projects, community engagement evidence *in the paper* may be the primary signal, and that a solo project which is otherwise clearly open "will not be rejected solely for lacking a PR workflow". What the repository has: a written collaborator route (`docs/COLLABORATOR_DATA_INTERFACE.md`, the local-cohort adapter, the notice in `README.md`/`TEAM.md`) that is designed but not yet exercised by another group, and an author list of one. JOSS also notes that a single developer working closely with a research group is "meaningfully different from a solo project with no broader community", and that a short author list for a project of this scope is worth raising. Whether a collaborator is credited — and the CUHK connection is documented in `wiki/Jiang-4mer-Validation.md` — is the author's call, not the repository's |

### PLOS Computational Biology / Bioinformatics Advances

| Requirement | State |
|---|---|
| Methods manuscript | ✅ `paper/METHODS_PAPER.md` (all sections filled, numbers regenerated) |
| Data availability statement | ✅ §5 |
| Code availability | ✅ MIT, public, versioned |
| **Acknowledgements section** | ⚠️ **the section now exists in both manuscripts** — `paper/paper.tex` (§ Acknowledgements) and `paper/METHODS_PAPER.md` (§6) thank the public data providers, worded as in the JOSS manuscript's existing section (`paper.md`). **Still open:** whether to add a sentence disclosing that the multi-persona audits were AI-assisted. The section is present so the submission is not desk-reject-incomplete; the wording of that one sentence is the author's call, and two ready-to-paste variants are in the local checklist `paper/BIORXIV_SUBMISSION.md` (gitignored — the choice is not published with the repo) |
| **Generative-AI disclosure** | ⬜ **required by policy — and the section it names now exists.** PLOS's Ethical Publishing Practice (applies to all PLOS journals; https://journals.plos.org/plosone/s/ethical-publishing-practice, checked 2026-11-17): "Contributions by artificial intelligence (AI) tools and technologies to a study or to an article's contents must be clearly reported in a dedicated section of the Methods, or in the Acknowledgements section for article types lacking a Methods section. This section should include the name(s) of any tools used, a description of how the authors used the tool(s) and evaluated the validity of the tool's outputs, and a clear statement of which aspects of the study, article contents, data, or supporting files were affected/generated by AI tool usage." Non-compliance "will be considered misrepresentation of methods, contributions, and/or results", and "PLOS may notify the authors' institution(s)" — which is why the affiliation question above has to be settled *before* a PLOS submission, not after. Three components are required, and this project can answer the third one honestly and specifically: every number in the manuscript is reproduced by the committed code from the committed artifact, and every fix the audits produced is verified by a test. Wording: local checklist §7 |
| Figure files | ✅ `docs/figures/fig1..fig4*.png`, rendered from the committed artifacts |
| **Statistical reviewer concerns** | ⚠️ expect a request for the per-cancer denominators, the single-cohort HCC_J reading and **what stage the held-out cohort is at** — the first two are already stated in §3.3 and §4.2, the third in §3.2.2: the held-out cohorts are advanced-stage where the public records characterise them (publication 1 is Stage IV throughout, publication 9 is metastatic), so the external AUCs are a late-stage transfer result and not an early-detection one |
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
