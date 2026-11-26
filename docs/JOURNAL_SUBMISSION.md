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
        arXiv preprint  ──►  JOSS  (software paper, DOI via Zenodo)
              │                    └─ OR ─►  PLOS Computational Biology
              │                              / Bioinformatics Advances
              └─►  (data + DOI) available for whichever journal follows
```

The preprint is not optional: arXiv assigns the DOI and the priority, it
moderates rather than peer-reviews, and it is what makes the *journal*
submission credible for an author without an institutional affiliation.
Everything in the manuscript is already inside it, so the preprint costs
one upload — the source bundle `dist/deepcatch_v2_3_0.tar.gz`, built and
verified to compile standalone, and the author presses **Submit Article**.

**The route is arXiv.** bioRxiv was the first choice and it is closed: it
rejected this manuscript on 2026-11-25 because it requires an
institutional affiliation, and `Independent Researcher` is not an entity
it will open for (the rejection, in bioRxiv's own words, is in §3's
bioRxiv table). arXiv asks for no institutional affiliation at all — its
own directions say *"If you are unaffiliated you may enter 'Independent'
as your organization"*.

---

## 2. Venue comparison

| | **arXiv** (preprint) | **JOSS** | **PLOS Computational Biology** | **Bioinformatics Advances** | **GigaScience** |
|---|---|---|---|---|---|
| Reviews | moderation, not peer review | software + docs | science | science | science + data |
| Fits this work? | **yes** — priority + DOI, and the only route that needs no institutional affiliation | **yes** — the contribution *is* the reproducible pipeline | yes — the honest external validation and the nulls | yes — methods emphasis | partly — wants more biological interpretation |
| Length | no limit | **750–1,750 words** | no limit | no limit | no limit |
| Abstract length | **≤ 1,920 characters** on the metadata field (the committed abridged text is 1,915); the manuscript's own abstract is unconstrained | no limit | **≤ 300 words** | no limit stated | no limit |
| Data/code required | TeX source bundle + code URL in the manuscript | **public repo, installable, versioned** | data availability statement | data availability | data + workflow |
| Needs clinical cohort | no | no | no (but reviewers will ask) | no | likely |
| Time to decision | 1–4 days of moderation, announced Sun–Thu | 1–3 months | 2–4 months | 1–3 months | 2–4 months |
| DOI | yes (auto, DataCite: 10.48550/arXiv.NNNNNN) | **yes (Zenodo: 10.5281/zenodo.24453044)** | yes | yes | yes |

**bioRxiv was the first choice, and it is closed.** Its own submission
guide requires an affiliation so that "there is an entity providing
oversight we can contact in the event concerns are raised about research
misconduct", and on 2026-11-25 it **rejected this manuscript** for exactly
that reason — "we require submissions to be associated with an organization
that can provide oversight of research activities" (§3's bioRxiv table
carries the rejection in full). `Independent Researcher` is not such an
entity, so the route moved to arXiv, whose registration offers the
unaffiliated author a way in by name.

**Recommendation.** Submit the preprint to arXiv now, then target
**JOSS** with the software framing (the reproducibility guarantees,
the artifact-as-tests design and the honest-null discipline are the
reviewable contribution), and keep **PLOS Computational Biology** as the
parallel science submission if a second journal is wanted.

What the route change costs: bioRxiv could transmit the manuscript files
and its metadata straight to a partner journal ("direct transfer to the
journal", B2J), and **arXiv has no equivalent** — its submission workflow
ends at "Submit Article" (https://info.arxiv.org/help/submit/index.html,
checked 2026-11-26). Each journal form is therefore filled in by hand.
Both are short: the JOSS form asks for the repository URL, the archive DOI
and the paper's Summary section, and the PLOS form asks for the manuscript
re-formatted on PLOS's own template.

Rationale for JOSS specifically: it is the only one of these venues that
*reviews the software*, which is exactly the axis the work is strongest
on. It also requires exactly the things this repository now has — a
versioned public repo, an installable package, a test suite, and a
DOI — and it does not ask for prospective clinical data that no
independent researcher can produce.

---

## 3. Readiness, venue by venue

### arXiv — the route: the affiliation-free server, and what it actually costs

**The route: arXiv** (the author's choice, 2026-11-25). It is the only
venue that does not ask for an institutional affiliation, so it is the
route that does not depend on the affiliation question — which has since
been settled: the declaration the author chose for it is **Independent
Researcher** (see the affiliation rows above). Its requirements are quoted
from arXiv's own submission documentation
(https://info.arxiv.org/help/submit/index.html,
https://info.arxiv.org/help/prep.html, checked 2026-11-25).

| Requirement | State |
|---|---|
| **Account** | ⬜ author action — arXiv's registration asks for email, username, password and a **captcha** (the captcha is the author's to solve). **The registration's own field list, from arXiv's directions (checked 2026-11-25):** username ("must be unique and **cannot be changed**"), name, **organization affiliation**, "**Email affiliated with the organization**", the archives to submit to (this filters the category list later), and a password. The live form's email field is **empty** — the author types his own address, and with the organization set to "Independent" it should be a **personal** one, because arXiv names the e-mail address as one of the things that make a claimed affiliation current. Two details that decide whether this route is open to an unaffiliated author at all: arXiv says in as many words "**If you are unaffiliated you may enter 'Independent' as your organization**", and the password rule is "8 characters or more… must contain at least one underscore (`_`)… **may not contain any other special characters besides underscore**" (arXiv's own example: `Pass_Word_123`) — a password with a hyphen or a dot is rejected. The email has to be **verified** from the account page afterwards, so this is a third step after the form and the captcha |
| **Endorsement** | ⚠️ **possible.** "We only accept submissions from registered authors. If you are a new user or are submitting to a new category, you may be required to find endorsements." A first-time submitter may be asked for an endorsement; a moderator can auto-endorse a submission that is plainly in a category's scope |
| **Abstract length** | ✅ **the limit is on the form field, not on the paper.** "Keep it short — abstracts longer than **1920 characters** will not be accepted; abridge your abstract if necessary" — and that rule sits under arXiv's "Title and Abstract Fields ... for inclusion in metadata" (https://info.arxiv.org/help/prep.html), the text used on the abstract page, in the announcement and in the RSS feed. arXiv does not truncate the manuscript, so **the paper keeps its own 2,940-character abstract** (which bioRxiv also accepts, under its 3,000-character limit) and the *field* gets an abridged one: **1,915 characters**, generated from the same manuscript by `python3 scripts/generate_preprint_abstract.py --variant arxiv` into `docs/ARXIV_ABSTRACT.md`. The generator refuses to build it if it states a number the manuscript's abstract does not contain, so the field cannot drift into a second, different paper. **This corrects an earlier reading of the rule here**, which said the manuscript's abstract "will be rejected" and would have to be abridged |
| **Self-submission** | ⬜ **the author must press "Submit Article".** "Authors are expected to self-submit. Submissions by a third party are only accepted under limited conditions... Prior authorization is required for such submissions. If you are not the author of this work you must contact arXiv administrators for permission. Submissions without prior authorization will be removed." So on arXiv the author submits; what an assistant can do is prepare every field and the source archive |
| **Affiliation must be current** | ✅ **satisfied by the same string.** arXiv's authors-field policy: "It is a violation of our policies to misrepresent your identity or organizational affiliation. Claimed affiliation should be current in the conventional sense: e.g., physical presence, funding, e-mail address, etc. Misrepresentation of identity or affiliation, for any reason, is possible grounds for immediate and permanent suspension." The declared affiliation is `Independent Researcher` and the arXiv account's organization is set to "Independent", so the two agree. Two things follow for the account: its **organization** field is "Independent" (arXiv's own wording: "If you are unaffiliated you may enter 'Independent' as your organization"), and its **email affiliated with the organization** should be a personal address rather than an institutional one, because arXiv names the e-mail address as one of the things that make a claimed affiliation current |
| **Files** | ✅ the source archive `dist/deepcatch_v2_3_0.tar.gz` (flattened TeX plus the three figures), built and verified standalone by `python3 scripts/make_arxiv_bundle.py --verify` and uploaded by the paper workflow as the `deepcatch-arxiv-source` artifact — arXiv takes TeX source, not a PDF. **The artefact CI actually publishes was checked end to end on 2026-11-25:** the archive from the `90df1d9` run was downloaded, unpacked (one `.tex` plus the three flattened figures, no build leftovers), and compiled standalone to 9 pages, byte-for-byte the same size as the committed `paper/biorxiv_submission_v2.3.0.pdf`; its TeX carries the AI-usage disclosure and the current test count |
| **Licence** | ⚠️ **the author's choice, and it is irrevocable.** arXiv offers CC BY 4.0, CC BY-SA 4.0, CC BY-NC-SA 4.0, CC BY-NC-ND 4.0, the arXiv perpetual non-exclusive licence and CC0, and states that "the license chosen is irrevocable and cannot be changed" — and that for every option except CC0 "the original copyright holder retains ownership after posting". Two things to weigh. First, arXiv's own description of the perpetual licence is that it "limits re-use of any type from other entities or individuals", which is the opposite of what this project is for; **CC BY 4.0** is the recommendation, and it is the licence bioRxiv posted under, so the paper's licence does not change with the server. Second, a copyright statement belongs **on the first page of the text**, not in the metadata ("Copyright notices should not be included in the separate metadata and will be removed"). `paper/paper.tex` carries none today, so that line is optional for a plain CC BY 4.0 submission — but it is what arXiv asks for when the paper's own licence is to differ from the dropdown's choice ("authors can select the arXiv license and then indicate the desired license in the first page of the article") or when a funder requires one, and then it means a line in the title block, a rebuilt PDF and a rebuilt bundle |
| **Categories** | ✅ **verified against arXiv's own taxonomy on 2026-11-25** (https://arxiv.org/category_taxonomy): **q-bio.GN** (Genomics) as the primary category — it covers "genomic structure and processes … mutational processes" — and **q-bio.QM** (Quantitative Methods) as the secondary, which covers "All experimental, numerical, statistical and mathematical contributions of value to biology". No computer-science category: the contribution is a biological benchmark and the software is the instrument, not the subject |
| **DOI** | ✅ automatic: "Articles submitted to arXiv are automatically assigned DOIs that correspond to their arXiv ID, and the associated article metadata is submitted to **DataCite**." So an arXiv preprint's DOI is a **DataCite** DOI (`10.48550/arXiv.NNNNNN`), not a Crossref one — which changes the registry to pick when the work is added to the ORCID record (see `docs/SCHOLAR_VISIBILITY.md` §3.2). **And it can reach the ORCID record by itself**: arXiv's own ORCID page says that linking the iD to the account makes your arXiv works "unambiguously connected to your works in other systems", and that "arXiv will use ORCID iDs in preference to the internal arXiv author identifiers" — so with the iD linked to the account **and** supplied in the submission's author field, the DataCite deposit carries it and ORCID adds the DOI-published work on its own. The manual "Add DOI" in §3.2 is then the fallback rather than the plan |
| **Timeline** | ℹ️ moderation's quality-assurance checks "can take between one to four days"; submissions are announced Sunday through Thursday |
| **Generative AI** | ✅ handled, and now checked rather than assumed. "Generative AI language tools should not be listed as an author; instead authors should refer to moderation related policy" — and arXiv's moderation policy asks for the use of text-to-text generative AI to be reported **in the work itself**, so a disclosure kept in a file beside the manuscript would not satisfy it. All four manuscripts carry the disclosure now (`paper/paper.tex`, `paper/METHODS_PAPER.md`, `paper/PAPER.md` and the JOSS `paper.md`), each naming the tools used and confirming that the author reviewed the result, and `test/test_author_affiliation_is_consistent.py` fails if one of them loses it. `paper/PAPER.md` did not carry one until 2026-11-25, and nothing had checked the other three either |

### bioRxiv — closed: the manuscript was rejected on 2026-11-25

**This is the record of the attempt that closed, not an action list.** It is
kept because it is *why* the route is arXiv: bioRxiv requires an
institutional affiliation, and `Independent Researcher` is not an entity it
will open for. Nothing in the table below is pending — the "⬜ author
action" rows are the steps that were planned, and the account was never
created. If the affiliation question is ever reopened, this table is what
to read first; until then, the live route is the arXiv table above.

| Requirement | State |
|---|---|
| Compiled PDF | ✅ `paper/biorxiv_submission_v2.3.0.pdf` (9 pages, 3 figures embedded) — **now the compiled reference copy rather than the uploaded file**, because arXiv takes the TeX source bundle and the PDF is what a reader and the Zenodo release open. **The manuscript itself has not changed since the attempt:** the submission that bioRxiv rejected on 2026-11-25 carried the v2.0-era title "DeepCatch: Performance-Weighted Multi-Modal Fusion with Cumulative Evidence Tracking for Pan-Cancer Detection from Cell-Free DNA", so it was not this manuscript. The current title is "An Open-Source Benchmark for Panel-Based Ultra-Sensitive Detection of Circulating Tumor DNA from Real Tumor Mutations" |
| Abstract within limits | ✅ |
| Submission metadata | ✅ `paper/BIORXIV_SUBMISSION.md` (gitignored: it carries personal details) |
| ORCID | ✅ `0009-0008-9113-769X`, also in `CITATION.cff` |
| Code URL in the manuscript | ✅ |
| **Affiliation consistent with the ORCID record** | ✅ **settled 2026-11-25 — `Independent Researcher`.** The author first stated that he had graduated and was studying genomics and bioinformatics at CUHK, so every surface was changed to that institution on the same day; he then chose to present the work as independent research, and every surface went back to `Independent Researcher` on the same day. The surfaces that declare it: the manuscript's author block (and so the title page of the uploaded PDF), `CITATION.cff`, the JOSS paper, the model card, the two markdown manuscripts, both rows of the submission sheet's form values, and the bylines — and `test/test_author_affiliation_is_consistent.py` holds them to that one string (it failed until every one of them moved). **The ORCID record is the one surface this repository cannot change:** its single employment is "University of Michigan" with no start and no end date, which ORCID renders as *current*, and it belongs to the M.S. he has finished — so it needs an end date (or deletion), and the CUHK study belongs in the record's **educations** section. ORCID needs his sign-in |
| **Affiliation is a registration requirement, not a formality** | ⚠️ **the rule stands, and the manuscript cannot pass it.** bioRxiv's own submission guide says so: "Authors wishing to deposit manuscripts must first register on the site and supply an affiliation so there is an entity providing oversight we can contact in the event concerns are raised about research misconduct" (https://www.biorxiv.org/submit-a-manuscript, checked 2026-11-17). openRxiv's 2026-10-20 statement repeats it: "All bioRxiv and medRxiv authors must therefore be verifiable human beings and affiliated with an organization to which we can turn if there are allegations of fraud or misconduct." **How hard it is in practice — answered, and the earlier answer was wrong:** counted from bioRxiv's own API (`http://api.biorxiv.org/details/biorxiv/2026-01-01/2026-01-31/0/json`, all 5,480 preprints posted in January 2026), **115 of them, 2.1%, are single-authored**, and the earlier version of this row read that as evidence that the rule is not enforced. That inference was wrong, because it measured the wrong population: a single-authored preprint can perfectly well carry an institutional affiliation (one PhD student is a single author), so "single-authored" and "affiliation-less" are not the same thing. **A submission with no institutional affiliation was rejected on 2026-11-25**, and bioRxiv's own words were: *"your manuscript is inappropriate for bioRxiv because bioRxiv requires authors to add an institutional affiliation. As a channel for rapid communication of science that has not undergone formal peer review, we require submissions to be associated with an organization that can provide oversight of research activities."* (manuscript ID BIORXIV/2026/725086; nothing was posted, no DOI was minted, and the record does not appear in bioRxiv's API). The rule is therefore a hard gate — and `Independent Researcher` is not an entity it will open for. **The route is arXiv**, whose registration accepts "Independent" by name; bioRxiv can take the same paper later only if the declaration names an institution |
| **Disclosure of AI use** | ⚠️ **expected, and cheap to satisfy.** bioRxiv's FAQ: "Authors are responsible for the content they submit, including material produced by generative AI and related technologies. As with other computational methodologies, use of these tools should be detailed in manuscripts, and the output must not plagiarize other work or breach copyright." openRxiv's 2026-10-20 statement adds that an author who uses AI "must shape, review, stand behind and have intellectual input into the interpretation of the work they submit", and that a submission whose AI-generated material shows no critical review "may receive further review". The multi-persona audits this project ran are documented in the repository itself (`docs/`, `AUDIT_2_FINDINGS.md`, `CHANGELOG.md`), so a one-sentence disclosure in the acknowledgements is the consistent choice. Ready-to-paste wording: the local (gitignored) checklist `paper/BIORXIV_SUBMISSION.md` §7 |
| **Upload** | ⬜ **closed with the route — the account was never created.** The steps that were planned, in order, verified against the live site on 2026-11-24: **1)** open `submit.biorxiv.org/cgi/registration` and enter the corresponding author's email — bioRxiv checks whether an account exists for that address, so this one field is the whole of step 1; **2)** if no account exists, create one (name, password, affiliation); if one exists, sign in with it — **Bench>Press shares ONE registration across ~70 journals** (medRxiv, the CSH journals, *Genetics*, *Genome Research*, *The Plant Cell* and ~65 more), so an address already registered with any of them is already a bioRxiv account: do not repeat the registration, just sign in with the same email and password and complete the bioRxiv profile page; **3)** in the Personal Info area, associate the ORCID iD `0009-0008-9113-769X` with the account — **until that is done, ORCID sign-in cannot work at all**: bioRxiv's callback answers *"No account is associated with orcid, please login and associate orcid via the edit personal information page"*, which is what the `orcidcode=0003` in the returned URL records, so the ORCID side can succeed while the bioRxiv side still refuses the author; **4)** the author area opens and the manuscript can be submitted. **The root cause is in the ORCID record itself**: it exposes **no email address** (ORCID keeps email addresses private by default — read from the public API on 2026-11-24), so "Create an account with ORCiD" has no address to register with, and the iD has to be attached to an account made with an email address instead. The affiliation given in step 2 is the affiliation the preprint must declare, so it is the same decision as the row above and the repository has to be changed to match it on the same day |
| **The arXiv route, source bundle prepared** | ✅ **built and verified — and it is the route now, not the fallback.** arXiv takes a LaTeX source archive rather than a PDF, and this manuscript reaches outside the repository for three figures (`results/real_tcga_performance.png`, `docs/figures/fig1_pooled_roc.png`, `docs/figures/fig4_sens_operating_points.png`), so a bundle made by copying `paper/paper.tex` alone would compile with three missing figures. `python scripts/make_arxiv_bundle.py --verify` flattens them into `dist/deepcatch_v2_3_0.tar.gz`, compiles it standalone and compares the result with the committed PDF (9 pages, the same). It needs no institutional affiliation — which is exactly the field bioRxiv's registration requires. The paper workflow now uploads the bundle as the `deepcatch-arxiv-source` artifact, so the archive that is actually submitted is the one built from the committed manuscript rather than a local copy that can be stale |
| **Direct transfer to the journal (B2J)** | ⚠️ **this advantage is lost with the route change.** bioRxiv's submission guide: "bioRxiv can save authors time in submitting papers to journals or peer review services by transmitting their manuscript files and metadata directly from the bioRxiv site." The partner list includes **PLOS Computational Biology, PLOS Genetics, PLOS Biology, PLOS Medicine, PLOS One, Bioinformatics, Bioinformatics Advances, GigaScience, eLife, PNAS, Nature and Science**. On bioRxiv the preprint would therefore have been not just a priority record but the submission form, already filled in, for every venue on the list. **arXiv offers no equivalent** — its submission workflow runs from "Prepare Files" through the metadata page to "Submit Article" and stops there (https://info.arxiv.org/help/submit/index.html, checked 2026-11-26) — so each journal form is now filled in by hand |
| **A posted preprint is permanent** | ℹ️ "After an article is deposited in bioRxiv, it is citable and therefore cannot be removed." Revisions keep the same DOI and must carry their supplemental material again, so the PDF uploaded is the version that stays — which is the reason the uploaded PDF is the committed, CI-built `paper/biorxiv_submission_v2.3.0.pdf` and not a locally rebuilt one |
| **Article category** | ✅ **New Results.** bioRxiv's submission guide: "Articles in bioRxiv are categorized as New Results, Confirmatory Results, or Contradictory Results… New Results is most appropriate for articles that describe an advance in a field." This paper is a benchmark that advances the field, so it is New Results — and not Confirmatory (it does not "largely replicate and confirm previously published work") |
| **Licence** | ✅ **the author's choice, and the options are broader than arXiv's.** bioRxiv's terms of distribution: "The copyright remains with the author, who is asked to select from several distribution/reuse options under which to make the article available (CC BY, CC BY-NC, CC BY-ND, CC BY-NC-ND, CC0, or no reuse)." Two further terms in the same paragraph matter for this project's visibility goal: "By posting on bioRxiv, authors explicitly consent to text mining of their work (e.g., by search engines or researchers)" — which is exactly how a preprint gets indexed — and the article "is citable and therefore cannot be removed" |
| **Corresponding-author email** | ⬜ **closed with the route** — and deliberately not in the repo either way |
| **Discoverability (ORCID, Google Scholar)** | ✅ the metadata is consistent and guarded; the three account-side steps are in [`docs/SCHOLAR_VISIBILITY.md`](SCHOLAR_VISIBILITY.md) (post the preprint, add the two works to the ORCID record, claim one Scholar profile) |

### JOSS — one structural gap

| Requirement | State |
|---|---|
| Public, versioned repository | ✅ MIT-licensed, tagged releases |
| Installable package | ✅ `pyproject.toml` |
| Test suite | ✅ 1040 tests (`pytest src test -q`), split across two CI jobs |
| Documentation | ✅ `README.md`, `MODEL.md`, `USAGE.md`, `docs/` |
| Citation metadata | ✅ `CITATION.cff` |
| Contribution guidelines | ✅ `CONTRIBUTING.md` (setup, the drift-guard rules, how to bring a cohort in) |
| Author list | ✅ `AUTHORS` |
| Code of conduct | ✅ `CODE_OF_CONDUCT.md` |
| **DOI** | ✅ **10.5281/zenodo.24453044** (v2.3.0, archived and published; declared in `CITATION.cff`, `MODEL.md` and both manuscripts) |
| **Paper file discoverable from the repo root** | ✅ `paper.md` in the repo root, in the **current** JOSS format: summary, statement of need, **state of the field**, **software design**, features, quality/reproducibility, documentation, code availability, **research impact statement**, limitations, **AI usage disclosure**, acknowledgements. The four bold sections are new: JOSS changed its submission format in late 2025 and its own review bot now checks for each of them, so a paper in the old format is rejected in pre-review. `test/test_author_affiliation_is_consistent.py` requires all twelve headings, and checks that the AI disclosure is a real disclosure (it must name a model provider and confirm the human review), because JOSS's AI usage policy treats an incomplete one as a research-integrity matter |
| **Paper length** | ⚠️ **the paper is 2,509 words and JOSS's current guidance is 750–1,750.** The guidance (https://joss.readthedocs.io/en/latest/paper.html, checked 2026-11-26): "The paper should be between 750-1750 words. Authors submitting papers significantly longer than 1750 words may be asked to reduce the length of their paper." The venue table above recorded "~2,500 words", which was neither the current guidance nor the older one (250–1,000); `paper.md` measures **2,509 words** — front matter, tables and markup excluded — which is **759 words over** the upper bound. It is a request rather than a rejection, but it would cost a revision round on a review that takes 1–3 months. Where the weight sits: Features 399, Software design 351, Summary 340, Research impact statement 265, State of the field 236, References 235, AI usage disclosure 140, Quality and reproducibility 139, Statement of need 126, Limitations 108, Code availability 87. Which 759 words go is the author's |
| **Short description** (the JOSS submission form's one-paragraph field) | ✅ the Summary section of `paper.md` is the one-paragraph description |
| **Submission** | ⬜ author action — the JOSS form, and it is **not** blocked on the affiliation: JOSS asks for no institutional affiliation, so this submission can go ahead on its own terms. The flow is one ORCID consent click (ORCID asks "Authorize access" for the client named **Open Journals**, which asks only for the ORCID iD) and then the form: the repository URL, the archive DOI `10.5281/zenodo.24453044` (which already points at the `v2.3.0` tag) and the Summary section as the description; the paper file is `paper.md` in the repository root, and the review takes 1–3 months |
| **Pre-review gate 1 — public development history** | ✅ **met.** JOSS requires the repository to have been public for **more than six months** with development spanning that period: the GitHub repository was created on **2026-04-28**, so the six-month threshold was cleared on 2026-10-28, and the commit history spans April, May, August, September, October and November 2026 (358 commits at the time of writing). JOSS's own check is an automated commit-distribution analysis, so the gap between June and July is visible in it and is worth being able to explain |
| **Pre-review gate 2 — demonstrated research impact** | ⚠️ **partly met, and it moves with the preprint.** JOSS asks for evidence that the software is being used for research, "at minimum by the developers themselves, and ideally by others", and says aspirational statements about future use are not sufficient. What exists today: the developers' own use (the package's artifacts are the evidence base of the methods manuscript), a working integration with the companion pipeline (the extractor is bit-identical, so the numbers describe the same features of the same samples), and public artifacts with DOIs on both ends. What does not: any documented external research use. The paper's research impact statement says exactly that rather than implying otherwise, and the preprint is named there as in preparation. JOSS also asks whether related publications exist — a preprint server is explicitly **not** a previous publication, so the preprint can be posted before, during or after the JOSS review |
| **Pre-review gate 3 — open source practices** | ✅ **met.** For a single-author project JOSS wants *multiple* indicators present at submission time: a meaningful public commit history over time ✅, tagged releases and a changelog ✅ (`v2.3.0`, `CHANGELOG.md`), tests and CI ✅ (1040 tests, two CI jobs), clear documentation ✅ (`README.md`, `MODEL.md`, `USAGE.md`, `docs/`), a `CONTRIBUTING.md` ✅, and stated support or governance expectations ✅ (`CONTRIBUTING.md` § "Support and governance"). JOSS states that a solo project which is otherwise clearly open will not be rejected for lacking a PR workflow |
| **Pre-review gate 4 — iterative development over time** | ✅ **met.** The history shows continued iteration rather than one burst: 31 commits in April, 30 in May, 39 in August, 61 in September, 8 in October and 189 in November 2026, with a tagged release in between (v2.2.0 on 2026-08-21, v2.3.0 on 2026-11-23) |
| **Scope — "must not focus on new research results"** | ✅ **addressed.** JOSS requires the paper not to focus on new research results accomplished with the software. The Summary now says explicitly that the measurements are what the package produces on its own benchmark, reported to show what the software does and how it behaves under its controls, and that the research findings themselves are discussed in the separate methods manuscript |
| **Community engagement — the weakest signal, and it is not a gate** | ⚠️ **no external use yet, and JOSS says so plainly.** The collaborative-effort criterion is scored "Good / OK / Not acceptable", and for a single-author project with no evidence of community engagement the answer is "not acceptable" — although the same paragraph says that for single-author projects, community engagement evidence *in the paper* may be the primary signal, and that a solo project which is otherwise clearly open "will not be rejected solely for lacking a PR workflow". What the repository has: a written collaborator route (`docs/COLLABORATOR_DATA_INTERFACE.md`, the local-cohort adapter, the notice in `README.md`/`TEAM.md`) that is designed but not yet exercised by another group, and an author list of one. JOSS also notes that a single developer working closely with a research group is "meaningfully different from a solo project with no broader community", and that a short author list for a project of this scope is worth raising. Whether a collaborator is credited — and the CUHK connection is documented in `wiki/Jiang-4mer-Validation.md` — is the author's call, not the repository's |

### PLOS Computational Biology / Bioinformatics Advances

| Requirement | State |
|---|---|
| Methods manuscript | ✅ `paper/METHODS_PAPER.md` (all sections filled, numbers regenerated) |
| **Manuscript format** | ⚠️ **the preprint's PDF is not a PLOS manuscript, and the repository had not recorded that.** PLOS's own LaTeX page (https://journals.plos.org/ploscompbiol/s/latex, checked 2026-11-25): "PLOS provides a template for LaTeX submissions… Follow guidelines in the template to make sure it is properly formatted for the PLOS system", and "**If your submission does not follow the template format, your manuscript may be sent back to you for revisions.**" The same page fixes the shape of the initial submission: "upload the PDF as your manuscript file, containing your manuscript text, **figure legends**, tables, and references. **Do not include figures in your PDF; they should be uploaded as separate files.**" So the PLOS manuscript needs its own variant of `paper/paper.tex` on the PLOS template, with the figure **legends as paragraphs in the text** — they are `\caption{}`s under embedded figures today — and the three figures uploaded separately as TIFF/EPS (the format row below). That work is queued behind the preprint, not done |
| **Bioinformatics Advances' own format rules** | ⚠️ **three requirements the generated `.tex` does not yet meet.** The journal's author guidelines (https://academic.oup.com/bioinformaticsadvances/pages/author-guidelines, checked 2026-11-26, read in a browser because the page refuses a plain fetch): "**You must place all tables at the end of the main text**"; "**You must include figure titles and legends within the manuscript file — they should not be included in the image file**"; and the initial submission is "the PDF in the first instance for review featuring any tables, figures and references alongside the main text", with the figures submitted as individual image files at the production stage. The generated manuscript has its tables inline rather than at the end, and its figures are embedded rather than submitted separately. **The abstract has no word limit** on that page, so the 300-word PLOS cap binds one target venue and not the other |
| Data availability statement | ✅ §5 |
| Code availability | ✅ MIT, public, versioned |
| **Abstract length** | ⚠️ **the methods manuscript's abstract is 414 words, and PLOS's limit is 300.** PLOS's own submission guidelines (https://journals.plos.org/ploscompbiol/s/submission-guidelines, checked 2026-11-26): "The Abstract should be succinct; **it must not exceed 300 words**." The methods manuscript's abstract measures **436 words**, so on that ground alone it would come back for revision. The same page is explicit that the *manuscript* is unconstrained -- "Manuscripts can be any length. There are no restrictions on word count, number of figures, or amount of supporting information" -- so the 300 words binds the abstract alone, not the paper. Abridging it belongs to the PLOS-variant work that is queued behind the preprint, and which words to cut is the author's. The count is **measured, not remembered**: `python3 scripts/count_manuscript_abstract.py` reports it and `--check` fails if any document states a different one — the figure had been recorded as 414 words, which was wrong when it was written (the abstract has not changed since), and the gap is **136 words**, not 114. **The preprint's own abstract is unaffected**: arXiv's limit is on the metadata field's characters (1,920, with a committed 1,915-character text ✓) and the manuscript's own abstract is unconstrained |
| **Acknowledgements section** | ⚠️ **the section now exists in both manuscripts** — `paper/paper.tex` (§ Acknowledgements) and `paper/METHODS_PAPER.md` (§6) thank the public data providers, worded as in the JOSS manuscript's existing section (`paper.md`). **Still open:** whether to add a sentence disclosing that the multi-persona audits were AI-assisted. The section is present so the submission is not desk-reject-incomplete; the wording of that one sentence is the author's call, and two ready-to-paste variants are in the local checklist `paper/BIORXIV_SUBMISSION.md` (gitignored — the choice is not published with the repo) |
| **Generative-AI disclosure** | ⬜ **required by policy — and the section it names now exists.** PLOS's Ethical Publishing Practice (applies to all PLOS journals; https://journals.plos.org/plosone/s/ethical-publishing-practice, checked 2026-11-17): "Contributions by artificial intelligence (AI) tools and technologies to a study or to an article's contents must be clearly reported in a dedicated section of the Methods, or in the Acknowledgements section for article types lacking a Methods section. This section should include the name(s) of any tools used, a description of how the authors used the tool(s) and evaluated the validity of the tool's outputs, and a clear statement of which aspects of the study, article contents, data, or supporting files were affected/generated by AI tool usage." Non-compliance "will be considered misrepresentation of methods, contributions, and/or results", and "PLOS may notify the authors' institution(s)" — which is why the affiliation question above has to be settled *before* a PLOS submission, not after. Three components are required, and this project can answer the third one honestly and specifically: every number in the manuscript is reproduced by the committed code from the committed artifact, and every fix the audits produced is verified by a test. Wording: local checklist §7 |
| Figure files | ⚠️ **the submission is fine; three production details are not.** PLOS's own figures page (https://journals.plos.org/ploscompbiol/s/figures, checked 2026-11-25): "**PLOS Computational Biology waives all formatting requirements until your manuscript has received a provisional Editorial Accept decision.**" So the PNGs are acceptable at submission, and the figures exist and are rendered from the committed artifacts ✓ (`results/real_tcga_performance.png`, `docs/figures/fig1_pooled_roc.png`, `docs/figures/fig4_sens_operating_points.png`). Measured against PLOS's Figure File Requirements: the pixel dimensions are **inside** the spec (widths 2085, 1198 and 1088 px against 789–2250 ✓, heights ≤2625 ✓, 53–64 KB against <10 MB ✓), and as placed in the manuscript — at 0.48 and 0.62 of a 6.5-inch text width — the effective resolution is **349–517 dpi**, above PLOS's 300 ✓. Three things still differ: the format must become **TIFF or EPS** (these are PNG), the file metadata must say **300 dpi** (the PNGs tag 110 and 150 — the pixels are fine, the tag is not), and text inside a figure must be **Arial, Times or Symbol, 8–12 pt** (the plots use matplotlib's default DejaVu Sans at 9–10 pt). The captions are already in the manuscript rather than in the image files ✓. None of the three is a submission gate |
| **Statistical reviewer concerns** | ⚠️ expect a request for the per-cancer denominators, the single-cohort HCC_J reading and **what stage the held-out cohort is at** — the first two are already stated in §3.3 and §4.2, the third in §3.2.2: the held-out cohorts are advanced-stage where the public records characterise them (publication 1 is Stage IV throughout, publication 9 is metastatic), so the external AUCs are a late-stage transfer result and not an early-detection one |
| **DOI** | ⬜ same Zenodo dependency |

---

## 4. The software DOI — minted, and what it still carries

**Zenodo's GitHub integration.** Every venue above wants a DOI for the
software, and Zenodo only mints one automatically if the repository is
linked to the account. That is done: `v2.3.0` is archived and published
as **10.5281/zenodo.24453044**, and the DOI is declared in `CITATION.cff`,
`MODEL.md` and both manuscripts.

One detail in that archive is out of step with the repository today.
Zenodo's GitHub integration takes the author block from `CITATION.cff` at
tag time — the `v2.3.0` record's creator block matches `CITATION.cff`
**at that tag** exactly, in name, affiliation and ORCID — so the archive
carries the affiliation as it stood then: *"Independent Researcher (no
institutional affiliation; verified email domain umich.edu)"*. Every
surface today declares the shorter **Independent Researcher**, and that
parenthetical names a domain which has nothing to do with the
declaration any more.

There are two ways to close that, and either is the author's account
action:

1. **Edit the record's creator metadata in place.** Zenodo's own help
   says it plainly: "You can edit the metadata (title, creators, etc) of
   a published record at any time… This does not affect your DOI"
   (https://help.zenodo.org/docs/deposit/manage-records/, checked
   2026-11-26). The DOI stays `10.5281/zenodo.24453044` and every
   citation of it stays correct.
2. **Cut a new tag** and let Zenodo archive a new version, whose DOI
   differs; the new DOI then has to replace the old one in
   `CITATION.cff`, `MODEL.md` and the manuscripts. The steps are in
   `ZENODO_SETUP.md`, and the version number is the author's call.

Neither blocks the preprint: arXiv issues its own DOI automatically (a
DataCite DOI, `10.48550/arXiv.NNNNNN`). It matters for the journal
submissions, which cite the software's permanent identifier.

---

## 5. What a reviewer will press on, and the answer

| Likely objection | Where the answer is |
|---|---|
| "The pooled AUC is measured inside the cohort that trained it." | §3.2.1 and the external validation in §3.2.2: 0.9444 on a held-out whole-genome publication, 0.708 on a second, 0.6805 in the single-type direction, with the reason spelled out. |
| "Is the methylation channel real?" | §3.2.4: it is a **null**, with the mechanism (size selection at 300 bp, 1.09 observations per CpG) and a control (depth ablation on a genuine MspI library). |
| "Are the numbers reproducible?" | §2.5: the PCA is seeded, and `test/test_benchmark_reproducibility.py` asserts that two runs agree exactly. |
| "Do the artifacts match the text?" | `scripts/benchmark_numbers.py` prints every quoted figure straight from the committed JSON. |
| "No clinical cohort." | §4.2 and §4.3: this is an open-data methods benchmark, not a clinical validation, and it says so. |
