# Making the work findable: Google Scholar, ORCID and the preprint

**Status (2026-11-23):** the metadata side is done and guarded; the
three steps that need the author's own accounts are listed in §3.

This is the part of publishing that is not the paper: a reader who
searches for the work has to be able to *find* it, see it attributed to
one author, and cite it with a permanent identifier. Google Scholar is
the index most readers of a cfDNA paper will use, and it is
unforgiving about the details.

## 1. How Google Scholar actually finds an author

Scholar does not read an ORCID to work out who wrote a paper. It reads
the **name string**, the **affiliation** and the **email** on the
record, and it uses those three to decide whether a paper belongs to a
profile. The ORCID is the second half of the story: it links a work to
a person once the two are known to be the same.

That order matters, because it means:

- A name that differs by a hyphen, an initial or a particle
  ("Yu-Ching Lam" against "Yu Ching Lam") puts half the works in one
  profile and half in another, and the split is invisible until
  somebody merges them by hand.
- The affiliation on the profile and the affiliation on the paper have
  to agree, or the paper is attributed to somebody else with the same
  name.
- An unclaimed profile (or a second one created later) splits the
  record again.

So the repository's job is to make sure that every surface an index
reads says the *same* name, the *same* affiliation and the *same*
ORCID. That is now enforced:
`test/test_author_affiliation_is_consistent.py` checks the name
string and the ORCID across the LaTeX manuscript, the JOSS paper, both
markdown manuscripts, `CITATION.cff` and `AUTHORS`.

## 2. What each source holds today (read from the public APIs)

| Source | State |
|---|---|
| ORCID `0009-0008-9113-769X` | name **Yu Ching Lam** (given "Yu Ching", family "Lam") — identical to every surface above. **The works list is empty**: the Zenodo DOI is not linked to the record, so ORCID-based discovery finds nothing. The record also still carries one public, open-ended employment (University of Michigan), which is the affiliation question in §4. |
| Zenodo `10.5281/zenodo.24453044` | published, v2.3.0, MIT, creator "Lam, Yu Ching" **with the ORCID**, related identifier pointing at the `v2.3.0` tag. This is the citable software identifier. Its one stale detail is the creator's affiliation, which still reads as it stood at that tag — editable in place, DOI unchanged (`docs/JOURNAL_SUBMISSION.md` §4) |
| The preprint | **not posted yet** — no arXiv DOI exists, so the paper itself is not in any index. Everything below depends on this. |
| The repository | the name, the affiliation, the ORCID, the title and the DOI are consistent across the surfaces an index reads, and a guard keeps them that way. |

## 3. The three steps that need the author's accounts

In this order; each one is quick.

### 3.1 Post the preprint

**arXiv — the route.** bioRxiv needs an institutional affiliation (see
§4) and rejected this manuscript on 2026-11-25 under that rule; arXiv
does not, and the source archive is already built and checked
(`python3 scripts/make_arxiv_bundle.py --verify`). It feeds Google
Scholar automatically once the record is public, and it takes the ORCID
from the manuscript, which `paper/paper.tex` already carries in the
author block. This is the step that creates the paper's citation, and it
is the reason the manuscript's affiliation has to be decided first (it
is copied onto the record).

**The arXiv account comes before the ORCID link, not after.** arXiv's own
account page offers "Link my arXiv account with ORCID" — "This linkage
will allow your works on arXiv to be unambiguously connected to your
works in other systems", and "arXiv will use ORCID iDs in preference to
the internal arXiv author identifiers". Until the link exists, an ORCID
authorization alone changes nothing: the iD is what the DataCite
deposit carries, so without it the work reaches the record without an
ORCID. The working order is: register (email, username `rollroyces`,
password, captcha) → verify the email from the message arXiv sends →
link the iD → submit. The affiliation the account carries is the one the
preprint must declare, so that field and §4 are the same decision.

### 3.2 Add the works to the ORCID record

The record has **zero works** at the moment, so the Zenodo DOI is
invisible to anything that discovers through ORCID — including Google
Scholar's "publications from ORCID" panel.

ORCID's own method for this is **Works → +Add → Search & Link**, and
the source differs per work, because the two DOIs come from different
registries:

| Work | Source to pick in Search & Link | Why |
|---|---|---|
| the release (Zenodo) | **DataCite** | DataCite issued the DOI, and ORCID's DataCite integration covers software and datasets |
| the preprint (arXiv) | **DataCite** | arXiv assigns its own DOIs (`10.48550/arXiv.NNNNNN`) and "the associated article metadata is submitted to DataCite", not Crossref. (bioRxiv would have been Crossref, prefix 10.1101 — but it rejected this manuscript, so no bioRxiv DOI was ever minted) |

Search for the DOI, confirm the record, and grant the source
permission. Two works, two clicks each.

**One of the two works can be added today.** The Zenodo record does not
wait for the preprint: v2.3.0 is archived and citable as
`10.5281/zenodo.24453044`, the DOI resolves, and the ORCID record
carries zero works — so the software is invisible to anything that
discovers through ORCID, including Scholar's "publications from ORCID"
panel. Adding it is the cheapest step on this page: one DOI, one "Add
DOI", no trusted-party grant. The preprint's DOI is added the same way
once the server has minted it.

There is a second route, and for the preprint it may be the cheaper one.
ORCID's **Add DOI** ("Works → Add work → Add DOI") accepts a DOI, a
PubMed ID or an arXiv ID, fills the record in from the registry and lets
the author review the imported data before saving it — and it needs **no
trusted-party grant**, so nothing is delegated to a registry. And if the
ORCID iD is supplied *at submission*, the publisher passes it to the
registry, so the work arrives by itself: "all your publications with a
DOI will automatically be added to ORCID as they are published (as long
as you have submitted your ORCID to the publisher during the submission
process)". The submission forms in §3.1 and §3.2 both carry the iD, so
the preprint may not need a manual step at all.

### 3.3 Claim the Google Scholar profile

Create or claim **one** profile for "Yu Ching Lam", put the ORCID in
it (Scholar shows it and links the works), and use the affiliation
from §4 — the same string the preprint carries. Then check that both
works appear under it; if one does not, it is a name/affiliation
mismatch and the fix is in the metadata, not in Scholar.

Do not create a second profile under a variant of the name: two
profiles split the citations, and the merged one has to be curated by
hand afterwards.

## 4. The affiliation, which is the one thing still undecided

Two strings have to agree: the affiliation the preprint declares and
the one the Scholar profile carries. Here is what each surface carries
today, exactly.

**Surfaces that DECLARE the affiliation** — the ones an index reads, and
the ones a reviewer compares against each other:

| Surface | The affiliation it declares | Guarded? |
|---|---|---|
| `paper/paper.tex` (the manuscript's author block, and so the title page of the uploaded PDF) | `Independent Researcher` | yes, exactly |
| `CITATION.cff` (and the Zenodo record, which archived the same declaration) | `Independent Researcher` — and the Zenodo v2.3.0 record carries that same string with a qualifier, `Independent Researcher (no institutional affiliation; verified email domain umich.edu)`, read back from the published record on 2026-11-25; Zenodo lets a published record's metadata be edited without changing its DOI ("You can edit the metadata (title, creators, etc) of a published record at any time… This does not affect your DOI"), so the archive can be brought into line in place — or, if a new release is wanted instead, as a new version with a new DOI | yes, exactly |
| `docs/PREPRINT_ABSTRACT.md` (generated from the manuscript) | `Independent Researcher` | yes, exactly |
| `paper.md` (the JOSS paper's author field) | `Independent Researcher` | yes, exactly |
| `MODEL.md` (the model card's author row) | `Independent Researcher` | yes, exactly |
| `paper/METHODS_PAPER.md` (the long-form methods manuscript) | `Independent Researcher` (a parenthetical on the author's footnote line) | yes, since 2026-11-24 |
| `paper/PAPER.md` (the earlier panel-only revision) | `Independent Researcher` | yes, since 2026-11-24 |
| the bioRxiv / arXiv form values in `paper/BIORXIV_SUBMISSION.md` | `Independent Researcher` | yes, exactly (two rows) |
| the docs-site bylines (`docs/README.md`, `docs/site/index.html`) | `Independent Researcher` | no — a byline may add a locality |

The two markdown manuscripts put the affiliation in prose, so their check
is that the declared string *appears* in the author block, not that the
prose parses to it exactly.

**Surfaces that describe the author's SITUATION** rather than declare an
affiliation on a record. The guard does not compare these, but they must
not contradict the declaration either — a reviewer who reads the paper's
title page and then the landing page must not find two different
institutions:

| Surface | What it says today |
|---|---|
| the public GitHub profile (`rollroyces`) | "ML & bioinformatics engineer — production AI for mRNA design and cfDNA cancer detection. M.S. ADS @ Michigan · bioinformatics @ CUHK · HKT, remote-OK" (read 2026-11-25) — **the one public surface that names an institution in the present tense, and it is not Michigan.** CUHK is where this project's *collaborating* laboratory is (`wiki/Jiang-4mer-Validation.md`: Professor Jiang's laboratory, the source of the 129 real plasma samples), so "bioinformatics @ CUHK" reads as an employment to anyone who does not know that; the Michigan half is the M.S. ADS programme, which the ORCID record carries as an open-dated employment. **No guard compares this** — a GitHub profile is not a tracked file, and it is a public statement of the author's situation |
| `TEAM.md` §1.1 | "Lead — Yu Ching Lam (Independent Researcher)", and "Location: Hong Kong SAR, Independent Researcher" |
| `.github/SPONSORS.md` | "Independent, reproducible, open-source cancer detection. Built in Hong Kong. No paywall. No lock-in." and "I'm an independent researcher in Hong Kong. I don't have a lab of my own or a clinical cohort of my own, and DeepCatch runs entirely on open data." |
| `docs/site/index.html` | "a genomics and bioinformatics student at Independent Researcher. No lab of my own." |
| `docs/PATH_TO_IMPACT.md` | "Affiliation: \"Independent Researcher\"" |
| `METHYLATION_PROJECT.md` | "**Author:** Yu Ching Lam (Independent Researcher)" — its earlier "no institutional affiliation" line was the author's own statement of 2026-09-10 and has been corrected |
| `docs/JOURNAL_SUBMISSION.md` | the gate row that records the decision itself |

The citation metadata's parenthetical is a qualifier *about* the
declaration, not a different declaration, so the declared affiliation is
one string — `Independent Researcher` — everywhere it is declared. The
bylines add a locality because a byline is contact information; a
declared affiliation is the field an index matches on. That is now
enforced rather than asserted:
`test/test_author_affiliation_is_consistent.py` compares the manuscript,
`CITATION.cff`, the generated abstract, the JOSS paper and the model card,
requires the declared string to appear in the two markdown manuscripts'
author blocks, and reads the submission sheet's two form rows and
compares them with the manuscript.
It ignores a trailing parenthetical and nothing else — a locality added
to one surface only is the drift it exists to catch.

bioRxiv's submission guide asks for an affiliation so that there is an entity to contact about research misconduct, and it does not consider a job title to be one. **`Independent Researcher` is therefore not an entity bioRxiv can accept**, so the preprint goes to arXiv, whose registration asks for the same thing and offers "Independent" by name: *"If you are unaffiliated you may enter 'Independent' as your organization."* The preprint's own declaration is `Independent Researcher` on every surface, so the form's organization and the uploaded PDF's author block agree — and arXiv's currency rule ("physical presence, funding, e-mail address") is why the account's email should be a personal address rather than an institutional one.

The ORCID record is the one surface that still says something else: it shows one public employment (University of Michigan, no start or end date), which is the M.S. he has finished and which needs an end date, with the CUHK study belonging in the record's educations section. See §7 of `paper/BIORXIV_SUBMISSION.md` for what is left to do.

**Three public surfaces have to be reconciled, and a fourth is already
fixed.** The ORCID record carries an open-dated University of Michigan
employment; the public GitHub profile reads "M.S. ADS @ Michigan ·
bioinformatics @ CUHK" in the present tense, where CUHK is the
*collaborating* laboratory rather than an employer; and the papers
declare no institution. Only the author can reconcile those three — no
guard can infer it from a document, and a reviewer who compares the
paper's title page with the ORCID record has to find them agreeing.

The fourth is already fixed, and it is the archive. The Zenodo v2.3.0
record is public and immutable, and its creator field carries the
affiliation string quoted above, so the DOI's landing page states the
declaration as well as the paper's title page does. A change of
affiliation after v2.3.0 therefore does not rewrite that page: it reaches
the archive as a new release with a new DOI, and the old DOI keeps the
old string as the historical record — the same rule the repository
already applies to its changelog and its dated audit reports.

## 5. What is *not* needed

- A second ORCID, a second Scholar profile, or a variant name: see
  §1.
- A separate "software paper" DOI: the release already has one, and it
  is the same DOI the citation metadata declares.
- Manual entry of the paper into Scholar: the preprint server does
  that when the record goes public, and a manual entry can only
  duplicate it.
