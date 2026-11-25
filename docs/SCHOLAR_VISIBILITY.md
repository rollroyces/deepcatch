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
| Zenodo `10.5281/zenodo.24453044` | published, v2.3.0, MIT, creator "Lam, Yu Ching" **with the ORCID**, related identifier pointing at the `v2.3.0` tag. This is the citable software identifier, and its metadata is complete. |
| The preprint | **not posted yet** — no bioRxiv or arXiv DOI exists, so the paper itself is not in any index. Everything below depends on this. |
| The repository | the name, the affiliation, the ORCID, the title and the DOI are consistent across the surfaces an index reads, and a guard keeps them that way. |

## 3. The three steps that need the author's accounts

In this order; each one is quick.

### 3.1 Post the preprint

bioRxiv or arXiv — both feed Google Scholar automatically once the
record is public, and both take the ORCID from the manuscript, which
`paper/paper.tex` already carries in the author block. This is the
step that creates the paper's citation, and it is the reason the
manuscript's affiliation has to be decided first (it is copied onto
the record).

bioRxiv needs an institutional affiliation (see §4); arXiv does not,
and the source archive is already built and checked
(`python3 scripts/make_arxiv_bundle.py --verify`).

**The bioRxiv account comes before the ORCID, not after.** ORCID sign-in
only works once an account exists and the iD is associated with it from
the account's Personal Info area: until then the callback answers "No
account is associated with orcid, please login and associate orcid via
the edit personal information page", and the returned URL carries
`orcidcode=0003`. The ORCID side can therefore succeed while the
bioRxiv side still refuses the author — which is what an ORCID
authorization alone looks like. The reason is in the record: the ORCID
exposes **no email address** (ORCID keeps email addresses private by
default, and this one does), so "Create an account with ORCiD" has no
address to register with. The working order is: email → "Check for
Account" → create or sign in → associate the iD → submit. The
affiliation the account carries is the one the preprint must declare,
so that field and §4 are the same decision.

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
| the preprint (bioRxiv) | **Crossref** | bioRxiv mints Crossref DOIs (prefix 10.1101) |
| the preprint (arXiv, if that route is used instead) | **DataCite** | arXiv assigns its own DOIs (`10.48550/arXiv.NNNNNN`) and "the associated article metadata is submitted to DataCite", not Crossref — so the source to pick depends on which server posted the preprint, and the arXiv case is DataCite |

Search for the DOI, confirm the record, and grant the source
permission. Two works, two clicks each.

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
| `CITATION.cff` (and the Zenodo record, which archived the same value) | `Independent Researcher (no institutional affiliation; verified email domain umich.edu)` | yes, exactly |
| `docs/PREPRINT_ABSTRACT.md` (generated from the manuscript) | `Independent Researcher` | yes, exactly |
| `paper.md` (the JOSS paper's author field) | `Independent Researcher` | yes, exactly |
| `MODEL.md` (the model card's author row) | `Independent Researcher` | yes, exactly |
| `paper/METHODS_PAPER.md` (the long-form methods manuscript) | `Independent researcher` (a parenthetical on the author's footnote line) | yes, since 2026-11-24 |
| `paper/PAPER.md` (the earlier panel-only revision) | `Independent researcher` | yes, since 2026-11-24 |
| the bioRxiv / arXiv form values in `paper/BIORXIV_SUBMISSION.md` | `Independent Researcher` | yes, exactly (two rows) |
| the docs-site bylines (`docs/README.md`, `docs/site/index.html`) | `Independent Researcher, Hong Kong` | no — a byline may add a locality |

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
| `TEAM.md` §1.1 | "Lead — Yu Ching Lam (Independent Researcher, no affiliation)", and, under what the author does not have, "the umich.edu verified email on the ORCID profile is from a previous .edu address, not active enrollment" |
| `.github/SPONSORS.md` | "No institution. No lab. No academic email." and "I don't have an institutional affiliation, an academic email, or lab access" |
| `docs/site/index.html` | "an independent researcher in Hong Kong. No institution. No lab. No academic email." |
| `docs/PATH_TO_IMPACT.md` | "Affiliation: \"Independent Researcher\"" |
| `METHYLATION_PROJECT.md` | "**No institutional affiliation** (per user, 2026-09-10)" |
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

bioRxiv's submission guide asks for an affiliation so that there is
an entity to contact about research misconduct, and it does not
consider a job title to be one. The ORCID record still shows one
public employment (University of Michigan, no start or end date), so
"who is the author, and under which affiliation" is a real question
rather than a formatting one — and it is the author's to answer. See
§7 of `paper/BIORXIV_SUBMISSION.md` for the options and what
each one costs.

## 5. What is *not* needed

- A second ORCID, a second Scholar profile, or a variant name: see
  §1.
- A separate "software paper" DOI: the release already has one, and it
  is the same DOI the citation metadata declares.
- Manual entry of the paper into Scholar: the preprint server does
  that when the record goes public, and a manual entry can only
  duplicate it.
