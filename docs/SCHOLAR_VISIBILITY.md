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

Search for the DOI, confirm the record, and grant the source
permission. Two works, two clicks each.

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
today, exactly:

| Surface | The affiliation it declares |
|---|---|
| `paper/paper.tex` (the manuscript's author block, and so the title page of the uploaded PDF) | `Independent Researcher` |
| `CITATION.cff` (and the Zenodo record, which archived the same value) | `Independent Researcher (no institutional affiliation; verified email domain umich.edu)` |
| `docs/PREPRINT_ABSTRACT.md` (generated from the manuscript) | `Independent Researcher` |
| `paper.md` (the JOSS paper's author field) | `Independent Researcher` |
| the bioRxiv / arXiv form values in `paper/BIORXIV_SUBMISSION.md` | `Independent Researcher` |
| the docs-site bylines (`docs/README.md`, `docs/site/index.html`) | `Independent Researcher, Hong Kong` |

The citation metadata's parenthetical is a qualifier *about* the
declaration, not a different declaration, so the declared affiliation is
one string — `Independent Researcher` — everywhere it is declared. The
bylines add a locality because a byline is contact information; a
declared affiliation is the field an index matches on. That is now
enforced rather than asserted:
`test/test_author_affiliation_is_consistent.py` compares the manuscript,
`CITATION.cff`, the generated abstract and the JOSS paper, and reads the
submission sheet's two form rows and compares them with the manuscript.
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
