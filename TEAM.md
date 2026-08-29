# Project Team — DeepCatch + cfdna-fragmentomics-pipeline

This document describes the **current team structure**, the **open
collaboration opportunities**, and the **governance model** for the
DeepCatch + cfdna-fragmentomics-pipeline project. It is intentionally
honest: as of August 2026, this is a **single-contributor project**
with explicit open roles. The plan is to grow it into a small
collaborative team as the right people appear.

> **Not affiliated with any institution, company, or grant.** The
> project is independent, MIT-licensed, open-source, and depends on
> volunteer contributors + public data.

---

## 1. Current team (as of 2026-08-28)

### 1.1 Lead — Royce (Independent Researcher, Hong Kong SAR)

**Role**: project founder, sole author, sole maintainer, sole
contributor of every commit across both repos as of 2026-08-28.

**What Royce does**:
- Writes all code (Python, shell, LaTeX)
- Designs the experimental methodology
- Runs the benchmarks
- Writes the documentation (READMEs, RESULTS, AUDIT_REPORT_2, USAGE)
- Maintains the bioRxiv submission package
- Manages the GitHub issues, PRs, releases

**What Royce does NOT have**:
- Wet-lab access (this is a computational-only project)
- Clinical collaborators with patient samples
- An ORCID (TODO — 5 min on https://orcid.org/register)
- Institutional affiliation, funding, or computing cluster
- A co-author for clinical-impact journals

**Contact** (placeholder — fill in before submission):
- bioRxiv profile: pending
- ORCID: pending
- Email: `[your email]`
- GitHub: `@rollroyces`
- Country: Hong Kong SAR

### 1.2 Indirect contributors (acknowledged but not yet authors)

These people / projects have enabled this work but have **not
contributed code, design, or analysis**:

- **Cristiano et al. 2019** (Nature, GSE317007) — 537 of the 627
  headline samples. Their preprocessing pipeline is reused (via
  FinaleDB).
- **Jiang et al. 2015** (PNAS) — 121 of the 627 headline samples
  (HCC + healthy low-pass WGS).
- **FinaleDB team** (CCHMC, Liu lab) — public host of uniformly-
  preprocessed cfDNA WGS data, source of all 627 samples.
- **TCGA-LUAD consortium** (NCI) — source of 20 patients' somatic
  mutations used in the panel LLR demo.
- **OpenAI / Anthropic / Google** — the LLM providers behind the
  AI agents that have been doing 4+ rounds of multi-persona audits.

These are acknowledged in [`paper/PAPER.md`](paper/PAPER.md) and
[`BENCHMARK.md`](BENCHMARK.md) citations, not as authors.

---

## 2. Open collaboration roles (the gaps)

The project is at a stage where **specific skills are missing** for it
to grow beyond a single-author benchmark. The roles below are
**what we need**, not what we have. Each role has a one-paragraph
description of the contribution, a time estimate, and a contact
mechanism.

### 2.1 Clinical / wet-lab co-author — **HIGHEST PRIORITY**

**What we'd want from you**: real plasma cfDNA data from patients
with confirmed cancer diagnoses (preferably pan-cancer, but any
tumor type useful). The data should be:
- WGS or targeted panel sequencing
- Matched tumor + plasma (for the fusion story to be honest)
- 30-100 samples (any size is useful; we just need held-out clinical
  data to upgrade from "honest benchmark" to "validated assay")
- Public-accessible OR with patient consent for our re-analysis

**Why this is highest priority**: The round-4 journal reviewer
explicitly said the work needs "real held-out clinical data" before
Nature Medicine / Cancer Discovery will accept it. Single-author
benchmarks on public data are appropriate for *Bioinformatics* or
*PLOS Comp Bio* but not for clinical-impact journals.

**Time commitment**: a one-time data drop + one Skype call to align
on the re-analysis.

**Compensation**: co-authorship on the next published version of the
work. (If we get to a clinical journal, this is a significant
co-authorship.) No financial compensation — this is a volunteer
research project.

**Contact**: open a GitHub Issue tagged `clinical-data` or email the
address in `paper/BIORXIV_SUBMISSION.md`.

### 2.2 Methylation / Galleri-style expert — **HIGH PRIORITY**

**What we'd want from you**: help interpreting where the
fragmentomics-only approach (this work) sits relative to the
methylation-based MCED work (Galleri, PATHFINDER, Liu 2020). The
round-4 statistical reviewer noted that "methylation is the largest
single signal in any MCED context" — and our framework already has
a [`src/methylation_gnn/`](../../deepcatch/src/methylation_gnn/)
module that is smoke-tested but not validated.

**Time commitment**: 2-4 hours per month for 3 months, to:
- Audit the existing methylation GNN scaffolding
- Specify what a head-to-head comparison would require
- Co-author the comparison paper (if/when we do it)

**Contact**: open a GitHub Issue tagged `methylation` or comment on
any of the `src/methylation_gnn/*.py` files.

### 2.3 Statistician / methodological reviewer — **MEDIUM PRIORITY**

**What we'd want from you**: peer-review the statistical methodology
in [`RESULTS.md`](RESULTS.md) and [`AUDIT_REPORT_2.md`](AUDIT_REPORT_2.md),
particularly:
- The 4 round-4 statistical caveats (Q1 Bonferroni, Q2 winner's-curse,
  Q3 fusion calibration, Q5 PPV prevalence interpretation)
- The DeLong CI computation (`deepcatch/validation/delong_test.py`)
- The per-fold NaN imputation in `cfdna-fragmentomics-pipeline/scripts/train_classifier.py`

**Time commitment**: 1-2 hours per quarter, async via GitHub PR
reviews. No required meetings.

**Contact**: comment on any `AUDIT_REPORT_2.md` section or open an
issue tagged `stats-review`.

### 2.4 Software engineer (Python / sklearn / pandas) — **MEDIUM PRIORITY**

**What we'd want from you**: code-review the next 100 PRs. The
project is mostly < 2000 LoC per script but it has accumulated
technical debt from being a one-person project for 8+ months:
- 5 scripts still use the old `argparse` style (no `--help`)
- 2 scripts have hardcoded paths (mostly fixed in round-3 audit)
- Test coverage of script-level integration is 29% (12 of 17 scripts
  in pipeline; 3 of 12 in deepcatch/src/fragmentomics)
- The `model_ablation.py` script in pipeline has a module-level
  execution bug (pre-existing, deferred)

**Time commitment**: 1-2 hours per week, async.

**Contact**: open an issue tagged `code-review` or comment on a PR.

### 2.5 Wet-lab / sequencing collaborator — **LOW PRIORITY (but transformative)**

**What we'd want from you**: real WGS or panel sequencing of cfDNA
samples. This is different from §2.1 (which is about analyzing
existing data); §2.5 is about **producing new data**.

**Why this is transformative**: the project currently relies on
public FinaleDB data. To upgrade the panel LLR demo from
"simulated plasma from TCGA mutations" to "real plasma with
real variant calling," we need real WGS data from at least 5-10
patient samples.

**Time commitment**: a one-time sequencing run (~1-2 weeks of lab
time + ~$5,000-10,000 USD in reagents, depending on platform).

**Compensation**: senior co-authorship on the next paper; you
retain all rights to your raw data and may publish it elsewhere
in parallel.

**Contact**: open an issue tagged `wet-lab` or email the
corresponding-author address.

### 2.6 Documentation / technical writer — **LOW PRIORITY**

**What we'd want from you**: copy-edit the READMEs, USAGE.md, and
bioRxiv submission package for grammar, clarity, and tone. The
project currently has 3500+ lines of user-facing documentation
written by one person (Royce); an outside editor would catch
issues the author is too close to see.

**Time commitment**: 1-2 days of async review.

**Contact**: open a PR against `USAGE.md` or `paper/BIORXIV_SUBMISSION.md`.

### 2.7 Bioinformatician student / postdoc — **LOW PRIORITY**

**What we'd want from you**: a 3-6 month project using this framework
to answer a research question that Royce hasn't had time to address.
Suggested topics:
- Per-cancer-type performance (the §6.1 deferred item)
- Methylation GNN head-to-head vs the fragmentomics baseline
- ComBat / limma-style harmonization comparison
- Held-out validation on a real cohort (requires data access)

**Time commitment**: 3-6 months, part-time or full-time.

**Compensation**: first-author on the resulting paper. (Royce
co-authors; your advisor co-supervises.)

**Contact**: open an issue tagged `student-project` or email the
corresponding-author address with a CV and a 1-page project
proposal.

---

## 3. Contribution model

The project follows a **fork-and-pull-request** model on GitHub:

1. **Fork** the relevant repo (`rollroyces/deepcatch` or
   `rollroyces/cfdna-fragmentomics-pipeline`).
2. **Make your changes** in a topic branch.
3. **Open a Pull Request** with a clear description of:
   - What changed
   - Why it changed
   - How you tested it (commands, expected output, actual output)
4. **CI must pass** (`.github/workflows/*.yml` runs automatically on
   every PR).
5. **A maintainer will review** within 1-2 weeks.
6. **Once approved**, the PR is squash-merged to `main`.

There is no CLA (Contributor License Agreement). All contributions
are MIT-licensed under the same terms as the rest of the project.

### 3.1 What kind of contributions are welcome

- Bug fixes
- Documentation improvements
- New feature extractors
- New model architectures (especially non-linear: RF, LightGBM,
  XGBoost, kernel SVM)
- Held-out validation experiments
- Reproducibility improvements (CI smoke tests, dependency pinning)
- Methylation GNN scaffolding completion
- Cross-platform fixes (Windows, ARM)
- Translations of USAGE.md / README.md

### 3.2 What is **NOT** accepted

- Anything that requires private data the contributor cannot share
  (the project's reproducibility depends on public data only)
- Anything that introduces proprietary dependencies
- Anything that requires GPU/TPU compute (the framework is designed
  to run on CPU; deep-learning code paths are smoke-tested only)
- Rebranding or forking under a different name
- Adding new top-level repositories without discussion

### 3.3 What is **ACCEPTED BUT REQUIRES DISCUSSION**

- Changing the headline numbers (any change to `RESULTS.md` TL;DR)
- Adding new headline claims (e.g., claiming a new best-AUC)
- Changing the license
- Changing the maintainer structure (this document)

These need an issue with the `governance` tag, and 2 weeks of
discussion before any change.

---

## 4. Governance

The project is a **single-maintainer project with a community-input
governance model**. Specifically:

### 4.1 Decision authority

- **Code changes** (anything that doesn't touch the headline
  numbers, license, or maintainer): Royce has final say after
  one round of community review.
- **Headline changes** (changes to `RESULTS.md` TL;DR, `BENCHMARK.md`
  top section, `paper/PAPER.md` abstract, `paper/biorxiv_submission_v2.2.0.pdf`):
  requires Royce approval + 2 weeks of public discussion.
- **License changes**: requires Royce approval + 1 month of
  public discussion + 2 community-member agreement.
- **Maintainer changes**: requires Royce approval + a clear
  handoff plan documented in `TEAM.md`.

### 4.2 What is NOT in the governance

The project is **not**:
- A foundation or non-profit
- A working group or consortium
- Affiliated with any institution
- Bound by any external code-of-conduct (we follow the GitHub
  community guidelines: <https://docs.github.com/en/site-policy/github-terms/github-community-guidelines>)

If the project grows to need any of these, the governance model
will be revisited.

### 4.3 Conflict resolution

If Royce and a contributor disagree on a change:
1. Try to resolve in the PR comments.
2. If unresolved, open a `governance` issue with both positions.
3. After 2 weeks of public discussion, Royce has tiebreaker authority.
4. If the contributor feels Royce is being unreasonable, they
   may fork the project under the MIT license (which they are
   free to do).

---

## 5. Roadmap (6 months)

If the right collaborators appear, here's what could be done
in the next 6 months:

| Quarter | Goal | Required collaboration |
|---|---|---|
| **Q1 (now)** | bioRxiv preprint (current state) | None — Royce is preparing |
| Q1 | ORCID + Zenodo deposit (manual) | None — Royce |
| Q2 | Per-cancer-type AUC table | §2.4 software engineer (data prep) |
| Q2 | Head-to-head vs Galleri (methylation) | §2.2 methylation expert |
| Q3 | Held-out clinical validation | §2.1 clinical co-author |
| Q3 | ComBat / limma harmonization comparison | §2.3 statistician |
| Q4 | Methods-journal submission (Bioinformatics) | All of the above |
| Q4 (stretch) | Clinical-journal submission | §2.5 wet-lab collaborator |

The roadmap is **aspirational, not committed**. If the
collaborators don't appear, the project will continue as-is
(single-author benchmark on public data) and re-target to
Bioinformatics or PLOS Computational Biology.

---

## 6. How to engage

### 6.1 If you want to ask a question

- **GitHub Issues** (preferred): open an issue on either repo with
  the `question` tag.
- **GitHub Discussions** (if enabled): for open-ended conversation.
- **Email**: for sensitive or off-the-record questions, use the
  address in `paper/BIORXIV_SUBMISSION.md` (after you fill in
  the placeholder).

### 6.2 If you want to report a bug

- **GitHub Issues** (preferred): open an issue on the relevant repo
  with the `bug` tag. Include:
  - The command you ran
  - The expected output
  - The actual output
  - The Python version (`python --version`)
  - The OS (e.g., `Darwin 24.0.0`)

### 6.3 If you want to contribute code

- See [§3 Contribution model](#3-contribution-model) above.
- First step: open an issue with the `contribution` tag describing
  what you want to do. This avoids wasted work on a PR that gets
  rejected for design reasons.

### 6.4 If you want to be a co-author (not just a contributor)

- See [§2 Open collaboration roles](#2-open-collaboration-roles)
  above. The 7 roles there are the only paths to co-authorship.
- Co-authorship is granted after substantive contribution (≥ 1
  research-question-level contribution, not just bug fixes).

---

## 7. Current state (audit summary, August 2026)

| Asset | Status | Documented in |
|---|---|---|
| Headline tumor-naive AUC 0.974-0.978 | ✓ Verified, 4 audit rounds | `RESULTS.md` |
| Headline fusion AUC 0.989 | ✓ Verified, partly synthetic | `RESULTS.md` §3 |
| bioRxiv preprint | ⏳ Awaiting user action (ORCID + email) | `paper/BIORXIV_SUBMISSION.md` |
| 85 unit tests passing | ✓ | `test/` in both repos |
| 9/9 CI jobs green | ✓ | `.github/workflows/` |
| 4 audit rounds completed | ✓ | `~/JOURNAL_REVIEW_*`, `~/ROUND4_*`, repo `AUDIT_REPORT*.md` |
| 3500+ lines of user documentation | ✓ | `USAGE.md` + `README.md` + `RESULTS.md` + audit reports |
| Cross-repo fusion (DeepCatch + pipeline) | ✓ Working | `src/fragmentomics/fusion_ablation.py` |
| Tumor-naive adapter (file-format contract) | ✓ Tested in CI | `src/fragmentomics/tumor_naive_adapter.py` |

### 7.1 What is **NOT** done (the honest list)

- ❌ Held-out clinical validation (needs clinical co-author)
- ❌ Per-cancer-type AUC (needs FinaleDB metadata or new cohort)
- ❌ Methylation channel (scaffolding only)
- ❌ Real mutation channel (synthetic only)
- ❌ ComBat / limma harmonization comparison
- ❌ Production deployment (regulatory framing not addressed)
- ❌ Single-source funding
- ❌ Second author

---

## 8. History

| Date | Event |
|---|---|
| 2026-04-28 | `data/deepcatch_data.xlsx` removed from repo for privacy (commit 8c812c0) |
| 2026-08-18 | FinaleDB 627-sample feature set extracted (407 MB, gitignored) |
| 2026-08-22 | First bioRxiv submission attempt blocked on user action |
| 2026-08-23 | Round-1 audit (4 reviewers): RESULTS.md caveats |
| 2026-08-25 | Round-2 audit: PPV, 8ch, --help bug, file-handling |
| 2026-08-26 | Round-3 audit: E1 (nuc_ablation dup), E4 (NaN), E6 (Gemma), E8 (paths) |
| 2026-08-28 | Round-4 audit: Q4 (decision-curve), Q5 (PPV prevalence), B3 (README), B5 (traceability), S4 (LICENSE), S7 (workflow) |
| 2026-08-28 | USAGE.md (510-520 lines) added to both repos |
| **2026-08-28** | **This TEAM.md (project team structure documented)** |
| 2026-08-28 | bioRxiv preprint submission: **2 user actions pending (ORCID + email)** |
| 2026-Q3 (planned) | If collaborators appear: per-cancer-type AUC, head-to-head vs Galleri, held-out validation |
| 2026-Q4 (planned) | Methods-journal submission to Bioinformatics or PLOS Comp Bio |

---

*Last updated: 2026-08-28* (round-4 audit consolidation + USAGE.md + TEAM.md)

*See git log for the precise commit.*
