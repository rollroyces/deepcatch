#!/usr/bin/env python3
"""Print the benchmark numbers that the manuscript and docs quote.

One place to read them from, so the prose and the artifacts cannot drift
apart: the script reads the committed ``results/*.json`` and prints every
figure the documents quote, with the artifact it came from.

``QUOTED_NUMBERS`` is the authoritative list -- artifact, the path to the
number inside it, and the documents that must contain it.
``test/test_paper_numbers_match_artifacts.py`` imports that list and
fails the build when an artifact is regenerated and the prose is not
updated, so the two never disagree about what is quotable.

What this list does not cover, established by an audit on 2026-11-25 that
checked every four-decimal figure in the tracked documents against the
numbers the committed artifacts hold: the figures it does not pin are the
companion repositories' and an external deposit's, or derived values (a
difference between two runs, a DeLong bound) -- traceable, but not from
this tree. Two blind spots are worth naming. Some figures live in binary
artifacts (``results/alphamissense_real_index_full.pkl`` and
``results/alphamissense_index_cohort.pkl`` hold the AlphaMissense figures
``docs/ALPHAMISSENSE_REAL_LLR.md`` quotes), which no diff shows and no
check here reads. And a figure that no artifact holds at all is not
caught by anything in this file: the audit found two such figures on the
project's wiki -- a nested-CV AUC and an optimal-k AUC whose run's input
and output were never committed -- and they are now labelled as such on
the page rather than presented as reproducible.
"""
from __future__ import annotations

import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Where the artifacts live. Overridable so that a reproduction run writing
# into a scratch directory can be checked against the committed artifacts
# (see ``scripts/reproduce_check.py``) without touching the repo.
RESULTS_DIR = os.environ.get("DEEPCATCH_RESULTS_DIR") or os.path.join(
    ROOT, "results")

# (artifact file, path into the artifact, documents that must quote it)
QUOTED_NUMBERS = [
    ("external_validation_finaledb.json",
     ("per_publication_auc", "7", "auc"),
     ["paper/paper.tex", "paper.md", "README.md", "README.zh-CN.md",
      "MODEL.md", "RESULTS.md", "TEAM.md", "CITATION.cff",
      "docs/EXTERNAL_VALIDATION.md", "docs/JOURNAL_SUBMISSION.md",
      "docs/PREPRINT_ABSTRACT.md", "docs/site/index.html",
      "paper/README.md", "results/README.md"]),
    # The second held-out whole-genome publication, and the case-only
    # one reported as sensitivity at the training threshold.
    # The second held-out whole-genome publication. Every document that
    # quotes it is named, including the JOSS manuscript `paper.md`, whose
    # summary quotes 0.708 -- it was pinned for the methods paper and
    # not for the one a JOSS reviewer reads.
    ("external_validation_finaledb.json",
     ("per_publication_auc", "1", "auc"),
     ["paper/METHODS_PAPER.md", "paper.md", "paper/paper.tex",
      "paper/README.md", "paper/REPRODUCE.md", "MODEL.md", "README.md",
      "README.zh-CN.md", "RESULTS.md", "REVIEWERS.md", "TEAM.md",
      "CITATION.cff", "docs/ARXIV_ABSTRACT.md", "docs/EXTERNAL_VALIDATION.md",
      "docs/JOURNAL_SUBMISSION.md", "docs/PUBLIC_DATA_SOURCES.md",
      "docs/PREPRINT_ABSTRACT.md", "docs/site/index.html",
      "results/README.md"]),
    # The lower bound of that publication's DeLong interval. The AUC is
    # pinned; the bound was not, so a document could carry a stale
    # interval and nothing would notice.
    ("external_validation_finaledb.json",
     ("per_publication_auc", "7", "ci_lo"),
     ["paper/METHODS_PAPER.md", "paper.md", "paper/paper.tex",
      "paper/README.md", "paper/REPRODUCE.md", "MODEL.md", "README.md",
      "README.zh-CN.md", "RESULTS.md", "REVIEWERS.md", "CITATION.cff",
      "docs/ARXIV_ABSTRACT.md", "docs/EXTERNAL_VALIDATION.md",
      # The ccRCC section now quotes the held-out publication's DeLong
      # lower bound, to say which interval belongs to which cohort.
      "docs/PUBLIC_DATA_SOURCES.md",
      "docs/PREPRINT_ABSTRACT.md", "docs/site/index.html", "results/README.md"]),
    ("external_validation_finaledb.json",
     ("pooled_external_auc_matched_publications_only", "auc"),
     ["paper/METHODS_PAPER.md", "MODEL.md", "README.md", "README.zh-CN.md",
      "RESULTS.md", "CITATION.cff", "docs/ARXIV_ABSTRACT.md",
      "docs/EXTERNAL_VALIDATION.md",
      "docs/PREPRINT_ABSTRACT.md", "results/README.md"]),
    # The case-only publication is reported as a sensitivity at the
    # threshold fixed in training, so the threshold and that sensitivity
    # are quoted numbers too.
    ("external_validation_finaledb.json",
     ("training_derived_threshold", "threshold"),
     ["paper/METHODS_PAPER.md", "docs/EXTERNAL_VALIDATION.md"]),
    ("external_validation_finaledb.json",
     ("per_publication_auc", "9", "sensitivity_at_train_threshold_99spec"),
     ["paper/METHODS_PAPER.md", "docs/EXTERNAL_VALIDATION.md"]),
    # The naive pooled external AUC -- the figure that shows what
    # folding a case-only publication into an AUC denominator would
    # produce. It was corrupted once by a test-count string replacement
    # (0.7097 -> 0.7207 -> 0.7237), so it is quoted explicitly now.
    ("external_validation_finaledb.json",
     ("pooled_external_auc", "auc"),
     ["paper/METHODS_PAPER.md", "docs/EXTERNAL_VALIDATION.md"]),
    ("cross_study_finallydb_gc_corrected.json",
     ("pooled", "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md", "paper.md", "README.md",
      "README.zh-CN.md", "MODEL.md", "RESULTS.md", "TEAM.md", "REVIEWERS.md",
      "CITATION.cff", "results/README.md",
      "ULTRA_EARLY_READINESS.md", "USAGE.md",
      "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md",
      "docs/CROSS_PLATFORM_BENCHMARK.md", "docs/HARMONIZATION.md",
      "docs/site/index.html", "paper/README.md", "paper/REPRODUCE.md",
      # The EGA/dbGaP access request states the pooled AUC to justify the
      # request ("reports AUC 0.968 +/- 0.003"), so it is live in it.
      "data/EGA_DBGAP_DATA_ACCESS_TEMPLATE.md"]),
    # The per-cohort rows of the same table. They were quoted from an
    # earlier, unseeded run for a long time, so they are on the list too.
    ("cross_study_finallydb_gc_corrected.json",
     ("per_publication", "6", "auc_mean"),
     ["paper/METHODS_PAPER.md", "paper/REPRODUCE.md", "docs/HARMONIZATION.md",
      "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md"]),
    ("cross_study_finallydb_gc_corrected.json",
     ("per_publication", "8", "auc_mean"),
     ["paper/METHODS_PAPER.md", "paper/REPRODUCE.md", "docs/HARMONIZATION.md",
      "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md"]),
    # The true-confound control, which the manuscript quotes to argue the
    # pooled AUC is a cancer signal and not a study batch.
    ("cross_study_finallydb_gc_corrected.json",
     ("true_confound_control", "cancer_jiang_healthy_cristiano",
      "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md"]),
    # The legacy (no-GC) artifact carries the same control; the
    # harmonisation narrative quotes that one.
    ("cross_study_finallydb.json",
     ("true_confound_control", "cancer_jiang_healthy_cristiano",
      "harmonized", "auc_mean"),
     ["docs/HARMONIZATION.md"]),
    ("cross_study_finallydb.json",
     ("pooled", "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md", "README.md", "README.zh-CN.md",
      "ULTRA_EARLY_READINESS.md", "MODEL.md", "RESULTS.md", "TEAM.md",
      "results/README.md", "docs/CROSS_STUDY_BENCHMARK.md"]),
    ("methylation_cohort_channel.json",
     ("beta_summary_auc_loo",),
     ["paper/paper.tex", "paper.md", "README.md", "MODEL.md",
      "ULTRA_EARLY_READINESS.md", "docs/METHYLATION_PIPELINE.md",
      # These four documents quote the same 0.9444 but named
      # ``methylation_channel_diagnosis.json`` as their source -- and
      # that artifact does NOT hold the figure; it holds the MECHANISM
      # (the depth figures, the coverage baselines, the decoder
      # degeneracy). Naming the documents here is what makes the
      # misattribution fail the guard.
      "paper/METHODS_PAPER.md", "paper/README.md",
      "docs/PORTFOLIO_ROADMAP.md", "docs/METHYLATION_PIPELINE.md"]),
]
# The matched-locus test the same null rests on: 103,863 loci, a median
# cancer-minus-control difference of exactly 0.0, and a permutation
# p-value of 1.0. These live in the COHORT artifact, not the diagnosis
# one.
QUOTED_NUMBERS.append(
    ("methylation_cohort_channel.json",
     ("common_cpg_set", "per_cpg_test", "n_cpgs_tested"),
     ["paper/METHODS_PAPER.md", "paper/paper.tex", "README.md",
      "docs/METHYLATION_PIPELINE.md"]))

# The top-5 per-cancer rows of the §3.3 table. Generated rather than
# written out, so a new cancer type joins the guard automatically.
_TOP5_CANCERS = ("LUAD", "BRCA", "OV", "PAAD", "HCC_J")
for _cancer in _TOP5_CANCERS:
    QUOTED_NUMBERS.append((
        "cross_study_finallydb_gc_corrected.json",
        ("per_cancer", _cancer, "auc_mean"),
        # Every document that quotes the row is named, so a copy of one
        # cannot drift unnoticed.
        ["paper/METHODS_PAPER.md", "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md"]
        + (["docs/SENS_AT_SPEC_ABLATION.md"] if _cancer == "BRCA" else []),
    ))
    QUOTED_NUMBERS.append((
        "per_cancer_sens_at_spec.json",
        ("per_cancer", _cancer, "sens_at_99"),
        ["paper/METHODS_PAPER.md", "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md",
         "ULTRA_EARLY_READINESS.md", "docs/ULTRA_EARLY_READINESS_AUTO.md"],
    ))

# Numbers that need a lookup rather than a plain path.
# Documents that quote a pinned figure but were not named in the entry
# above. The guard checks that every named document quotes the figure, so
# a document that is not named is a copy that can drift silently. This
# audit found them by scanning every tracked document with the guard's
# own matcher; the historical records (the changelog, the audit reports,
# the deprecated v1.0 draft, the dated status note) are excluded, because
# naming them would require them to keep quoting figures forever.
EXTRA_DOCS_ENTRIES = [
    # The entries below are PAIRS, not dict keys. A dict literal
    # with a repeated key silently keeps only the LAST list, so
    # naming a document for a figure that already had an entry
    # dropped the earlier documents from the guard without any
    # error. Seven keys here were repeated, which had already
    # removed 31 document/figure protections -- including the
    # documents the superseded-figure guard is scoped to. The
    # pairs are merged by union below, so a repeated key adds to
    # the protection instead of replacing it.
    (("external_validation_finaledb.json", ('per_publication_auc', '7', 'auc')), [
        "NEXT_STEPS.md",
        "REVIEWERS.md",
        "ULTRA_EARLY_READINESS.md",
        "ZENODO_SETUP.md",
        "docs/METHYLATION_PIPELINE.md",
        "docs/PORTFOLIO_ROADMAP.md",
        "docs/PUBLIC_DATA_SOURCES.md",
        "paper/METHODS_PAPER.md",
        "paper/REPRODUCE.md",
    ]),
    (("external_validation_finaledb.json", ('per_publication_auc', '7', 'ci_lo')), [
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
    ]),
    (("external_validation_finaledb.json", ('pooled_external_auc_matched_publications_only', 'auc')), [
        "REVIEWERS.md",
        "TEAM.md",
        "docs/CADD_PER_SUBGROUP_LLR.md",
        "docs/CROSS_STUDY_BENCHMARK.md",
        # The public-data survey quotes the pooled figure when it says
        # what the ccRCC cohort would change about it (§6).
        "docs/PUBLIC_DATA_SOURCES.md",
        "docs/site/index.html",
        "paper/README.md",
        "paper/REPRODUCE.md",
        "paper/paper.tex",
    ]),
    (("external_validation_finaledb.json", ('training_derived_threshold', 'threshold')), [
        "RESULTS.md",
        "paper/REPRODUCE.md",
    ]),
    (("external_validation_finaledb.json", ('per_publication_auc', '9', 'sensitivity_at_train_threshold_99spec')), [
        "MODEL.md",
        "README.md",
        "README.zh-CN.md",
        "TEAM.md",
        "paper/REPRODUCE.md",
        "results/README.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('pooled', 'harmonized', 'auc_mean')), [
        "docs/CADD_PER_SUBGROUP_LLR.md",
        "docs/ULTRA_EARLY_READINESS_AUTO.md",
        "docs/USER_GUIDE.md",
        # The portfolio roadmap quotes the current headline in its
        # repository table, and it is a state document: if the number
        # moves, the roadmap has to move with it.
        "docs/PORTFOLIO_ROADMAP.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('per_publication', '6', 'auc_mean')), [
        "MODEL.md",
        "docs/CROSS_STUDY_BENCHMARK.md",
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('per_publication', '8', 'auc_mean')), [
        "docs/site/index.html",
        "paper/PAPER.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('true_confound_control', 'cancer_jiang_healthy_cristiano', 'harmonized', 'auc_mean')), [
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
        "paper/METHODS_PAPER.md",
    ]),
    (("cross_study_finallydb.json", ('true_confound_control', 'cancer_jiang_healthy_cristiano', 'harmonized', 'auc_mean')), [
        "docs/HARMONIZATION.md",
    ]),
    (("cross_study_finallydb.json", ('pooled', 'harmonized', 'auc_mean')), [
        "ZENODO_SETUP.md",
        "docs/CROSS_PLATFORM_BENCHMARK.md",
        "docs/HARMONIZATION.md",
        "docs/site/index.html",
        "paper/README.md",
        "paper/REPRODUCE.md",
    ]),
    (("methylation_cohort_channel.json", ('beta_summary_auc_loo',)), [
        "CITATION.cff",
        "NEXT_STEPS.md",
        "README.zh-CN.md",
        "RESULTS.md",
        "REVIEWERS.md",
        "TEAM.md",
        "ZENODO_SETUP.md",
        "docs/EXTERNAL_VALIDATION.md",
        "docs/JOURNAL_SUBMISSION.md",
        "docs/PREPRINT_ABSTRACT.md",
        "docs/PUBLIC_DATA_SOURCES.md",
        "docs/site/index.html",
        "paper/REPRODUCE.md",
        "results/README.md",
    ]),
    (("methylation_cohort_channel.json", ('common_cpg_set', 'per_cpg_test', 'n_cpgs_tested')), [
        "MODEL.md",
        "README.zh-CN.md",
        "RESULTS.md",
        "REVIEWERS.md",
        "TEAM.md",
        "ULTRA_EARLY_READINESS.md",
        "docs/PORTFOLIO_ROADMAP.md",
        "docs/site/index.html",
        "paper.md",
        "results/README.md",
    ]),
    (("per_cancer_sens_at_spec.json", ('per_cancer', 'OV', 'sens_at_99')), [
        "docs/USER_GUIDE.md",
    ]),
    (("per_cancer_sens_at_spec.json", ('per_cancer', 'PAAD', 'sens_at_99')), [
        "docs/CADD_GDC_VALIDATION.md",
    ]),
    (("per_cancer_sens_at_spec.json", ('per_cancer', 'HCC_J', 'sens_at_99')), [
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
    ]),
    (("cross_study_finallydb.json", ('per_publication', '6', 'auc_mean')), [
        "MODEL.md",
        "TEAM.md",
        "docs/CADD_ALPHAMISSENSE_COMBINED.md",
        "docs/CROSS_STUDY_LOSO.md",
        "paper/METHODS_PAPER.md",
        "results/README.md",
    ]),
    (("cross_study_finallydb.json", ('per_publication', '8', 'auc_mean')), [
        "docs/site/index.html",
        "paper/PAPER.md",
    ]),
    (("real_tcga_validation.json", ('panel_detection', 'panel_llr', 16, 'mean')), [
        "RESULTS.md",
        "docs/ALPHAMISSENSE_REAL_LLR.md",
        "docs/ALPHAMISSENSE_WEIGHTED_LLR.md",
        "docs/ARXIV_ABSTRACT.md",
        "docs/CADD_PER_SUBGROUP_LLR.md",
        "docs/CADD_WEIGHTED_LLR.md",
        "docs/PORTFOLIO_ROADMAP.md",
        "docs/PREPRINT_ABSTRACT.md",
        "paper/README.md",
        "paper/REPRODUCE.md",
    ]),
    (("real_tcga_validation.json", ('panel_detection', 'panel_llr', 18, 'mean')), [
        "RESULTS.md",
        "docs/CADD_PER_SUBGROUP_LLR.md",
        "docs/CADD_WEIGHTED_LLR.md",
    ]),
    (("fusion_regularized_smoke.json", ('summary', 'regularized_late_fusion', 'auc_mean')), [
        "RESULTS.md",
        "ZENODO_SETUP.md",
        "docs/CADD_PER_SUBGROUP_LLR.md",
        "docs/CROSS_STUDY_BENCHMARK.md",
        "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md",
        "paper/REPRODUCE.md",
    ]),
    (("fusion_regularized_smoke.json", ('summary', 'bayesian_model_average', 'auc_mean')), [
        "docs/CROSS_STUDY_BENCHMARK.md",
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
    ]),
    (("cross_study_finallydb_with_motifs.json", ('pooled', 'harmonized', 'auc_mean')), [
        "docs/ALPHAMISSENSE_REAL_LLR.md",
    ]),
    (("cross_study_finallydb_with_motifs.json", ('per_publication', '6', 'auc_mean')), [
        "RESULTS.md",
        "USAGE.md",
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
        "docs/PORTFOLIO_ROADMAP.md",
    ]),
    (("cross_study_finallydb_with_motifs.json", ('per_publication', '8', 'auc_mean')), [
        "USAGE.md",
        "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
    ]),
    (("decision_curve.json", ('per_strategy', 'tumor_naive', 'auc')), [
        "paper/METHODS_PAPER.md",
        "results/README.md",
    ]),
    # Three more documents that quote the headline figures and were not
    # named: the methods-paper status file, the wiki's Jiang page and the
    # V3 design note. They quote the numbers, so a regenerated artifact
    # could have moved them without anything noticing.
    #
    # Each pin was read in its document before it was added: the status
    # file restates the methods paper's §3.2 rows (Jiang 0.9789,
    # Cristiano 0.9704, pooled 0.9770, and the GC-corrected 0.9661 /
    # 0.9698 / 0.9684 with sens@99 0.700), the wiki's LC-vs-Control row
    # is the Jiang publication's per-publication AUC, and the V3 design
    # note quotes the panel LLR anchor (0.921). A document that merely
    # rounds to the same number -- a range's endpoint, or a threshold
    # written as a target -- is NOT named here: the quote test matches
    # at three decimals, so such a document would pass the check without
    # ever quoting the figure.
    (("cross_study_finallydb.json", ('per_publication', '6', 'auc_mean')), [
        "paper/METHODS_PAPER_STATUS.md",
        "wiki/Jiang-4mer-Validation.md",
    ]),
    (("cross_study_finallydb.json", ('per_publication', '8', 'auc_mean')), [
        "paper/METHODS_PAPER_STATUS.md",
    ]),
    (("cross_study_finallydb.json", ('pooled', 'harmonized', 'auc_mean')), [
        "paper/METHODS_PAPER_STATUS.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('per_publication', '6', 'auc_mean')), [
        "paper/METHODS_PAPER_STATUS.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('per_publication', '8', 'auc_mean')), [
        "paper/METHODS_PAPER_STATUS.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('pooled', 'harmonized', 'auc_mean')), [
        "paper/METHODS_PAPER_STATUS.md",
    ]),
    (("cross_study_finallydb_gc_corrected.json", ('pooled', 'harmonized', 'sens_at_spec_99')), [
        "paper/METHODS_PAPER_STATUS.md",
    ]),
    (("real_tcga_validation.json", ('panel_detection', 'panel_llr', 16, 'mean')), [
        "docs/V3_DESIGN.md",
    ]),
]

# Merge the pairs by union. ``EXTRA_DOCS`` keeps the mapping the
# rest of the module (and the guard) reads; a repeated key now
# ADDS its documents to the protection instead of replacing it.
EXTRA_DOCS: dict[tuple, list] = {}
for _extra_key, _extra_docs in EXTRA_DOCS_ENTRIES:
    _known = EXTRA_DOCS.setdefault(_extra_key, [])
    _known.extend(_doc for _doc in _extra_docs if _doc not in _known)

# Documents that RECORD the project's history rather than state its
# current results: the changelog and the dated audit reports. They name
# superseded figures on purpose -- describing the old value is the point
# of those files -- so they are not held to the current artifacts. The
# changelog has no "unreleased" section: every entry in it is a released,
# dated record.
HISTORICAL_RECORDS = {
    "CHANGELOG.md",
    "AUDIT_REPORT.md",
    "AUDIT_REPORT_2.md",
    "AUDIT_2_FINDINGS.md",
}

# Documents whose text HAPPENS to contain the digits of a pinned figure
# without quoting that figure: a range that spans it, a target written in
# a design table, a dated review record, or a number that belongs to a
# different experiment. Each one is listed with the reason it is not a
# copy, and the guard checks the list in both directions -- an excuse
# that stops being needed fails, and an excuse with no reason fails, so
# the list cannot quietly become a place where drift is allowed.
#
# The figure is named by (artifact, path) so a document can be excused
# for one figure and stay protected for another.
EXCUSED_COPIES = {
    ("external_validation_finaledb.json",
     ("per_publication_auc", "7", "ci_lo")): {
        "docs/V3_DESIGN.md":
            "a design-table TARGET ('now >= 0.873'), not a quote of the "
            "held-out study's DeLong lower bound",
    },
    ("cross_study_finallydb_gc_corrected.json",
     ("per_cancer", "OV", "auc_mean")): {
        "results/prof_jiang_4mer_analysis/summary_for_professor_jiang.md":
            "HCC-vs-control AUCs from the 4-mer analysis (0.9845 / "
            "0.9885 / 0.9886); the match with the ovarian per-cancer AUC "
            "is a coincidence, not a quote of it",
        "docs/ALPHAGENOME_PROXY_RESULTS.md":
            "the AlphaGenome proxy methods' own AUCs (all eight at >= "
            "0.989), not the ovarian per-cancer AUC",
        "docs/CADD_PER_SUBGROUP_LLR.md":
            "the Top-K=200 AUC of the STK11-mutant subgroup (0.989 +/- "
            "0.015), not the ovarian per-cancer AUC",
    },
    ("cross_study_finallydb_gc_corrected.json",
     ("per_cancer", "PAAD", "auc_mean")): {
        "review/agent_review_2026-08-10.md":
            "a dated review record quoting a sensitivity table (0.969 / "
            "0.967), not the pancreatic per-cancer AUC",
    },
    ("fusion_regularized_smoke.json",
     ("summary", "regularized_late_fusion", "auc_mean")): {
        "review/agent_review_2026-08-10.md":
            "a dated review record listing the OLD README's 0.982, not "
            "the late-fusion AUC",
    },
    ("cross_study_finallydb_with_motifs.json",
     ("per_publication", "6", "auc_mean")): {
        "METHYLATION_PROJECT.md":
            "the companion repository's cross-study AUC range "
            "(0.974-0.978), not the with-motifs Jiang row",
        "docs/V3_DESIGN.md":
            "the v2.2 design baseline ('already at AUC 0.978'), not the "
            "with-motifs Jiang row",
    },
    ("cross_study_finallydb_with_motifs.json",
     ("per_publication", "8", "auc_mean")): {
        "METHYLATION_PROJECT.md":
            "the companion repository's cross-study AUC range "
            "(0.974-0.978), not the with-motifs Cristiano row",
    },
    ("decision_curve.json",
     "decision curve tumor_naive @0.95 threshold"): {
        "results/prof_jiang_4mer_analysis/summary_for_professor_jiang.md":
            "the LC-vs-Control sensitivity 0.700 in the 4-mer analysis "
            "table, not the decision curve's operating threshold",
        "paper/METHODS_PAPER_STATUS.md":
            "the pooled sensitivity at 99% specificity (sens@99 = 0.700), "
            "not the decision curve's operating threshold",
    },
    ("cross_study_loso.json", "in_domain cristiano"): {
        "review/agent_review_2026-08-10.md":
            "a dated review record quoting a variant-caller AUC table "
            "(0.969 / 0.967), not the within-publication reference AUC",
    },
}

QUOTED_SPECIAL = [
    ("cross_study_loso.json", "loso jiang -> cristiano",
     "paper/paper.tex"),
    ("cross_study_loso.json", "loso cristiano -> jiang",
     "paper/paper.tex"),
    # CITATION.cff is what GitHub's "Cite this repository" button hands a
    # reader, and it quotes the jiang -> cristiano transfer direction.
    ("cross_study_loso.json", "loso jiang -> cristiano",
     "CITATION.cff"),
    # The model card and the results index quote both transfer directions.
    ("cross_study_loso.json", "loso jiang -> cristiano",
     "MODEL.md"),
    ("cross_study_loso.json", "loso cristiano -> jiang",
     "MODEL.md"),
    ("cross_study_loso.json", "loso jiang -> cristiano",
     "RESULTS.md"),
    ("cross_study_loso.json", "loso cristiano -> jiang",
     "RESULTS.md"),
    ("fusion_ablation_isotonic.json", "naive_average",
     "paper/paper.tex"),
]

# Every document that quotes a transfer direction or an in-domain
# reference is named, so a copy of one cannot drift unnoticed. The
# documents are the ones that actually quote the figure today; the
# guard fails if one of them stops quoting it.
_LOSO_DOCS = [
    "NEXT_STEPS.md", "README.md", "README.zh-CN.md", "REVIEWERS.md",
    "TEAM.md", "docs/ARXIV_ABSTRACT.md", "docs/CROSS_STUDY_BENCHMARK.md",
    "docs/CROSS_STUDY_LOSO.md",
    "docs/JOURNAL_SUBMISSION.md", "docs/PREPRINT_ABSTRACT.md",
    "docs/site/index.html", "paper.md",
    "paper/METHODS_PAPER.md", "paper/README.md", "paper/REPRODUCE.md",
    "results/README.md",
]
for _doc in _LOSO_DOCS:
    QUOTED_SPECIAL.append(
        ("cross_study_loso.json", "loso jiang -> cristiano", _doc))
_LOSO_DOCS_CRISTIANO = [
    "NEXT_STEPS.md", "REVIEWERS.md", "TEAM.md", "docs/CROSS_STUDY_LOSO.md",
    "docs/site/index.html", "paper/METHODS_PAPER.md", "paper/README.md",
    "paper/REPRODUCE.md", "results/README.md",
]
for _doc in _LOSO_DOCS_CRISTIANO:
    QUOTED_SPECIAL.append(
        ("cross_study_loso.json", "loso cristiano -> jiang", _doc))
for _doc in ("docs/CROSS_STUDY_BENCHMARK.md", "paper/REPRODUCE.md"):
    QUOTED_SPECIAL.append(("cross_study_loso.json", "in_domain jiang", _doc))


# Figures that the v2.3.0 regeneration superseded. They must not reappear
# in a current-facing document; the CHANGELOG and the dated audit records
# are allowed to name them, because describing the old value is the
# point of those files.
#
# The check is a literal one, so a value that merely shares a prefix
# with a current figure does not trip it: the old GC-corrected pooled
# AUC "0.9670" is banned, but the per-cancer PAAD AUC "0.967" and the
# clean-panel LLR AUC "0.922" are different numbers and stay.
# Figures that the v2.3.0 regeneration superseded.
#
# A global literal ban is only safe while no CURRENT artifact renders the
# same digits -- and several of these do:
#   0.9915 is the decision curve's naive-average AUC (0.991505968778696)
#   0.9670 is the OV per-cancer DeLong CI upper bound
#   0.9701 is one of the pooled per-seed AUCs
#   0.9753 is the OTHER_C per-cancer AUC in the with-motifs artifact
#   0.499  is a decision-curve net-benefit value
#   0.9225, 0.9458, 0.9747 are rendered by other artifacts too
# Banning those everywhere would forbid a correct document, so each is
# banned only where it would be a regression: the documents that quote
# the quantity it replaced. The four values below are rendered by no
# current artifact and are banned everywhere.
#
# The check is a literal one, so a value that merely shares a prefix
# with a current figure does not trip it: the old GC-corrected pooled
# AUC "0.9670" is banned, but the per-cancer PAAD AUC "0.967" and the
# clean-panel LLR AUC "0.922" are different numbers and stay.
SUPERSEDED = [
    "0.4657",   # the old shuffled-label null
    "0.8819",   # the old pooled external AUC
    "0.4454",   # the old publication-9 external sensitivity
    # The old mutation-only fusion AUC, from the pre-regeneration
    # fusion artifact. It is the only one of the stale fusion figures
    # that can be banned literally: the others ("0.9887", "0.937",
    # "0.845") are prefixes of legitimate current numbers, so they are
    # pinned by the table-row checks below instead.
    "0.9242",
]

# Superseded values that a current artifact also renders, and therefore
# can only be banned where they would be a regression:
# value -> (artifact, path/what) whose quoting documents may not carry it.
SCOPED_SUPERSEDED = {
    # The old publication-7 external AUC. Replaced by the corrected
    # hold-out AUC.
    "0.9915": ("external_validation_finaledb.json",
               ("per_publication_auc", "7", "auc")),
    # The old GC-corrected pooled AUC.
    "0.9670": ("cross_study_finallydb_gc_corrected.json",
               ("pooled", "harmonized", "auc_mean")),
    # The old in-domain Cristiano AUC (the seeded LOSO run replaced it).
    "0.9701": ("cross_study_loso.json", "in_domain cristiano"),
    # The old legacy pooled AUC.
    "0.9747": ("cross_study_finallydb.json",
               ("pooled", "harmonized", "auc_mean")),
    # The companion pipeline's unseeded tumor-naive figure.
    "0.9753": ("cross_study_finallydb.json",
               ("pooled", "harmonized", "auc_mean")),
    # The old PCA-20 pooled AUC.
    "0.9225": ("cross_study_finallydb_gc_corrected.json",
               ("pooled", "harmonized", "auc_mean")),
    # The old PCA-20 per-publication AUC.
    "0.9458": ("cross_study_finallydb_gc_corrected.json",
               ("per_publication", "6", "auc_mean")),
    # The with-motifs pooled AUC of the LEGACY-footed run, superseded when
    # the GC correction was actually applied to that path. Scoped: 0.9785
    # is also the CADD Top-K=200 AUC in docs/ALPHAMISSENSE_REAL_LLR.md.
    "0.9785": ("cross_study_finallydb_with_motifs.json",
               ("pooled", "harmonized", "auc_mean")),
    # The readiness generator's hard-coded true-confound value, superseded
    # by the GC-corrected control.
    "0.499": ("cross_study_finallydb_gc_corrected.json",
              ("true_confound_control", "cancer_jiang_healthy_cristiano",
               "harmonized", "auc_mean")),
}

# Documents where a scoped value is legitimate for another reason: the
# decision curve's naive-average AUC renders as 0.9915, so that document
# may quote it; and the AlphaMissense document quotes 0.9785 as the
# CADD Top-K=200 AUC while also quoting the with-motifs pooled AUC, so
# it may quote that too.
SCOPED_ALLOWED = {
    "0.9915": {"results/README.md"},
    "0.9785": {"docs/ALPHAMISSENSE_REAL_LLR.md"},
}

# The special (label-addressed) figures have the same problem: a document
# that quotes one was not named for it. These rows are the missing names.
QUOTED_SPECIAL_EXTRA = [
    ("fusion_ablation_isotonic.json", "naive_average", "MODEL.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "README.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "README.zh-CN.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "RESULTS.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "TEAM.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "USAGE.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "docs/CADD_PER_SUBGROUP_LLR.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "docs/CROSS_STUDY_BENCHMARK_with_motifs.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "docs/FUSION_ISOTONIC.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "docs/ARXIV_ABSTRACT.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "docs/PREPRINT_ABSTRACT.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "docs/site/index.html"),
    ("fusion_ablation_isotonic.json", "naive_average", "paper/PAPER.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "paper/README.md"),
    ("fusion_ablation_isotonic.json", "naive_average", "results/README.md"),
    # The project note quotes the same fusion headline (it used to carry
    # the older run's 0.989, which no artifact holds).
    ("fusion_ablation_isotonic.json", "naive_average", "METHYLATION_PROJECT.md"),
    ("cross_study_loso.json", "in_domain jiang", "docs/CADD_ALPHAMISSENSE_COMBINED.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "MODEL.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "README.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "README.zh-CN.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "ULTRA_EARLY_READINESS.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "docs/CADD_PER_SUBGROUP_LLR.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "docs/ULTRA_EARLY_READINESS_AUTO.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "paper/METHODS_PAPER.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.95 threshold", "paper/REPRODUCE.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.99 threshold", "docs/EXTERNAL_VALIDATION.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.99 threshold", "paper/METHODS_PAPER.md"),
    ("decision_curve.json", "decision curve tumor_naive @0.99 threshold", "paper/REPRODUCE.md"),
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_95", "MODEL.md"),
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_95", "docs/FUSION_ISOTONIC.md"),
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_99", "MODEL.md"),
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_99", "docs/CROSS_STUDY_BENCHMARK.md"),
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_99", "docs/FUSION_ISOTONIC.md"),
]


def documents_for(artifact: str, path) -> set[str]:
    """The documents that quote a given (artifact, path) or (artifact, what)."""
    docs: set[str] = set()
    for name, quoted_path, documents in QUOTED_NUMBERS:
        if name == artifact and tuple(quoted_path) == tuple(path):
            docs.update(documents)
    for (name, quoted_path), extra in EXTRA_DOCS.items():
        if name == artifact and tuple(quoted_path) == tuple(path):
            docs.update(extra)
    for name, what, document in QUOTED_SPECIAL:
        if name == artifact and what == path:
            docs.add(document)
    for name, what, document in QUOTED_SPECIAL_EXTRA:
        if name == artifact and what == path:
            docs.add(document)
    return docs


def scoped_superseded_documents() -> dict[str, list[str]]:
    """{superseded value: the documents it may not reappear in}."""
    out: dict[str, list[str]] = {}
    for value, (artifact, target) in SCOPED_SUPERSEDED.items():
        docs = documents_for(artifact, target)
        docs -= SCOPED_ALLOWED.get(value, set())
        out[value] = sorted(docs)
    return out


# The legacy (no-GC) per-publication rows. They are a separate benchmark
# from the GC-corrected one and are quoted by the legacy narrative, so
# they are guarded separately -- otherwise a document could satisfy the
# guard with the GC row while quoting the legacy one.
QUOTED_NUMBERS.append(
    ("cross_study_finallydb.json",
     ("per_publication", "6", "auc_mean"),
     ["docs/CROSS_STUDY_BENCHMARK.md", "paper/REPRODUCE.md"]))
QUOTED_NUMBERS.append(
    ("cross_study_finallydb.json",
     ("per_publication", "8", "auc_mean"),
     ["docs/CROSS_STUDY_BENCHMARK.md", "paper/REPRODUCE.md"]))

# The in-domain references of the leave-one-study-out benchmark. They
# are the within-publication CV the transfer numbers are compared
# against, and they live in their own artifact.
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain jiang",
     "paper/METHODS_PAPER.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain cristiano",
     "paper/METHODS_PAPER.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain jiang",
     "docs/CROSS_STUDY_LOSO.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain cristiano",
     "docs/CROSS_STUDY_LOSO.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain jiang",
     "MODEL.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain cristiano",
     "MODEL.md"))
# The team page and the results index quote the same in-domain pair
# against the transfer numbers.
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain jiang",
     "TEAM.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain cristiano",
     "TEAM.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain jiang",
     "results/README.md"))
QUOTED_SPECIAL.append(
    ("cross_study_loso.json", "in_domain cristiano",
     "results/README.md"))

# Documents that quote the true-confound control. The GC-corrected
# value is the one the headline benchmark uses; the legacy value is
# quoted by the legacy narrative. Both are guarded, because a document
# that quotes the wrong one is mixing artifacts.
QUOTED_NUMBERS.append(
    ("cross_study_finallydb_gc_corrected.json",
     ("true_confound_control", "cancer_jiang_healthy_cristiano",
      "harmonized", "auc_mean"),
     ["README.md", "README.zh-CN.md", "RESULTS.md", "ULTRA_EARLY_READINESS.md",
      "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md",
      "docs/ULTRA_EARLY_READINESS_AUTO.md"]))
QUOTED_NUMBERS.append(
    ("cross_study_finallydb.json",
     ("true_confound_control", "cancer_jiang_healthy_cristiano",
      "harmonized", "auc_mean"),
     ["docs/CROSS_STUDY_BENCHMARK.md"]))

# The panel-LLR headline of the mutation-informed pipeline -- the first
# figure in the preprint's abstract. The driver regenerates the artifact
# (step 11), so it belongs in the quotable-number list like the rest.
QUOTED_NUMBERS.append(
    ("real_tcga_validation.json", ("panel_detection", "panel_llr", 16, "mean"),
     ["paper/paper.tex", "paper/PAPER.md", "paper.md", "MODEL.md",
      "REVIEWERS.md", "README.md", "README.zh-CN.md", "NEXT_STEPS.md",
      "ULTRA_EARLY_READINESS.md", "results/README.md",
      "docs/site/index.html"]))

# The panel's sensitivity at 99% specificity, quoted next to its AUC.
# It was unguarded, and one document carried it as "0.49" -- a two-decimal
# rendering the guard cannot see, because a number written with fewer than
# three decimals cannot identify a value. Written as 0.460 and named here,
# a later rewrite of that document fails the guard instead of drifting.
QUOTED_NUMBERS.append(
    ("real_tcga_validation.json", ("panel_detection", "panel_llr", 18, "mean"),
     ["paper/paper.tex", "paper/PAPER.md", "MODEL.md",
      "ULTRA_EARLY_READINESS.md", "docs/PREPRINT_ABSTRACT.md"]))

# The goal-2 headline: the regularized late fusion and the Bayesian
# model average that replaced the overfitting neural encoder. It was
# quoted by ten documents and pinned by none, and the driver did not
# regenerate it -- so the central claim of the release was both
# unguarded and unreproducible.
QUOTED_NUMBERS.append(
    ("fusion_regularized_smoke.json",
     ("summary", "regularized_late_fusion", "auc_mean"),
     ["MODEL.md", "README.md", "README.zh-CN.md", "ULTRA_EARLY_READINESS.md",
      "docs/ULTRA_EARLY_READINESS_AUTO.md", "paper/METHODS_PAPER.md",
      "paper/README.md", "CITATION.cff",
      "results/README.md", "TEAM.md", "REVIEWERS.md", "NEXT_STEPS.md"]))
QUOTED_NUMBERS.append(
    ("fusion_regularized_smoke.json",
     ("summary", "bayesian_model_average", "auc_mean"),
     # paper/METHODS_PAPER_STATUS.md was removed from this list: it quoted
     # 0.9650 only as a superseded legacy no-harmonize figure, so the
     # requirement was satisfied by a coincidence and would have blocked
     # the correction of that stale number.
     ["MODEL.md", "README.md", "README.zh-CN.md", "ULTRA_EARLY_READINESS.md",
      "paper/METHODS_PAPER.md", "paper/REPRODUCE.md", "CITATION.cff"]))
# The POOLED measured null on the paired cohort -- the number the
# documents read the control against, because on this design the null
# sits below 0.5 rather than at it.
#
# The documents used to say "0.42 +/- 0.12", which is the HIGHEST
# seed's mean and std rounded to two decimals; the cohort's pooled null
# is 0.4063 +/- 0.1143. The difference is under a rounding unit, so
# nothing flagged it -- and at two decimals the number guard cannot
# see it at all, because it requires at least three decimals to
# identify a value. The artifact now carries the pooled figure at
# `protocol.measured_null`, and these two pins hold the documents to
# it at the precision they write it.
QUOTED_NUMBERS.append(
    ("fusion_regularized_smoke.json",
     ("protocol", "measured_null", "mean"),
     ["MODEL.md", "README.md", "README.zh-CN.md", "ULTRA_EARLY_READINESS.md",
      "paper/METHODS_PAPER.md", "CHANGELOG.md"]))
QUOTED_NUMBERS.append(
    ("fusion_regularized_smoke.json",
     ("protocol", "measured_null", "std"),
     ["MODEL.md", "README.md", "README.zh-CN.md", "ULTRA_EARLY_READINESS.md",
      "paper/METHODS_PAPER.md", "CHANGELOG.md"]))
# The shuffled-label null. Six documents quote 0.4460 as "the signal is
# not in the labels"; nothing pinned it.
QUOTED_NUMBERS.append(
    ("cross_study_finallydb_shuffled_control.json",
     ("pooled_auc_mean",),
     ["README.md", "README.zh-CN.md", "ULTRA_EARLY_READINESS.md",
      "paper/METHODS_PAPER.md",
      # The EGA/dbGaP access request states the headline AUC and the
      # shuffled-label null side by side ("AUC 0.968 +/- 0.003 with
      # shuffled-label null 0.446") to justify the request, so both
      # figures are live in it.
      "data/EGA_DBGAP_DATA_ACCESS_TEMPLATE.md"]))
# The with-motifs benchmark: a third tumor-naive run, quoted by eight
# documents as "the same pipeline with 4-mer end motifs".
QUOTED_NUMBERS.append(
    ("cross_study_finallydb_with_motifs.json",
     ("pooled", "harmonized", "auc_mean"),
     ["README.md", "README.zh-CN.md", "RESULTS.md", "TEAM.md",
      "ULTRA_EARLY_READINESS.md", "paper/METHODS_PAPER.md",
      "results/README.md", "docs/CROSS_STUDY_BENCHMARK_with_motifs.md"]))
# The per-publication rows of the same run, which the methods paper
# quotes as the "both publications improve" claim.
for _pub in ("6", "8"):
    QUOTED_NUMBERS.append(
        ("cross_study_finallydb_with_motifs.json",
         ("per_publication", _pub, "auc_mean"),
         ["paper/METHODS_PAPER.md"]))

# The decision-curve operating points of RESULTS.md §4. They are a
# separate benchmark from the fusion one -- the decision curve averages
# each sample's out-of-fold score across seeds, the fusion benchmark
# averages the per-seed metric -- so they are guarded against their own
# artifact, not against the fusion one.
QUOTED_NUMBERS.append(
    ("decision_curve.json",
     ("per_strategy", "tumor_naive", "auc"),
     ["RESULTS.md"]))
for _strategy in ("tumor_naive", "naive_average"):
    for _spec in (0.95, 0.99):
        QUOTED_SPECIAL.append(
            ("decision_curve.json",
             f"decision curve {_strategy} @{_spec} threshold",
             "RESULTS.md"))
# The lr_fusion operating points in the same table come from the fusion
# benchmark, because the decision-curve CLI does not compute that
# strategy.
QUOTED_SPECIAL.append(
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_95",
     "RESULTS.md"))
QUOTED_SPECIAL.append(
    ("fusion_ablation_isotonic.json", "fusion lr_fusion sens_at_99",
     "RESULTS.md"))

# ─────────────────────────────────────────────────────────────────────
# Table rows that must be internally consistent
# ─────────────────────────────────────────────────────────────────────
# A document can satisfy the document-level checks above and still be
# wrong: if a table row quotes numbers from two different runs, every
# number in it exists somewhere in the repo and the guard stays green.
# That is exactly what happened in the methods paper's ablation table,
# where the harmonization row mixed the legacy pooled AUC with the
# GC-corrected true-confound value and pointed at the legacy artifact.
#
# Each case below pins one table row: the row that contains ``anchor``
# must quote every required number, must NOT quote any forbidden one
# (the same quantity from a different run), and -- when ``source`` is
# given -- must name that artifact in its source column.
TABLE_ROWS = [
    # The GC-corrected pooled row: the AUC *and* its standard
    # deviation, with the legacy run's pair forbidden. The AUC was
    # pinned, but the std was not -- and the std is the number a
    # reviewer checks against the "+/-" in the abstract. Because the two
    # runs' stds render differently (0.0027 vs 0.0013), a row that mixed
    # them would now fail.
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| **pooled (harmonized, GC corrected)** |",
        "required": [
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "auc_std")),
        ],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "auc_std")),
        ],
        "source": None,
    },
    # The methods paper's section 3.1 table: the paired synthetic-plasma
    # honest baseline. Its figures were an audit-2 record that no
    # committed artifact held, so the artifact is now regenerated by the
    # documented command and each row is pinned to it.
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `panel_only AUC` |",
        "required": [("sens_at_spec_ce_n20.json", ("panel_only_auc_mean",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `frag_only AUC` (real mutation-derived) |",
        "required": [("sens_at_spec_ce_n20.json", ("frag_only_auc_mean",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `foundation AUC` |",
        "required": [("sens_at_spec_ce_n20.json", ("foundation_auc_mean",)),
                     ("sens_at_spec_ce_n20.json", ("foundation_auc_std",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `lr_baseline AUC` (sklearn LR on [panel, frag]) |",
        "required": [("sens_at_spec_ce_n20.json", ("lr_baseline_auc_mean",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `naive_avg AUC` |",
        "required": [("sens_at_spec_ce_n20.json", ("naive_avg_auc_mean",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `shuffled_foundation AUC` (permuted-label null) |",
        "required": [("sens_at_spec_ce_n20.json",
                       ("shuffled_foundation_auc_mean",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        # The raw gap between the honest AUC and its permuted-label
        # null. It is pinned separately from ``delta_auc_normalized``
        # because the normalised form divides by ``(real - 0.5)``,
        # which assumes a null at 0.5 -- and on this paired design the
        # null sits at ~0.35-0.42, below chance. The gap is the
        # interpretable quantity.
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `foundation_auc_minus_null` (honest − null) |",
        "required": [("sens_at_spec_ce_n20.json",
                       ("foundation_auc_minus_null",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| `delta_auc_normalized` |",
        "required": [("sens_at_spec_ce_n20.json", ("delta_auc_normalized",))],
        "forbidden": [],
        "source": "results/sens_at_spec_ce_n20.json",
    },
    # The readiness document's panel-sweep table. The document-level check
    # only asks whether 0.460 appears SOMEWHERE in it -- and it did, in
    # another table -- while this table carried 0.490 for the same
    # quantity. A row pin asks the ROW to quote the artifact's value.
    {
        "doc": "ULTRA_EARLY_READINESS.md",
        "anchor": "**Sens @ 99% spec**",
        "required": [("real_tcga_validation.json",
                       ("panel_detection", "panel_llr", 18, "mean"))],
        "forbidden": [],
    },
    {
        "doc": "ULTRA_EARLY_READINESS.md",
        "anchor": "**Panel AUC**",
        "required": [("real_tcga_validation.json",
                       ("panel_detection", "panel_llr", 16, "mean")),
                     ("real_tcga_validation.json",
                       ("panel_detection", "panel_llr", 16, "std"))],
        "forbidden": [],
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "Can the features separate cancer from healthy",
        "required": [("cross_study_finallydb_gc_corrected.json",
                       ("pooled", "harmonized", "auc_mean"))],
        "forbidden": [("cross_study_finallydb.json",
                       ("pooled", "harmonized", "auc_mean"))],
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "Is that signal cancer rather than study identity?",
        "required": [("cross_study_finallydb_gc_corrected.json",
                       ("true_confound_control",
                        "cancer_jiang_healthy_cristiano",
                        "harmonized", "auc_mean"))],
        "forbidden": [("cross_study_finallydb.json",
                       ("true_confound_control",
                        "cancer_jiang_healthy_cristiano",
                        "harmonized", "auc_mean"))],
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "Does the classifier transfer to an unseen publication?",
        "required": [
            ("external_validation_finaledb.json",
             ("per_publication_auc", "7", "auc")),
            ("external_validation_finaledb.json",
             ("pooled_external_auc_matched_publications_only", "auc")),
            ("external_validation_finaledb.json",
             ("per_publication_auc", "9",
              "sensitivity_at_train_threshold_99spec")),
        ],
        "forbidden": [],
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "Does it transfer to an unseen publication in the "
                  "single-type direction?",
        "required": [("cross_study_loso.json", "loso jiang -> cristiano")],
        "forbidden": [],
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "**harmonization** (per-cohort batch-effect check)",
        "required": [
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "no_harmonize", "auc_mean")),
        ],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "harmonized", "auc_mean")),
        ],
        "source": "cross_study_finallydb_gc_corrected.json",
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| cancer=Jiang, healthy=Cristiano |",
        "required": [
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "no_harmonize", "auc_mean")),
        ],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "harmonized", "auc_mean")),
        ],
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| cancer=Cristiano, healthy=Jiang |",
        "required": [
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_cristiano_healthy_jiang",
              "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_cristiano_healthy_jiang",
              "no_harmonize", "auc_mean")),
        ],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("true_confound_control", "cancer_cristiano_healthy_jiang",
              "harmonized", "auc_mean")),
        ],
    },
    # The fusion table. Every row comes from the same artifact, so each
    # row is pinned individually: the table had drifted to a mix of two
    # runs (a stale tumor-naive AUC, a stale mutation-only AUC and a
    # stale LR-fusion AUC) while every one of those numbers still
    # existed somewhere in the repository.
    {
        "doc": "MODEL.md",
        "anchor": "| Tumor-naive only |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "tumor_naive", "auc_mean"))],
        "forbidden": [],
    },
    # The model card's tumor-naive table. It used to carry the sibling
    # repository's unseeded numbers (0.9716 and 0.9746 -- the second is
    # one seed's value out of the legacy artifact's five, not the pooled
    # mean), which no committed artifact held. Each row is now pinned to
    # the artifact it comes from, and the two cross-study rows are kept
    # apart because they are different benchmarks.
    {
        "doc": "MODEL.md",
        "anchor": "| Single-study Jiang 2015 (121 samples) |",
        "required": [("cross_study_finallydb_gc_corrected.json",
                       ("per_publication", "6", "auc_mean"))],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("per_publication", "6", "auc_mean")),
        ],
    },
    {
        "doc": "MODEL.md",
        "anchor": "| Cross-study pan-cancer (627 samples, harmonized) |",
        "required": [
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "sens_at_spec_95")),
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "sens_at_spec_99")),
        ],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "auc_mean")),
        ],
    },
    {
        "doc": "MODEL.md",
        "anchor": "| Cross-study pan-cancer, legacy no-GC correction |",
        "required": [
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "sens_at_spec_95")),
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "sens_at_spec_99")),
        ],
        "forbidden": [
            ("cross_study_finallydb_gc_corrected.json",
             ("pooled", "harmonized", "auc_mean")),
        ],
    },
    {
        "doc": "MODEL.md",
        "anchor": "| Mutation-only (synthetic, AUC 0.92) |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "mutation_only", "auc_mean"))],
        "forbidden": [],
    },
    {
        "doc": "MODEL.md",
        "anchor": "| **Naive average** |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "naive_average", "auc_mean"))],
        "forbidden": [],
    },
    {
        "doc": "MODEL.md",
        "anchor": "| **LR fusion** |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "lr_fusion", "auc_mean"))],
        "forbidden": [],
    },
    {
        "doc": "docs/FUSION_ISOTONIC.md",
        "anchor": "| Tumor-naive only |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "tumor_naive", "auc_mean"))],
        "forbidden": [],
    },
    {
        "doc": "docs/FUSION_ISOTONIC.md",
        "anchor": "| Mutation-only (calibrated AUC 0.92) |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "mutation_only", "auc_mean"))],
        "forbidden": [],
    },
    {
        "doc": "docs/FUSION_ISOTONIC.md",
        "anchor": "| Naive average | **",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "naive_average", "auc_mean"))],
        "forbidden": [],
    },
    {
        "doc": "docs/FUSION_ISOTONIC.md",
        "anchor": "| LR fusion (uncalibrated) |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "lr_fusion", "auc_mean"))],
        "forbidden": [],
    },
    # The Simplified-Chinese README quotes the same figures as the
    # English one; the translation-drift check in the test file
    # requires it to carry every figure the original states.
    {
        "doc": "README.zh-CN.md",
        "anchor": "| **+ 4-mer 末端基序**（`--include-motifs`） |",
        "required": [("cross_study_finallydb_with_motifs.json",
                       ("pooled", "harmonized", "auc_mean"))],
        "forbidden": [],
    },
    # The readiness document's per-cancer sens@99 row, which had been
    # left quoting a pre-regeneration run of the per-cancer table.
    {
        "doc": "ULTRA_EARLY_READINESS.md",
        "anchor": "Per-cancer sens@99% spec on the same cohort",
        "required": [("per_cancer_sens_at_spec.json",
                       ("pooled", "sens_at_99"))],
        "forbidden": [],
    },
    # The Fragmentomics row of the cross-platform table quoted the HCC_J
    # OvR AUC as 0.7532 -- the value the per-cancer artifact held on
    # 2026-09-23. The artifact was regenerated on 2026-11-16 and the
    # row's own source (`per_cancer_sens_at_spec.json`, HCC_J) now reads
    # 0.7221. Nothing caught it because that row was not guarded: the
    # document carried a figure from a run the artifact no longer
    # described. The pin holds the row to the artifact it names, and
    # FORBIDS the superseded value so the old figure cannot come back.
    {
        "doc": "docs/CROSS_PLATFORM_BENCHMARK.md",
        "anchor": "| Fragmentomics |",
        "required": [("per_cancer_sens_at_spec.json",
                       ("per_cancer", "HCC_J", "auc_mean"))],
        # The superseded CI bound from the 2026-09-23 run. It is not in
        # the artifact any more, so it is quoted as a literal -- a
        # forbidden value is compared as a NUMBER, and 0.688 was the
        # old auc_ci[0]. If the row ever quotes it again, the row has
        # gone back to a run the artifact does not describe.
        "forbidden": [("per_cancer_sens_at_spec.json", ("__literal__", 0.688))],
        "source": "results/per_cancer_sens_at_spec.json",
    },
    # The n=5 pair's paired deltas (rows 3 and 4) are DERIVED from the
    # two n=5 per-run arrays -- the n=5 ablation artifact was
    # overwritten by the n=20 run, so no artifact stores them. They
    # are therefore guarded by a recomputation test rather than by a
    # row pin: `test/test_label_control.py` recomputes both deltas
    # from the arrays and requires the document to quote them.
    #
    # The n=5 per-run artifacts themselves are the historical
    # leaky-CV record: the leaky StratifiedKFold they were produced
    # with is gone from the code, so the driver cannot regenerate
    # them and must not be asked to. Pinning them here would put the
    # reproduction driver in the impossible position of having to
    # re-create a superseded protocol, which is why the n=5 rows are
    # guarded by recomputation instead.
    # The decision-curve table in RESULTS.md §4. Its rows mix two
    # artifacts on purpose -- the decision-curve CLI does not compute
    # the lr_fusion strategy, so that column comes from the fusion
    # benchmark -- and the text says so. What must not happen is the row
    # carrying a threshold from a different run of the decision curve.
    {
        "doc": "RESULTS.md",
        "anchor": "| **95%** |",
        "required": [
            ("decision_curve.json",
             ("per_strategy", "tumor_naive", "threshold@0.95")),
            ("decision_curve.json",
             ("per_strategy", "naive_average", "threshold@0.95")),
        ],
        "forbidden": [],
    },
    {
        "doc": "RESULTS.md",
        "anchor": "| **99%** |",
        "required": [
            ("decision_curve.json",
             ("per_strategy", "tumor_naive", "threshold@0.99")),
            ("decision_curve.json",
             ("per_strategy", "naive_average", "threshold@0.99")),
        ],
        "forbidden": [],
    },
    # The TL;DR rows of RESULTS.md. Three tumor-naive benchmarks exist
    # (0.9684 GC-corrected, 0.9770 legacy, 0.9785 with motifs) and the
    # headline row used to quote one of the others -- the row pin is
    # what keeps it on the headline one.
    {
        "doc": "RESULTS.md",
        "anchor": "| Tumor-naive fragmentomics (627 FinaleDB cohort) |",
        "required": [("cross_study_finallydb_gc_corrected.json",
                       ("pooled", "harmonized", "auc_mean"))],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("pooled", "harmonized", "auc_mean")),
            ("cross_study_finallydb_with_motifs.json",
             ("pooled", "harmonized", "auc_mean")),
        ],
    },
    {
        "doc": "RESULTS.md",
        "anchor": "| Mutation + tumor-naive naive-average |",
        "required": [("fusion_ablation_isotonic.json",
                       ("per_strategy", "naive_average", "auc_mean"))],
        "forbidden": [],
    },
    # The harmonization row quotes the true-confound control's two
    # states. Both come from the GC-corrected artifact; the legacy
    # artifact's pair is a different run.
    {
        "doc": "RESULTS.md",
        "anchor": "| Cross-study harmonization effect |",
        "required": [
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "harmonized", "auc_mean")),
            ("cross_study_finallydb_gc_corrected.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "no_harmonize", "auc_mean")),
        ],
        "forbidden": [
            ("cross_study_finallydb.json",
             ("true_confound_control", "cancer_jiang_healthy_cristiano",
              "harmonized", "auc_mean")),
        ],
    },
    # The methods paper's section 3.3 per-cancer rows. Their three
    # sensitivity columns are pinned to the committed per-cancer artifact,
    # which carries them for every cancer. They were not pinned before,
    # and the gap had teeth: a blanket search-and-replace of the suite's
    # test count rewrote OV's sens@95 from 0.964 to 0.967 -- a different
    # quantity -- and every guard stayed green.
    #
    # The AUC column and the bootstrap CI are NOT pinned here: they come
    # from the paper's 5-seed run, and the committed artifact holds a
    # single-seed AUC (auc_std 0) and its own bootstrap draw. Pinning them
    # to that artifact would fail on numbers that are not wrong, so they
    # stay outside this guard.
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| HCC_J | 89 | 353 (jiang only) |",
        "required": [
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "HCC_J", "sens_at_95")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "HCC_J", "sens_at_98")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "HCC_J", "sens_at_99")),
        ],
        "forbidden": [],
        "source": None,
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| LUAD | 79 | 343 |",
        "required": [
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "LUAD", "sens_at_95")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "LUAD", "sens_at_98")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "LUAD", "sens_at_99")),
        ],
        "forbidden": [],
        "source": None,
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| PAAD | 60 | 324 |",
        "required": [
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "PAAD", "sens_at_95")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "PAAD", "sens_at_98")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "PAAD", "sens_at_99")),
        ],
        "forbidden": [],
        "source": None,
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| BRCA | 53 | 317 |",
        "required": [
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "BRCA", "sens_at_95")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "BRCA", "sens_at_98")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "BRCA", "sens_at_99")),
        ],
        "forbidden": [],
        "source": None,
    },
    {
        "doc": "paper/METHODS_PAPER.md",
        "anchor": "| OV | 28 | 292 |",
        "required": [
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "OV", "sens_at_95")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "OV", "sens_at_98")),
            ("per_cancer_sens_at_spec.json",
             ("per_cancer", "OV", "sens_at_99")),
        ],
        "forbidden": [],
        "source": None,
    },
]


def load(name: str, results_dir: str | None = None):
    path = os.path.join(results_dir or RESULTS_DIR, name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


# Every decimal number in a text, with the precision it was written at.
_NUMBER_RE = re.compile(r"(?<![\d.])[-−]?\d+\.\d+(?![\d])")
# An integer written with a thousands separator, e.g. "103,863". A
# count like this is written EXACTLY -- there is no rounding to allow
# for -- so it is matched exactly rather than at a decimal precision.
# Without this, a figure the documents quote in that form is invisible
# to the guard: the decimal-only pattern skips it, so it could drift
# freely while the guard stayed green.
_THOUSANDS_RE = re.compile(r"(?<![\d.,])[-−]?\d{1,3}(?:,\d{3})+(?![\d,.])")


def numbers_in(text: str):
    """Yield (value, decimals_written) for each decimal number in ``text``.

    A document may quote a benchmark value at the precision its table
    uses, so a bare substring search is too blunt: "0.970" cannot
    distinguish the GC-corrected publication-8 AUC (0.9698) from the
    legacy one (0.9704), because both round to it. Comparing each
    written number against the target at the precision it was written
    at separates them: "0.9698" and "0.9704" are distinct, and only
    the exact rendering matches.

    A sign is part of the number, and the manuscripts write a negative
    figure with the typographic minus (U+2212) as well as with the
    ASCII hyphen -- so both are accepted. Without the sign, a negative
    figure parsed as its positive magnitude would never match the
    artifact it came from, and the row pin would fail on a document
    that quotes it correctly.

    A thousands-separated integer yields ``decimals = 0`` and is
    matched EXACTLY by the caller: "103,863" names one value, so a
    tolerance of half a unit would let "103,864" pass.
    """
    for match in _NUMBER_RE.finditer(text):
        raw = match.group(0)
        decimals = len(raw.split(".")[1])
        yield (float(raw.replace("−", "-")), decimals)
    for match in _THOUSANDS_RE.finditer(text):
        raw = match.group(0)
        yield (float(raw.replace(",", "").replace("−", "-")), 0)


def quotes(text: str, value: float, tol: float = 1e-9) -> bool:
    """Whether ``text`` quotes ``value``, at the precision it was written.

    A written number matches when it is the target rounded to the
    written precision -- i.e. within half a unit of the last written
    digit. A number written to more digits than the target carries is
    accepted only if it rounds to the target at that precision, so
    "0.48360" matches 0.4836 but "0.9704" does not match 0.9698.

    At least three decimals are required (the repository's convention
    for quoting a benchmark figure). With fewer, the rendering is too
    coarse to identify a value: "0.50" cannot tell 0.4958 from 0.5033,
    so accepting it would make both the presence check and the
    forbidden-value check unreliable.

    A thousands-separated integer is the exception: it is written
    EXACTLY, so it is matched exactly rather than at a decimal
    precision.
    """
    for written, decimals in numbers_in(text):
        if decimals == 0:
            if abs(written - value) <= tol:
                return True
            continue
        if decimals < 3:
            continue
        if abs(written - value) <= 0.5 * (10.0 ** -decimals) + tol:
            return True
    return False


def line(label, value, source):
    print(f"  {label:<46} {value}   [{source}]")


def dig(artifact, path):
    node = artifact
    for key in path:
        if key == "__literal__":
            # A literal value, not a key to walk. Used for a FORBIDDEN
            # number that is no longer in the artifact -- the superseded
            # figure a row must not quote again.
            continue
        node = node[key]
    return node


def dig_number(artifact, path):
    """Resolve a path to a number, including the fusion per_strategy shape.

    A path of ``("per_strategy", "<strategy>", "auc_mean")`` reaches into
    the fusion artifact's list-of-metric-rows shape, which a plain key
    walk cannot do. ``("per_strategy", "<strategy>", "threshold@<spec>")``
    reaches into the decision-curve artifact's list-of-operating-points
    shape the same way.

    A path of ``("__literal__", <number>)`` returns the number itself.
    That is how a forbidden value is expressed when it is NOT in the
    artifact any more -- the superseded figure a row must not quote
    again.
    """
    if len(path) == 2 and path[0] == "__literal__":
        return float(path[1])
    if (len(path) == 3 and path[0] == "per_strategy"):
        if path[2] == "auc_mean":
            return fusion_strategy_metric(artifact, path[1], "auc")
        if path[2].startswith("threshold@"):
            spec = float(path[2].split("@", 1)[1])
            return float(decision_curve_point(
                artifact, path[1], spec)["operating_threshold"])
        if path[2] in ("sens_at_95", "sens_at_99"):
            return fusion_strategy_metric(artifact, path[1], path[2])
    return float(dig(artifact, path))


def decision_curve_point(artifact, strategy, spec):
    """One operating point of the decision-curve artifact.

    ``per_specificity`` is a list of
    ``{"specificity", "sensitivity", "operating_threshold"}`` rows, so a
    plain path lookup cannot reach into it.
    """
    for row in artifact["per_strategy"][strategy]["per_specificity"]:
        if abs(row["specificity"] - spec) < 1e-9:
            return row
    raise KeyError(f"no specificity {spec} for {strategy}")


def loso_direction(artifact, train, test):
    for entry in artifact["loso"]:
        if (entry["train_study"] == train and entry["test_study"] == test
                and entry["harmonization"] == "self_zscore"):
            return float(entry["pooled_test_auc"]["auc"])
    raise KeyError(f"{train} -> {test} not in the LOSO artifact")


def fusion_auc(artifact, strategy):
    for row in artifact["per_strategy"][strategy]:
        if row["metric"] == "auc":
            return float(row["mean"])
    raise KeyError(f"no auc row for {strategy}")


def fusion_strategy_metric(artifact, strategy, metric):
    """One metric of one fusion strategy, from the per_strategy list.

    ``per_strategy`` maps a strategy name to a list of
    ``{"metric": ..., "mean": ..., "std": ...}`` rows, so a plain
    path lookup cannot reach into it.
    """
    for row in artifact["per_strategy"][strategy]:
        if row["metric"] == metric:
            return float(row["mean"])
    raise KeyError(f"no {metric} row for {strategy}")


def loso_in_domain(artifact, study):
    """The within-publication CV AUC for one study, from the LOSO artifact."""
    for row in artifact["in_domain"]:
        if row["study"] == study:
            return float(row["cv_auc_mean"])
    raise KeyError(f"no in_domain row for {study}")


def special_value(artifact_name, artifact, what):
    """Resolve one ``QUOTED_SPECIAL`` entry to (value, label)."""
    if what.startswith("loso "):
        train, test = (part.strip() for part in what[5:].split("->"))
        return (loso_direction(artifact, train, test),
                f"{artifact_name} {train} -> {test}")
    if what.startswith("in_domain "):
        study = what[len("in_domain "):].strip()
        return (loso_in_domain(artifact, study),
                f"{artifact_name} in_domain {study}")
    if what.startswith("decision curve "):
        rest = what[len("decision curve "):].strip()
        metric = "threshold"
        if rest.endswith(" threshold"):
            rest = rest[:-len(" threshold")].strip()
        strategy, _, spec = rest.rpartition(" @")
        strategy = strategy.strip()
        row = decision_curve_point(artifact, strategy, float(spec))
        if metric == "threshold":
            return (float(row["operating_threshold"]),
                    f"{artifact_name} {strategy} @ {spec} threshold")
        return (float(row["sensitivity"]),
                f"{artifact_name} {strategy} @ {spec} sensitivity")
    if what.startswith("fusion "):
        strategy, _, metric = what[len("fusion "):].strip().rpartition(" ")
        return (fusion_strategy_metric(artifact, strategy, metric),
                f"{artifact_name} {strategy} {metric}")
    return fusion_auc(artifact, what), what


def quoted_artifacts():
    """The artifact file names the documents quote numbers from.

    ``paper/REPRODUCE.sh`` regenerates every one of these; the test
    ``test/test_paper_reproduce_driver.py`` asserts the driver keeps up
    with this list, so a newly quoted artifact cannot be left behind by
    the reproduction driver.

    The table-row pins count too: an artifact pinned only by a row
    (the with-motifs benchmark was) is quoted by the documents just as
    much as one pinned by number, and the driver has to regenerate it
    either way.
    """
    names = [name for name, _path, _docs in QUOTED_NUMBERS]
    names += [name for name, _what, _doc in QUOTED_SPECIAL]
    for row in TABLE_ROWS:
        for artifact_name, _path in list(row["required"]) + list(row["forbidden"]):
            names.append(artifact_name)
    return names


def quoted_entries(results_dir: str | None = None):
    """Yield (label, value, artifact, key, documents) for every quoted figure.

    ``key`` is the pair the figure is filed under -- ``(artifact, path)``
    for a number and ``(artifact, what)`` for a special lookup -- which
    is the key ``EXCUSED_COPIES`` uses. It is what lets a document be
    excused from one figure and stay protected for another.
    """
    for artifact_name, path, documents in QUOTED_NUMBERS:
        artifact = load(artifact_name, results_dir)
        if artifact is None:
            continue
        named = list(documents)
        for doc in EXTRA_DOCS.get((artifact_name, tuple(path)), []):
            if doc not in named:
                named.append(doc)
        yield (f"{artifact_name}:{'.'.join(map(str, path))}",
               dig_number(artifact, path), artifact_name,
               (artifact_name, tuple(path)), named)
    for artifact_name, what, document in QUOTED_SPECIAL:
        artifact = load(artifact_name, results_dir)
        if artifact is None:
            continue
        value, label = special_value(artifact_name, artifact, what)
        yield (label, value, artifact_name, (artifact_name, what), [document])
    for artifact_name, what, document in QUOTED_SPECIAL_EXTRA:
        artifact = load(artifact_name, results_dir)
        if artifact is None:
            continue
        value, label = special_value(artifact_name, artifact, what)
        yield (label, value, artifact_name, (artifact_name, what), [document])


def quoted_values(results_dir: str | None = None):
    """Yield (label, value, artifact, documents) for every quoted number."""
    for label, value, artifact_name, _key, named in quoted_entries(results_dir):
        yield (label, value, artifact_name, named)


def table_row_values(results_dir: str | None = None):
    """Yield (doc, anchor, required, forbidden, source) for each table row.

    ``required`` and ``forbidden`` are lists of (artifact_name, value)
    pairs: the numbers the row must quote, and the numbers it must not
    (the same quantity measured by a different run).
    """
    for case in TABLE_ROWS:
        required = []
        for artifact_name, path in case["required"]:
            artifact = load(artifact_name, results_dir)
            if artifact is None:
                continue
            if isinstance(path, str) and path.startswith("loso "):
                value, _ = special_value(artifact_name, artifact, path)
            else:
                value = dig_number(artifact, path)
            required.append((artifact_name, value))
        forbidden = []
        for artifact_name, path in case["forbidden"]:
            artifact = load(artifact_name, results_dir)
            if artifact is None:
                continue
            forbidden.append((artifact_name, dig_number(artifact, path)))
        yield (case["doc"], case["anchor"], required, forbidden,
               case.get("source"))


def main() -> int:
    print("== numbers the documents quote (authoritative list)")
    for label, value, artifact_name, documents in quoted_values():
        line(label, f"{value:.4f}", artifact_name)
        print(f"      quoted in: {', '.join(documents)}")
    print()

    gc = load("cross_study_finallydb_gc_corrected.json")
    legacy = load("cross_study_finallydb.json")
    pc = load("per_cancer_sens_at_spec.json")

    for tag, art, src in (("GC-corrected (headline)", gc,
                           "cross_study_finallydb_gc_corrected.json"),
                          ("legacy (no GC correction)", legacy,
                           "cross_study_finallydb.json")):
        if not art:
            continue
        print(f"== {tag}")
        pooled = art["pooled"]["harmonized"]
        line("pooled harmonized AUC",
             f"{pooled['auc_mean']:.4f} +- {pooled['auc_std']:.4f}", src)
        line("pooled no_harmonize AUC",
             f"{art['pooled']['no_harmonize']['auc_mean']:.4f}", src)
        line("sens@95 / @98 / @99",
             " / ".join(f"{pooled[f'sens_at_spec_{s}']:.3f}"
                        for s in (95, 98, 99)), src)
        for pub, row in (art.get("per_publication")
                         or art.get("per_cohort", {})).items():
            name = row.get("study", pub)
            line(f"per-publication {name} (pub {pub})",
                 f"{row['auc_mean']:.4f} +- {row['auc_std']:.4f}", src)
        tc = (art.get("true_confound_control", {})
              .get("cancer_jiang_healthy_cristiano", {}))
        if tc and isinstance(tc.get("harmonized"), dict):
            line("true-confound harm / no-harm",
                 f"{tc['harmonized']['auc_mean']:.4f} / "
                 f"{tc['no_harmonize']['auc_mean']:.4f}", src)
        print()

    if pc:
        print("== per-cancer standalone")
        line("pooled AUC (DeLong CI)",
             f"{pc['pooled']['auc_mean']:.4f} "
             f"[{pc['pooled']['auc_ci'][0]:.4f}, {pc['pooled']['auc_ci'][1]:.4f}]",
             "per_cancer_sens_at_spec.json")
        for cancer, row in sorted(pc.get("per_cancer", {}).items()):
            line(f"{cancer} AUC",
                 f"{row['auc_mean']:.4f} "
                 f"[{row['auc_ci'][0]:.4f}, {row['auc_ci'][1]:.4f}]",
                 "per_cancer_sens_at_spec.json")
            s = row.get("sens_at_spec", {})
            if s:
                line(f"{cancer} sens@95/98/99",
                     " / ".join(f"{s[k].get('sensitivity', float('nan')):.3f}"
                                for k in ("0.95", "0.98", "0.99") if k in s),
                     "per_cancer_sens_at_spec.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
