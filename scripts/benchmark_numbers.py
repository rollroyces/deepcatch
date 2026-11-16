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
"""
from __future__ import annotations

import json
import os
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
     ["paper/paper.tex", "paper.md", "README.md",
      "docs/EXTERNAL_VALIDATION.md"]),
    # The case-only publication is reported as a sensitivity at the
    # threshold fixed in training, so the threshold and that sensitivity
    # are quoted numbers too.
    ("external_validation_finaledb.json",
     ("training_derived_threshold", "threshold"),
     ["paper/METHODS_PAPER.md", "docs/EXTERNAL_VALIDATION.md"]),
    ("external_validation_finaledb.json",
     ("per_publication_auc", "9", "sensitivity_at_train_threshold_99spec"),
     ["paper/METHODS_PAPER.md", "docs/EXTERNAL_VALIDATION.md"]),
    ("cross_study_finallydb_gc_corrected.json",
     ("pooled", "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md", "paper.md", "README.md", "MODEL.md",
      "ULTRA_EARLY_READINESS.md", "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md"]),
    ("cross_study_finallydb.json",
     ("pooled", "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md", "README.md", "ULTRA_EARLY_READINESS.md",
      "docs/CROSS_STUDY_BENCHMARK.md"]),
    ("methylation_cohort_channel.json",
     ("beta_summary_auc_loo",),
     ["paper/paper.tex", "paper.md", "README.md", "MODEL.md",
      "ULTRA_EARLY_READINESS.md", "docs/METHYLATION_PIPELINE.md"]),
]

# Numbers that need a lookup rather than a plain path.
QUOTED_SPECIAL = [
    ("cross_study_loso.json", "loso jiang -> cristiano",
     "paper/paper.tex"),
    ("cross_study_loso.json", "loso cristiano -> jiang",
     "paper/paper.tex"),
    ("fusion_ablation_isotonic.json", "naive_average",
     "paper/paper.tex"),
]


def load(name: str, results_dir: str | None = None):
    path = os.path.join(results_dir or RESULTS_DIR, name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def line(label, value, source):
    print(f"  {label:<46} {value}   [{source}]")


def dig(artifact, path):
    node = artifact
    for key in path:
        node = node[key]
    return node


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


def quoted_artifacts():
    """The artifact file names the documents quote numbers from.

    ``paper/REPRODUCE.sh`` regenerates every one of these; the test
    ``test/test_paper_reproduce_driver.py`` asserts the driver keeps up
    with this list, so a newly quoted artifact cannot be left behind by
    the reproduction driver.
    """
    names = [name for name, _path, _docs in QUOTED_NUMBERS]
    names += [name for name, _what, _doc in QUOTED_SPECIAL]
    return names


def quoted_values(results_dir: str | None = None):
    """Yield (label, value, artifact, documents) for every quoted number."""
    for artifact_name, path, documents in QUOTED_NUMBERS:
        artifact = load(artifact_name, results_dir)
        if artifact is None:
            continue
        yield (f"{artifact_name}:{'.'.join(map(str, path))}",
               float(dig(artifact, path)), artifact_name, documents)
    for artifact_name, what, document in QUOTED_SPECIAL:
        artifact = load(artifact_name, results_dir)
        if artifact is None:
            continue
        if what.startswith("loso "):
            train, test = (part.strip() for part in what[5:].split("->"))
            value = loso_direction(artifact, train, test)
            label = f"{artifact_name} {train} -> {test}"
        else:
            value = fusion_auc(artifact, what)
            label = f"{artifact_name} {what}"
        yield (label, value, artifact_name, [document])


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
