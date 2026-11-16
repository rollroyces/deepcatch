"""Real methylation cohorts for the GNN branch — no synthetic fallback.

Why this module exists
----------------------
``src/methylation_gnn/`` was architecture-only: its graph builder filled
the node set with randomly generated regions and its data module with
randomly generated beta values, so the "methylation" the GNN saw had no
relationship to any real sample. The V3 pipeline (``scripts/
methylation_cohort_channel.py``) now runs FinaleMe end to end on real
public FinaleDB fragments and writes per-CpG decoded calls, so the GNN can
be pointed at real data.

This module is the seam. It reads FinaleMe's decoded output -- the same
files the cohort channel script reads -- and returns the structure
``RegulatoryGraphBuilder.build_graph`` expects, plus the real CpG
coordinates so the graph's nodes are real loci rather than synthetic
placeholders.

What is and is not real
-----------------------
Real: the beta values, the CpG coordinates, the per-region methylation
variance across samples, the labels, the sample identities.

Not real: ``cpg_density``, ``gc_content`` and ``cpg_obs_exp``. The
FinaleMe output does not carry region-level annotations, so those node
features are left out of the returned dict and the graph builder's
existing "missing data is filled with zeros and logged" contract applies.
They are NOT filled with random numbers, which is what the synthetic path
did -- a fabricated GC content would make a real-data run look better
calibrated than it is.

Provenance is explicit: the returned dict carries ``data_source`` and
``is_synthetic``, so a downstream artifact can state where its numbers
came from.
"""
from __future__ import annotations

import gzip
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# The decoded files carry the methylation call as an integer percentage in
# column 4 (0-based index 3).
_BETA_COLUMN = 3
_BETA_SCALE = 100.0

# The shared-locus coverage threshold. Matches
# ``scripts/methylation_cohort_channel.py`` (COVER_FRACTION), so the
# cohort this module loads is the same cohort the channel was measured
# on -- the two must not disagree about which loci are "shared".
DEFAULT_COVERAGE_FRACTION = 0.8
DEFAULT_MIN_SHARED_CPGS = 50


class RealDataError(RuntimeError):
    """Raised when the real cohort cannot be loaded.

    Deliberately an error and not a warning: a GNN run that silently
    falls back to synthetic regions would produce numbers that look real
    and are not.
    """


def _read_decoded(path: str) -> Dict[Tuple[str, str], float]:
    """One sample's decoded CpG calls: ``{(chrom, start): beta}``."""
    out: Dict[Tuple[str, str], float] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= _BETA_COLUMN:
                continue
            try:
                beta = float(parts[_BETA_COLUMN]) / _BETA_SCALE
            except ValueError:
                continue
            out[(parts[0], parts[1])] = beta
    return out


def _read_labels(labels_tsv: str) -> Dict[str, Tuple[int, str]]:
    """``{sample: (label, disease)}`` from the cohort labels TSV."""
    labels: Dict[str, Tuple[int, str]] = {}
    with open(labels_tsv, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "sample":
                continue
            try:
                label = int(parts[1])
            except ValueError:
                continue
            disease = parts[2] if len(parts) > 2 else "unknown"
            labels[parts[0]] = (label, disease)
    return labels


def _decoded_paths(decoded_dir: str) -> List[str]:
    if not os.path.isdir(decoded_dir):
        raise RealDataError(
            f"decoded directory does not exist: {decoded_dir}")
    paths = [
        os.path.join(decoded_dir, name)
        for name in sorted(os.listdir(decoded_dir))
        if name.endswith(".decoded.txt.gz")
    ]
    if not paths:
        raise RealDataError(
            f"no *.decoded.txt.gz files in {decoded_dir}; run "
            f"scripts/methylation_cohort_channel.py first")
    return paths


def load_finaleme_cohort(
    decoded_dir: str,
    labels_tsv: str,
    *,
    coverage_fraction: float = DEFAULT_COVERAGE_FRACTION,
    min_shared_cpgs: int = DEFAULT_MIN_SHARED_CPGS,
) -> Dict[str, Any]:
    """Load a real methylation cohort from FinaleMe's decoded output.

    Parameters
    ----------
    decoded_dir : str
        Directory of ``<sample>.decoded.txt.gz`` files.
    labels_tsv : str
        TSV with ``sample<TAB>label(0/1)<TAB>disease``.
    coverage_fraction : float
        A locus joins the shared set when at least this fraction of the
        samples cover it. The default matches the cohort channel script.
    min_shared_cpgs : int
        Refuse to build a cohort with fewer shared loci than this.

    Returns
    -------
    dict with:
        ``beta_values``       (n_regions, n_samples) float32, NaN where a
                              sample does not cover the locus.
        ``covered``           (n_regions, n_samples) bool.
        ``labels``            (n_samples,) int64.
        ``sample_names``      list[str].
        ``disease``           list[str].
        ``regions``           list[GenomicRegion] -- one per shared locus.
        ``cpg_keys``          list[(chrom, start)].
        ``methylation_variance`` (n_regions,) float32, variance across the
                              samples that cover the locus.
        ``data_source``       "finaleme_decoded".
        ``is_synthetic``      False.
        ``region_annotations`` note about which node features are zero-filled.

    Raises
    ------
    RealDataError
        If the cohort cannot be loaded. Never a silent fallback.
    """
    # Imported here so this module can be used for the data-only path
    # without pulling in the graph builder's dependencies.
    from .graph_builder import GenomicRegion

    paths = _decoded_paths(decoded_dir)
    labels = _read_labels(labels_tsv)
    if not labels:
        raise RealDataError(f"no labels parsed from {labels_tsv}")

    decoded: Dict[str, Dict[Tuple[str, str], float]] = {}
    for path in paths:
        sample = os.path.basename(path)[: -len(".decoded.txt.gz")]
        calls = _read_decoded(path)
        if calls:
            decoded[sample] = calls

    # Keep only the samples the labels file knows about: an unlabelled
    # sample cannot be scored and would silently shift the label vector.
    samples = sorted(s for s in decoded if s in labels)
    dropped = sorted(set(decoded) - set(samples))
    if dropped:
        logger.warning("dropping %d decoded samples with no label: %s",
                       len(dropped), ", ".join(dropped[:5]))
    if len(samples) < 2:
        raise RealDataError(
            f"only {len(samples)} labelled samples in {decoded_dir}; "
            f"need at least 2")

    # The shared loci: covered by at least `coverage_fraction` of the
    # samples.
    coverage: Dict[Tuple[str, str], int] = {}
    for sample in samples:
        for key in decoded[sample]:
            coverage[key] = coverage.get(key, 0) + 1
    need = max(2, int(np.ceil(coverage_fraction * len(samples))))
    threshold_loci = sorted(k for k, c in coverage.items() if c >= need)
    selection = "coverage_threshold"
    shared = threshold_loci
    if len(shared) < min_shared_cpgs:
        # Same fallback as ``scripts/methylation_cohort_channel.py``: the
        # most widely covered loci, when too few reach the coverage
        # threshold. On the published cohort the threshold yields only 16
        # loci, so the 1,550 the paper quotes is THIS pool -- the most
        # widely covered loci, each covered by a different subset of
        # samples, not 1,550 loci covered by 80% of them.
        #
        # The two paths must select the same loci, or the GNN's cohort and
        # the channel's cohort would not be the same cohort.
        widest = max(50, len(samples) * 50)
        shared = [k for k, _ in
                  sorted(coverage.items(), key=lambda kv: (-kv[1], kv[0]))
                  [:widest]]
        selection = "most_widely_covered_fallback"
        logger.warning(
            "only %d loci reach the %d/%d coverage threshold; using the "
            "%d most widely covered loci instead (same rule as "
            "scripts/methylation_cohort_channel.py)",
            len(threshold_loci), need, len(samples), len(shared))
    if len(shared) < min_shared_cpgs:
        raise RealDataError(
            f"only {len(shared)} loci could be selected from {decoded_dir}; "
            f"the loader needs at least {min_shared_cpgs} to build a graph")

    n_regions = len(shared)
    n_samples = len(samples)
    beta = np.full((n_regions, n_samples), np.nan, dtype=np.float32)
    covered = np.zeros((n_regions, n_samples), dtype=bool)
    for j, sample in enumerate(samples):
        calls = decoded[sample]
        for i, key in enumerate(shared):
            if key in calls:
                beta[i, j] = calls[key]
                covered[i, j] = True

    # Per-locus variance across the samples that actually cover it. A
    # locus covered by one sample has no variance, which is 0 here and
    # flagged by `covered` rather than by a magic number.
    counts = covered.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.nanmean(np.where(covered, beta, np.nan), axis=1)
        variance = np.nanvar(
            np.where(covered, beta - means[:, None], np.nan), axis=1)
    variance = np.nan_to_num(variance, nan=0.0).astype(np.float32)

    regions = [
        GenomicRegion(
            chrom=chrom,
            start=int(start),
            end=int(start) + 1,
            region_type="cpg",
        )
        for chrom, start in shared
    ]

    cohort = {
        "beta_values": beta,
        "covered": covered,
        "labels": np.array([labels[s][0] for s in samples], dtype=np.int64),
        "sample_names": samples,
        "disease": [labels[s][1] for s in samples],
        "regions": regions,
        "cpg_keys": shared,
        "methylation_variance": variance,
        "data_source": "finaleme_decoded",
        "is_synthetic": False,
        "n_locus_observations": int(counts.sum()),
        "locus_selection": selection,
        "coverage_threshold": need,
        "n_threshold_loci": len(threshold_loci),
        "region_annotations": (
            "cpg_density, gc_content and cpg_obs_exp are not carried by "
            "the FinaleMe output and are therefore absent here; the graph "
            "builder zero-fills them and logs it. They are not filled "
            "with random values, which is what the synthetic path did."),
    }
    logger.info(
        "loaded a REAL methylation cohort: %d samples, %d shared loci, "
        "%d observations", n_samples, n_regions, int(counts.sum()))
    return cohort


def per_sample_methylation_data(
    cohort: Dict[str, Any],
    sample_index: int,
) -> Dict[str, np.ndarray]:
    """The ``methylation_data`` dict for one sample of a real cohort.

    Same contract as ``build_graph`` expects: one beta value per region,
    plus the per-region methylation variance. The uncovered loci are
    reported as NaN, which the graph builder's feature assembly turns
    into the missing-data fill it already documents.
    """
    beta = cohort["beta_values"][:, sample_index]
    return {
        "beta_values": np.where(np.isnan(beta), 0.0, beta).astype(np.float32),
        "methylation_variance": cohort["methylation_variance"],
        "beta_is_observed": cohort["covered"][:, sample_index].astype(
            np.float32),
    }


def build_real_cohort_graphs(
    builder,
    cohort: Dict[str, Any],
    *,
    max_samples: Optional[int] = None,
) -> List[Any]:
    """Build one PyG graph per sample of a real cohort.

    ``builder`` is a ``RegulatoryGraphBuilder`` whose regions have already
    been set from the cohort's real loci (see ``load_real_regions``).
    """
    graphs = []
    n = len(cohort["sample_names"])
    if max_samples is not None:
        n = min(n, max_samples)
    for index in range(n):
        methylation = per_sample_methylation_data(cohort, index)
        graphs.append(
            builder.build_graph(
                sample_name=cohort["sample_names"][index],
                methylation_data=methylation,
                label=int(cohort["labels"][index]),
                cancer_type=cohort["disease"][index],
            )
        )
    return graphs


def load_real_regions(builder, cohort: Dict[str, Any]) -> None:
    """Point a graph builder at the cohort's real CpG loci.

    After this call the builder will NOT fall back to synthetic regions:
    the loci are real, so a later "no regions loaded" state is a bug, not
    a default.
    """
    builder.load_regions_from_list(cohort["regions"])
    builder.require_real_regions = True
