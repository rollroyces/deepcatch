# Cross-Platform FinaleMe Pipeline — Operational Guide

**TL;DR.** When the FinaleDB 5-channel fragmentomics cache is present
AND the FinaleMe Java tool is fully installed (JAR + pretrained models
+ 1.4 GB reference files), running
`scripts/cross_platform_finaledb_validate.py` emits a real
cross-platform AUC using the methylation β-values imputed by FinaleMe
on top of the existing fragmentomics baseline. **When any input is
missing, the pipeline refuses to print a cross-platform AUC** —
the JSON carries `data_source` and a `refusal_reason` instead.

This is the **scaffolding** for the cross-platform pipeline. The
two preamble scripts (`scripts/finaleme_pipeline.py`,
`scripts/finaleme_to_deepcatch_bridge.py`,
`scripts/cross_platform_finaledb_validate.py`) are honest by
construction: they probe inputs, surface what's missing in plain
language, and never synthesize predictions.

---

## What it WILL do (when fully installed)

```bash
env -u PYTHONPATH ./.venv/bin/python scripts/cross_platform_finaledb_validate.py \
    --finaleme-dir ~/.hermes/.local/finaleme/output \
    --features-dir /Users/hermes/cfdna-fragmentomics-pipeline/data/features \
    --labels-tsv /Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv \
    --output results/cross_platform_readiness.json
```

Emits a JSON with:

```json
{
  "schema_version": "1.0",
  "data_source": "both",
  "methylation_available": true,
  "fragmentomics_available": true,
  "n_overlap_samples": N,
  "methylation_channel_auc": <float>,
  "fragmentomics_channel_auc": <float>,
  "cross_platform_auc": <float>,
  "fusion_strategies": {
    "naive_average": {"auc_mean": …, "auc_std": …, "per_seed_auc": […]},
    "logit_average": {"auc_mean": …, …},
    "lr_fusion":     {"auc_mean": …, …}
  },
  …
}
```

The output is the runnable cross-platform fusion that the
methodology demo (`cross_platform_methylation_fusion.json`) called
out as "a true cross-platform fusion would require paired
fragmentomics + methylation measurements on the SAME patients".
The structural pipeline is now in place; running it on real data
requires the FinaleMe installation below.

## What it WON'T do

When **any** required input is missing, the validator emits a
not-runnable JSON and **does not print** a `cross_platform_auc`:

```json
{
  "data_source": "methylation_only",  // or fragmentomics_only / no_data
  "cross_platform_auc": null,
  "refusal_reason": "Methylation channel is missing — install FinaleMe JAR + pretrained models + reference files."
}
```

This is the load-bearing honesty gate. The pipeline never
synthesizes fake per-sample β-values from the published AUCs to
make a number look good — `scripts/finaleme_to_deepcatch_bridge.py`
validates β-values are in [0, 1] and that sample IDs overlap with
the FinaleDB features cache.

## Current state (as of the last `probe` run)

Run the always-safe diagnostic:

```bash
env -u PYTHONPATH ./.venv/bin/python scripts/_cross_platform_readiness_probe.py
```

Expected output on this machine **today** (Sept 2026):

```
[OK]        Java runtime           openjdk 21.x.x
[MISSING]   FinaleMe JAR            no JAR under ~/.hermes/.local/finaleme/
[MISSING]   Pretrained HMM models   ~/.hermes/.local/finaleme/models/ missing
[MISSING]   Reference files         hg19.2bit, methylation prior bw, CpG motif, mappability BED, chrom sizes
[MISSING]   FinaleMe output         ~/.hermes/.local/finaleme/output/ missing
[OK]        FinaleDB features       627 samples × 5 channels

data_source: fragmentomics_only
cross-platform AUC will NOT be printed
```

This is the honest verdict. Installing the FinaleMe suite (below)
moves the verdict to `data_source: both`.

## Prerequisites (~3 hours of one-time setup)

### 1. Java 21 (tarball, no `brew install`)

`brew install --cask zulu@21` fails on this M4 because
`/opt/homebrew` is root-owned. Use the tarball recipe:

```bash
JDK_DIR=~/.hermes/.local/jdk
mkdir -p "$JDK_DIR"
curl -L -o /tmp/zulu21.tar.gz \
  "https://cdn.azul.com/zulu/bin/zulu21.34.18-ca-jdk21.0.4-macosx_aarch64.tar.gz"
tar -xzf /tmp/zulu21.tar.gz -C "$JDK_DIR"
export JAVA_HOME="$JDK_DIR/zulu21.34.18-ca-jdk21.0.4-macosx_aarch64"
export PATH="$JAVA_HOME/bin:$PATH"
java -version   # openjdk version "21.0.4"
```

Persist by adding those two exports to `~/.zshrc`.

### 2. FinaleMe JAR (hybrid recipe, ~30 min)

The Zenodo pretrained models serialize against a legacy Java
package (`main.java.edu.mit.compbio.ccinference.*`) that no public
FinaleMe JAR's `LegacyPackageObjectInputStream` remaps. The hybrid
recipe combines **v0.61's streaming decoder** with **v0.58.1's
TreeMap-based HMM classes** (renamed under v0.61's package to
avoid split-package conflicts) plus an additional legacy-package
remap.

Five steps (full recipe in the cfdna-fragmentomics skill's
references/honest-benchmarking.md FinaleMe pitfall):

1. Audit legacy packages with `javap -p` on both the model file
   and the candidate JARs.
2. Patch `LegacyPackageObjectInputStream.remapLegacyClassName()`
   to add a remap for `main.java.edu.mit.compbio.ccinference.*` →
   `edu.northwestern.epifluidlab.finaleme.hmm.*`.
3. Check field-layout compatibility (pi/a as
   `TreeMap<Integer, Double[][]>` vs primitive `double[][]`) —
   v0.58.1 matches; v0.61 doesn't.
4. Build a hybrid JAR with v0.61's streaming decoder + v0.58.1's
   HMM classes (renamed). Compile together with `mvn clean package`.
5. Verify biological direction: cancer model gives mean β ~2-3 pp
   LOWER than healthy on real cfDNA (global hypomethylation).

Drop the resulting JAR at:
```bash
mkdir -p ~/.hermes/.local/finaleme
cp FinaleMe-hybrid-jar-with-dependencies.jar ~/.hermes/.local/finaleme/
```

**Do NOT ship the 50 MB hybrid JAR to a public release** — commit
the patched sources instead and rebuild locally.

### 3. Pretrained HMM models from Zenodo (~150 kB total)

```bash
mkdir -p ~/.hermes/.local/finaleme/models
cd ~/.hermes/.local/finaleme/models
curl -L -O "https://zenodo.org/records/14013719/files/healthy_WGS.mincg7.example.hmm_model?download=1"
curl -L -O "https://zenodo.org/records/14013719/files/cancer_WGS.mincg7.example.hmm_model?download=1"
ls -lh   # ~148 kB total
```

### 4. Reference files (~1.4 GB, ~90 sec parallel download)

```
~/.hermes/.local/finaleme/hg19.2bit                                      778 MB
~/.hermes/.local/finaleme/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw     324 MB
~/.hermes/.local/finaleme/CG_motif.hg19.common_chr.pos_only.bedgraph.gz 148 MB
~/.hermes/.local/finaleme/mappability.bed                                ~3 MB
~/.hermes/.local/finaleme/hg19.chrom.sizes                               <1 MB
```

The cfdna-fragmentomics skill's references/finaledb-api.md has
the parallel-download recipe. Verify SHA-256 against the Zenodo
manifest before use (the methylation prior bw is the most likely
to silently mismatch).

## Running the pipeline (when prerequisites are met)

```bash
# Optional: decode per-sample β-values for a single sample (32 sec/sample)
env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_pipeline.py \
    features --input BH01.chr22.frag.bed.gz --output-dir ~/.hermes/.local/finaleme/output/
env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_pipeline.py \
    decode \
    --features-file ~/.hermes/.local/finaleme/output/BH01_chr22_cpg_features.hg19.bed.gz \
    --output-dir ~/.hermes/.local/finaleme/output/ \
    --model both

# Validate cross-platform readiness (the canonical entry point)
env -u PYTHONPATH ./.venv/bin/python scripts/cross_platform_finaledb_validate.py
# → results/cross_platform_readiness.json with cross_platform_auc when both channels present
```

## Honest ceilings

Even with a perfect installation, the realistic upper-bound AUC
addition from methylation on low-pass cfDNA is **≤ +0.005 AUC**
above the fragmentomics baseline (skill: ceiling is intrinsic
biology). The methylation channel is most useful for **tissue-of-origin
classification**, not cancer-vs-healthy — and there's no public
cohort of meaningful size with paired fragmentomics + methylation
on the **same** patients as of 2026-Q3.

## Cross-platform fusion contract (3 strategies)

When `data_source == 'both'`, the validator computes three fusion
strategies (matching the methodology demo's
`cross_platform_methylation_fusion.json`):

1. **Naive average**: `0.5 · (p_meth + p_frag)` per sample.
2. **Logit average**: `sigmoid(0.5 · (logit p_meth + logit p_frag))`.
3. **LR fusion**: sklearn `LogisticRegression(C=1.0)` head on the
   concatenated logit-transformed per-channel scores, evaluated via
   **5-fold StratifiedKFold pooled OOF across 5 seeds**.

Same metric, same input, same labels — only the aggregation
differs. Report all three, not just the headline.

## Relationships

| Tool                                              | Purpose                              |
| ------------------------------------------------- | ------------------------------------ |
| `scripts/_cross_platform_readiness_probe.py`        | Always-safe diagnostic; exits 0.     |
| `scripts/finaleme_pipeline.py`                     | Orchestrates FinaleMe Steps 1 + 3.   |
| `scripts/finaleme_to_deepcatch_bridge.py`          | Converts per-sample β-values → JSON. |
| `scripts/cross_platform_finaledb_validate.py`      | Cross-platform orchestrator.         |
| `scripts/_finaledb_feature_loader.py`              | 5-channel feature loader helper.     |
| `results/cross_platform_readiness.json`            | Default output (operator-checked).   |
| `results/cross_platform_finaledb_bridge.json`      | Intermediate bridge output.          |
| `docs/CROSS_PLATFORM_FINALEME.md` (this file)      | Operator guide.                      |

The methodology demo's `results/cross_platform_methylation_fusion.json`
is **not** touched by the new pipeline (it uses disjoint cohorts and
synthesized per-sample scores — see its `provenance.honest_caveats`).
