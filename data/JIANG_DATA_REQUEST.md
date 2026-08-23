# Jiang Lab — Data Request Template (for re-running nested-CV)

**Recipient:** Prof. Dennis Lo / Prof. Peiyong Jiang, CUHK
**Dataset needed:** Table S1 — 4-mer end-motif frequencies, 129 plasma DNA samples
(38 Control + 91 Cancer: HCC 34, LC 10, HNSCC 10, CRC 10, NPC 10, + HBV carriers)
**Paper:** Jiang et al., "Plasma DNA End-Motif Profiling as a Fragmentomic Marker
in Cancer, Pregnancy, and Transplantation", *Cancer Discovery* 2020,
DOI 10.1158/2159-8290.CD-19-0622

---

## Email Draft

> Subject: DeepCatch re-analysis of 4-mer end-motif data (Table S1) — data access request

> Dear Prof. Jiang and colleagues,

> You may recall the DeepCatch analysis of your 4-mer end-motif dataset that we
> shared with you earlier this year (`summary_for_professor_jiang.md` — HCC
> AUC 0.9845 with the two-stage CET approach). We have since completed a
> thorough methodological review of our pipeline, and found that the original
> analysis selected the top-50 discriminative motifs **before** cross-validation
> splitting — a form of feature-selection leakage that inflates the reported
> AUC slightly.

> We have now fixed the pipeline with proper nested cross-validation
> (motif selection inside each training fold), and would like to re-run the
> analysis on the same 129-sample Table S1 dataset to obtain an honest,
> leakage-free performance estimate. This is a quick re-analysis (~1 hour of
> compute) and we will share the updated summary with you before any
> publication.

> Could you share the Table S1 file (4-mer end-motif frequency matrix, 129
> samples × 256 motifs, with sample labels) via a secure channel? We will use
> it solely for this validation and will not redistribute it.

> Thank you for your continued collaboration.

> Best regards,
> Royce
> DeepCatch Project — github.com/rollroyces/deepcatch

---

## Once the file arrives (one command)

```bash
# Place the file at:
#   data/deepcatch_data.xlsx        (or set DEEPCATCH_DATA_DIR)

cd ~/deepcatch
source .venv/bin/activate
python run_jiang_analysis.py \
    -i data/deepcatch_data.xlsx \
    -o results/jiang_nested_cv/ \
    --top-k 50 --seed 42 --nested-cv --report
```

**What to expect:** per-cancer-type nested-CV AUCs in
`results/jiang_nested_cv/summary_report.md`. The old (leaky) HCC AUC was
0.9845; the honest nested-CV value will be somewhat lower — that is the
expected and correct outcome.
