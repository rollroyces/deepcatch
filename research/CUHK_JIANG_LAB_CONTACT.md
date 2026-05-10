# CUHK Jiang Lab — Contact & Collaboration Prep

**Generated:** 2026-05-10  
**Purpose:** Contact Peiyong Jiang's lab at CUHK for 4-mer end motif frequency data collaboration (DeepCatch FragmentoSign validation)  
**Status:** DRAFT — Ready for Royce review before sending

---

## 1. ⚠️ Citation Correction (IMPORTANT)

DeepCatch code (multiple files) incorrectly cites:

> "Jiang, P. et al. (2020). Nature Genetics 52:712-719. PMID: 32514122"

**This is WRONG on two counts:**
1. PMID 32514122 is a Japanese GWAS paper (Ishigaki et al.), NOT a Jiang paper
2. The 4-mer end motif paper is in **Cancer Discovery**, NOT Nature Genetics

**CORRECT citation:**
> Jiang P\*, Sun K\*, Peng W\*, Cheng SH, Ni M, Yeung PC, Heung MMS, Xie T, Shang H, Zhou Z, Chan RWY, Wong J, Wong VWS, Poon LC, Leung TY, Lam WKJ, Chan JYK, Chan HLY, Chan KCA, Chiu RWK, Lo YMD. **Plasma DNA End-Motif Profiling as a Fragmentomic Marker in Cancer, Pregnancy, and Transplantation.** *Cancer Discov.* 2020;10(5):664-673. DOI: [10.1158/2159-8290.CD-19-0622](https://doi.org/10.1158/2159-8290.CD-19-0622). PMID: [32111602](https://pubmed.ncbi.nlm.nih.gov/32111602/).

**Files to fix in DeepCatch:**
- `src/fragmentomics/bam_motif_extractor.py` (lines 50, 112, 182, 393)
- `src/fragmentomics/fragment_gmm.py` (line 59)
- `src/fragmentomics/__init__.py` (line 38)
- `src/fragmentomics/normalization.py` (line 302)
- `src/fragmentomics/themis_features.py` (line 178)

Recommended action: Fix the cite key and PMID before contacting the Jiang lab — citing the wrong journal would look sloppy.

---

## 2. Lab & PI Information

| Field | Detail |
|---|---|
| **PI Name** | **Peiyong Jiang (江培勇)** — NOT "Peng Jiang" |
| **Position** | Professor, Department of Chemical Pathology |
| | Director, Biomedical Computing Centre |
| | Li Ka Shing Institute of Health Sciences (LiHS) |
| **Institution** | The Chinese University of Hong Kong (CUHK) |
| **Address** | Dept of Chemical Pathology, Prince of Wales Hospital, 30-32 Ngan Shing Street, Shatin, New Territories, Hong Kong SAR |
| **Email** | **jiangpeiyong@cuhk.edu.hk** |
| **Phone** | Not publicly listed — check lab page |
| **Lab Page** | https://www.cuhk.edu.hk/med/cpy/Research/PeiyongJiang.htm |
| **Research Profile** | https://research.cuhk.edu.hk/en/persons/peiyong-jiang/ |
| **Google Scholar** | https://scholar.google.com/citations?user=RIC8c6cAAAAJ |
| **Group** | CNARG (Circulating Nucleic Acids Research Group) |
| **Lab Head** | **Y.M. Dennis Lo (盧煜明)** — Corresponding author on most papers |
| **Key Collaborator** | Wenlei Peng (彭文磊) — co-first author on end motif paper |
| **Notable Recognition** | Top 20 Translational Researchers of 2019 (Nature Biotechnology, 2020) |
| | Senior Member, National Academy of Inventors (first HK scholar) |

### Lab Research Focus
1. Bioinformatics software for NGS & third-generation sequencing data analysis
2. AI/deep learning for cfDNA fragmentomics
3. Noninvasive prenatal testing (NIPT) & cancer detection from plasma DNA
4. cfDNA fragmentation biology (end motifs, nucleosome positioning, methylation)
5. Single-cell transcriptomics + plasma cfRNA

### Most Relevant Publications (chronological)
1. Jiang P et al. (2020) — **Plasma DNA End-Motif Profiling** — *Cancer Discov* 10:664-673 ← primary target
2. Jiang P et al. (2026) — **Holistic determination of ends of cfDNA molecules** — *Cell Genomics* ← latest update
3. Jiang P et al. (2018) — **Preferred end coordinates in HCC** — *PNAS* 115:E10925
4. Jiang P et al. (2015) — **Lengthening/shortening of plasma DNA in HCC** — *PNAS* 112:E1317
5. Cristiano S et al. (2019) — DELFI genome-wide fragmentation — *Nature* 570:385 (collaborator group)

---

## 3. Data Availability Assessment

### 3.1 Primary Dataset

| Detail | Info |
|---|---|
| **EGA Accession** | **EGAD00001005093** |
| **Title** | "Plasma DNA motif analysis" |
| **Samples** | 118 samples |
| **Platform** | Illumina HiSeq 4000 |
| **DAC** | EGAC00001000078 (CNARG) |
| **Contact** | jiangpeiyong@cuhk.edu.hk |
| **Data Type** | Likely FASTQ (raw sequencing data) |
| **Format** | Not specified — expect FASTQ or BAM |
| **Contains** | HCC vs non-HCC plasma cfDNA sequencing data |

### 3.2 Related Datasets

| EGA ID | Title | Samples | Tech | Relevance |
|---|---|---|---|---|
| EGAD00001006054 | plasma dna fragmentations | 29 | HiSeq 4000 | HCC + pregnant women, includes buffy coat & tumor tissue |
| EGAD00001005088 | (EBV methylation) | 236 | NextSeq 500 | Nasopharyngeal carcinoma — different cancer type, shared methodology |
| EGAD00001005071 | DNASE1L3 deletion & end motifs | 41 | NextSeq 500 | Mouse model — validates end motif biology |
| EGAD00001004561 | Plasma DNA without PCR amplification | 169 | HiSeq 2000 | Library-prep-free cfDNA — gold standard for fragmentomics |

### 3.3 Access Process

1. **Data is CONTROLLED ACCESS** — NOT publicly downloadable
2. Must sign **CNARG Data Access Agreement (DAA)** via EGA
3. DAC reviews application — Jiang Peiyong is the DAC contact
4. Processing time: typically 2-4 weeks for EGA DAC applications
5. Requirements:
   - Institutional affiliation (or independent researcher justification)
   - Non-commercial research purpose
   - Specific project description
   - Agreement to acknowledge CNARG/CUHK in publications

### 3.4 What's NOT Publicly Available

- ❌ No processed 4-mer frequency matrices in supplementary materials
- ❌ No GEO/SRA deposit — data is EGA-only
- ❌ No pre-computed MDS values per sample (must compute from raw data)
- ❌ No 256-dim motif frequency vectors shared as supplementary tables

### 3.5 Alternative Approach (Preferable)

Instead of requesting raw FASTQ/BAM data (which requires formal DAA + weeks of review), we could ask for:

1. **Processed 4-mer frequency vectors**: 256-dim vectors per sample (HCC vs control labels)
   - Much smaller data transfer
   - No IRB/consent issues (aggregate, de-identified frequency data)
   - Faster turnaround
2. **MDS values per sample**: Even simpler — just one scalar per sample with clinical metadata
3. **Summary statistics**: Mean ± SD per motif across cohorts (published in paper figures)

This approach aligns with the paper's results (they show MDS and motif rankings) and would be sufficient to validate DeepCatch's FragmentoSign pipeline against their published findings.

---

## 4. Technical Integration Plan

### 4.1 Jiang Lab Methodology → DeepCatch Mapping

| Jiang Lab Protocol | DeepCatch FragmentoSign | Status |
|---|---|---|
| Extract 4-mer from cfDNA fragments (5' end) | `bam_motif_extractor.extract_4mer_end_motifs()` | ✅ Implemented |
| 256 possible 4-mers (AAAA → TTTT) | `ALL_4MERS` array, `MOTIF_TO_IDX` mapping | ✅ Identical |
| Strand orientation (to + strand) | `stranded=True` with `_reverse_complement()` | ✅ Implemented |
| cfDNA fragment length filter (90-250 bp) | `fragment_length_min=90`, `fragment_length_max=250` | ✅ Matched |
| MAPQ ≥ 30 | `mapq_threshold=30` | ✅ Matched |
| Motif Diversity Score (MDS) | `compute_MDS()` — normalized Simpson index | ✅ Identical formula |
| HCC → CCCA depletion signature | MDS captures diversity shift; CCCA can be isolated | ✅ Compatible |
| Tissue-of-origin clustering by motif profiles | TOO classifier uses motif enrichment scores | ✅ Compatible |

### 4.2 Expected Data Format

Assuming Jiang lab can share processed frequency data:

```python
# Expected format (per sample)
{
    'sample_id': 'HCC_001',
    'group': 'HCC',  # or 'healthy', 'HBV_carrier', 'cirrhosis'
    'motif_frequencies': np.array([256])  # normalized 4-mer frequencies
}

# Or as a table:
# sample_id | group | AAAA | AAAC | AAAG | AAAT | ... | TTTT
# HCC_001   | HCC   | 0.0042 | 0.0031 | ...
```

### 4.3 Validation Pipeline

1. **Load Jiang data** → 256-dim vectors per sample
2. **Run through FragmentoSign**:
   - Compute MDS from Jiang's frequency vectors
   - Compare MDS distributions (HCC vs control) against published values
3. **Run our BAM extractor on Jiang's raw data** → generate our own frequency vectors
4. **Compare**:
   - Correlation of motif frequencies (Pearson r per sample)
   - Concordance of MDS values
   - CCCA depletion signal in HCC
5. **Integrate into fusion classifier**:
   - Use Jiang's validated motif features as one modality
   - Compare single-modality vs fused performance

### 4.4 Key Metrics to Validate

| Metric | Jiang et al. 2020 (Cancer Discov) | DeepCatch to Validate |
|---|---|---|
| MDS in HCC vs non-HCC | Significant increase in HCC | Reproduce direction & magnitude |
| CCCA abundance | Much lower in HCC | Reproduce fold-change |
| Top discriminatory motifs | CCCA, others reported | Rank correlation |
| AUC for HCC detection | Reported in paper (need exact value) | Compare with FragmentoSign single-modality AUC |

### 4.5 Integration with Fusion Classifier

```python
# Proposed integration:
from fragmentomics import extract_4mer_end_motifs, compute_MDS

# Option A: Use Jiang lab's pre-computed frequencies as validation
jiang_frequencies = load_jiang_data()  # 256-dim per sample
our_frequencies = extract_4mer_end_motifs(bam, ref)['motif_frequencies']

# Option B: Use motif features as a modality in DeepCatch fusion
motif_features = {
    'mds': compute_MDS(our_frequencies),
    'motif_256': our_frequencies,  # full vector for classifier
    'top_discriminatory': our_frequencies[top_motifs],  # reduced set
}
# Feed into performance_weighted_fusion()
```

---

## 5. Draft Email (English) — DEPRECATED (replaced by Section 9)

See Section 9 for the final polished version.

---

## 6. Draft Email (Chinese Version) — DEPRECATED (replaced by Section 9)

See Section 9 for the final polished version.

---

## 7. Concrete Next Steps

### Before Contacting (Recommended)

- [x] **Fix citation in DeepCatch codebase** — replace "Nature Genetics 52:712-719" with "Cancer Discovery 10:664-673" and correct PMID to 32111602 across all 7 files ✅ (fixed 2026-05-10)
- [ ] **Update paper bibliography** — ensure `jiang2026holistic` reference is correct; add the Cancer Discovery 2020 paper as a separate citation if not already present
- [ ] Fix Issue #1 from validation report (BAM FASTA caching) — improves credibility if Jiang lab reviews our code
- [ ] Prepare a 1-page technical summary of FragmentoSign (ready to attach/share)

### Contact Strategy

| Step | Action | Timeline |
|---|---|---|
| 1 | Send email to jiangpeiyong@cuhk.edu.hk | After citation fix |
| 2 | If no reply in 2 weeks, follow up + CC Wenlei Peng (co-first author) | +2 weeks |
| 3 | Alternative: contact via CUHK Chemical Pathology dept office | +4 weeks |
| 4 | Consider reaching out through Dennis Lo's lab (corresponding author on most papers) | Fallback |
| 5 | If formal DAA needed: apply through EGA (EGAC00001000078) | 2-4 week process |

### Alternative Data Sources (If Jiang Lab Cannot Share)

| Source | What | Access |
|---|---|---|
| **ICGC-ARGO** | cfDNA fragmentomics from HCC cohorts | May have public motif data |
| **GRAIL CCGA** | Multi-cancer cfDNA methylation + fragmentomics | Controlled access only |
| **TCGA** | Tissue-derived fragmentomics proxies | Public (but not cfDNA) |
| **Snyder et al. 2016 (Cell)** | cfDNA nucleosome positioning | GEO GSE71378 (public) |
| **Cristiano et al. 2019 (Nature)** | DELFI fragmentation patterns | EGA EGAS00001003266 (controlled) |
| **Mock data approach** | Generate synthetic 4-mer frequencies based on published MDS distributions | Immediate (but lower impact) |

### Risk Assessment

| Risk | Likelihood | Mitigation |
|---|---|---|
| No reply from cold email | Medium-High | Follow up; try different channels |
| Data cannot be shared (IRB/consent) | Medium | Ask for summary statistics instead |
| Requires institutional affiliation | Medium | Emphasize open-source, non-commercial nature |
| Jiang has left CUHK | Low | Lab page shows active professor status; recent 2026 publications |
| Time zone / language barrier | Low | English email is standard; lab publishes in English |

---

## 8. Appendix: Lab Context & Organizational Structure

Jiang Peiyong works within the **Circulating Nucleic Acids Research Group (CNARG)** led by **Prof. Y.M. Dennis Lo (盧煜明)** — the pioneer of NIPT and cfDNA diagnostics. The group is based at:

- **Department of Chemical Pathology**, CUHK
- **Li Ka Shing Institute of Health Sciences (LiHS)**
- **Prince of Wales Hospital**, Shatin, Hong Kong

The lab's work spans the full cfDNA pipeline: wet-lab protocols → sequencing → bioinformatics → clinical validation. Jiang leads the bioinformatics/computational arm. Co-first author Wenlei Peng (彭文磊) handles much of the motif analysis pipeline. Dennis Lo is the senior author and decision-maker on data sharing.

If direct communication with Jiang does not yield results, approaching through the Lo lab infrastructure (e.g., contacting the department's administrative office) may be more effective.

---

## 9. Final Email (Ready to Send)

**Status:** ✅ REVIEWED & POLISHED  
**Date:** 2026-05-10  
**Reviewer:** Ironman (subagent)

---

### 9.1 Deliverability Check

| Check | Result |
|---|---|
| **Email** `jiangpeiyong@cuhk.edu.hk` | ✅ Confirmed — consistent with CUHK format, listed on EGA DAC EGAC00001000078, appears on multiple data access pages |
| **PI still at CUHK?** | ✅ Yes — active in 2025-2026 (Cell Genomics 2026, Cancer Cell 2025, Communications Biology 2025) |
| **Lab page up?** | ✅ https://www.cuhk.edu.hk/med/cpy/Research/PeiyongJiang.htm — live, lists correct publications |
| **Research profile?** | ✅ https://research.cuhk.edu.hk/en/persons/peiyong-jiang/ — active, recent pubs listed |
| **Alternative contacts** | Wenlei Peng (co-first author, handles motif pipeline), Dennis Lo (lab head), CUHK Chemical Pathology dept office |

⚠️ **Flag:** `contact@deepcatch.org` — verify this email exists before sending. If not, use GitHub profile URL or a personal email Royce is comfortable with. A non-deliverable reply-to address hurts credibility.

---

### 9.2 Subject Line Options (A/B Test)

| Option | Subject | Rationale |
|---|---|---|
| **A (Recommended)** | `DeepCatch: Validating 4-mer cfDNA End Motif Pipeline Against Your Cancer Discovery 2020 Data` | Specific, references their paper, action-oriented. Best for academic cold email. |
| **B** | `Question about 4-mer End Motif Data — DeepCatch Fragmentomics Validation` | Softer open, "Question" piques curiosity. Slightly less specific. |
| **C (Short)** | `DeepCatch × Jiang Lab: 4-mer cfDNA Motif Validation` | Brief, clear collaboration framing. Good if Jiang's inbox is high-volume. |

**Recommendation:** Use **Option A** for first send. If no reply in 2 weeks, use **Option B** for follow-up (it looks like a different email, reducing spam-filter grouping).

---

### 9.3 Final Polished Email — English (Copy-Paste Ready)

**To:** jiangpeiyong@cuhk.edu.hk  
**Subject:** DeepCatch: Validating 4-mer cfDNA End Motif Pipeline Against Your Cancer Discovery 2020 Data

---

Dear Professor Jiang,

I'm Royce, an independent researcher building **DeepCatch** — an open-source multimodal framework for pan-cancer early detection from cfDNA (github.com/rollroyces/deepcatch, MIT license).

Your group's work on plasma DNA end motif profiling — from the Cancer Discovery 2020 paper establishing the 4-mer methodology through to your recent Cell Genomics 2026 paper on holistic determination of cfDNA ends — has been foundational to our fragmentomics module. We've implemented the extraction protocol (strand-aware orientation, normalized Simpson MDS, 90–250 bp filtering, MAPQ ≥ 30), validated it through unit testing, and now want to benchmark against your clinical data.

**Specifically, we'd be grateful for one of the following (in order of preference):**

1. **Processed 4-mer frequency vectors** (256-dim per sample) from your HCC vs non-HCC plasma samples with group labels — aggregate, de-identified frequencies; no raw sequencing data, minimal privacy concern, small transfer
2. **Per-sample MDS values** with clinical annotations — sufficient to validate the direction and magnitude of the cancer-specific diversity shift
3. **Guidance** on whether EGAD00001005093 is the most suitable dataset, or if a different repository better fits our purpose

**Our commitment:** Fully open-source (MIT), non-commercial. Your work will be prominently cited (targeting Bioinformatics or PLOS Computational Biology). We'll share all validation results before publication, and we're happy to sign any data access agreements.

If data sharing isn't feasible, I'd still greatly appreciate your advice on publicly available 4-mer datasets or whether published MDS distributions could serve as validation targets.

Thank you for considering this — and for your foundational contributions to cfDNA fragmentomics.

Best regards,
Royce
DeepCatch Project
GitHub: github.com/rollroyces/deepcatch

---

### 9.4 Final Polished Email — 中文版 (Copy-Paste Ready)

**收件人：** jiangpeiyong@cuhk.edu.hk  
**主旨：** DeepCatch項目：請求4-mer cfDNA末端基序數據協作驗證

---

江教授您好：

我是DeepCatch的獨立開發者Royce。DeepCatch是一個開源的多模態cfDNA泛癌早篩框架（github.com/rollroyces/deepcatch，MIT許可證），目標是建立可重現、社區共同維護的液體活檢工具。

貴團隊在cfDNA末端基序分析領域的系列工作——從2020年Cancer Discovery確立4-mer分析方法，到2026年Cell Genomics對cfDNA末端的全面鑑定——是我們FragmentoSign模塊的核心基礎。我們已完整實現了4-mer提取流程（鏈定向處理、歸一化Simpson MDS、90–250 bp片段過濾、MAPQ ≥ 30），並通過單元測試驗證，現希望與貴團隊的臨床數據進行基準對比。

**我們希望獲取以下之一（按優先級排列）：**

1. HCC vs 對照組血漿樣本的**處理後4-mer頻率向量**（每樣本256維，附分組標籤）——僅需聚合頻率，無原始測序數據，私隱風險極低，數據量小
2. 或**每樣本MDS值**及臨床標註——已足以驗證癌症特異性多樣性變化的方向與幅度
3. 關於EGAD00001005093是否為最合適數據集的建議

**我們的承諾：** 完全開源（MIT）、非商業用途。論文將顯著引用貴團隊工作（目標期刊Bioinformatics或PLOS Computational Biology）。發表前與您分享所有驗證結果。樂意簽署任何數據訪問協議。

如數據無法共享，也懇請指點公開可用的4-mer數據集，或建議替代驗證方案。

期待您的回覆。感謝您在cfDNA fragmentomics領域的卓越貢獻！

Royce 敬上
DeepCatch項目
GitHub: github.com/rollroyces/deepcatch

---

### 9.5 Key Improvements from Draft

| Area | Draft | Final |
|---|---|---|
| **Length** | ~400 words (too long for cold email) | ~250 words (optimal for academic cold email) |
| **Opening** | "I am writing as an independent researcher developing…" (stiff) | "I'm Royce, an independent researcher building…" (warm, direct) |
| **Recent work** | Only referenced 2020 paper | References both 2020 Cancer Discovery AND 2026 Cell Genomics (shows diligence) |
| **Ask structure** | Unordered list, unclear priority | Numbered by preference, each option self-contained |
| **Data justification** | Separate paragraph on "why processed data" | Inline justification per option ("no raw sequencing, minimal privacy") |
| **"Why this matters"** | Separate vague paragraph | Removed — let the ask speak for itself |
| **Contact info** | `contact@deepcatch.org` (⚠️ unverified) | GitHub profile (reliable) — with flag to verify email |
| **Chinese version** | Machine-translated feel, 沒有引用2026論文 | Natural flow, 加入Cell Genomics 2026, 語氣更地道 |

---

### 9.6 When to Send

**Best window:** HK business hours (HKT = UTC+8)

| Day | Window | Rationale |
|---|---|---|
| **Tue–Thu** | **9:00–11:00 HKT** (01:00–03:00 UTC) | Best time — academics clear inbox in the morning, mid-week avoids Monday rush and Friday wind-down |
| **Mon** | 10:00–11:30 HKT | OK, but Monday inboxes are crowded |
| **Fri** | Before 10:00 HKT | Avoid — emails sent Friday afternoon get buried over weekend |
| **Sat–Sun** | ❌ Avoid | Weekend email risks being ignored or perceived as unprofessional |

**Optimal send time:** **Tuesday or Wednesday, 9:30 AM HKT**

**Timezone conversion for Royce (HKT):** Just send during your morning. No conversion needed.

---

### 9.7 Follow-Up Plan

| Timeline | Action | Details |
|---|---|---|
| **Day 0** | Send email | Use Subject A, Tue/Wed 9:30 AM HKT |
| **+14 days** | **Follow-up #1** — Reply to original thread | Keep it brief: "Dear Professor Jiang, I'm following up on the email below regarding DeepCatch's 4-mer end motif validation. I understand you're very busy — if you could point me toward the most appropriate channel for this request, I'd be very grateful. Best, Royce" |
| **+21 days** | **Follow-up #2** — New email to Jiang, **CC Wenlei Peng** (co-first author, handles motif pipeline) | Slightly different framing: "I wanted to follow up on my earlier email regarding DeepCatch. Since Wenlei Peng was co-first author on the end motif methodology, I'm CC'ing in case the data query falls within his domain." |
| **+35 days** | **Escalate** — Contact CUHK Chemical Pathology department office | Phone or general department email (cpy@cuhk.edu.hk). Ask for guidance on reaching Prof. Jiang's group for a research collaboration inquiry. |
| **+42 days** | **Last resort** — Contact Dennis Lo's lab | Only if all above fails. Frame as: "I've been trying to reach Prof. Jiang regarding a non-commercial, open-source validation study referencing your group's end motif work. Could you advise on the best contact channel?" |

**Important:** At each follow-up, reference the original email date. Don't send more than 3 follow-ups total — after that, accept it as a "no" and pivot to alternative data sources (Section 7).

---

### 9.8 Pre-Send Checklist

Before hitting send, Royce should verify:

- [x] **Fix DeepCatch citation** — Change "Nature Genetics 52:712-719" → "Cancer Discovery 10:664-673" and PMID → 32111602 across all source files ✅ (fixed 2026-05-10)
- [ ] **Verify reply-to email** — `contact@deepcatch.org` may not exist. Use a real, deliverable email or just GitHub URL
- [ ] **Check GitHub repo is public** — Jiang will likely visit github.com/rollroyces/deepcatch
- [ ] **README looks professional** — Ensure the repo README explains DeepCatch clearly for an academic visitor
- [ ] **FragmentoSign code is clean** — Consider fixing Issue #1 (BAM FASTA caching) before sending; Jiang's lab may review the code
- [ ] **Prepare 1-page technical summary** — Have it ready if Jiang responds asking for more details
- [ ] **Send at right time** — Tue/Wed, 9:00–11:00 HKT
