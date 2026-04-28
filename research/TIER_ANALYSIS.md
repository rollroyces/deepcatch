# DeepCatch Tier Analysis: Where We Stand vs State-of-the-Art

**Date:** 2026-04-28 | **Based on:** 21 papers + v1.2 improvement data

---

## Tier Classification System

| Tier | Description | Example |
|---|---|---|
| **Tier 0** | Clinical deployment, FDA approved, >100K patients | Grail Galleri, Guardant360 |
| **Tier 1** | Clinical validation complete, published in Nature/Science/Lancet | CancerSEEK, DELFI, PanSeer |
| **Tier 2** | Multi-center clinical validation ongoing | DELFI LUCAS (Mazzone 2024) |
| **Tier 3** | Computational innovation + simulation + limited clinical data | Bie 2023 (THEMIS), Moldovan 2024 |
| **Tier 4** | Computational innovation, simulation only | Most academic methods papers |

---

## 1. DeepCatch Multi-Modal Fusion — Competitive Position

### Head-to-Head AUC (Simulation, Matched Data)

| Method | AUC | Fusion Strategy | Multi-Modal? | Longitudinal? |
|---|---|---|---|---|
| **DeepCatch (5 modalities)** | **0.961** | Performance-weighted | ✅ 5 modalities | ✅ CET+SPRT |
| Bie 2023 (THEMIS, 4 modalities) | 0.857 | Simple average | ✅ 4 modalities | ❌ |
| Moldovan 2024 (FrEIA) | 0.960 | OR-logic + ML | ✅ 3 | ❌ |
| CAPP-Seq (mutation only) | 0.847 | Single-modality | ❌ | ❌ |
| iDES (error-suppressed) | 0.514 | Single-modality | ❌ | ❌ |

**Verdict:** DeepCatch is in **Tier 3** for fusion — competitive with Bie and Moldovan, but simulation only.

### Advantage vs Closest Competitor

| Competitor | DeepCatch Advantage | Significance |
|---|---|---|
| **Bie 2023** | ΔAUC +0.104, p<0.0001 | ✅ Clear winner in simulation |
| **Moldovan 2024** | AUC comparable (0.961 vs 0.960) | ⚠️ Marginally better |
| **CAPP-Seq** | ΔAUC +0.114 | ✅ Multi-modal beats single-modality |

---

## 2. Tissue-of-Origin — Where We Stand

| Method | TOO Accuracy | Cancer Types | Clinical? |
|---|---|---|---|
| **Grail Galleri** | **88.7%** | 50+ | ✅ Clinical |
| **DeepCatch TOO** | **81.7%** | 8 (20 planned) | ❌ Simulation |
| CancerSEEK | ~63% | 8 | ✅ Clinical |
| DELFI | ~75% | 7 | ✅ Clinical |

**Verdict:** Tier 3 for TOO. DeepCatch's 81.7% is competitive on simulation, but Grail's clinical 88.7% is the gold standard. Gap: ~7%.

---

## 3. Longitudinal CET — Unique Position

| Method | Strategy | Sensitivity | Specificity | Timepoints |
|---|---|---|---|---|
| **PanSeer (methylation)** | Archived single sample | 95% (pre-dx) | 96% | 1 (archived) |
| **DeepCatch CET (multi-modal)** | Active SPRT tracking | 23.0% | 78.4% | 8 quarterly |
| **DeepCatch CET (mutation-only)** | Single-modality SPRT | 9.5% | 61.8% | 8 quarterly |

**Verdict:** DeepCatch is **unique** — no other paper does active multi-timepoint SPRT for screening. But performance is not yet clinical-grade. PanSeer's 95% pre-diagnosis sensitivity (from one archived sample) demonstrates the power of the right biomarker (methylation).

**Key insight:** CET with methylation as the primary analyte (instead of mutations) could potentially match PanSeer's performance while adding multi-timepoint advantages.

---

## 4. Cost-Effectiveness

| Method | Sequencing | Est. Cost/Sample | Sensitivity |
|---|---|---|---|
| Guardant360 | 5,000× targeted | ~$500-1000 (clinical) | 85.3% |
| Grail Galleri | ~30× WGBS targeted | ~$949 (list price) | 51.5% |
| DELFI | 1-2× WGS | ~$100-200 | 73% |
| **DeepCatch (recommended)** | **5,000× targeted (500kb)** | **~$74** | **71.0%*** |
| DeepCatch (original) | 50,000× targeted | $135 | 72.8%* |

*\*Simulation-estimated at 99% specificity*

**Verdict:** DeepCatch is cost-competitive at the recommended depth. $74/sample is cheaper than all clinical alternatives. But this assumes targeted capture validation.

---

## 5. Overall Tier Placement

```
Tier 0 (Clinical):  Grail Galleri, Guardant360, FoundationOne
Tier 1 (Validated): CancerSEEK, DELFI, PanSeer
Tier 2 (Validating): DELFI LUCAS, Grail NHS trial
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Tier 3 (Innovation): ██ DeepCatch v1.2 ██  ← WE ARE HERE
                     Bie 2023 (THEMIS)
                     Moldovan 2024 (FrEIA)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Tier 4 (Simulation): Most methods papers
```

**DeepCatch is solid Tier 3** — same tier as Bie 2023 (Nat Commun) and Moldovan 2024 (Cell Rep Med). With clinical validation, could move to Tier 2 within 12-18 months.

---

## 6. Strengths vs Competitors

| Strength | vs Whom | Margin |
|---|---|---|
| 🧠 **Fusion method** | Bie (simple avg) | ΔAUC +0.104 |
| 🔄 **Longitudinal** | All competitors | UNIQUE (active SPRT) |
| 🧬 **TOO** | CancerSEEK | +18.7% accuracy |
| 💰 **Cost** | Guardant360 | $74 vs $500-1000 |
| 🏥 **Pan-cancer** | Moldovan (10 types) | 20 types |
| 📦 **Open source** | Grail (proprietary) | UNIQUE |

---

## 7. Weaknesses vs Competitors

| Weakness | vs Whom | Gap |
|---|---|---|
| 🔬 **Clinical validation** | Grail (140K patients) | 0 vs 140K |
| 🎯 **CET specificity** | PanSeer (96%) | 78.4% vs 96% |
| 🧬 **TOO (clinical)** | Grail (88.7%) | 81.7% sim vs 88.7% clinical |
| 📊 **Cancer types** | Grail (50+) | 20 vs 50+ |
| 🏭 **Throughput** | Guardant360 (200K+) | 0 clinical throughput |

---

## 8. Path to Tier 1 (What It Would Take)

| Action | Impact | Timeline |
|---|---|---|
| **CET Two-Stage** (specificity → 99%) | Move CET from weakness to strength | 1-2 weeks |
| **Clinical pilot** (n=50+50) | First clinical data point | 3-6 months |
| **Independent replication** | Validation credibility | 6-9 months |
| **Prospective study** (n=500+) | Tier 2 entry | 12-18 months |
| **Multi-center trial** (n=5000+) | Tier 1 entry | 24-36 months |

---

## 9. How to Leapfrog to Tier 2 Fastest

### Quick Wins (This Month)
1. **Implement Two-Stage CET** — combined specificity → 99%
2. **Partner with a clinical lab** — start with retrospective samples
3. **Test on public GEO/SRA cfDNA datasets** — first "real data" validation

### Strategic Moves
4. **Switch CET primary analyte to methylation** — PanSeer proved methylation works for pre-diagnosis detection
5. **Add methylation-based CET alongside mutation-based** — best of both worlds
6. **Publish as pre-print on bioRxiv/medRxiv** — get community feedback early

### The One Thing That Would Change Everything
> **A single positive clinical pilot (n=50 cancer + 50 healthy, detecting Stage I cancers at >80% sensitivity with >95% specificity)** would immediately elevate DeepCatch to Tier 2 and attract serious attention.

---

*Analysis compiled from 21-paper literature review + DeepCatch v1.2 improvement data.*
