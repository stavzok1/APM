# CPTAC protein validation of the miRNA-Hallmark pressure spine — session summary

**Session date:** 2026-06-23/24 · **Registry:** MH-33 → MH-36 · **Modules:**
`cptac_validation.py`, `cptac_resid_survivors.py`, `cptac_orphan_discovery.py` ·
**Data:** `scripts/cptac/{fetch,ingest}_cptac2_breast_mirna.py`

> One-line: the subproject's miRNA-pressure spine is validated **at the protein level** in an
> **independent** breast cohort, the repression is shown to be **predominantly mRNA-mediated**,
> the small genuinely-translational residual is carried by **known, subtype-structured oncomiR
> edges (incl. the basal miR-17~92→PTEN hub)**, and an orphan-edge scan nominates **93 novel
> protein-level candidate edges** with the method self-validated by recovering miR-29→collagen.

---

## 1. What was done (in order)

| # | Analysis | Module | Cohort |
|---|----------|--------|--------|
| 1 | Two CPTAC BRCA cohorts disambiguated (TCGA-105 iTRAQ vs prospective-122) | — | — |
| 2 | Pressure → protein validation, 4 pressure variants × 3 layers | `cptac_validation` | TCGA-105 (same patients) |
| 3 | Built CPTAC-2 prospective miRNA-seq from GDC (104 cases, open-access) → arm matrix | `scripts/cptac/*` | prospective |
| 4 | Same validation, **self-contained** in the independent prospective cohort | `cptac_validation --cohort prospective` | prospective (101) |
| 5 | `protein_resid` layer added — protein residualized on its own mRNA (slope-free "beyond mRNA") | `cptac_validation` | both |
| 6 | Single-edge resolution of the residual survivors (which miRNA / subtype / literature) | `cptac_resid_survivors` | prospective |
| 7 | Orphan-edge discovery scan (TargetScan ∪ ENCORI, miRTarBase-absent) → protein anti-corr | `cptac_orphan_discovery` | both |

### The two CPTACs (the framing that made this work)
- **CPTAC TCGA-105** (Mertins 2016, iTRAQ): the *same* 105 TCGA patients re-assayed for proteome. Deterministic join.
- **CPTAC prospective-122** (Krug 2020, TMT): *independent* patients, NOT TCGA. Its own miRNA-seq lives on **GDC CPTAC-2** (open-access, 104 breast cases; `01BR001` ↔ LinkedOmics `X01BR001`) — so it becomes a **self-contained** miRNA+RNA+protein cohort with no cross-cohort matching.

---

## 2. Prominent results

### 2.1 The spine validates at protein — but only in the independent cohort
| | TCGA-105 (same patients, iTRAQ) | Prospective (independent, TMT) |
|---|---|---|
| Genes FDR-neg (protein ~ pressure) | **2 / 1006** | **53 / 1143** |
| Hallmarks FDR-neg | MYOGENESIS tail | **6 / 50**, EMT top (−0.32, q=0.011) |
| `protein_resid` survivors | **0** | **44** |

The noisier same-patient iTRAQ proteome cannot see the coupling; the cleaner prospective TMT proteome (minimized ischemic time) recovers it. **ZEB1** — the canonical miR-200 target — leads at ρ=−0.55 (q=1.4×10⁻⁵).

### 2.2 The repression is predominantly mRNA-mediated (your mRNA snapshot model, confirmed)
Residualizing protein on its own mRNA (`protein_resid`) collapses most coupling — e.g. **ZEB1 −0.55 → −0.24 (NS)** — leaving only **11/1132** genes with a genuine protein-beyond-mRNA component. miRNA pressure acts mainly via mRNA destabilization (already captured by RNA-seq); the translational residual is small. Consistent with Bartel/Eichhorn canonical mammalian miRNA biology.

### 2.3 The translational survivors are known, subtype-structured edges — incl. the basal hub
Single-edge resolution: **17/23 survivors have an FDR driver arm; ALL 17 are literature-validated** (15 miRTarBase functional-MTI, 2 breast-PMID). **PTEN's residual is carried by the miR-17~92 cluster** (miR-17-5p ρ=−0.34 q=0.008; miR-19a/18a/20a) — exactly the basal hub the subproject identified (MH-18/20/21), now at the protein layer. Subtype split: **8 Basal** (PTEN, KRT8, ZYX, HMGN2, SLPI) / **7 LumA** (DDAH1, PEBP1, ARHGDIA) / 2 LumB.

### 2.4 Orphan discovery: 93 novel protein-level candidate edges
Scanning 20,663 orphan edges (TargetScan ∪ ENCORI, miRTarBase-functional-MTI-absent): **539 protein / 182 translational / 93 translational+TCGA-replicated** candidates (**76 fully uncurated, 29 dual-prior**). **Method self-validated** by de-novo recovery of the **miR-29→collagen/ECM axis** (miR-29a→COL11A1/COL6A2/BMP1). Top fully-uncurated dual-prior hits: miR-30e→MINPP1/VIM, miR-29a→COL11A1, miR-17→NEDD4L, miR-15b→AP2B1. **Wet-lab nominations, not validated edges.**

---

## 3. Actual anti-correlations: TCGA RNA vs independent CPTAC protein

Per-gene partial Spearman (pressure vs target). **TCGA RNA** = `target_combined_anticorr` (gated, CPE+HRD-adjusted, n≈1077). **CPTAC protein / resid** = prospective `cptac_validation` (highev|gated, purity+CIN-adjusted, n≈101). The genes split into **three mechanistic classes**:

### Class A — mRNA-dominant (strong RNA + protein, residual ~0): repression visible in RNA, propagates to protein
| gene | TCGA RNA ρ | CPTAC protein ρ | protein_resid ρ | note |
|------|-----------:|----------------:|----------------:|------|
| ZEB1 | −0.24 | **−0.55** | −0.24 (NS) | canonical miR-200 target; protein amplifies, residual gone |
| CALD1 | −0.33 | −0.51 | −0.23 (NS) | |
| COL5A2 | −0.40 | −0.32 | −0.24 (NS) | |
| MMP2 | −0.34 | −0.28 | −0.10 (NS) | |
| FBN1 | −0.32 | −0.19 | +0.01 (NS) | ECM; purely mRNA-level |
| TGFB3 | −0.35 | −0.20 | +0.04 (NS) | ECM; purely mRNA-level |

### Class B — protein/translational-dominant (weak RNA, strong protein **and** residual): regulation RNA misses
| gene | TCGA RNA ρ | CPTAC protein ρ | protein_resid ρ | driver miRNA (subtype) |
|------|-----------:|----------------:|----------------:|------------------------|
| **PTEN** | −0.10 | −0.30 | **−0.38 (q=0.011)** | **miR-17-5p (Basal)** — the hub |
| DDAH1 | −0.07 (NS) | −0.50 | −0.40 (q=0.007) | miR-210-3p (LumA) |
| ARHGDIA | −0.05 (NS) | −0.36 | −0.45 (q=0.002) | miR-151a-5p (LumA) |
| PEBP1 (RKIP) | −0.11 | −0.45 | −0.43 (q=0.002) | miR-155-5p (LumA) |
| EPHB2 | +0.10 | −0.34 | −0.40 (q=0.046) | miR-128-3p (LumA) |
| MTCH2 | +0.09 | −0.26 | −0.35 (q=0.035) | miR-135b-5p (Basal) |

### Class C — mixed (strong RNA + protein, residual borderline)
| gene | TCGA RNA ρ | CPTAC protein ρ | protein_resid ρ | driver miRNA (subtype) |
|------|-----------:|----------------:|----------------:|------------------------|
| KRT8 | −0.22 | −0.39 | −0.43 (q=0.002) | miR-186-5p (Basal) |
| BMP1 | −0.26 | −0.47 | −0.40 (q=0.007) | miR-194-5p (Basal) |
| ZYX | −0.31 | −0.35 | −0.30 (q=0.11) | miR-16-5p (Basal) |
| HMGN2 | −0.10 | −0.33 | −0.33 (q=0.054) | miR-23a-3p (Basal) |
| MINPP1 | −0.13 | −0.31 | −0.35 (q=0.037) | miR-30b-5p (Basal) |
| SLPI | +0.06 (NS) | −0.15 (NS) | −0.33 (q=0.051) | miR-146a-5p (Basal) |

**Reading:** Class A confirms your "mRNA is the intermediate snapshot" model — the effect is already in the RNA and rides through to protein, with no residual. Class B is the payoff: **PTEN, DDAH1, ARHGDIA, PEBP1 are barely visible in TCGA RNA but strongly repressed at protein with a genuine translational residual** — regulation the mRNA layer alone would have missed. PTEN is the standout: weak in RNA (−0.10), genuine at protein-residual (−0.38) via the basal miR-17~92 hub.

---

## 4. Comparison with literature

| Result | Literature standing |
|--------|---------------------|
| ZEB1 ← miR-200 at protein | **Confirms canonical** miR-200/ZEB1 EMT double-negative feedback — here at endogenous tumor protein |
| PTEN ← miR-17~92 (Basal, protein-residual) | **Confirms** the textbook miR-17~92→PTEN axis; adds endogenous-protein, in-tumor, translational-specific evidence |
| PEBP1(RKIP) ← miR-155 | **Confirms** known miR-155→RKIP repression |
| Repression mostly mRNA-mediated | **Consistent** with Eichhorn/Bartel 2014 (mRNA destabilization dominates steady-state miRNA effect) |
| 17 driver edges: all single-study, 15/17 non-breast, 4/17 reporter-only | Your contribution = **evidence-tier upgrade** (endogenous protein, breast tumor, translational, subtype-resolved), not new edges |
| miR-29 → COL11A1/COL6A2/BMP1 (orphan, de novo) | **Recovers** the established miR-29-collagen/ECM axis at uncurated specific targets — validates the orphan method |
| Prospective TMT > same-patient iTRAQ | **Consistent** with the Mertins(iTRAQ)→Krug(TMT) quality progression (reduced ischemic time) |

**Honest ceilings:** protein anti-correlation ≠ direct edge (CLIP/luciferase needed); per-subtype n is small (descriptive); "literature-validated" = the repo's curated miRTarBase/breast-PMID layer, not a fresh per-edge search.

---

## 5. Registry / outputs

| ID | Claim | Output |
|----|-------|--------|
| MH-33 | TCGA-105 (same patients) protein null — reproduces parent D63 | `cptac_validation/` |
| MH-34 | Independent prospective cohort validates at protein (ZEB1/EMT); mostly mRNA-mediated | `cptac_validation/prospective/` |
| MH-35 | Translational survivors = known subtype-structured arms; PTEN←miR-17~92 hub at protein | `cptac_validation/prospective/resid_survivors/` |
| MH-36 | Orphan scan → 93 novel candidate edges; miR-29→collagen recovered de novo | `cptac_validation/orphan_discovery/` |

**Reproduce:** `scripts/cptac/fetch_cptac2_breast_mirna.py` → `ingest_cptac2_mirna_arms.py` →
`python -m mirna_hallmark.cptac_validation --cohort {tcga105,prospective}` →
`cptac_resid_survivors` → `cptac_orphan_discovery`.
