# 08 — HLA per-allele layer

HLA typing is available for the full cohort. This promotes a whole class of hypotheses from postdoc-scale to lead-set and adds a per-allele resolution tier to the platform.

---

## Why this is a resolution step up

Every other feature in the platform treats HLA genes as three loci (A, B, C). With typing, each locus becomes **two haplotypes**, and outcome variables can be quantified per-allele:

- `expr[HLA-A*02:01]` vs `expr[HLA-A*24:02]` separately, not `expr[HLA-A]` pooled.
- `β_prom[HLA-A allele 1]` vs `β_prom[HLA-A allele 2]` if haplotype-resolved bisulfite reads are available (and frequently they are, because HM450/EPIC probes are polymorphic in MHC region).
- `CRML` computed per-allele, because SNVs in a regulatory element are haplotype-phased by the sample's typing.

This means the platform's entire per-element / per-variant attribution machinery runs at **per-allele resolution on HLA loci**. Protein-coding driver analysis has long had per-residue resolution. HLA regulatory analysis at per-allele × per-element resolution at cohort scale is rare to absent in the public literature.

---

## Required ingredients

| Ingredient | Status | Notes |
|---|---|---|
| Class-I HLA typing (A, B, C) | **Available, full cohort** | User-provided |
| Class-II typing (DR, DQ, DP) | Check availability | Optional but valuable for CIITA / MHC-II claims (H25, H49) |
| Per-allele RNA quantification | To build | `arcasHLA` or `HLApers`-style allele-aware alignment from TCGA BAMs, or `PHLAT`-derived per-allele counts |
| LOHHLA (genomic HLA copy / allele loss) | To build or check | Needed to separate regulatory ASE from genomic LOH; McGranahan lab's tool, re-implementable |
| Haplotype-phased methylation at HLA | Approximate via overlapping-SNP probe filtering | HM450/EPIC probes with polymorphic sites must be flagged (Atlas table 6a already has `contains_snp`) |

---

## New L1 feature builders

Added to the contract in `03_architecture.md` §L1:

| Builder | Returns |
|---|---|
| `build_hla_allele_typing()` | per-sample 6-allele genotype (A1/A2, B1/B2, C1/C2), normalized to IMGT nomenclature |
| `build_hla_allele_expression()` | per-sample × per-allele log2(TPM+1) at HLA-A/B/C; total loci = 6 columns per sample |
| `build_hla_allele_ase()` | per-sample × per-locus allele imbalance score and flag (significant ASE present) |
| `build_hla_lohhla()` | per-sample × per-allele genomic copy status, `loh_flag` per allele |
| `build_hla_allele_methylation()` | per-sample × per-allele promoter β (where haplotype-phasing possible; else NA flag) |
| `build_hla_allele_crml()` | per-sample × per-allele CRML using haplotype-phased SNVs at HLA enhancers |
| `build_neoantigen_load_per_allele()` | per-sample × per-allele predicted neoantigen count (from existing MHC affinity prediction ingestion; requires upstream SNV → peptide → netMHCpan pipeline) |

## New derived phenotypes

| Phenotype | Definition | Hypothesis |
|---|---|---|
| `ASE_reg` | Allele expression drop **without** genomic LOH → regulatory loss | H46, H63 |
| `ASE_LOH` | Allele expression drop **with** genomic LOH → genomic loss | H46 (reference class) |
| `editing_residual_per_allele` | Per-allele neoantigen load vs presentation competence; residual | H65 |

---

## Hypotheses (Family M)

### H46 (promoted) — HLA allele-specific expression loss is regulatory, not genomic
- **Statement**: A substantial fraction of per-allele HLA expression loss in BRCA occurs without genomic LOH, i.e., is regulatory (allele-specific methylation, haplotype-phased enhancer damage, 3'UTR/miR asymmetry).
- **Features**: `ASE_reg` vs `ASE_LOH` classification; enhancer β per allele; per-allele CRML.
- **Tier**: T2 | **Priority**: 1 | **Lead-set**: yes | **Claim**: C1, C2

### H63 ✱ — Regulatory ASE is mechanistically attributable to per-allele CRML + per-allele β
- **Statement**: Within the `ASE_reg` class, which allele drops is predictable from per-allele CRML and per-allele promoter/enhancer β on that haplotype.
- **Features**: Haplotype-resolved CRML + β on the deflated vs retained allele.
- **Tier**: T2 | **Priority**: 1 | **Lead-set**: yes | **Claim**: C1, C2

### H64 — HLA supertype × APM structural-integrity interaction predicts ICB response
- **Statement**: Chowell supertype classification × four-class hot/cold decomposition (H59) better predicts ICB response than either alone; specific supertypes are more vulnerable to specific mechanism classes.
- **Features**: HLA supertype (Sidney) + four-class label + ICB response.
- **Tier**: T2 | **Priority**: 2 | **Claim**: C1, C5

### H65 — Per-allele immunoediting signature at element resolution
- **Statement**: Per-allele neoantigen load is anti-correlated with per-allele presentation competence within the same tumor — a per-allele-resolved immunoediting fingerprint, attributable to specific CRML / methylation lesions on the more-presenting allele.
- **Features**: `neoantigen_load[allele]`, `editing_residual_per_allele`, per-allele CRML + β.
- **Method**: mixed model with patient-level random effect; both alleles in-tumor as paired observations.
- **Tier**: T2 | **Priority**: 2 | **Claim**: C1, C2, C5

---

## Integration with other families

- Family G (methylation): `COMB_HLA` extends to per-allele PC1 where probe-level haplotype resolution allows.
- Family J (within-panel): H38 bottleneck and H39 paralog compensation refined by per-allele expression.
- Family N (TNBC hot/cold, next file): four-class decomposition × per-allele ASE gives a finer sub-classification; "hot-discordant" tumors may turn out to be driven by mono-allelic presentation loss with the retained allele's neoantigens well-presented.

---

## Caveats, honestly

- Haplotype-phased methylation at HLA is partial: EPIC/HM450 probes with polymorphic sites are unreliable; many MHC-region probes must be filtered. The `contains_snp` flag in the methylation probe reference (Atlas 6a) already supports this.
- LOHHLA inference from BAMs is not trivial and depends on coverage. Budget 1 month for a pipeline stand-up or reuse McGranahan's code.
- Per-allele neoantigen prediction (H65) depends on an SNV → peptide → netMHCpan pipeline. Check whether this is already in the analysis/ layer or needs to be added.

---

## Paper arc

- **Methods paper** (standalone, short): "Regulatory loss distinct from genomic loss at HLA in BRCA", built on H46 + H63 + 13-HiChIP calibration. Genome Biology / Nat Commun band.
- **Mechanism paper** (integrated with Aim 2): H65 immunoediting signature, threaded through the four-class hot/cold TNBC story.
