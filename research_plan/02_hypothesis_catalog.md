# 02 — Hypothesis catalog (H1–H65)

Every hypothesis in one place with a consistent schema. This is the file the agent (or a new collaborator) reads to pick the next test.

## Schema

Each hypothesis has:

- **ID** (H1–H65)
- **Title** — one-liner
- **Biological statement** — mechanism being claimed
- **Atlas dependencies** — which tables from `pipeline/md/DATA_STRUCTURE_ATLAS.md` it needs
- **Derived features** — L1/L2 quantities it consumes (see `03_architecture.md`)
- **Automation tier**: T1 (fully automatable, minutes–hours), T2 (semi-automatable, days), T3 (design-heavy, weeks)
- **Lead-set flag** ✱ — pre-registered for A1/A2 committee defense (tight subset)
- **Priority** (1 = must-run, 2 = plan-run, 3 = exploratory)
- **Claim(s) it supports** — C1..C5

Lead pre-registered set (15 hypotheses): **H1, H4, H9, H12, H14, H18, H19, H38, H41, H46, H51, H53, H56, H59, H63**.

Cohort facts: HLA typing (class-I A/B/C) is available for all samples, enabling per-allele outcome variables across Family M. 13 per-tumor HiChIP samples calibrate the conservation prior (`04_conservation_prior.md`). Wet-lab validation (`07_wet_lab.md`) is available.

Automation tiers and ranking live in `06_automation_plan.md`.

---

## Family A — Cis-regulatory mutational load and element-resolved causal claims (H1–H3, H9, H16)

### H1 ✱ — CRML at IRF1/STAT1 enhancers predicts APM expression loss beyond methyl + CNV
- **Statement**: Σ |Δ motif score| × 1[TF bound by ChIP] at ChIP-bound, ABC-strong enhancers of each APM gene predicts expression loss after conditioning on CNV, promoter/enhancer β, purity, PAM50, IFN-γ.
- **Atlas**: Element Focus (`gene_links`, `chip_hits`), SNV (`motif_hits`, `cCRE_hits`), SV (`elem_hits.motif_hits`, `flank_motif_hits`, `chip_hits`), Methylation, CNV, RNA.
- **Features**: `CRML_snv_g`, `CRML_sv_g`, promoter β, enhancer β, CNV dosage, IFN-γ ssGSEA.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C1, C2

### H2 — CRML explains residually cold tumors
- **Statement**: In samples with IFN-γ low but no structural APM lesion, elevated CRML at STAT1/IRF1-bound sites drives the cold phenotype.
- **Atlas**: SNV, SV, Element Focus, RNA.
- **Features**: `CRML`, APM_intrinsic (IFN-γ residualized).
- **Tier**: T1 | **Priority**: 2 | **Claim**: C2

### H3 — Subtype-skewed coding vs cis damage
- **Statement**: Basal-like carries more coding APM damage (nonsense/frameshift in `has_*_canonical`); luminal carries more cis-regulatory damage (CRML).
- **Features**: `CRML_snv`, `CRML_sv`, PAM50, canonical coding flags.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H9 ✱ — Conservation-weighted boundary-disruption (BCD) predicts expression where raw SV count does not
- **Statement**: SV disruption of TAD boundaries weighted by their across-24-biosample conservation (and calibrated to 13 per-tumor HiChIP where available) predicts gene dysregulation; raw SV counts do not.
- **Atlas**: Element Focus (`TAD_domains`×24), SV, 13 per-tumor HiChIP.
- **Features**: `BCD_g` with conservation weights.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C2

### H16 — Serial mediation: CRML → CXCL9/10/11 → activated-CD8 → CYT → survival
- **Statement**: A four-link chain from cis-regulatory IRF1-driven CXCL expression to infiltration to cytolytic activity to outcome has significant cumulative indirect effect.
- **Atlas**: SNV, SV, RNA, Thorsson, Clinical.
- **Features**: `CRML`, `activated_CD8`, `CYT_score`, `OS`.
- **Tier**: T3 | **Priority**: 2 | **Claim**: C2, C5

---

## Family B — Proteostatic visibility (protein-layer claims; H4–H6, H15)

### H4 ✱ — mRNA-competent, protein-deficient APM class exists
- **Statement**: A subset of tumors show `HLA-A`/`B2M` mRNA normal but protein low (RPPA `protein_rna_discordance`), driven by 3'UTR SNVs, miR pressure, and/or NMD substrates.
- **Atlas**: RPPA, RNA, SNV, miRTarBase.
- **Features**: `PRD_g`, `MIR_g`, 3'UTR variants.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C3

### H5 — IFN+ / pSTAT1− signaling-refractory class
- **Statement**: Samples with high IFN-γ signature but low `p-STAT1/STAT1` ratio are signaling-refractory due to `SOCS1/3` up, `PTPN2` amplification, or JAK1/2 damaging variants.
- **Atlas**: RPPA, RNA, CNV, SNV.
- **Features**: `p-STAT1/STAT1` ratio, IFN-γ ssGSEA, SOCS1/3 expression, CNV/SNV of JAK1/2, PTPN2.
- **Tier**: T2 | **Priority**: 2 | **Claim**: C3

### H6 — Protein APC outperforms mRNA APM for ICB prediction
- **Statement**: In an external ICB cohort, RPPA-estimated APM capacity stratifies response better than mRNA APM.
- **Atlas**: RPPA, Clinical, external ICB cohort.
- **Tier**: T3 | **Priority**: 2 (pending cohort access) | **Claim**: C3, C5

### H15 — Proteasome-primed / MHC-lost subgroup
- **Statement**: Tumors with high `PSMB8/9/10` but low `HLA-A/B/C` define a proteasome-primed / MHC-lost subgroup with distinct therapeutic handles (oncolytic, MHC-restoring).
- **Features**: `PSMB` vs `HLA` expression and protein.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1, C3

---

## Family C — Enhancer redundancy, robustness, and buffering (H7, H8)

### H7 — R_g × burden interaction (buffering)
- **Statement**: For equal CNV/SV/meth burden, genes with high enhancer redundancy R_g show smaller expression change; interaction term R_g × burden has negative sign.
- **Atlas**: Element Focus (`gene_links`: SCREEN, ABC, HiChIP conservation), ATAC, CNV, SNV.
- **Features**: `REDUND_g`, aggregate burden.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H8 — R_g rank explains B2M fragility
- **Statement**: `B2M` has few strong enhancers → low R_g → observed fragility in ICB-resistant tumors; `HLA-A/B/C` high R_g → robustness.
- **Features**: `REDUND_g`.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

---

## Family D — Deconvolution-state × structural intactness (H12, H13, H14, H50)

### H12 ✱ — Activated CD8 drives APM only when JAK–STAT structurally intact (moderation)
- **Statement**: `APM_g ~ activated_CD8 × JSI + covariates`, where `JSI` integrates CNV, SV, damaging SNV, and promoter methylation of JAK1/2 and STAT1. Interaction is positive.
- **Atlas**: Thorsson, CNV, SV, SNV, Methylation, RNA.
- **Features**: `AST_CD8a`, `JSI`.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C2

### H13 — Activated NK tracks NK-ligand minus shedding, not mRNA alone
- **Statement**: Thorsson activated-NK fraction correlates with `(MICA/B mRNA) − shedding index (ADAM10/17)` better than mRNA alone.
- **Features**: `DAV` (dual-axis visibility), shedding score.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1, C3

### H14 ✱ — Four-quadrant MHC-I / NK-ligand escape map
- **Statement**: Tumors segregate into MHCI±/NKL± quadrants with distinct infiltration patterns and survival.
- **Atlas**: RNA, RPPA, Clinical, Thorsson.
- **Features**: `QUAD` (from `DAV`).
- **Tier**: T1 | **Priority**: 1 | **Claim**: C1, C5

### H50 — ERAP1/2 × HLA-E × NKG2A inhibitory axis
- **Statement**: Low `ERAP1/2` expression disrupts VL9 presentation on HLA-E → NKG2A inhibition lifted → NK-cell infiltration up, even when MHC-I expression is high.
- **Features**: ERAP1/2 expression, HLA-E expression, activated NK.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

---

## Family E — RPPA archetype / signaling-block phenotypes (H10, H11)

### H10 — Canonical escape archetypes from RPPA block vector
- **Statement**: Cluster samples on the 7-bit `signaling_blocks` vector → canonical archetypes; each has a distinct genetic/epigenetic signature.
- **Atlas**: RPPA.
- **Features**: `BLK` (7-bit vector), `ARCH` (cluster label).
- **Tier**: T2 | **Priority**: 2 | **Claim**: C3

### H11 — RPPA archetype vs Thorsson C1–C6 disagreement is informative
- **Statement**: Where archetype and Thorsson disagree, the disagreement predicts non-response better than either alone.
- **Tier**: T2 | **Priority**: 3 | **Claim**: C3

---

## Family F — Platform and prognosis (H17, H18)

### H17 — Mechanism fingerprint discriminates pathways
- **Statement**: `fingerprint(gene_set)` run on APM, DDR, EMT yields distinct dominant-mechanism profiles.
- **Tier**: T1–T2 (needs second/third panel) | **Priority**: 2 | **Claim**: C4

### H18 ✱ — Intrinsic Visibility Index predicts OS beyond stage + PAM50
- **Statement**: IVI composite predicts OS/PFS in TCGA-BRCA adjusted for stage, grade, PAM50, age; replicates in METABRIC.
- **Atlas**: all regulatory + Clinical.
- **Features**: `IVI`.
- **Tier**: T2 | **Priority**: 1 | **Claim**: C5

---

## Family G — Methylation-centric (H19–H25)

### H19 ✱ — Coordinated HLA-locus co-methylation block
- **Statement**: Promoter β of `HLA-A/B/C`, `HCP5`, `HCG18` co-vary → a single PC1 captures more APM variance than per-gene β; predicted by DNMT3A/B expression and local CNV.
- **Atlas**: Methylation (per-sample gene + cCRE), RNA.
- **Features**: `COMB_HLA` = PC1 of locus β.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C1

### H20 — Promoter β vs enhancer β differential effect
- **Statement**: Enhancer β controls inducibility (`β_enh × IFN-γ` interaction), promoter β controls baseline.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H21 — PRC2-silenced APM subset
- **Statement**: Subset of tumors carry bivalent chromatin → APM silenced despite normal CNV and normal promoter β; `EZH2/SUZ12/EED` expression elevated.
- **Features**: PRC2 core expression, APM residuals.
- **Tier**: T2 | **Priority**: 2 | **Claim**: C1, C3

### H22 — STAT1 methylation cascade
- **Statement**: Serial mediation `β_prom(STAT1) → STAT1 expr → IRF1/NLRC5 expr → APM expr`.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C2

### H23 — CpG-island translocations via SV
- **Statement**: SVs that relocate CpG-island regions into/out of enhancers cause locus-scale methylation reprogramming.
- **Atlas**: SV, Methylation probes with `overlapping_ccre_ids`, `overlapping_atac_ids`.
- **Tier**: T2 | **Priority**: 3 | **Claim**: C2

### H24 — mQTL interaction: HLA genotype × methylation
- **Statement**: HLA-A/B/C alleles associate with differential methylation at MHC-I promoters and with ASE loss.
- **Tier**: T3 (needs HLA typing) | **Priority**: 3 | **Claim**: C1

### H25 — CIITA methylation as MHC-II biomarker
- **Statement**: `β_prom(CIITA)` correlates with activated-CD4 Thorsson fraction and survival.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

---

## Family H — miRNA-centric (H26–H31)

### H26 — miR-148a / miR-152 → HLA-G
- **Statement**: Evidence-weighted miR-148a/152 pressure correlates inversely with HLA-G after controlling for CNV and β; segregates by Thorsson subtype.
- **Atlas**: miRTarBase summary (`evidence_score`), RNA.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1, C3

### H27 — miR-125a/b → B2M/MHC-I coordinated repression
- **Statement**: High miR-125a/b co-occurs with B2M/MHC-I drop independent of CNV and β.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C3

### H28 — miR-155 → SOCS1 disinhibits JAK-STAT (hyper-responsive APM)
- **Statement**: Interaction `IFN-γ × miR-155` predicts APM higher than IFN-γ alone because SOCS1 is suppressed.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

### H29 — Network-level miR master regulators (data-driven)
- **Statement**: Enrichment of predicted targets (weighted by `evidence_score`) among APM residuals identifies data-driven master-regulator miRs.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H30 — miR-346 → TAP1
- **Statement**: High miR-346 predicts TAP1 loss independent of CNV/β.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C3

### H31 — Coordinated miR-pressure phenotype
- **Statement**: A subset of tumors co-elevate multiple APM-repressive miRs (148a, 9, 125, 346, 34a); `miR_panel_pressure` score defines a distinct escape class.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C3

---

## Family I — lncRNA-centric (H32–H37)

### H32 — HCP5 ceRNA → PD-L2 axis
- **Statement**: HCP5 sponges miR-106b-5p / miR-17-5p → derepresses `PDCD1LG2`. HCP5 expression correlates with PD-L2 beyond CNV/methylation.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1, C3

### H33 — HCG18 ceRNA → PD-L1 / FGFR1 axis
- **Statement**: HCG18 expression correlates with `CD274` and oncogenic targets; mediated by miRNA sponging.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

### H34 — HLA-locus lncRNAs as CTCF / boundary maintainers
- **Statement**: HCP5/HCG18 expression correlates with MHC-locus TAD boundary strength (`TAD_boundary_overlaps`, 24 biosamples + 13 tumor HiChIP). Loss predicts topology disruption.
- **Tier**: T2 | **Priority**: 2 | **Claim**: C1, C2

### H35 — lncRNA–gene TAD colocalization predicts cis effect
- **Statement**: For each lncRNA–APM pair within the same TAD, lncRNA expression correlates with APM more than size-matched trans-lncRNAs.
- **Atlas**: `lncRNAs_within_1000kb`, `TAD_domains`.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H36 — lncRNA-mediated DNMT/PRC2 recruitment (mechanism for H19)
- **Statement**: Chromatin-modifying lncRNAs elevated in samples with coordinated HLA hypermethylation — the molecular mechanism for H19.
- **Tier**: T2 | **Priority**: 3 | **Claim**: C1, C3

### H37 — lncRNA expression as locus-competence proxy for samples lacking ChIP
- **Statement**: lncRNA expression at a locus predicts local transcriptional activity; validate by correlating with ChIP strength at nearby enhancers in same-subtype biosamples.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C4

---

## Family J — Within-panel gene-gene relations (H38–H50)

### H38 ✱ — Bottleneck principle (min-rule)
- **Statement**: APM output is capped by `min(JAK1, JAK2, STAT1, IRF1, NLRC5)` better than by their sum/mean.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C1

### H39 — Paralog compensation
- **Statement**: Loss of one paralog (e.g. HLA-A) predicts upregulation of the others (HLA-B/C) — regression of residuals has a negative sign when compensation exists.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H40 — STAT1/STAT3 antagonist ratio
- **Statement**: `log(STAT1/STAT3)` or `p-STAT1/p-STAT3` outperforms either alone in predicting APM and ICB response.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H41 ✱ — IRF1/IRF2 antagonist ratio
- **Statement**: (Lit-backed, PMC6761035.) `IRF1/IRF2` ratio predicts APM expression, TAP2/ERAP1 regulation, and PD-L1 expression.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C1

### H42 — Negative-feedback loop intactness (SOCS1 residual)
- **Statement**: Define `feedback tension` = residual of SOCS1 given predicted IFN response. Leaky feedback → hyper-responsive APM; over-damped → dampened.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H43 — EP300/CREBBP coactivator rate-limit
- **Statement**: Interaction `STAT1 × (EP300 + CREBBP)` on APM induction; coactivator-depleted tumors fail to induce APM despite STAT1/IRF1 competence.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

### H44 — Pathway coherence as sample phenotype
- **Statement**: Per-sample coherence = `var(PC1) / total_var` across 66 APM genes. Low coherence = decoherent pathway, a distinct phenotype from low mean.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H45 — Transmission-bottleneck detection
- **Statement**: Per-sample "where the signal breaks" map by regressing residuals of downstream on immediate upstream.
- **Tier**: T2 | **Priority**: 3 | **Claim**: C2

### H46 ✱ — HLA allele-specific expression (ASE) loss is regulatory, not only genomic
- **Statement**: (Lit: PMC9912643; McGranahan 2017 LOHHLA.) A substantial fraction of HLA-A/B/C per-allele expression loss occurs **without** genomic LOH — driven by allele-specific methylation, haplotype-phased enhancer damage, and 3'UTR/miR asymmetry. Predict per-allele imbalance from local enhancer β, per-allele CRML, and SV, over and above CNV.
- **Atlas**: per-sample HLA class-I typing (available), RNA (allele-aware), CNV, SV, Methyl, SNV. LOHHLA to be layered on.
- **Tier**: T2 | **Priority**: 1 | **Claim**: C1, C2.
- See Family M below for the detailed allele-resolved sub-catalog.

### H47 — Peptide-pipeline coherence (proteasome ↔ transport ↔ loading)
- **Statement**: Composite coherence across the three APM stations outperforms any single-station score.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H48 — cGAS–STING axis intactness (parallel)
- **Statement**: Tumors with intact cGAS–STING but blocked IFN-γ retain APM induction via type-I IFN; two-axis integrity score improves APM prediction.
- **Tier**: T2 (requires Tier-1 panel extension) | **Priority**: 2 | **Claim**: C1

### H49 — NLRC5 vs CIITA coordination (MHC-I vs MHC-II)
- **Statement**: Four-quadrant classification by NLRC5/CIITA expression predicts visibility phenotype for CD8 vs CD4.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

---

## Family K — Conservation-as-prior (H51–H55) — the methodological backbone

All details in `04_conservation_prior.md`. The 13 per-tumor HiChIP samples anchor calibration.

### H51 ✱ — Subtype-consensus HiChIP edges predict expression, calibrated by 13 tumor HiChIP
- **Statement**: Same-subtype cell-line HiChIP edge strength predicts per-tumor edge strength in the 13 HiChIP samples; calibrated prior then applies to all other BRCA tumors.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C2, C4

### H52 — Subtype-consensus ATAC predicts expression where per-tumor ATAC missing
- **Statement**: Consensus ATAC (from the ~74-sample BRCA ATAC cohort) + per-tumor RNA improves prediction over either alone.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C4

### H53 ✱ — "Deviation from subtype prior" is a better feature than absolute
- **Statement**: For each feature with a subtype prior, residual (observed − prior) is more predictive of expression and outcome than raw value.
- **Tier**: T1 | **Priority**: 1 | **Claim**: C4

### H54 — Core APM enhancers (defined by biosample conservation) have largest effects
- **Statement**: Enhancers in ≥k breast biosamples at strong level have ×N the per-sample expression effect of opportunistic enhancers.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

### H55 — Subtype-dependent TAD boundary effect magnitudes
- **Statement**: SV disruption of subtype-conserved boundaries has larger per-sample effects than non-conserved; effect sizes differ by PAM50.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C2

---

## Family L — Innate sensing extension (H56–H58)

### H56 ✱ — Dual-axis IFN intactness predicts APM induction
- **Statement**: Composite `IFN-γ intactness × cGAS–STING intactness` predicts APM; tumors with only one axis retain partial induction.
- **Tier**: T2 (requires Tier-1 panel extension) | **Priority**: 1 | **Claim**: C1, C2

### H57 — GAS vs ISRE mode preference
- **Statement**: Per-tumor preference for STAT1 homodimer (GAS) vs ISGF3 (ISRE) mode, determined by STAT2/IRF9 availability.
- **Tier**: T2 | **Priority**: 3 | **Claim**: C1

### H58 — Type-I vs type-II IFN receptor dosage ratio
- **Statement**: `log(IFNGR1·IFNGR2 / IFNAR1·IFNAR2)` is a tumor-intrinsic bias indicating which IFN therapy will work.
- **Tier**: T1 | **Priority**: 3 | **Claim**: C1

---

## Family M — HLA per-allele layer (H46, H63–H65)

All details in `08_hla_allele_layer.md`. Requires per-sample HLA class-I typing (**available**) plus an allele-aware RNA quantification builder and LOHHLA-style genomic copy inference.

H46 (above) is the entry point; H63–H65 extend it. This family is the single biggest upgrade of the plan since HLA typing arrived.

### H63 ✱ — Regulatory ASE is mechanistically attributable to per-allele CRML + per-allele β
- **Statement**: Within the `ASE_reg` class (allele expression drop without genomic LOH), which allele drops is predictable from per-allele CRML and per-allele promoter/enhancer β on that haplotype.
- **Features**: Haplotype-resolved CRML + β + SV on the deflated vs retained allele.
- **Tier**: T2 | **Priority**: 1 | **Lead-set**: yes | **Claim**: C1, C2

### H64 — HLA supertype × structural-integrity interaction on ICB response
- **Statement**: Chowell supertype × four-class hot/cold decomposition (H59) better predicts ICB response than either alone; specific supertypes are more vulnerable to specific mechanism classes.
- **Features**: Sidney supertype + four-class label + ICB response (external cohort).
- **Tier**: T2 | **Priority**: 2 | **Claim**: C1, C5

### H65 — Per-allele immunoediting signature at element resolution
- **Statement**: Per-allele neoantigen load is anti-correlated with per-allele presentation competence within the same tumor — a per-allele immunoediting fingerprint, mechanistically attributable to specific CRML / β / SV lesions on the more-presenting allele.
- **Features**: `neoantigen_load[allele]`, `editing_residual_per_allele`, per-allele CRML + β + SV.
- **Method**: mixed model, patient random effect, both alleles as paired observations.
- **Tier**: T2 | **Priority**: 2 | **Claim**: C1, C2, C5

---

## Family N — TNBC hot/cold lead clinical case (H59–H62)

Lead clinical application. Four-class decomposition using the platform's composite scores, tested in TNBC as the decisive subtype, with external-cohort and wet-lab tie-ins (see `07_wet_lab.md`).

### H59 ✱ — Four-class hot/cold decomposition
- **Statement**: TNBC (and broader BRCA) segregates into four mechanism-attributed classes — **hot-competent**, **hot-discordant**, **cold-primed-not-triggered**, **cold-intrinsically-invisible** — defined by `{IFN-γ, CYT, APC, JSI, COHER, CRML, COMB_HLA}`. Classes differ in ICB response, in tumor-intrinsic reversibility, and in optimal therapeutic stratum.
- **Atlas**: RNA, RPPA, SNV, SV, CNV, Methyl, Thorsson, Clinical; external ICB cohort for validation.
- **Features**: All six L2 composite scores; class assignment via threshold or clustering.
- **Tier**: T1 (class assignment) + T2 (external ICB validation) | **Priority**: 1 | **Lead-set**: yes | **Claim**: C1, C3, C5

### H60 — "Hot-discordant" class is enriched for `PRD + MIR + 3'UTR-SNV` signature
- **Statement**: The hot-discordant class (high infiltrate, low presentation competence) is mechanistically driven by proteostatic loss: high `PRD_g` + high `MIR_g` + enriched 3'UTR variants on APM genes.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C3

### H61 — "Cold-primed" class is STING-responsive, demethylation-non-responsive
- **Statement**: Cold-primed tumors (`JSI` high, `COHER` high, `CRML` low) retain inducibility and respond to cGAS–STING agonists; cold-intrinsic tumors do not. Wet-lab W4 tests this directly.
- **Tier**: T2 (computational) + wet-lab W4 | **Priority**: 2 | **Claim**: C1, C3

### H62 — TNBC subtype (Lehmann or Burstein) × four-class distribution is non-uniform
- **Statement**: Lehmann / Burstein TNBC sub-subtypes map non-trivially onto the four-class decomposition; specific sub-subtypes concentrate in specific mechanism classes (e.g., LAR enriched for cold-intrinsic; BLIA enriched for hot-competent).
- **Features**: Lehmann/Burstein TNBC sub-subtype + class assignment.
- **Tier**: T1 | **Priority**: 2 | **Claim**: C1

---

## Cross-family summary table

| Claim | Hypotheses supporting |
|---|---|
| **C1** — gene- and subtype-specific mechanism | H1, H3, H7, H8, H14, H19–H22, H25, H26, H29, H31, H32, H35–H37, H38–H44, H46, H47–H50, H54, H56–H58, H59, H61, H62, H63–H65 |
| **C2** — element-resolved causal inference | H1, H2, H9, H16, H22, H23, H34, H45, H46, H51, H55, H56, H63, H65 |
| **C3** — proteostatic visibility | H4, H5, H6, H10, H11, H13, H15, H21, H26, H27, H30, H31, H32, H36, H59, H60, H61 |
| **C4** — platform generalizes | H17, H37, H51–H53 |
| **C5** — intrinsic visibility is prognostic | H6, H14, H16, H18, H59, H64, H65 |
