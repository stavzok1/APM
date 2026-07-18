# The miRNA progression landscape — consolidated characterization report

> **Goal:** the readable, standalone synthesis of the GTEx-healthy → NAT → tumour miRNA trajectory biology —
> the 17-thread descriptive map produced by `analyses/progression/landscape_characterization.py`.
> **What belongs here:** the consolidated NARRATIVE + headline numbers of the progression *characterization*
> (structure · mechanism · players · validation). NOT the atomic record (that is the sync-partner) and NOT the
> null-tested realization findings (those are MH-158/160/162/163 in the registry).
> **Update trigger:** the characterization module gains/changes a thread, or a thread's verdict is corrected.
> **Sync-partner:** `BIOLOGICAL_SYNTHESIS.md §12` (the atomic biology home — edit both together).

> ⚠ **READ FIRST — this is a DESCRIPTION, not a null-tested claim set.** It reads structure off annotations with
> known caveats: (1) `shift_class` per-edge FDR labels rest on the MH-123/124 null measured **3–4× too narrow**
> (`LANDSCAPE_REPORT.md`); (2) gene-role annotations are incomplete (many targets `unknown`). Read the *patterns*
> (cluster co-movement, convergence, the mechanism splits), not per-edge significance. The **null-TESTED** progression
> results live in the registry: **MH-158** (paired realization), **MH-160** (class→realization, repression-liveness),
> **MH-162/163** (within-sample + program-specificity negatives).
> ⚠ **Cross-platform (GTEx TPM → TCGA RPM):** trust the RANK axis (percentile deltas dHT/grank) + Shapley over the
> QN magnitude bridge (`log2fc`); QN is correct+settled but a softer assumption. The **paired NAT→tumour** dose is
> TCGA same-platform ⇒ not a QN concern, and the cohort (~1000) view agrees with it at ρ=0.94.

---

## Headline

The malignant miRNA trajectory in TCGA-BRCA is **cluster-coordinated, convergent, two-step (field-dominated),
regulatory-rewiring, and mostly buffered** — and its acquisition mechanism is **heterogeneous** (promoter
hypomethylation + copy-number amplification each explain only a minority; the rest is transcriptional). It is
textbook-concordant (miR-200↔ZEB1, miR-17~92, 14q32/DLK1–DIO3, miR-21→PTEN, C19MC).

---

## 1. STRUCTURE — how the trajectory is organized

- **Cluster-coordinated (co-transcription drives it).** Genomic clusters move as *units* (dose-concordance 0.9–1.0):
  oncomiR polycistrons **ACQUIRE** — miR-17~92 (MIR17HG), its paralogues miR-106a~363 / miR-106b~25, and
  miR-183/96/182; tumour-suppressor clusters **LOSE** — the **14q32/DLK1–DIO3 imprinted locus** (MEG8 concordance
  1.00, 12 arms; MEG9 21 arms), miR-30, let-7/miR-100/125.
- **Convergent (cancer genes as multi-family sinks).** **PTEN is the supersink — 28 acquired repressive edges from
  20 distinct seed families**; then CCND1 (17/14), CDKN1A (15/10), ZEB1 (dose 2.10), EZH2, STAT3, BCL2, FOXO1,
  TGFBR2, PDCD4, RB1, FBXW7.
- **Two INDEPENDENT waves (field-dominated).** Decomposing the healthy→tumour rank move: the **field step**
  (healthy→NAT, 150 arms) is larger than **and uncorrelated with** (ρ(dHN,dNT)=0.003) the **malignant step**
  (NAT→tumour, 32 arms). The DIRECT GTEx→tumour difference tracks the field step (ρ=0.84) far more than the
  malignant step (ρ=0.45) ⇒ *what's different in cancer vs truly-healthy is mostly established pre-malignantly.*
- **Regulatory rewiring.** 333 genes switch their **dominant regulator** healthy→tumour, typically an oncomiR taking
  the seat (PTEN miR-148a→miR-21, AGO2 miR-99a→miR-375, HIF1A miR-199a→miR-30a, MYB miR-155→miR-16).
- **Mostly buffered.** Functional 2×2: 455 genes acquire dose but DON'T realize vs 284 DRIVERS (ESR1/SMAD3/BCL2).
  Even PTEN — the 28-edge supersink — is buffered at the gene-aggregate level (dose ≫ effective repression).

## 2. MECHANISM — what drives acquisition and loss

- **Copy number: a real SECONDARY contributor** (⚠ corrected from an initial coarse-metric null). Amplified-locus
  arms (CN>2.5) acquire more than diploid (+0.127 vs −0.044, MWU p=0.035); key oncomiRs **miR-96/182/21 are
  CN-amplified**. But only ~11% of acquirers are CN-driven ⇒ not dominant genome-wide.
- **AGO/RISC is NOT the realization bottleneck** — edge realization is flat across `ago_loading` (ρ=0.08).
- **Epigenetics: two PARTIAL routes + a refutation** (methylation Δβ gate, imprinting-aware; positive controls
  miR-124/9/129/137 pass at Δβ +0.18). Acquisition splits into **promoter-hypomethylation** (miR-21/141/200c/17,
  Δβ −0.14…−0.18) and **CN-amplification** (miR-96/182) — but together only ~35% of top acquirers; the majority are
  TF/enhancer-driven. **⛔ The 14q32 coordinated LOSS is NOT hyper-methylation silencing** — the imprinted locus
  (baseline β 0.87) *hypo*-methylates in tumour (73% of arms) ⇒ hypothesis refuted, loss mechanism still open.
- **Program effect needs ARCHITECTURE, not a flat mean** (`geneset_architecture`, signed + master-regulator +
  gene-role). Acquired miRNA pressure is **pro-tumour by repressing P53/apoptosis** (+18/+3) and **anti-tumour by
  damaging G2M/glycolysis** (−11) — a sign a bag-of-genes average completely misses.

## 3. THE PLAYERS

- **Acquired oncomiRs → tumour suppressors:** miR-200 family (patient-specific, own-frac 0.71–0.75) + miR-96/183 +
  miR-21 → PTEN, FOXO1 (miR-96→FOXO1 ρ −0.32), CDKN1A, TXNIP, DLC1, PTPN14.
- **Silenced TSG-miRs → oncogenes:** miR-145/486/451/335/15-16-195-497 → CDK4, AKT1, BIRC5, IGF1R, E2F3, HER2.
- **Trajectory archetypes:** most change is field-effect (91 arms) ≫ tumour-acquired (12); acquisition is magnitude
  more than rank (⚠ the magnitude-only set rests on the softer QN bridge — the rank-supported set of 101 is the
  confident one).
- **Subtype:** subtype-specific acquisition concentrates in Basal/TNBC; distinct clusters per lineage — LumA→C19MC
  (chr19 miR-524/525/526), Her2→miR-371~373.
- **Genomic context predicts direction:** intergenic arms ACQUIRE (dose +0.47), the imprinted 14q32 locus LOSES
  hardest (68% of its arms losing).

## 4. VALIDATION & DISCIPLINE

- **Paired (103) ≈ cohort (~1000).** The 103-pair `mean_own_shift` and the full-cohort `arm_lfc_NAT_TUM` agree at
  Spearman ρ=0.94 (zero big sign-divergences) ⇒ the paired-based map generalises; it is not a matched-subset
  artifact. `coupling_tum` (~1000-tumour) is the well-powered cohort coupling.
- **Trust rank + Shapley over QN** for anything crossing GTEx→TCGA (above).
- **What is DESCRIPTIVE vs TESTED.** Everything in this report is descriptive. The null-tested spine is the registry:
  MH-158 (paired realization is real, weak, set-level), MH-160 (the cross-state taxonomy predicts within-patient
  realization by repression-liveness, A), MH-162/163 (within-sample biomarker + program-specificity are negatives).

## 5. OPEN QUESTIONS

- **The 14q32/DLK1–DIO3 coordinated-loss mechanism.** Three drivers are now RULED OUT — not CN (§2), not AGO (§2),
  not promoter hyper-methylation (§2, imprinting-aware). Remaining candidates: the **imprinting control region
  (IG-DMR / MEG3-DMR)** or a transcriptional/`MEG3`-lncRNA effect.
- **What drives the ~60% of acquirers with no CN/hypomethylation route** (TF/super-enhancer hypotheses).
- **Cross-cohort replication** — everything is single-cohort TCGA; the cap on strength (A not S) is external
  replication, which is data-gated.

---

## Provenance

- **Module:** `mirna_hallmark/analyses/progression/landscape_characterization.py` (17 threads; `characterize()`).
- **Outputs:** `output/learned/realization/landscape_*.tsv` (convergence_hubs, cluster_units, driver/loser_roster,
  subtype_*, trajectory_archetypes, regulatory_handoffs, dose_realization_quadrants, genomic_context_census,
  program_acquired_dose, trajectory_two_step, healthy_anchored_drivers, cohort_view, direct_gtex_tumor, cn_layer,
  ago_layer, epigenetic_mechanism/routes).
- **Reused (not reinvented):** the integrated `progression_edge_card`/`progression_gene_card` (MH-158/160),
  `mirna_state_class` (rank trajectory + PAM50), `genomic_context`, `mirna_locus_cnv`, `learned.methylation`
  (imprinting-aware Δβ gate), `geneset_architecture` (signed architecture-weighted program pressure).
- **Atomic home:** `BIOLOGICAL_SYNTHESIS.md §12`. **Tested findings:** registry MH-158/160/162/163.
