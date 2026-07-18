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

The malignant miRNA trajectory in TCGA-BRCA is **cluster-coordinated (incl. real across-loci trans-coordination),
convergent (PTEN a genuine outlier), and two-step (field-dominated)** — with **acquisition ≫ realization** (a
continuum, not a bimodal split). Its acquisition mechanism is **heterogeneous** (promoter hypomethylation +
copy-number amplification each explain only a minority; the rest transcriptional). *(Rigor-audited: the
"regulatory-handoff" and "buffered/driver 2×2" claims were downgraded/retracted — see §1.)* Textbook-concordant
(miR-200↔ZEB1, miR-17~92, 14q32/DLK1–DIO3, miR-21→PTEN, C19MC).

---

## 1. STRUCTURE — how the trajectory is organized

*(Rigor-audited 2026-07-18 — `landscape_hardening_audit.tsv`. ✅ = survives a proper test; ⚠ = soft/downgraded.)*

- ✅ **Cluster-coordinated — and REAL beyond co-transcription.** Same-host (co-transcribed) arms co-vary
  trivially (paired Δx-corr +0.57, shared primary transcript). The non-trivial finding: **across-loci, same-family
  arms co-vary at +0.29 vs random +0.12** (excess +0.17) — genuine *trans*-coordination of paralogous loci
  (miR-200 chr1↔chr12, the miR-17~92 super-family). OncomiR polycistrons ACQUIRE; the 14q32/DLK1–DIO3 imprinted
  locus, miR-30, and let-7/100/125 LOSE.
- ✅ **Convergent — PTEN is a genuine outlier** (not just "many regulators"). Acquired-edge count partly tracks a
  gene's total regulator count (ρ=0.61), but **PTEN has +11.8 acquired edges beyond its regulator-count
  expectation — the top residual** (then WEE1, CDKN1A, TGFBR2, XIAP, TGFB1). *(Report the residual, not the raw
  count of 28 — PTEN has 90 total regulators.)*
- ✅ **Two INDEPENDENT waves (field-dominated).** field step (healthy→NAT, 150 arms) larger than **and uncorrelated
  with** (ρ(dHN,dNT)=0.003) the malignant step (NAT→tumour, 32 arms); direct GTEx→tumour tracks the field (ρ=0.84)
  » malignant (ρ=0.45). Rank-based, trusted.
- ⚠ **Regulatory rewiring — DOWNGRADED (mostly argmax-noise).** Of 292 genes whose *dominant* regulator (by
  abundance rank) differs healthy→tumour, only **~64 (22%) are DECISIVE**; the rest are close-call flips (median
  tumour margin **2 percentile points** — a coin-flip). ⚠ Even the flagship **PTEN miR-148a→miR-21 is a TIE in
  tumour** (margin 0.0). The *dose* rise of miR-21 is real (§3); the "dominant-regulator handoff" framing is weak.
- ⛔ **Functional 2×2 — RETRACTED (axiom-5 threshold artifact).** The driver/buffered counts are entirely
  threshold-dependent (driver fraction 0.54→0.38→0.25→0.14→0.06 as the cut moves −0.05→−0.25; 28% of genes sit
  within ±0.05 of the −0.1 boundary; adjacent-cut overlap 0.64). There is **no stable bimodal split** — realization
  is a *continuum* of acquired dose. The "455 buffered / 284 driver" headline is void. *(That acquired dose ≫
  realized dose overall is the real, tested statement — see MH-160/163, not this 2×2.)*

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
