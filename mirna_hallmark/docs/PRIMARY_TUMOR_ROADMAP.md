# Primary-tumor leverage roadmap — `mirna_hallmark`

**Goal.** Take the rigor/compartment lens proven on the DCIS/EV arc (MH-48..56; `apm-rigor-protocol`)
to **primary breast tumours** and the project's **core framework** (MH-1..47). The core is built on
**bulk TCGA** — so it carries the same composition confound that flipped miR-29c (MH-54) and that
decoupled the protein programs from miRNA pressure (MH-56d). The roadmap tests, gradually, **whether
the miRNA-pressure framework holds up cell-resolved in primary tumours.** Execute phase-by-phase;
each phase reuses existing machinery and ends with a registry MH-id.

## The binding constraint (read first)
**No measured single-cell / spatial miRNA exists for breast TUMOURS — primary or in-situ** (the field
is scRNA = mRNA). BUT measured single-cell miRNA **does exist for CELL LINES** (the method papers:
GSE114071 K562 miRNA+mRNA co-seq; MCF7 RNAome / method-comparison) — so the **non-circular per-cell
miRNA↔target test is possible in vitro** (Phase 4c). For tumours, every per-cell-miRNA route is
**(i) deconvolution** of bulk miRNA via mRNA cell-fractions, **(ii) inference** (SiCmiR/STmiR —
circular for miRNA→target, hypothesis-only), or **(iii) sorted/LCM bulk** small-RNA. The
mRNA/target/compartment side is fully cell-resolved; the tumour miRNA side is not. Plan accordingly.

## Reference data to stage
- **Primary scRNA:** GSE176078 (Wu/Swarbrick, 26 primary + spatial), GSE161529 (68 specimens), Pal/Visvader.
- **Primary spatial:** Visium/Xenium/MERFISH breast; MIBI/CODEX/IMC primary (Risom is DCIS).
- **Core/bulk (have):** TCGA-BRCA miRNA/mRNA/CNV/methylation/RPPA; CPTAC primary proteome+miRNA.
- **METABRIC (have, UNDER-USED):** `data/METABRIC/` (cBioPortal: mRNA microarray + CNA + clinical
  **with long-term outcomes**, ~2,000 tumours; currently only `edge_prior_refinement` touches it) +
  GSE22216 METABRIC-miRNA (used by `buffa_validation`). A large outcome-bearing validation/retest set.
- **Controlled DCIS outcome:** HTAN Breast PreCancer Atlas (Strand et al., Cancer Cell 2022,
  TBCRC038+RAHBT; 677 DCIS, recurrence outcomes, 812-gene classifier) — Level-1/2 via **dbGaP
  phs002371**; derived (Strand bulk/spatial + Risom MIBI, already used) largely open on the HTAN portal.
- **Progression (have + scarce):** GSE59247/162670/93740/297448/METABRIC; new EV candidate GSE293809
  (nipple-discharge/fluid/tumour EV miRNA); controlled HTAN DCIS/TBCRC038 for outcome power.
- **Inferred sc-miRNA:** SiCmiR (github Cristinex/SiCmiR), STmiR — circular caveat.

---

## Phase 0 — Foundation (deconvolution reference) — ✓ COMPLETE → MH-70 (BayesPrism); MH-58 (NNLS partial)
Built cell-type signature from GSE176078 (streamed 2.4 GB mtx → 9-type pseudobulk) + NNLS-deconvolved
1,095 TCGA samples (`brca_deconvolution.py`). **Validated (vs metagenes): immune ρ=0.85, CAF ρ=0.57 —
USABLE (the stromal/immune confounders); epithelial ρ=0.39 + endothelial over-fit — UNRELIABLE.**
The quick NNLS suffices for CAF/immune; epithelial needs a **proper tool (CIBERSORTx/BayesPrism,
gene-length-matched)** — that's the remaining Phase-0 upgrade. Until then the marker-metagene proxy is
primary for epithelial (so MH-57a's 93% stands). Output: `tcga_celltype_fractions.tsv`.

## Phase 1 — Does the core framework survive composition? **(the headline test) — ✓ DONE (full scale) → MH-57**
**FULL-SCALE result (`core_coupling_composition_retest`, ALL 1,013 edges, CPE+HRD+batch baseline):** **93% of headline BY-neg couplings stay negative** after epi/immune/stroma+proliferation adjustment (median ρ −0.148→−0.124;
774/1013 still FDR-neg) → the framework is **largely composition-robust**; the attenuated ~7% is precisely the **miR-29→collagen/ECM** class (miR-29b-3p→COL5A1/2/3, miR-29c-3p→PDGFRB). **Refined version DONE (MH-71): with BayesPrism deconvolution fractions, 96% stay negative (>93% metagene) — framework MORE robust under gold-standard composition adjustment; miR-29→ECM the known exception.** Companion: `program_pressure_by_role` (MH-57b, gene-role pressure).
Re-run the framework's headline miRNA-Hallmark **couplings (MH-1..47) composition-adjusted** with the
Phase-0 fractions (partial Spearman, adjusted-primary; `apm-rigor-protocol` B). 
*Establishes:* which couplings are tumour-cell-intrinsic vs stromal/immune-composition — the
primary-tumour version of the miR-29c lesson. *Caveat:* deconvolution noise; only the target/mRNA
side is cell-resolved (miRNA bulk). *Deliverable:* MH-id "core couplings, composition-retested".

## Phase 2 — Compartment-annotate the network + definitive orphan retrofit — ✓ 2a DONE → MH-59
**2a DONE (`compartment_annotate_network.py`):** all 1,013 headline edges annotated by target-gene compartment (Wu scRNA, peak-celltype). miR-29→ECM targets CAF-enriched (38% vs 19%) = gene-resolved confirmation of MH-57; only ~19% of targets epithelial-specific. **2b (pending):** definitive orphan retrofit on the TCGA spine.
(a) For every framework edge, label **which compartment expresses the target** (scRNA) → tumour /
CAF / immune / endothelial / cross-compartment(paracrine). (b) Run the **definitive orphan retrofit**
(MH-55 was n=40 DCIS screening) on the TCGA spine + composition metagenes. *Establishes:* a
compartment-annotated regulatory network; which orphans survive composition at primary scale.

## Phase 3 — Spatial localization + outcome — ✓ OUTCOME (MH-60/61/62/63) + ✓ SPATIAL (MH-64) DONE
**Outcome DONE (`outcome_survival.py`):** TCGA OS + GSE22216 DRFS, composition-adjusted Cox on 145 headline arms. No single arm FDR-sig for OS; 3 replicate cross-cohort (protective: miR-342-3p/miR-30a-5p/miR-181c-5p). Follow-up: aggregate/pressure-signature outcome; METABRIC mRNA outcome when EGA miRNA lands. **Pressure×resolution + per-subtype DONE (MH-61):** pressure > expression for outcome (AGO-capacity→OS LumA; gene-pressure MKNK2→Basal; METABRIC miR-210-3p hypoxia RISK pan-ER). **Independent Buffa/METABRIC retest (MH-62/63):** TCGA pressure-outcome does NOT cross-cohort replicate; durable anchor = miR-210/hypoxia *expression* (known).
**Spatial component DONE (MH-64) — `spatial_*` suite, 4 MEASURED resolution layers:** single-cell PROTEIN (Risom DCIS MIBI), Visium SPOT primary (Wu/Zenodo, 6 primary tumours), Visium SPOT public (10x block, replication), single-cell RNA (Xenium, 168k cells). **Localization unanimous across all 4** (proliferation/glycolysis→epithelial, hypoxia/EMT→CAF, immune→immune; programs spatially structured into niches, median Moran's I 0.45–0.60). **MH-56(d) decoupling reproduces in situ** on the MIBI normal→DCIS→IBC axis (SLC2A1/FAP brake-release, ERBB2/ESR1/CDH1 concordant-repression, HIF1A/VIM discordant-rise), and the **brake-release genes spatially track their programs** in primary Visium (SLC2A1→glycolysis, FAP→EMT, 7 sections). Binding constraint held: no measured spatial miRNA → miRNA side is route-(i) composition projection only; inferred sc-miRNA (route ii / SiCmiR) deferred to Phase 4a.
Place the miRNA-Hallmark **programs** (proliferation/hypoxia/glycolysis/EMT/immune) in primary-tumour
compartments and niches (Visium/Xenium/MIBI), with **survival/metastasis outcome**. Re-test the
**MH-56(d) program↔pressure DECOUPLING** spatially (programs vs the deconvolved miRNA pressure):
are the programs driven independently of miRNA pressure in primary tumours too? *Establishes:* where
the programs live + whether pressure is a brake or a driver, with outcomes.

## Phase 4 — The miRNA side, bounded (incl. the MEASURED in-vitro route) — ✓ 4c DONE → MH-65
**4c DONE (`sc_mirna_target_k562.py`):** acquired **GSE114071** (Wang 2019 half-cell co-seq) — **19 K562 cells**,
each physically split into a small-RNA half (miRNA) and an RNA-seq half (mRNA): the *only* data where a
miRNA and its target mRNA are MEASURED in the same single cell. Pairing validated by the shared physical
**HalfCell** id (intersection = exactly 19, matches the paper). Tested the 1,013 headline BY-neg edges as
per-cell Spearman. **Result = honest INCONCLUSIVE:** headline edges null (229 testable, 51% neg, median
sc_ρ −0.009, no decoy specificity), BUT the **max-power positive control** (abundant miRNA → broad curated
targets vs decoy) is only weakly directional (cognate −0.016 vs decoy −0.004, MWU p=0.17) → the binding
limit is the **data** (19 cells × single-half-cell dropout × log2-floored miRNA), not the framework. The one
arm with real dynamic range (**miR-92a-3p**) carries the expected sign (median target ρ=−0.040, n=377).
Net: executing the only measured per-cell route returns neither confirmation nor refutation at this
resolution, and **empirically motivates the 4b sorted/LCM ask**.
**4c breast-retarget DONE (`ccle_breast_target_anticorr.py` → MH-67):** to answer "why leukemia?", retested the
edges **across 50 CCLE breast lines** (subtype-resolved Basal/HER2/Luminal, n=50 ≫ K562's 19, matched
NanoString miRNA + DepMap RNA-seq). Edges are **positively** correlated (median ρ +0.046, subtype-robust,
0 FDR-neg), and it is **DOUBLY inconclusive**: (i) the **positive control fails** (gold-standard
miR-200→ZEB1/let-7→HMGA2 not recovered, median ρ +0.07) ⇒ cross-line axis insensitive to repression
(co-regulation dominates); (ii) the **miRNA data-validity gate FAILS** — the 2018 NanoString miR-200c/141
don't even mark known epithelial identity (epi 6.16 ≤ mes 7.38) while the DepMap mRNA EMT markers are
textbook ⇒ the NanoString miRNA is the broken link. A cleaner miRNA-**seq** swap would help but **no large
breast-line small-RNA-seq matched to mRNA exists** (CCLE miRNA = the NanoString panel; GEO panels ~9 lines).
Same honest verdict as K562: **INCONCLUSIVE**, the measured route can't deliver the repression test. **4a/4b still open.**
(a) **Inferred sc-miRNA** (SiCmiR) on primary scRNA → per-cell-type miRNA *activity*, hypothesis-only
(never causal for target edges; `apm-rigor-protocol` D). (b) **Scan/acquire sorted-or-LCM compartment
small-RNA** (non-circular per-compartment miRNA) for tumours; else the wet-lab ask (FACS
tumour/CAF/immune → miRNA-seq). **(c) MEASURED matched sc-miRNA+mRNA in CELL LINES — the non-circular
per-cell test, in vitro:** GSE114071 (K562 co-seq), MCF7 RNAome / method-comparison (Nat Commun 2021)
→ test the framework's edges directly as **per-cell miRNA↔target anti-correlation** (measured, not
inferred). Cell-line not tumour + small-n, but the **only measured per-cell miRNA↔target data that
exists** — the cleanest possible causation check for the pressure axis. *Establishes:* whether the
framework's miRNA→target edges hold at true single-cell resolution where both are measured.

## Phase 5 — EV / liquid-biopsy + generalization
(a) Export/biomarker signature in primary-patient plasma-EV cohorts with staging/outcome (+ GSE293809).
(b) Pan-cancer: the rigor protocol + composition lens are tumour-agnostic; port to other types.

---

## Parallel track — progression power (small, ongoing)
Public DCIS→invasive miRNA cohorts are **scarce** (we hold most). Add **GSE293809** (EV) for the EV
angle; pursue **controlled HTAN DCIS/TBCRC038** (dbGaP) for outcome-powered progression. Don't expect
many more public miRNA progression cohorts — the scarcity is why MH-48/50 are set-level/underpowered.

## Sequencing recommendation
Do **Phase 0 → Phase 1 first** (the headline test, reuses everything). Phases 2–5 are independent and
can follow in any order as appetite allows. Each phase is one analysis module + an MH-id; nothing
requires controlled access except the optional HTAN progression-power track.
