# What's next — `mirna_hallmark`

Forward-looking extensions, the discoveries each could yield, and the open
questions / artifacts they would resolve. Ordered roughly by expected
payoff-to-effort. Anchors: findings in [`BIOLOGICAL_SYNTHESIS.md`](BIOLOGICAL_SYNTHESIS.md)
(MH-/D- ids), methods in [`METHODS.md`](METHODS.md).

Legend: **[discovery]** = new biology it could surface; **[resolution]** = a
current caveat/artifact it would settle; effort ≈ S/M/L.

---

## Tier 1 — highest payoff, builds directly on current results

### 1.-1 Aggregate "raw force" vs abundance-sum baseline, tumor↔GTEx Δρ (M) — **DESIGN**
Full plan in [`method_dev/aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md`](../method_dev/aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md); results in `method_dev/aggregate_pressure/AGGREGATE_PRESSURE_FINDINGS.md`.
Gene-rung test (aggregation lemma regime): does promiscuity-attenuated force
`Σ a(m,s)·w_eff/D(m)` detect net repression better than `Σ a(m,s)`, and is the
tumor→healthy Δρ stronger in the force form? Healthy = **GTEx breast (~300 paired,
already ingested)**, not NAT. Key new artifact = per-gene **contribution effective-N**
(`1/Σp_m²`) as moderator + falsification, wired into the existing `:473` crowding
stratification in `coupling_predictor_comparison.py`. Guardrails: partial-ρ on
CPE+HRD+prolif+target-CN/methyl, no circular promiscuity tuning, functional-abundance
floor, independent ±validation gene set. See skill `apm-gene-question`.

### 1.0 Proliferation + curation-bias robustness (M) — **DONE (`robustness_checks.py`)**
Pre-registered reviewer-objection controls (grant Aims 1–2), Basal n=197.
- **[resolution delivered]** The basal pressure↔repression coupling **survives and
  strengthens** under a proliferation covariate added to CPE+HRD — proliferation is a
  *suppressor* confound (MH-17): within Basal **8/8** key programs stay negative & p<0.05
  under E2F+G2M; **6/8** under MKI67 alone.
  **Over-adjustment control done:** because the E2F+G2M metagene shares the hub's
  E2F/MYC driver, we repeated under **MKI67 alone (6/8)** and an **E2F/MYC-stripped
  G2M+mitotic metagene (**8/8**, 271 genes)** — strengthening reproduces, so it is not a
  shared-E2F artifact. `prolif_sign_structure.tsv` shows the suppressor sign
  structure (prolif positively associated with both pressure and target expression).
  PTEN/CDKN1A/TGFBR2/VIM routes survive (MH-19).
- **[resolution delivered]** The hub is **not** a study-count artifact (MH-18):
  miR-17-5p ranks #1 under evidence, binary/degree, **and** TargetScan weighted-
  context++; family tops all three (rank-corr ρ=0.65 evidence-vs-TargetScan).
- **[caveat surfaced]** The **IRF1 route fails** proliferation adjustment within
  Basal and is reassigned to miR-130b/301 by sequence prediction → demoted to
  exploratory; this is what §1.2 must now resolve before any causal claim.

### 1.1 Pin the cold-Basal axis to the *archetype*, not bulk PAM50 (M) — **DONE (Hallmark + gene panel)**
Implemented in `visibility_archetype_contrasts.py`: all 50 Hallmarks + 10-gene
route panel (IRF1, PTEN, FOXO1, CDKN1A, BCL2L11, TGFBR2, TP53, STAT1, NLRC5,
CIITA) for cold_Basal / hot_Basal / other_Basal / hot_Luminal / other_Luminal.
- **[discovery delivered]** Apoptosis, p53, PI3K-AKT-mTOR, Notch join immune
  programs in cold_Basal > hot_Luminal pressure gap; per-gene PTEN/BIM/FOXO1
  resolution like IRF1.
- **[still open]** Within-basal cold vs other_Basal Hallmark-pressure deltas are
  ~0 — refine with per-sample attribution (§1.2) or module-level APM scores.

### 1.2 Per-sample hub-route OLS explainability (M) — **in progress (`gene_expression_explainability.py`)**
**Context:** Partial ρ / antagonism % are cohort association summaries. The next
step is **within-sample × gene variance decomposition**: nested OLS on
`log2(TPM+1) ~ CN + CPE + PAM50 + promoter meth [+ route miRNA block]`, mirroring
APM `cn_adjusted_rna_states` and `per_gene_explainability`.

**Module:** `.venv/bin/python3 -m mirna_hallmark.analyses.misc.gene_expression_explainability`
(cohort + `--pam50 Basal`). miRNA blocks: route HE pressure, family means, per-arm
expression. Outputs: `output/gene_expression_explainability/{cohort,pam50_Basal}/`.

**Early read (Basal, `core_plus_route_arms` vs CN+CPE+prolif+meth+SV core):**
CDKN1A ΔR²≈0.16 (all) / 0.05 (rna_low); PTEN ≈0.11 / 0.05; IRF1 **curated route ≈0.006**
but **TargetScan alt (miR-130b/301) ≈0.086 (all) / 0.30 (rna_low)** — rescues IRF1 at
the OLS layer with the correct arms. Per-subtype table:
`incremental_r2_by_subtype_summary.tsv`.

- **[discovery]** Which hub genes have miRNA-explained variance beyond cis layers, **by subtype**?
- **[resolution]** IRF1 grant route should use **`curated_plus_targetscan`** (union of
  miR-23a/106b + miR-130b/301); combined Basal ΔR²≈0.10 (all) / 0.32 (rna_low), vs
  alt-only 0.086 / 0.30 and curated-only 0.006 / 0.17.

### 1.3 Stratum-resolved anti-correlation done right (S)
Current cohort negatives are partly Normal-like-driven (MH-8). Recompute the full
anti-correlation **within each subtype group with PAM50/purity already removed**,
and report only within-stratum effects.
- **[resolution]** Removes the Normal-like admixture confound from the headline
  repression claims; tells us which repression is luminal- vs basal-intrinsic.

---

## Tier 2 — sharpen the pressure model

### 2.1 Family-arm-aware pressure (S/M)
The miR-29a/b vs miR-29c opposite-direction split (§7) shows arm-level structure.
Add a per-seed-family pressure roll-up and flag intra-family discordant arms
cohort-wide.
- **[discovery]** Systematic catalog of "tug-of-war" families where arms oppose.
- **[resolution]** Avoids cancelling opposing arms inside a family-level score.

### 2.2 Degree-normalized differential pressure (S)
Raw differential favors hub genes; the z-version ties single-miRNA genes (§6).
Add a **per-edge** differential (Δ realized pressure per (miRNA,gene) by group)
and a permutation-null effect size.
- **[resolution]** A cross-gene-comparable differential-pressure ranking without
  either artifact.

### 2.3 Replace the heuristic AGO gate with a fitted interaction term (M)
Instead of a fixed logistic gate, fit `expression ~ pressure * AGO_capacity`
(+ confounders) per Hallmark and test the **interaction** coefficient.
- **[discovery]** Whether RISC capacity *actually* modifies the pressure→expression
  slope (would upgrade AGO from "motivated" to "demonstrated", MH-6/limitation).
- **[resolution]** Settles whether the bounded gate's "rescale-not-reorder"
  behavior is a modeling artifact or a real null.

---

## Tier 3 — new evidence layers (reuse parent APM data)

### 3.1 miRNA-locus CNV → target de-repression coupling (M)
Join miRNA-universe locus CNV (`mirna_locus_cnv`) to target-gene expression: do
**lost** high-evidence miRNA loci coincide with **up-regulated** Hallmark targets?
- **[discovery]** Loss-of-repression events (structural, causal-leaning) per
  Hallmark — the natural follow-up to MH-2/MH-4 (registry hypothesis MH-H2).

### 3.2 Bring in the parent miRTarBase *isoform/3'UTR* and SV layers (L)
Cross the EMT/NF-κB pressured genes with APM SV/3'UTR-disruption to find tumors
where the miRNA binding site itself is altered.
- **[discovery]** Sample-level escape from miRNA repression (binding-site SVs).

### 3.3 Protein-level check via CPTAC (M) — **DONE (`cptac_validation.py`, MH-33 / MH-34)**
Tested whether **protein** (CPTAC) tracks predicted miRNA pressure, on two cohorts:
same-patient **TCGA-105** iTRAQ (MH-33) and **independent CPTAC-2 prospective** breast
(MH-34, self-contained pressure from CPTAC-2's own GDC miRNA-seq). Result: the spine
**validates at protein in the independent cohort** (53 genes FDR, ZEB1/EMT led) but not
in the noisier same-patient iTRAQ; a slope-free **`protein_resid`** layer shows the
repression is **mostly mRNA-mediated** (ZEB1 −0.55→−0.24 NS) — translational component small.
- **Next refinement (open): single-(miRNA→gene)-edge pressure vs protein.** MH-33/34
  used gene-level *aggregate* pressure (all regulators summed), which can dilute a strong
  arm sitting among weak co-regulators. Resolve per-edge pressure via
  `compute_gene_pressure_contributions` and test the specific arm→target protein coupling
  (e.g. miR-200→ZEB1) — sharper than the aggregate, and the natural test of the hub arms.
- **Optional:** extend `protein_resid` + single-edge to the phospho/acetylome layers (deferred per scope).

### 3.4 External replication in METABRIC (M)
Re-run the targeting map + basal immune-pressure axis on METABRIC (Illumina
expression; no miRNA-seq, so use target-expression + edge structure only).
- **[resolution]** Cohort-independent confirmation of the subtype contrasts.

---

## Tier 4 — scope / infrastructure

### 4.0 Paralog MIMAT CN aggregation (S) — **DONE (`dosage_landscape_cnv.assemble_mirna_cnv_long_table`)**
MIMAT arm CN is aggregated from hairpin loci (not a single host-gene overlap row).
Methods: `median` | `mean` | `max` | **`expression_weighted`** (cohort-median
locus×MIMAT reference RPM — build with `build_mirna_locus_mimat_expression.py`;
same weights for all samples, not per-sample expression).
- **[resolution]** Fixes triplicated arm rows and expression↔CN mismatch for paralog arms.

### 4.1a Extend MirGeneDB coords to high-abundance expr-only MIMATs (M) — **DEFERRED**
~1,418 GDC MIMATs lack MirGeneDB loci; ~24 have mean log2RPM>2 (lead gap:
**miR-625-3p**, rank ~94). Tier-B coordinate extension from miRBase/miRGenDB for
that tail only — not full 2k expansion.
- **[discovery]** CN↔expression for dosage-visible arms currently absent from CNV tables.

### 4.1b Coordinate overlap audit + bridge layers (M)

Three independent overlap joins still have gaps:

| Join | Mapped | Gap | Next step |
|------|--------|-----|-----------|
| MI* → GDC hairpin (expr weights) | 348/506 | 158 zero-bp (MirGeneDB v3 vs miRBase v21 spans) | **A.** Curated `MI*↔hsa-mir-*` table from miRBase registry; **B.** Use GDC hairpin coords as fallback locus definition for unmapped high-RPM arms |
| MI* → ASCAT3 segment (CN) | ~506 loci queried | loci with no overlapping segment → `cn_state=NA` | Report per-locus `pct_mapped` in `entity_segment_reference_map.tsv`; optional ±N bp flank (`mirna_cnv_segment_flank_bp`) |
| MIMAT expr-only (no MI*) | ~810 in coord catalog | ~1,418 expr-only | §4.1a — extend coords for high-abundance tail only |

Deliverables when prioritized:
1. ~~`scripts/mirna/audit_mirna_locus_maps.py`~~ — **DONE** (`data/miRNA/audit/`)
2. `data/miRNA/mirbase_hairpin_hg38.bed` + name bridge (optional canonical layer between MirGeneDB and GDC)
3. Sensitivity run: segment vs gene CN concordance on paralog arms

### 4.1 Vectorize per-sample segment×locus overlap (S)
Segment-mode `mirna_locus_cnv` loops samples × loci. An interval-join (pyranges /
sorted merge) per chromosome would cut wall-clock substantially.
- **[resolution]** Removes the main slow step from `run_all`.

### 4.1c Legacy entity→ENSG gene path (S) — **optional / deprecated default**
`mirna_cnv_source=ascat3_gene` retained for comparison. ENSG vectorization (old §4.1) only applies to that path.

### 4.2 Multi-set governance (S)
Genes in many sets are counted per context (dedup on (miRNA,gene) for totals).
Add an option for set-size-weighted / Jaccard-deduplicated Hallmark scoring.
- **[resolution]** Quantifies how much of a Hallmark's signal is shared vs unique.

### 4.3 Figures (S)
Heatmaps for the immune-pressure axis (subtype × immune-Hallmark) and the IRF1
route; volcano for differential gene pressure. Currently tables only.

---

## 5. DCIS / pre-malignant / EV directions (after MH-48..51)

These follow from the in-situ→invasive + EV work. Implemented so far: timing (MH-48),
within-patient mRNA coupling (MH-49), EV export (MH-50), confounder-hardening +
cellular/LCM triangulation + direct export (MH-51). Forward, **not yet built**:

### 5.1 EV recipient-cell delivery — GSE297447 (P)
GSE297447 is the EV **mRNA** companion to GSE297448. Test whether the targets of the
*exported* miRNAs (miR-140-3p, miR-205) are de-repressed in EV-recipient cells / show up
as EV cargo — i.e. does export have a delivery consequence, not just shedding.
- **[mechanism]** Closes the export→delivery loop the MCF10A EVs only half-answer.

### 5.2 Rab27A / secretion-machinery axis (P)
The GSE297448 source paper attributes export to **Rab27A**. Test whether the
export-selectivity (cellular-down→EV-up) tracks Rab27A / secretory-pathway expression
across the progression (where mRNA is available), to separate "selective packaging" from
"bulk shedding".

### 5.3 Plasma-EV biomarker bridge (S)
Connect the tissue/cell-line export signature to the existing **plasma-EV cohorts**
(`ev_mirna_multi_cohort`): do the directly-exported arms (miR-140-3p, miR-205) appear
elevated in patient plasma EVs vs normal? Tissue export → circulating biomarker.

### 5.4 Duke / Strand DCIS atlas — RNA + outcome complement (verdict: YES, two uses)
The Duke Breast PreCancer Atlas (Strand et al.; TBCRC 038 + RAHBT; ~677–774 DCIS, ~481–542 pts,
7+yr follow-up; bulk + spatial RNA + methylation; 812-gene recurrence classifier) is the **RNA +
outcome instrument** Risom's MIBI protein can't be:
1. **Gene-level miR-29 targets** (COL1A1/3A1/4A1/5A1, **LOX/LOXL2**, SPARC, FBN1) across
   normal→DCIS→IBC **and by recurrence outcome** — direct effector readout at RNA resolution; if it
   has **GeoMx epithelial-vs-stromal ROIs**, that is compartment-resolved RNA (complements Risom).
2. **Outcome prediction** — the MH-48 abundance-asymmetry signature + **miR-185-5p** as a DCIS→
   recurrence predictor; the one question the public cohorts (3 events) cannot answer.
**Access reality:** via HTAN Data Portal / Synapse — patient bulk sequencing is likely
dbGaP/EGA-controlled (DAR), but **derived L3/L4 (expression summaries, spatial ROI tables) are
often open** — confirm the specific open files on the portal before committing (heavier lift than
Risom's one-click Mendeley).

### 5.5 miR-29c compartment dissection via HTAN — **✓ DONE 2026-06-26 → MH-52** (S)
MH-51 left miR-29c as cell-intrinsic-capable but **stroma-dominant** in tissue. Resolved on Risom
MIBI (69,672 cells, normal/DCIS/IBC): the answer is **scenario (i) regulatory CAF ACTIVATION, not
(ii) fraction** — total fibroblasts *decrease* but the CAF-activation ratio rises **0.03→0.44**
(JT p=3e-4) with within-fibroblast COLI/FAP/SMA up and **CAF≫NORMFIBRO**; myoepithelial fraction
collapses 0.24→0 (the miR-145/-143 program). Corroborated by the authors' Table_S4 (miR-29 target
collagens CAF-high). The MH-51 "stromal confound" = real desmoplastic CAF biology. **Caveat held:**
MIBI is protein → confirms the miR-29 *targets* + CAF compartment, not miR-29c (see §5.8). Original
scoping retained below for provenance.

- **Datasets (two, complementary; both HTAN, derived data open):**
  - **PRIMARY — Risom et al. 2022 *Cell* (MIBI-TOF single-cell spatial proteomics).** 37-plex,
    **79 matched normal/DCIS/IBC** resections; TME states defined by **myoepithelium + fibroblasts +
    immune**; key result = myoepithelial disruption + stromal transitions at DCIS→IBC. Near-perfect
    for the **myoepithelial-loss (miR-145/-143)** and **CAF/stroma (miR-29c)** programs at
    single-cell + spatial resolution; panel carries SMA/FAP/Vimentin/collagen + myoepithelial
    markers (p63/CK5/calponin). PMID 35063072. **Data location (corrected): Mendeley Data
    `d87vg86zd8` (v3)** — the single-cell table (cells × markers + FlowSOM/manual lineage + tissue
    type) is the L4-equivalent we need. *(NB the Zenodo `10.5281/zenodo.5945388` that a search
    suggested is a DIFFERENT paper — the Liu/Bosse CIMAC MIBI methods set, phenotypes B_cell/DC/
    PAX5 — not Risom DCIS; do not use it.)* **Access caveat:** open but the Mendeley **file API is
    OAuth/SPA-gated** (root listing returns empty, folders 404) → **one manual browser download** of
    the single-cell CSV into `data/external/dcis_mibi/`; no headless pull.
    **Module BUILT (2026-06-26): `dcis_htan_compartment.py`** — schema-adaptive (auto-detects
    lineage/tissue/marker columns), tests fraction-shift (ii) vs within-CAF ECM trend (i) +
    myoepithelial fraction (Kruskal + Jonckheere); prints the download instruction and exits cleanly
    until the CSV is present, then runs end-to-end.
  - **RNA complement — Strand et al. 2022 (TBCRC 038 + RAHBT; 677 DCIS / 481 pts; outcomes).**
    bulk + spatial RNA → gives the **miR-29 RNA target set at gene level** (COL/LOX/SPARC) +
    recurrence outcomes that MIBI's protein panel can't.
- **Levels:** **L3 (expression / per-cell or per-ROI matrices) + L4 (cell-type / ROI annotations)
  are SUFFICIENT — L1/L2 (raw FASTQ / aligned BAM / registered images) are NOT needed** (we consume
  derived expression + labels; no re-alignment or re-segmentation). The open L3/L4 is the right
  layer; controlled dbGaP sequencing (L1/L2) is not required.
- **Portal filters to add to the Level-3/4 filter:** Atlas = *HTAN Duke (Breast PreCancer)*;
  Assay ∈ {spatial transcriptomics (GeoMx DSP — epithelial-vs-stromal ROIs, ideal), scRNA-seq,
  bulk RNA-seq}; Diagnosis ∈ {normal, DCIS, invasive}. Access via HTAN Data Portal / Synapse / CDS.
- **The test (new module `dcis_htan_compartment.py`):** (ii) stromal/CAF fraction across
  normal→DCIS→IDC (and by outcome); (i) within the stromal/CAF compartment, expression of the
  **miR-29 ECM target set** (COL1A1/3A1/4A1/5A1, **LOX/LOXL2**, SPARC, FBN1) across stage — rising
  in CAFs at invasion ⇒ consistent with per-CAF miR-29 loss. GeoMx ROI = region-resolved; spatial
  localises to the invasive front.
- **miRNA caveat (decisive for what to expect):** HTAN profiles **mRNA / protein / spatial, NOT
  miRNA** → this tests the miR-29 **targets/effectors + CAF compartment**, *not* miR-29c itself.
  Per-cell miR-29c would need small-RNA scSeq / miR-ISH (separate modality, not in HTAN).

### 5.6 Stroma-aware compartment extension of the pressure framework (S)
The framework (MH-1..47) is epithelial/tumor-cell-centric; miR-29c→ECM shows the invasion biology
lives in non-epithelial compartments. Build a **myoepithelial/CAF/endothelial-aware** layer; the
**myoepithelial-loss program (miR-145/-143/-126)** is the LATE invasion leader (MH-48) and deserves
its own module rather than being adjusted away.

### 5.7 miR-185-5p — the composition-independent epithelial grade decliner (S)
The one arm that survives composition adjustment as an epithelial grade-decliner (MH-51, adj
ρ=−0.75, q_adj=0.026). Cheap in existing data: its targets, Hallmark links, and whether it is a
bona fide epithelial progression suppressor (vs the stroma-driven miR-29c).

### 5.8 miRNA single-cell / spatial data search (PARKED — for later)
HTAN (5.5) resolves the compartment question only at the **target/effector + CAF** level, never
miR-29c itself (no miRNA modality). To get miR-29c **per-compartment directly**, scout dedicated
small-RNA / miRNA spatial resources: **single-cell small-RNA-seq** (e.g. Smart-seq-total,
half-seq, scsRNA-seq), **miRNA in-situ hybridisation / miRNAscope** on DCIS sections, or
**spatial small-RNA** methods. Deliverable when picked up: a `data/external` candidate list (like
the DCIS GEO corpus) of breast/DCIS miRNA sc/spatial datasets + a feasibility note. This is the
*definitive* answer to "which compartment loses miR-29c"; HTAN is the effector-level proxy until then.
- **Lead evaluated — SiCmiR Atlas** (Cai et al., *Adv. Sci.* 2026; arXiv 2508.05692; code
  github.com/Cristinex/SiCmiR; webpage awi.cuhk.edu.cn/~SiCmiR; 9.36M cells / 362 datasets, breast
  covered). A 2-layer NN predicting **1,298 miRNAs from 977 L1000 genes** (trained on 6,462 TCGA;
  mean Pearson **0.67**). **Usable** as a trained model to apply to a breast/DCIS scRNA dataset for
  per-cell-type *inferred* miR-29c. **BUT a decisive circularity caveat for OUR question:** it
  infers miRNA activity *from target-gene expression*, so "miR-29c-low in CAFs" would be predicted
  *from* "collagen/ECM-high in CAFs" — i.e. it would **restate the effector observation, not
  independently confirm the miRNA**. ⇒ good for hub-miRNA landscape / hypothesis-generation on
  *other* arms and as a cross-check, **not** as causal evidence for the miR-29c→ECM compartment
  claim. True small-RNA sc/spatial (or miR-ISH) remains the only non-circular answer.
- **Landscape search (2026-06-26) — the field is INFERENCE-dominated, all circular for our axis.**
  The public "single-cell miRNA" resources are predictors from mRNA, not measurements: **SiCmiR**
  (above), **STmiR** (XGBoost spatial miRNA activity, bioRxiv 2025.03.18.644021), and **"Inferring
  single-cell and spatial microRNA activity from transcriptomics"** (*Comms Biol* 2025,
  s42003-025-07454-9). All infer miRNA from target expression → **all circular** for miR-29c→ECM.
  **No public *measured* single-cell/spatial small-RNA breast/DCIS dataset was found** (the area is
  emerging; spatial = Xenium/MERFISH mRNA panels, not miRNA). ⇒ the only non-circular routes are
  **(a) compartment-SORTED bulk small-RNA-seq** (FACS/LCM CAF vs epithelial → miRNA-seq; some such
  exist for CAFs generally — worth a targeted GEO scan) or **(b) miR-ISH / miRNAscope** for miR-29c
  on DCIS sections (targeted wet-lab, not a public dataset). Recommend (a) as the realistic
  data-only next step; inference tools stay as landscape/hypothesis only.
- **GEO scan done (2026-06-26) — compartment-sorted miRNA exists, but does NOT cleanly resolve
  miR-29c.** Strong convergent lead: **miR-29b is published DOWN in breast CAFs** (Zhang et al.,
  PMC5503632 — 5 ductal-carcinoma CAF vs NF cultures, α-SMA/FSP1-confirmed, qRT-PCR) → independent,
  *measured* (non-circular) support that the miR-29 family falls in the CAF compartment, upgrading
  MH-52 beyond the target/effector level. Also: **miR-26b down in ER+ breast CAFs** (PMC4030585);
  CAF-vs-NF miRNA profiling studies; and **LCM epithelium-vs-stroma breast miRNA** (~51 DE arms).
  **Gaps:** (1) the miR-29b CAF paper has **no GEO deposit** + is array; (2) it measured **miR-29b,
  not miR-29c** — and MH-51 showed the *paralogs dissociate* (cell-line miR-29b ↑ while miR-29c ↓),
  so miR-29b-down-in-CAF does **not** settle miR-29c; (3) no deposited LCM/sorted **small-RNA-seq**
  resolving miR-29c per compartment surfaced. **Realistic next:** scan GEO directly for a *deposited*
  LCM-stroma-vs-epithelium small-RNA-seq accession (the ~51-DE-miRNA LCM studies), else miR-ISH for
  miR-29c. Net: family-level CAF-down is published; **miR-29c-specific compartment is still the open,
  paralog-resolved question.**
- **RESOLVED 2026-06-26 → MH-54.** Found + tested the deposited compartment dataset **GSE37527** (6
  paired breast CAF vs NF, miRNA array, resolves paralogs) — the non-circular measurement.
  **miR-29c is NOT down in CAFs** (p=0.69; 29a/b flat/up) → CAF desmoplasia (MH-52) is
  **miR-29-independent**; bulk miR-29c signal was composition-correlational, not causal. CAF-lost
  arms = **miR-143/145/126/205/185**. Generalized small-RNA-**seq** scan: deposited human
  small-RNA-seq is **EV/biofluid-dominated**, not sorted tissue → array (GSE37527, GSE162670) is the
  best compartment route. **Remaining open:** (i) the desmoplastic driver if not miR-29 (TGFβ-axis on
  CAF mRNA); (ii) a **DCIS-staged** CAF/epithelial compartment set (GSE37527 is cultured tumour-CAF
  generally) — scan persists for a deposited LCM/sorted small-RNA-seq with stage.

### 5.9 Data-acquisition stance for more small-RNA-seq (assessed 2026-06-26)
**Verdict: NO broad sweep — targeted only.** We already hold 6+ plasma/serum EV miRNA cohorts
(`ev_mirna_multi_cohort`) + the MCF10A EV set; the compartment question is answered by array
(GSE37527/162670). More small-RNA-seq has diminishing returns because the two real gaps are *modality*
(true single-cell/spatial miRNA — doesn't exist for breast/DCIS) and *compartment-with-stage*
(no deposited DCIS-staged sorted/LCM small-RNA-seq), neither of which "more EV/biofluid seq" fills.
**Acquire only if a specific aim is prioritised:** (a) **EV-biomarker aim** → one *independent*
plasma-EV breast cohort with normal/DCIS/invasive or outcome, to test the miR-205/miR-140-3p export
signature as a circulating marker; (b) **desmoplastic-driver aim** → a **CAF mRNA** set (e.g.
GSE215307) to test the TGFβ-axis directly (GSE37527 is miRNA-only). General EV/biofluid trawling is
not worth it.

**Both aims spot-checked with data in hand (2026-06-26):**
- *EV-biomarker (a):* miR-205-5p / miR-140-3p tested across the 6 existing plasma-EV cohorts —
  **inconsistent**: miR-205-5p case-higher only in the one clean cancer-vs-normal cohort (GSE197020,
  AUC 0.71, q=0.009), miR-140-3p only in the screening cohort (GSE270497, AUC 0.72); not replicated
  elsewhere, and miR-29c also elevated despite not being a tissue-export candidate → plasma EV is
  **multi-source**, the tissue-export signature does NOT translate to a robust circulating marker.
  Conclusion: a new EV cohort is **not** worth it for this.
- *Desmoplastic driver (b):* tested single-cell in Risom (no download) — within fibroblasts the
  myCAF/TGFβ markers carry the collagen: **FAP↔COLI ρ=+0.38 (p=2.7e-255)**, SMA↔COLI +0.15, and
  CAF≫NORMFIBRO for COLI/FAP/SMA (all p≪1e-20). ⇒ the desmoplastic ECM is an **activated FAP⁺/αSMA⁺
  myCAF program** (TGFβ-induced state), miR-29-independent — confirms MH-54's "look to TGFβ" without a
  new dataset. (Caveat: FAP/αSMA are TGFβ-myCAF *proxies*; TGFβ signalling not directly measured.)
- *Desmoplastic driver — direct mRNA test (GSE196354, 7 paired breast NF/CAF RNA-seq, `dcis_caf_tgfb_gse196354.py`):*
  **inconclusive due to CULTURE DRIFT** — positive control FAP/ACTA2 only +0.15 (NS), collagens flat,
  TGFβ targets mixed (TGFBI/POSTN borderline up, SERPINE1/CTGF down). Cultured primary CAFs lose the
  in-vivo myCAF/desmoplastic program → **cultured-CAF datasets (this + GSE37527) under-represent the
  in-vivo state**; the TGFβ-desmoplasia conclusion rests on the **in-vivo single-cell Risom**. ⇒ the
  remaining clean test is **in-situ** (spatial/single-cell CAF), not another cultured CAF set. The
  rigor-protocol null+positive controls flagged the unsuitability (apm-rigor-protocol move I).

## Open questions worth a dedicated study

- Is the basal **immune-pressure tone** a *cause* (miRNA-driven IFN silencing) or
  a *consequence* (low infiltration → low immune-gene/miRNA co-program)? §1.1 + a
  TME-deconvolution adjustment would arbitrate.
- Does **8q24 amplification** (miR-151a, AGO2) create a coordinated
  "amplicon-driven RISC + guide" state that amplifies pressure non-linearly?
- Are the **differentiation-program** repressions (estrogen/adipogenesis/myogenesis)
  a generic loss-of-identity signal shared with other carcinomas (pan-cancer test)?
