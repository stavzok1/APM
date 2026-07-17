# miRNA CNV Dosage Report - TCGA-BRCA

Run focus: 2026-06-14 segment-based miRNA CN refresh; **2026-06-17 pressure-model
refresh** (`softmax_z_logrpm` + `evidence_mass` spine; M7 hybrid uses
`softmax_z` + `combined_mass`). This report covers the miRNA-universe copy-number
work plus how **miRNA pressure → target repression** is scored when linking dosage
to target RNA.

Primary outputs live under:

- `mirna_hallmark/output/mirna_locus_cnv/`
- `mirna_hallmark/output/mirna_locus_cnv/maps/`
- `mirna_hallmark/output/ccle_mirna_cn_concordance/` — CCLE cell-line CN→miRNA expr replication (§7)
- `mirna_hallmark/output/target_combined_anticorr/`
- `mirna_hallmark/output/hybrid_pressure/` — M0/M7/M8/M11 Hallmark coupling
- `mirna_hallmark/output/pressure_sensitivity/` — weighting grid vs legacy z
- `data/miRNA/audit/`

The cohort is TCGA-BRCA primary tumors, participant-level collapsed. The current
segment-based cache contains 1,049 primary participants.

---

## 1. What Changed In This Refresh

The miRNA CNV layer was rebuilt around physical miRNA loci rather than gene
proxies. Each MirGeneDB hairpin locus (`MI*`) is assigned copy number from the
best-overlapping ASCAT3 allelic segment. MIMAT arms are then reconstructed from
their contributing hairpin loci.

Important changes:

1. **Segment-based locus CN is now the default.** The previous gene-level ASCAT3
   proxy has been replaced by direct overlap of `cnv_miRNA.csv` hairpin
   intervals against per-sample ASCAT3 segments.
2. **Paralog MIMAT CN is explicit.** A mature arm can come from one or more
   hairpin loci. The MIMAT layer now keeps the locus structure rather than
   pretending every MIMAT is a single genomic site.
3. **Expression weights are cohort-reference weights.** For multi-locus MIMATs,
   `w_norm` is derived from cohort-median locus x MIMAT RPM. It is not
   sample-specific, avoiding circular CN -> expression circularity.
4. **The healthy reference is defined explicitly** (full scheme: §2).
   - Locus: diploid CN = 2.
   - MIMAT: weighted composite on the same scale (diploid composite = 2), not the
     sum of paralog copy numbers.
5. **Subtype summaries are participant-deduped.** The earlier inflated
   stratum counts came from duplicated `(participant, entity_id, entity_level)`
   rows in an old cache. The current outputs report the expected scale:
   approximately 1,049 primary participants, LumA approximately 482.
6. **Target-repression pressure uses the adopted dual model** (2026-06-17).
   **Spine (M0 miRTar):** `softmax_z_logrpm` + `evidence_mass` via
   `pressure_build` — `hallmark_interaction`, `target_combined_anticorr`,
   subtype contrasts, robustness, explainability. **Hybrid extended claims (M7):**
   `softmax_z` + `combined_mass` on tiered miRTar+TS edges (`hybrid_pressure`).
   Legacy `softmax_z+degree` archived in `pressure_sensitivity/`.

Core generated files:

| File | Purpose |
|------|---------|
| `tables/sample_entity_cnv.tsv.gz` | Participant x entity CN for loci, MIMAT arms, and families |
| `mirna_cnv_by_stratum.tsv` | Locus/arm/family CN summaries across strata |
| `mirna_cnv_expr_concordance.tsv` | CN -> mature-arm expression Spearman + partial (CPE, HRD) |
| `mirna_cnv_expr_concordance_scope_audit.tsv` | Full cohort vs CN>2 / gain\|amp sensitivity per arm |
| `figures/mirna_cnv_expr_concordance_rho_boxplot.png` | Boxplot of marginal vs partial ρ across tested arms |
| `concordance_top_arm_target_join.tsv` | miRTarBase target checks for top CN-concordant arms |
| `maps/mirna_cnv_locus_genome_map.tsv` | Locus-level genome map, including PAM50 deltas |
| `maps/mirna_cnv_chromosome_summary.tsv` | Per-chromosome hairpin/mature-arm counts, hg38 size, mean Δ |
| `maps/mirna_cnv_chromosome_arm_summary.tsv` | Per cytoband arm (chrNp/chrNq) counts, arm size, mean Δ |
| `maps/mirna_cnv_locus_cluster_summary.tsv` | Dense hairpin clusters (≥3 loci within 100 kb) |
| `maps/mirna_cnv_cluster_segment_focality.tsv` | Per-cluster ASCAT3 segment span vs cluster window |
| `sv_overlap/mirna_locus_sv_hits.tsv.gz` | Long-form DEL/DUP SV × hairpin locus overlaps |
| `sv_overlap/mirna_locus_sv_recurrence.tsv` | Per-locus SV recurrence (participants with DEL/DUP) |
| `sv_overlap/mirna_cluster_sv_recurrence.tsv` | Dense-cluster SV rollup |
| `sv_overlap/mirna_locus_cnv_sv_layers.tsv.gz` | Participant × locus: cnv_only / sv_only / both / neither |
| `sv_overlap/mirna_locus_sv_both_patterns.tsv` | SV direction concordance within CNV-altered loci |
| `sv_overlap/mirna_locus_sv_only_patterns.tsv` | sv_only focality and SVTYPE mix by stratum |
| `sv_overlap/mirna_locus_sv_both_concordance_by_cn_state.tsv` | both-layer DUP/DEL vs amp/gain/loss/deep_del |
| `sv_overlap/mirna_locus_sv_only_top_loci.tsv` | Top recurrent sv_only loci with CNV context |
| `sv_overlap/review_queues/review_queue_*.tsv.gz` | IGV review queues (focal sv_only, discordant/mixed both) |
| `sv_overlap/review_queues/review_queue_summary.tsv` | Row/participant/locus counts per queue |
| `sv_overlap/chr19_megacluster/chr19_megacluster_loci.tsv` | Per-hairpin stats inside chr19:53.7 megacluster |
| `sv_overlap/chr19_megacluster/chr19_megacluster_participant_locus.tsv.gz` | Participant × locus layers inside megacluster |
| `sv_overlap/chr19_megacluster/chr19_megacluster_summary.json` | Megacluster rollup (layer mix, SV n) |
| `sv_overlap/chr19_megacluster/gallery/index.html` | All-participant chr19 megacluster schematics (filterable) |
| `sv_overlap/chr19_megacluster/gallery/chr19_gallery_manifest.tsv` | One row per participant plot |
| `maps/mirna_locus_dispersion_*.tsv` | Dispersion tables (chromosome, arm, 5 Mb bins, clusters) |
| `maps/figures/mirna_locus_genome_dispersion_multires.png` | Five-panel genome dispersion (genome → chr19 zoom) |
| `maps/figures/mirna_locus_dispersion_by_chromosome_facets.png` | Per-chromosome locus scatter (within-chr Mb) |
| `maps/figures/by_pam50/{PAM50}/mirna_locus_genome_dispersion_multires.png` | Five-panel dispersion using mean locus Δ within that PAM50 stratum |
| `maps/figures/by_pam50/{PAM50}/mirna_locus_dispersion_by_chromosome_facets.png` | Per-chromosome facets by PAM50 stratum |
| `maps/by_pam50/{PAM50}/mirna_locus_dispersion_*.tsv` | Per-subtype dispersion tables (chromosome, arm, 5 Mb bins, clusters) |
| `sv_overlap/review_queues/igv/README_IGV.txt` | IGV desktop load order (one tumor sample per review row) |
| `sv_overlap/review_queues/igv/igv_session_guide.tsv` | Combined session guide: participant, sample_vial, igv_window, SV intervals |
| `sv_overlap/review_queues/igv/mirna_hairpin_loci.bed` | All 506 hairpins as an IGV annotation track |
| `sv_overlap/review_queues/igv/review_queue_sv_intervals.bed` | Queued DEL/DUP intervals (filter by participant prefix in IGV) |
| `sv_overlap/review_queues/case_plots/index.html` | BAM-free schematic gallery (segment + hairpin + SV) |
| `sv_overlap/review_queues/case_plots/case_plot_manifest.tsv` | One row per generated schematic PNG |
| `maps/mirna_cnv_mimat_paralog_map.tsv` | MIMAT x locus long map with `w_norm` and weighted contributions |
| `maps/mirna_cnv_mimat_summary.tsv` | One row per MIMAT weighted CN summary |
| `maps/mirna_cnv_focus_subtype_distribution.tsv` | Per-focus-arm CN distributions by PAM50 |
| `maps/mirna_cnv_neighborhood_by_subtype.tsv` | Same-segment and +/-2 Mb neighborhood co-alteration |
| `target_combined_anticorr/target_combined_anticorr.tsv` | AGO-gated combined miRNA pressure vs target RNA (cohort) |
| `target_combined_anticorr/target_combined_anticorr_by_pam50.tsv` | Within PAM50: raw + CPE/HRD + target-CN partial Spearman |
| `hybrid_pressure/{m0,m7,m8,m11}/hallmark_coupling_by_pam50_prolif.tsv` | Hybrid edge modes × PAM50 (prolif-adjusted Hallmark coupling) |
| `pressure_sensitivity/coupling_summary_by_spec.tsv` | Weighting grid (z / softmax_z / logrpm × degree / evidence_mass / combined_mass) |

---

## 1b. Canonical miRNA Pressure Scoring (Decided 2026-06-17)

This subproject separates **three decisions** that are easy to conflate:

1. **Per-miRNA / per-sample expression weighting** — how much each mature arm
   contributes to pressure on a target gene in tumor *s*.
2. **Edge-weight normalization** — how edge evidence is scaled per miRNA (degree,
   evidence mass, TS mass, combined mass).
3. **Edge universe** — which `(miRNA, gene)` pairs enter the sum (miRTarBase only
   vs TargetScan fill-in vs pure TS).

### A. Per-arm expression weighting — spine (M0 modules)

**Use `softmax_z_logrpm` + `target_norm=evidence_mass` + `aggregate=sum`.**

Configured in `config.PRESSURE_*`; implemented in `pressure_engine.py` /
`pressure_build.py`.

Per gene *g*, sample *s*, targeting arms *m*:

```text
edge_w(m→g) = log1p(evidence) / log1p(Σ_g' log1p(evidence_m→g'))   # evidence_mass
logit_m(s)  = log2(RPM+1)_m,s − cohort_median_m
share_m(s)  = softmax(logit_m) among arms targeting g only
z_m(s)      = cohort z-score of log2(RPM+1) for arm m
pressure(g,s) = Σ_m  share_m(s) · z_m(s) · log2(RPM+1)_m,s · edge_w(m→g)
```

| Component | Role |
|-----------|------|
| **softmax share** | Relative abundance among regulators of *g* |
| **× z** | Down-weights always-on oncomiRs not elevated in *s* |
| **× log2(RPM+1)** | Absolute abundance anchor (avoids proliferation co-expression trap) |
| **evidence_mass** | Down-weights miRNAs with heavy total miRTar evidence mass (vs raw target count) |
| **sum (not mean)** | Mean-over-arms weakens signal |

**Sensitivity winner (M0 edges, Basal prolif-adj):** `softmax_z_logrpm` +
`evidence_mass` → **42/50** neg-sig Hallmarks (8/8 key), vs prior
`softmax_z+degree` **35/50**. Pure `softmax` or `z_softmax` without × z still
flip positive — see `pressure_sensitivity/`.

**Hybrid modules (M7/M8/M11)** use a separate profile: `softmax_z` +
`combined_mass` (`config.PRESSURE_HYBRID_*`) — see §1b.B–D.

**Abundance floor:** TS-only / orphan arms require cohort-median log2(RPM+1) ≥
**1.0** (`PRESSURE_TS_ABUNDANCE_FLOOR`); miRTarBase arms use floor **0.0**.

**Legacy reference:** `softmax_z+degree` and cohort-z-only sums archived under
`pressure_sensitivity/`.

### B. Edge universe — miRTarBase vs TargetScan (decided roles)

| Mode | Code | Edge set | Use for |
|------|------|----------|---------|
| **Spine / headline** | **M0** | miRTarBase, evidence ≥ 2 | Targeting map, `target_combined_anticorr`, subtype contrasts, OLS explainability, **this report §8** — weighting: **logrpm + evidence_mass** |
| **Hybrid primary** | **M7** | HE edges TS-boosted + capped orphan TS fill-in | Extended Hallmark coupling + exploratory IRF1 hub partial — weighting: **softmax_z + combined_mass** |
| **Hybrid sensitivity** | **M8** | `P_mirtar + λ·P_orphan` (separate tracks) | Robustness (39/50 Basal); 7/8 key programs |
| **TS-only lane** | **M11** | All TS pairs ≥ 0.25, top-k/gene, TS floor | Curation-gap / sequence-only sensitivity; not headline |
| **Orphan coupling** | `targetscan_orphan` | TS pairs **not** in miRTarBase HE | Hub-level curation gaps (miR-130b/301 → IRF1); separate from M0 pressure |

**Decided policy:** report **M0 logrpm+evidence_mass** for dosage→repression links in this
document; cite **M7 combined_mass** when arguing TS-supported extension or IRF1 hub
partial; **M8/M11** as robustness only. Do not merge TS orphans into M0 without caps.

### C. Mature-arm name harmonization (miRTar → GDC expression)

Before pressure is computed, miRTarBase arm strings are mapped onto GDC miRNA-seq
row names (`mirna_arm_resolve.py`; full policy in **`docs/METHODS.md` §2b**).

**Applied to:** combined pressure, hybrid modes, `target_combined_anticorr`, hub OLS
explainability — any module that scores **miRTar edges × arm abundance**.

**Not applied to:** primary **CN→miRNA expression** concordance
(`mirna_cnv_expr_concordance.tsv`), which joins the miRNA entity catalog directly
to the expression matrix. CCLE replication uses its own alias map.

**Default rules (Tier A+B only):** isoform suffixes (`.1`/`.2`/`.3`), loci/MIMAT
map, star strip; precursor `-3p`/`-5p` append only when unambiguous. **Excluded:**
opposite-arm flip (`3p`↔`5p`), letter slim, sibling collapse — different molecules
or paralog risk.

**Impact on Hallmark M0 edges:** 245/730 previously dropped rows recovered; audit:
`output/logs/mirna_arm_resolve_audit.tsv`.

### D. Proliferation-adjusted Hallmark coupling (Basal, M0 vs hybrids)

CPE + HRD + E2F/G2M + mean target CN; within-subtype z on pressure & expression.
Source: `hybrid_pressure/` refresh 2026-06-17 (prolif-adj; within-subtype z).

| Mode | Weighting | Basal neg-sig Hallmarks | Key programs neg-sig |
|------|-----------|------------------------:|-------------------:|
| M0 (miRTar spine) | logrpm + evidence_mass | **42/50** | **8/8** |
| M7 (tiered fusion) | softmax_z + combined_mass | **41/50** | **8/8** |
| M8 (sum tracks) | softmax_z + combined_mass | **39/50** | 7/8 |
| M11 (pure TS capped) | softmax_z + combined_mass | **33/50** | 6/8 |
| Legacy softmax_z+degree (sensitivity) | degree norm | 35/50 | 8/8 |
| Legacy z-only (sensitivity) | none | 27/50 | 7/8 |

### E. Where each module gets its pressure

| Analysis | Edge universe | Expression weighting |
|----------|---------------|---------------------|
| `target_combined_anticorr` (§8) | M0 miRTarBase | **logrpm + evidence_mass** |
| `hallmark_interaction` | M0 | same (spine) |
| `robustness_checks` / hub routes | M0 (+ family filter) | same (spine) |
| `hybrid_pressure` M0 baseline | M0 | spine |
| `hybrid_pressure` M7/M8/M11 | mode-specific | **softmax_z + combined_mass** |
| `targetscan_orphan_coupling` | TS orphan only | spine + TS abundance floor |
| `gene_expression_explainability` | Route-specific edges | spine (route pressure term) |

Formal methods: `docs/METHODS.md` §3; registry: `DISCOVERY_REGISTRY.md` MH-30.

---

## 2. Healthy Reference, MIMAT Weighting, And Interpretation

### 2.1 Hairpin locus (MI*) — direct diploid reference

ASCAT3 reports **total copy number per segment**. For a single MirGeneDB hairpin
locus, the theoretical healthy reference is **diploid CN = 2** (one copy per
homolog on autosomes):

```text
locus_delta = observed_locus_CN − 2
```

This is the primary layer: segment overlap → `copy_number` per `MI*` per
participant (`tables/sample_entity_cnv.tsv.gz`, `entity_level = locus`).

### 2.2 MIMAT arm — a weighted composite, not a genomic interval

A **MIMAT** (mature arm) has no single genomic CN. It is a mature product that
can be produced from **one or more paralog hairpin loci** (`mature_accessions` on
`cnv_miRNA.csv` → MIMAT map).

**Important:** MIMAT CN in this report is a **synthetic dosage index** on the
**same scale as per-locus ASCAT CN** (diploid ≈ 2). It is **not**:

- the **sum** of paralog copy numbers (three diploid loci would be 6, not 2);
- measured mature-RNA abundance (RPM);
- a claim that one abstract “MIMAT locus” is diploid.

It answers: *“Given how this arm is usually split across paralogs in the
cohort, what is the arm-level copy-number dosage index for this tumor?”*

Single-locus arms (e.g. miR-21) are a special case: one `MI*`, `w_norm = 1`,
MIMAT CN = locus CN, reference = 2.

### 2.3 Weight construction (`w_norm`)

Weights are **cohort-reference** — identical for every participant and tumor
state. They are **not** re-estimated from each sample’s expression or CN (that
would circularly couple CN→expression concordance tests to the same readout).

**Step 1 — locus × MIMAT expression reference (`weight_ref`):**

From GDC miRNA-seq, build per-sample RPM for each `(locus_id, mimat)` pair
(`data/miRNA/locus_mimat_rpm.tsv.gz`; hairpin coords mapped via
`locus_gdc_hairpin_map.tsv`). Take the **cohort median RPM** per pair →
`weight_ref` (`build_locus_mimat_weight_reference` in
`pipeline/miRNA_exp/locus_mimat_expression.py`).

**Step 2 — normalize within each MIMAT (paralog shares):**

```text
w_norm(locus | mimat) = weight_ref(locus, mimat) / Σ weight_ref(all paralogs of that mimat)
```

So `Σ w_norm = 1` per MIMAT. `w_norm` is the cohort-typical **fraction of that
mature arm’s expression pool** attributed to each hairpin locus — a fixed prior
on which paralog’s CN move should matter most for arm-level dosage.

Stored in `maps/mirna_cnv_mimat_paralog_map.tsv` (`w_norm`, `weight_ref_rpm`).

### 2.4 Tumor MIMAT CN (per participant)

For participant `p` and MIMAT `m`, with paralog loci `i` and observed locus CN
`CN_i(p)` from ASCAT3:

```text
MIMAT_CN(p) = Σ_i  w_norm(i | m) × CN_i(p)
```

Implemented in `mirna_cnv_genome_maps.py`, `mirna_cnv_subtype_depth.py`, and
`dosage_landscape_cnv.aggregate_mimat_cnv_from_loci` when
`mimat_cn_aggregate = expression_weighted` (focus-arm and paralog tables use this
scheme). If all weights for a MIMAT are zero/missing, fall back to the median of
paralog locus CNs.

**Rationale:** If the highly expressed paralog sits on an amplicon and a minor
paralog stays diploid, the arm-level index should move with the **expressed**
locus, not with an unweighted median (which can hide focal amp on the dominant
paralog). See flagged “paralog-spread” cases in §6 (miR-9, miR-320b, …).

### 2.5 Healthy reference for MIMAT (why 2, not “more”)

Theoretical healthy means **every contributing locus is diploid** (`CN_i = 2`).
Apply the **same** `w_norm` as for tumors:

```text
healthy_MIMAT_CN = Σ_i  w_norm(i | m) × 2  =  2 × Σ_i w_norm  =  2
```

So the reference is **2 by definition** of this index: when all paralogs are
diploid, the weighted composite equals single-locus diploid scale. It does **not**
mean “only two DNA copies exist genome-wide” — a three-locus arm still has
3 × 2 = **6 physical precursor gene copies** in the genome. The index **averages**
(with expression weights), it does not **sum**, so arms remain comparable to
single-locus MIMATs (miR-21 diploid ≈ 2, let-7a diploid ≈ 2).

**Tumor shift vs healthy:**

```text
MIMAT_delta(p) = MIMAT_CN(p) − healthy_MIMAT_CN
               = Σ_i  w_norm(i) × (CN_i(p) − 2)
```

Healthy and tumor use **identical weights**; the reference is not a separate
weighting scheme.

### 2.6 Paralog decomposition columns (cohort summaries)

`maps/mirna_cnv_mimat_paralog_map.tsv` exposes the accounting per locus (example:
let-7a-5p in §5.5):

| Column | Meaning |
|--------|---------|
| `w_norm` | Paralog share of mature-arm pool (§2.3) |
| `healthy_cn_contribution` | `2 × w_norm` — that locus’s share of the **reference composite 2** |
| `mean_locus_cn` | Cohort-mean ASCAT CN at that `MI*` |
| `delta_locus_vs_diploid` | `mean_locus_cn − 2` |
| `weighted_delta_contribution` | `w_norm × (mean_locus_cn − 2)` — locus’s contribution to cohort MIMAT delta |

Summing `healthy_cn_contribution` over paralogs = **2**. Summing
`weighted_delta_contribution` = cohort `mimat_delta_vs_healthy_weighted`.

These columns are for **interpretation and QC**, not a second reference model.

### 2.7 Alternatives considered (not the default here)

| Aggregation | Healthy ref | Issue for multi-paralog arms |
|-------------|------------|------------------------------|
| **Expression-weighted mean** (default) | 2 | Keeps per-arm scale; weights dose by expressed paralog |
| Unweighted **median** of paralog CNs | 2 (if all diploid) | Can hide amp on the main-expressed locus |
| **Sum** of paralog CNs | 2 × n_loci | Different scale per arm; not comparable to single-locus MIMATs |
| **Max** paralog CN | 2 | Pessimistic / amp-driven; ignores minor diploid paralogs |

### 2.8 Caveats

1. **Linear CN → dosage proxy** — processing, stability, and post-transcriptional
   regulation are not modeled; weights are a pragmatic prior, not fitted per tumor.
2. **Fixed cohort weights** — a tumor with unusual paralog usage still uses
   cohort `w_norm`; per-sample weights would confound CN→expression tests.
3. **Unmapped hairpins** — loci with zero GDC hairpin overlap get `weight_ref = 0`
   and do not contribute to weighted MIMAT CN (locus-level CN is still reported).
4. **Composite ≠ total paralog burden** — for “how many precursor copies exist
   across the genome for this arm?”, read **per-locus CN** or sum paralog CNs
   explicitly; do not read MIMAT CN = 2 as that sum.

### 2.9 Empirical baseline (Normal-like)

The empirical PAM50 **Normal-like** group is useful as an in-cohort anchor, but
it is not the same as healthy normal breast tissue. Normal-like samples still
show mean miRNA CN above 2 and approximately 24% altered participant-arm
observations. Therefore the report uses both references:

- **Theoretical healthy:** diploid composite = 2 (§2.4–2.5).
- **Empirical low-alteration baseline:** PAM50 Normal-like.

Formal methods cross-ref: `docs/METHODS.md` §8; implementation:
`mirna_cnv_genome_maps._load_weight_reference`,
`mirna_cnv_subtype_depth._participant_mimat_cn`.

---

## 3. Global miRNA CN Landscape

Across the mature-arm layer, the miRNA CN landscape is strongly gain-biased.

| Quantity | Value |
|----------|------:|
| Primary participants | 1,049 |
| Mature arms in current concordance table | 756 tested |
| MIMATs in genome map | 605 |
| MIMATs with >1 paralog locus | 36 |
| Hairpin loci in genome map | 506 |
| Mean locus delta vs diploid | +0.859 |
| Hairpin loci with mean delta > +1 | 138 / 506 |

The broad pattern is not loss of miRNA loci. It is amplification and arm-level
copy gain. Cohort-wide, no MIMAT arms have mean CN below 2 in the segment-based
summary. The strongest chromosome-level locus deltas are (hairpin loci vs
unique mature arms; chromosome sizes are hg38):

| Chromosome | Chr size (Mb) | Hairpin loci | Mature arms | Mean locus Δ vs CN=2 | Loci with Δ > +1 |
|------------|-------------:|-------------:|------------:|---------------------:|-----------------:|
| chr8 | 145 | 15 | 15 | +1.711 | 10 / 15 |
| chr1 | 249 | 35 | 51 | +1.593 | 20 / 35 |
| chr20 | 64 | 9 | 4 | +1.534 | 9 / 9 |
| chr19 | 59 | 66 | 31 | +1.021 | 53 / 66 |
| chr7 | 159 | 25 | 42 | +1.003 | 9 / 25 |
| chr16 | 90 | 11 | 9 | +0.897 | 6 / 11 |
| chr12 | 133 | 12 | 23 | +0.875 | 2 / 12 |
| chr21 | 47 | 5 | 8 | +0.861 | 0 / 5 |
| chr17 | 83 | 26 | 43 | +0.837 | 13 / 26 |

**Arm-level split (same loci, cytoband p/q).** Gain is concentrated on long q
arms and chr1q/chr8q rather than spread uniformly:

| Cytoband arm | Arm size (Mb) | Hairpin loci | Mature arms | Mean locus Δ vs CN=2 |
|--------------|-------------:|-------------:|------------:|---------------------:|
| chr8q | 100 | 9 | 10 | +2.526 |
| chr1q | 127 | 20 | 30 | +2.352 |
| chr20q | 38 | 8 | 2 | +1.590 |
| chr16p | 36 | 6 | 5 | +1.537 |
| chr17q | 58 | 16 | 25 | +1.261 |
| chr19q | 32 | 57 | 15 | +1.050 |
| chr7p | 60 | 8 | 11 | +1.097 |
| chr8p | 45 | 6 | 5 | +0.489 |
| chr1p | 122 | 15 | 22 | +0.580 |

For contrast, chr8 carries only **15 hairpin loci** on a **145 Mb** chromosome,
but **chr8q alone** (100 Mb) accounts for the entire chr8 signal (mean Δ +2.53
vs +0.49 on chr8p). chr1 shows the same q-arm skew (+2.35 on chr1q vs +0.58 on
chr1p). chr19 is the opposite pattern: **66 hairpins** packed into **59 Mb**,
with a dense **53.7 Mb cluster** of 46 loci spanning ~96 kb (mean Δ +1.06) —
a polycistronic miRNA cluster, not a megabase amplicon.

**Segment focality vs dense clusters.** For each dense cluster (≥3 hairpins
within 100 kb), we asked whether the ASCAT3 segment carrying amp/del is confined
to a small window around the cluster, using per-participant `overlap_segment`
labels from the segment-based CNV cache. Focality thresholds: segment overhang
beyond the cluster interval ≤ **5 Mb**; or total segment span ≤ **10 Mb**; or
segment span ≤ **10×** the cluster genomic span (minimum 500 kb).

| Cluster (chr:Mb) | Hairpins | Cluster span | Same segment | Median segment | Median overhang | Focal (≤10× cluster) | Focal amp (≤10×) |
|------------------|---------:|-------------:|-------------:|-------------:|----------------:|---------------------:|-----------------:|
| chr19:53.7–53.8 | 46 | 96 kb | 99.6% | 58 Mb | 58 Mb | 0.7% | 1.1% |
| chr19:53.8 | 3 | 1 kb | 100% | 58 Mb | 58 Mb | 0.2% | 0.4% |
| chr7:100.1 (miR-17~92) | 3 | 0.5 kb | 100% | 136 Mb | 136 Mb | 0% | 0% |
| chr19:13.8 | 5 | 39 kb | 99.7% | 58 Mb | 58 Mb | 0.4% | 0.8% |
| chr21:16.5 | 3 | 51 kb | 100% | 36 Mb | 36 Mb | 0.1% | 0.2% |

Read: cluster hairpins are almost always on **one shared ASCAT3 segment**
(>99%), and that segment **contains** the cluster interval (>99%). But the
segment itself is **not** micro-focal — median spans are **36–136 Mb** with
overhang comparable to full segment length. Only **≤1%** of participants (and
≤1% of amp/gain participants) have a segment within 10× the cluster span. The
chr19 polycistronic block co-amplifies because hairpins sit on a **chr19-scale
ASCAT3 segment**, not because amp/del is confined to the ~96 kb cluster window.
The miR-17~92 cluster on chr7 behaves the same way (median segment **136 Mb**).
Scattered loci on chr8 (no dense cluster) share one segment in only **~30%** of
participants (median segment **68 Mb**), consistent with arm-level co-amp rather
than a focal cluster lesion.

Canonical tables: `maps/mirna_cnv_cluster_segment_focality.tsv`,
`maps/mirna_cnv_chromosome_summary.tsv`, `maps/mirna_cnv_chromosome_arm_summary.tsv`,
`maps/mirna_cnv_locus_cluster_summary.tsv`.

**SV DUP/DEL at miRNA loci (structural vs segment CNV).** The SV strict-set
(930 samples, `07_final_sv_with_fimo`) was scanned for **DEL/DUP** with direct
pre-miRNA overlap (`mir_hits` ``region_hit=overlaps`` or ``overlap_bp>0``). All
**506** hairpin loci receive at least one overlapping DEL/DUP in **828**
participants (163,842 hit rows). This is common at the locus level, but it is
**not** the primary dosage layer relative to ASCAT3 segment CNV.

Participant × locus layer assignment (diploid + no SV = ``neither``; ASCAT3
``cn_state`` gain/amp/loss/deep_del = CNV altered):

| Layer | Participant × locus rows | Share |
|-------|------------------------:|------:|
| cnv_only | 253,852 | 49.8% |
| neither | 246,912 | 48.5% |
| both | 20,667 | 4.1% |
| sv_only | 9,363 | 1.8% |

Among the **138 hairpin loci** with mean CNV Δ > +1, **sv_only** is rare
(2,129 rows, **1.5%**); **cnv_only** dominates (78,231 rows). The chr19
46-locus polycistronic cluster illustrates the split: **528** participants carry
CNV amp/gain at those loci, but only **42** have an overlapping SV DEL/DUP;
within-cluster layers are **cnv_only 50%**, **neither 46%**, **both 2.8%**,
**sv_only 0.9%**. The chr7 miR-17~92 cluster (3 hairpins) has **100**
participants with SV overlap but the same CNV-heavy layer mix.

Top recurrent SV loci (participants with any overlapping DEL/DUP) include
miR-139 (184), miR-486 (183), miR-194/192 paralogs on chr11, miR-21 (169), and
chr1q amplicon neighbors (miR-214, miR-199, miR-1295). Dense clusters with the
most SV participants: chr7:100.1 miR-17~92 (100), chr13:50 (81), chr19:13.8
(79), chr4:112.6 (78) — not the chr19:53.7 megacluster (42).

**What ``both`` looks like (SV inside CNV-altered loci).** At participant × locus
grain (**20,667** rows), **77%** carry at least one **direction-concordant** SV
(DUP on amp/gain, DEL on loss/deep_del). Concordance is highest on **loss/deep_del**
(89–96%) and lowest on **gain** (68%). **37%** of ``both`` rows have **both DUP and
DEL** overlapping the same hairpin in the same sample — complex rearrangement on
an already-altered background, not a clean single SV call. Her2/LumB show the most
mixed DUP+DEL (**~43%**). Median overlapping SV span is **~29 Mb** (amp **~26 Mb**,
gain/loss **~32 Mb**), i.e. the SV layer inside CNV is still **arm-scale**, not a
focal lesion nested inside a focal CNV.

| CNV state | n (participant × locus) | Any concordant SV | Discordant SV only | Mixed DUP+DEL | Median SV span |
|-----------|------------------------:|------------------:|-------------------:|--------------:|---------------:|
| amp | 11,475 | 78% | 22% | 40% | 26 Mb |
| gain | 5,582 | 68% | 32% | 34% | 32 Mb |
| loss | 3,529 | 89% | 11% | 30% | 33 Mb |
| deep_del | 81 | 96% | 4% | 10% | 9 Mb |

**What ``sv_only`` looks like (SV without ASCAT3 change).** All **9,363** rows are
strictly **CN=2 / neutral** — by construction. This is **not** the high-Δ amplicon
story: top recurrent loci are chr11 paralogs (miR-192/194, **71** participants each;
cohort mean Δ **+0.75**) and chr3 let-7/miR-30 neighbors (mean Δ **+0.46–0.47**).
**Her2** has the highest ``sv_only`` rate (**3.2%** of participant × locus rows vs
**1.1%** in LumA). Overlapping SVs are **not focal**: median span **~36 Mb**; only
**2.5%** are ≤1 Mb. SVTYPE mix is heterogeneous (**41%** DEL-only, **28%** DUP-only,
**31%** both DEL and DUP at the same locus). Median tumor alt support **~39** —
these are real calls, but they are typically **large-span SVs overlapping a diploid
segment**, consistent with subclonal structural events or ASCAT3 smoothing rather
than a pure miRNA micro-deletion/duplication invisible to segments. The chr19:53.7
megacluster contributes **427** ``sv_only`` rows (median SV span **~7 Mb** — smaller
than the cohort median but still megabase-scale).

**Interpretation for later adjudication:** overlapping SV DEL/DUP at miRNA loci
is widespread but usually **redundant with segment CNV** (``both`` or
``cnv_only``). In ``both``, SV direction **usually matches** CNV state (~77%) but
**mixed DUP+DEL** is common (~37%) on arm-scale events. The **sv_only** tail
(~1.8%; ~1.5% at high-Δ loci) is enriched in **moderate-Δ paralog clusters** and
**Her2**, with **megabase SV spans on CN=2 segments** — prioritize for
subclonal/SV-vs-segment review, not for treating as independent focal miRNA amp/del.
The natural next step is a merged dosage call that down-weights redundant ``both``
events and flags ``sv_only`` for manual or purity-aware review.

**IGV review queues (exported).** Three pre-filtered tables under
``sv_overlap/review_queues/`` — each row is one participant × locus with
``igv_locus``, CNV state, overlapping SV ids/types, and queue-specific flags:

| Queue | Rows | Participants | Loci | Selection rule |
|-------|-----:|-------------:|-----:|----------------|
| ``review_queue_sv_only_focal_le_1mb`` | 235 | 112 | 151 | ``sv_only`` and median overlapping SV span ≤1 Mb |
| ``review_queue_both_discordant_sv_only`` | 4,674 | 522 | 461 | ``both`` with **no** direction-concordant DEL/DUP |
| ``review_queue_both_mixed_dup_del`` | 7,548 | 580 | 502 | ``both`` with both DEL and DUP at the same locus |

The focal ``sv_only`` queue is the smallest and most IGV-friendly (~2.5% of all
``sv_only``). The ``both`` queues are larger because arm-scale mixed events are
common; sort by ``median_sv_span_mb`` or ``PAM50_final`` before manual review.

Regenerate queues only (no SV rescan):
``.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_overlap --queues-only``

**SV VAF / clonality columns (Manta + CPE).** Overlap hits and review queues carry
tumor read VAF from Manta strict-set CSVs and purity-adjusted CCF from ASCAT3 ``CPE``:

| Column | Meaning |
|--------|---------|
| ``tumor_vaf`` | Manta read VAF: ``tumor_alt / (TUMOR_PR_ref + TUMOR_SR_ref + tumor_alt)`` |
| ``tumor_vaf_adj`` | CPE-scaled cancer-cell fraction (CCF) via ``purity_adjusted_vaf`` |
| ``tumor_vaf_expected_locus`` | Expected VAF from locus ``copy_number`` + ``CPE`` (not SV midpoint segment) |
| ``tumor_vaf_ratio_adj`` | ``tumor_vaf_adj / tumor_vaf_expected_locus`` — subclonal vs locus-expected |
| ``median_tumor_vaf`` / ``median_vaf_ratio_adj`` | Per participant×locus medians in review queues |

Case schematics label each SV bar with ``VAF=… · CCF=… · ratio=…``; ASCAT3 segment
track shows ``CPE=…``. Refresh without full SV rescan:
``.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_overlap --enrich-vaf-only``
(rebuilds hits VAF, review queues, case plots, chr19 gallery).

**chr19:53.7 megacluster dossier.** Densest hairpin block in the panel:
**46** loci in **96 kb** (``chr19:53.667-53.762``), cohort mean Δ **+1.06**.
Layer mix across 48,254 participant × locus rows: **cnv_only 50%**, **neither 46%**,
**both 3%**, **sv_only 0.9%** (427 rows). Only **42** participants have any SV
overlap here vs **528** with CNV amp/gain on these loci — SV is a thin tail on a
predominantly segment-amplified polycistronic block. Per-hairpin breakdown:
``chr19_megacluster/chr19_megacluster_loci.tsv``; full participant table:
``chr19_megacluster_participant_locus.tsv.gz``.

**chr19 megacluster HTML gallery (all 1,049 participants).** One schematic per participant
over the concentrated window (~256 kb around the 46-locus block). Hairpin ticks are colored
by dosage layer; SV DUP/DEL and ASCAT3 segment included. Filter by layer / PAM50 / SV in
``chr19_megacluster/gallery/index.html``. Regenerate:
``python -m mirna_hallmark.analyses.cnv_locus.mirna_locus_chr19_megacluster_gallery`` or ``--chr19-gallery-only``.
Open via local server (same as review case plots):
``cd …/chr19_megacluster/gallery && python3 -m http.server 8766`` →
``http://localhost:8766/index.html`` in Cursor Simple Browser.

**Genome dispersion (multi-resolution).** Figures and tables summarize where the
506 hairpins sit relative to aneuploidy:

- **Genome-wide Manhattan** — all loci on a concatenated hg38 x-axis, colored by
  mean Δ vs diploid.
- **Chromosome bar chart** — locus counts per chr (chr19 leads with 59 hairpins).
- **Top cytoband arms** — arms ranked by locus density (chr19q, chr1q, …).
- **5 Mb bins** — fixed-width bins highlight local pile-ups vs sparse arms.
- **chr19 zoom (53.4–54.0 Mb)** — megacluster band vs flanking loci.

Outputs: ``maps/figures/mirna_locus_genome_dispersion_multires.png``,
``mirna_locus_dispersion_by_chromosome_facets.png``, and
``maps/mirna_locus_dispersion_{chromosome,cytoband_arm,genome_bins,dense_clusters}.tsv``.

Regenerate: ``.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_genome_dispersion``

**Per-PAM50 dispersion.** The cohort-wide figures use all-participant mean Δ. For
subtype-specific **mean locus Δ within each PAM50 stratum** (same 506 loci, different
vertical scale — e.g. Her2/LumB sit higher), see:

- Figures: ``maps/figures/by_pam50/{Normal,LumA,LumB,Basal,Her2}/``
- Tables: ``maps/by_pam50/{PAM50}/mirna_locus_dispersion_*.tsv``

Locus **counts** per chromosome/arm/bin are unchanged across subtypes; only the mean Δ
(and thus Manhattan/facet **heights**) reflect subtype-specific aneuploidy.

**IGV desktop review (per sample, optional BAM).** Each review-queue row is **one participant × one
hairpin** — IGV with BAM/CRAM is optional. **Without BAM**, use:

1. **Case schematic gallery** — ``review_queues/case_plots/index.html`` (open in a browser).
   Each PNG shows ASCAT3 segment (gray), hairpin (red), and overlapping DEL/DUP (blue/red).
   Each SV bar is labeled with span, **VAF / CCF / ratio**, and ``sv_id``; segment track shows **CPE**.
   Default: all 235 focal ``sv_only`` cases + top 40 smallest-span cases per ``both_*`` queue.
   Regenerate: ``python -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_case_plots`` or ``--case-plots-only``.
2. **UCSC Genome Browser** — ``ucsc_url`` column in ``igv_session_guide.tsv`` (hg38, no login).
   Add ``mirna_hairpin_loci.bed`` / ``review_queue_sv_intervals.bed`` as custom tracks.
3. **Tables** — ``sv_intervals``, ``overlap_segment``, ``cn_state``, ``dosage_layer`` in the queue TSVs.

With BAM/CRAM (IGV):
2. Pick a row; note ``participant`` (case ID, e.g. ``TCGA-3C-AALI``), ``sample_vial`` (tumor
   aliquot barcode when present), and ``igv_window`` (recommended viewport).
3. In IGV: **File → Load from URL / local** → primary tumor **BAM/CRAM** for that case
   (GDC Portal → search ``participant`` → open the ``-01A`` WGS/RNA aliquot).
4. Load annotation tracks from ``review_queues/igv/``:
   - ``mirna_hairpin_loci.bed`` — hairpin intervals (``igv_locus`` is the tight locus).
   - ``review_queue_sv_intervals.bed`` — DEL (red) / DUP (blue); filter to rows whose name
     starts with the current ``participant``.
5. Optional: ASCAT3 segment track or ``overlap_segment`` from the table; per-sample SV CSV
   under ``data/SV/pipeline_output/07_final_sv_with_fimo/{sample_vial}_strict_sv_set.csv``
   (match ``sv_ids`` in the row).
6. Navigate: paste ``igv_window`` into the IGV locus box (padded around the hairpin and
   overlapping SV span).

Full step-by-step: ``sv_overlap/review_queues/igv/README_IGV.txt``. Rebuild after queue
refresh: ``--queues-only`` (includes IGV bundle).

Regenerate full SV scan: `.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_overlap --workers 8`

Biological read: the miRNA layer mirrors breast-cancer aneuploidy and focal
amplicon structure. The important question is not "are miRNAs altered?" but
"which mature arms are carried by subtype-specific amplification blocks, and do
those dosage shifts propagate to expression or target repression?"

---

## 4. PAM50 Copy-Number Hierarchy

At the mature-arm level, the average CN hierarchy is:

```text
LumB ~= Her2 > Basal > LumA >> Normal-like
```

Representative participant-deduped arm-level summary:

| PAM50 | n participants | Mean CN | Median CN | Mean delta vs 2 | % diploid | % altered |
|-------|---------------:|--------:|----------:|----------------:|----------:|----------:|
| Normal | 35 | 2.314 | 2 | +0.314 | 75.3 | 24.3 |
| LumA | 482 | 2.718 | 2 | +0.718 | 55.6 | 43.5 |
| Basal | 190 | 2.738 | 2 | +0.738 | 41.5 | 57.3 |
| Her2 | 90 | 3.095 | 3 | +1.095 | 38.8 | 60.0 |
| LumB | 250 | 3.122 | 3 | +1.122 | 38.1 | 60.8 |

Interpretation:

- **Normal-like** remains close to diploid but is not perfectly neutral.
- **LumA** has many near-diploid samples with a long amplicon tail.
- **LumB and Her2** are the major miRNA-amplification states.
- **Basal** has substantial gain burden but its strongest biology is more
  selective, especially around the miR-200 family.

---

## 5. Case Presentations

**Median CN [95% CI]** in the focus-arm tables below: one expression-weighted
MIMAT CN per participant (cohort `w_norm` for multi-locus arms; §2); within each
PAM50 stratum, a percentile bootstrap CI on the participant median (2,000
resamples, all finite participants retained, equal weight per participant — no
outlier trimming and no CN-weighting inside the bootstrap;
`mirna_cnv_subtype_depth.py` → `maps/mirna_cnv_focus_subtype_distribution.tsv`).
Read **Range**, **% amp**, and **% diploid** alongside the CI when a diploid mass
and amp tail coexist; cross-subtype claims use Kruskal–Wallis / Mann–Whitney, not
CI overlap.

### 5.1 miR-21 - Her2/LumB High-Level Amplification

The miR-21 hairpin is one of the clearest cases where the subtype distribution
is not subtle. The median shifts from diploid in Normal-like and LumA to CN=4 in
LumB and Her2, with high-level tails up to CN 37-40.

| PAM50 | n | Median CN [95% CI] | Range | % amp | % diploid |
|-------|--:|--------------------|-------|------:|----------:|
| Normal | 35 | 2 [2, 2] | 1-4 | 6 | 83 |
| LumA | 482 | 2 [2, 3] | 1-34 | 35 | 50 |
| LumB | 250 | 4 [4, 4] | 1-37 | 59 | 25 |
| Basal | 189 | 3 [2, 3] | 1-9 | 25 | 36 |
| Her2 | 90 | 4 [3, 4] | 1-40 | 57 | 24 |

Cross-subtype evidence is strong (Kruskal-Wallis p approximately 9.3e-23).
LumB vs LumA has median difference +2 CN copies with FDR-significant
Mann-Whitney support.

Neighborhood read: miR-21 is not a hairpin-isolated event. In LumB, same-segment
neighbors such as miR-130, miR-142, miR-4524, and miR-3615 co-amplify with
miR-21 in approximately 51-57% of samples, with segment-level CN correlations
around rho 0.73-0.92. The typical anchor segment is large, approximately 83 Mb,
so this is an arm-level 17q amplification context.

Interpretation: miR-21 dosage gain is robust in Her2/LumB. However, dosage gain
does not automatically prove target repression in bulk RNA. The target-resolved
join should be read with the repressor prior: miRNA CN/expression should
anti-correlate with target RNA if repression is visible at bulk level.

### 5.2 chr8 miRNA Amplicon - miR-151, miR-4661, miR-30d/b

The chr8 block is the strongest genome-wide miRNA CN neighborhood. It contains
miR-151, miR-4661, miR-4662, miR-4664, miR-30d/b paralogs, miR-599, and miR-875
among other loci.

Representative distributions:

| Arm | PAM50 peak | Median at peak | Range at peak | % amp at peak | KW p |
|-----|------------|---------------:|---------------|--------------:|-----:|
| miR-151a-3p | LumB | 5 | 1-19 | 66 | 3.6e-19 |
| miR-4661-5p | LumB | 5 | 1-26 | 70 | 3.0e-20 |
| miR-30d-5p | LumB/Basal | 5 | 1-22 / 1-15 | 69 / 63 | 2.5e-21 |

miR-151a-3p by subtype:

| PAM50 | n | Median CN [95% CI] | Range | % amp | % diploid |
|-------|--:|--------------------|-------|------:|----------:|
| Normal | 35 | 2 [2, 2] | 2-7 | 26 | 69 |
| LumA | 482 | 3 [3, 3] | 1-20 | 43 | 42 |
| LumB | 250 | 5 [4, 5] | 1-19 | 66 | 19 |
| Basal | 190 | 4 [4, 5] | 1-14 | 63 | 18 |
| Her2 | 90 | 4 [3, 4] | 1-13 | 57 | 29 |

Neighborhood read: the chr8 events are highly synchronized. In LumB:

- miR-151 co-amplifies with miR-30-P1a/P2a and miR-4664 in approximately 66% of
  samples.
- miR-4661 co-amplifies with miR-875, miR-599, and miR-4662 in approximately
  69-70% of samples.
- Same-segment CN correlations are very high, often rho 0.85-0.97.
- Typical same-segment span is large, approximately 145 Mb.

Interpretation: chr8 miRNA dosage is a megabase/arm-level subtype event,
especially in LumB. It should be analyzed as a neighborhood amplification block,
not as independent single-miRNA lesions.

### 5.3 chr1 Amplicon - miR-135b, miR-205, miR-29, miR-181

The chr1 block is another high-confidence amplification neighborhood, strongest
in LumB but already visible in LumA.

| Arm | Normal median | LumA median | LumB median | Her2 median | % amp LumB |
|-----|--------------:|------------:|------------:|------------:|-----------:|
| miR-135b-5p | 2 | 4 | 5 | 4 | 77 |
| miR-205-5p | 2 | 4 | 5 | 4 | 75 |

Both have strong subtype differences:

- miR-135b-5p: Kruskal-Wallis p approximately 4.0e-15.
- miR-205-5p: Kruskal-Wallis p approximately 2.0e-13.

Neighborhood read:

- miR-135b and miR-205 co-amplify with miR-29-P2d and miR-181-P1a/P2a in
  approximately 74-75% of LumB samples.
- Same-segment correlations are high, around rho 0.90-0.96.
- Typical same-segment span is large, approximately 106 Mb.

Interpretation: this is a luminal/LumB amplification neighborhood, not a
single-miRNA event. It is biologically interesting because it groups mature
arms with different regulatory reputations: miR-135b, miR-205, miR-29, and
miR-181. Any downstream expression or target-pressure claim needs to respect
this shared CN background.

### 5.4 miR-200c - Basal/EMT-Axis Gain, Not A Hard Amplicon

miR-200c is a different type of case. It is not a high-level amplicon like
miR-21 or chr8. It shows a smaller but coherent Basal-skewed gain.

| PAM50 | n | Median CN [95% CI] | Range | % amp | % diploid |
|-------|--:|--------------------|-------|------:|----------:|
| Normal | 35 | 2 [2, 2] | 1-5 | 14 | 77 |
| LumA | 482 | 2 [2, 2] | 1-8 | 28 | 58 |
| LumB | 250 | 3 [3, 3] | 1-8 | 43 | 34 |
| Basal | 190 | 3 [3, 4] | 1-18 | 46 | 27 |
| Her2 | 90 | 3 [2, 3] | 1-7 | 36 | 43 |

Basal vs LumA median difference is +1 CN copy, FDR-significant. The effect is
smaller than chr8/LumB events but fits the biology: miR-200 family members sit
near epithelial/mesenchymal state control.

Neighborhood read: the strongest same-segment neighbor is the adjacent miR-8
paralog, with approximately 46% Basal co-amplification and rho approximately 1
when both are on the same segment. This is local/co-segment behavior, not a
massive focal amplification signature.

Interpretation: miR-200c should be treated as a Basal-enriched CN gain case,
not a universal high-level amplicon. It is relevant for EMT-state hypotheses,
but the CN effect size is moderate.

### 5.5 let-7a - Multi-Paralog, Mild, Useful Negative Control

let-7a-5p has three paralog loci contributing to the same mature arm. Column
definitions for the decomposition table: §2.6.

| Locus | Chrom | w_norm | Healthy contribution | Mean locus CN | Weighted delta contribution |
|-------|-------|-------:|---------------------:|--------------:|----------------------------:|
| MI0000062 | chr22 | 0.335 | 0.669 | 2.37 | +0.12 |
| MI0000060 | chr9 | 0.333 | 0.666 | 2.63 | +0.21 |
| MI0000061 | chr11 | 0.332 | 0.664 | 2.24 | +0.08 |

Subtype medians stay near diploid:

| PAM50 | Median CN | Range | % amp |
|-------|----------:|-------|------:|
| Normal | 2.0 | 1-4 | 6 |
| LumA | 2.0 | 1-5 | 6 |
| LumB | 2.33 | 1-6 | 6 |
| Basal | 2.33 | 1-5 | 8 |
| Her2 | 2.67 | 1-5 | 16 |

Interpretation: let-7a is a useful control case. It illustrates why the
weighted MIMAT reference is methodologically correct, while also showing that
not every multi-locus MIMAT is a high-amplitude CN event.

---

## 6. Paralog-Specific Warnings

Most high-amplitude cases above are single-locus MIMATs or tightly co-amplified
blocks. But some mature arms have paralog loci with distinct CN histories. These
should not be interpreted using one interval or one unweighted median alone.

**`n_paralog_loci`** = MirGeneDB hairpin loci (`MI*`) that produce that mature
arm (same map as §2.3; cohort summary: `maps/mirna_cnv_mimat_summary.tsv`). Of
**605** MIMATs in the genome map, **36** have `n_paralog_loci > 1` (§3).

Flagged paralog-spread examples (cohort locus-Δ spread ≥ 1.0 among paralogs of
the same arm):

| MIMAT arm | n paralog loci | Cohort locus delta spread | High locus | Low locus |
|-----------|---------------:|--------------------------:|------------|-----------|
| miR-9-5p | 3 | 1.66 | MI0000466 delta +2.30 | MI0000467 delta +0.63 |
| miR-9-3p | 3 | 1.66 | MI0000466 delta +2.30 | MI0000467 delta +0.63 |
| miR-320b | 2 | 1.63 | MI0003839 delta +2.38 | MI0003776 delta +0.75 |
| miR-181a-5p | 2 | 1.62 | MI0000289 delta +2.36 | MI0000269 delta +0.73 |
| miR-181b-5p | 2 | 1.62 | MI0000270 delta +2.36 | MI0000683 delta +0.73 |
| miR-194-5p | 2 | 1.60 | MI0000488 delta +2.35 | MI0000732 delta +0.75 |
| miR-29b-3p | 2 | 1.46 | MI0000107 delta +2.42 | MI0000105 delta +0.96 |

(miR-9-5p and miR-9-3p share the same three precursor loci but are distinct
mature arms with separate `w_norm` splits per locus×MIMAT pair.)

Interpretation: for multi-paralog arms, the key biological unit is the mature
arm, but the genetic perturbation can be carried by only one locus. The
weighted map reports both the parent MIMAT delta and the locus-level
contributions so that an apparent mature-arm dosage shift can be traced back to
its genomic source.

---

## 7. CN -> Expression Concordance And Target Direction

The CN -> miRNA expression question and the miRNA -> target repression question
have opposite sign expectations.

| Link | Expected sign if model holds |
|------|------------------------------|
| miRNA CN vs miRNA expression | Positive |
| miRNA CN/expression vs target RNA | Negative |

### 7.1 TCGA CN → miRNA expression (primary: **partial Spearman**)

The concordance table tests **756** mature arms in TCGA-BRCA primary tumors
(n ≈ 1,032 with both CN and expression per arm). **Primary inferential read:**
**partial Spearman ρ** (`partial_rho`, `partial_q` in
`mirna_cnv_expr_concordance.tsv`) — not bare Pearson/Spearman on raw variables.

**Method.** Per arm, across participants with both measurements:

1. **X** = continuous arm `copy_number` (ASCAT; one value per participant).
2. **Y** = continuous arm expression **log2(RPM+1)**.
3. OLS-residualize **both** X and Y on **CPE** and **thornsson_hrd_score**.
4. **Partial ρ** = Spearman correlation of the residuals
   (`analysis.utils.common.loaders.partial_spearman`).

All participants with data are included (diploid + altered; **no gain-only
filter**). **Marginal** `spearman_rho` / `spearman_q` (no confounder adjustment)
are retained in the TSV as a sensitivity column only.

**Sample scope check.** `concordance_sample_scope_audit()` in
`mirna_locus_cnv.py` recomputes ρ on three participant sets per arm:

| Subset | Rule | Role |
|--------|------|------|
| **Full** (primary) | all participants with CN + expression | default concordance |
| CN > 2 | `copy_number > 2` | sensitivity — amp/gain-enriched |
| gain \| amp | `cn_state ∈ {gain, amp}` | sensitivity — discrete altered call |

Across **756** arms (median **43%** of participants per arm have CN > 2; **730**
arms have ≥20 CN>2 samples for sensitivity). Full-cohort **partial ρ** correlates
with CN>2-only partial ρ at **r ≈ 0.60** (730 arms); **217** arms flip sign
between full and CN>2 — so gain-only analysis is **not** equivalent and is **not**
used for primary calls. Example: miR-151a-3p partial ρ **0.42** (full n=1,032) vs
**0.39** (CN>2 n=681); miR-21-5p **0.19** full vs **0.21** CN>2.

**ρ distribution (boxplot).** Per-arm marginal vs partial ρ across all tested
arms (`figures/mirna_cnv_expr_concordance_rho_boxplot.png`):

![TCGA CN→expression concordance: marginal vs partial ρ](../output/mirna_locus_cnv/figures/mirna_cnv_expr_concordance_rho_boxplot.png)

| Statistic | Marginal Spearman ρ | Partial ρ (CPE + HRD) |
|-----------|--------------------:|----------------------:|
| Median | 0.059 | 0.046 |
| IQR | 0.002 – 0.134 | −0.003 – 0.113 |
| FDR+ positive (q < 0.05, ρ > 0) | 237 / 756 | **196 / 756** |

Red points = arms passing the primary partial FDR filter. Most arms sit near
ρ ≈ 0; the **196** significant positive partial hits are the right-tail outliers
(chr8q, miR-21 neighborhood, miR-200c, etc.).

| Filter | Arms passing |
|--------|-------------:|
| **Partial** partial_q < 0.05 and partial_rho > 0 (**primary**) | **196 / 756** |
| Marginal spearman_q < 0.05 and spearman_rho > 0 (sensitivity) | 237 / 756 |

Top **partial** examples (CPE + HRD adjusted):

| Arm | partial ρ | partial q | n |
|-----|----------:|----------:|--:|
| miR-151a-3p | 0.419 | 9.5e-40 | 1,032 |
| miR-30d-5p | 0.390 | 4.3e-34 | 1,032 |
| miR-30b-5p | 0.385 | 2.3e-33 | 1,031 |
| miR-151a-5p | 0.353 | 7.8e-28 | 1,032 |
| miR-200c-3p | 0.304 | 2.3e-20 | 1,032 |

These support dosage propagation for selected miRNAs. They do not by themselves
prove target repression. Within-PAM50 **partial ρ** for focus arms: Section
8.3.

### 7.2 CCLE cell-line replication (DepMap + CCLE NanoString)

External check: do **PureCN absolute** segment CN at pre-miRNA loci predict
**log2(50 + nSolver)** mature-arm expression in CCLE? Module:
`mirna_hallmark/ccle_mirna_cn_concordance.py`; breast/TAD depth:
`ccle_cn_expr_subtype_depth.py`. Outputs:
`mirna_hallmark/output/ccle_mirna_cn_concordance/` (shared locus-CN cache under
`tables/sample_locus_cn.tsv.gz`).

**Method (CCLE).** Locus CN = max bp overlap of DepMap
`OmicsAbsoluteCNSegmentsProfile` on MirGeneDB hairpins (WES profile when default
WGS is absent). Expression = Broad CCLE miRNA GCT (2018); mirBase arms mapped to
NanoString probe names via alias rules (e.g. miR-151a-3p → `hsa-miR-151-3p`,
30d-5p/3p → precursor `hsa-miR-30d`). **Partial Spearman** on CCLE uses
**OncotreeLineage** dummies (all-lines strata, n ≈ 867–874) or **PAM50 group**
dummies (breast stratum, n ≈ 48 with known subtype from
`pipeline.cell_line_subtype_map` / TADs labels) — not CPE/HRD (unavailable in
local DepMap drop).

**Coverage.** 1,767 cell lines with locus CN; 949 with matched miRNA expression;
96 DepMap breast lines (12 TAD reference lines; **11/12** join CN + expr —
MCF10A lacks both GCT column and DNA profile).

**Arm panels tested** (alias match to CCLE probes):

| Panel | Arms in panel | CCLE expr mapped | Notes |
|-------|-------------:|-----------------:|-------|
| `focus_8q` | **22** (focus amplicon table) | **20** | Same arms as Section 8.3 / neighborhood focus |
| `tcga_sig` | **196** (TCGA partial FDR+) | **156** | Not an intersection with CCLE — TCGA filter only |
| `high_evidence` | 791 (Hallmark edge miRNAs) | 569 | Exploratory superset |

**Results (partial, 2026-06-16 refresh).**

1. **Genome-wide CCLE does not replicate TCGA FDR calls.** Across all lineages
   (lineage-adjusted partial, n ≈ 867), **0 / 20** focus arms and **0 / 156**
   `tcga_sig`-mapped arms reach partial_q < 0.05 with positive rho. Effect sizes
   are much smaller than TCGA (typical partial rho 0.04–0.08 vs 0.3–0.4 in
   tumors).

2. **Breast-only CCLE (PAM50-adjusted, n ≈ 48)** raises partial ρ for several focus
   arms but still **0 / 20** FDR passes. Top breast **partial** ρ among focus
   arms: let-7a-5p **0.28**, miR-454 **0.22**, miR-205 **0.22**, miR-30b **0.20**,
   miR-301a **0.19**, miR-30d **0.18**, miR-151a-3p **0.14**. This is **not**
   restricted to chr8q: chr1 (miR-205), chr17 (miR-21 weak), and let-7a show
   comparable or stronger breast-lineage partial coupling than some 8q arms.
   Chr8q amplicon miRNAs (151a, 30b, 30d, 301a) remain the **best TCGA sign
   match** (same positive partial direction in both cohorts) but with ~5–10×
   attenuated magnitude in CCLE.

3. **TAD reference panel — the only CCLE stratum with visible marginal ρ (see below).**
   Do **not** read this as “CCLE replicates TCGA globally.” It is an **11-line
   breast model subset** (Kim/TAD biosamples mapped via
   `pipeline.cell_line_subtype_map`) with **n too small for partial adjustment**
   (PAM50 partial is undefined at n = 11). All TAD **marginal** ρ values below are
   **marginal Spearman** with **one point per line per precursor probe** (12
   focus precursors from the 20 mapped focus arms).

##### TAD panel composition (11 / 12 reference lines)

| PAM50 (TADs map) | n | Lines | Role in CCLE test |
|------------------|--:|-------|-------------------|
| LumA | **4** | MCF7, T47D, ZR751, ZR7530 | Only stratum with ≥4 lines → within-subtype marginal ρ possible |
| Basal | **5** | BT549, HCC70, MDA-MB-231, MDA-MB-468, HCC1954 | Within-subtype marginal ρ possible |
| LumB | **1** | BT474 | Descriptive CN/expr only (no correlation) |
| HER2 | **1** | SKBR3 | Descriptive only |
| Normal_like | **0** | MCF10A excluded (no GCT, no CN) | — |

##### Three ways to read the same 11 lines (global situation)

**A. Pooled across all 11 lines (what looked like “signal”).**  
Each precursor: correlate CN vs expression across the 11 lines (mixes LumA, Basal,
LumB, HER2).

| Precursor | marginal ρ (11 lines) | p | Positive? |
|-----------|---------------------:|---:|:---------:|
| miR-30d | **0.61** | 0.048 | yes (only borderline p) |
| miR-205 | 0.40 | 0.22 | yes |
| miR-30b | 0.37 | 0.27 | yes |
| miR-151-3p | 0.32 | 0.33 | yes |
| miR-21 | 0.31 | 0.36 | yes |
| miR-454 | 0.26 | 0.44 | yes |
| miR-200c | 0.01 | 0.97 | flat |
| let-7a | −0.27 | 0.43 | no |

**9 / 12** precursors have marginal ρ > 0 pooled; **5 / 12** have marginal ρ > 0.3. This is the
**only** CCLE stratum approaching TCGA effect sizes — but it **confounds subtype
with dosage**: LumA lines (especially MCF7, ZR7530) carry higher chr8q CN on
several precursors while Basal lines cluster lower, and expression shifts
partly follow that split (see miR-30d medians: Basal CN 4.0 / expr 7.9 vs LumA
CN 5.0 / expr 8.0 vs LumB CN 6.0 / expr 8.2). Correlating across all 11 lines
partly measures **“Lum-like vs Basal-like line identity”**, not the same
within-tumor partial relationship as TCGA Section 8.3.

**B. Within PAM50 group only (subtype-stratified line-level marginal ρ).**  
Here each stratum uses only lines of one PAM50 label — the fair CCLE analogue of
TCGA Section 8.3 (but n = 4 or 5, **marginal not partial**, not FDR-tested).

| Precursor | LumA marginal ρ (n=4) | Basal marginal ρ (n=5) | Same sign in both? |
|-----------|----------------------:|-------------------------:|:------------------:|
| miR-30d | **0.74** | 0.21 | yes, weak Basal |
| miR-151-3p | **−0.63** | **0.56** | **no — opposite** |
| miR-151-5p | −0.32 | 0.56 | **no** |
| miR-200c | −0.20 | 0.45 | no |
| miR-205 | 0.80 | 0.00 | LumA only |
| miR-454 | 1.00 | 0.53 | yes |
| miR-301a | 0.80 | 0.53 | yes |
| miR-21 | 0.40 | −0.53 | no |

**Within-subtype read:** dosage-like CN→expr coupling is **not uniform**. It is
**split by PAM50 bucket**:

- **LumA (4 lines):** chr8q precursors **miR-30d / 30b / 301a / 454** show the
  strongest *within-LumA* positive marginal ρ (miR-30d marginal ρ = 0.74), to be
  compared with TCGA **partial** ρ ≈ 0.36–0.39 for 30d/30b (Section 8.3).
  **miR-151a is negative within LumA** (marginal ρ ≈ −0.63) despite strong TCGA
  LumA/LumB/Her2 **partial** — the pooled 11-line 151 signal is driven by
  **Basal lines**, not LumA.
- **Basal (5 lines):** **miR-151a** shows the best within-stratum marginal ρ
  (≈ 0.56), closer to TCGA Basal **partial** (151a partial ρ ≈ 0.36).
  **miR-200c** is modestly positive within Basal (marginal ρ ≈ 0.45), echoing
  TCGA Basal-strong 200c (partial ρ ≈ 0.45). **miR-30d** is weak within Basal
  (marginal ρ ≈ 0.21) even though TCGA Basal **partial** for 30d is 0.46 —
  LumA carries the 30d line-level story in CCLE, not Basal.
- **LumB (BT474) and HER2 (SKBR3):** single lines — inspect per-line CN/expr
  only. SKBR3 has **low** miR-151 CN (2.0) and low expr vs amp LumA lines;
  BT474 has high miR-21 / miR-301a CN (17q context) with moderate expr.

**C. Subtype centroids (4 PAM50 medians — not independent samples).**  
If one collapses to median CN/expr per PAM50 group (4 points), miR-30d shows
marginal ρ ≈ 0.95 — but that is **between-subtype ordering** (Basal < LumA < LumB for CN),
not replication of within-subtype TCGA **partial** ρ. Do not treat this as confirmatory.

##### Per-line snapshot (why pooled marginal ρ can mislead)

**miR-30d** (TCGA partial ρ ≈ 0.39, broad across subtypes): CN ranges 3–12
(ZR7530 amp) and expr 6.1–10.3. High-CN LumA lines (MCF7, ZR7530) have high
expr; Basal lines sit mid-range — pooled marginal ρ ≈ 0.61 is partly **ZR7530/MCF7
leverage**.

**miR-151-3p** (TCGA partial ρ ≈ 0.42): SKBR3 HER2 CN = **2.0** (low); LumA
MCF7/ZR7530 CN = **6–12**; Basal lines CN = **3–6**. Pooled marginal ρ > 0 comes from
**Basal + low-HER2** vs **high-LumA** spread; **within LumA marginal ρ is negative**
because ZR7530 has high CN but **low** 151 expr (6.9) vs MCF7/ZR751.

##### High-evidence panel on TAD 11 lines (569 mapped arms)

The focus-8q precursor tables above use **12 precursors** from 20 mapped focus
arms. The complementary test uses the full **Hallmark high-evidence** arm panel
(**791** arms → **569** CCLE-mapped) restricted to the same **11 TAD breast
lines** (MCF10A excluded). Two correlation designs:

1. **Panel-wide (per arm):** marginal Spearman(CN, expr) for each arm across the
   11 lines (plus within LumA n=4 and Basal n=5). **387 / 569** arms have
   complete CN+expr on ≥3 lines for the pooled stratum.
2. **Per line (across arms):** within each TAD line, Spearman(CN, expr) across
   all mapped arms with ≥30 non-missing pairs (~384–387 arms per line).

**Panel-wide (pooled 11 lines):** median marginal ρ = **0.06**; **89** arms have marginal ρ > 0.3
(e.g. miR-197-3p marginal ρ = 0.67, miR-106b marginal ρ = 0.66, miR-30d marginal ρ = 0.61, miR-205 marginal ρ =
0.40). **0 / 387** survive BH-FDR at q < 0.05 — power is too low to declare
genome-wide replication, but the **chr8q / 17q cluster** from TCGA **partial**
concordance still ranks among the top TAD pooled hits (30d TCGA partial ρ =
0.39 → TAD pooled marginal ρ = 0.61, p = 0.048 uncorrected).

**Within PAM50 on high-evidence arms:** LumA median marginal ρ = **0.00** (139 arms marginal ρ >
0.3 — driven by a few chr8q outliers); Basal median marginal ρ = **0.05** (119 arms marginal ρ >
0.3). Subtype-split sign flips persist at scale (not only focus precursors).

**Per-line across all high-evidence arms:**

| Line | PAM50 | marginal ρ (across ~387 arms) | q |
|------|-------|------------------------------:|--:|
| SKBR3 | HER2 | **+0.15** | **0.012** |
| BT474 | LumB | +0.10 | 0.11 |
| HCC70 | Basal | +0.08 | 0.18 |
| BT549 | Basal | +0.08 | 0.18 |
| ZR751 | LumA | +0.07 | 0.23 |
| HCC1954 | Basal | +0.05 | 0.38 |
| MDA-MB-468 | Basal | +0.04 | 0.42 |
| T47D | LumA | −0.05 | 0.38 |
| MCF7 | LumA | −0.10 | 0.11 |
| **ZR7530** | LumA | **−0.22** | **1.1×10⁻⁴** |
| **MDA-MB-231** | Basal | **−0.24** | **1.6×10⁻⁵** |

Only **SKBR3** shows FDR-positive *global* CN→expr concordance across the
high-evidence arm set; **ZR7530** and **MDA-MB-231** are strongly
**anti-concordant** (high CN does not track high expression genome-wide on those
lines). This matches the per-precursor leverage story: ZR7530 amps chr8q CN but
not uniformly chr8q expression.

**TCGA partial vs TAD pooled marginal (130 TCGA partial-q<0.05 arms in join):** **80 / 130**
have the same positive sign in TAD pooled marginal ρ; best TAD matches include miR-197,
106b cluster, **miR-30d** (TCGA partial ρ 0.39 → TAD marginal ρ 0.61), miR-93/25/324.

Outputs: `breast_tad_lines/high_evidence/tad_panel_arm_concordance.tsv`,
`per_line_across_arms_concordance.tsv`, `per_line_arm_cn_expr.tsv.gz`,
`tcga_partial_vs_tad_panel_compare.tsv`.

##### How this compares to TCGA Section 8.3 (same focus arms)

| Arm / precursor | TCGA partial ρ (strongest stratum) | CCLE TAD pooled marginal ρ (11) | CCLE within LumA marginal ρ | CCLE within Basal marginal ρ |
|-----------------|-----------------------------------:|--------------------------------:|----------------------------:|-----------------------------:|
| miR-151a | LumB/Her2 partial ρ ≈ 0.53–0.56 | marginal ρ ≈ 0.32 | **marginal ρ ≈ −0.63** | **marginal ρ ≈ 0.56** |
| miR-30d | broad partial ρ ≈ 0.36–0.46 | **marginal ρ ≈ 0.61** | **marginal ρ ≈ 0.74** | marginal ρ ≈ 0.21 |
| miR-200c | **Basal partial ρ ≈ 0.45** | marginal ρ ≈ 0.01 | marginal ρ ≈ −0.20 | marginal ρ ≈ 0.45 |
| miR-205 | Her2 partial ρ ≈ 0.52 | marginal ρ ≈ 0.40 | marginal ρ ≈ 0.80 | marginal ρ ≈ 0.00 |

**Bottom line for TAD-only CCLE:** the panel is **informative but not
confirmatory**. Visible marginal ρ appears only when **(i)** restricting to TAD breast
lines and **(ii)** pooling across subtypes, or when inspecting **within-LumA**
chr8q precursors / **within-Basal** 151a and 200c with **n = 4–5** and wide
p-values. It does **not** reproduce TCGA’s genome-wide **partial** FDR block (196
arms), and **subtype-stratified CCLE contradicts a simple global dosage rule**
(e.g. 151a sign flips LumA vs Basal). Treat TAD CCLE as **line-biology context**
for the same biosamples used in TAD/HiChIP — not as a second cohort validating
tumor **partial** concordance.

Outputs: `breast_tad_lines/tad_line_precursor_cn_expr.tsv` (per-line CN/expr),
`breast_tad_lines/breast_line_level_pam50_correlation.tsv` (within-PAM50 line marginal ρ),
`breast_tad_lines/tad_line_coverage_audit.tsv`,
`breast_tad_lines/high_evidence/` (full high-evidence panel-wide + per-line).

**CCLE summary (all strata).**

| CCLE stratum | n | Partial adjustable? | FDR replication? | Typical partial ρ vs TCGA partial ρ |
|--------------|--:|:-------------------:|:----------------:|:------------------------------------|
| All lineages | ≈867 | Yes (lineage dummies) | **No** (0 / 20 focus) | 0.04–0.08 vs 0.3–0.4 |
| Breast (DepMap) | ≈48 | Yes (PAM50 dummies) | **No** | 0.14–0.28 vs 0.3–0.4 |
| **TAD panel pooled** | **11** | **No** | **No** (marginal only) | marginal 0.3–0.6 vs TCGA partial 0.3–0.4 (some precursors) |
| TAD within LumA / Basal | 4 / 5 | No | No | marginal, subtype-split (tables above) |
| TAD per-line (HE across arms) | 1 line × ~387 arms | No | **1 / 11** (SKBR3 q=0.012) | ZR7530 / MDA-MB-231 **anti**-concordant (marginal) |

Also: `focus_8q/ccle_mirna_cn_expr_concordance.tsv`,
`focus_8q/tcga_ccle_concordance_compare.tsv`.

The target-resolved join over top CN-concordant arms (selected on **marginal**
`spearman_q` in `concordance_target_join.py`; see note below) contains 167
miRNA-target rows across 21 arms. Repressor-consistent negative links exist, but
they are mostly modest. **Target-join ρ values are marginal Spearman** (not
partial-adjusted):

| Link | Stronger large-n examples |
|------|---------------------------|
| miRNA CN vs target RNA | miR-92b-3p -> SMAD3 marginal ρ -0.19; miR-21-3p -> EGFR marginal ρ -0.18; miR-4661-5p -> IL6ST marginal ρ -0.17 |
| miRNA expression vs target RNA | miR-17-5p -> ERAP1/JAK1/TGFB1/SMAD3/STAT3 marginal ρ -0.20 to -0.28; miR-200c-3p -> EGFR marginal ρ -0.26 |

Interpretation:

- CN gain clearly raises dosage for some miRNA arms.
- Target repression is much harder to see in bulk RNA.
- Positive miRNA CN-target RNA associations are not "expected direction"; they
  can reflect shared CN blocks, purity, cell-state coupling, or weak target
  regulation in bulk.
- Stronger target claims should use the **canonical pressure framework** (§1b,
  Section 8): **logrpm + evidence_mass** on miRTarBase edges, AGO-gated, with subtype and
  target-CN controls. Hybrid **M7** (`softmax_z` + `combined_mass`) extends the edge
  universe when TS agreement is needed (exploratory IRF1 hub partial).

---

## 8. Combined Target Repression And Subtype Stratification

Pressure for this section: **M0 miRTarBase edges**, **`softmax_z_logrpm` +
`evidence_mass` + `sum`**, via `pressure_build.compute_gene_pressure` (see §1b).
AGO gate applied by default. TargetScan hybrid modes (M7/M8/M11) use
`PRESSURE_HYBRID_*` weighting — reported under `hybrid_pressure/` for
Hallmark-program coupling and exploratory hub partials only.

### 8.1 Cohort-wide combined pressure anti-correlation

For each Hallmark-universe target gene, all targeting MIMAT arms fold into one
combined pressure vector per sample, then Spearman-tested against target RNA.
Confounder-adjusted partials: CPE + HRD + target gene CN.

Cohort results (AGO-gated view, n = 1,077 participants, 1,349 genes tested).
Primary inferential read: **partial Spearman** (+ target CN).

| Metric | Count |
|--------|------:|
| Negative FDR q < 0.05 (partial rho + target CN) | **373 / 1,349** |
| Negative FDR q < 0.05 (partial rho, CPE + HRD only) | 385 / 1,349 |
| Negative FDR q < 0.05 (raw Spearman; sensitivity) | **480 / 1,349** |

Top repressor-consistent hits (gated, partial rho + target CN):

| Gene | partial rho | q |
|------|------------:|--:|
| COL5A2 | -0.42 | 7.3e-40 |
| ZFPM2 | -0.41 | 4.2e-39 |
| SPARC | -0.39 | 5.9e-35 |
| GSN | -0.38 | 3.5e-33 |
| COL1A1 | -0.36 | 1.0e-29 |
| LOX | -0.36 | 1.1e-29 |
| ADAM12 | -0.35 | 1.2e-27 |
| ESR1 | -0.34 | 7.8e-27 |
| KLF2 | -0.34 | 4.9e-26 |
| TGFB3 | -0.34 | 6.8e-26 |

These are bulk-RNA repressor signatures after confounder adjustment. Counts rose
vs the prior softmax_z+degree model (370 partial hits) because **logrpm +
evidence_mass** strengthens repressor-consistent signal; top ECM/stromal genes
are preserved.

### 8.2 PAM50-stratified repression: where subtype matters

Within-subtype tests reuse the same combined pressure but restrict samples to
each PAM50 group (minimum n = 20). All gene-level tables below use **partial
rho (+ target CN)**; FDR is computed separately per subtype. LumA is the largest
stratum (n approximately 503), so it contributes the most FDR passes; compare
effect sizes across subtypes, not raw hit counts alone.

Negative FDR q < 0.05 by PAM50 (gated, partial rho + target CN):

| PAM50 | n (approx.) | genes with sig negative partial rho |
|-------|------------:|------------------------------------:|
| LumA | 503 | **228** |
| LumB | 258 | **139** |
| Basal | 188 | **94** |
| Her2 | 92 | **28** |
| Normal | 36 | **1** |

Compared with raw Spearman within subtype (311 / 168 / 107 / 45 / 0 in the same
strata), partial adjustment removes confounded hits while preserving core
stromal/ECM repressor genes.

Pattern summary (partial + target CN):

- **187 genes** are significant in exactly one subtype; **zero genes** pass FDR
  in all five strata on partials.
- **32 genes** are Basal-only (partial sig in Basal but not LumA); **33 genes**
  are LumB-only (partial sig in LumB but not LumA).

**Basal-enriched repressor signal** (stronger partial anti-correlation in Basal
than cohort-wide, FDR sig within Basal):

| Gene | Cohort partial | Basal partial | Notes |
|------|---------------:|--------------:|-------|
| DGUOK | +0.01 | -0.24* | Basal-only partial FDR |
| CD4 | +0.01 | -0.22* | immune; Basal-only partial FDR |
| TNFRSF12A | -0.07 | -0.32* | stronger within Basal |
| FCGR2B | -0.11 | -0.21* | immune receptor; Basal-enriched |
| SSRP1 | +0.05 | -0.18 | directional in Basal only (q > 0.05) |
| CDKN2A | -0.11 | -0.12 | LumA partial stronger (-0.23*) |

**LumB-enriched repressor signal** (partial):

| Gene | Cohort partial | LumB partial | Notes |
|------|---------------:|-------------:|-------|
| VIM | -0.03 | -0.22* | cohort null; LumB EMT repressor read |
| RAB17 | -0.13 | -0.34* | stronger within LumB amp context |
| NOS2 | -0.13 | -0.32* | LumB >> cohort |
| VCAN | -0.07 | -0.20* | stromal ECM |
| ERN1 | +0.03 | -0.17 | UPR axis; LumB directional only (q > 0.05) |

**Route genes and panel anchors** — partial rho (+ target CN) within PAM50;
* = FDR q < 0.05:

| Gene | LumA | LumB | Basal | Her2 | Interpretation |
|------|-----:|-----:|------:|-----:|----------------|
| LOX | -0.32* | -0.39* | -0.37* | -0.43* | robust pan-subtype |
| GSN | -0.37* | -0.37* | -0.39* | -0.36* | same |
| COL1A1 | -0.30* | -0.31* | -0.32* | -0.25 | similar LumA/LumB/Basal; Her2 weakest |
| ZFPM2 | -0.38* | -0.38* | -0.30* | -0.48* | Her2 partial peak |
| VIM | -0.17* | -0.22* | -0.16 | -0.18 | LumB strongest partial |
| IRF1 | -0.08 | -0.16 | +0.05 | -0.06 | no subtype partial sig |
| TGFBR2 | -0.03 | -0.07 | -0.18 | +0.01 | no subtype partial sig |
| PTEN | +0.07 | +0.07 | -0.13 | +0.14 | no within-subtype partial sig |
| PGR | -0.21* | +0.04 | -0.06 | -0.12 | LumA-specific partial |
| EGFR | -0.07 | +0.14 | -0.11 | +0.06 | cohort positive (+0.11*); not repressor |
| MYC | -0.13* | -0.03 | +0.25 | -0.31 | LumA negative; Basal positive (not FDR) |
| CDKN2A | -0.23* | -0.16 | -0.12 | -0.03 | LumA partial sig |
| CD4 | +0.06 | +0.02 | -0.22* | +0.05 | Basal-only partial FDR |
| SSRP1 | +0.06 | +0.07 | -0.18 | +0.03 | Basal directional only (q > 0.05) |

Takeaways for interpretation:

1. **Stromal/ECM targets (LOX, GSN, ZFPM2, KLF2)** survive within-subtype
   partial adjustment across multiple strata — the most credible bulk repressor
   signatures in this layer under the canonical pressure model.
2. **Basal partial hits are distinct from cohort mixing.** CD4 and DGUOK pass
   FDR only within Basal (cohort partial null); TNFRSF12A and FCGR2B are
   Basal-enriched relative to cohort-wide reads.
3. **LumB amp context + repression:** VIM, RAB17, NOS2, and VCAN show LumB
   partial FDR passes stronger than the cohort partial column.
4. **PTEN cohort partial (-0.09) is weak and does not localize** — no
   within-subtype partial FDR pass — consistent with cross-stratum composition
   rather than LumA-intrinsic repression.
5. **EGFR (cohort +0.11*) and MYC (Basal +0.25)** stay positive or weak after
   partial — not miRNA-repressor stories under M0 logrpm+evidence_mass pressure.

### 8.3 PAM50-stratified CN to expression (focus arms)

CN to miRNA expression concordance also varies by subtype for the focus amplicon
arms. Values below are **partial Spearman rho** (CN vs mature-arm RPM), adjusting
for CPE and HRD within each PAM50 stratum:

| Arm | Basal partial ρ | LumA partial ρ | LumB partial ρ | Her2 partial ρ | Notable pattern |
|-----|----------------:|---------------:|---------------:|---------------:|-----------------|
| miR-151a-3p | 0.36 | 0.32 | 0.56 | 0.63 | strongest in Her2/LumB amp strata |
| miR-200c-3p | **0.45** | 0.16 | 0.26 | 0.33 | **Basal-strongest** EMT-axis dosage |
| miR-30d-5p | 0.46 | 0.36 | 0.41 | 0.48 | broad concordance |
| miR-17-5p | 0.33 | 0.25 | 0.31 | 0.44 | chr13 neighborhood |
| miR-21-5p | 0.16 | 0.12 | **0.30** | 0.10 | LumB amp context helps |
| miR-205-5p | 0.11 | 0.12 | 0.13 | **0.52** | Her2-local chr1 block |
| let-7a-5p | 0.15 | 0.09 | 0.18 | 0.09 | weak everywhere |

miR-200c is the standout subtype-specific dosage case: partial CN to expression
coupling is approximately 3x stronger in Basal than in LumA, matching the
Basal-enriched miR-200c co-gain pattern in Section 9. miR-151 and miR-30d show
the expected positive CN-expression link across strata, with the highest partial
rho in the most amp-burdened subtypes (LumB, Her2).

### 8.4 Method note

Combined pressure uses **`pressure_build.compute_gene_pressure`**: **`softmax_z_logrpm`**
× **`evidence_mass`** on **M0 miRTarBase** edges (evidence ≥ 2). This replaces
the pre-2026-06 parent `compute_pressure` cohort-z-only sum and the interim
softmax_z+degree default. **All Section 8 gene-level tables report partial Spearman:**
target-repression tests adjust for CPE + HRD + target gene CN; CN→expression
tests in Section 8.3 adjust for CPE + HRD within each PAM50 stratum (PAM50
excluded because stratum-fixed). FDR is computed separately per (view, PAM50)
stratum for within-subtype repression. Normal-like (n ≈ 36) has very low partial
power; treat its single partial hit as descriptive only. Raw Spearman counts
(480 cohort / 631 within-subtype gene×PAM50 pairs) are in the Section 8.1 table. For TS
hybrid edge modes and Hallmark-level coupling, see §1b.D and `hybrid_pressure/`.

---

## 9. Neighborhood Amplification And Deletion Context

Neighborhood analysis was computed two ways:

1. **Same ASCAT3 segment:** participant-specific loci sharing the anchor's
   segment.
2. **+/-2 Mb window:** static genomic proximity around the anchor locus.

Same-segment analysis is the more informative read for TCGA ASCAT3 here. Many
events span tens to hundreds of Mb, so the correct interpretation is often
chromosome-arm co-amplification.

Key subtype-stratified co-amplification patterns:

| Anchor | Peak subtype | Neighborhood behavior |
|--------|--------------|-----------------------|
| miR-21 | LumB / Her2 | 17q same-segment co-amp with miR-130/142/4524/3615, approximately 51-57% LumB co-amp |
| miR-151 | LumB | chr8 same-segment co-amp with miR-30, miR-4662/4664, miR-599, approximately 65-66% LumB co-amp |
| miR-4661 | LumB | chr8 co-amp with miR-875, miR-599, miR-4662, approximately 69-70% LumB co-amp |
| miR-135b | LumB | chr1 co-amp with miR-29, miR-181, miR-205, approximately 74-75% LumB co-amp |
| miR-205 | LumB | same chr1 block as miR-135b, miR-29, miR-181 |
| miR-200c | Basal / LumB | local miR-8 family segment co-gain, approximately 46% Basal co-amp |
| let-7a | Her2 weak | modest multi-locus gain, not a focal amplicon |

Deletion context: for the focus cases, coupled "anchor amplified, neighbor
deleted" patterns are near 0%. This refresh does not support a gain/loss
neighborhood grammar for miRNA loci. The dominant pattern is co-gain and
co-amplification.

---

## 10. Main Biological Conclusions

1. **The miRNA CNV layer in BRCA is amplification-dominated.** The global
   mature-arm and locus landscapes are shifted above diploid, with the strongest
   burden in LumB and Her2.

2. **Most strong miRNA CN effects are neighborhood effects.** miR-21, miR-151,
   miR-4661, miR-30d/b, miR-135b, and miR-205 are embedded in large same-segment
   events. Treating them as isolated hairpin lesions would overstate locus-level
   specificity.

3. **Subtype specificity is real but case-dependent.** LumB/Her2 carry the hard
   amp phenotype. Basal shows a more selective miR-200c/EMT-axis gain. LumA is
   less CNV-burdened at the median but still has amplicon-positive tails.

4. **Weighted MIMAT handling is essential but not always transformative.** It is
   crucial for paralog arms with discordant locus CN. For single-locus or
   co-amplified loci, weighted and median summaries agree.

5. **CN -> expression is easier to establish than expression -> target
   repression.** The strongest CN -> expression arms are credible dosage
   propagation events (Section 8.3 partial CN→expr). Target repression requires
   anti-correlation and is modest in bulk RNA after partial adjustment: **373 / 1,349**
   genes cohort-wide (partial + target CN, logrpm+evidence_mass M0 pressure), with
   zero passing FDR in all five PAM50 strata; Basal and LumB carry distinct
   partial signatures (CD4/DGUOK vs VIM/RAB17). Hallmark-level coupling under
   the spine weighting reaches **42/50** Basal neg-sig programs (M0, prolif-adj).

---

## 11. Caveats And Next Analysis Priorities

1. **Normal-like is not healthy normal.** Use it as a low-alteration empirical
   anchor, not a biological normal breast reference.
2. **ASCAT3 segments are often broad.** Same-segment co-amplification can reflect
   arm-level events, not focal miRNA-specific selection.
3. **Bulk RNA hides repression.** Target suppression should be tested with
   AGO-gated pressure, subtype stratification, and target CN controls. Section 8
   reports partial Spearman throughout (CPE + HRD + target CN for repression;
   CPE + HRD for CN→expression).
4. **Some labels are missing for rare MIMATs.** Unlabeled MIMATs are retained in
   tables but should be resolved before figure presentation.
5. **Paralog weights are cohort-reference estimates.** They are deliberately not
   sample-specific. That avoids circularity, but it cannot capture tumor-specific
   locus usage.

Recommended next steps:

- Add a compact figure panel for the focus arms: PAM50 box/violin distributions
  plus same-segment co-amp bars.
- Convert same-segment neighborhoods into named genomic blocks (chr8, chr1,
  chr17, chr12/miR-200 family) with coordinates and carrier counts.
- Manually review unlabeled top MIMATs from the chr8 block before presentation.

---

## 12. CN → Regulatory-Network Attribution (End-To-End)

Sections 7–8 measure the chain as two *separate* halves — CN→miRNA-expression
concordance (§7/§8.3) and pressure→target repression (§8). This section closes the
loop: **how much of the bulk repressive network is actually driven by the genetic
(copy-number) layer, in which context, and at which genomic granularity.** Module:
`cn_dosage_attribution.py`; outputs `output/cn_dosage_attribution/`. Pressure is the
spine (`softmax_z_logrpm` + `evidence_mass`, M0 miRTarBase, AGO-gated); all couplings
are partial Spearman over a **cumulative covariate stack**: baseline (CPE + HRD +
target CN) → `*_aneu` (+ broad genome-wide aneuploidy, Taylor `thornsson_aneuploidy_score`)
→ `*_full` (+ the E2F/G2M **proliferation** metagene, as `robustness_checks`) — to
separate a locus-specific CN effect from the genome-wide CIN tide *and* the
proliferation/metabolic state.

### 12.1 Method

Each regulator arm's log2(RPM+1) is split per-arm by OLS on its own copy number into a
**CN-driven** component (`cn_hat = a + b·CN`) and a **CN-residual** (CN-independent)
component. Three pressure tracks are built on the *same* edges and tested against
target RNA: `total` (real expression), `cn` (CN-driven only), `resid` (residual only).
The share of the total repressive coupling reproduced by the `cn` track =
`cn_attribution = ρ_cn / ρ_total` (evaluated only where ρ_total < −0.05).

### 12.2 How much (cohort, n = 1,408 genes tested, AGO-gated)

Every CN coupling is reported under a **cumulative covariate stack**: baseline
(CPE+HRD+target-CN) → `+aneu` (broad genome-wide aneuploidy, `thornsson_aneuploidy_score`)
→ `+full` (also the E2F/G2M **proliferation** metagene, as `robustness_checks`).

| Track | Neg-FDR (baseline) | + aneuploidy | + aneuploidy + proliferation (**full**) |
|-------|-------------------:|-------------:|----------------------------------------:|
| `total` (real expression) | **383** | — | — |
| `cn` (CN-driven only) | **744** | **665** | **675** |
| `resid` (CN-independent only) | 353 | — | — |

| Median `cn_attribution` among total-neg genes | baseline | + aneu | + **full** |
|-----------------------------------------------|---------:|-------:|-----------:|
| (share of total repressive ρ reproduced by `cn`) | **0.70** | 0.54 | **0.51** |

- The `total` count (383) reproduces the canonical `target_combined_anticorr`
  headline (373/1,349) — sanity check that the spine is unchanged.
- **The CN-driven dosage component reproduces ~half-to-two-thirds of the genuine
  repressive network and HOLDS under the full control stack**: median attribution
  0.70 → 0.54 (aneuploidy) → **0.51 (also proliferation)**; **57.7 %** of total-neg
  genes are CN-driven. So ~**half the bulk repressive network is attributable to the
  copy-number layer net of purity, HRD, target-CN, aneuploidy, and proliferation.**
- **Proliferation is a *suppressor*, not the confound.** Adding it *increases* the
  CN-driven count (cohort 665 → **675**; Basal **127 → 262**, Her2 **17 → 72**) rather
  than collapsing it — because proliferation is positively tied to both pressure and
  target expression, it was *masking* the negative coupling. This is the same
  suppressor structure `robustness_checks` found for the basal hub (**MH-17**), now
  reproduced for the CN layer — and it is the decisive answer to "isn't this just a
  proliferation/metabolic shadow?": no, proliferation *hides* CN's effect.

> **Note on the raw 744.** The `cn` track surfaces ~2× more neg-FDR genes than `total`
> partly because `cn_hat` is a **denoised dosage instrument** (it strips CN-independent
> expression fluctuation, raising power), and partly because some surplus genes have a
> *flat* `total` pressure and are metabolic/mitochondrial (DLD, FH, FXN, MCL1, SMAD2).
> Quote the **attribution among genuine total-neg genes (0.70 → 0.51 full)**, not the
> raw `cn` count. The proliferation control (which only strengthens it) settles the
> earlier metabolic-shadow caveat.

### 12.3 In which context (CN's grip on the network is luminal / immune-hot)

Within-stratum CN attribution (partial, AGO-gated; `cn_attribution_summary.tsv`):

| Stratum | n_neg_total | n_neg_cn | +aneu | +**full** | median CN-attribution (raw / +aneu / +full) |
|---------|------------:|---------:|------:|----------:|:-------------------------------------------:|
| **LumA** | 221 | 602 | 549 | 517 | **0.56 / 0.53 / 0.46** |
| LumB | 136 | 269 | 248 | 277 | 0.24 / 0.22 / 0.23 |
| Basal | 77 | 193 | 127 | **262** | 0.29 / 0.26 / 0.29 |
| Her2 | 28 | 40 | 17 | **72** | 0.37 / 0.39 / 0.38 |
| Normal | 1 | 10 | 19 | 3 | ~null |
| **Immune C2** (IFN-γ dominant) | 204 | 546 | 496 | 512 | **0.57 / 0.48 / 0.46** |
| Immune C1 (wound healing) | 175 | 480 | 353 | 455 | 0.56 / 0.42 / 0.49 |
| Immune C3 | 79 | 317 | 255 | 276 | 0.36 / 0.33 / 0.31 |
| Immune C4 (lymph-depleted) | 34 | 85 | 59 | 65 | 0.05 / 0.05 / 0.05 |
| TNBC subtypes (BL1/BL2/LAR/M) | ≤6 | ≤20 | ≤0 | ≤20 | underpowered |

- **CN drives the miRNA network most in LumA and the immune-hot C1/C2 strata** —
  and there it survives the full stack (LumA attribution 0.56→0.46; C2 0.57→0.46).
- **Basal and Her2 are where proliferation was *masking* CN**: their CN-driven counts
  **rise** under proliferation adjustment (Basal 127→**262**, Her2 17→**72**) — the
  suppressor structure of §12.2 / MH-17. So Basal *does* carry a real CN-driven
  component (previously hidden), though its per-gene attribution stays **modest (0.29)**:
  Basal repression is still comparatively **expression/state-driven** (selective
  miR-200/EMT, §5.4/§8.2) rather than amplicon-dosage-driven.
- **LumB** has the heaviest CN burden (§4) yet a **modest attribution (0.23)**: CN
  surfaces *many* neg genes but reproduces only a quarter of each — CN moves a broad
  set weakly rather than a focused set strongly.

### 12.4 At which granularity — does combining loci beat one locus?

Per target gene, a single per-sample scalar = the evidence-weighted summed CN dosage
(Δ vs diploid) of all its regulators, built at four resolutions and tested for
repressor-consistent (negative) coupling (`cn_load_granularity_summary.tsv`,
cohort n = 1,335):

| Granularity | Neg-FDR genes | median ρ (neg) | median ρ (all) |
|-------------|--------------:|---------------:|---------------:|
| `single_best` (most-amplified single regulator) | 842 | −0.141 | −0.098 |
| `locus` (combined, all regulators) | 900 | −0.156 | −0.108 |
| `cluster` (dense-cluster co-dosage) | 897 | −0.156 | −0.108 |
| **`chrom`** (broad chromosome-scale aneuploidy) | **941** | **−0.164** | −0.120 |
| `locus` + aneuploidy | 811 | −0.144 | −0.094 |
| `chrom` + aneuploidy | 834 | −0.154 | −0.105 |
| `locus` + aneuploidy + proliferation (**full**) | 852 | −0.160 | −0.104 |
| `chrom` + aneuploidy + proliferation (**full**) | 886 | −0.170 | −0.115 |

- **Combining loci beats any single locus**: the combined multi-locus load couples
  more (900 vs 842 genes; median ρ −0.156 vs −0.141) than the single most-amplified
  regulator. The regulatory CN effect is a **set property**, not a single-hairpin one.
- **Broad (chromosome-scale) load is the strongest unit** (941 / −0.164) > combined
  locus ≈ cluster — matching the landscape finding that miRNA CN is arm/chromosome
  scale, not focal (§3, §9).
- **Both survive the full control stack and even sharpen** (`locus_full` 852 / −0.160;
  `chrom_full` 886 / −0.170 — stronger median ρ than baseline, the same proliferation
  *suppressor* effect). The broad edge over locus persists (886 vs 852) but
  locus-combined load is robust throughout — a genuine locus-specific multi-locus
  dosage component beneath the broad aneuploidy tide.

### 12.5 Bottom line

CN is **not** a minor contributor to the miRNA regulatory network in BRCA. The
CN-driven component reproduces ~**half** of the genuine bulk repressive signal
(attribution 0.51) **even after removing purity, HRD, target-CN, broad aneuploidy, and
proliferation together**; it is a **set / arm-scale** effect (combined > single,
broad ≥ focal, all surviving the full stack); and it is **luminal- and
immune-hot-concentrated** while **Basal's network is comparatively expression-driven**.
Crucially, **proliferation is a *suppressor* (it masks CN, not mimics it)** — adjusting
for it *raises* CN-driven detection (cohort 665→675; Basal 127→262; Her2 17→72), the
same MH-17 structure — so the earlier "proliferation/metabolic shadow" caveat is
resolved: the copy-number layer is a genuine upstream driver of the network.

---

## 13. Reproducibility

Key commands:

```bash
.venv/bin/python3 -m mirna_hallmark.mirna_locus_cnv --force-rebuild
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_cnv_genome_maps
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_overlap --workers 8
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_overlap --queues-only
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_overlap --case-plots-only
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_chr19_megacluster_gallery
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_locus_genome_dispersion
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.mirna_cnv_subtype_depth
.venv/bin/python3 -m mirna_hallmark.analyses.pressure_dev.pressure_sensitivity
.venv/bin/python3 -m mirna_hallmark.hybrid_pressure
.venv/bin/python3 -m mirna_hallmark.analyses.misc.target_combined_anticorr
.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.cn_dosage_attribution   # §12 CN→network attribution
```

The full rebuild extracts ASCAT3 segment CN for 1,049 primary participants.
Subsequent map and depth analyses reuse
`mirna_hallmark/output/mirna_locus_cnv/tables/sample_entity_cnv.tsv.gz`.
