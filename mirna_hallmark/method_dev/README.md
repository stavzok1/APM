# `method_dev/` — method & conceptual development

Working area for methodological/conceptual development in `mirna_hallmark` (distinct from the
result-reporting docs in `docs/`). Organized into **four arcs**, each in its own subfolder:

## Layout

```
method_dev/
  aggregate_pressure/   force-vs-abundance arc — design log, findings card, followups (docs)
  arm_expression/       canonical detectability floor — which arms are too silent to act (docs)
  site_ladder/          target-SITE confidence: seed-scan → filter ladder → experimental validation (code+data+doc)
  landscape/            genome loci + family / hub / competition figures (code+data)
  figures/              all rendered figures (+ figures/families/) — shared output dir
  README.md             this index
```

Run any module from **repo root**: `.venv/bin/python3 -m mirna_hallmark.method_dev.<subfolder>.<module>`.
Subfolder docs: `aggregate_pressure/AGGREGATE_PRESSURE_FINDINGS.md` (reference card, start here) +
`AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md` (design log) + `FOLLOWUPS.md`; `arm_expression/ARM_EXPRESSION_FLOOR.md`
(+ `SILENT_ARM_REMOVAL.md`); `site_ladder/SITE_FILTER_LADDER_PLAN.md`.

## Cohort overview — "what we deal with" (`landscape/cohort_overview.py` → `figures/cohort_overview.png`)

Orientation figure, 9 panels. **miRNA expression** (top): **A** rank-abundance — all 2,236 measured arms vs
the 716 HE-edge arms we model (both shown; **top 10 arms = 74% of the RISC pool**, top 50 = 97%; the modelled
arms span the whole range); **B** tier *composition* (100%-stacked) of all arms vs HE arms vs HE **edges** —
**silent is 39% of HE arms but only 11% of edges** (robust arms carry 62% of edges; silent arms are
low-degree); **C** per-arm median-abundance histogram (huge silent spike near 0, thin expressed tail above the
RPM-10 floor). **Gene regulation** (bottom, HE edges over 1,424 genes): **D** regulators/gene (in-degree;
median 2, mean 3.7, max 91=PTEN); **E** expressed regulators/gene (silent removed; mean 3.2) with the
convergence hubs labelled (PTEN, BCL2, CCND1, IGF1R, VEGFA, CDK6); **F** seed *diversity*/gene (# distinct
seeds; max 64); **G** site *load*/gene (total functional 7mer+ sites summed over the gene's regulators — counts
**sites not seeds**, cf. F; mean 4.2, max 151=PTEN). All four medians are exactly 2 (most genes are lightly
regulated; the means diverge in the tails). **Across-sample spread** (**H**, full-width): boxplots of the **top
50 most-abundant arms** over all 1,079 tumors with a dot per tumor — even the dominant arms swing **~10–100×
across tumors** (some, e.g. miR-375 / miR-205 / miR-141, far more), so abundance *rank* is a cohort-level
statistic, not fixed per tumor. **Seed-site 3′UTR position map** (**I**, full-width): family-collapsed
seed-site **instances** (same-seed arms counted once, each seed at every position) across all **1,181 sited
genes' 3′UTRs**, UTR-length-normalised (stop codon 0% → poly-A 100%) — sites **concentrate near the stop codon**
(~0.17/gene at the 5′ end) and dip mid-UTR (~0.09), the classic positional enrichment. Take-home: miRNA
abundance is extremely skewed and silent arms barely carry edges.

## Gene-regulator expression shape (`landscape/gene_regulator_expression.py` → `figures/gene_regulator_expression.png`)

How abundance distributes *within* a gene's regulator set (not just how many): **A** within-gene rank-abundance
(each gene's regulators normalized to its top one, median+IQR) — the **#2 regulator is typically ~15% of the
#1, #3 ~3%**; **B** dominance — the top regulator's share of total regulator abundance (**median 84%; 92% of
genes >50%-dominated, 56% >80%**); **C** effective vs nominal regulators (Simpson 1/Σpᵢ²) — **effective-N stays
~1.4 however large the in-degree** (a gene with 50 nominal regulators still has ~2 effective). Take-home:
a gene's incoming pressure is almost always concentrated in **one dominant regulator** — most nominal
regulators are abundance-negligible. Run `-m mirna_hallmark.method_dev.landscape.gene_regulator_expression`.

## Cross-state top-50 expression (`landscape/cross_state_top50.py` → `figures/cross_state_top50.png`)

**Four** stacked boxplot bands (dot per sample) comparing arm expression healthy→cancer, reusing
`cross_state_landscape._state_matrices()` (tumour / NAT / GTEx, GTEx already in TCGA-arm space):
- **Tumour** (TCGA primary, n=1,060) — top-50 arms, **dot per tumour coloured by PAM50 subtype** (LumA/LumB/Her2/Basal).
- **NAT** (tumour-adjacent normal, n=104) — the **same 50 arms** (tumour order), dot by matched-patient subtype;
  labels carry **rank-Δ and log2FC vs tumour** (NAT is the same RPM platform, so FC is meaningful) — 21/50 up in cancer.
- **GTEx ∥ tumour-50** (n=346) — the **same 50 tumour arms** shown in healthy GTEx, **raw log2(TPM+1)** (own scale,
  compare by rank); labels = rank-Δ vs tumour → which tumour-abundant arms are healthy-low (**31/50 cancer-elevated**,
  e.g. miR-21 Δr+406, miR-10b Δr+48).
- **GTEx OWN top-50** (n=346) — GTEx's most-abundant arms, raw TPM; labels green=more abundant in cancer / red=lowered.

GTEx bands are **raw TPM** (heights **not** comparable to the RPM panels — compare by **rank**); QN-onto-TCGA was
tried and **rejected** because it imposes TCGA's distribution *shape* (the GTEx rank-1 arm would inherit miR-21's
extreme dominance — GTEx's pool is genuinely flatter). All green/red up-down calls and Δ's use **platform-invariant
ranks** (GTEx vs tumour over 783 shared arms; NAT vs tumour over the full 2,236-arm TCGA platform).

**Finding:** of GTEx's 50 most-abundant healthy arms, **41 fall in rank in cancer, only 9 rise** — the known global
miRNA down-regulation in tumours. Reds are the tumour-suppressor families (let-7, miR-99/100, miR-30, miR-126,
miR-143/145, miR-148); greens include oncomiRs (miR-10b, miR-22, miR-200c). Run
`-m mirna_hallmark.method_dev.landscape.cross_state_top50` (~2–3 min; loads the cross-state matrices).

## Genome landscape (`landscape/he_edge_arm_landscape.py`)

Genome ideogram of HE-edge miRNA-arm loci: dot size/colour ∝ HE-edge degree (colorbar), **p/q-arms separated**
by a centromere dot, a right-margin column giving per-chromosome **miRNA-arm and locus counts**, top loci
labelled with non-overlapping leader lines. `--families top` (or `--families let-7,miR-17`) → one ideogram per
seed family (expressed arms filled, colour∝degree; **silent arms = grey hollow rings**, e.g. miR-302/520 = 6
expressed + 6 silent). `--bipartite "miR-29-3p"` → **bipartite gene-mirror** (miRNA arm loci right ↔ HE-edge
target-gene loci left, both at true positive genomic position; e.g. miR-29 → its ECM/collagen targets).
`--cluster "hsa-miR-17-5p"` → **cluster-COORDINATION** bipartite: a polycistron's member arms (coloured by
seed) → genes co-targeted by ≥2 *distinct* seeds. Figures → `figures/` & `figures/families/`; data
(`he_edge_arm_loci.tsv`, `gene_seed_diversity.tsv`) stay in `landscape/`. (Cluster tags kept in the TSV.)

**Coordination finding (from `--cluster`):** scanning all polycistronic clusters for genes co-targeted by ≥2
*distinct-seed* members, the cleanest coordination is the **miR-17~92 cluster (chr13:91,350,618–91,351,361)** —
10 arms / **8 distinct seeds** converging on 18 genes, with **PTEN, HIF1A, DNMT1 each hit by 5 different seeds**
and TGFBR2, BCL2L11/BIM by 3 — the canonical coordinated oncogenic program, recovered from the HE edges.
(Runner-up coordinated clusters: chr14 DLK1-DIO3 megacluster — diffuse; miR-143/145 chr5; miR-23a~27a~24-2 chr19; miR-200 chr1.)

**Non-textbook coordinated cluster:** `--cluster "hsa-miR-424-5p"` → the **miR-424/503/542 chrX cluster**
(chrX:134.54–134.55 Mb, 8 arms/8 seeds) converging on a coherent **G1/S cell-cycle module** — WEE1,
CHEK1, CDC25A, CCNE1, CCND1, CCND3, CCNF, ANLN, FGF2, FGFR1. (Caveat: convergence is dominated by the
miR-424-5p/miR-503-5p sister pair, near-identical seeds — a coherent-program example rather than the
seed-diverse convergence of miR-17~92.) Other interesting non-textbook hits: miR-183/96/182 chr7 →
FOXO1/FOXO3/PTEN/GSK3B (3 seeds each); miR-34b/c chr11 → BCL2/CDK4/CDK6/MET/MYC.

**Multi-seed 3'UTR hubs** (`--seed-hubs` → `gene_seed_diversity.tsv`): genes whose 3'UTR integrates the
most *distinct seeds* genome-wide — **PTEN (64 seeds / 91 miRNAs)**, BCL2 (52), CCND1 (48), IGF1R (42),
EGFR/STAT3 (35), VEGFA/CDK6 (34), HIF1A (33), CDKN1A (31)… the canonical convergence hubs.

**Site-count edge weighting** (`_site_counts`): coordination-figure edge width ∝ # predicted 3'UTR sites
(TargetScan rows per arm-gene pair; 484 HE edges are multi-site, e.g. **miR-200b/c/429 → ZEB1 = 5 sites**).

**Gene convergence hub** (`--gene-hub PTEN` → `figures/hub_PTEN.png`): the distinct seeds targeting a gene,
ordered by # 3'UTR sites; node size ∝ #arms/seed, edge width/colour ∝ #TargetScan sites. PTEN: 64 seeds /
91 miRNAs; top multi-site = miR-130a-3p (12), miR-26a-5p (~11), miR-106a/17-family (7).

## Gene competition figure (`landscape/pten_competition.py` → `figures/<GENE>_competition.png`)

Layered competition mechanism for one gene (default PTEN). **A** expression hierarchy of all regulators
(cohort-median RPM + cross-sample IQR, tier-coloured; PTEN: top 2 miRNAs = 50% of the RISC-competing pool,
81 expressed / 10 silent). **B** seed-structure bars (sites/miRNA, miRNAs/seed, site-type counts). **C** 3′UTR
competition map — expressed regulators in multi-member seed families placed on the 6,461-nt PTEN 3′UTR at
their seed-site positions. **Colour = seed** (same colour = direct seed competition); **dot size = seed-match
type** (8mer > 7mer-m8 > 7mer-A1, legend in panel). Site-pair regimes are distinguished: **red bands = S2
overlap ≤8 nt** (mutually exclusive — true competition), **blue bands = S1 proximal 8–40 nt** (cooperative
AGO). Shows the miR-17 family co-piling on one ~300 nt site, miR-25/32/92/363, miR-130/301, miR-29, miR-200
competing for shared loci. **C₀** (thin track above C) **seed-site density along the 3′UTR** — family-collapsed
**instances** (each seed counted at *every* position it matches; same-seed arms once), sliding 150 nt window —
peaks = binding-site hotspots, x-aligned with C. **D** position-ordered zoom on the densest hotspot (auto-detected: PTEN ~2255–2345
nt = **5 distinct seeds / 12 miRNAs**) — seed-match footprints as bars (overlap shaded), revealing the 3
offset-but-overlapping seeds (CAGUGCA/AGUGCAA/GUGCAAA: miR-148/152, miR-130/301, miR-19 = 8 miRNAs on one
footprint) that C scatters across family blocks. Run `-m mirna_hallmark.method_dev.landscape.pten_competition [GENE]`.

**Dot accounting:** in C each dot = one **site** (one seed-match occurrence) for one miRNA, coloured by its
seed; same-seed miRNAs at the same position stack on different rows (visible, not hidden — 0 same-row
coincidences for PTEN). "single-seed regulators omitted" = miRNAs that are the **sole expressed carrier of
their seed** (seed-family size 1, no same-seed competitor), not single-location.

## Site-pair regimes (`site_ladder/site_interaction_schematic.py` → `figures/site_interaction_regimes.png`)

Conceptual schematic of two miRNA sites on one 3'UTR, ordered by **Axis S (increasing overlap)**:
**S0** disjoint-distal (≳40nt, independent, log-additive `F=F_i·F_j`) · **S1** disjoint-proximal (8–40nt,
cooperating AGO, super-additive `F<F_i·F_j`) · **S2** overlapping (mutually exclusive, competition
`F≈min(F_i,F_j)`).

## Arm expression-detectability floor — CANONICAL (`arm_expression/ARM_EXPRESSION_FLOOR.md`)

`mirna_hallmark/arm_expression.py` (config `ARM_EXPRESSED_MIN_RPM=10` / `ARM_EXPRESSED_MIN_FRAC=0.01`;
`output/matrices/arm_expression_tiers.tsv`). Per-arm DETECTABILITY tiers (robust/conditional/silent),
**max-based**: KEEP if RPM≥10 in ≥1 tumor (not cohort-median, to keep context-specific induction — the 111
median-silent *spiker* arms are kept); REMOVE only if never reaches the floor. Of 789 HE-edge arms: **482
expressed / 307 silent**; **edges 5,219→4,628 kept** (591 dropped, 11% — silent arms are low-degree). Silent
arms can't occupy RISC (≤2% of pool) or show coupling (≈0 variance); a **noise filter, not a functional
verdict**. Headline unchanged on the kept universe. `filter_expressed_edges()` to apply; landscape tags every
arm's tier (silent = grey rings). Concise rationale: `SILENT_ARM_REMOVAL.md`.

## Site filter ladder (`site_ladder/SITE_FILTER_LADDER_PLAN.md`) — BUILT + TRIPLE-VALIDATED

Graded-confidence miRNA target sites via published filters (Bartel/Grimson/Friedman/Agarwal) — all code +
`.tsv.gz` data in `site_ladder/`: `site_filter_ladder.py` (L1 type · L2 AU+position · L3 3′-supp · L4
conservation PCT · APA-aware) → `utr_site_ladder.tsv.gz`; `site_genomic_l5.py` (UTR→GRCh38 lift + Manakov
site overlap) → `utr_site_ladder_genomic.tsv.gz`; `validate_ladder_experimental.py` (ENCORI+TarBase
v9+Manakov) → `ladder_experimental_support.tsv.gz`. **Validated:** experimental support rises monotonically up every rung
across 3 independent assays (edge-level ENCORI 10→78% / TarBase 49→88% / Manakov 9→35%) and at single-site
resolution (Manakov duplex coincidence 6→24%). TarBase v9 (`data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz`,
4.72M rows) ∩ HE = 57%. (POSTAR3 unreachable here; ENCORI+Manakov cover CLIP.)

## Universal 3'UTR seed-scan (`site_ladder/utr_seed_scan.py` → `site_ladder/utr_seed_sites.tsv.gz`) — BUILT

Fixes TargetScan's coverage gap (TS gives positions but only for **321 conserved-family arms**; ENCORI is
broader but interaction-level, no positions). The seed-scan reads each HE-edge gene's **MANE-Select 3'UTR**
(GENCODE annotation + GRCh38 genome FASTA, via a built-in faidx reader — no pysam) and scans it for each
miRNA's seed (canonical **6mer / 7mer-A1 / 7mer-m8 / 8mer**, from `mature.fa`), giving site **counts AND
positions** for every edge.

- **Coverage: 767 arms** (vs TargetScan's 321, **2.4×**); 3'UTRs for 1,420/1,424 genes (MANE).
- All **5,219 HE edges scanned**; **3,795 (73%) have a functional 7mer+ site** (1,429 multi-site, max 27);
  **348 carry only a 6mer**, and **1,076 (21%) have NO canonical seed site at all** on the MANE 3′UTR — these
  miRTarBase-functional edges must act via non-canonical / ORF / non-MANE-isoform sites or indirectly.
- **Validated:** Spearman ρ=0.62 vs TargetScan site counts (2,139 shared edges); confirms miR-200b/c/429 →
  ZEB1 = 6 sites (3× 8mer). (Counts differ from TS — TS pools 6mers + UTR isoforms; we count 7mer+ on MANE.)
- Output columns: `miRNA, gene, utr_len, n_6mer, n_7mer_plus, n_8mer, site_pos`. The **`site_pos`** (UTR
  coordinates of functional sites) is what enables real **S0/S1/S2** classification (site-pair separation/overlap).

## Seed-less edges — the 1,076 with no canonical site (`site_ladder/no_seed_edges.py` → `figures/no_seed_edges.png`)

Characterises the **1,076 HE edges whose miRNA seed matches nowhere on the target MANE 3′UTR** — real or
artefactual, and does an alt 3′UTR isoform rescue them? 6 panels: **A** experimental support
(ENCORI/TarBase/Manakov/any) vs sited; **E** **region rescue**; **F** support by rescue region; **B** MANE-UTR
length; **C/D** top miRNAs/genes. **Verdict — far weaker:** **ENCORI AGO-CLIP just 3% vs 36%** (no physical site →
no AGO occupancy), and only **37% have ANY experimental support** (vs 70%). **Region rescue (GENCODE — every
transcript's 3′UTR isoform + MANE CDS + 5′UTR):** **32% (348) gain a canonical seed in another region** — 15%
(160) in a **non-MANE 3′UTR isoform** (APA), **9% (94) in the CDS/ORF**, 2% (17) in the 5′UTR (+7% only a 6mer).
The alt-3′UTR and CDS rescues are **better-supported (49% any evidence each, vs 33%** for seedless) → they find
real sites; the 5′UTR rescue (29%, n=17) is marginal. **68% (728) stay seedless in every region** — the
genuinely non-canonical / indirect set. The no-seed edges sit on **shorter MANE UTRs** (median 1,176 vs 2,071 nt) and concentrate on **atypical-seed
miRNAs** (miR-375/137/320a) and **studied hubs** (HIF1A/MYC/DNMT1/PTEN/EZH2). Caches: `_noseed_evidence.tsv.gz` +
`_noseed_region.tsv.gz` (`--force` recomputes); per-edge region + flags in `no_seed_edges.tsv.gz`.
| `figures/he_edge_arm_landscape.png`, `figures/families/*.png`, `he_edge_arm_loci.tsv` | Landscape + per-family figures + per-locus table. **773/789 arms mapped (98%) = 100% of total HE-edge degree**: MirGeneDB loci (`mature_loci.csv`, the curated subset) primary + **GENCODE v49 miRNA-gene fallback** (1,879 miRNA genes; recovers the 232 non-conserved arms MirGeneDB excludes — in-repo, no miRBase gff3 needed). The 16 still-unmapped are name artifacts (comma-joined resolver outputs, mouse-style/ambiguous names). |

Headline of the arc: at the gene aggregate, **edge weighting and promiscuity normalization are
inert/harmful**; the one construction that beats summed abundance is **acquired abundance**
`max(a − h_healthy, 0)` (robust, proliferation-survived). See `AGGREGATE_PRESSURE_FINDINGS.md`.
