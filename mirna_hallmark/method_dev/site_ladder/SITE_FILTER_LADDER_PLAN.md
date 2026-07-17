# miRNA target-site confidence — filter-ladder plan

**Goal.** Turn the universal seed-scan (`utr_seed_sites.tsv.gz`, 767 arms, all seed matches with positions)
into a **graded-confidence** site set by applying a **ladder of published filters**, each rung trading recall
for precision. Report site/edge counts **and a validation metric at every rung** so we pick the operating
point empirically rather than by fiat.

**Operating principle (how we "believe it").** Each rung is a filter/score layered on the seed-scan sites.
Climbing = higher precision. We validate by recomputing, at each rung, the **functional signal** on the
curated ± sets — the oncomiR-TSG positive set should separate from the miRNA-poor negatives *more* as we
climb (and acquired-ρ / net-repression should strengthen). If precision rises monotonically up the ladder,
the filters are doing real work; the rung where it plateaus is the operating point.

All gold-standard methods are cited; thresholds follow the source papers.

## Inputs (availability)

| input | source | status |
|-------|--------|--------|
| seed-scan sites (type + UTR position) | `utr_seed_sites.tsv.gz` | ✓ built |
| MANE 3'UTR sequence | GENCODE v49 + GRCh38 FASTA | ✓ (faidx reader) |
| miRNA full mature seq (for 3'-supp) | `mature.fa` | ✓ |
| **conservation (PCT)** | `Conserved_Family_Info.txt` (col 11 `PCT`, Friedman 2009) | ✓ (conserved families only) |
| **APA / UTR shortening** | APAatlas `breast_pdui_per_gene.tsv` | ✓ |
| AGO-CLIP peaks / chimeric eCLIP | POSTAR3 / ENCORI CLIP; Manakov chimeric (Zenodo 14730307) | ⤓ download (L5) |
| TarBase v9 (edge support) | DIANA `Homo_sapiens_TarBase-v9.tsv.gz` | ⤓ downloading (orthogonal) |

## The ladder (rungs, increasing precision)

**L0 — Seed-match catalog (built).** All canonical sites: 8mer > 7mer-m8 > 7mer-A1 > 6mer
(Bartel 2009 site taxonomy). Recall ceiling; `n_6mer/n_7mer_plus/n_8mer` + positions already emitted.

**L1 — Site-type stringency + weighting** (Grimson 2007; Agarwal 2015).
- *Filter:* drop **6mers** (~50% chance background in a 6 kb UTR; weakly functional). Keep 7mer-A1/7mer-m8/8mer.
- *Weight:* per-site type coefficient = the mean context++ repression by type from Agarwal 2015
  (8mer ≫ 7mer-m8 > 7mer-A1). Store `type_weight`.

**L2 — Local context** (Grimson 2007 "AU / position"; folded into Agarwal 2015 context++).
- *Local AU content:* position-weighted AU fraction in the 30 nt up- and down-stream of the site
  (Grimson 2007 weighting ∝ 1/distance). Higher AU ⇒ more accessible ⇒ functional. Store `local_AU`.
- *Position-in-UTR:* (a) **exclude sites within 15 nt of the stop codon** (ribosome displacement, Grimson 2007);
  (b) `min_dist` = min(distance to UTR 5′ end, distance to 3′ end) — sites near UTR ends score higher; deep-interior
  sites in long UTRs are down-weighted. *Filter:* drop the <15 nt sites; *weight:* by AU and min_dist.

**L3 — 3′-supplementary pairing** (Grimson 2007).
- Check Watson–Crick pairing of miRNA **nt 13–16** against the UTR region 5′ of the seed match at the canonical
  offset (seed-match 5′ end − ~1..5 nt window). Flag/up-weight sites with ≥3 contiguous supplementary pairs
  (strongest for 8mer). Store `supp_pairs` (0–4) and `has_3p_supp`.

**L4 — Conservation** (Friedman 2009 PCT; Lewis 2005).
- Join each site to `Conserved_Family_Info.txt` on (miR family, gene, UTR position) → **`PCT`** (probability of
  conserved targeting). *Filter:* `PCT ≥ 0.5` = conserved (tunable; report the PCT distribution). Defined only
  for conserved families (~321 arms); other arms get `PCT = NA` → flagged **unconserved** (lower confidence,
  not hard-excluded, since non-conserved functional sites exist).

**L5 — Experimental occupancy (gold standard).**
- *AGO-CLIP:* intersect `site_pos` (lifted to genomic coords) with AGO-CLIP/eCLIP **peaks** (POSTAR3 or ENCORI
  CLIP clusters). Site `clip_supported` if it falls under a peak.
- *Chimeric eCLIP:* a site with a direct miRNA–target **chimeric duplex** read (Manakov) = `chimeric_supported`
  — the strongest single-site evidence ([[orphan-chimeric-eclip-validation]]).
- This rung is occupancy truth, not prediction; reported as the ceiling of belief.

## Definitions & clarifications (what each thing actually means)

**L1 `type_weight`.** *Not* a confidence — it's an **efficacy coefficient**: the mean repression a site of that
type confers, taken from Agarwal 2015 (TargetScan7 context++ regressions): **8mer 0.31, 7mer-m8 0.161,
7mer-A1 0.099** (log-fold-change magnitude). So one 8mer ≈ the repression of ~2 × 7mer-m8 ≈ ~3 × 7mer-A1. It
lets us *weight* sites by expected potency instead of counting them equally (an 8mer is worth ~3 weak sites).

**Why drop sites <15 nt from the stop codon (L2).** The terminating ribosome covers ~the first ~15 nt of the
3′UTR; a seed match sitting in that "ribosome shadow" is **sterically occluded** — RISC can't load there while
the ribosome terminates/recycles. Grimson 2007 measured it: sites within 15 nt of the stop are essentially
non-functional. So it's a mechanistic accessibility filter, not statistical.

**What L3 is (3′-supplementary pairing).** The seed (miRNA nt 2–7/8) gives specificity; some sites get *extra*
stability from the miRNA's **3′ end (nt 13–16)** base-pairing to the UTR just 5′ of the seed match. We revcomp
nt 13–16 and look for ≥3 contiguous matches in a 17-nt window upstream of the seed. Sites with it are stronger
(especially 8mers). It's the "does the back half of the miRNA also grab on" check.

**What PCT is.** **P**robability of **C**onserved **T**argeting (Friedman 2009; TargetScan's
`Conserved_Family_Info` col). For a site, the probability that its conservation across vertebrates is due to
**selection to maintain targeting** rather than background sequence conservation. PCT≥0.5 ≈ "conserved on
purpose." Defined only for the ~321 conserved miRNA families; non-conserved arms get PCT = NA (flagged, not
excluded — non-conserved functional sites exist).

**AGO-CLIP vs chimeric eCLIP.** Both crosslink Argonaute to RNA, but:
- **AGO-CLIP** (POSTAR3, ENCORI CLIP) sequences the *target* fragment under AGO → tells you **where AGO sits**
  (an occupancy peak) but **not which miRNA** guided it; you *infer* the miRNA from the seed match under the peak.
- **Chimeric eCLIP / CLASH** (Manakov) ligates the miRNA and its target into **one chimeric read** → you read
  **both partners directly**, so it names the exact miRNA–site duplex. Strongest single-site evidence, but sparser.

**Third validation table vs L5 — different resolution.**
- **Third validation** (`validate_ladder_experimental.py`, edge-level): per **(miRNA, gene) edge**, is the
  *interaction* supported *anywhere*? Columns = ENCORI (AGO-CLIP interactions), TarBase v9 (validated
  interactions), Manakov (any chimeric duplex for the pair). Answers "is this edge real."
- **L5** (`site_genomic_l5.py`, site-level): per **site position** (lifted to GRCh38), does the *predicted
  nucleotide location* fall under experimental binding? Columns = Manakov duplex, TarBase-CLIP coords,
  POSTAR3 AGO peak. Answers "is this exact site occupied."

**Why ENCORI and POSTAR-AGO are kept mutually exclusive.** They are the **same assay class (AGO-CLIP)** — partly
the same underlying experiments — so they are **not independent evidence**. The whole validation claim is that
the three columns are *orthogonal* assays. So each table uses exactly **one** AGO-CLIP representative: ENCORI at
the edge level (it ships interaction calls), POSTAR at the site level (it ships peak coordinates). Putting both
in one table would double-count AGO-CLIP and fake a third independent confirmation.

## Orthogonal modifiers (not rungs)

- **APA-aware UTR** (improvement #4; Mayr & Bartel 2009; APAatlas Hong 2020). Applied **before L0** as a UTR-mask:
  for each gene use the APAatlas **breast PDUI** (proximal-polyA usage). In genes with strong proximal APA
  (low PDUI) the distal 3′UTR is **not expressed in breast** → mask sites distal to the proximal polyA. Produces
  a "breast-APA-aware" site set alongside the MANE set; compare.
- **TarBase v9 edge support** (interaction-level, not site-level). Adds, per edge, whether an experiment supports
  it at all (and in which cell type/assay) — an orthogonal confidence axis to intersect with the site ladder.

## Output

`utr_site_ladder.tsv.gz` — one row per **site** with: `miRNA, gene, utr_pos, type, type_weight, local_AU,
dist_to_stop, min_dist, position_ok, supp_pairs, PCT, conserved, clip_supported, apa_retained`, and a derived
**`max_rung`** (highest ladder rung the site survives). Plus a per-edge rollup (`n_sites_Lk` for k=1..5) and a
per-rung summary table.

## Validation plan (the belief test)

For k = L1…L5, take the edges retaining ≥1 site at rung k and recompute on the curated ± sets:
1. **positive-vs-negative separation** of acquired-ρ / net-repression (Mann-Whitney; should ↑ with k);
2. **# net-repressed positives recovered** (should hold/↑ while negatives fall);
3. **internal consistency:** CLIP-support rate (L5) should rise across L1<L2<L3<L4 (higher rungs enrich for
   occupied sites);
4. **APA cross-check:** APA-masked sites should be the ones with weaker breast coupling.
A rung that improves (1)+(2) and is corroborated by (3) is real; the plateau = operating point.

## Build order / status

1. ✓ **L1+L2+L3 built** (`site_filter_ladder.py` → `utr_site_ladder.tsv.gz`): 6,030 functional (7mer+) sites /
   3,795 edges; 1,793 8mer / 2,631 7mer-m8 / 1,606 7mer-A1; 36% with 3′-supp pairing.
2. ✓ **L4 (PCT) built**: 4,100 sites carry PCT (conserved-family arms); 2,542 conserved (PCT≥0.5), median PCT 0.68.
3. ✓ **APA modifier built**: 988 sites in breast-APA-shortened genes; 222 distal sites dropped from the
   parallel breast-APA-aware set.
4. ⤓ L5 — download AGO-CLIP (POSTAR3/ENCORI) + Manakov chimeric; genomic intersect. (remaining, bonus)
5. ⤓ TarBase v9 — downloading (`Homo_sapiens_TarBase-v9.tsv.gz`); ingest as orthogonal edge support + validation.
6. Validation harness — per-rung ± separation (pending) + TarBase enrichment.

## First validation — internal coherence (no external data)

The independent rungs **agree on site strength**, which they wouldn't if they were noise:

| site type | n **sites** | % conserved | median PCT | % 3′-supp |
|-----------|---|-------------|-----------|-----------|
| 8mer | 1,793 | **50%** | **0.75** | 37% |
| 7mer-m8 | 2,631 | 40% | 0.64 | 36% |
| 7mer-A1 | 1,606 | 36% | 0.62 | 34% |

Rows = the **6,030 functional sites**, partitioned by seed-match type. The **`% conserved` column** is where
the monotonic claim lives: 8mer 50% > 7mer-m8 40% > 7mer-A1 36% (L1 strength ↔ L4 conservation agree). The
"35% vs 26%" is the **reverse cut of the same counts — not a table row**: of the 2,542 conserved sites **35%
are 8mers**, vs only **26%** of the 3,488 non-conserved (8mer share 1.35× higher among conserved).

**Definition of "conserved"** (used throughout): a site is conserved iff its (miR-family, gene) carries
TargetScan **PCT ≥ 0.5** (probability of conserved targeting, Friedman 2009). The 3,488 non-conserved are sites
with PCT < 0.5 **or** — the majority — whose miRNA family isn't in TargetScan's conserved-family set at all
(no PCT defined → counted as non-conserved, *not* dropped). So "conserved" is a positive-selection flag, not
mere sequence identity.

## Second validation — experimental AGO-CLIP (ENCORI) — PASSED, strongly

Belief test: do higher ladder rungs **predict experimental binding**? Per-edge ladder features vs ENCORI
AGO-CLIP support (P(CLIP-supported | feature)):

| ladder feature | n edges | % AGO-CLIP supported |
|----------------|---------|----------------------|
| weak (no 8mer / cons / 3′-supp) | 888 | **9%** |
| 8mer (L1) | 1,551 | 46% |
| 3′-supp (L3) | 1,748 | 37% |
| conserved PCT≥0.5 (L4) | 1,513 | **65%** |
| 8mer + conserved + 3′-supp (top rung) | 391 | **74%** |

**~8× enrichment from weak → top rung.** The ladder rungs predict experimental occupancy.

*Reading the table:* **n = edges** (unique miRNA–gene pairs; 3,795 total), **not sites** — site features are
rolled up per edge (`has_8mer` = the edge has ≥1 8mer site, etc.). Rows **overlap** (an edge that has an 8mer
*and* is conserved is counted in both) — they are feature subsets, **not a partition**; **weak** = edges with
none of {8mer, conserved, 3′-supp}. The "top rung" row is the intersection (8mer ∧ conserved ∧ 3′-supp).
*Why 3,795 and not 5,219?* the ladder only contains the **HE edges that have ≥1 functional 7mer+ site** to
score; the other **1,424** HE edges have no canonical 7mer+ site at all — a **site-availability** gap, **not** a
statement that they lack CLIP support. Breaking down those 1,424: only **348 (24%) even have a 6mer**; the
other **1,076 have no canonical seed match at all** on the MANE 3′UTR (likely non-canonical / ORF /
non-MANE-isoform sites, or indirect regulation). **Would adding 6mers help?** It would recover those 348 — but
we deliberately don't: a 6mer occurs **~50% by chance in a multi-kb UTR** (Bartel 2009), so as a *ladder rung*
it adds recall at a steep precision cost and would *flatten* the very enrichment the ladder demonstrates.
Reporting the 6mer count here (as context) is harmless; using 6mers as a confidence tier is not advised.

## Third validation — TRIPLE experimental cross-check (ENCORI + TarBase v9 + Manakov) — PASSED

`validate_ladder_experimental.py` → `ladder_experimental_support.tsv.gz`. Three **independent** experimental
layers, all monotonic up the ladder:

| ladder feature | n edges | ENCORI AGO-CLIP | TarBase v9 | Manakov chimeric |
|----------------|---|-----------------|-----------|------------------|
| weak (no 8mer/cons/3′-supp) | 888 | 10% | 49% | 9% |
| 8mer (L1) | 1,551 | 49% | 69% | 22% |
| 3′-supp (L3) | 1,748 | 39% | 67% | 20% |
| conserved PCT≥0.5 (L4) | 1,513 | 69% | 87% | 31% |
| **top rung (8mer+cons+3′-supp)** | 391 | **78%** | **88%** | **35%** |

Three orthogonal assays (AGO-CLIP occupancy, validated interactions, direct chimeric duplexes) **all** rise
monotonically with the ladder. The filters are calibrated against experimental truth, not just internally coherent.

*Same row convention as the second table:* **n = edges** (miRNA–gene; the same 3,795-edge rollup, overlapping
feature subsets — not a partition). **Resolution note:** this table and the second are **per-EDGE** ("is the
interaction supported anywhere?"); the **first-validation** table and **L5** below are **per-SITE** ("is this
exact site conserved / occupied?"). Same ladder, two resolutions.

## TarBase v9 — union / intersection with the HE edge set

A **pair** = one distinct **(miRNA, gene) interaction** — collapsing TarBase's many per-experiment rows down to
the unique miRNA→target relationships; this is the **same unit as one of our HE edges**, so pair counts are
directly comparable to edge counts.

TarBase v9: 4.72M rows → **1.31M distinct (miRNA,gene)** pairs (398k breast-tissue). Of our **5,219 HE edges**:
- **HE ∩ TarBase = 2,968 (57%)** experimentally validated; HE ∩ ENCORI 27%; HE ∩ Manakov 15%.
- HE-only (no TarBase) = 2,251; **HE ∪ TarBase (within our genes) = 171,157** — TarBase adds ~33× more
  candidate edges on our genes (the edge-universe expansion; ingest as its own evidence tier before using).

## L5 site-level — BUILT (`site_genomic_l5.py` → `utr_site_ladder_genomic.tsv.gz`)

Lifted all 6,030 MANE-3'UTR sites to **GRCh38 genomic** (UTR exon-block walk; coordinates verified — Manakov's
PTEN 3'UTR duplexes fall inside PTEN's GENCODE 3'UTR) and overlapped with **three independent CLIP datasets**:
Manakov chimeric (miRNA-specific duplexes), TarBase-CLIP site coordinates (miRNA-specific, HITS/PAR-CLIP/CLASH),
and POSTAR3 AGO-CLIP peaks (2.36M, AGO occupancy). *Does the predicted site coincide with experimental binding?*
— all three climb the ladder:

| rung | n | Manakov | TarBase-CLIP | POSTAR-AGO | any |
|------|---|---------|--------------|-----------|-----|
| 7mer-A1 | 1,606 | 6% | 24% | 48% | 55% |
| 7mer-m8 | 2,631 | 8% | 35% | 50% | 62% |
| 8mer | 1,793 | 16% | 46% | 55% | 68% |
| conserved | 2,542 | 18% | 59% | 63% | 78% |
| **8mer+conserved+3′-supp** | 352 | **24%** | **71%** | **70%** | **86%** |

**86% of top-rung sites coincide with experimental AGO/CLIP binding** — the ladder is validated at single-site
resolution against three orthogonal occupancy datasets. (POSTAR3 = `data/external/POSTAR/human (1).txt.gz`.)

## TarBase v9 at the miRTarBase HIGH-EVIDENCE level (apples-to-apples)

The 171k union was the raw (HTP-CLIP-dominated) set. Filtered to the **same evidence bar as our miRTarBase HE**
(low-throughput functional: Luciferase Reporter Assay + pSILAC):

| TarBase tier | rows | pairs | HE ∩ | HE ∪ (our genes) |
|--------------|------|-------|------|------------------|
| raw (all methods) | 4.72M | 1.31M | 2,968 (57%) | 171,157 |
| non-CLIP functional | 167k | 117k | 1,599 (31%) | 21,134 |
| **LTP functional (≈ miRTarBase HE)** | 9,819 | **5,256** | **1,135 (22%)** | **5,914** |

At matched evidence, **TarBase has 5,256 LTP-functional pairs genome-wide** (from 9,819 luciferase+pSILAC
rows) — comparable in *size* to our 5,219 HE edges, but **genome-wide** (all genes). The funnel to our overlap:
**5,256 genome-wide → 1,830 fall on our Hallmark genes → 1,135 match our exact arm-level edges** (the 22%). So
the modest overlap is partly our Hallmark-gene scoping (most of TarBase's 5,256 target non-Hallmark genes) and
partly different literature curation. The union (within our genes) expands only to **5,914** (+~700), not 171k.

**What about the other 78%? — database incompleteness, not contradiction.** Of our 5,219 HE edges: **22%
(1,135)** are in TarBase at matched LTP evidence; **57% (2,968)** are in TarBase at *some* evidence tier (the
extra ~35% TarBase carries only via CLIP/HTP, not LTP); and **43% (2,251) are in TarBase not at all** — these
come from miRTarBase's curated LTP experiments that TarBase's curators simply didn't include. Both
miRTarBase-HE and TarBase-LTP are **sparse samples of the same true LTP-validated interaction space** (each
curates a different slice of the literature), so 22% mutual agreement is concordance between two *independent*
curations — not 78% failure. The unmatched HE edges are "edges TarBase hasn't tested", not "edges TarBase refutes".
