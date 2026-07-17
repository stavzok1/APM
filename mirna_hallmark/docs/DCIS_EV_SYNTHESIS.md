# DCIS / pre-malignant / EV synthesis (session 2026-06-26)

Self-contained biology-first synthesis of the DCIS→invasive + EV arc (registry **MH-48..51**).
Written fresh (does **not** assume `BIOLOGICAL_SYNTHESIS.md` is current). Mechanical records:
`DISCOVERY_REGISTRY.md` (MH-48..51), `ANALYSES_CATALOG.md`, `ANALYSIS_RUN_LEDGER.md`,
`DATA_SOURCES.md` (cohorts + platform/annotation/overcome). Cohorts: GSE59247 (DCIS+IBC miRNA),
GSE59246 (paired mRNA), GSE297448 (MCF10A EV miRNA), GSE93740 (MCF10A **cellular** miRNA),
GSE162670 (LCM TNBC miRNA). All external; integrated by **rank/direction**, never absolute scale.

---

## Headline: a two-compartment, two-phase model of the in-situ→invasive transition

> **The epithelial oncomiR program is an early, pre-invasive lock-in; invasion is a late,
> microenvironmental, stroma-led *loss* program. The pressure framework has been reading the
> invasion biology in the wrong compartment.**

- **Phase 1 — pre-invasive, epithelial (done by DCIS):** oncomiR GAINS / brake-release are
  complete at the in-situ stage and are **robust to composition adjustment**.
- **Phase 2 — invasion, microenvironmental (continues past DCIS):** the miRNA **LOSS** axis is
  what keeps moving across the in-situ→invasive boundary — and it is dominated by the
  microenvironment: myoepithelial-barrier loss (miR-145/-143), endothelial/angiogenic
  (miR-126), stromal ECM/desmoplasia (miR-29c). Invasion is a tissue-architecture event
  (breach myoepithelium, remodel matrix), not a cell-autonomous epithelial step.

---

## Findings, each with its rigor status (raw → composition+site-adjusted)

### 1. Temporal asymmetry of gains vs losses — **ROBUST** (MH-48)
Same-platform DCIS-vs-IBC abundance, acquired direction from the GTEx-anchored TCGA axis,
set-level Wilcoxon (per-arm BH underpowered at n_IBC=14, so set-level is the statement).
- Acquired **GAINERS** (n=351): oriented rank-biserial −0.077 → **−0.089 adjusted**, signed-rank
  p 4e-6 → **3e-15** → *pre-invasively complete* (not composition-driven).
- Acquired **LOSERS** (n=57): +0.077 → **+0.170 adjusted**, p 6.9e-5 → **9.5e-7** → *loss continues
  to invasion*, and the adjustment **sharpens** it.
- Strongest nominal invasion movers are all downward, myoepithelial/stromal/immune (miR-133b,
  -145, -126, -140, -223) = the histological myoepithelial-loss signature of invasion.

### 2. Within-patient miRNA↔mRNA coupling — **broad replication survives, loss-leader enrichment is composition-driven** (MH-49, MH-51)
40 fully-paired DCIS+IBC samples (GSE59247×GSE59246), within-sample partial Spearman.
- Broad curated-edge replication vs TCGA: 59% → **54% adjusted** (still p=7.7e-4); FDR-neg 97→16.
- Loss-leader anti-correlation enrichment: p=0.0099 → **p=0.81 adjusted** → **collapses**
  (it was a stromal-composition signal — the MH-39 caveat, realized).
- Edges that **survive** adjustment: **miR-145-5p→COL5A1** (ρ_adj −0.53, q=0.047, strengthens),
  **miR-29c→SERPINH1** (−0.46, q≈0.10) and **→CDK6** (−0.42); basement-membrane collagens
  (COL4A1/2, COL15A1, COL7A1) are composition-attenuated.

### 3. miR-29c→ECM: a *compartment* result, triangulated (MH-51)
| System | stromal content | miR-29c | reads as |
|---|---|---|---|
| MCF10A line (GSE93740) | **zero** (pure epithelial culture) | ρ=−0.88 (q<0.001); −5p −0.75 | epithelium *can* do it |
| LCM epithelium (GSE162670) | reduced (microdissected) | ρ=−0.14 (q=0.70, n.s.) | not robust in real epithelium |
| Bulk DCIS (GSE59247) | full | raw grade −0.71 → **adj +0.25** | carried by composition |

Arm-specific: in the cell line miR-29**a/b rise** while miR-29**c falls**. ⇒ the in-vivo
miR-29c→ECM loss is **substantially a stromal-compartment phenomenon**, with epithelial-intrinsic
capacity visible only in the artificial line. **This is a relocation of the biology, not a
dismissal** (see "stroma reframe" below). The clean composition-independent *epithelial* grade
decliner is **miR-185-5p** (raw −0.70 → adj −0.75, q_adj=0.026) — underappreciated, worth chasing.

### 4. EV export is selective but a minority route, with two distinct "loss" mechanisms (MH-50, MH-51)
- Raw EV readcounts carry a **global cargo increase** (95% of miRNA rise; median trend ρ=0.60) →
  must analyse compositionally (CLR EV-share). The source paper reports the same global increase.
- Compositional selectivity: cellular LOSERS gain EV-share more than gainers (MWU p=0.027;
  25/55 losers FDR-confirmed exporters). **Direct** test (matched cellular GSE93740 + EV
  GSE297448, same lineage, no TCGA borrow): cellular & EV trends *positively* correlate
  (ρ=+0.15, p=0.002) → **co-accumulation dominates; only 43 miRNA directly exported**.
- **Two mechanisms of "loss":** **miR-205-5p / miR-140-3p are vesicularly exported** (miR-205
  shows a weak monotone trend but a large invasion-specific EV-share surge, DCIS→invasive step
  **+3.75**, recovered only by the 4-stage step decomposition; matches the paper's Rab27A story);
  **miR-29c is transcriptionally depleted** (down in cells AND EVs — *not* exported). The route
  matters mechanistically and for biomarkers.

---

## The stroma reframe (why "composition-driven" ≠ "uninteresting")

"Composition-driven" is statistical (the signal co-varies with cell fraction), not a biological
verdict. For miR-29c→ECM the stromal reading is arguably the **central** invasion mechanism:
- miR-29 is **the anti-fibrotic miRNA** — represses collagens (COL1/3/4/5), **LOX/LOXL2**
  (crosslinking → matrix stiffening, causally pro-invasion), SPARC, MMPs, fibrillins.
- The cells that build the **desmoplastic** invasive stroma are **CAFs**. miR-29-low CAFs →
  collagen/LOX-high → stiff, invasion-permissive matrix and basement-membrane breach.
- The pressure framework (MH-1..47) is implicitly **epithelial/tumor-cell-centric**; miR-29c→ECM
  is where the data points to a **compartment the framework never modelled**. The "confound" is a
  signpost, not noise.
- **Per-compartment question — RESOLVED (MH-52, Risom 2022 MIBI, 69,672 single cells):** it is
  **(i) regulatory CAF activation**, not (ii) a passive fraction shift. Total fibroblasts *decrease*
  normal→IBC, but the **CAF-activation ratio rises 0.03→0.44** (JT p=3e-4) and within fibroblasts
  COLI/FAP/SMA intensify with stage and are **CAF≫NORMFIBRO** (VIM flat = control); myoepithelial
  fraction collapses 0.24→0 (the miR-145/-143 program). Corroborated by the authors' Table_S4
  (miR-29 target collagens COL1A1/3A1/4A2/5A2 + SPARC CAF-high). The "stromal confound" was real
  desmoplastic CAF biology — **but the protein readout cannot say it is miR-29c-driven.**
- **miR-29c attribution RETRACTED (MH-54, GSE37527 — the non-circular test):** in 6 paired breast
  CAF vs NF (miRNA array, paralog-resolved), **miR-29c is NOT down in CAFs** (+12, p=0.69; miR-29a/b
  also flat/up). So the desmoplastic ECM up-regulation is **miR-29-INDEPENDENT**; the bulk miR-29c
  signal (composition-driven, MH-49/51) was a *correlation* (stromal fraction ↔ collagen), not a
  CAF-miR-29c mechanism. **The arms genuinely lost on CAF activation are miR-143/145/126/205/185**
  (p=0.031) — the canonical myofibroblast **miR-143/145 cluster**. ⇒ Net honest picture: the
  *desmoplasia* is real (MH-52), the *miR-29c→ECM-in-CAF mechanism* is not (MH-54); look to the
  TGFβ-axis for the driver. **miR-143/145 is the consolidated through-line** (MH-48 invasion
  loss-leader + MH-52 myoepithelial program + MH-54 CAF-activation-lost).

---

## Specific-resolution nominations (arm / edge) and where they sit on the timing

The broad biology is mostly confirmatory (next section), but the rigor pass produced concrete,
*specific* nominations — several genuinely under-studied — placed here on the in-situ→invasive
axis. Timing axis: **normal → (intra-DCIS grade) → DCIS → (in-situ→invasive step) → IBC.**

**A. EARLY / intra-DCIS epithelial severity (composition-INDEPENDENT grade decliners).** These
fall with DCIS grade *within* in-situ disease and survive composition+site adjustment — the
earliest epithelial-intrinsic progression signal we can see. Notably, adjustment did not only
*kill* signals (miR-29c) — it **revealed** two that composition was masking:
| arm | raw ρ(grade) | adj ρ | q_adj | note |
|---|---|---|---|---|
| **miR-185-5p** | −0.70 | **−0.75** | 0.026 | robust both; under-studied in DCIS |
| **miR-627-5p** | −0.27 | **−0.79** | 0.021 | composition-**masked** → revealed by adjustment (novel) |
| **miR-421** | −0.08 | **−0.76** | 0.026 | composition-**masked** → revealed (novel) |

*Orphan/target status (2026-06-26 check):* **miR-627-5p** is the cleanest orphan — only **1** curated
Hallmark target (→KDM3A, HYPOXIA), not on the acquired axis, GTEx-undetected, **no breast/DCIS
literature** (a colorectal miRNA elsewhere), and **ZERO TargetScan predictions** (absent from
TargetScan's miRNA set entirely — non-conserved; cf. miR-145-5p 1027, miR-29c-3p 1463 targets). So
it is *mechanistically a blank*: a genuinely novel but un-anchored **phenotypic** nomination (could
equally be a real dark-matter miRNA or a noisy low-expressed probe — treat with caution; n≈18 DCIS).
**miR-421** is less orphan — 11 curated targets incl.
**CDH1 (E-cadherin)/CASP3/FOXO4** (APOPTOSIS/APICAL_JUNCTION/TGFβ), broader pan-cancer literature
but under-characterised in DCIS; the CDH1 link is apt for an epithelial grade-decliner. Both are
nominations from ~18 graded DCIS (small) — replicate before foregrounding.

**B. Composition-robust within-patient regulatory edges (MH-49, stage-agnostic backbone).** The
edges that survive within-sample + composition adjustment — the "real" regulatory core in this
independent cohort: **miR-21-5p→EGFR** (ρ_adj −0.69, q=0.001, strongest), **miR-182-5p→FOXO1 /
SNAI2 / NDRG1**, **miR-330-3p→PRRX1** (EMT TF), **miR-26a-5p→TGFBR1**, **miR-92a-3p→E2F1**, and the
one loss-leader survivor **miR-145-5p→COL5A1** (q=0.047). Most are *oncomiR/broad* edges, not the
stromal loss-leaders (which were composition-driven) — i.e. the curated network that replicates is
epithelial-regulatory, not the desmoplastic axis.

**C. LATE / invasion-coupled losses (the microenvironmental program).** Strongest DCIS→IBC
down-movers + the compartment events: **miR-145 / miR-143** (myoepithelial), **miR-126** (endothelial),
**miR-140-3p**, **miR-29c** (CAF/ECM), plus nominal **miR-133b** (rb −0.67), **miR-338-3p**, **miR-223-3p**.
These map onto the MH-52 single-cell events (myoepithelial loss; CAF activation).

**D. Two mechanisms of "loss" at invasion (MH-50/51).** *Vesicular export* — **miR-205-5p**
(invasion-specific EV-share surge, DCIS→invasive step +3.75) and **miR-140-3p** (the one directly
exported loss-leader); vs *transcriptional / stromal* — **miR-29c** (down in cells AND EVs; lost in
the CAF compartment, not exported).

**E. Paralog dissociation within a seed family (MH-51, MCF10A).** miR-29 is not monolithic: during
transformation **miR-29c-3p/-5p fall** (−0.88 / −0.75, FDR-sig) while **miR-29a/b-3p rise**
(+0.53…+0.61) — a within-family arm/paralog split (ties the seed-family non-identifiability theme,
MH-45: a shared seed does not mean a shared trajectory).

**Where on the timing, in one line:** oncomiR *gains* and the epithelial regulatory backbone (B) are
**set by DCIS**; the composition-independent epithelial decliners (A: miR-185/627/421) move **with
DCIS grade (intra-DCIS)**; the suppressor/stromal *losses* (C/D) are **invasion-coupled (latest)**,
split into exported (miR-205/140-3p) vs stromal-transcriptional (miR-29c).

## Known vs novel — and what dataset availability allowed

**Most of the *biology* is well established; the contribution is the miRNA-lens re-derivation,
the rigor, and the temporal/compartment framing.** Honest standing per finding:
- **Confirmatory (textbook):** "most dysregulation already present at DCIS" (Ma 2003 / Volinia 2012);
  miR-29→collagen/ECM/fibrosis (van Rooij, Roderburg); myoepithelial loss = the histological
  definition of invasion; CAF activation/desmoplasia at invasion (FAP/αSMA/collagen/LOX); selective
  exosomal TS-miRNA export incl. miR-205 (the GSE297448 source paper's own claim).
- **Genuinely new here:** (1) the **gain-early / loss-late temporal asymmetry** framing; (2) the
  **rigorous composition attribution** showing the invasion miRNA signal is *stromal* not epithelial,
  with adjustment both killing (miR-29c) and **revealing** (miR-627/miR-421) specific arms; (3)
  **export-is-minority** (matched cellular↔EV) + the miR-29c-not-exported vs miR-205/140-exported
  mechanism split; (4) the specific composition-robust edges (miR-21→EGFR etc.) and the miR-29
  paralog dissociation. Strength = **triangulated consistency**, not one new mechanism.

**What the data conceptually + practically allowed (and didn't):**
- *Allowed:* a full arc — timing (GSE59247) → within-patient coupling (GSE59247+59246 companion) →
  EV mechanism (GSE297448+93740 isogenic) → single-cell compartment (Risom MIBI) — each
  independent of the TCGA spine, integrated by rank/direction.
- *Did not allow:* per-patient DCIS↔IBC pairing (none exist); well-powered per-arm DCIS tests
  (n_IBC=14 → set-level only); and — the binding limit — **miRNA at single-cell/compartment
  resolution** (every compartment result is at the protein/target/effector level). The miR-29c-in-CAF
  claim is therefore *supported* (target-in-CAF + known miR-29 biology + composition triangulation)
  but **not directly measured**; that needs compartment-sorted small-RNA-seq or miR-ISH (WHATS_NEXT §5.8).

---

## Honest caveats (carry these with the findings)
- External cohorts, small n (DCIS n_IBC=14; isogenic 3/stage) → **direction, not population power**.
- Cross-platform: TCGA used for *direction* only (rank); v16→modern arm bridge (GSE59247 ships
  miRBase **r14** IDs); feature→gene map for GSE59246.
- **No matched synchronous DCIS+IBC pairs** (confirmed: title 0-overlap + GEO "48 arrays/41
  tissues, independent groups"); the only within-patient layer is miRNA↔mRNA.
- LCM is morphology- not compartment-resolved (composition reduced, not eliminated); "normal" is
  tumor-adjacent (field effect). Grade is confounded with subtype (grade-1=LumA, grade-3=Basal/Her2).
- Cell-line decline may be a culture/Ras-transformation program, not the in-vivo epithelial program.

---

## What to pursue (priority order)
1. **✓ DONE → MH-52.** Single-cell spatial DCIS (Risom MIBI) resolved the miR-29c→ECM compartment:
   regulatory CAF activation + myoepithelial loss (effector level). Next within this thread: the
   **non-circular miR-29c** measurement via small-RNA sc/spatial (WHATS_NEXT §5.8) — SiCmiR
   inference is circular for this specific axis.
2. **Stroma-aware extension of the pressure framework** — model CAF/myoepithelial/endothelial
   compartments instead of adjusting them away; the myoepithelial-loss program (miR-145/-143/-126)
   is the LATE invasion leader and deserves its own module.
3. **miR-185-5p** — the composition-independent epithelial grade decliner: targets, Hallmark links,
   is it a bona fide progression suppressor? Cheap in existing data.
4. **EV export mechanism + liquid biopsy** — Rab27A axis (GSE297447 EV-mRNA); do miR-205/miR-140-3p
   export signals appear in the plasma-EV cohorts as an invasion biomarker?
