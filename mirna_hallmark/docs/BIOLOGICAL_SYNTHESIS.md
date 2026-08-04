# Biological synthesis — `mirna_hallmark`

The distilled, biology-first surface of what this subproject has established so
far, and how each claim stands against the literature. This is the high-level
companion to the granular [`REPORT.md`](reports/REPORT.md) and the claim-by-claim
[`DISCOVERY_REGISTRY.md`](DISCOVERY_REGISTRY.md); forward-looking work is in
[`WHATS_NEXT.md`](PROGRAM_FORWARD_BOARD.md).

> **Refresh status (2026-06-26 — Normal-like-excluded + batch + TNRC6 gate).** The cohort is now
> PAM50-Normal-like-excluded (**n=1065**; cross-state tumour 1060) with sequencing-batch covariates;
> the cross-state landscape, coupling, CPTAC and Buffa numbers are settled in
> [`LANDSCAPE_REPORT.md`](reports/LANDSCAPE_REPORT.md) and **DISCOVERY_REGISTRY MH-53**. The per-PAM50 program
> coupling below (Basal 42/50 etc.) is the **M0 proliferation-adjusted `pressure_sensitivity`** metric
> (MH-30, within-subtype-z'd + CPE/HRD/proliferation/target-CN **+ TCGA sequencing-batch**-adjusted on
the M0 spine — via `robustness_checks._partial_batch`, same batch covariate as the rest of the spine), **re-run
> 2026-06-26 on the Normal-excluded cohort: Basal 42/50 holds, LumA/LumB rise 18/17→22/22** (key 8/8,
> 5/8, 5/8, 0/8). This is a *distinct* metric from the rerun's gated CPE+HRD per-PAM50 coupling
> (`hallmark_anticorrelation_by_pam50`: **LumA 17 / Basal 5 / LumB 5 / Her2 0**; raw spearman 9/7/7/0) —
> the two differ because `pressure_sensitivity` is within-subtype-z'd and proliferation-adjusted.

> **Methodological note on framing.**
> - **Pressure weighting (2026-06-17):** **Spine (M0)** is **`softmax_z_logrpm` +
>   `evidence_mass`** (`pressure_build`) — abundance share × cohort z × log2(RPM+1),
>   edge weight ÷ evidence-mass per miRNA. **Hybrid headline:** **M7** tiered fusion
>   with **`softmax_z` + `combined_mass`** (`PRESSURE_HYBRID_*`). **Sensitivity:** M8,
>   M11. Legacy grids in ``pressure_sensitivity/``.
> - **Enrichment (limitation §8 item 2):** read as a **ranking**, not absolute fold
>   magnitudes — miRTarBase curation bias inflates well-studied pathways.
> - **Multi-set overlap (limitation §8 item 6):** genes in many Hallmarks appear
>   in multiple `(miRNA, gene, hallmark_set)` rows; cohort totals dedup on
>   `(miRNA, gene)` but per-Hallmark tables are not mutually independent.
> - Anti-correlations: confounder-adjusted (CPE + HRD). **Normal-like set aside**
>   in PAM50 contrasts. The primary subtype readout is now **actual PAM50 tumor
>   cohorts** (LumA/LumB/Her2/Basal); cold/hot visibility archetypes are secondary.

---

## Headline

In TCGA-BRCA, high-evidence miRNA regulation concentrates on **signaling and
cell-state-decision programs** (TGF-β, WNT, NF-κB, EMT rank highest — by
**ranking**, not absolute enrichment folds; see limitation §8.2). AGO-gated
miRNA pressure tracks **intrinsic subtype / lineage**, not stage.

The most useful subtype result is now per-gene, actual-PAM50 resolution:
**Basal** concentrates miRNA pressure on EMT, immune, apoptosis and
tumor-suppressor nodes (CCN2/SDC4/VIM/BCL2L10/PTEN/BIM/IRF1/CDKN1A), driven by
miR-17/20a/106b/23a/130b-like regulators. **Her2** has a weaker TLR3/WNT5A/CENPF
/IL15 pattern. **LumB** leans E2F/MYC/WNT/NF-kB pressure (HUS1, HNRNPA1, AXIN2,
RELA), while **LumA** carries a luminal/differentiation pressure signature
(TFAP2C, WNT2, VCAN/NT5E). Normal-like is excluded.

Cold/hot visibility archetypes were checked as a secondary layer. They do **not**
produce a strong cold-vs-other-Basal Hallmark-pressure split; the more stable
biology is the PAM50 subtype structure above.

---

## 1. What miRNAs regulate (the targeting map) — by rank

Hallmark programs **rank** (high → low) for high-evidence miRNA targeting as:
TGF-β ≳ WNT/β-catenin > apoptosis ≈ NF-κB ≈ EMT > angiogenesis ≈ UV-down ≈
Notch ≈ G2M ≈ IL6/JAK/STAT3 ≈ PI3K-AKT-mTOR, with **OXPHOS, bile-acid metabolism
and DNA-repair at the bottom**. Interpretation: miRNAs preferentially buffer
plasticity/signaling networks rather than constitutive metabolic machinery.
**Literature: concordant** with the canonical miRNA-regulation map; under-coverage
of OXPHOS/DNA-repair matches the thinner direct-targeting literature.

## 2. The hub miRNAs and the program-specific regulators

The broadest high-evidence regulators are the textbook BRCA miRNAs — **miR-21,
miR-155, miR-34a, let-7, miR-125b, miR-16, miR-29, miR-124**. Unbiased
target-count attribution recovers the *mechanistically correct* program leaders:
**miR-29 family → EMT**, **miR-34a → WNT/β-catenin (and estrogen/NF-κB)**,
**miR-21 + miR-155 → NF-κB/inflammation**, **miR-17/miR-181a → TGF-β**.
**Literature: strongly concordant** — this convergence is the project's main
internal-validity check.

## 3. Does pressure repress? Yes for differentiation/metabolic programs; the rest is confound

23/50 Hallmarks show a significant **negative** pressure↔target-expression
coupling that **survives** CPE+HRD adjustment — cleanest for differentiation /
metabolic-identity programs (bile-acid, estrogen/androgen response, adipogenesis,
myogenesis, peroxisome, angiogenesis, coagulation). Positive couplings split:
**proliferation** (G2M, E2F) persists after adjustment (co-regulation, *not*
de-repression), while **immune** positives (IFN-γ, IL6/JAK/STAT3, allograft)
**collapse** under purity/PAM50 adjustment (infiltration artifact). This is the
central caveat for using the pressure score as a repression readout.

## 4. Subtype contrast (Normal-like set aside): luminal vs non-luminal

Grouping HER2+Basal as *non-luminal* vs LumA+LumB as *luminal*:

- **Differentiation/hormone programs are more repressively coupled in
  non-luminal** (Δrho most negative for peroxisome, estrogen-late/early,
  xenobiotic, androgen, adipogenesis) — consistent with active miRNA-mediated
  suppression of luminal/differentiation identity in basal-like disease.
- **Immune + proliferation couplings are more positive in non-luminal**
  (allograft, IL6/JAK/STAT3, G2M, E2F) — i.e. the co-expression/infiltration
  confounds are stronger where immune and proliferative content is higher.

## 5. ★ Actual PAM50 per-gene resolution (LumA / LumB / Her2 / Basal)

Module: `pam50_gene_resolution.py`. This is the primary subtype-resolution layer:
each subtype is compared against the complement of the other three tumor PAM50
groups (Normal excluded), with raw pressure, per-gene z-pressure, target RNA, and
top realized miRNA drivers.

### 5.1 Basal: EMT + apoptosis/tumor-suppressor pressure

The strongest Basal-vs-complement **standardized** pressure gains are EMT /
stromal-interface and apoptosis-related genes: MTCH2, MID1, EVI5, FCGR2B, CCN2,
TNFSF11, SDC4, BCL2L10, SOCS5, VIM. Raw hub-gene pressure also highlights the
tumor-suppressor/apoptosis panel (PTEN, BCL2L11/BIM, FOXO1, CDKN1A, IRF1).

Driver examples:
- **PTEN:** Basal drivers miR-20a-5p, miR-17-5p, miR-18a-5p; LumA drivers are
  negative-realized (low expression), explaining the large subtype split.
- **BCL2L11/BIM:** Basal drivers miR-17-5p, miR-92a-3p, miR-20a-5p.
- **IRF1:** Basal drivers miR-106b-5p, miR-130b-3p, miR-23a-3p — matching the
  APM D76 route.
- **CDKN1A/p21:** Basal drivers miR-17-5p, miR-20a-5p, miR-106b-5p.

Interpretation: Basal tumors carry a broad miRNA pressure pattern against
immune priming, apoptosis, cell-cycle brakes and tumor-suppressor nodes. This is
more informative than the cold/hot archetype comparison, because it follows the
actual intrinsic subtype biology.

### 5.2 Her2: weaker, focused immune/developmental pressure

Her2-specific pressure is more modest. Top standardized genes include **TLR3**,
FXN, CENPF, DLL1, WNT5A, COL6A2, IL15, VAV2, CDK5R1, RPS6KA5. This suggests a
focused immune/developmental and mitotic interface, not the broad basal
tumor-suppressor pressure.

### 5.3 LumB and LumA: proliferative vs luminal/differentiation pressure

LumB gains concentrate on E2F/MYC/WNT/NF-kB-adjacent genes: **HUS1**, HNRNPA1,
PTGER2, TSPAN1, AKR1C2, HDAC3, SHMT2, ENG, AXIN2, RELA. LumA gains include
**TFAP2C** and WNT2 plus stromal/luminal identity genes (VCAN, NT5E, CD99).

Interpretation: LumB pressure follows proliferation and stress adaptation,
whereas LumA pressure tracks luminal/differentiation identity. This split is
more specific than a pooled luminal group.

## 6. Visibility-archetype contrasts (secondary; cold/hot within subtype)

Module: `visibility_archetype_contrasts.py` — joins APM
`visibility_archetypes/archetype_assignments.tsv` (cold_silenced / hot_IFN_MHCI).

### 5.1 All 50 Hallmarks — cold_Basal vs hot_Luminal (pressure delta)

Largest **cold_Basal > hot_Luminal** gated-pressure gaps (target mRNA often
*similar or higher* in hot_Luminal → repressive-load interpretation):

| Hallmark | cold_Basal | hot_Luminal | Δ pressure | target cold | target hot |
|----------|----------:|------------:|-----------:|------------:|-----------:|
| PI3K_AKT_MTOR | +3.51 | −0.31 | **+3.82** | 2.20 | 2.42 |
| NOTCH | +3.26 | −0.13 | +3.40 | 2.23 | 2.37 |
| **APOPTOSIS** | +3.16 | −0.19 | **+3.36** | 2.29 | 2.53 |
| TGF_BETA | +3.02 | −0.29 | +3.31 | 2.40 | 2.59 |
| UV_RESPONSE_DN | +3.19 | −0.12 | +3.31 | 2.13 | 2.39 |
| IFN_GAMMA | +2.98 | −0.14 | +3.12 | 2.09 | 2.51 |
| **P53_PATHWAY** | +2.53 | −0.14 | +2.66 | 2.23 | 2.42 |
| TNFA_NFKB | +2.50 | −0.12 | +2.63 | 2.20 | 2.44 |

**Observation.** The cold-Basal miRNA tone is **not immune-only** — apoptosis,
p53, PI3K-AKT-mTOR, and Notch rank alongside interferon/inflammatory programs.
This is **novel resolution** beyond the bulk basal/luminal split.

**Within basal only** (cold_Basal vs other_Basal): Hallmark-pressure deltas are
**near zero** (max Δ≈+0.03). The archetype signal is therefore a
**cold-Basal-vs-luminal** (especially vs hot_Luminal) axis, not a fine split
among all basal tumors at the Hallmark-pressure level.

### 5.2 Per-gene route panel (IRF1-style resolution)

| gene | cold_Basal pressure | hot_Luminal pressure | cold target RNA | hot target RNA |
|------|--------------------:|---------------------:|----------------:|---------------:|
| **IRF1** | +5.58 | −0.83 | 2.00 | 2.55 |
| **PTEN** | +47.2 | −4.1 | 2.15 | 2.50 |
| **FOXO1** | +20.4 | −2.9 | 1.93 | 2.32 |
| **BCL2L11/BIM** | +24.7 | −3.3 | 2.16 | 2.53 |
| **CDKN1A/p21** | +17.9 | −3.0 | 2.70 | 2.96 |
| **TP53** | +6.7 | −0.2 | 2.54 | 2.72 |

**Interpretation.** Cold-Basal carries high repressive miRNA load on
tumor-suppressor/apoptosis nodes **with lower target expression** for several
(IRF1, FOXO1, BIM) — consistent with active miRNA-mediated suppression, not
merely co-expression. **Yes, there is per-gene resolution like IRF1** for the
whole panel; see `gene_route_focus.tsv` + `gene_route_top_mirnas.tsv`.

**Driver miRNAs (realized load, cold_Basal vs hot_Luminal):**
- **IRF1:** miR-130b-3p, **miR-106b-5p**, **miR-23a-3p** (APM D76 route confirmed;
  miR-23a realized load +9.6 cold vs −2.8 hot).
- **PTEN / BCL2L11:** **miR-106b-5p**, miR-17-5p, miR-20a-5p dominate in cold_Basal;
  same miRNAs are **negative realized load** in hot_Luminal.
- **FOXO1:** miR-1271-5p, miR-96-5p, miR-15a-5p (weaker edges; smaller structural weights).

**Literature:** PTEN/BIM/FOXO1 miRNA repression in TNBC/basal is known; tying it
specifically to the **cold_Basal archetype** (not all basal) is the novel piece.

### 5.3 Bulk PAM50 contrast (subtype_contrasts) — cohort-wide lineage view

Still available: luminal (771) vs non-luminal (290), Normal aside. This captures
**lineage-wide** basal skew but **confounds** cold/hot within subtype. Use it for
cohort-wide differential gene pressure; use §5.1–5.2 for translational archetypes.

---

## 7. Differential gene pressure (bulk basal signature)

The genes most heavily miRNA-pressured in basal vs luminal (raw, high-confidence
multi-miRNA hubs) form a coherent program: **PTEN, FOXO1, CDKN1A (p21), BCL2L11
(BIM), TGFBR2** (tumor-suppressor / pro-apoptotic / differentiation), alongside
**CCND1, HIF1A, STAT3, VEGFA, BCL2** (proliferation / hypoxia / survival /
angiogenesis). Standardized (per-gene SD) contrasts add an estrogen-axis signal
**more pressured in luminal** (TFAP2C, AGR2) vs **EMT/immune genes more pressured
in non-luminal** (CCN2/CTGF, SDC4, CASP7, FCGR2B).

**Interpretation.** Basal tumors concentrate miRNA pressure on the very nodes
whose loss defines basal biology — PTEN/FOXO1/p21/BIM suppression is a recognized
basal/TNBC hallmark. **Literature: concordant** (PTEN-low, p21/FOXO dysregulated,
BIM-suppressed TNBC). Caveat: raw pressure scales with miRNA in-degree (hub genes
dominate); the standardized table corrects this but ties low-degree single-miRNA
genes — use both views together.

## 8. miRNA-level attribution per Hallmark (and how it flips by subtype)

Decomposing each Hallmark's pressure into per-miRNA structural load + group-realized
pressure surfaces the dominant regulator and its subtype behavior, e.g.:

- **IFN-γ response:** miR-146a-5p, miR-17-5p, **miR-106b-5p**, miR-221-3p —
  all realize far more pressure in non-luminal (the immune-pressure axis above).
- **EMT:** miR-29b/a-3p realize *more* in non-luminal, but **miR-29c-3p realizes
  much less in basal** (lower miR-29c in basal) — a within-family split.
- **NF-κB (TNFA):** miR-27a-3p strongly non-luminal-skewed, while **miR-199a-5p
  is luminal-skewed** — opposite-direction regulators within one program.

**Interpretation.** The *identity of the most impactful miRNA on a program is
subtype-dependent*, and within a seed family different arms can move in opposite
directions (miR-29a/b vs miR-29c). This is the granular structure the cohort-level
view hides.

## 9. Copy-number / dosage and AGO capacity (supporting layer)

- ~31% of mature arms are dosage-coupled (locus CN → expression); the most
  coupled (**miR-151a, miR-30b/d, miR-320a**) sit on **8q24**, so part of the
  "pressure" signal in aggressive subtypes is a structural CNV readout.
  **Literature: concordant** (recurrent 8q gain).
- **AGO2** is amplified across aggressive subtypes (LumB/Basal CN ~5/4.9, 80%
  gain; Normal-like near-diploid). The bounded RISC gate **rescales but rarely
  reorders** results, so AGO availability is a **motivated modulator**, not a
  demonstrated rate-limiter here.

## 10. ★ Robustness: proliferation + curation-bias controls (grant Aims 1–2)

Module: `robustness_checks.py`. Two pre-registered reviewer objections, both
tested on the real cohort (Basal n=197).

### 10.1 Proliferation is a *suppressor* confound, not the driver — and not an over-adjustment artifact (MH-17, MH-19)

The hub is the proliferation machinery (miR-17~92 / miR-106b~25 are MYC/E2F-driven
and Basal is the most proliferative subtype), so the obvious objection is that
elevated hub-miRNA + depressed tumor suppressors are both just downstream of
proliferation. We added a proliferation covariate to the CPE+HRD ladder. The
negative coupling does **not** weaken — it **strengthens**, because proliferation is
a positive confound masking repression. Because our first proxy (E2F+G2M metagene)
shares its upstream driver (E2F/MYC) with the hub, we ran the test under **three
proxies** to rule out over-adjustment:

| program (Basal) | CPE+HRD | +E2F/G2M | +MKI67 alone | +ortho (E2F/MYC-stripped) |
|-----------------|--------:|---------:|-------------:|--------------------------:|
| APOPTOSIS | −0.16 | −0.29 | −0.24 | **−0.33** |
| P53_PATHWAY | −0.18 | −0.26 | −0.23 | **−0.30** |
| EMT | −0.19 | −0.24 | −0.22 | **−0.26** |
| TNFα/NF-κB | −0.16 | −0.26 | −0.24 | **−0.25** |
| PI3K-AKT-mTOR | −0.03 | −0.23 | −0.19 | **−0.22** |
| IFN-γ | −0.06 | −0.22 | −0.19 | −0.15 |

Key-program survival (neg & p<0.05): **8/8 (E2F/G2M), 6/8 (MKI67), 8/8 (ortho)**.
The **independence claim leads on MKI67-alone** — a single histological proliferation
marker that does not partial out a 200-gene E2F-target program. The ortho metagene
(mean z of G2M_CHECKPOINT + MITOTIC_SPINDLE after removing all E2F/MYC targets, 271
genes) is a *largely* E2F/MYC-independent cross-check (G2M genes remain partly
E2F-influenced), not a fully orthogonal one. The strengthening reproduces under both,
so the effect is not a shared-E2F over-adjustment artifact.

> **Magnitude caveat (read the numbers as illustrative, not effect sizes).** Part of
> the deepening is **variance deflation**: when target expression co-varies strongly
> with proliferation (e.g. PI3K target +0.89), partialling proliferation removes most
> of the target's variance and the partial-r is read off a thin residual, which
> inflates |partial-r| through the denominator. The deepened values are therefore
> **proxy-dependent** (IFN-γ −0.22→−0.15 between proxies is the visible tell). The
> **stable, reportable** quantities are **sign, FDR-survival, and survival count**
> (8/8 key programs under E2F+G2M and ortho; 6/8 under MKI67); the magnitudes are illustrative of direction only.

- **Sign structure — program level (`prolif_sign_structure.tsv`, ortho proxy, Basal):**
  proliferation is positively associated with **both** pressure (≈+0.05…+0.13) and
  target expression (**+0.36…+0.89**; PI3K +0.89, TGF-β +0.80, NOTCH +0.76). High-
  proliferation basal tumors co-express these target genes, masking the miRNA
  repression — exactly the suppressor structure required for adjustment to deepen the
  negative coupling. This is *shown*, not asserted.
- **Sign structure — gene level (the headline robust claim, confirmed individually):**
  the suppressor sign (prolif positively associated with the hub arm **and** the target)
  holds for the robust genes, but is weaker at single-gene than program resolution, as
  expected — `rho(prolif, target_expr)`: PTEN **+0.54** (ortho) / +0.24 (MKI67), FOXO1
  +0.34/+0.05, TGFBR2 +0.37/+0.07, CDKN1A +0.27/+0.08, VIM +0.24/**−0.02** (MKI67); the
  hub arm side (miR-17-5p/20a-5p) is +0.16…+0.38 across proxies. So PTEN shows the clean
  suppressor sign on its own; p21/TGFBR2/VIM are weak-positive (and VIM is ~0 under
  MKI67) — single tumor-suppressor genes track proliferation less than a growth-program
  score, which is exactly why the robust claim rests on **sign + survival**, not on the
  deepened magnitude.
- **Per-gene routes, within Basal (ortho proxy):** PTEN (all three miR-17/20a/106b
  arms), CDKN1A (miR-106b-5p, miR-17-5p), TGFBR2 (both arms), VIM survive; **8/14
  routes (ortho); 9/14 (E2F/G2M and MKI67); 13/14 cohort-wide.**
- **Casualties (honest):** the **IRF1** route fails within Basal under every proxy
  (miR-23a-3p p≈0.4–0.6) and **BCL2L11/BIM** fails. The TSG/cell-state core is
  robust; the IRF1 immune route is not.

### 10.2 The hub is not a miRTarBase curation artifact — at the cluster level (MH-18, MH-19)

Evidence-weighting upweights well-studied oncomiRs, so we recomputed each miRNA's
basal regulatory load under three weightings of the **same** high-evidence edges:
evidence-weighted (current), binary/degree, and **sequence-based TargetScan
weighted-context++** (study-count-independent). The hub holds at the family level:

- **miR-17-5p ranks #1 under all three** weightings; miR-20a-5p #3/#3/#2;
  miR-106b-5p #9/#10/#5; miR-23a-3p #15/#15/#13. The miR-17~92/miR-106b~25 family
  (incl. miR-19a-3p, miR-92a-3p) tops the TargetScan ranking too.
- Load rank-correlation: evidence-vs-TargetScan ρ=0.65; binary≈evidence ρ=1.00
  (so evidence-weighting adds little beyond degree).
- **Caveat (per-gene):** TargetScan reassigns the *specific arm* per target —
  PTEN's strongest sequence-predicted repressors are miR-19a/b-3p (the miR-17~92
  PTEN arm) + miR-29; IRF1's are **miR-130b-3p / miR-301a/b**, not the nominated
  miR-23a-3p/miR-106b-5p. So the **cluster** is robust; specific single-edge
  nominations (especially IRF1) are predictor-dependent.

### 10.3 Per-family decomposition: how much, and is it basal-*specific*? (MH-20, MH-21)

We decomposed the basal pressure by miRNA family — both its **share of load** and,
critically, its **within-subtype coupling** (family-restricted per-sample pressure,
partial Spearman | CPE+HRD+proliferation, FDR within (family, subtype)).

- **Load share** (`family_load_share.tsv`): the **miR-17~92/106b~25/106a~363**
  family carries **22.1%** of total positive basal load (28.2% under TargetScan) and
  supplies **4 of the top-10 / 7 of the top-20** miRNAs. Next: **miR-23~27~24**
  (9.1%; contains the IRF1-route arm miR-23a-3p) and **miR-221/222** (4.0%).
- **Coupling** (`family_coupling_by_pam50.tsv`, `family_specificity_summary.tsv`):
  - The miR-17~92 family represses its targets in **every** subtype but **far most
    deeply in Basal** — median ρ Basal **−0.45** (44/50 programs) vs Her2 −0.31
    (37/50), LumB −0.22 (39/50), LumA −0.19 (40/50). So it is **basal-MAXIMAL, not
    basal-exclusive**: a constitutive oncomiR program **amplified** in Basal (deepest
    coupling + highest cluster expression, basal z ≈ +1.3).
  - **miR-23~27~24 and miR-221/222 are basal-*enriched* in expression but
    functionally *uncoupled*** — 0–1 of ~40–48 programs negative-significant in any
    subtype (median ρ ≈ 0, even +0.19 in Her2). So among the top-load families, the
    **miR-17~92 cluster is the only one whose load actually translates to repression.**
- The **all-miRNA aggregate** coupling is genuinely **basal-concentrated** (Basal
  **42/50**, LumB **22/50**, LumA **22/50**, Her2 **0/50** under the adopted
  logrpm+evidence_mass spine; key programs **8/8**, **5/8**, **5/8**, **0/8**) —
  but **LumB is non-trivial**, so do not call luminal coupling "absent."

**Interpretation.** Two distinct claims must not be conflated: (1) the *aggregate*
miRNA→repression coupling concentrates in Basal (true, Fig. 1A); (2) the *cluster's
own* coupling is basal-restricted (**false** — it is pan-subtype, ~2× deeper in
Basal). The defensible basal story is **amplification + functional uniqueness**: the
cluster is the single most concentrated *and* the only functionally-coupled top-load
family, and its repression is deepest where the cluster is highest (Basal).

**Bottom line for the grant.** The defensible, proliferation- and weighting-robust
proof of concept is a **miR-17~92 / miR-106b~25 program repressing a tumor-suppressor
/ cell-state node set (PTEN, p21, TGFBR2, VIM)** — multi-program, **basal-maximal**
(deepest coupling in Basal; the only top-load family that is functionally coupled),
and not explained by proliferation or study counts. The **IRF1 immune-priming route
is exploratory** (fragile to proliferation and reassigned by sequence prediction) and
should be validated, not assumed.

### 10.4 TargetScan orphans — breadth without a new Basal program layer (MH-23, MH-24)

We tested whether miRNA→target pairs **predicted by TargetScan but absent from the
miRTarBase high-evidence graph** show independent repression coupling.

- **Program level (MH-23):** under the **prior** softmax_z+degree panel, orphan
  TS-only pressure raised neg-sig counts (Basal 24→38/50) without adding Basal key
  programs beyond miRTar HE. Under the **adopted logrpm+evidence_mass spine**, full
  M0 and orphan program coupling both reach **42/50 Basal** (8/8 key) — orphan breadth
  inflation is **not** additive vs the refreshed spine. Orphan signal at program
  resolution still **tracks** curated pressure rather than a parallel Basal layer.
- **Gene level (MH-24):** aggregating orphan arms per hub within Basal, signal
  **clusters on a few nodes**: **CDKN1A** (8/14 CN-surviving arms; miR-301/130/454/499
  + let-7/98 family), **PTEN** (miR-148b), **TGFBR2** (miR-520d). Seed-family
  heatmaps show coherent blocks, not diffuse TS noise. Most zero-study,
  basal-concentrated arms sit on **p21/PTEN** — exploratory sequence-only candidates.
  **IRF1/BIM orphan arms do not couple.**

- **Subtype intersection (MH-25):** the **curated HE hub is constitutive** — only 4
  active miRNAs, 3 of 4 in all PAM50 subtypes, Jaccard ≥0.75, **zero Basal-only**. The
  **orphan TS layer** is what introduces subtype-private structure (52 active miRNAs;
  **5 Basal-only, 5 LumB-only, 3 pan-subtype**; Jaccard 0.27–0.47, luminal pair
  highest, Basal most divergent). The apparent subtype specificity of orphans is
  therefore **private low-evidence arms**, not a coherent shared Basal program.
- **miRTarBase-HE curation-gap coupling queue (MH-26) — NOT a novelty claim.** Using HE
  membership as a *negative* filter on TargetScan (zero curated studies) + CN+prolif-
  adjusted coupling: **38/62** orphan hub relations are absent from miRTarBase HE; **19**
  couple. Cohort-wide the strongest sit **outside Basal**: **miR-183→FOXO1** (ρ −0.48),
  **miR-124→VIM** (ρ −0.46), **miR-383→IRF1** (ρ −0.42, Her2), plus the Basal **p21/PTEN**
  block (miR-499a/301a/301b/454/130b→CDKN1A; miR-148b→PTEN).
  **CORRECTION (all 19 edges PubMed-triaged):** miRTar-HE absence is a **curation gap,
  not literature novelty.** **14/19 are luciferase-validated elsewhere** (miR-183/135a→FOXO1,
  miR-301a/499/130b→CDKN1A, miR-301a→TGFBR2, miR-383/17/124→IRF1, miR-148b/454→PTEN,
  miR-30→VIM, miR-181b→BIM); **1 is likely indirect** (miR-124→VIM via EMT); **only 4 are
  genuine direct-edge gaps** — **miR-454→CDKN1A**, **miR-301b→CDKN1A** (no p21 luciferase;
  both validated on other targets), and **miR-24/miR-519d→IRF1** (prediction/network only).
  The queue's value is **breast-context / in-tumor confirmation** of under-curated edges,
  with only a small truly-uncharted tail (`mirtar_gap_coupling_queue.tsv`,
  `literature_status` column).

**Interpretation.** Orphans are worth mining for **follow-up validation** at specific
hub targets, but they do **not** overturn the curated-graph story or rescue the IRF1
route.

### 10.5 Combining miRTarBase HE + TargetScan into a pressure set — recommendation

Given MH-23/24/25, a **flat union** (HE ∪ TS-orphan) is **not** advisable: it inflates
program breadth and recruits subtype-idiosyncratic noise without adding a Basal-specific
layer. The principled combination is **agreement-weighted / confidence-tiered** pressure,
using TargetScan as an *orthogonal corroboration prior* on the curated graph rather than
a peer source:

- **Tier 1 — canonical (HE):** current evidence-weighted, AGO-gated pressure. Unchanged.
- **Tier 2 — cross-validated (HE ∩ TS):** pairs supported by **both** curated HE evidence
  and TargetScan context-score get **up-weighted** (two independent evidence types). This
  is the highest-confidence set and concentrates on the hub (miR-17~92 → p21/PTEN/TGFBR2),
  which is both curated and sequence-predicted. **Testable claim:** HE∩TS pressure should
  couple to target RNA *more deeply / more robustly* within Basal than HE-alone.
- **Tier 3 — exploratory (TS-only orphans):** kept **separate**, down-weighted, used for
  hypothesis generation only (luciferase/CLIP follow-up per MH-24), never merged into the
  primary pressure score.

The inputs already exist (`ts_weight` per pair from the orphan module; evidence weights
from the HE graph), so building an agreement-weighted Tier-2 set is feasible without new
data. The added value is **robustness on the existing hub**, not new breadth.

---

## 11. ★ Primary-tumour roadmap — composition-robustness + the outcome verdict (MH-57→72)

This arc took the framework to **primary TCGA tumours** and asked two questions: (i) does the coupling
survive *proper* composition control, and (ii) is the pressure *prognostic*. The answers are clean and
opposite — **robust regulatory biology, but not a prognostic classifier.**

### 11.1 Deconvolution is validated by two gold-standard methods (Phase 0; MH-70, MH-72)
TCGA bulk was deconvolved against the Wu primary-breast scRNA atlas. Naive NNLS and standalone CIBERSORT.R
**failed on epithelial** (5–24 % vs the expected ~60 %) — the platform gap, not the algorithm. **BayesPrism**
(tumour-aware malignant `key`, 2-round Gibbs) gave **63 % Cancer-Epithelial** and **CIBERSORTx S-mode**
(batch-corrected, user's Docker) gave a strong fit (R=0.805). The two **converge** on the trustworthy
compartments — **CAF cross-method ρ=0.77 (both ≈0.85 vs metagene); both pass the NAT positive control**
(NAT ↓Cancer-Epithelial, ↑Normal-Epithelial); both flag the rare types (PVL/Plasmablasts) as noisy. The
double-log bug in `load_rna` (it logs an already-log file) was found and fixed for the linear-scale inputs.

### 11.2 The framework's coupling is composition-ROBUST (Phase 1-refined; MH-71) — method-robust
Re-running the headline composition-retest with the **real deconvolution fractions** as covariates:
**96 % of the 1,013 BY-neg couplings stay negative (BayesPrism) / 94 % (CIBERSORTx) / 93 % (metagene, MH-57)**
— all three converge, so it is **not a deconvolution-method artifact**. The only attenuated class is the
**same miR-29→collagen/ECM edges** every time (COL5A1/2/3, COL1A1, FBN1, BMP1, ADAM12) — the known
stromal-composition confound (MH-54), now gene-resolved: **miR-29-family targets are 2× CAF-enriched**
(MH-59; COL3A1 is CAF, whereas miR-29→**LAMC2** is *epithelial* basement-membrane — one family, two
compartments). With proper composition control the framework is **more** robust, not less.

### 11.3 The outcome verdict — pressure is NOT prognostic; functional repression > abundance (MH-60→69)
Exhaustively tested across **every construct × resolution × gene-role × topology + multivariate**, in
**TCGA (DFI+OS) and METABRIC-miRNA/Buffa (DRFS)**:
- **Pressure magnitude is not prognostic** — 0 recurrence FDR across 3,368 features; the multivariate
  elastic-net and Random-Survival-Forest signatures **do not beat clinical** out-of-sample (DFI 0.59 vs
  0.64). miRNA-CNV dosage, pressure residuals, acquired pressure, program×role×topology — all null.
- **The durable outcome signals are KNOWN biology** — proliferation programs → RFS (METABRIC mRNA,
  q<0.001) and **miR-210 (hypoxia) → DRFS** (= Buffa 2011), all on the expression/mRNA side.
- **Functional repression beats abundance** (MH-67): whether a miRNA *actually represses its targets*
  (anti-corr) out-prognoses how much miRNA is present (beats abundance 62 %); the **functional TS-miRNA
  panel beats the abundance panel** (FDR-sig); **escape** (abundant-but-not-repressing) → worse outcome
  and opens the only TCGA crack (miR-93/103a).
- **miR-200/EMT deep-dive (MH-69):** the signal runs through **LOX** (hypoxia-induced ECM crosslinker,
  q=0.022) and **BCL2**, **not** the canonical ZEB1/2 — matrix remodeling, converging with the
  miR-210/hypoxia axis and the DCIS miR-29→collagen theme. Target-set role-weighting and per-miRNA
  pressure↔expression coherence did **not** add.
- **One lead for the powered cohort:** regulatory **coherence → OS** (more-intact repression → worse OS,
  q=0.004, TCGA) — OS-specific (recurrence null in both cohorts), untestable in Buffa (DRFS only) →
  needs **METABRIC-full (EGA-pending)**. **CPTAC ruled out** (13 deaths / 0 recurrence — power floor).

### 11.4 Framework identity (the bottom line of this arc)
The miRNA→Hallmark **coupling is real and composition-robust** (survives gold-standard deconvolution,
two methods) — a **regulatory-biology** tool. It is **not a prognostic classifier**: pressure adds no
out-of-sample survival information beyond clinical, and the durable miRNA-prognosis (miR-210, miR-200→LOX,
proliferation) is known hypoxia/ECM biology found on the *functional* (not abundance) axis. The honest
remaining lead is coherence→OS, for METABRIC-full. **A target-aware *functional* pressure construct (MH-73/74/75) does beat full clinical *within* METABRIC (+0.06, unit-robust) — but it does NOT generalize to TCGA when frozen (MH-76: +0.002/+0.006), so it is METABRIC-specific (likely platform), not yet a generalizable biomarker.**

### 11.5 Pressure earns its complexity by ATTRIBUTION, on robust textbook genes (MH-80)
"Does the constructed pressure beat just summing the regulators' abundance?" — on aggregate the gain is
modest (MH-79, +12.6% from a NAT reference), but the *value* is concrete and per-gene. The softmax-**share
attribution** re-ranks which arm dominates a gene's regulation (459/758 genes flip), and on **94 genes** the
arm it **surfaces** — masked behind a more-abundant arm by a naive abundance-sum — is a **literature-backed
regulator**, on a gene whose net-repression **survives composition+proliferation adjustment (75/94)**. The
clean exemplar is the user's reference case: **ABCA1**, where the share surfaces **miR-33b-5p** (the textbook
cholesterol-efflux regulator) over the more-abundant miR-93-5p (ρ −0.21 composition-adjusted, p=1e-10) — and
miR-33b has *lower* miRTarBase evidence than the masked arm, so it is the **share, not the evidence score**,
that recovers the canonical edge. Robust survivors are TSG/apoptosis-loaded (FOXO1, BTG2, SPRY2, CASP8/7,
SMAD5/7, ESR1, NCOA3). **Caveats:** the share is the attribution axis (⊥ realized coupling, FORMULAS §5a), so
each case rests on the *coupling* for "is repressed" and the *literature* for "which arm" — never the share
alone; and "canonical" is a noisy automated proxy, so the 75 are candidate-grade and the headline claims are
the curated, recognizable ones. **Bottom line:** abundance-sum tells you *how much* miRNA; pressure tells you
*which* miRNA — and on these genes the "which" is the textbook regulator.

## 12. ★ The progression landscape — cross-state trajectory biology (characterization, 2026-07-18)

The descriptive biological map of the GTEx→NAT→tumour landscape, read off the integrated `edge_card`
(MH-158/160) + arm trajectory (`mirna_state_class`) + genomic context. **⚠ DESCRIPTIVE** — `shift_class` per-edge
FDR labels rest on the MH-124 null (3–4× too narrow); role annotations are incomplete. Read the *structural
patterns* (cluster co-movement, convergence), not per-edge significance. Module:
`analyses/progression/landscape_characterization.py` → `output/learned/realization/landscape_*.tsv`.

> **⚠ RIGOR-AUDITED 2026-07-18 (`landscape_hardening_audit.tsv`) — this map was hardened after an initial too-clean
> pass.** ✅ threads ①②③④⑤⑨⑩⑬⑭⑯⑰ survive a proper test; ⚠ **⑥ downgraded** (mostly argmax-noise); ⛔ **⑦ retracted**
> (axiom-5 threshold artifact). Read the per-thread ✅/⚠/⛔ tags.

**① CONVERGENCE HUBS — cancer genes as multi-family repression sinks (✅ hardened).** Acquired repressive edges
converge on a who's-who of BRCA genes. ✅ **Nulled against regulator count:** acquired-edge count partly tracks a
gene's total regulators (ρ=0.61), but **PTEN is a genuine outlier — +11.8 acquired edges beyond its regulator-count
expectation, the TOP residual** (then WEE1, CDKN1A, TGFBR2, XIAP, TGFB1). *(Report the residual — the raw "28 from
20 families" is inflated because PTEN has 90 total regulators.)* Also **ZEB1** (acquired dose 2.10), **EZH2** (dose
1.32, realized ρ −0.19), CCND1, CDKN1A, STAT3, BCL2, FOXO1, TGFBR2. The canonical **miR-200 ↔ ZEB1/EZH2** axis
appears directly as coordinated acquired repression.

**② CLUSTERS MOVE AS UNITS — and it's REAL beyond co-transcription (✅ hardened).** ⚠ Same-host dose-concordance
0.9–1.0 is **near-trivial** (co-transcribed arms share a primary transcript ⇒ paired Δx-corr +0.57 by construction).
✅ **The non-trivial, hardened finding:** across-loci, same-family arms co-vary at **+0.29 vs random +0.12** (excess
+0.17) — genuine *trans*-coordination of paralogous loci (miR-200 chr1↔chr12; the miR-17~92 super-family). OncomiR
clusters ACQUIRE (miR-17~92/MIR17HG, paralogues miR-106a~363 / miR-106b~25, miR-183/96/182); the **14q32/DLK1–DIO3
imprinted locus** (MEG8/MEG9/MIR493HG, ~47 arms), **miR-30**, and **let-7/miR-100/125** LOSE.

**③ DRIVER / LOSER rosters (named).** *Acquired oncomiRs → TSGs:* the miR-200 family (patient-specific, own-frac
0.71–0.75) + miR-96/183 cluster + miR-21 → **PTEN, FOXO1 (miR-96→FOXO1 ρ −0.32), CDKN1A, TXNIP, DLC1, PTPN14**.
*Silenced TSG-miRs → oncogenes:* miR-486→CDK4, miR-451→AKT1, miR-145→CDK4/MUC1, miR-335→BIRC5, miR-15/16/195/497→
E2F3/IGF1R/RET, miR-205→ERBB2/HER2.

**④ SUBTYPE.** Subtype-specific acquisition concentrates in **Basal/TNBC** (miR-599/802/137/551a/516a/519a, all
subtype_tau>0.85 Basal) — but magnitude is small (rank-only), so a lineage flavour on the shared trajectory, not a
separate program. (Complements §5's bulk PAM50 pressure view.)

**⑤ MOST trajectory change is a pre-malignant FIELD EFFECT, not tumour-specific acquisition.** Per-arm trajectory
census (`mirna_state_class.primary_class`, built on within-state percentile-**rank** deltas dHN/dNT/dHT — the
cross-platform-robust axis, trusted over QN magnitude): of the arms that move, **field-effect 91** (58 gain + 33
loss — already changed in histologically-normal NAT) vs **truly tumour-acquired only 12** (miR-449a/b cluster,
miR-802, miR-33b); plus 467 stable. ⇒ the malignant miRNA program is largely laid down in the field.
**⚠ RANK vs QN-magnitude (trust the rank axis):** the confident acquired-gainer set is the **rank-supported 101**
(69 rank-only dHT>0.15 + 32 rank+magnitude); a further **141 are magnitude-only**, i.e. defined *solely* by the
QN'd `log2fc_tumor_vs_healthy` — the cross-platform bridge is a softer assumption than ranks, so these are treated
as the lower-confidence supporting layer, not headlined. (The paired NAT→tumour dose used in ①–④ is TCGA-only,
same-platform, so not a QN concern at all.)

**⑥ REGULATORY HANDOFFS — ⚠ DOWNGRADED (mostly argmax-noise).** Of 292 genes whose *dominant* regulator (by
abundance rank) differs healthy→tumour, only **~64 (22%) are DECISIVE** (both-state margins >5 pctile); the rest are
close-call flips — median tumour margin **2 percentile points** (a coin-flip). ⚠ **Even the flagship PTEN
miR-148a→miR-21 is a TIE in tumour** (margin 0.0). The *dose* rise of miR-21 (§3) is real; the "dominant-regulator
handoff" framing by abundance-argmax is not. Use the ~64 decisive switches only, not "333."

**⑦ THE FUNCTIONAL 2×2 (acquired dose × realized) — ⛔ RETRACTED (axiom-5 threshold artifact).** The driver/buffered
counts are entirely threshold-dependent: driver fraction **0.54→0.38→0.25→0.14→0.06** as the realization cut moves
−0.05→−0.25, 28% of genes sit within ±0.05 of the −0.1 boundary, adjacent-cut set overlap 0.64. **No stable bimodal
split — realization is a CONTINUUM.** The "455 buffered / 284 driver" headline is void. *(That acquired dose ≫
realized dose overall is the real, TESTED statement — MH-160/163, not this 2×2.)*

**⑧ SUBTYPE-SPECIFIC CLUSTER ACQUISITION.** Basal broadest (20 subtype-specific gainers); **LumA acquires C19MC**
(chr19 miR-524/525/526); **Her2 acquires miR-371~373** (miR-371a/302b); LumB miR-325/124. Distinct clusters per lineage.

**⑨ GENOMIC CONTEXT predicts trajectory direction** (paired NAT→tumour dose, same-platform — trusted). **Intergenic
arms ACQUIRE** (mean dose +0.47, dominant `acquired_realized` — the oncomiR clusters, e.g. miR-200a/429); the
**imprinted 14q32/DLK1–DIO3 locus LOSES hardest** (50 arms, mean −0.44, 68% losing); lncRNA-hosted arms lean loss
(many are the 14q32/MEG lncRNAs); sense-coding-host mixed (+0.16). Where a miRNA lives predicts gained-vs-silenced.

**⑩ WHICH PROGRAMS THE ACQUIRED PRESSURE LANDS ON.** ⚠ **A FLAT bag-of-genes MEAN is inadequate — retracted as a
real assessment.** Averaging gene `acquired_vs_nat` over a program ignores how the program is BUILT: master
regulators vs peripheral effectors, and the SIGN (repressing a pathway activator vs a repressor has opposite effects
on program output). The flat mean merely ranked IL6/JAK/STAT3, Hedgehog, OXPHOS, cholesterol, complement top — a topology-blind first
pass that MISSES THE SIGN (repressing a suppressor program is pro-tumour; repressing a proliferation engine is
anti-tumour). **The correct assessment is `geneset_architecture`'s architecture-weighted ACQUIRED pressure**
(master-regulator reverse-PageRank + signed effect-on-program-output + gene-role malignancy weighting;
`total_mal_pro_tumor_acquired`, `architecture_all_set_summary.tsv`), which gives a QUALITATIVELY DIFFERENT and
correct answer: **acquired miRNA pressure pushes toward malignancy by REPRESSING tumour-suppressor programs —
P53_PATHWAY (+18.0) and APOPTOSIS (+3.0) top the pro-tumour-acquired list** (brake-release on TSG programs), while
it **DAMAGES proliferation engines — G2M_CHECKPOINT (−11.2)**, glycolysis, angiogenesis (net anti-tumour). ⇒ the
*where the dose lands* (flat) and *what it does to program output* (signed architecture) are different questions,
and only the latter is the real program assessment. (Architecture-acquired here is GTEx-anchored `c_tumor−c_gtex`;
the topology/sign is the point, not the dose source.) A lesson in *read the doctrine before inventing* — the
flat-mean was a reinvention of a tool that already does this correctly.

**⑪ TWO INDEPENDENT WAVES (healthy anchor, RANK-based).** Decomposing each arm's healthy→tumour rank move into the
**FIELD step dHN** (healthy→NAT) and the **MALIGNANT step dNT** (NAT→tumour): **field_only 150 arms ≫ malignant_only
32**, and the two steps are **UNCORRELATED (ρ(dHN,dNT)=0.003)**, field-dominant in 68% of arms. ⇒ the pre-malignant
field wave is both larger AND a *distinct set of miRNAs* from the malignant-transition wave — not one continuous
progression. malignant-only = miR-449a/b, miR-802, miR-33b (truly tumour-specific).

**⑫ Rank-acquired gainers are mostly DE NOVO** (dHT>0.15, trusted rank axis): **83 acquired-de-novo** (silent in
healthy, rank_gtex<0.5) vs 18 amplified-from-healthy — e.g. miR-486-5p/105/155-3p/124 (de novo). ⚠ **MEASURE TENSION
(flagged, not reconciled):** miR-486-5p reads *gained* by rank (dHT +0.75) but *lost* by paired dose (mean_own_shift
−3.94, loser roster) — different sample sets (full-cohort median rank vs 103 paired) + scales (rank vs log2
abundance); a genuine divergence point between the rank and paired views, per-arm flag not blind reconciliation.

**⑬ COHORT vs PAIRED (all ~1000 tumours vs ~104 NAT, unpaired, via `arm_lfc_NAT_TUM`).** The cohort NAT→tumour dose
AGREES strongly with the 103-pair `mean_own_shift` (**Spearman ρ=0.94, sign-agreement 0.87, ZERO big
sign-divergences** over 526 arms) ⇒ **the paired view generalises — ①–④ are not a matched-subset artifact.** The
cohort is better-powered (keeps all ~900 unmatched tumours) and `coupling_tum` (~1000-tumour) is the well-powered
cohort coupling vs the n=103 paired `edge_rho_adj`. This also RESOLVES the miR-486-5p tension (⑫): both paired
(−3.94) and cohort (−4.06) agree it is LOST NAT→tumour; the "gain" was the healthy→tumour RANK (dHT +0.75) ⇒
miR-486 **rose in the field then fell in the malignant step** — real two-step biology (⑪), not a measure artifact.

**⑭ DIRECT GTEx→TUMOUR (rank dHT, bypasses the field-contaminated NAT intermediate).** The direct healthy→tumour
rank shift tracks the **FIELD step dHN at ρ=0.84** but the **malignant step dNT only at ρ=0.45** ⇒ *what's different
in cancer vs TRULY healthy is predominantly field-established* — NAT already carries most of it. 101 direct gainers,
44 losers. (Confirms ⑪/⑬ from the healthy-anchored direction; miR-486 dHT +0.75 = field dHN +0.80, malignant dNT ~0.)

**⑮ CN LAYER — copy number is a real SECONDARY contributor to acquisition, not the dominant mechanism** (⚠ CORRECTED
2026-07-18 — the initial cross-arm ρ=0.10 was a coarse-metric null: cohort-median CN is near-integer 2/3/4, so a
Spearman *correlation* is insensitive; the right test is dichotomized). **Amplified-locus arms (CN>2.5, n=116)
acquire +0.127 vs diploid (n=356) −0.044, MWU p=0.035** — a significant CN effect. And it operates for KEY oncomiR
acquirers: **miR-96-5p/182-5p/21-5p/21-3p sit at amplified loci (CN=3) with CN→expr concordance ρ≈0.19–0.27**;
miR-200c/141 are diploid but CN-concordant (0.23–0.30). But it is NOT dominant genome-wide — only **8%** of acquiring
arms have CN→expr partial ρ>0.2 (median 0.046). ⇒ acquisition is **predominantly transcriptional/epigenetic, WITH a
real copy-number-amplification contribution for a subset that includes textbook oncomiRs (miR-21/96/182)**.

**⑯ AGO/RISC LAYER — realization is NOT AGO-loading-gated.** Edge realization (`edge_rho_adj`) is flat across
`ago_loading` tertiles (−0.042 low → −0.050 mid → −0.042 high; ρ=0.08) ⇒ RISC-loading capacity is not the realization
bottleneck.

**⑰ EPIGENETIC MECHANISM (methylation Δβ gate, imprinting-aware — reuses `learned.methylation`; Δβ=tumour−normal
sees through the 14q32 imprinting that absolute β cannot).** ✅ **Positive controls validate the detector**:
miR-124/9/129/137 hyper-methylate (Δβ **+0.18**, gate 100%, normal β 0.17→tumour 0.34). ⛔ **The 14q32 coordinated
LOSS is NOT hyper-methylation silencing** (hypothesis REFUTED): the locus has high *imprinting* baseline β (0.87) and
in tumour **73% of its arms HYPO-methylate** (mean Δβ −0.12) — the opposite of silencing; the loss mechanism is NOT
promoter methylation (unexplained — imprinting-control-region / transcriptional). **Acquisition has two PARTIAL
routes**: promoter-hypomethylation (miR-21/141/200c/miR-17, Δβ −0.14 to −0.18) and CN-amplification (miR-96/182, from
⑮) — but together they explain only a **MINORITY (~35%, ~13 hypometh + ~12 CN of ~60 top acquirers)**; the majority
(39/60) have neither identifiable route (likely TF/enhancer-driven). ⇒ acquisition is mechanistically HETEROGENEOUS,
not one epigenetic switch.

**Bottom line:** the malignant miRNA trajectory is *cluster-coordinated* (oncomiR polycistrons acquired, the 14q32
TS-cluster and miR-30/let-7 lost, as co-transcribed units), *convergent* (a handful of cancer genes — PTEN above
all — absorb pressure from many independent families), *field-established and TWO-STEP* (the healthy→NAT field wave
is larger than and UNCORRELATED with the NAT→tumour malignant wave, ρ=0.003; the direct GTEx→tumour difference is
field-dominated, dHT~dHN ρ=0.84), *rewiring* (~64 DECISIVE dominant-regulator switches — not the 333 argmax flips),
and *acquisition ≫ realization* (a continuum, NOT a bimodal buffered/driver split — the 2×2 counts are retracted; dose acquired
≫ dose realized). MECHANISM: acquisition is predominantly transcriptional/epigenetic with a real SECONDARY
copy-number contribution (amplified loci acquire more, MWU p=0.035; key oncomiRs miR-21/96/182 are CN-amplified —
CORRECTED from an initial coarse-metric null), and realization is NOT AGO-gated (ρ=0.08). Acquisition splits into
promoter-HYPOmethylation (miR-21/200) + CN-amplification (miR-96/182) routes but only for a MINORITY (~35%); the
14q32 coordinated LOSS is NOT hyper-methylation (it hypo-methylates — imprinting-aware Δβ refutes the silencing
hypothesis; positive controls miR-124/9/129/137 pass at +0.18). Program effect requires ARCHITECTURE not a flat mean (acquired pressure is pro-tumour by
repressing P53/apoptosis, anti-tumour by damaging G2M — `geneset_architecture`). Paired (103) and cohort (~1000)
agree (ρ=0.94). Cross-platform GTEx: trust the rank axis (dHT/grank) + Shapley over the QN magnitude bridge.
Textbook-concordant (miR-200/ZEB1, miR-17~92, 14q32, miR-21→PTEN, C19MC) — a *description*, not a null-tested claim.

## Literature standing — one-line summary

| theme | standing |
|-------|----------|
| Targeting map (signaling/plasticity > metabolism) | concordant |
| Hub + program-specific miRNAs (miR-29→EMT, miR-34a→WNT, miR-21/155→NF-κB) | strongly concordant |
| Differentiation/hormone programs repressed by miRNA pressure | plausible, partly novel framing |
| Proliferation/immune positive couplings are confounds | expected |
| **Basal immune-pressure axis (cold-Basal embodiment)** | **novel-in-this-frame, concordant with APM** |
| **Basal TSG/apoptotic/cell-state hub (PTEN/p21/TGFBR2/VIM via miR-17~92/106b)** | **proliferation- and weighting-robust (MH-17/18/19)** |
| **IRF1 route (miR-23a-3p, miR-106b-5p) basal-specific** | **exploratory — fails proliferation adjustment within Basal; sequence prediction favors miR-130b/301 (MH-19)** |
| **TargetScan orphans (TS-only, not miRTar HE)** | **program breadth inflation; hub clusters on p21/PTEN/TGFBR2 seed families; subtype-private but idiosyncratic — exploratory (MH-23/24/25)** |
| **Curated HE hub across subtypes** | **constitutive — 3/4 miRNAs pan-subtype, Jaccard ≥0.75, zero Basal-only (MH-25)** |
| Basal tumor-suppressor/apoptosis pressure signature (PTEN/FOXO1/p21/BIM) | concordant with TNBC biology |
| 8q24 dosage coupling (miR-151a) + AGO2 amplification | concordant |

## Confidence ladder

- **Strong (multi-line / adjusted / cross-pipeline):** targeting map ranking;
  hub/program miRNAs; **basal TSG/apoptotic/cell-state repression hub
  (PTEN/p21/TGFBR2/VIM via miR-17~92/106b), robust to proliferation and to
  degree-/sequence-based weighting (MH-17/18/19)**.
- **Moderate (adjusted but single-pipeline):** differentiation-program repression;
  broad basal multi-program pressure (immune programs at the *program* level survive
  proliferation but the specific immune-gene routes are weaker).
- **Exploratory (descriptive / confounded / small-n):** the **IRF1 immune-priming
  route** (miR-23a-3p/miR-106b-5p) — fails proliferation adjustment within Basal and
  is reassigned to miR-130b/301 by sequence prediction; BCL2L11/BIM route; Normal-like
  -driven pooled negatives; AGO rate-limiting role; absolute enrichment folds.
