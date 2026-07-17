# Cross-state miRNA pressure landscape

> ⚠ **READ FIRST — the per-edge FDR class counts below rest on an UNCALIBRATED NULL (MH-123/124, added 2026-07-17).**
> The theoretical t-null used throughout this report is **3–4× too narrow**: site-free pairs — which *cannot*
> repress — have an observed SD of **0.083–0.132** against the theoretical **0.031**, and **25.0–35.3% of those
> impossible pairs pass "FDR q<0.05"**. Under an honest empirical null, the fraction of curated HE edges surviving
> per-edge drops from a nominal **31.3% to 0.0%**.
> **This does NOT mean the edges are fake or that this report is void.** A typical HE edge is ρ≈−0.15 against
> σ₀≈0.09 (z≈−1.4) — the signal is a **set-level distributional shift, not per-edge significance**. So:
> **the class STRUCTURE and the top-movers are readable; the per-edge FDR COUNTS are not.** Do not re-quote the
> `R00`/`00R`/neg-sig counts as significance claims.
> ⚠ Also: this report is **abundance-level**. It is **not** the state *channel* — that was measured and
> **CANCELLED** (τ≈0, info 0.6%; `learned/channel_state.py` was never built). Current state:
> [`../STATE_OF_PLAY.md`](../STATE_OF_PLAY.md) Axis 6.

**Scope.** GTEx-healthy breast → TCGA-NAT → TCGA-tumor, 721 miRNA arms / 5,108 miRNA→gene edges over
the MSigDB Hallmark gene universe (1,410 target genes; `mature.fa` MIMAT mapping → 2,221 expressible
`hsa-` arms in the expression matrix). Built by `mirna_state_class` (Stage 1), `mirna_comovement`
(Stage 2), `geneset_architecture` (Stage 3); outputs under
`mirna_hallmark/output/tissue_reference/{mirna_state_class,mirna_comovement,geneset_architecture,decoupling_validation}/`.
This report states the **preferred/default** readout for each question and only references the
alternative methods where they act as a sensitivity check (with the reason the default was chosen).
Formula definitions live in `FORMULAS.md`; per-module I/O in `ANALYSES_CATALOG.md`; run history in
`ANALYSIS_RUN_LEDGER.md` (this is not a changelog).

> **Refresh status (2026-06-26 — Normal-like-excluded + batch-corrected + TNRC6 co-limitation gate).**
> Updated to the settled rerun (cross-state tumor n=1060): §1 acquired-axis split; §3a/§3b/§3d
> edge-class and three-state pattern tables; §4 de-coupling end-to-end (conditioning masks, arm-/gene-side
> mechanism splits, and the symmetric gained view, all on the synced 210-edge `lost`/`nat_decoupled`
> set after re-running `mirna_comovement`→`decoupling_validation` against the final `mirna_state_class`);
> §4c gene-level net-repression; and §6a module composite. **§6 NMF / architecture factor tables are
> NOT individually refreshed:** NMF factor IDs are permutation-arbitrary (they renumber on every
> rebuild), so the per-factor cell values and the `/1096` cohort denominators in §6 reflect the
> pre-rerun fit — the stable, refreshed takeaway is *which subtype runs which program* (Basal-concentrated
> TNFA/KRAS/E2F gain, LumB immune/ALLOGRAFT), not the exact loadings.

## 1. Method spine

Every edge carries **two** measurements, and the headline class is their join:

- **Pressure (potential):** `c(m,g,X) = w(m,g) · f_expr(m;X)` — static mirTarBase/TargetScan edge
  weight × abundance multiplier (`softmax_logrpm`, cross-state comparable). Direction inherited from
  the miRNA's abundance trajectory. Edge weight `w(m,g)` uses the **S1 scorer** (`confidence_logclass`:
  assay-directness × log1p × Functional MTI cross-counts, binding-only CLIP down-weighted 0.5) plus
  ENCORI AGO-CLIP depth boost (α=0.5) on the ~25% of pairs with an ENCORI entry; see METHODS.md §2a.
- **Coupling (realized):** ESTIMATE-adjusted partial Spearman ρ between the miRNA's abundance and the
  target's expression, computed **per state** (composition + proliferation controlled). Negative ρ =
  repression realized. Rolled up per miRNA as the Stage-2 module composite.

Coupling everywhere is the **composition-adjusted partial-ρ** per state: **ESTIMATE**
immune + stromal + **epithelial metagene** + proliferation (`estimate_epi`, the default),
**plus per-state sequencing-batch dummies** (tumour/NAT analyte-plate, GTEx SMNABTCH·SMGEBTCH;
covariates, not ComBat — added since the 2026-06-26 batch refresh).
Plain `estimate` (no epithelial axis) and full `metagene` remain sensitivity audits — the
latter agrees with `estimate_epi` on ~84% of focus edges; `estimate` under-adjusts epithelial
lineage co-expression (miR-200 family).

**Acquired signal is read on two complementary axes** (both in `mirna_state_class.tsv`):
- **rank** (`dHT` = healthy→tumor percentile-rank delta on the global pressure budget): cross-platform
  robust, but **saturates** — an already-abundant arm cannot gain rank (max dHT ≈ 0.08 in the 80–95%
  bin; Spearman(dHT,log2fc)=0.53), so **81% of >4-fold-up arms are rank-invisible** (140/173 magnitude-only);
- **magnitude** (`log2fc_tumor_vs_healthy`, same QN'd TCGA scale): recovers exactly those arms.

`acquired_gainer` = rank-gain OR magnitude-gain; `acquired_axis` labels which fired. **242/721 arms
gain on ≥1 axis** (69 rank-only, 141 magnitude-only, 32 both) vs ~101 by rank alone. The default
acquired readout is **both axes together**; rank alone is retained only to show what saturates.

**GTEx collapse-repair.** GTEx's uniquely-mappable pipeline zeros ~278 canonical arms (let-7a, miR-30a,
miR-182, miR-17…) by multi-mapping collapse, which would floor their GTEx rank and inflate
`dHT`/`acquired_vs_gtex`. The GTEx **abundance** matrix is therefore repaired at source
(`impute_gtex=True`): collapsed arms take their healthy-anchor `gtex_median` (same log2(TPM+1) scale;
miTED rank-transfer for `mited_anchor`) before any contribution, so `c_gtex`/`dHT`/`acquired_vs_gtex`
are collapse-free — the let-7a/miR-30a/miR-182 collapse artifacts do not register as gainers (dHT≈0),
and `frac_gain_imputable_artifact` = **0/1,410 genes**. **Coupling stays on the raw GTEx matrix** (it
re-derives its own bundles), so partial-ρ is untouched.

**Healthy-reference coverage.** GTEx is tried first; only arms missing/collapsed in GTEx (median ≤0.05)
go to the miTED-anchor rank-transfer, then GTEx seed-family — flooring at 0 is the last resort. Of the
721 edge-universe arms: **323 direct GTEx, 234 imputed** (216 miTED-anchor + 18 seed-family),
**164 `floor0`** (no reference in GTEx *or* miTED). At the edge level (5,108 edges): 65% direct GTEx,
27% imputed, 8% floor0. floor0 arms are left un-imputed — flagged (`low_healthy_confidence` covers
**408/721 arms** = 164 floor0 + 234 imputed + blood), not imputed, because imputing a truly-absent arm
would erase a real acquired gain — and only **22** floor0 arms are flagged acquired gainers (the sole
residual false-acquired risk). The imputed baseline feeds **only the acquired/dev axis**
(`softmax_devhealthy_logrpm`), never the within-tumor z-spine or coupling, so headline coupling findings
are baseline-independent. `--no-impute-gtex` reverts for audit.

### 1a. The edge-weight prior is pan-context — breast-context re-ranking (`edge_breast_context`)

The edge weight `w(m,g)` (and the HE filter that selects which edges enter the landscape at all) is
**tissue-agnostic**: miRTarBase carries no tissue field, so an edge is "high-evidence" if its functional
repression was validated in *any* cell line / tissue (HeLa, HEK293, liver, …), and `w` is the same in a
breast cohort as anywhere else. This is a real limitation of the prior, and it is quantifiable. Mining the
**PMIDs** behind each HE edge's Functional-MTI rows (PubMed MeSH breast-neoplasm terms OR a breast
cell-line / mammary mention in title+abstract; cached under `data/external_cache/pubmed/`) shows that only
**652 / 5,293 HE edges (12.3%)** have *any* breast-specific validating paper. The rest are credited to the
landscape purely on out-of-tissue evidence.

As a **non-circular literature prior** (independent of TCGA expression, so it cannot leak into the coupling
findings), re-ranking each gene's regulators by
`breast_weight = log1p(evidence) + β·log1p(n_breast_pmids)` changes *which arm is credited as the gene's
top regulator* for **90 / 783** multi-regulator genes at β=1. A β-sweep {0.5, 1, 2} → {37, 90, 94} flips;
**37 are robust** (already flip at the most conservative β and keep the same credited arm), of which **9
are "strong"** (≥2 breast PMIDs). The strong reassignments recover textbook breast biology that the
pan-context weight had ranked below a more heavily (but non-breast) studied arm:

| gene | regulators | pan-context top arm (ev) | breast-credited arm (ev) | breast PMIDs | breast evidence |
|------|-----------:|--------------------------|--------------------------|-------------:|-----------------|
| PTEN  | 91 | miR-20a-5p (15) | **miR-21-5p** (6)  | 7 | MCF-7, MDA-MB-231, MeSH Breast/TNBC |
| ERBB2 | 22 | miR-375 (9)     | **miR-125b-5p** (6)| 3 | MCF7, SKBR3, MCF10A, MeSH Breast Neopl. |
| EGFR  | 37 | miR-218-5p (9)  | **miR-146a-5p** (6)| 3 | MDA-MB-231/435, MeSH Breast Neopl. |
| RAF1  | 8  | miR-7-5p (9)    | **miR-195-5p** (6) | 2 | MCF-7, MeSH Breast Neopl. |
| SNAI2 | 20 | miR-429 (9)     | **miR-506-3p** (6) | 2 | MDA-MB-231, MeSH Breast Neopl. |
| FOXO3 | 19 | miR-21-5p (9)   | **miR-96-5p** (6)  | 2 | MeSH Breast / Breast Neopl. |
| ESR1  | 18 | miR-9-5p (9)    | **miR-18a-5p** (6) | 2 | MeSH Breast Neopl. (ERα repressor) |
| STMN1 | 9  | miR-221-3p (9)  | **miR-101-3p** (6) | 2 | MeSH Breast Neopl. |
| BCL2  | 69 | miR-9-5p (10)   | **miR-16-5p** (9)  | 2 | MCF-7, MeSH Breast Neopl. |

The hub itself does **not** move (consistent with MH-18/MH-19, where weighting scheme reorders specific
arms per target but not the family-level hub) — the value is strictly at **per-target arm resolution**:
*which* arm to nominate for a given gene in a breast context. **Caveats:** this is a *prior*, never fed
back into the TCGA coupling; the remaining 28 robust flips rest on a single breast PMID (thin); and there
is no external breast-miRNA cohort to validate utility (METABRIC has no miRNA; CPTAC miRNA is
TCGA-reprocessed). Full per-edge scores and the β-sweep robustness table live under
`tissue_reference/edge_breast_context/`; see **MH-31** and METHODS §2a. The landscape spine still runs on
the pan-context S1 weight — breast re-ranking is an auditing/nomination overlay, not the default budget.

Breast context is **one** axis. `edge_prior_refinement` composes it with the other defensible,
non-circular axes into a single transparent nomination score, each axis a **separate visible term** so
nothing is hidden:
`refined = log1p(evidence) + w_ts·log1p(ts_weight) + w_brst·log1p(n_breast_pmids) + log(n_func/(n_func+λ·n_nonfunc))`.

| axis | what it adds | coverage over 5,300 HE edges | honest handling |
|------|--------------|------------------------------|-----------------|
| **base** | literature study count (the spine anchor) | 100% | unchanged S1 anchor |
| **A — sequence** (TargetScan Σ\|weighted context++\|, folds site-type + conservation/PCT) | binding-site strength, study-count-independent | **40%** have a predicted site | missing-site edges kept **neutral, not zeroed** — real edges can lack a conserved site (miR-21→PTEN has none) |
| **breast literature** (Functional-MTI breast PMIDs) | tissue relevance | **12%** | absolute count, not fraction (avoids 1/1-paper inflation) |
| **B — negative evidence** (Non-Functional MTI contradictions) | reliability down-weight | **25.9%** of HE edges are contradicted ≥once | ratio `n_func/(n_func+λ·n_nonfunc)` self-protects: miR-21→PTEN stays 0.86 despite 10 contradictions; only genuinely contested (1 yes/1 no) edges are demoted |
| **C — guide/passenger** (miRGeneDB `mat=1`/`star=1`; 3-tier map: exact-seq → seed+arm-end → hairpin propagation + legacy bare-name→guide) | passenger arm down-weight (+ guide = a breadth signal) | **91% of edges** classified | external curation, **non-circular** (not TCGA arm ratios); only confirmed passenger arms penalised, rest neutral; all passenger-driven flips demote a star arm to a guide (e.g. IRAK2 miR-15a-3p→miR-16-5p) |
| **D — APA / 3′UTR shortening** (Xia breast-TUMOUR ΔPDUI + APAatlas GTEx-breast + PolyASite) | **diagnostic flags only** | Xia breast-**tumour** ΔPDUI: 45 HE targets, 40 shortened (2% of edges, directional); APAatlas normal-breast 16%; PolyASite 95% | gene-level → **cancels in within-gene arm re-ranking**, flags not weights; Xia (user-supplied `XIA_APA.xlsx`) is the headline tumour-specific signal but small coverage; *site-resolved arm-discriminating version deferred* (needs TargetScan-site→genomic) |

At unit weights this reassigns the top regulator for **382/783** multi-regulator genes (high and
weight-dependent, stated plainly). Every reassignment is attributed to the axis that caused it (breast
146, negative-evidence 111, sequence 77, guide/passenger 48). Notably **a quarter of HE edges carry a
contradicting non-functional assay**, and **16% of edges target a 3′UTR that is shortened in breast
tissue** — signals the spine currently ignores. Two follow-on artefacts: a **miRGeneDB-validated
universe variant** (4709/5300 edges, high-confidence-miRNA subset) and an **exoneration** check showing
the arms miRGeneDB discards couple far weaker in TCGA-BRCA (neg-sig 10.6% vs 23.6%) — empirically
supporting the filter, with 40 couple-anyway candidates. Outputs under
`tissue_reference/edge_prior_refinement/`; see **MH-32**. Overlay, not the budget; **no external-cohort
validation**.

---

## 2. Per-gene total pressure and the acquired program

`compute_gene_pressure` returns the **total** repressive pressure on each gene summed over every
targeting miRNA, per sample; `gene_acquired_pressure.tsv` (1,410 genes) gives per-gene pressure in
each state, the acquired gain (`acquired_vs_gtex` healthy-anchored; `acquired_vs_nat` collapse-free),
the **top gaining arms**, and Hallmark membership. The per-miRNA decomposition of any one gene is the
edge table grouped by `gene`.

The imputable-artifact inflation is absent (`frac_gain_imputable_artifact = 0/1,410`) and
`acquired_vs_gtex_corrected ≈ acquired_vs_gtex` (the post-hoc correction is redundant once the matrix
is repaired, kept only as a diagnostic). **1,242/1,410 genes are fully clean** (`frac_gain_from_floor0
= 0`); the residual **70 `floor0=1.0` genes** (CANX, DHRS3, ARPP19, ARF1, RAB14… — gains carried only
by arms with no healthy reference, e.g. miR-320a/miR-2355/miR-3074) are a genuine reference-coverage
gap (low/absent in healthy breast), not a collapse artifact. The top acquired genes are a coherent
**miR-21-5p program**: PCBP1 (+7.56), PYCR1 (+7.28, miR-2355-5p), PTPN14 (+7.21), IL12A (+7.00),
TPM1 (+7.00), WWP1 (+6.33), EIF4A2/DNM1L/STRN/DOCK4 (+6.22), CCL20 (+5.85), PLAT (+4.98). The
collapse-free NAT-anchored readout is a **same-platform cross-check of the acquired-pressure claim**,
not a gene-for-gene repeat of the GTEx list: `acquired_vs_nat = pressure_tumor - pressure_nat` ranks
the final NAT→tumor step and brings up HUS1 (+2.79), ENG (+2.30), YWHAG (+1.91), PTEN (+1.89),
INHBA (+1.73), GSN (+1.44). Thus GTEx→tumor highlights the large miR-21 healthy-baseline jump, while
NAT→tumor shows a complementary malignant-step ramp. The **magnitude_only** arms the rank was blind to
are the BC oncomiR surge (miR-183/182/375/93/141, +9 to +13 log2; see §1).

---

## 3. Edge classes — biology and novelty

### 3a. `acquired_realized` (640 edges, 129 arms) — pressure rises **and** coupling turns repressive
The newly-switched-on repressive edges (no significant normal repression → significant tumor
repression; coupling pattern `00R`/`P0R`) — the genuinely-acquired-and-coupled subset. Biological
standouts:

| edge | edge share (tumor) | ρ gtex→tumor | note |
|---|---|---|---|
| miR-590-5p → RBPJ | 0.86 | +0.02 → −0.11 | acquired Notch-effector repression |
| **miR-744-3p → HLA-G** | 0.76 | +0.01 → −0.10 | **acquired repression of a non-classical MHC immune-tolerance ligand** — novel |
| miR-582-3p → CXCL14 | 0.67 | +0.05 → −0.12 | acquired chemokine repression |
| miR-301a-3p → BTG1 | 0.66 | −0.11 → −0.20 | strengthening repression of an anti-proliferative TS |
| miR-129-2-3p → CCP110, TIMELESS | 0.42 | ~0 → −0.15/−0.14 | acquired centriole / circadian-replication |
| miR-21-3p → RBPMS | 0.42 | −0.02 → −0.21 | miR-21 *passenger* −3p arm; RBPMS = TS splicing factor |
| miR-744-3p → TGFB1 | 0.24 | +0.01 → −0.24 | acquired TGFβ-ligand repression |

**Novelty.** The standout is **miR-744-3p acquiring coherent repression of both HLA-G and TGFB1** — an
immune-axis finding. miR-21-3p→RBPMS and miR-582-3p→CXCL14 are under-described (the *passenger* −3p arm
of miR-21 acquiring its own realized target is notable). The well-known miR-21 PTEN/PDCD4 axis appears
as **pressure** but not as a top *realized* edge — the realized edges land on less-canonical genes, the
platform's distinctive contribution.

### 3b. Released brakes — repression present in healthy, gone in tumor (`R00` coupling pattern)

The released-brake biology is read off the **`R00` coupling pattern** (FDR-significant repression in
GTEx, gone in both NAT and tumor): **189 edges / 109 arms**. Note the joint `lost` class is the
**stricter subset** (209 edges / 89 arms) where the arm's abundance is *also* flat; where the arm
additionally **gained abundance**, the same brake-release edge is joint-classed `acquired_unrealized`
(pressure up, coupling never realized) — so the example brakes below span both joint classes while
sharing the `R00` pattern.

| edge | ρ gtex → tumor | pattern (joint class) | literature |
|---|---|---|---|
| **miR-193a-5p → ERBB2** | −0.28 → +0.07 | R00 (acquired_unrealized) | miR-193a a reported ERBB2/HER2-axis TS — concordant |
| miR-331-3p → ERBB2 | −0.19 → −0.09 | R00 (acquired_unrealized) | second ERBB2 brake |
| miR-224-5p → CXCR4 | −0.25 → +0.02 | R00 (lost) | CXCR4 metastatic homing; brake loss is pro-metastatic |
| miR-22-3p → MMP14, miR-708-5p → MMP2 | −0.30 → +0.19 / −0.31 → +0.03 | R0P / R00 (acquired_unrealized) | MMP de-repression → invasion; both documented BC TS |
| miR-326 / miR-23b-3p → NOTCH2 | −0.33/−0.30 → ~0 | R00 (acquired_unrealized / lost) | Notch de-repression |
| miR-96-5p → FHL1 | −0.16 → −0.07 | R00 (acquired_unrealized) | (miR-96→PTEN is `000` — below FDR in GTEx) |
| miR-363-3p → PDPN, BCL2L11(BIM) | −0.36/−0.28 → ~0 | R00 (lost) | apoptosis/podoplanin brakes |

**Reading.** The theme — *healthy breast actively holds oncogenic/invasion targets (ERBB2, CXCR4,
MMP2/14, NOTCH2) under miRNA repression and tumor releases the brake* — is biologically coherent and,
being coupling-derived, immune to the GTEx-collapse caveat. **miR-96-5p→PTEN** is below FDR
significance in GTEx (ρ≈−0.08, pattern `000`), so it is not a brake; miR-96-5p→FHL1 is R00.

**Novelty — be precise.** Two separable claims: (1) the miRNA→target *pairs* are individually
documented (miR-193a/ERBB2, miR-22/708/MMP, miR-326/23b/NOTCH2, miR-363/BIM) and these miRNAs are
reported down in BC — the *direction* is textbook (low novelty); (2) the **state-resolved loss of
realized coupling** (ρ ≈ −0.35 in healthy → ≈ 0 in tumor, confounder-adjusted, on a healthy→tumor
trajectory) is our contribution (strength **P/H**) — the literature asserts the pair and the miRNA's
down-regulation but rarely measures the *brake being engaged in normal tissue and released in tumor*.
Strongest exactly for the thinly-reported pairs (miR-708→MMP2, miR-363→BIM). See §4 for whether these
breaks are genuine loss of repression.

### 3c. `nat_decoupled` (1 edge) — NAT-only coupling sign-flip
The strict `nat_decoupled` joint class contains a **single** edge, miR-183-5p.1→PTPA (pattern `0P0`:
≈0 in GTEx, positive/co-expressed in NAT, ≈0 in tumor). No multi-edge NAT-only sign-flip set survives
the strict FDR sign-flip-vs-both-states criterion, so the NAT-only sign-flip phenomenon — the diagnostic
the three-state separation was built to surface — is essentially **absent** as an FDR-significant class;
the few directional NAT-only oddities survive only in `coupling_dir_pattern`/`nat_underpowered_repressive`
(§3d), consistent with NAT being underpowered (n≈104) rather than carrying a robust decoupling signal.

### 3d. Three-state coupling trajectories: where the switch occurs
Use this as the explicit NAT-aware readout behind §3a–c. A state is called **R** when adjusted coupling
is significantly repressive (`ρ < -0.10`, FDR `q < 0.05`), **P** when significantly positive/co-expressed
(`ρ > +0.10`, FDR `q < 0.05`), and **0** otherwise. The three letters are **GTEx / NAT / tumor**. This is
a threshold-crossing trajectory, not a formal test that two independent correlations differ.

**NAT is underpowered (n≈104).** Report two layers in `mirna_gene_edge_class.tsv`:
- `coupling_sig_pattern` — strict FDR letters above (`P0R`, `0RR`, …);
- `coupling_dir_pattern` — directional only (`r`/`p`/`0` at |ρ|>0.10, no FDR);
- `nat_underpowered_repressive` — NAT ρ<−0.10 but NAT q≥0.05 while tumor is significantly repressive
  (keeps BNIP3L-style cases visible: dir `PrR`, sig `P0R`).

**Joint classes (pressure × coupling):** arm acquired = `arm_acquired_gainer` (`dHT>0.15` OR log2fc>2.0),
not `dHT` alone. **`field_established_realized`** = arm acquired + `0RR` + **`gtex_coupling_scored`**
(GTEx partial ρ scored on the MIMAT-mapped TCGA-arm matrix; unscored healthy → `?` in
`coupling_sig_pattern`, not `0`). Default coupling covariates: **`estimate_epi`** (ESTIMATE + epithelial sig).

**Repressive coupling gained** (GTEx MIMAT coupling + `gtex_coupling_scored` gate):

| pattern | interpretation | all edges | miRNA arms | joint class | examples |
|---|---|---:|---:|---|---|
| `00R` | tumor-step only: no significant normal repression, significant tumor repression | 798 | 195 | `acquired_realized` | miR-590→RBPJ, miR-744→HLA-G/TGFB1 |
| `0RR` | NAT-established and retained in tumor | 7 | 6 | `field_established_realized` (if arm acquired + GTEx scored) | miR-21-5p→TGFBR3, miR-375→AHR, miR-25→PTEN |
| `P0R` / `PPR` | normal co-expression replaced by tumor repression | 28 / 0 | 23 / 0 | `acquired_realized` | miR-30d→BNIP3L, miR-324→SOCS3, miR-106b-3p→CCND1 |
| `PrR` (dir) | NAT directionally repressive but FDR-fails | see `nat_underpowered_repressive` | — | often `acquired_realized` + `nat_underpowered` | BNIP3L-style cases |

**Repressive coupling lost:**

| pattern | interpretation | all edges | miRNA arms | examples |
|---|---|---:|---:|---|
| `R00` | healthy-only brake: significant in GTEx, already gone in NAT and tumor | 189 | 109 | miR-363→PDPN, miR-224→CXCR4, miR-96→FHL1, miR-628-3p→TP53 |
| `RR0` | tumor-step loss: retained in NAT, lost only in tumor | 7 | 6 | miR-206→SFRP1, miR-133b→FGFR1, miR-485-3p→NTRK3 |
| `0R0` | NAT-only brake: not seen in GTEx, present in NAT, gone in tumor | 14 | 14 | miR-499a-5p→PDCD4, miR-26b-5p→FH, miR-99a-5p→FGFR3 |
| `R0R` | **not lost** — constitutive; GTEx+tumor sig R, NAT FDR dip | 114 | 68 | miR-887-3p→DNMT1, miR-362-5p→CYLD, miR-199b-3p→HK2 (class `constitutive`, not `lost`) |

Edge-level pressure sensitivity (not a hard joint-class gate): `edge_pressure_gain_gtex` /
`edge_pressure_gain_nat` = per-edge `c_tumor > c_gtex` / `c_nat`.

---

## 4. De-coupling — is it real loss of repression, and what actually breaks?

A broken coupling is a within-state *rank correlation* lost, so it is **level-independent by
construction** — we never use raw expression direction (that re-introduces the confound the partial-ρ
design removes). De-coupling is validated by **conditioning and platform-matching**, in `decoupling_validation.tsv`
(**210** `lost`/`nat_decoupled` edges; the strict `lost` joint class is the abundance-flat subset of
the R00 brake-release pattern, §3b).

### 4a. Confidence and what is ruled out as the driver
- **Platform-matched confidence (`nat_anchored_confident`, 14/210).** NAT and tumor are the same
  platform (TCGA-seq) and both ESTIMATE-adjusted, so a partial-ρ significantly negative in NAT and ≈0
  in tumor is not a cross-platform/composition artifact: miR-206→SFRP1 (−0.54→+0.12), miR-499a-5p→PDCD4
  (−0.51→−0.02), miR-26b-5p→FH (−0.40→0.00), miR-133a-3p.2→ANXA2 (−0.38→−0.02), let-7b-5p→CPEB3
  (−0.35→+0.27), miR-543→PTK2, miR-99a-5p→FGFR3, miR-1-3p→TKT. NAT is *conservative* (field
  cancerization + small biased n≈100 + tumor-epithelium admixture all mask early loss), so 14 is a
  **lower bound**, not the full set.
- **Conditioning (all three layers, fixed pipeline).** Re-fitting each lost edge's tumor partial-ρ
  with one extra covariate and flagging re-emergence (Δ≥0.10 and conditioned ρ<−0.10):
  - **Copy number** (target ASCAT3): **1/210 masked** — de-coupling is essentially **not** a target-amplification artifact.
  - **Promoter methylation** (per-sample β over decoupled-target Hallmark genes,
    `hallmark_gene_methylation.tsv.gz`): **182 tested, 2 masked** — effectively a
    real null; de-coupling not driven by the target's own promoter methylation.
  - **Chromatin accessibility** — *imputed cohort-wide* (≈subtype-imputed, dilutes signal): 178 tested,
    **1 masked**; *strict directly-assayed* (ATAC-measured tumors, no imputation): 178 tested,
    **7 `atac_direct_masked`** — miR-7-5p→PTK2, miR-652-3p→FGFR1, miR-335-5p→IGF1R, miR-205-5p→EGLN2,
    miR-144-3p→TLR2, miR-381-3p→STC2, miR-140-5p→DCX (growth-factor / hypoxia loci) — a small
    accessibility-linked minority the imputation washes out.

  Net: de-coupling is **not driven by copy number (1/210) or promoter methylation (2/182)**; a small
  ATAC-direct minority at growth-factor loci (7/178) is accessibility-linked. The dominant mechanism
  is therefore post-transcriptional/competitive (3′UTR APA/shortening, ceRNA, RBP competition, AGO/RISC),
  with a small chromatin contribution at specific growth-factor loci.

### 4b. What breaks — arm-side mechanism and gene-side scope
The same break reads two ways, and we report both:

- **Arm-side (`decoupling_mechanism`, from the Stage-2 module composite):** **78 `arm_wide_failure`**
  (the arm's *whole* target module also decouples — the arm broadly loses realized control),
  **33 `target_specific_escape`** (module stays repressed, one edge breaks — the target's variance is
  taken over by a tumor-specific driver while the miRNA's pressure *share* stays flat), 16 intermediate,
  83 unknown.
- **Gene-side (`escape_scope`, from gene-level realized repression, §4c):** **46 `edge_local_escape`**
  (gene still net-repressed by its other regulators), **63 `gene_wide_derepression`** (the gene lost net
  miRNA control; largely single-/few-regulator genes where the lost edge *was* the brake),
  **101 `gene_not_repressed`** (gene never under aggregate realized repression — coupling modelled but
  not gene-decisive). So most lost couplings sit on genes miRNAs never bulk-controlled; of the rest,
  gene-wide de-repression outnumbers local escape and concentrates on lightly-regulated genes.

The **symmetric (gained) view** uses the same machinery (`gene_repression_scope`, all classes): of the
640 `acquired_realized` edges, **143 `gene_wide_acquired_repression`** (the gene newly comes under net
miRNA repression and this edge is part of it — miR-590-5p→RBPJ, miR-582-3p→CXCL14, miR-26b-3p→POLR1H),
196 `reinforces_constitutive`, 301 `edge_specific_acquired`. So a gained coupling can be read as
gene-wide acquired repression vs an edge-specific acquired hit, exactly mirroring the loss side.

### 4c. Gene-side share and gene-level realized repression
Two gene-centric measures complement the arm-centric ones:
- **Gene-side share (`gene_share_{gtex,nat,tumor}` = c/Σ_regulators c).** The complement of the
  arm-side `edge_share` (the arm's fraction of *its own* budget on the gene): `gene_share` is this
  arm's fraction of the **gene's total incoming pressure**, with `gene_share_rank_tumor` (1 = dominant
  regulator) and `gene_share_delta_tumor_gtex` (grip change). They answer opposite questions — e.g.
  **miR-9-5p is PTEN's #1 regulator** (gene_share 0.080, rising from 0.048 → tightening grip) yet PTEN
  is a trivial slice of miR-9's portfolio (edge_share 0.006).
- **Gene-level realized net-repression (`gene_corepression.tsv`).** Gene-side mirror of the arm module
  composite: per gene, the ESTIMATE-adjusted partial-ρ of its **own expression** vs its **aggregate
  incoming pressure** (summed over all regulators), per state. **282/1,424 genes net-repressed in
  tumor** — led by luminal TFs **GATA3, ESR1**, GATA1, and **AURKB/CCNA2** (representative
  exemplars; per-gene ρ re-derive each rebuild). Trajectory classes: 115 `constitutive_repressed`,
  **258 `lost_repression`** (gene-wide de-repression), 167 `gained_repression`, 867 `never_repressed`
  (+17 `na`). These classes feed `escape_scope`/`gene_repression_scope` above.
- **Co-regulator context.** The informative co-regulator number is **13/210** where the focal arm
  is top-quartile regulator share on the gene yet still lost coupling (the higher-confidence escapes).
  The "112/210 have a coupled competitor" figure is ≈ the base rate (many genes have ≥1 coupled arm)
  and is not a signal by itself — a heavily-contested gene like PTEN keeps a coupled regulator while
  individual arms decouple, explaining why no single arm "controls" it.

---

## 5. Subtype specificity

**Tau** indexes how concentrated an arm's acquired pressure is across the four PAM50 subtypes
(0 = pan, 1 = single-subtype) — it is the *confidence* of the specificity call. Counts
(`subtype_specific`, τ≥0.5; 250 arms total): Basal **141**, LumB 41, Her2 34, LumA 34.

τ is a Yanai index, so it saturates at ≈1 for any arm present in **one** subtype regardless of
magnitude — the very-highest-τ arms (miR-2467-3p, miR-302d-3p) carry near-zero pressure. **So the
default ranking is τ≥0.5 then by share in the dominant subtype**, which also reveals a structural fact:
subtype-specific arms are intrinsically **low-abundance specialists** (dominant share ≤0.005), because
the high-abundance pressure backbone is pan-subtype (the static-c NMF first factor is subtype-invariant,
§6).

This makes the two subtype resolutions answer different questions — only the first is abundance-limited:
- **Per-arm τ specificity is a low-magnitude story** — "miR-X is Basal-specific" is real but small.
- **NMF subtype *programs* are the opposite — strongly, significantly subtype-divergent**, driven by
  *abundant* arms whose program *activity* differs by subtype (within-tumor T3 hypoxia/miR-17~92 and
  gene-pressure P4 Basal TNFA/NFKB both reach Kruskal–Wallis q≈0). Subtype structure is **not**
  negligible — it lives in *how strongly each subtype runs the shared abundant programs*, a
  concentration, not in low-abundance exclusive arms.

**Top subtype-specific arms (τ≥0.5, by dominant-subtype share):**

| subtype | arms (share, τ) | note |
|---|---|---|
| **Basal** | miR-18a-5p (.004, .68), miR-19a-3p (.004, .51), miR-135b-5p (.003, .79), miR-106a-5p (.003, .52), miR-31-5p (.003, .60) | the **miR-17~92 / miR-106a cluster arms** plus miR-135b/31 — reported basal/TNBC, concordant |
| **Her2** | miR-184 (.005, .59), miR-224-5p (.003, .51), miR-187-3p (.003, .54), miR-452-5p (.001, .58) | miR-184 & miR-187 Her2-acquired — candidate-novel |
| **LumB** | miR-124-3p.1 (.005, .59), miR-1251-5p (.0005, .57), miR-489-3p (.0002, .75), miR-34b-3p | miR-124 LumB-localized acquired gain |
| **LumA** | miR-153-3p/5p (.0006/.0004, .59/.67), miR-520c/524/525-3p (C19MC 19q13.42) | the **C19MC cluster is LumA-concentrated** — locus-level signal |

Also locus-coherent: the **DLK1-DIO3 14q32 cluster loss is Basal-concentrated** (miR-543/329/487a/376a).
Edge-level (`edge_subtype_tau`): the sharpest subtype-localized edges (τ≈0.98) are
**miR-138-5p → CDK6 / CCND1, Basal** (basal-restricted cell-cycle repression), miR-302d-3p → CXCL8/
MAP3K1 (Basal), miR-490-3p → TGFBR1 (Her2).

---

## 6. Co-movement and pressure programs (Stage 2)

### 6a. Target-module co-repression (the composite)
For each miRNA we build **one composite target signal** per sample = the mean within-state z-score
across its target genes (`z(g,s)=(x−mean_s)/sd_s`; we z *before* averaging because raw target
expression spans orders of magnitude and an unstandardised mean would track the single
highest-expressed gene), then take **one ESTIMATE-adjusted partial Spearman** of the arm's abundance
vs that composite, per state. So it is the **Spearman of the mean-z** — one ρ per miRNA per state,
asking whether the arm tracks the *coordinated* down-shift of its *whole* module (which single-edge ρ
cannot see); `Δ = ρ_tumor − ρ_gtex` is acquired module repression.

The strongest **confidence-independent** signal in the run is acquired module co-repression (composite
ρ strongly negative in tumor, q<0.05; 103 arms): miR-18a-5p (19 targets, −0.51), miR-106b-5p (29t, −0.51),
miR-17-5p (48t, −0.42), miR-301b-3p (−0.40), miR-19a-3p (27t, −0.39), miR-15b-3p (−0.39), miR-130b-3p (19t, −0.38),
**miR-93-5p (18t, −0.38)**, miR-301a-3p (−0.37). The module-level realized
onset for the miR-17~92 / miR-106b~25 clusters and miR-130b/301 is a distinctive platform output.

*Defaults & alternatives.* An **edge-weighted** composite (`rho_composite_w`, weight = evidence_score)
correlates r=0.97–0.99 with the equal-weight one (targets carry near-even evidence mass), so the
equal-weight composite is the default and weighting is shipped as a robustness column. A
**role-stratified** composite (`rho_composite_{tsg,onco}`, same ρ over the arm's TSG vs oncogene
targets) is added as an *interpretive overlay, not a reweight* — folding role into one number would
conflate "is it repressed" with "is it pro/anti-tumor" (§7). It is informative: validated oncomiRs
concentrate coordinated repression on **suppressor** targets — miR-93-5p TSG −0.37 (onco below
coverage), miR-106b-5p −0.49 vs −0.24, miR-21-5p −0.25 vs −0.07 (coverage partial: 33 arms reach the
TSG stratum, 71 the oncogene stratum). TF status is not added — it is captured upstream by the architecture
master-regulator centrality (§7).

### 6b. NMF pressure programs
The default abundance pressure (`softmax_logrpm`) drives the NMFs; the dev/healthy-anchored mode enters
only the acquired axis, so the NMFs are mechanically unaffected by the GTEx repair (which edits only
the GTEx state matrix), and the imputation never double-counts the log2(RPM+1) magnitude the pressure
formula already encodes.

> **NMF factor IDs are permutation-arbitrary.** The factor indices (T3, P4, S7, D1, …) **renumber on
> every rebuild** — the stable handle is the content-derived `factor_label` (e.g. `Basal:TNFA_SIGNALING_VIA_NFKB`).
> Read the tables for *which subtype runs which program*, not for the specific factor index. The stable
> structure is the Basal-concentrated TNFA/EMT/E2F gain programs (KW q ≈ 1e-79 to 1e-140), the LumB
> immune/ALLOGRAFT and LumA EMT/complement factors, the subtype-invariant abundance backbone, and the
> signed-NMF `gain_mass_frac` ≈ 0.97–0.99 (no organised pressure-relief program).

- **Aggregate (`c_tumor`) NMF:** top tumor-amplified program **F1 = miR-21/210/let-7e/155/29a**
  (hypoxia–inflammation), F3 = let-7a/let-7b/miR-582/196a/30b.
- **Within-tumor per-sample NMF** (programs co-varying across patients, PAM50-tested) — the per-PAM50
  activity fraction makes the subtype loading concrete and shows it is **graded, not exclusive**:

| factor | program (top arms) | LumA | LumB | Her2 | **Basal** | dominant (KW q) |
|---|---|---|---|---|---|---|
| **T3** | miR-9 / miR-17~92 / miR-210 hypoxia-oncomiR | .10 | .18 | .25 | **.48** | Basal (q≈0) |
| T8 | immune miR-142/155/146a/223 | .24 | .19 | **.29** | .28 | Her2 |
| T5 | miR-375/196a/149 luminal | .23 | **.32** | .30 | .16 | LumB |
| T2 | let-7/miR-30/miR-143 | **.34** | .22 | .22 | .23 | LumA |

  (The abundant let-7/miR-30 backbone factor is the one that is *flat* across subtypes, KW q=0.29.)
- **Per-subtype abundance NMF** (NMF run *inside* each PAM50 group; `nmf_subtype_{static,persample}.tsv`).
  The **static** first factor (C1) is near-identical across all four subtypes — the abundance backbone
  is subtype-invariant (an informative negative) — while the **per-sample** secondary factors (S2/S3)
  are genuinely divergent:

| subtype | static C1 backbone (invariant) | divergent per-sample secondary |
|---|---|---|
| LumA | miR-21/210/let-7e/155 | S2 = miR-375/103a/182/93 (luminal-proliferative) |
| LumB | miR-21/210/let-7e/155 | S3 = miR-205/196a/20a/146b |
| Her2 | miR-21/210/155/let-7e | **S2 = immune miR-224/452/150/155/146a/223**; S3 = miR-200/375/149 |
| Basal | miR-21/210/155/miR-9 | **S2 = full miR-17~92 (miR-17/19a/19b/20a/106a)** |

  So the *target-sharing* backbone does not differ by subtype; the subtype signal lives entirely in the
  per-sample secondary decomposition (Her2 immune split, Basal miR-17~92).
- **Gene-pressure per-sample NMF** (gene×tumor-sample total-pressure) — strongly PAM50-divergent gene
  programs (`nmf_gene_pressure_factor_subtype.tsv`):

| factor | dominant Hallmark | LumA | LumB | Her2 | **Basal** | dominant (KW q) | top driver arms |
|---|---|---|---|---|---|---|---|
| **P4** | TNFA_SIGNALING_VIA_NFKB | .06 | .16 | .24 | **.54** | Basal (5e-139) | miR-9/92a/378a/660 |
| P2 | COMPLEMENT | **.35** | .31 | .25 | .09 | LumA (2e-63) | miR-10a/10b/143/let-7a |
| P3 | ALLOGRAFT_REJECTION | .28 | **.36** | .30 | .06 | LumB (6e-62) | miR-375/182/342/200a |
| P9 | EMT | **.38** | .27 | .09 | .25 | LumA (1e-39) | miR-30a/30c/let-7b |
| P10 | TNFA_SIGNALING_VIA_NFKB | .19 | .19 | .27 | **.35** | Basal (3e-22) | miR-30e/150/142/126 |

  Also run **within each subtype** (`nmf_gene_pressure_within_subtype`); leading within-subtype program:
  Basal HYPOXIA, LumB G2M_CHECKPOINT, Her2 APOPTOSIS, LumA ESTROGEN_RESPONSE_EARLY.

  *How signatures distribute across samples* (`nmf_sample_signatures.py` →
  `nmf_sample_signatures/summary_all.tsv` + per-analysis heatmaps under
  `…/mirna_comovement/nmf_sample_signatures/figures/`). For every NMF that has a
  **sample axis**, the module refits the decomposition on its appropriate cohort, captures the
  factor×sample activation matrix H, column-normalises to a per-tumour signature composition, and
  scores dominance (top share ≥50% = "dominant"). Cohort-wide heatmaps order samples by PAM50 block
  then dominant signature (colour bar = subtype). **Not included:** aggregate `c_tumor` NMF
  (`nmf_programs`), share NMF (`nmf_share`), and per-subtype *static* arm NMF (miRNA×gene cells, no
  sample axis).

  **Within-subtype gene-pressure** (same refit as before; k=12 inside each PAM50):

| subtype | n | % dominant | median top-share | median eff. #sigs | leading signature |
|---|--:|--:|--:|--:|---|
| LumA | 512 | 2% | 0.32 | 6.6 | ESTROGEN_RESPONSE_EARLY |
| LumB | 258 | 14% | 0.34 | 6.3 | G2M_CHECKPOINT |
| Her2 | 93 | 12% | 0.33 | 5.9 | APOPTOSIS |
| Basal | 195 | 8% | 0.32 | 6.6 | HYPOXIA |

  **Cohort-wide (all ~1096 tumors)** — mixing is even flatter than within-subtype; subtype blocks
  visible in heatmaps but no clean global sub-subtypes:

| NMF type | n | % dominant | median top-share | median eff. #sigs | most-common signature |
|---|---:|--:|--:|--:|---|
| Gene-pressure (abundance) | 1096 | 4% | 0.32 | 6.5 | ESTROGEN_RESPONSE_EARLY (517) |
| Gene-pressure (acquired gain) | 1096 | 6% | 0.32 | 6.2 | ALLOGRAFT_REJECTION (269) |
| Gene-pressure (signed ±) | 1096 | **45%** | 0.48 | 4.6 | S3 RASGRP1/MSH2/CCL20 [flat `Her2:KRAS`] (725) |
| miRNA within-tumor (arm abundance) | 1096 | **0%** | 0.29 | 7.1 | miR-196a/199a/205 (486) |

  The cohort-wide gene-pressure view confirms the PAM50-divergent factors from the KW table (Basal
  TNFA, LumB immune) are **population enrichments**, not discrete patient classes — only ~4% of tumours
  are dominated by one *abundance* program. The within-tumor miRNA NMF is the flattest (0% dominant,
  median ~7 active programs): arm co-variation is almost always a blend across patients. The one
  exception is the **signed (gain+loss) gene-pressure NMF — 45% dominant**: its decomposition is
  dominated by a single shared acquired-gain backbone factor (top signature in 725/1096 tumours), so
  "dominance" there reflects that one universal gain program, not discrete patient subtypes.

  **Per-subtype acquired gene-pressure** (dev-anchored gain, k=12): LumA 6% / Her2 12% / LumB 19%
  dominant, but **Basal 68% dominant** — within Basal the acquired-gain decomposition is carried by one
  shared backbone gain factor (the same single-program-dominance seen in the cohort-wide signed NMF),
  so Basal's acquired pressure is far less heterogeneous than its abundance pressure.

  **Per-subtype miRNA within-subtype NMF** (k=8, arm×sample inside each PAM50 — the divergent
  secondary factors from the static/per-sample table above): **more discrete than gene-pressure** in
  Basal and LumA:

| subtype | n | % dominant | median top-share | median eff. #sigs | leading signature |
|---|--:|--:|--:|--:|---|
| LumA | 512 | **16%** | 0.39 | 5.2 | S3 miR-150/142/223 immune-pan (342) |
| LumB | 258 | 12% | 0.38 | 5.1 | S3 miR-375/99b/9 luminal (122) |
| Her2 | 93 | 1% | 0.36 | 5.4 | S3 miR-628/142 immune-pan (50) |
| Basal | 195 | **26%** | 0.42 | 4.9 | S4 miR-150/142 immune-pan (141) |

  Basal is the one subtype where a quarter of tumours are single-arm-program-dominant (S4 =
  miR-150/142/223 immune module) — closer to a sub-subtype split on the **miRNA heterogeneity** axis,
  while the **abundance gene-pressure** axis within Basal stays graded (8% dominant — though its
  *acquired-gain* gene-pressure does not, 68%, §ACQ above). Read as: subtype-specific miRNA co-variation
  modules can partition a minority of patients, but abundance-driven gene programs co-occur continuously
  inside a PAM50 class.

  For the **abundance** decompositions, only **4–26%** of tumours are dominated by a single signature
  (4% cohort gene-pressure → 26% Basal miRNA-within), median ~5–7 active signatures.

  **Divergent baseline, convergent acquisition (the key dominance contrast).** The abundance and gain
  decompositions answer two *different* questions and give opposite answers:
  - The **abundance** pressure NMF (raw `c_tumor` / arm abundance) factorises the *steady-state*
    pressure each tumour runs. Here tumours are a **continuous blend** — 4–26% dominant, top signature
    only ~⅓ of the mass, ~5–7 active programs per tumour. No discrete miRNA-pressure sub-subtypes exist;
    PAM50-divergent factors (Basal TNFA, LumB immune) are **population enrichments**, not patient classes.
  - The **gain** NMF (acquired = `c_tumor − c_gtex`; signed = ±-channel gain/loss) factorises the
    *increment the tumour added over healthy*. Here a **single shared backbone factor dominates** —
    45% of tumours cohort-wide (signed), 68% within Basal (acquired). That backbone is **subtype-flat**
    — the dominant factor is **S3** (mean H-activation ≈2× any other factor; ≈0.27/0.27/0.30/0.15
    across LumA/LumB/Her2/Basal, so its `Her2:KRAS_SIGNALING_UP` label is a thin-margin argmax over a
    near-flat profile, **not** real Her2 specificity; cf. the divergent secondaries that carry the KW
    signal). It is a broad **miR-21/26a/let-7-driven** gain program led by RASGRP1/MSH2/CCL20/VHL/HPGD.
    The §2 acquired oncomiR genes (LIFR/ABCE1/IGFBP5; IL12A/PCBP1/TPM1) sit in *distinct minor* factors
    (S9 `Basal:MYC`; S2/S8 `ESTROGEN`), **not** in the dominant backbone — a correction from the prior
    report version, where a stale standalone sample-signature fit (since unified with the comovement
    fit — see the §6b note below) mislabelled the backbone as the LIFR/ABCE1/IGFBP5 factor.

  **What it means.** *Which* miRNA programs a tumour runs at steady state is heterogeneous and continuous
  (every tumour mixes many, in subtype-skewed proportions). But *what a tumour acquires relative to
  healthy breast* is largely **one common gain program shared across tumours and subtypes** — breast
  tumours converge on the same acquired-pressure backbone rather than each acquiring a private one. So
  "dominance" on the gain axis is **not** evidence of discrete patient subtypes (the opposite of what
  46–68% would suggest at face value); it is the signature of a **universal, co-directional acquisition**.
  Two caveats keep this honest: (i) it is partly methodological — the gain matrix is overwhelmingly
  one-signed (`gain_mass_frac` 0.97–0.99: almost everything moves *up* vs healthy), so NMF's leading
  factor necessarily captures that shared "everything-up" direction; (ii) the **subtype-divergent**
  acquired structure still exists, but it lives in the *secondary* factors (Basal:KRAS, LumB:ALLOGRAFT,
  §ACQ table), not in the dominant backbone.

  Heatmaps echo this: within each subtype the **abundance** programs are a continuous mixture with only a
  minority "leading edge" (within-subtype leaders: Basal HYPOXIA, LumB G2M, Her2 APOPTOSIS, LumA
  ESTROGEN_RESPONSE_EARLY), behaving like **continuous axes that co-occur in varying proportions**, not
  discrete sub-subtypes — only the tail of each program is private to a tumour subset, and the
  miRNA-within-subtype NMF is modestly more discrete for Basal/LumA on the arm axis.
- **Acquired (dev-anchored, clipped ≥0) gene-pressure NMF** — which gene programs *co-gain* pressure vs
  healthy; even more Basal/LumB-concentrated than the abundance view (`nmf_gene_pressure_acquired_factor_subtype.tsv`):

(keyed on stable `factor_label` per §6b note; F-index in parentheses)

| factor_label | LumA | LumB | Her2 | **Basal** | dominant (KW q) | top driver arms |
|---|---|---|---|---|---|---|
| **Basal:KRAS_SIGNALING_UP** (P7) | .08 | .17 | .22 | **.53** | Basal (7e-82) | miR-9/22/103a/27b |
| LumB:COMPLEMENT (P1) | .33 | **.38** | .22 | .08 | LumB (9e-52) | miR-21/10a/30d/10b |
| LumB:ALLOGRAFT_REJECTION (P8) | .28 | **.34** | .31 | .06 | LumB (1e-45) | miR-375/625/10a/210 |
| **Basal:APOPTOSIS** (P3) | .13 | .19 | .26 | **.41** | Basal (2e-44) | miR-21/183/30e/361 |
| LumA:APOPTOSIS (P2) | **.36** | .18 | .27 | .19 | LumA (8e-41) | miR-21/183/10b/127 |

- **Signed (gain+loss in one decomposition) gene-pressure NMF** — ±-channel stack, direction preserved
  (`nmf_gene_pressure_signed_factor_subtype.tsv`). **Every program is gain-dominated (gain-mass
  0.97–0.99)** → there is no organised pressure-*relief* program, so the gain-only acquired NMF discards
  nothing; subtype divergence persists (KW q ≤ 1e-44):

(keyed on stable `factor_label`; F-index in parentheses. **All 12 factors**, ordered by KW q —
i.e. most subtype-divergent first; the near-flat **backbone S3 sits near the bottom**. Subtype fracs
are LumA/LumB/Her2/Basal. The two `Basal:MYC_TARGETS_V1` rows are *distinct* factors — S12 is the
sharp, genuinely Basal-divergent MYC program; S9 a flat Basal-leaning secondary.)

| factor_label (F-idx) | gain frac | LumA/LumB/Her2/Basal | dominant (KW q) | top gain genes | distinctive arm(s) |
|---|--:|---|---|---|---|
| **Basal:MYC_TARGETS_V1** (S12) | .99 | .10/.21/.17/**.51** | Basal (8e-95) | KDM6B;ALOX15;DAB2;ARHGDIA;KPNB1 | miR-21/93/143/660 |
| **Basal:E2F_TARGETS** (S1) | .98 | .08/.17/.24/**.51** | Basal (8e-76) | CCNG1;SOCS5;CUL4A;STMN1;RAB34 | miR-9/99a/92a |
| LumA:E2F_TARGETS (S10) | .94 | **.43**/.28/.08/.21 | LumA (1e-53) | DONSON;WNT2;APBB2;CBX3;NT5E | miR-30a/30a-3p/155 |
| LumB:ALLOGRAFT_REJECTION (S5) | .95 | .28/**.34**/.31/.06 | LumB (1e-47) | KLF5;EIF4G3;RPN1;JAK2;LDHB | miR-375/125b |
| LumB:ESTROGEN_RESPONSE_LATE (S6) | .96 | .33/**.36**/.24/.08 | LumB (5e-40) | BACH1;ANXA1;S100A9;KRT5;CDKN1B | miR-196a/let-7a-d |
| LumA:ESTROGEN_RESPONSE_EARLY (S2) | .97 | **.35**/.18/.26/.20 | LumA (5e-32) | PCBP1;PTPN14;IL12A;TPM1;CCL20 | miR-21/127/379/183 |
| Basal:ESTROGEN_RESPONSE_LATE (S8) | .98 | .16/.18/.28/**.38** | Basal (1e-23) | IL1A;PCBP1;PTPN14;TPM1;IL12A | miR-21/142-3p/191/183 |
| Basal:MYC_TARGETS_V1 (S9) | .98 | .17/.17/.28/**.38** | Basal (1e-19) | LIFR;ABCE1;IGFBP5;ABL1;RAN | miR-203a/197/125a |
| **Her2:KRAS_SIGNALING_UP (S3) — BACKBONE** | .93 | .27/.27/**.30**/.15 | Her2 (3e-17) | RASGRP1;MSH2;CCL20;VHL;HPGD | miR-21/155 (+ universal miR-26a/let-7) |
| Her2:EMT (S7) | .98 | .17/.28/**.29**/.26 | Her2 (5e-16) | EFNA3;RAD52;SDHD;PCBP1;PTPN14 | miR-210/127/218 |
| LumB:E2F_TARGETS (S4) | .98 | .28/**.34**/.15/.23 | LumB (4e-08) | MIF;CAB39;EPOR;HELLS;DNAJC3 | miR-451a/144/92b |
| LumB:IL2_STAT5_SIGNALING (S11) | .94 | .26/**.28**/.27/.19 | LumB (3e-07) | FGF9;AKAP12;GPX4;CADM1;KIF2A | miR-183/182/186/324 |

  **Reading the full catalog (the secondary structure under the backbone).** Three facts the 5-row
  excerpt hid:
  1. **The backbone is *not* the most divergent factor — it is one of the least** (S3, KW q 3e-17, the
     flattest non-trivial profile). The subtype signal lives in the *high-KW top* of the table; the
     dominant-in-725 factor is precisely the one that does **not** separate subtypes — exactly what
     "shared everything-up gain" predicts. So "most divergent" and "most dominant" are *opposite* ends.
  2. **Proliferation is the divergent Basal axis, twice over:** the two sharpest factors are Basal MYC
     (S12) and Basal E2F (S1), both ≈.51 Basal, KW ≤8e-76 — Basal tumours concentrate acquired pressure
     on MYC/E2F proliferation targets (STMN1, CUL4A, CCNG1, DAB2). The Luminal/Her2 secondaries are
     hormone + immune (ESTROGEN_LATE/EARLY, ALLOGRAFT, IL2_STAT5), never proliferation.
  3. **Driver arms split into a universal core + a distinctive tail.** `miR-26a-5p` and `let-7i-5p`
     are top-8 drivers of **all 12** factors (`let-7g-5p` 11/12) — the carriers of the shared backbone.
     Each factor's *identity* is its distinctive tail (right column): miR-21/93/143 (Basal-MYC S12),
     miR-375/125b (LumB-ALLOGRAFT S5), miR-183/182 (LumB-IL2_STAT5 S11), miR-30a (LumA-E2F S10). No
     factor is loss-loaded (gain frac .93–.99) — there is no organised pressure-relief program anywhere.

> **Single-source-of-truth note (signed NMF).** The per-sample dominance table
> (`nmf_sample_signatures/`) and the factor tables above (`nmf_gene_pressure_signed_{loadings,
> factor_subtype}.tsv`) are now derived from the **one** comovement fit: `_nmf_signed_gene_pressure`
> persists the full `W`/`H` (`nmf_gene_pressure_signed_{W,H}.tsv`) and `nmf_sample_signatures` reuses
> them instead of re-fitting an independent copy. Two fits returned permutation-arbitrary indices
> (the backbone landed on S5 in one, S9-vs-S3 in the other), which previously mislabelled the
> cohort-wide backbone. Factor indices are still permutation-arbitrary across *rebuilds* — always join
> by `factor_label`, never the bare F-index.

- **Share NMF** (`_nmf_share`) — NMF on the *regulatory-share* matrices instead of raw pressure (the
  `c_tumor` NMF is dominated by the abundance backbone, so normalising exposes complementary structure).
  - **Dominant-regulator programs** (`nmf_share_dominant`, NMF on `gene_share_tumor`). The top factor D1
    is the same abundance backbone (loading corr ≈1 with the `c_tumor` NMF); the **secondaries are
    uncorrelated (≈0)** and isolate arms that *own* their few targets' budgets despite low global
    abundance:

| factor | top arms | top genes owned |
|---|---|---|
| D1 | miR-21/9/29a/let-7e/183 (backbone) | TPM1;BASP1;STRN;IL12A;DOCK4;DNM1L;PCBP1 |
| D2 | miR-125b/125a/221/139/199a | DRAM2;ENPEP;ICAM2;IKZF4;ACADS;SCNN1A |
| D3 | miR-24/192/191/23a/22 | NCAN;CTSD;FEN1;INSIG1 |
| D4 | miR-27a/27b/26b/29c/181a | ZBTB10;SLC6A8;SYK;CD55;LAMP2 |
| **D5** | **miR-200c/200b/200a/23a/148a** | **ETV5;ITGB4;FLNA;FBLN5;EFNA1** (EMT) |

  - **Grip-change programs** (`nmf_share_grip_change`, signed NMF on `gene_share_tumor − gene_share_gtex`,
    gain/loss channels): which arms *co-gain* vs *co-lose* regulatory grip across genes, independent of
    pressure inflation — the share analogue of the signed gene-pressure NMF and the only NMF reading
    *competitive restructuring*:

| factor | direction | top arms | genes (grip moved) |
|---|---|---|---|
| **G1** | **gain** | **miR-21/17-3p/217/323b/21-3p** | **VHL;MSH2;RASGRP1;HPGD;DDAH1** (TSG/repair) |
| G4 | gain | miR-210/375/135b/449a | DDAH1;FOXN3;ISCU;IGFBP3;BNIP3 |
| G5 | gain | miR-105/184/203a/503 | ADAMTS5;TJP1;PDK1;MYB;CDK6 |
| G3 | loss | miR-7/103a/29a/34c | HELLS;ARRB1;RELA;SERPINB5 |
| G6 | loss | miR-494/204-5p/335 | ARHGAP5;HNRNPA3;ATF3;HDAC2 |

**Subtype interpretation.** The *heterogeneity* axis (per-sample) carries the subtype signal; the
*target-sharing* axis (static-c) does not. miR-375-luminal, miR-9/210-basal-hypoxia, miR-142/155-immune
and the basal miR-17~92 program are literature-concordant; the **within-subtype factor decomposition**
(e.g. Her2 splitting into epithelial S1 + immune S2) is the novel resolution this platform adds.

### 6c. Within-pathway NMF — sub-modules *inside* one Hallmark program (`within_pathway_nmf`)
Everything above factorises the **whole** ~1,410-gene universe and labels each factor by its dominant
Hallmark. `within_pathway_nmf` does the opposite: it restricts the matrix to **one set's genes** and
decomposes *that* (33/50 sets with ≥40 targetable genes; per-set k tuned by a reconstruction-error
elbow), surfacing co-regulation **sub-modules inside a program** that the cross-set NMF averages away.
Two complementary families:

**Family A — gene×sample pressure (who is co-pressured across patients), 5 pressure flavours.** Run for
abundance (`softmax_logrpm`), the canonical **spine** (`softmax_z_logrpm`, the pressure §1–5 couple on),
**hybrid** (`softmax_z`+combined_mass), **acquired** (dev-anchored), **signed** (gain|loss); each at
cohort + within-PAM50 scope. The flavours are **not interchangeable** — the per-set elbow assigns each a
different complexity, and *the chosen k is itself the signal of resolvable sub-structure*:

| flavour | median k | median eff. #sub-modules | reading |
|---|--:|--:|---|
| abundance (no-z) | **2** | 1.7 | collapses to ~one module — the few high-abundance arms own the set's pressure |
| **spine** (`softmax_z_logrpm`, canonical) | **5** | 2.6 | **richest** — z removes the abundance backbone, exposing balanced sub-modules |
| hybrid (`softmax_z`) | 5 | 2.9 | richest (relative-share) |
| acquired (dev-anchored) | 4 | 2.7 | moderate; co-gain structure |
| signed (gain\|loss) | 4 | 2.9 | moderate; gain/relief structure |

So the **z-spine is the resolving lens, not abundance** — the opposite of what a naive "% dominant" read
(k-confounded) suggests, and the reason all flavours are worth keeping (abundance alone would have hidden
the sub-structure). Under the spine flavour the sub-modules are sharp and subtype-split:

| set (spine, k) | sub-module (top genes; lead arms; dominant subtype, KW q) |
|---|---|
| **EMT** (k=4) | **f2 inflammatory-contractile** MYLK/INHBA/VEGFA/TIMP3/IL6 (miR-130b/103a/93; **Basal .49, q1e-109**) · **f3 collagen** COL4A1/COL1A2/SPARC/COL3A1/LOXL2/FBN1/COL5A2 (**miR-29a/b/c**; pan-subtype, **q≈1**) · f4 mesenchymal-marker VCAN/VIM/ITGB3/THBS1 (miR-30a/182; LumA) · f1 SDC1/IL6/CD44/CDH11 (miR-10a/143/200c) |
| **E2F** (k=5) | f2 replication RAD21/MTHFD2/STMN1/SHMT1/CKS1B (miR-92a/197/9; **Basal .56, q7e-104**) · f3 mitotic RRM2/TOP2A/AURKB (let-7; LumA) · f4 splicing/checkpoint SRSF1/TRA2B/CDKN2A/PLK1/TP53 (miR-10a/b; LumA) · f5 MYBL2/DNMT1/CHEK1 (miR-30e/16/196a; Basal) |

The **miR-29a/b/c→collagen** EMT sub-module (f3) independently recovers the orphan-discovery miR-29→ECM
axis (MH-36/38) — internal validation that the decomposition is real — and is correctly flagged
**pan-subtype / KW-non-significant**, consistent with collagen co-varying with stroma (MH-39 caveat).
**Within-subtype** re-fitting adds structure for some programs: APOPTOSIS goes from ~1.5 effective
sub-modules cohort-wide to **3.2 inside Her2 / 2.3 inside Basal**, HYPOXIA to **3.6 inside LumB** — i.e.
a program that looks like one module cohort-wide resolves into several when the subtype is fixed.

**The recurring grammar (what is actually there).** Across all 11 focal programs the same **three
sub-module archetypes** recur, each driven by a distinct arm family — i.e. the platform finds the *same
regulators partitioned the same way* inside every program, which is the main descriptive result:

| archetype | driver arms | genes it lands on (examples by program) |
|---|---|---|
| **Basal replication / TSG-repression core** | miR-17~92 (17-5p/17-3p/19), miR-92a-3p, miR-25-3p, miR-197-3p, miR-203a-3p | RAD21/STMN1/RAN/TYMS (E2F/G2M/MYC) · PTEN/CDKN1A/CDKN1B/FOXO3/RB1/FBXW7 (PI3K/P53/HYPOXIA) · CD69/ICAM1 (IFN-γ) |
| **LumA mitotic / splicing** | let-7a/e, miR-30a/c, miR-10a/10b | TOP2A/AURKB/RRM2/KIF11 · SRSF1/TRA2B/PLK1 |
| **immune antigen-presentation / evasion** | **miR-375**, miR-200a-5p, miR-148a-3p | TAP1/STAT1/SOCS1/JAK2 (IFN-γ, ALLOGRAFT) |

So the **basal oncomiR hub is not one module** — at this resolution it appears *separately* as the
replication core in E2F/G2M/MYC **and** the TSG-repression core in PI3K/P53/HYPOXIA/IFN-γ (the same hub
MH-17/18/19 flagged, now split into its per-program effector vs suppressor arms); and **miR-375 carries a
distinct antigen-presentation/immune-evasion sub-module** in the immune programs (the §7 miR-375→ALLOGRAFT
pro-tumor call, localised to TAP1/STAT1).

**How the sub-modules distribute across tumours (sample signatures).** From `per_sample.tsv.gz` (spine,
cohort): **60–71% of tumours are dominated by a single sub-module** (median top-share 0.54–0.60) — much
more discrete than the cross-set NMF (4–45%, §6b). But the shape is a **fractal of the §6b backbone
story**: each program = **one large pan-subtype "backbone" sub-module that most tumours share + a few
subtype-specialist sub-modules that cleanly partition** a minority. The *biologically interesting* modules
are the specialists, not the backbone:

| program | backbone sub-module (spans subtypes) | subtype-specialist (cleanly partitions) |
|---|---|---|
| HYPOXIA | f2 CDKN1B/FOXO3/VEGFA — **685/1096**, all subtypes | f3 LumA SERPINE1/RORA (186/248 LumA) |
| PI3K/AKT | f2 MAPK8/ATF1 — **697/1096** | f6/f3 **PTEN/CDKN1A** (Basal/LumB-enriched, 25+132) |
| IFN-γ | f4 TAP1/STAT1 — **519/1096** | f1 **Basal CD69/ICAM1 — 8/8 Basal** (sharpest) |
| EMT | f2 inflammatory (340) + f3 collagen (303) | f1 LumA SDC1 (222/318 LumA), f4 LumA VCAN (78/97) |

i.e. PTEN-repression, the basal CD69 immune-restriction module, and the LumA hormone/ECM modules are
**minority signatures dominant in a focused tumour subset**, while the program-wide pressure backbone is
shared and subtype-flat.

**Family B — edge-level miRNA×gene lenses (who *co-regulates / owns / restructures* the program).** A
genuinely **different matrix** (not a pressure multiplier; no sample axis — the within-pathway analogue
of §6a's `nmf_share_*`): each factor pairs an **arm module** with a **gene module**. Three lenses:
- **share** (`gene_share_tumor`, who *owns* the budget): EMT splits the program's regulatory ownership
  into a **miR-29b/c/a** module (the collagen genes), a **miR-21/9/223/210** module, and a
  **miR-200c/127/27b** module — ownership tracks the same sub-modules Family A found, from the regulator
  side.
- **grip_change** (`share_tumor − share_gtex`, competitive restructuring): the clean tumor story — in
  EMT, **miR-21/miR-17-3p *gain* grip** on RHOB/TGFBI/COL4A1/VEGFA/TIMP3 while **miR-29a/b/c *lose* grip**
  on FBN1/COL4A1/COL5A2/SPARC/collagen; in E2F, **miR-21/449a/b gain grip** on MSH2/CDC25A/MYC/CHEK1
  while miR-7/miR-1246 lose grip on HELLS/CCNE1/CDK4/EZH2. So the program's *control* is being handed from
  the suppressor families (miR-29) to the oncomiR (miR-21) — a restructuring neither the pressure NMF nor
  the topology layer (§7) shows.
- **aggregate** (`c_tumor`): co-regulating arm-module + co-pressured gene-module, the within-set analogue
  of the §6b aggregate NMF.

**Novelty — be precise.** Three separable claims, only the middle one is the platform's contribution:
1. *The individual arm→gene modules are textbook* (LOW novelty — they function here as positive
   controls / internal validation, not discovery): miR-17~92→PTEN/replication, **miR-29→collagen**
   (the canonical ECM-miRNA axis), let-7/miR-30→mitotic, miR-375 as an oncomiR. If these had *not*
   appeared the decomposition would be suspect; that they do is the evidence it is real.
2. *The within-program **resolution** is new* (MODERATE, methodological): standard Hallmark/GSEA gives
   **one score per set** and the cross-set NMF (§6b) gives **one factor per program** — neither shows that
   a program is internally a **shared pan-subtype backbone + named subtype-specialist sub-modules**, nor
   *which arm drives each*. "Same known players, newly resolved organisation," not new players. The
   sharpest specialist (Basal CD69/ICAM1 immune-restriction in IFN-γ, 8/8 Basal) and the localisation of
   PTEN-repression to a *minority* PI3K/P53 sub-module are the most useful new reads.
3. *The grip_change handover is a state-resolved competitive read* (P/H, like the §3b brake-release): the
   pairs (miR-21→EMT effectors, miR-29→collagen) are known; the **engaged-in-healthy → handed-to-oncomiR
   framing** is the contribution, and it is descriptive (a share-Δ, **not** an FDR-tested coupling).
Plus one genuine **methods** finding: under a principled per-set elbow the canonical z-spine resolves
sub-structure (k≈5) while the no-z abundance attribution does not (k=2) — relevant to anyone running
gene-pressure NMF.

**What this is / is NOT (so "what's there" is unambiguous).** IS: a reproducible **co-pressure**
decomposition (static edge weight × abundance, summed per gene) → which arms pressure which genes of a
program together, how that splits by subtype, and how tumours distribute over the sub-modules. Is **NOT**:
(i) a coupling/repression claim — there is **no measured anti-correlation and no FDR test on the modules
themselves** (the only test is KW for subtype divergence of a module's activity); (ii) causal or
externally validated; (iii) composition-free — abundance/acquired pressure partly tracks stroma/immune
content (the collagen and some immune modules especially, MH-39). Also: per-set k is elbow-chosen not
cross-validated, and factor indices are permutation-arbitrary across rebuilds (join by content, not
`f#`). Outputs under `tissue_reference/within_pathway_nmf/{factor_summary,set_dominance,per_sample.tsv.gz,
lens_factors}.tsv` + per-focal-set heatmaps. See **MH-42**.

---

## 7. Gene-set architecture (Stage 3)

Overlays each Hallmark set's regulatory **topology** (live OmniPath directed+signed interactions,
cached) onto otherwise-flat pressure. Per-gene weight = reverse-PageRank (master-regulator influence)
÷ redundancy-group size. The default effect model is **signed-path propagation with a calibrated
continuous edge sign**, then a tumor-direction prior; the original immediate-neighbour binary sign is
kept only as a robustness sibling (below).

**Default direction is pathway *damage*, not activation.** Most program nodes are activators, so a
repressing miRNA *reduces* program output: across all 50 sets, **7,840** (set,miRNA) pairs are net
pathway-repressing vs only **114** net-activating. The largest effects are damaging
(COMPLEMENT×miR-21-5p, GLYCOLYSIS×miR-183-5p, EMT×miR-30a-5p). The net-*activating* cases
(an arm that raises output by repressing **inhibitor** nodes) are the rare, counter-intuitive calls
worth highlighting — only **11/50 sets** carry ≥5 inhibitor nodes (APOPTOSIS 13, COAGULATION 12, E2F
12, G2M/PI3K 9, IFNγ/P53 8); the rest are effector bags where `arch ≈ flat`.

**Per-set aggregate** (`architecture_all_set_summary.tsv`) — the overall miRNA pressure on each whole
program, the complement to the across-set (per-arm) and within-set (per-gene) views. `total_program_
output_change = Σ_m program_output_change(m)` is the summed effect of the *entire* miRNA layer on the
set; **all 50 sets are net `damaged`** (no program is net-amplified by miRNAs). `total_pro_tumor_pressure
= prior · total` is positive *only* for the suppressive (−1) programs — the layer's malignant net effect
is **brake-release on suppression, not engine amplification** — and is `ambig` for the 27 prior=0 sets
(direction kept, polarity not forced).

| Hallmark | prior | inhib nodes | total output Δ | net pro-tumor (acq) | top pro-tumor arm |
|---|--:|--:|--:|--:|---|
| APICAL_JUNCTION | 0 | 4 | −607 | ambig | — |
| G2M_CHECKPOINT | +1 | 9 | −559 | −559 | (engine damaged) |
| EMT | +1 | 1 | −508 | −508 | (engine damaged) |
| TNFA_NFKB | 0 | 5 | −452 | ambig | — |
| **P53_PATHWAY** | −1 | 8 | −400 | **+400 (98)** | **miR-21-5p** |
| **APOPTOSIS** | −1 | 13 | −390 | **+390 (87)** | **miR-21-5p** |
| E2F_TARGETS | +1 | 12 | −383 | −383 | (engine damaged) |
| **ALLOGRAFT_REJECTION** | −1 | 4 | −381 | **+381 (74)** | **miR-375** |
| **INFLAMMATORY_RESP** | −1 | 2 | −351 | **+351 (53)** | miR-21-5p |
| **INTERFERON_GAMMA** | −1 | 8 | −251 | **+251 (47)** | miR-22-3p |

The +1 proliferation/metabolism engines (G2M/EMT/E2F) carry the *largest raw damage* but **negative**
pro-tumor totals (miRNAs net-damage them); the positive pro-tumor totals are exclusively the −1
suppressive programs — the aggregate-set restatement of the brake-release thesis.

**Gene-role malignancy per-set aggregate** (`total_mal_pro_tumor = Σ_m mal_pro_tumor`, + `total_tsg_credit`
/ `total_onco_collateral` / `malignancy_direction`). The total above uses the **coarse prior** (one
hand-assigned sign per set, `ambig` for the 27 prior=0 sets). This second aggregate is **prior-free** —
it sums the curated TSG/oncogene-role pressure (§7b) over the whole set, so it is defined for **all 50
sets** and resolves the ambiguous ones. It validates against the prior (prior=−1 suppressive sets are
malignancy-pro-tumor in **5/6**) and ranks **P53_PATHWAY #1** (net +49.1; TSG credit 163.5 ≫ onco
collateral 57.5):

| Hallmark | prior | **mal total** (tsg − onco) | acquired | malignancy dir |
|---|--:|--:|--:|---|
| **P53_PATHWAY** | −1 | **+49.1** (163.5 / 57.5) | +18.0 | pro_tumor |
| **TNFA_NFKB** | **0** | **+17.9** (38.8 / 32.9) | −0.6 | pro_tumor (prior was `ambig`) |
| TGF_BETA | +1 | +16.1 (11.8 / 8.8) | +1.7 | pro_tumor |
| INTERFERON_GAMMA | −1 | +12.0 (48.7 / 32.9) | −0.6 | pro_tumor |
| **HEME_METABOLISM** | **0** | **+8.1** (17.8 / 6.4) | +3.8 | pro_tumor (prior was `ambig`) |
| **ADIPOGENESIS** | **0** | **+7.0** (7.9 / 1.1) | +2.0 | pro_tumor (prior was `ambig`) |

Two caveats make the two aggregates **complementary, not interchangeable**: (i) the malignancy total only
counts the ~232 curated drivers, so it is **high-specificity / low-coverage** — ALLOGRAFT scores +381 by
prior (whole-set immune-evasion brake-release) but only +8 by malignancy (few of its genes are curated
drivers), so a small `total_mal_pro_tumor` means "few driver genes hit," not "not oncogenic"; (ii) most
sets are net **anti-tumor** by gene-role (31/50) — the layer mostly lands collateral on oncogenes — but
19 are net pro-tumor where TSG repression dominates.

### 7a. Calibrated sign + propagation + prior → pro-tumor pressure (the default readout)
**Why keep the calibrated continuous sign at all** (vs just the gene-role malignancy score of §7b)?
Because it is the **only** readout that captures *network-mediated, indirect* effects: an arm
repressing an upstream hub changes the whole downstream program even if it never touches an oncogene/TSG
directly. Flat pressure, anti-correlation, and the gene-role score are all **local** (per-edge or
per-gene); the signed propagation is what turns "this miRNA hits these nodes" into "this is the net
effect on *program output*". It is therefore retained as the **topology/pathway lens**, explicitly
paired with — not replacing — the gene-role malignancy call (§7b): the two answer different questions and
miR-21 is the canonical case where they must both be read. The *continuous* (vs binary) calibration is
kept because the hard ±1/0 discretisation demonstrably over-counted dual-edge scores (below).
- **Calibrated continuous sign.** Each ordered gene→gene pair aggregates its OmniPath records to
  `net_sign = (n_stim − n_inh)/(n_stim + n_inh) ∈ [−1,+1]`, so dual/conflicting edges contribute a
  *graded* value instead of being flattened to 0. It inherits OmniPath's annotation completeness, so the
  **binary `effect_sign` / `arch_signed_pressure` are retained as a sensitivity check** — they agree on
  the TS leaders and family ordering, and the one place they disagree loudly (miR-21-5p) is exactly where
  the graded encoding is more defensible.
- **Signed-path propagation.** A finite-horizon (Katz, α=0.25, 6 hops) signed-path sum on the `net_sign`
  adjacency gives each gene `output_influence(g)` = its net effect on the rest of the program; an arm
  repressing `g` changes program output by `−c_tumor·output_influence(g)`, summed = `program_output_change(m)`.
- **Tumor-direction prior** (coarse BC polarity): +1 pro-tumor (EMT/HYPOXIA/GLYCOLYSIS/PI3K/E2F/G2M/
  MYC/…), −1 suppressive (APOPTOSIS/P53/IFN/INFLAMMATORY/ALLOGRAFT), 0 ambiguous. **`pro_tumor_pressure
  = prior · program_output_change`** (`sum_pro_tumor_pressure` per arm in `architecture_mirna_cross_set.tsv`).

**This recovers textbook biology (validation).** Net **anti-tumor** leaders are the canonical TS
miRNAs — miR-30a-5p (−116.9), miR-24-3p (−78.7), miR-29b-3p (−63.0), miR-203a-3p (−57.2), let-7a-5p
(−53.5), miR-29c-3p (−49.4) — damaging pro-tumor engines across many sets each. Net **pro-tumor**
leaders are known oncomiRs acting by **brake-release**: miR-200a-5p (+26.7, top P53), miR-375 (+23.8,
top ALLOGRAFT = immune-evasion), let-7f-5p (+14.6), miR-222-3p (+9.8), miR-130a-3p (+9.2), let-7d-5p
(+9.0).

**The key mechanistic readout:** for the +1 proliferation/metabolism engines (E2F, G2M, PI3K, HYPOXIA,
MYC) the top `pro_tumor_pressure` is ≈0/negative — miRNAs **damage** those engines, they don't amplify
them. Positive pro-tumor signal concentrates on the −1 **suppressive** programs, i.e. malignancy via
this layer is driven by miRNAs **releasing brakes** (apoptosis, p53, interferon, immune rejection), not
revving engines. Reinforced by the ambiguous-set decomposition (§7c): **all 27 ambiguous sets are
net-damaged too** — the whole miRNA layer is net-repressive on every program; only the *targeted*
sub-pattern is oncogenic.

### 7b. Gene-role malignancy score (and the miR-21 resolution)
The topology score has a blind spot for canonical oncomiRs: it scores miR-21-5p as net-neutral
(`pro_tumor_pressure` −0.6) because miR-21 mostly *damages* the +1 engines it broadly represses (it has
76 targets), and its canonical oncogenic targets are mis-handled (PTEN acts on a **lipid (PIP3)**, not a
gene node — a true ceiling of gene-gene topology; PDCD4 surfaces only via another arm). The fix is a
curated **gene-role** layer (`gene_roles.py`, 232 high-confidence drivers: 114 oncogene, 106 TSG, 12
dual; COSMIC-CGC/OncoKB consensus) giving each gene a `malignancy_sign`. The **default malignancy score
is `mal_pro_tumor` = Σ −malignancy_sign · c_tumor** (credits repressing a TSG as pro-tumor, independent
of the graph), with the architecture-weighted **`mal_pro_tumor_arch` = Σ −malignancy_sign · c_tumor ·
w_arch** (gene-role direction × edge strength × master-regulator centrality; `w_arch` always positive,
so `tsg_credit` is non-negative and `onco_collateral` clean).

| miRNA | topology `pro_tumor` | **`mal_pro_tumor`** | `mal_pro_tumor_arch` | **acquired** | TSG / onco targets |
|---|--:|--:|--:|--:|---|
| **miR-21-5p** | −0.6 (neutral) | **+34.2 (#1)** | +28.5 | **+31.8 (93%)** | 4 TSG (PTEN/PDCD4/FOXO axis) / 3 |
| miR-182-5p | −28.0 | +25.8 | +32.2 (top) | +11.6 (45%) | 3 / 2 |
| miR-374a-5p | +1.0 | +18.2 | +17.5 | ~0 (constitutive) | 1 / 1 |
| miR-10b-5p | −17.9 | +16.2 | +21.7 | +4.4 (27%) | 4 / 0 |
| miR-93-5p | −2.3 | +14.8 | +15.2 | +4.6 (31%) | 3 / 2 |

The malignancy ranking is a who's-who of validated BC oncomiRs (miR-21/182/10b/93/106b/25/221/9/10a) —
external validation the gene-role layer captures what topology missed. **Both scores are kept and
answer different questions:** topology = "what the pathway wiring says," gene-role = "what the driver
annotation says." **Rule: for high-breadth abundant arms (n_targets ≫ median 3) read the gene-role /
coupling / acquired axes, not the topology `pro_tumor` sign.**

**Acquired malignancy** (`sum_mal_pro_tumor_acquired`, weight `max(c_tumor − c_gtex, 0)`) is computed
for every arm and set (`total_mal_pro_tumor_acquired`). It sharpens the call: **miR-21's TSG repression
is 93 % tumor-acquired** (+31.8 of +34.2) — not a constitutive trait — whereas miR-374a-5p (#3 by total
abundance, +18.2) is **constitutive** (drops out of the acquired top), as do miR-10a/452/210/150. So
"miR-21 acquires TSG repression in tumor" is a stronger, collapse-safe statement than the abundance score.

> **Why malignancy direction is binary + driver-only (and what already weights genes continuously).**
> The all-genes, *non-binary, downstream-weighted* score is the **topology** score, not this one:
> `output_influence(g)` = the finite-horizon signed-path (Katz, α=0.25, 6 hops) sum on the calibrated
> `net_sign` adjacency, so every gene is weighted by the *propagated* effect it has on the rest of the
> program (a gene driving many downstream effectors gets a large value, a leaf ≈1, an inhibitor flips
> sign) — that is precisely "how it knows the weighted downstream effect," and the calibrated continuous
> sign feeds it as the edge weights. What topology lacks is a *per-gene malignancy direction*; it uses a
> coarse per-set prior instead. The gene-role score is the opposite trade — per-gene curated direction,
> but it is binary ±1 and only for the ~232 drivers **because per-gene tumor-direction is only known for
> those genes**; we deliberately do not fabricate a polarity for the ~unlabeled set members (averaging a
> meaningless set polarity is the documented "honest limit", §7c). The *magnitude* is already graded for
> drivers via `mal_pro_tumor_arch` (× `w_arch` centrality) — only the *direction* is categorical. We do
> **not** multiply `malignancy_sign` by signed `output_influence` (the rejected original
> `mal_pro_tumor_arch`): it double-counts direction (a TSG with negative downstream influence flips the
> TSG credit) and lets one hub edge dominate. A genuine all-genes *graded role weight* is now wired (see
> the continuous score below).

#### Continuous, all-gene role weight (`mal_pro_tumor_cont`) — DepMap dependency
The binary `malignancy_sign` covers only the 232 curated drivers. The graded counterpart uses **DepMap
26Q1 CRISPR (Chronos) gene-effect over the 96-line breast panel** (`gene_roles.load_gene_dependency` →
`data/CCLE/depmap_gene_effect_summary.tsv`, builder `_build_depmap_dependency.py`): `dep_role_weight(g) =
gene_effect(g)`, oriented like `−malignancy_sign` — strongly **negative** gene-effect = essential =
tumour *dependency* (oncogene-like, repressing it is anti-tumor), `≳0` = TSG-like (repressing it is
pro-tumor). `mal_pro_tumor_cont = Σ_g dep_role_weight·c_tumor` over **every gene in the set** (714/721
arms get a non-zero score vs the driver-restricted binary). Key results:

| view | binary (curated) | **continuous (all-gene DepMap)** |
|---|---|---|
| top anti-tumor arm | miR-30c-2-3p −39.3 | **miR-30c-2-3p −93.4, miR-24-3p −91.2** (binary near-0 — flipped) |
| arm coverage | curated-driver targets only | **714/721 arms** |
| most-hit sets | (driver-sparse) | **G2M −317, E2F −287, MYC −168, MITOTIC_SPINDLE −149** (proliferation dependencies) |

Two things this adds: (i) it **resolves the binary coverage blind spot** — miR-24-3p / miR-30c / let-7a
repress a large mass of breast-cancer **dependency** genes that aren't curated drivers, so binary scored
them ≈0 while continuous ranks them the top anti-tumor arms; (ii) at the **engine** sets (G2M/E2F/MYC)
the continuous dependency view shows the miRNA layer lands hard on proliferation dependencies — the
concrete, all-gene restatement of "miRNAs damage the engines" (§7a). **One-sided caveat**: CRISPR
resolves the dependency (anti-tumor, negative) tail well but the TSG (pro-tumor, positive) tail weakly
(max +1.6), so the continuous score is the **all-gene graded onco-collateral / anti-tumor axis** and the
curated `tsg_credit` remains the better *pro-tumor* detector — they are complementary
(corr 0.29). An acquired variant (`mal_pro_tumor_cont_acquired`) is computed in parallel.

### 7c. Family/total rollups, ambiguous sets, acquired-weighted
- **Family** (`architecture_family_*`, name-stem proxy — not exact TargetScan seed): net **anti-tumor
  families** dominate (miR-30-5p −193, miR-29-3p −149, let-7-5p −119, miR-24-3p −79, miR-26-5p −71); net
  **pro-tumor** (miR-200-5p +26.7, miR-375 +23.8, miR-130-3p +10.6, miR-222-3p +9.8, miR-196-5p +9.7).
- **Total layer** per set — the per-set aggregate table (`total_program_output_change` / `net_program_
  direction` / `total_pro_tumor_*`) is given in the §7 overview above; it is the denominator for "what
  fraction does arm/family X carry."
- **Ambiguous (prior=0) sets** are not forced a polarity but keep the directional readout: all 27 are
  net-damaged, and per-(set,miRNA) `program_output_change` is available for each. Gene-level
  decomposition (`prog_change(g) = −Σ_m c_tumor·output_influence`) localizes the damage to interpretable
  hubs — APICAL_JUNCTION onto **SRC/FAK/EGFR/ITGB4** (invasion hubs → plausibly anti-invasion);
  TNFA_NFKB onto **NFKBIA (the IκB inhibitor) *and* IL6/CCL20/SERPINE1** (a genuine tug-of-war, *why*
  it is ambiguous). The honest limit: turning that into "the arm's effect" needs **per-gene** polarity
  (now available via `gene_roles` for the 232 drivers; the rest stay unlabelled rather than averaging a
  meaningless set polarity).
- **Acquired-weighted parallel** (`pro_tumor_acquired`, weight = `max(c_tumor − c_gtex, 0)`) isolates
  the pro-tumor effect *gained vs healthy*. Top acquired pro-tumor arms — miR-375 (+11.2), miR-17-3p
  (+6.7, P53), miR-200a-5p (+5.7), miR-455-3p (+4.0) — are all strong magnitude gainers (+6.7 to +10.8
  log2, mostly `magnitude_only`); the abundance-axis leaders let-7d/f, miR-125a, miR-130a **drop out**
  (their pro-tumor pressure is constitutive). So the acquired axis is the one to use for "what the
  *tumor* switched on."

### 7d. Per-set gene-level hits — hubs, specific genes, and malignancy (`architecture_set_gene_hits.tsv`)
The set-level scores roll up to a per-(set, gene) "what is hit" table that localizes each program's
incoming miRNA pressure to specific genes, reading **both** directional lenses at gene resolution:
`incoming_pressure_tumor` (Σ_m c_tumor on the gene) × `output_influence` = `prog_change_gene` (the
topology/hub view of program-output damage), and `mal_hit = −malignancy_sign · incoming_pressure` (the
gene-role view — pressure on a TSG scores pro-tumor). `is_hub` flags the set's top-decile centrality
genes. This makes the engine concrete:

| set | top hub genes hit (arm; role) | top malignancy (TSG) hits |
|---|---|---|
| **P53_PATHWAY** | **TP53** (miR-375; TSG, influence 10.4 → prog_change 57) ≫ WWP1 (miR-21), PLK2, FOXO3 (TSG) | PTPN14, CDKN2A, BTG2 |
| **PI3K_AKT_MTOR** | RAC1 (miR-224; onco), HRAS (miR-143; onco), MAPK8, TIAM1 | CDKN1B, MAP3K7, **PTEN** |
| **E2F_TARGETS** | AURKA / AURKB / PLK1 (let-7; onco mitotic hubs) | CDKN1B, CDKN2A, MSH2 |
| **APOPTOSIS** | CASP3 / CASP8 (miR-224/181c), ERBB2 (miR-375; onco) | CDKN1B, BTG2, TIMP3 |
| **EMT** | COL1A2 / FN1 (miR-200c) / TNC / SPP1 (ECM effector hubs) | DAB2, TIMP3, GADD45A |

The two lenses cohere with the §7a/§7b story: the **hub** hits are mostly program effectors (so
repressing them *damages* the program — let-7→AURKA/PLK1 damages E2F, the anti-tumor case), while the
**malignancy** hits are the TSGs whose repression is pro-tumor regardless of network position. The
single most striking hit is **TP53 itself being the dominant P53-program node hit, driven by miR-375**
(the immune-evasion oncomiR) — a hub-and-malignancy convergence the set-level score alone hid.

**Within-set *functional role* — topology + curated TF identity.** The gene's role *inside the program*
(master regulator vs effector vs inhibitor) is the topology, not the oncogenic role: `w_arch`
(reverse-PageRank upstreamness), `output_influence` (propagated downstream effect) and `effect_sign_soft`
(activator/inhibitor) are computed per set from the OmniPath-induced signed subgraph — so the
*continuous, set-specific* role is already there (MSigDB Hallmark itself gives only membership, no role).
What that topology **misses** is curated TF identity — only 10% of `w_arch` hubs are curated TFs and
**~48% of curated TFs fall below median `w_arch`** (transcriptional regulation is under-covered in the
signed-PPI network). So a curated **human-TF census** (Lambert 2018, 1,639 TFs; `gene_roles.load_tf_census`)
is added as an independent identity flag `is_tf` + `hub_tf` (= curated TF ∧ topology hub = high-confidence
master regulator), plus a manually tunable hierarchy weight `w_hier = w_arch · (1 + TF_HIER_BOOST·is_tf)`
(default boost 0.5; `w_arch` left untouched so the pure-topology scores remain). **`w_hier` is wired into
scored columns**: `mal_pro_tumor_hier` (= Σ −malignancy_sign·c·w_hier, the TF-up-weighted analogue of
`mal_pro_tumor_arch`; +acquired), `program_output_change_hier` (topology damage with master-regulator
TFs up-weighted), and `tf_pressure`/`n_tf_targets` (raw pressure landing on curated TFs) — at per-(set,
miRNA), cross-set (`sum_*`) and per-set (`total_*`) levels. It stays close to `arch` where no TF is hit
(488/721 arms identical) and diverges for the ~233 arms whose targets include curated TFs
(e.g. miR-10b/150/25/106b climb: more of their pressure routes through master regulators). At the set
level the TF-up-weighted damage re-ranks the engines — **TNFA_NFKB jumps −452→−516** (TF pressure 116,
the highest, via NFKB1/RELA/JUN/FOS) and G2M −559→−610 (via E2F/MYC). `hub_tf` isolates
the **transcriptional master regulators** hit by miRNA pressure — a sharper list than topology-hub alone
(which mixes in signaling hubs like RAC1/CASP3):

| set | top master-regulator (`hub_tf`) hits | driving arm |
|---|---|---|
| P53_PATHWAY | **TP53** (prog 57) | miR-375 |
| G2M_CHECKPOINT | **E2F1 / E2F2 / E2F3 / MYC / MYBL2** | miR-205/let-7a/miR-184/449a/30e |
| INTERFERON_GAMMA | **STAT1** | miR-338-3p |
| MYC_TARGETS / E2F | MYC, MYBL2 | miR-184/30e |

### 7e. Are the anti-tumor brakes themselves released?
The net anti-tumor leaders are the canonical TS miRNAs — a natural worry is that they are lost in tumor.
They are **not**: every top anti-tumor arm is rank-flat (`own_dHT ≈ 0`, `pressure_move = flat`) but
**rises in absolute abundance** (`own_log2fc +3 to +9`, `magnitude_only`) — the tumor expresses them
*more*, not less. So "brake release" is **not** a TS-miRNA-silencing phenomenon at the abundance level
(and where rank ever looked dramatic, that was the now-repaired GTEx-collapse artifact, which never
touched the tumor-side pressure). The genuine brake release lives in the **coupling layer**
(released brakes / R00 pattern, §3b): healthy realizes repression of ERBB2/MMP/NOTCH2, tumor loses the co-variation
while the TS miRNA is still present — redirecting mechanism to target buffering / AGO saturation / 3′UTR
shortening, exactly where §4 points and where a functional follow-up should look.

---

## 8. Limitations and next steps
- **Gene-gene topology cannot see lipid/metabolite-mediated inhibition** (PTEN→PIP3) — a true ceiling;
  the gene-role layer (§7b) is the partial workaround for canonical drivers.
- **De-coupling mechanism is not yet pinned.** CN (1/210) and promoter methylation (2/182) are ruled
  out, and ATAC is mostly null (§4a); the leading remaining candidate is 3′UTR APA/shortening +
  ceRNA/RBP competition, which needs the APM motif/SV/APA layer per sample.
- **ATAC is sample-limited, not gene-limited** (after the rebuild): only the ATAC-assayed tumors carry a
  direct signal. The strict direct-only test (§4a) is the proof-of-concept — it surfaces a real **4/134
  accessibility-linked minority at growth/cell-cycle loci** (miR-7-5p→PTK2, miR-495-3p→FOXC1,
  miR-652-3p→FGFR1, miR-486-3p→CCND1) that the cohort-wide subtype-imputation washes out — but it is
  underpowered; a larger per-tumor ATAC cohort would make both the null and this signal conclusive.
- **164 `floor0` arms** (no GTEx *or* miTED reference; miR-320a-led) are a healthy-reference coverage
  gap, not true acquisition — only ~22 are flagged acquired gainers (the residual false-acquired risk).
  Resolve by extending GTEx/miTED arm coverage, not by imputing absent arms (which would erase real
  gains). The healthy-baseline cache must be force-rebuilt whenever the arm universe changes.
- **Family key is a name-stem proxy**; an exact TargetScan seed-family map would tighten miR-200a vs
  miR-200b/c-type splits.
- **Prior is coarse** BC-context literature polarity; ambiguous sets are handled per §7c rather than
  guessed.
