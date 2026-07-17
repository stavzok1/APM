# miRNA locus GENOMIC-CONTEXT axis (`learned/genomic_context.py`)

> **Goal:** classify every HE arm by the genomic context of its locus (host-coupled vs autonomous, strand-aware) and
> propagate that into the CN instrument (host pleiotropy) and the coupling confounder block (host-gene confounding) —
> the structural property that governs *what a regulator's abundance means*.
> **What belongs here:** the class taxonomy + Phase-0 distribution, the arm→locus provenance canonicalization, the
> universe-redefinition it exposed, and the DOWNSTREAM dependency map (CN / coupling / multi-arm). NOT the proliferation
> special case's own results (→ memory `proliferation-confound-vs-mechanism`, registry MH-100).
> **Update trigger:** when the classifier, the coordinate source, or a downstream wiring changes.
> **Sync-partner:** `DATA_SOURCES.md` (arm→locus canonicalization), memory `universe-redefinition-pending-refresh`
> (the `.N` universe fix), `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md §12 BP` (host-gene instrument validity).

**Status: ARC COMPLETE (2026-07-11). Phase 0 BUILT + provenance CANONICALIZED; annotation completed (FANTOM5 co-transcription
`promoter_class`); downstream (a) CN-exclusion `host_pleiotropy` + (b) coupling `host_confound_lens` + (c) family context
`family_context` ALL BUILT. ✅ PER-LOCUS HOST now BUILT too (`locus_host_map()` → `locus_context.tsv`; each CN locus classified
at its OWN coordinate; gate/diagnostic/dose-attribution/`host_target_relatedness` all condition each ACTIVE locus on its own
host — dropped 2 silent-locus spurious down-weights, 403→444 related). Open extension only: the `coding=True` full-universe scan.
NB the exclusion hot path is now 6.4× faster (batch ~24→2.4 min, pure caching; base table identical).**

---

## 1. Why it is an axis (not a footnote)
"Intronic vs intergenic" changes **what a regulator's abundance means**: a same-strand (sense) intronic miRNA is
Drosha-processed co-transcriptionally from its host's pre-mRNA, so its abundance ≈ host transcription × processing;
an antisense-overlap or intergenic miRNA has its own promoter. This cascades into (a) **confounding** (the host's
program confounds the miR→target coupling), (b) **CN-instrument validity** (locus CN moves the host too → BP
pleiotropy), (c) an **AL host-mRNA transcription proxy** (for intronic regulators the host mRNA is a directly-available
co-transcription readout — partially UN-blocking the design's "most load-bearing" transcription-rate latent).

## 2. Class taxonomy (STRAND-AWARE) + Phase-0 distribution (714 HE arms)
| class | arms | % | meaning |
|---|---|---|---|
| **sense_coding_host** | 337 | 47% | co-transcribed from a **protein-coding** host → host program **confounds** + host mRNA is an **AL proxy** |
| **sense_lncRNA_host** | 296 | 41% | host is a lncRNA — often the pri-miRNA's OWN unit (MIR17HG, DLEU2) → proxy, weak confound |
| **antisense_overlap** | 40 | 6% | opposite strand → **NOT** co-transcribed (autonomous-like) — strand matters |
| **intergenic** | 41 | 6% | autonomous, own promoter |

Flags: `clustered` (polycistronic), `mirgenedb` (MirGeneDB bona-fide quality), `n_loci`/`context_mixed` (multi-locus).
**Multi-locus** (72 arms; 20 context-mixed): every GENCODE copy classified; representative = majority, ties→coding-host
(abundance is the SUM over copies, so one host-coupled copy already injects the host program); dose-weighted host-coupled
fraction is a refinement. Validated on canonicals: miR-93→MCM7, miR-21→VMP1, miR-15a→DLEU2, miR-199a→DNM antisense.

## 3. Arm→locus provenance — CANONICALIZED (see DATA_SOURCES for full detail)
- **GENCODE `gene_type==miRNA` (1879 loci, strand) = the SOLE coordinate authority** (100% HE coverage; `_arm_gencode_loci`
  handles all copy conventions). **MirGeneDB (`cnv_miRNA.csv`, 506) is a QUALITY FLAG only** (`mirgenedb`), not a
  coordinate source. Same coordinate system as the host genes → no cross-source mismatch.
- **Universe-redefinition (GLOBAL, `.N` fix LIVE):** the arm-name `.N`/case/suffixless-guide normalization
  (`data._arm_name_map`) recovered **~18 HE arms (16 expressed) + ~300 edges** (incl. canonical miR-101/124/126-3p) ⇒
  corrected universe **730 arms / 494 expressed / 5799 edges**. Persisted outputs stay old-universe until a refresh.
  Edge-DB hygiene: unicode-dash / comma-joined rows cleaned at source (`ledger._clean_edge_arms`). No mouse leakage
  (all `hsa-`; 1 cosmetic mouse-style `miR-21a`). (Memory `universe-redefinition-pending-refresh`.)

## 4. Host-target relatedness (`host_target_relatedness`) — the gate
Host pleiotropy/confounding is real ONLY where the host is related to the target. Partial Spearman(host, target | CPE,
**miRNA**) over the 1928 coding-host edges: **21% (403) are host-related** ⇒ genuine risk; the rest = benign shared locus.
**The miRNA is partialled OUT (upgrade 2026-07-11, user-driven):** the gate must isolate the host's OWN target-association
(host protein → T), NOT the host merely PROXYING the miR's transcription (sense co-transcription, ρ≈0.29) — since the miR
is the very thing we study in both consumers, removing it sharpens rather than weakens. **Empirically it does NOT lose the
host signal** — mean |ρ| 0.184→0.190, net MORE inclusive (372→403 related): it reshuffles, DEMOTING co-transcription-
entangled edges (E2F5/BTG2@MCM7 drop <0.3) and REVEALING edges where the miR masked the host↔target link. Proliferation is
kept IN (it mediates for prolif hosts). Columns `r_host_target` (CPE-only, reference) + `r_host_target_mir` (the gate).
Top-risk arms: miR-26a→CTDSP2, miR-32→TMEM245, miR-197→GNAI3. Output `output/learned/host_target_relatedness.tsv`.

## 5. DOWNSTREAM DEPENDENCIES — (a) + (b) + (c) ALL BUILT 2026-07-11
The axis is upstream of three consumers; **(b) is the unification**, **(a) is the CN-exclusion wiring**, **(c) is the family mix**:

**(a) CN instrument / locus — ✅ BUILT (`instrument.host_pleiotropy`, 2026-07-11, axis BP).** For a coding-host arm's edge
where the host relates to the target, condition the exclusion `δ_s` on the **specific host gene** (host mRNA) and emit a
`pleio_down_weight`∈[0,1] that **inflates s²_π** (weaker CN pull, NOT a subtraction) — the relatedness-gated form of axis
**BP**. Host from annotation (sense_coding_host), gated by `host_target_relatedness` (partial ρ≥0.3); same
`condition-δ_s` + **ACT-ONLY-ON-REDUCTIONS** core as `coding_pleiotropy` (MH-99). **ONE JOINT GATE:** the host is FOLDED
INTO `coding_pleiotropy` as a MANDATORY GUARD-EXEMPT candidate (`host_genes=`/`scan=`), so `exclusion(coding=True, host=True)`
conditions {coding ∪ host} JOINTLY → one `pleio_down_weight`, RETIRING the earlier separate-gates-via-**MAX** (MAX assumed
same-slice; co-amplified confounders are collinear ⇒ the joint **de-double-counts**: ESR1 0.50 vs MAX 0.72, PIK3CG 0.67 vs
0.90). Guard-exempt = the scan's prolif-survival guard can drop MCM7 yet the host keeps it (E2F3/RBL1 = `MCM7*` alone).
**COMPLEMENTARY, not redundant (measured):** of 29 scoped host-confounds the generic scan MISSES 15 and lists the host
itself in only 5 — starkest for the MCM7 prolif hosts (⚠ cross-reference `prolif_verdict`, don't blind-apply). **δ_s-SIGN
CONFIDENCE (`pleio_confidence`):** a residual-miR leak has sign −sign(γ), so δ_s SAME-sign as γ **cannot** be lost-miR ⇒
`clean_antirepression` (unambiguous); opposite = `repression_direction` (residual-miR OR a co-targeting other family →
grade via miR-partial + `prolif_verdict`). **Batch (`run_exclusion(coding=True, host=True)`, full HE universe, 8-worker
BLAS-pinned + OPTIMIZED loop — precomputed lstsq projectors, vectorised co-amp, scan restricted to host families ⇒
now ~2.4 min after the hot-path caching speedup): base table now 45,187 (was 45,156; the +31 is the `.N` universe refresh
persisting, verified via a base-only run — NOT a regression); per-locus 36 host + 11 pure-coding well-instrumented (F>10)
down-weighted**: IL6/miR-26a@CTDSP2* 0.87, IGF1/miR-28@LPP* 0.84, HMGCS2/miR-425@DALRD3* 0.72, PIK3CG/miR-7@HNRNPK* 0.67,
CCNE1/miR-16@SMC4* 0.51, ESR1/miR-26a@CTDSP2* 0.50, FBP1/miR-21-3p@VMP1* 0.43, MCM7 cluster 0.30–0.36. **✅ MULTI-LOCUS
caveat RESOLVED (per-locus host, 2026-07-11):** `genomic_context.locus_host_map()` (→ `locus_context.tsv`, 506 CN loci)
classifies each CN locus at its OWN hg38 coordinate (bridge: `_precursor_loci` MI*→coords → `_classify_locus`, 506/506), and
`exclusion`/`host_pleiotropy`/`dose_attribute_host_downweights`/`host_target_relatedness` now condition each **ACTIVE** locus
on **its own** host. A coding host at a SILENT (non-instrumented) locus is no longer applied to the active instrument — miR-30c
→NFYC (active = LINC00472) and miR-194→BPNT1 (active = MIR194-2HG) correctly DROP (37→36 host down-weights); the 6 multi-coding-
host arms (miR-26a→CTDSP2+CTDSPL, miR-218→SLIT2+SLIT3, …) now score both hosts (`host_target_relatedness` 403→444, +`host_locus`).
Cols `pleio_down_weight`/`s2_pi_pleio`/`beta_weight_pleio`/`pleio_confidence`
/`pleio_sources` (host tagged `*`); `host_pleiotropy` remains the standalone host-only diagnostic. Also HARMONIZED the coding
scan's candidate filter to partial Spearman|C ≥0.3 (aligns cosmetics; conditioning sets stay different on purpose — coding
removes prolif, host keeps it). Extended CN matrix `output/matrices/mirna_locus_cn_gencode.tsv.gz` (1862 loci × 1060).

**(b) Coupling confounders — the generalization of MH-100 → ✅ BUILT (`host_confound_lens`).** For a coding-host arm
(host H)→T where H relates to T, H's program confounds the miR→T coupling (the miR abundance is host-transcription-driven).
**This is MH-100 generalized:** proliferation was the special case where H happened to be a proliferation gene (MCM7/SKA2).
So **the host gene is the a-priori, per-edge confounder that MH-99 (coding-gene gate) and MH-100 (proliferation) were both
approximating.** `host_confound_lens` runs `prolif_verdict`'s non-circular OOF 2×2 but conditions on the **specific host
gene** — the host axis `resid(H|C)` is added **IN ADDITION to** `mal_prolif` (kept in C) and orthogonalized to it, so it
tests the host confound *beyond* proliferation. The OOF coupling is **repression-directed −ρ, NOT |ρ|** (user-caught
fix, 2026-07-11): abs() falsely scored sign-FLIPS to anti-repression (ρ −0.04→+0.16) as "confound-strengthening"; −ρ
correctly routes those to over_control/fragile (miR-23b→PLAU confound→over_control; miR-224→KLK10→fragile) while genuine
repression-strengthening (miR-186→FOXO1 −0.35→−0.47) stays confound. **The same abs() flaw in MH-100's `prolif_verdict`
was FIXED + all three outputs REGENERATED under −ρ (2026-07-11; headline robust — GW over_control 389 vs confound 170,
panel confound 5→3).**
**RESULT (403 related coding-host edges under the miR-partial gate, signed −ρ): confound 77 · over_control 130 · fragile 32 · neutral 164; mean dC
−0.032** (net over-control, like MH-100 — mostly the host is the mechanism, per-gene mixed). **The value-add: 60 of the
77 confound edges (78%) are NOVEL non-proliferation hosts** (|ρ(host,prolif)|<0.4;
`host_prolif_corr`/`host_novel` in the output) — host confounds on a different axis than proliferation, INVISIBLE to
MH-100: miR-874→STAT3 (host KLHL3), miR-28→FOXO1 (LPP), miR-204→BIRC2 (TRPM3), miR-23b→PLAU (AOPEP). So the lens is
correctness/attribution (gene-specific, per-edge FLAG), not a global lever. Output `output/learned/host_confound_lens.tsv`.
(NB the per-edge `edge_verdict` column is noisier — tiny-β edges blow up the rel ratio; the gene-level OOF verdict is the robust signal.)
- **MECHANISM (antisense control, resolved 2026-07-11) — sense co-transcription IS real and STRAND-SPECIFIC.** Population-
  wide on EXPRESSED arms: **sense host↔miRNA mean ρ +0.30** (raw; +0.21 partial\|CPE; 25% ≥0.4 — miR-452/GABRE 0.81,
  miR-10a/HOXB3 0.73, miR-139/PDE2A 0.72, miR-93/MCM7 0.56) **vs antisense +0.05** (partial −0.04; only miR-199b/DNM1 0.37,
  the rest ~0). So sense-intronic co-transcription drives host↔miRNA; antisense neighbours essentially do NOT co-express
  (rare chromatin exceptions). *(An earlier draft wrongly concluded "local co-expression, not co-transcription" from the
  SELECTED 16-edge related-antisense subset — a selection artifact; the population/expression-filtered comparison above
  overturns it.)* ⇒ the confound is genuinely **sense-co-transcription-driven** (host supplies the miRNA); the strand
  distinction holds; antisense-overlap stays a distinct (autonomous-like) class. Caveat: host↔miRNA is 0.3 not ~1
  (splicing↔Drosha processing divergence), so the AL host-mRNA proxy is a transcription-rate proxy with processing noise.
- **✅ FOLDED INTO THE ANNOTATION (2026-07-11): `promoter` + `promoter_class` columns in `genomic_context.tsv`** — every HE
  arm now carries its FANTOM5 CAGE co-transcription call, queryable directly ("is miR-X co-transcribed?"). `promoter_class`
  ∈ {**host_shared** (promoter drives the host ⇒ co-transcribed), **shared_other** (a different gene's promoter), **own**
  (coordinate-only / own `MIR*` promoter ⇒ independent), **unknown** (absent from FANTOM5)}. Distribution: sense_coding_host
  **271 host_shared (80%)** / 33 shared_other / 32 own / 1 unknown; antisense & intergenic ~all `own` (correct — not
  host-co-transcribed); the `own`/`shared_other` sense-host arms are the **positionally-intronic-but-transcriptionally-
  INDEPENDENT** cases. **Classification = de Rie's LABEL is the discriminator + a TWO-PART ACCESSION RESOLVER** (the label
  can name the gene by a transcript accession instead of the symbol; both handled): (a) **`pN@ENST…` → gene via the GENCODE
  transcript map** (exact — let-7a→ENST00000533109→MIR100HG); (b) a **non-ENST / retired accession** (GenBank `AK…`,
  ENST-not-in-v49) is confirmed as the host when its FANTOM5 hg19 TSS **lifted to hg38** (`hg19ToHg38.over.chain`, pyliftover)
  lands **inside the host body** — applied ONLY to accession labels (miR-21→p1@AK310806=VMP1 ✓, miR-101→JAK1 ✓), NOT to
  coordinate labels (an intronic-OWN peak is also within the body → stays `own`). Recovered 9 host_shared that a label-only
  match under-called. (FANTOM5 hg19 vs GENCODE hg38, so a naïve TSS↔span match needs the liftover; the ENST map is exact.) (`_classify_promoter`.)
- **✅ INDEPENDENT COORDINATE annotation (2026-07-11, `_promoter_coord`/`_locus_at`) — a label-free cross-check + full detail.**
  Four more columns derived PURELY by lifting the FANTOM5 hg19 promoter TSS → hg38 and reading GENCODE there (no de Rie label):
  `promoter_gene` (which gene the TSS sits in, host-preferring for nested genes e.g. MCM7&AP4M1), `promoter_gene_type`
  (protein_coding / lncRNA / …), `promoter_feature` (**exon / intron / intergenic** — mostly exon 544, since a promoter TSS is
  in the driven gene's 1st exon; **47 intron** = independent-intronic promoters), `promoter_coord_class` (host_shared /
  shared_other / intergenic). **Label (transcriptional) vs coordinate (positional) agree 78%; the disagreements are the signal:
  de Rie's own `own` LABEL (its CAGE call) is what marks independence — the coordinate only ADDS position: 23 arms are `own`
  by label yet `host_shared` by coordinate = a miRNA sitting physically INSIDE the host (often an intron)** (miR-1205 in PVT1's intron, let-7e in SPACA6) — the "positionally-intronic-but-transcriptionally-
  INDEPENDENT" population flagged per-arm. (`promoter_feature=intron` and the label↔coord mismatch both mark it.) Intron NUMBER is a further refinement not built. **Also `promoter_host_tss_kb`** = distance from the (lifted) miRNA promoter TSS
  to the HOST's EXACT hg38 TSS (`_host_tss`, no liftover) — the direct host-promoter comparison: host_shared-label arms median **0.1 kb**
  (fire from the host promoter) vs own-label **2.3 kb**; it also reveals alt-promoter co-transcription (miR-21 host_shared but 130 kb from
  VMP1's canonical TSS = alternative VMP1 promoter). Liftover here is artifact-free (0% unmapped/multi-map on all 2,248 TSS); de Rie
  promoters are CAGE-**experimental**, not distance. So the table carries the de-Rie transcriptional call + an independent positional one + the host-TSS distance.
- **MECHANISTIC DIVE (2026-07-11) — DEFINITIVE via FANTOM5 promoter data (`data/external_cache/fantom5_human_miRNA_promoters.tsv`,
  de Rie 2017 CAGE).** The `promoter` field = the promoter each miRNA ACTUALLY uses (`p1@HOSTGENE` = host promoter). Of 211
  sense_coding_host arms: **74% HOST-SHARED promoter (co-transcribed, A) · 26% OWN promoter (independent, C).** So
  co-transcription is COMMON, not a minority. **KEY — mature-miRNA↔host-mRNA CORRELATION does NOT reflect promoter
  architecture** (host-shared mean ρ 0.29 ≈ own-promoter 0.28): POST-TRANSCRIPTIONAL PROCESSING (Drosha efficiency, LIN28,
  stability) decouples the mature level even when co-transcribed — e.g. **let-7f/HUWE1 (host-shared but ρ −0.09, LIN28
  processing block), miR-126/EGFL7 (host-shared, ρ 0.12)**; while miR-93/MCM7, miR-452/GABRE, miR-139/PDE2A are host-shared
  AND tightly coupled (0.56–0.81). **⇒ correlation is a POOR proxy for co-transcription (it measures processing coupling);
  an earlier correlation-based "31% independent" was mostly co-transcribed-but-processing-decoupled, NOT own-promoter —
  true own-promoter is 26%.** Implication: the host-confound premise (host drives the miRNA transcriptionally) holds for
  ~74%, BUT the host mRNA is a WEAK proxy for the mature miRNA (processing) — so conditioning on host mRNA captures only
  the transcriptional slice. The verdict still tracks coupling STRENGTH (strong→over_control=host supplies dose;
  weak→confound), but "strength" = processing-coupling, not promoter class. (miR-199b/DNM1 antisense = neighbourhood chromatin, distinct.)
- **RESOLUTION — the co-transcription framing is neither necessary nor sufficient; the lens is a host↔TARGET co-expression
  detector.** Partial ρ(arm,host|CPE,prolif) does NOT distinguish promoter class (host-shared 0.29 vs own 0.31); and
  **re-gating the lens on FANTOM5 promoter class does NOT separate the confound** (confound 16% for BOTH host-shared and
  own-promoter; over_control 33% vs 45%). So the host-program confound is driven by **host↔target co-expression**, with the
  host↔miRNA co-transcription a real-but-weak (processing-decoupled) SECONDARY factor, not the gate. Honest axis scope:
  *"an intragenic miRNA whose host co-expresses with the target introduces a host-program confound"* — independent of
  whether it is strictly co-transcribed. (Whole thread 2026-07-11, user-driven rigor: sign-fix, antisense control, independent
  hypothesis, FANTOM5 promoter data, partial-corr — each tightened the claim.)
- **What the novel host confounds ARE (diagnosis):** host↔target co-expression **robust to CPE+proliferation** (raw ρ 0.38
  → PARTIAL ρ(H,T|CPE,prolif) **0.40, all 46/46 ≥0.3, 45/46 POSITIVE**) but **NOT shared programs** (only 1/46 share a
  Hallmark) and **NO direct curated regulatory edge** — neither transcriptional (CollecTRI 0/46) nor signalling (full OmniPath
  fetched, 0/46; host node coverage 17/28, of which none regulate their target) ⇒ **host-gene-specific SHARED-UPSTREAM
  co-regulation** confounds: the host is a co-expressed **PROXY for a common upstream driver** of both host and target
  (metabolic PPARGC1B/CDK6; SMAD-phosphatase CTDSP2/ACVR1; RBP HNRNPK/IGF1R; lipid OSBP/CTNNB1), NOT a direct H→T
  regulator. Conditioning on the host still de-confounds (it proxies the shared driver). Structurally invisible to
  program-level controls (proliferation / C block) and to curated networks — which is why the specific host is required.

**(c) Multi-arm/family context — ✅ BUILT (`genomic_context.family_context` → `output/learned/family_context.tsv`, 2026-07-11).**
A seed family collapses arms of possibly-different contexts (some coding-host, some intergenic/antisense; some co-transcribed,
some independent) into ONE pooled `X_fam`, so the family's host-program confound / CN-instrument validity is a DOSE-WEIGHTED
MIX. `family_context` crosses the per-arm annotation with the family structure + member abundance: per family it emits the
member context **composition**, **dose-weighted coding-host & co-transcribed fractions** (Σ RPM-share over the sense_coding_host
/ host_shared members — the fraction of the family DOSE carrying a host confound), the **dominant (dose-max) arm** + its context,
and the distinct coding hosts. **573 families (85 multi-arm); 47/85 multi-arm families are context-MIXED; 53 families are
dose-co-transcribed≥0.5.** **KEY FINDING — the confound follows the DOSE-DOMINANT arm, not the composition:** miR-17~92
(17/20/**93**/106) is **88% co-transcribed** (dominant miR-93-5p@MCM7) and miR-15/16/195/424/497 is **71%** (miR-16-5p@SMC4) —
host-confounded — while **let-7/98 has FOUR coding-host members yet is only 11%** (dose-dominant let-7a-5p is lncRNA-hosted) and
miR-30-5p is 11% — **autonomous despite coding-host members** (they are dose-minor). This GENERALISES MH-100's arm-source result
("the miR-17~92 confound localises to miR-93-5p") to the whole family universe, derived purely from annotation + abundance; it
is **now WIRED into (a)'s output** (`instrument._attach_family_context` in `run_exclusion`): every CN-admission row carries `fam_dose_cotx_frac`
/ `fam_dominant_arm` / `fam_dominant_host` — arm-source attribution + a co-transcription-vs-independent QC on each host down-weight (e.g. SLC7A5/miR-744
tagged `fam_dose_cotx_frac`=0.0 = positionally-host-but-independent). INTERPRETATION layer (annotation), not an estimator change. **But `family_context` SHOULD feed the estimator — this is the
genuine open design, not a non-issue.** δ-pooling (`within_family.delta_pooling`, **MH-98**, BUILT) estimates WHICH member
carries the family dose (fusing mRNA + per-segment CN + chimeric into `M_member = M_fam + δ_m`); `family_context` annotates
the SAME dose-dominant arm from abundance, plus whether that arm is a coding-host. The intended pipeline (the reason (c) was
built — generalizing MH-100's "the confound follows the ONE dose-dominant arm") is the **dose-GATED host conditioning**:
down-weight only the member whose dose IS the coding-host, using δ-pooling's fused dose-share as the estimator-grade upgrade
of the abundance dominant-arm (= the `family_context` estimator-level use, §"Forward" item). They are the coarse map and its
estimator-grade refinement of one quantity, not inert parallels — currently un-wired (open), not "shouldn't meet." (Both the
CN between-family channel — `channel_cn`, MH-93/94 — and δ-pooling are already BUILT; see `CN_INSTRUMENT.md` §Bayes-unification.) Also the annotation is now COMPREHENSIVE: `host_region` (which NUMBERED intron/exon of the host the miRNA sits in — 479
intronic / 192 exonic) + `promoter_feature` numbered (e.g. miR-93 in MCM7 intron 13, promoter at MCM7 exon 1). **This closes the axis's
downstream — (a) + (b) + (c) all BUILT.**

## 6. Code + outputs
`learned/genomic_context.py` (`classify_he_arms`, `host_target_relatedness`); `data._arm_name_map` (universe fix);
`ledger._clean_edge_arms` (edge hygiene); `scripts/analysis/build_mirna_locus_cn.py`. Outputs:
`output/learned/{genomic_context.tsv, host_target_relatedness.tsv}`, `output/matrices/mirna_locus_cn_gencode.tsv.gz`.
