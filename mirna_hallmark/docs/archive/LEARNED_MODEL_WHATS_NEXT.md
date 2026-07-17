# Learned model — WHAT'S NEXT / TODO (living list, started 2026-07-05)

Forward list for the `mirna_hallmark/learned/` model. Companion to `LEARNED_MODEL_DESIGN_RESPONSE.md`
(design), `LEARNED_MODEL_BUILD_PLAN.md` (what's built), `LEARNED_MODEL_ESTIMATOR_MAP.md` (**which estimator
does which job & why over alternatives** — read this to orient), `LEARNED_MODEL_METHODS.md` (**formula-level
spec** of every estimator: objective + components + code ref), `LEARNED_MODEL_VALIDATION.md` (**the evidence
dossier**: gate numbers + positive controls + sharpening/novelty examples + honest negatives), `ATTRIBUTION_CONTEXT_AXIS.md` (the card
design), `EDGE_QUESTION_TAXONOMY.md`/`GENE_QUESTION_TAXONOMY.md` (E/G estimands). Status tags: ✅done · 🔨now · ⬜next · 🔬investigate.

## Where we are (2026-07-05)
The pipeline is built and confounder-honest: gene-focused NN adaptive-lasso; PMID-deduped method-centric
evidence ledger (transfection-calibrated); seed-family pooling; CN instrument (causal); CPTAC protein OOD;
genome-wide FDR-controlled **discovery** (+ deconv composition tag); program-wise M on 5 Hallmarks;
hierarchical Bayes (uncertainty); cross-state (NAT/GTEx) coupling, within-patient **realization** (E5),
within-gene **budget-shift** (E7/G4, GTEx→NAT→tumour). **Confounders (design-compliant, ablated):**
core C = CPE + target-CN + malignant-proliferation (in-fit); CIBERSORTx composition in-fit
program-conditional (+ retention tag); batch post-fit-only (≤0.03, negligible); TF-proxy diagnostic
(over-controls); HRD off. Everything re-run under corrected C (coupling strengthened; discovery 3034
candidates / 568 novel).

## ✅ OBSERVATION-MODEL DISTRIBUTION — RESOLVED: Student-t, not NB (MH-92, 2026-07-10)
**Ran the dense-mRNA residual diagnostic (the flagged concrete first step) → the dense §6b Gaussian-on-log-CPM likelihood
is meaningfully HEAVY-TAILED, and the fix is Student-t, not NB.** Diagnostic (85 genes × 1041 patients, exact `_gibbs_posterior`
π≡1 fit): pooled excess-kurtosis **3.5** (median **1.2** across the 82 well-expressed genes), Student-t MLE **df≈7**, tails
P(|z|>4)=**33×** Gaussian. Two mechanisms: **amplified-subset** tumours (+skew, ERBB2/ERBB3 — a patient mixture, axis Y3) and
**near-floor** genes (−skew, only 3, all essentially unexpressed). **NB rejected:** no low-expression count gradient
(Spearman(expr, kurtosis)=**+0.19**, *wrong* sign); the near-floor genes are a detectability-floor issue, not count structure.
- ✅ **BUILT — Student-t in `spike_slab._gibbs_posterior(nu=)`** via scale-mixture of normals (`ε_i~N(0,σ²/λ_i), λ_i~Gamma(ν/2,ν/2)`;
  per-obs λ_i Gibbs step) — **stays Gibbs, no HMC**; `nu=None` bit-identical to §6b. VALIDATED: nests Gaussian (Δβ≈0.001),
  **robust 1.76×** under 5% gross contamination, but **OOF coupling Δρ=+0.001 (n.s.)** and attribution unchanged (corr 0.999,
  5/5 drivers) ⇒ a **correctness/robustness upgrade, not a performance lever** at n≈1000 (the likelihood dominates — same shape
  as the CN-channel "real-but-weak" verdict). METHODS §3b (twins); CHANNEL_FUSION §AB; VALIDATION.
- ✅ **CN channel t** wired (`t['nu']` on the `π̂_s` pseudo-obs) — discounts a collider-outlier segment (β 0.535→0.508), consistent
  segment intact. ✅ **Exclusion-gate SE bootstrap** (`instrument.exclusion(n_boot=)`) — nonparametric bootstrap of the product
  `a·b`; **confirms the analytic Sobel SEs are adequate** (s²_π agrees 0.94–1.11×) because s²_π is **δ²-pleiotropy-floor-dominated**
  on the weak/pleiotropic segments (the δ² term, identical in both, swamps the delta-method sampling term).
- ⬜ **REMAINING per-channel (future channels only):** **protein** = log-normal + non-linear discordance link; **binding** =
  Poisson/NB (counts); **methylation** = Beta (bounded) — these break conjugacy → HMC (axis J), decide when each channel lands.
Raised 2026-07-10 (user: "very important, name it"); resolved same day.

## ✅ K_D SPECIFICITY — genome-wide, DISCARD OVERTURNED (MH-87, 2026-07-10) + 🔬 forward
**GENERAL identity/discovery infrastructure** — consumed by `discovery.py`, the selection prior (`priors.inclusion_prior`),
`attribution_identity`, AND (as the affinity leg) the CN instrument's Ring-1 (`CN_INSTRUMENT.md §7`) — but **not**
CN-specific. Lives here, not in the CN doc.
- ✅ **Discard OVERTURNED.** MH-86's "context++ > K_D" was an **HE-universe artifact**. All-site **genome-wide** scanMiR
  K_D (`kd.genome_affinity`, 746 detectable arms × 18,852 MANE genes; built by `scanmir_genomewide_par.R`, BiocParallel)
  recovers canonical specialists **0.89 vs 0.79** fair context++-genome (12 genes, `kd_fair_bench.tsv`); **strong-site
  thresholding HURTS**. It's **per-arm** — the resolution context++ (per-seed-family) structurally can't reach. Wired:
  `discovery.py` (biochemical support), `seq_specificity.affinity_percentile_kd`, `attribution_identity(spec_source="kd")`.
- 🔬 **context++ × K_D — revisit COMBINING.** context++ and K_D are **nearly orthogonal** — Spearman **+0.065** (78k
  (family,gene) pairs): K_D = scanMiR RBNS **binding**; context++ = TargetScan **repression-outcome + conservation +
  site-context**. Orthogonal predictors are exactly where a **fusion beats either** — MH-86's "don't combine" was
  HE-conditioned *and* on the worse (HE) K_D. Test a context++ × K_D fusion for specificity + discovery candidacy.
- ⬜ **Scan scope + floor sweep.** Current genome-wide K_D = **746 detectable arms** (RPM≥10) — the principled
  *discovery* universe (excludes unexpressed = the rarity-failure noise); per-arm *specificity* is universe-independent
  (self-contained denominator). Optional: (1) a **complete ~2,600-arm reference scan** (one `scanmir_genomewide_par.R`
  run) as a reusable `(arm,gene)` K_D resource; (2) a **discovery floor sensitivity sweep** (RPM ≥ 50/10/1 — a
  re-*filter* of candidates, **not** a re-scan) — discovery FDR is already null-controlled, so the floor is the
  sensitivity/precision lever, not an FDR lever.

## 🔨 IMMEDIATE — the per-edge (+ per-gene) attribution card
Generalises the standalone range-stat request. One row per edge, joining what we already compute:
1. **regime context (range stats)** — arm median RPM, % samples above the functional floor (RPM≥10),
   IQR, **spiker flag** (low median + high IQR = subset-driven, e.g. miR-135b 28%>floor); target median +
   IQR (a strong coupling on a low-IQR gene = miRNA explains most of its variance).
2. **budget share (E7/G4)** — `states.budget_shift`: rank + share of the gene's pressure, GTEx→NAT→tumour Δ.
3. **composition tag** — deconv retention: cell-intrinsic (≥0.7) / partial / composition-explained (<0.4).
4. **shift-class** — the precursor `mirna_state_class.joint_edge_class` (level dHN/dNT/dHT × realization)
   EXTENDED with our paired-E5 realization + low-variance/undetectable split. Build on the cell-intrinsic M.
Then the **per-gene card** = the G-series roll-up (net-repression G1, coherence+role G8, budget concentration
G4, composition fraction, shift-class G10).

## 🔬 INVESTIGATE — the top discoveries (deep literature + data)
This is a first-class task, not a footnote — the discovery lane's value is validated novel biology.
- **The miR-17~92 / miR-106b-25 target cluster** (miR-106b/93/17/19a → RABEP1, AHNAK, TGFBRAP1, IL6ST,
  AFF1, WWP1, M6PR): ubiquitously expressed arms (100%>floor), cell-intrinsic, strong, weakly-curated.
  → lit scan each (known target? breast context? tumour-suppressor targets AHNAK/WWP1?); data: CN-instrument
  causality, CPTAC protein, Manakov binding, within-patient realization; is this a coordinated cluster→program hit?
- **miR-135b-5p → GATA3** (fully novel, validated: CPTAC protein −0.49 + Cancer Cell Int 2021 ceRNA axis;
  a **spiker**, 28%>floor → coupling in the miR-135b-high subset). → deepen: which tumours (TNBC?), functional
  consequence (GATA3 loss → dedifferentiation), the circ_0044234 axis; is the subset a clinical stratum?
- **The 568 fully-novel candidates** → triage (cell-intrinsic + spiker/expressed + biochemically-supported),
  lit-novelty check, rank a shortlist for Bar-5 protein + Manakov validation.

## ⬜ EXPAND — new novel interactions (scale the discovery+validation lane)
- CLIP-source candidate expansion (ENCORI/POSTAR3 physical sites) beyond TargetScan nomination.
- Systematic Bar-5 (CPTAC protein) + Manakov (chimeric duplex) validation of the shortlist.
- Independent-cohort replication (METABRIC / Buffa) of the learned M / top edges.
- **State-stratified discovery (NAT / GTEx)** — run the *discovery lane per state* (its own permutation-FDR
  null + scanMiR biochemical + deconv composition, per state), then compare which edges are found where
  (found-in-tumour-not-NAT = candidate *acquired* regulator; found-in-all = constitutive; found-only-NAT =
  field/lost). This is the estimator-honest way to get **cross-state discovery** — it handles the *selection*
  instability via the FDR null (unlike a naïve wide per-state lasso). **DISTINCT from the wiring axis**: wiring
  (ΔM) uses a *fixed* support to *compare* known edges; this *discovers* new state-specific edges. **Power-
  bounded**, not meaningless: NAT n=104 / GTEx n=327 ≪ tumour ~1000 → fewer discoveries; and *undetectability*
  is worst in healthy (can't discover an edge where the miRNA/gene doesn't vary — miR-200/ZEB1) → "not found
  in healthy" ≠ "not a healthy edge". NAT (more variance than GTEx) is the more feasible end.

## ⬜ EPIGENETIC-REWIRING DISCOVERY — lost & gained specialists (the bidirectional methylation lane, planned 2026-07-06)
The silenced-specialist detector is built (`structural_identity` + `methylation`, both validated). Scale it to a
genome-wide **discovery lane** on the miRNA→gene edge, with the methylation axis read in **both directions** (the
gate is symmetric — validated: miR-124/9/129/137 hyper→LOST; miR-21/155 hypo→GAINED):

- **Axis A — LOST specialist (de-repression):** structural specialist for g, **baseline-active** (GTEx/NAT > floor),
  **tumour-silenced**, promoter **HYPER-methylated** (Δβ≥+0.15) → g's repressor is switched off → **g de-repressed
  (mRNA/protein ↑)**. Oncogenic-consequence edge (classic TSG-miRNA loss).
- **Axis B — GAINED specialist (acquired repression):** structural specialist for g, **baseline-silent**,
  **tumour-active**, promoter **HYPO-methylated** (Δβ≤−0.15) → repressor switched on → **g acquires repression
  (mRNA ↓)**. Tumour-suppressor-silencing-via-miRNA edge (epigenetically-activated oncomiR).

**Method (per axis, genome-wide over Hallmark genes × HE seed-families):**
1. **Nominate (prior, ⊥ coupling):** `structural_identity(g)` → families with `structural_share ≥ 0.10` ∧ `confidence ≥ 0.5`.
   *Potential×specificity* is the designed-owner score (potential = per-edge repression magnitude "how strong f on g";
   specificity = per-source concentration "fraction of f's targeting aimed at g"; product = strong AND concentrated =
   the specialist — see ATTRIBUTION_IDENTITY_VS_MAGNITUDE §0). Expression-free, so it survives silencing — that's why
   it, not M/budget/Shapley, seeds this lane.
2. **Expression switch:** family arm level in GTEx/NAT vs tumour state matrices → baseline-active∧tumour-silenced (A)
   / baseline-silent∧tumour-active (B), by the functional floor log2(11).
3. **Methylation mechanism:** `locus_methylation(arms).direction` = hyper (A) / hypo (B), carry Δβ + CGI context.
4. **Target consequence (the payoff, the ONE inferential claim):** does g move as predicted? (i) group contrast —
   g mRNA tumour-vs-normal (A: up, B: down); (ii) within-tumour, g mRNA vs promoter β across cases (A: g↑ with β;
   B: g↓ with β) — this closes methylation→miRNA→target into one coherent edge. Permutation/label-shuffle null on
   (ii); **BY-FDR** across the candidate set on the joint consequence.
5. **Rank** by (structural_share · |Δβ| · target-consequence effect), report q.

**Confounders (must condition, else Axis B especially is a global-hypomethylation artifact):**
- **Purity/composition** — β is bulk; tumour purity & epithelial fraction shift β globally. Condition on CPE + deconv
  fractions; require a *promoter-specific* Δβ that survives purity adjustment. Global tumour hypomethylation is the
  known trap for Axis B → demand CGI-promoter-local, not open-sea, signal.
- **CN vs methylation** — miRNA-locus deletion also silences (Axis A) / amplification activates (Axis B). Cross the
  `mirna_locus_cnv` locus CN: a LOST specialist that is hyper-methylated **and** CN-neutral = clean epigenetic; drop
  CN-confounded calls or report the mechanism split.
- **Imprinted / constitutively-methylated loci** (DLK1-DIO3 miR-127/379/370 cluster, miR-134…) — flag & exclude;
  their β is developmental, not tumour rewiring (miR-127 baseline β 0.92 in the control run = exactly this).

**Validation:** positive controls already pass (above); literature cross-check top hits; **CPTAC protein** (Axis A: g
protein ↑ — de-repression is often translational; Bar-5); Manakov chimeric-duplex for the physical edge.

**Prereqs / outputs / power:** needs the **UUID→TCGA-barcode map** for the within-tumour β↔mRNA leg (the group Δβ
needs none — 798 tumour / 97 normal, well-powered). Output `output/learned/epigenetic_rewiring_discovery.parquet`
(one row per axis×gene×family: structural_share, confidence, baseline/tumour arm level, Δβ, direction, CGI, locus_CN,
g mRNA shift, β↔mRNA ρ, q) + a ranked shortlist per axis. Power-bounded (normal n=97) but the tumour arm is deep.

## ⬜ PARALLEL VIEW & REGULATION-PATTERN CHARACTERISATION (from rank concordance, planned 2026-07-07)
`learned/parallel_view.py` is the design-model (learned M) vs abundance/rank-infra parallel view: (a)
`within_gene_ranks` — abundance-rank vs budget/identity/structural; (b) `change_trajectory` — trio
GTEx→NAT→tumour abundance-logFC vs wiring-ΔM (n-matched retention for the wiring verdict). Forward:

- **Pairwise-ranking gap as a per-gene regulation FINGERPRINT.** Each disagreement between two of the four
  rankings (abundance · budget · identity · structural) names a regulation MODE: budget≈abundance but
  identity≠budget → coupling owned by a non-dominant arm (quiet owner / collinearity); identity≠structural →
  the *realized* owner isn't the *designed* specialist (opportunistic vs designed); structural anti-abundance
  (PTEN ρ=−0.65) → the designed owner is silent (loss candidate). Tabulate the concordance matrix per gene →
  cluster genes into regulation ARCHETYPES (abundance-controlled / weight-controlled / designed-but-silenced /
  redundant-collinear). This is the systematic version of the ad-hoc PTEN/ESR1 reads.
- **Ground STRUCTURAL specificity in targetome RARITY (sequence, not just evidence).** Specificity is currently
  the evidence-mass fraction. Extend it to the *actual* uniqueness of the arm's SEED in the 3′UTR universe: how
  many other genes carry that arm's 7–8mer site genome-wide → a seed-rarity / effective-targetome-size statistic,
  tied back to the scanMiR K_D and TargetScan context++ reads (rare strong site = genuine specialist; ubiquitous
  seed = pseudo-specialist). Makes "designed to *specifically* repress g" a quantitative SEQUENCE claim. (This is
  also the principled replacement for "there is no biochemical specialist" — rarity is the biochemical axis that
  the flat per-arm promiscuity count missed.)
- **GLOBAL rank concordance** (not just within-gene) via `program_network` engagement (Σ_g M): global abundance
  rank vs total-regulatory-weight rank — the loud-but-unwired vs scarce-but-hub arms at genome scale.
- **TARGET-SIDE SNV mechanism leg (identity≫structural is not only "opportunistic").** Structural is REFERENCE-
  genome + curation based, so identity-high/structural-low splits into: (i) spurious (residual confound / shared
  credit — partly handled by C-residualisation + the collinear-pair flag), and (ii) REFERENCE-BLIND REAL — a
  coupling via a site the reference doesn't encode. The genetic case: a tumour **3′UTR SNV that CREATES a 7–8mer
  seed match** for a regulator arm = acquired target site → real coupling, zero reference structural. Symmetric:
  a **seed-DESTROYING SNV** removes a designed site → the target-side analog of the silenced specialist
  (structural high, identity/budget→0, but mechanism = site loss not miRNA silencing). Build a seed-SNV leg
  (APM has SNV/VEP data): scan each gene's 3′UTR tumour SNVs for gain/loss of its regulators' seed matches →
  a TARGET-side rewiring mechanism SIBLING to the miRNA-side methylation gate (§EPIGENETIC-REWIRING). Also in
  this family, reference-blind: A-to-I **editing** (creates/removes sites or edits the miRNA seed → retargeting)
  and **APA** (exposes/hides 3′UTR sites) — both make identity≫structural without being opportunistic. So the
  identity/structural quadrant is a MECHANISM router: loss side {methylation, seed-SNV-loss, CN-loss}; gain side
  {hypomethylation, seed-SNV-gain, editing/APA-exposed, uncurated-real, residual-confound}.
- **GTEx logFC gauge check**: `change_trajectory` now reports BOTH the naive RPM−TPM (`rawFC_GN`) and the
  QN-into-TCGA (`relFC_GN`) GTEx→NAT change — validate the miRNA RPM≈TPM assumption (near-uniform mature length,
  shared count-dominating arm sets) against the QN version to see where they diverge (housekeeping/spike arms).

## ✅ CANONICAL M — estimator unification (DONE 2026-07-05) + speed pass
**Estimator matched to estimand** (this is the resolution of "every fit?"):
- **Coefficient-reading fits (ATTRIBUTION)** → now the CANONICAL `states.canonical_M`: **bagged z-scored NNLS
  on the fixed HE seed-family support** (`states._bagged_nnls`), arm-broadcast (per-arm split = nomination).
  Verified reproducible (corr ~0.99) and recovers curated drivers (miR-21/103a/181/182/96 on PTEN). Wired into
  `states.budget_shift`, `card.decompose`, `card.state_wiring`, `card.stable_wiring`, and `card.gene_card`
  (budget via budget_shift). Added a **per-state variance floor** (std<0.1 arm → zero weight) so undetectable/
  near-constant arms in small states (NAT n=104) stop getting spurious M via z-scoring-by-~0 (miR-105-3p
  M_NAT 0.157→0). PTEN tumour drivers unchanged; only the small-state tail tightened.
- **Aggregate-prediction fits (VALIDATION/DISCOVERY)** = `mvp.oof_gate`, `discovery`, `_aggregate_coupling`,
  `realization` → **stay the adaptive lasso**. WHY (not inconsistency — estimator matched to estimand):
  (1) they read only the *aggregate* X·M (a prediction), estimator-robust under collinearity even when the
  coefficients aren't (verified ρ within 0.03, SD ~0.01) — swapping wouldn't change any verdict;
  (2) the gate's *purpose* is to certify the adaptive-lasso MODEL beats abundance / matches curated — it must
  fit the object it certifies; (3) the evidence PRIOR enters only through the lasso penalty, and DISCOVERY is
  a *selection* problem (which orphans enter, FDR-null) that NNLS-on-fixed-support structurally cannot do.
- **Speed pass (so the batch card is cheap):** cached the 913k-row edge TSV (`data_loaders.load_hallmark_edges`)
  → `assemble_gene` **742ms→82ms (9×)**, every hot path faster; `card.gene_card` assembles raw/deconv/GTEx
  **once** not per-arm → PTEN **20s→0.7s**. Batch = ~14s one-time warmup + ~0.7s/gene.
**Remaining:** one canonical rerun of the persisted attribution outputs (`output/learned/programs/*`, any cards)
now that the estimator is settled; and (optional) canonicalize the `shift_vs_weight` per-edge weight diagnostic.

## ⬜ REMAINING PHASES (from the design plan)
- Phase 4 cooperativity + ceRNA shared-pool (Decision G) — **caution**: sits on the saturation substrate
  that failed (occupancy); gate hard on held-out improvement.
- Decision I Shapley attribution (identity vs magnitude) — the formal per-edge share decomposition.
- State axis: abundance-vs-wiring decomposition (Δx·M vs x·ΔM); subtype interaction tests (Gelman–Stern).
- `priors.py` spike-and-slab π/μ/τ object.

## 🔬 VALIDATION GAPS (the real remaining validation)
- **CPTAC protein — SYSTEMATIC (Bar 5).** Only miR-135b→GATA3 validated at protein so far. The sharp test:
  does the mRNA-learned coupling / do the discovered edges hold at the **protein** level (functional readout,
  where mRNA-protein discordance lives)? Infra: `eval/cptac_validation.py`. **Highest-value item left.**
- **External-cohort replication (METABRIC / Buffa)** of the learned M / top edges. Buffa infra exists;
  METABRIC-miRNA is EGA-gated (partial subset only).

## ⬜ COMPOSITION / CONFOUNDER / INFRA TO-DOS
- **CIBERSORTx run — bundle PREPARED, awaiting external Docker run.** `output/cibersortx_transfer/`
  (`scripts/deconvolution/prepare_cibersortx_transfer.py`): state-matched refs — **HBCA** (healthy, level1+
  level2 = FB1-4/PV1-5/luminal) for **GTEx/NAT**, **Wu** (tumour, minor=myCAF/iCAF+cancer / major) for TCGA;
  + GTEx & TCGA linear-TPM mixtures. → run in Docker → drop `*_Fractions.txt` into `output/brca_deconvolution/`
  → wire `_DECONV_COLS` → **level1-vs-level2 over-control ablation** (42 fractions risk over-control; pick the
  coarsest sufficient). Caveat: no adipocytes in either atlas.
- **Proper batch** if ever needed: covariate-protected ComBat (preserve PAM50) — not naive dummies.
- Fold ENCORI/POSTAR3/Manakov into the ledger (union, not summed) for inclusion/discovery.
- `baselines/` re-export shims; `git add` the `learned/` tree (whole tree currently untracked).
- **DOC CONSOLIDATION (parked 2026-07-06, user deferred):** the learned model is production but its docs still
  read as a separate parked "cluster" with parallel/duplicate docs. Plan: (1) merge `LEARNED_MODEL_METHODS` →
  `FORMULAS.md` (one formula doc); (2) merge `MODULE_MAP` → `ANALYSES_CATALOG` (one component index); (3) archive
  `LEARNED_MODEL_{DESIGN_RESPONSE,DISCUSSION_PROMPT,BUILD_PLAN}` → `docs/archive/` (design-phase, now built);
  (4) promote `ESTIMATOR_MAP`/`VALIDATION`/`WHATS_NEXT` into the §1/§2 canonical set. Net ~13 model docs → ~6.
  Content already made current 2026-07-06 (pooled-HE, loss lens, un-STUB attribution); this is structure only.
  Leftover if NOT consolidating: one pooled-HE note in `MODULE_MAP`'s learned section + add structural_identity/
  methylation module rows.

## 📐 DESIGN vs AS-BUILT — the gaps, precisely (each defined so it's actionable)

> **⚠ PRE-CONVERGENCE (2026-07-05) — READ WITH THE §6b LENS.** This section is scored against the design-phase
> `LEARNED_MODEL_DESIGN_RESPONSE.md` decision log and predates the model converging to **one Bayesian posterior,
> lasso retired** ([`DISCOVERY_SYNTHESIS §6b`](LEARNED_MODEL_DISCOVERY_SYNTHESIS.md)). So "A Bayes-primary → PAYOFF
> LOW / lasso stays primary" and "F full family pooling → subsumed" describe **re-estimating the same observational
> data with a Bayesian prior** (verified −0.002 — no *new information*); they do **NOT** verdict against folding
> **exogenous** channels (CN dose, protein) or the causal/exclusion machinery **into** the one posterior, which is a
> different move (adds information / identifies flat directions). Don't quote these PAYOFF verdicts as "Bayes doesn't
> help" for the CN-anchor / unification work.

The design (`LEARNED_MODEL_DESIGN_RESPONSE.md`) = 10 decisions + 5 bars. **Done:** J (CN instrument) + Bars 1–4.
The rest, each with `DESIGN → BUILT → the GAP → APPLY (concrete build/test) → PAYOFF`:

**Deliberate skips (justified — leave unless a result demands it):**
- **C occupancy link.** DESIGN: saturating a/(a+K) + shared free-AGO pool (retire D(m)). BUILT: linear a·M.
  GAP: no saturation. Skip reason: the occupancy substrate **FAILED the held-out gate**.
- **G cooperativity/ceRNA.** Sits on C's failed substrate → skip until C earns its place.
- **I Shapley identity. ✅ BUILT (2026-07-06, `attribution.identity_vs_magnitude`)** — Shapley/LMG credit for
  the aggregate's R² across families = collinearity-fair coupling OWNERSHIP, vs budget = realized force.
  Delivers unique-owners vs magnitude-passengers (ESR1 miR-18 owns 0.54 on 0.20 budget; miR-22 budget 0.26 /
  identity 0.01). No longer a skip. (METHODS §15, ATTRIBUTION_IDENTITY_VS_MAGNITUDE §0.)

**Partials — the real menu:**
- **A Bayes-primary.** DESIGN: hierarchical Bayes as the *primary* coupling estimator. BUILT: lasso primary,
  Gibbs = uncertainty only. GAP: primary isn't Bayesian. APPLY: make the program Gibbs the primary coupling,
  compare OOF vs lasso. **PAYOFF: LOW** — pooling already shown not to improve coupling → recommend close.
- **B transcription-rate proxy** *(design's "most load-bearing")*. DESIGN: a proxy for the gene's TRANSCRIPTION
  RATE in C, so a positive M isn't a TF co-driving miRNA+target. BUILT: `mal_prolif` (proliferation program
  only); TF-regulon proxy OVER-controls (diagnostic). GAP: no orthogonal transcription-rate control beyond
  proliferation. APPLY: add an **intronic / pre-mRNA-read** proxy (nascent transcription, orthogonal to mature
  mRNA) as a C column; test it removes positive-weight artifacts WITHOUT over-controlling (the TF-proxy failure
  mode). **PAYOFF: HIGH but RISKY** (over-control) + needs intronic-read data (check availability first).
- **E direction channel.** DESIGN: fold TarBase-v9 Negative/Positive to EXCLUDE validated-ACTIVATING edges a
  count-based HE wrongly admits. BUILT: ledger+scanMiR+fused, no direction. GAP: activating "edges" not
  excluded. APPLY: ingest TarBase-v9 (Zenodo 14654626), map to HE, drop net-Positive edges, re-run gate +
  discovery. **PAYOFF: MEDIUM** (cleaner edge set) — data obtainable.
- **F full family pooling.** DESIGN: family μ + per-member δ (shrunk) + posterior width. BUILT: family collapse
  (∞-shrinkage) + conditioned nomination (§9, now **co-driver-aware**: single/co-drivers/family-only). GAP: no
  per-member δ with quantified uncertainty. APPLY: Phase-3b hierarchical family model; compare member
  resolution vs the nomination. **PAYOFF: LOW-MEDIUM** — the co-driver-aware nomination already does the job.
- **H Bayesian state nesting M_T=M_H+Δ. ✅ APPLIED (subtype flavour, 2026-07-05).** `subtype.subtype_wiring_pooled`
  = whole-cohort M_all as prior mean, per-subtype ridge-toward-M_all (regularised δ vs n-matched null). OUTCOME:
  it **corrected** the independent-fit over-call — PTEN wiring-remodeling is borderline (~1.35× across λ=1–30,
  not robust), only RB1-Her2 modest ⇒ **subtype wiring mostly conserved**; robust subtype signal = level + slope.
  REMAINING: the **state flavour** (M_H from GTEx as prior mean for M_T) once GTEx CIBERSORTx lands.
- **D spike-and-slab discovery.** DESIGN: soft π-inclusion alongside hard-HE. BUILT: lasso+permutation-FDR
  (`priors.py` stub unused). GAP: inclusion is frequentist FDR, not Bayesian π. APPLY: wire the spike-slab,
  compare to lasso-FDR. **PAYOFF: LOW** — lasso-FDR discovery works (FDR 0.001, validated).

**APPLIED OUTCOMES (2026-07-05 — applying the gaps *confirms* the design):**
- **E ✅** (`evidence/direction.py`) — TarBase-v9 covers 57% of HE; only **6 edges validated-activating (0.2%)**
  (miR-21→MMP9, miR-205→BCL2, miR-126→PITPNC1/MERTK, miR-663a→CDK1, miR-335→BRCA1) → HE set directionally
  clean; flag/drop the 6, negligible impact. (TarBase v9 is in-repo: `data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz`.)
- **A ✅** — pooling gate (Bayes-primary's test) Δρ(pooled−solo) = −0.002, hurts small-n ⇒ **lasso stays primary**.
- **D — subsumed** — discovery's permutation-FDR *is* the inclusion inference; spike-slab wouldn't change the
  FDR-0.001 edge set. **F — subsumed** — the co-driver-aware nomination already resolves members.
- **H ✅** (subtype flavour) — see above; **B** data-blocked (intronic reads).
**Remaining real work:** H **state** flavour (needs GTEx CIBERSORTx), and the **CPTAC-protein Bar-5** validation.
None of the partials beat their stand-ins; the design is robust to them.
