# Handoff — PROGRESSION / STATE axis session

> ⚠ **SUPERSEDED IN PART (2026-07-12) → read `../LEARNED_MODEL_STATE_CHANNEL_PLAN.md`.** Two premises below are
> now refuted by measurement: (1) *"the real first task is REFERENCE RECONCILIATION"* — **no**: `gtex_wu_major`
> exists (514, all 327 paired donors), and β_H is panel-invariant anyway (ρ=0.94 across atlases vs ρ=0.58 for
> no-C-vs-C). Zero CIBERSORTx runs are needed; the real first task is that **GTEx has no confounder block wired at
> all**. (2) *"`M_T = M_H + Δ`"* — the additive form is composition-attenuation-biased and would manufacture
> "acquired wiring"; use `M_T = a·M_H + Δ`. The inventory and reading list below remain valid.

> **You are one of TWO parallel sessions.** The other takes the **CPTAC / PROTEIN** channel
> (`HANDOFF_CPTAC_PROTEIN.md`). You share the fusion machinery, the inference engine, the gauge, and — most
> critically — the **CIBERSORTx / composition** dependency. Read the **Intersections** section and coordinate.
> **Prepare your OWN plan** (deliberately not written here). Orientation + reading list only.
>
> **Read two things before this handoff:** (1) **`HANDOFF_MODEL_INFRASTRUCTURE.md`** — the model object, its
> fundamental units (family-β, X_fam, the per-gene A1 fit unit), and the code spine + channel template you plug into;
> (2) **the WHOLE `PROGRAM_FORWARD_BOARD.md`** end-to-end — the consolidated open-work map (topics **A Fusion/Bayes ·
> B CN instrument · C Progression/state · D Discovery · E Validation · F Infra**) and the master index of everything
> standing and moving in parallel. Your axis is §C (+ §F for its CIBERSORTx blocker), but its intersections reach across
> the board — the protein session in §A③/§E, the discovery consumers in §D, the composition/infra in §F. Don't read only
> your slice.

---

## Your axis, in one line
Fold the **STATE / progression** axis into the model: **M_T = M_H + Δ** — healthy GTEx wiring `M_H` as the prior mean
for tumour wiring — the *acquired-vs-constitutive wiring* readout, plus the GTEx→NAT→tumour trajectory and
state-stratified discovery.

## Where it stands (verified)
- **PARTLY MATURE, partly blocked — READ the built infra before rebuilding.** The **trajectory / parallel-view arcs
  are BUILT + RUN**: `learned/parallel_view.py`, `learned/state.py`, ~25 tissue_reference lanes, the trio
  GTEx→NAT→tumour abundance-logFC vs wiring-ΔM. Memory `cn-and-trajectory-arcs-maturity` is explicit that these are
  mature and covered by canonical docs — **read them first, don't re-derive.**
- **The STATE CHANNEL (M_T = M_H + Δ) is NOT built — but the deconvolutions mostly ARE.** The home is
  **`data/external/cibersortx/`** (NOT `data/deconvolution/`). **Inventory (verified 2026-07-12):**
  - **References present:** **HBCA** (Human Breast Cell Atlas single-cell ref, ~11 breast-resolved types: luminal
    adaptive-secretory-precursor / luminal hormone-sensing / basal-myoepithelial / DDC1 / DDC2 / fibroblast / perivascular /
    vascular + lymphatic endothelial / myeloid / lymphoid) · **Wu-major** (9 lineages) · a single-cell **scref** panel.
  - **GTEx × HBCA — DONE = M_H:** `gtex_hbca_level1/CIBERSORTx_Adjusted.txt` (514 GTEx-breast samples, GTEX-* ids).
  - **NAT × HBCA — DONE = M_NAT:** `out_hbca_nat` (113) + `nat_hbca11` — the intermediate trajectory state, SAME HBCA panel.
  - **TCGA tumour × Wu-major — the model's CURRENT C:** `tcga_wu_major` / `out_NAT` (~1172) → surfaced as
    `output/brca_deconvolution/tcga_cibersortx_fractions.tsv`, of which `data._DECONV_COLS` uses **8** (Cancer Epithelial
    excluded = the malignant compartment, used separately for `mal_prolif`). Root `CIBERSORT-Results.txt` = a **redundant**
    earlier copy of this Wu panel (superseded — dedup, don't use).
  - **Other TCGA variants:** `tcga_tumor_no_nl` (scref, 1059, no-normal-like), `tcga_nat_alone`, `tcga_tumor_nat_pooled`.
  - **DEFERRED fine subtypes:** `data/deconvolution/cibersortx/hires_out/` = HiRes GEP **10 immune subtypes** + IFNG
    imputation — not wired; the level1(coarse)-vs-level2(fine) over-control ablation runs on this, no new run needed.
  - **⚠ THE ACTUAL FIRST TASK is REFERENCE RECONCILIATION, not a missing run:** GTEx + NAT are on **HBCA**, but the model's
    TCGA C is on **Wu-major**, and there is **no `tcga_hbca`**. For M_H↔M_T comparability the panels must match — either
    deconvolve TCGA on HBCA (the sigmatrix + refsample are already in the dir) or pick one panel across all three cohorts.
    Decide this before folding the state channel; it also decides which panel CPTAC's C should use.
- **Two homes for state (axis Q):** a **PRIOR** (`M_H` → Gaussian prior mean → Gibbs-cheap; the Decision-H flavour
  already validated) **OR** a **2nd likelihood** (healthy Y observed → healthy & tumour β co-estimated). The prior form
  is conjugate and stays on Gibbs.
- **Coverage:** GTEx-breast ~300 / NAT ~104 — power-bounded (state-stratified discovery is limited by NAT n=104).
- **Gauge (E) + the GTEx logFC gauge check** must hold before any cross-state β claim: validate the miRNA RPM≈TPM
  assumption across GTEx↔tumour (rawFC vs relFC) — the RPM↔TPM scale is *relative*, per memory `parallel-view-abundance-vs-wiring`.

## Intersections with the CPTAC / PROTEIN session — coordinate at these seams
1. **Composition / CIBERSORTx is the shared crux.** The **TCGA Wu-major fractions already feed the confounder block C for
   BOTH channels** (done + wired). Two shared decisions: (a) the **reference reconciliation** (GTEx/NAT are on HBCA, the
   model's TCGA C is on Wu-major) — whichever panel you standardize on becomes the C for BOTH channels, so decide it jointly
   with CPTAC, not solo; (b) the **fine (level2) immune panel** exists but is deferred — turning it on for the state
   over-control ablation also changes CPTAC's discordance model (protein discordance is often composition-driven), so
   coordinate. One panel, shared.
2. **Same fusion machinery + the engine (axis J).** Your state channel is Gaussian-conjugate (stays **J1/Gibbs**). The
   protein channel is the likely trigger for **J2/HMC** — and if it goes there, the whole posterior (incl. your state
   channel) comes along. Decide the engine jointly; if they build an HMC sampler, build your channel J2-compatible.
3. **State-conditional protein is where you two intersect.** Protein discordance may differ tumour vs NAT — a
   **(protein × state)** fusion is the natural joint object if both channels land. Flag it, don't build it solo.
4. **Shared gauge (E)** — both re-run the single-edge gauge for their own channel.

## What to read (in order)
0. **`HANDOFF_MODEL_INFRASTRUCTURE.md`** (model, units, code spine) **+ the WHOLE `PROGRAM_FORWARD_BOARD.md`** (all of
   §A–§F — the open-work landscape; your entries are §C + §F). Read both end-to-end before the spec docs.
1. **`LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md`** — **§Q** (state as a *dimension* of the latent vs a *channel*),
   **§I2** (state channel), **§7 channel catalog** (the **state / healthy** row is your spec — prior-mean nesting,
   Gaussian, coverage), **§E** (gauge), **§J** (engine). Plus §0/thesis + the status header.
2. **`LEARNED_MODEL_DESIGN_RESPONSE.md`** — **Decision H** (the state axis GTEx→NAT→tumour, §7) + **§2.7**
   (state / field-effect biology).
3. **`learned/parallel_view.py` + `learned/state.py`** (the BUILT trajectory infra — read, don't rebuild) +
   **`PRIMARY_TUMOR_ROADMAP.md`** (the phased state roadmap: composition-retest, spatial, orphan retrofit).
4. **`PROGRAM_FORWARD_BOARD.md`** — **§C** (progression/state — your board entries) + **§F** (the CIBERSORTx blocker) +
   `WHATS_NEXT.md` §PARALLEL-VIEW.
5. **Memory** — `cn-and-trajectory-arcs-maturity` (the arcs are mature+run — read the docs first),
   `parallel-view-abundance-vs-wiring` (design-model M vs abundance; the GTEx gauge caveats).
6. **`HANDOFF_CPTAC_PROTEIN.md`** — the parallel session, for the intersections above.

## Your job
Produce a plan: decide the **CIBERSORTx unblock path** (it gates you AND cleans C for CPTAC), **prior-form vs
2nd-likelihood** for `M_T = M_H + Δ`, how much of the trajectory infra to reuse vs fold into the posterior, and the
**composition + engine** coordination with the CPTAC session before committing shared machinery.
