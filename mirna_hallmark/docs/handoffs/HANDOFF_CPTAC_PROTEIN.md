# Handoff — CPTAC / PROTEIN channel session

> **You are one of TWO parallel sessions.** The other takes the **PROGRESSION / STATE** axis
> (`HANDOFF_PROGRESSION_STATE.md`). You share the fusion machinery, the inference engine, the gauge, and a
> composition/deconvolution dependency — so read the **Intersections** section and coordinate at those seams.
> **Prepare your OWN plan** (deliberately not written here). This handoff is orientation + a reading list only.
>
> **Read two things before this handoff:** (1) **`HANDOFF_MODEL_INFRASTRUCTURE.md`** — the model object, its
> fundamental units (family-β, X_fam, the per-gene A1 fit unit), and the code spine + channel template you plug into;
> (2) **the WHOLE `PROGRAM_FORWARD_BOARD.md`** end-to-end — it is the consolidated open-work map (topics **A Fusion/Bayes ·
> B CN instrument · C Progression/state · D Discovery · E Validation · F Infra**) and the master index of everything
> standing and everything moving in parallel. Your axis is a couple of entries on it, but its intersections reach across
> the board — the state session lives in §C/§F, your validation twin in §E, the discovery consumers in §D. Don't read only
> your slice.

---

## Your axis, in one line
Fold **PROTEIN (CPTAC)** into the model as a **second likelihood** — the *realized-at-protein* repression the
tumour-mRNA channel cannot see. Per the current program state this is **THE sole remaining fundamental lever**:
every internal refinement (CN channel, δ-pooling, coding/host gate, isomiR refit, Student-t likelihood, cross-gene
pooling) has landed *immaterial-at-n≈1000*, and the fit-unit A-axis is now closed to A1 (per-gene). Only an
**exogenous, non-redundant** channel can move past the mRNA floor — protein is it.

## Where it stands (verified)
- **NOT built as a channel.** It exists only as **validation** so far: `eval/cptac_validation.py` did
  miR-135b→GATA3 (registry **MH-33/34**). The *systematic* Bar-5 protein validation is flagged the
  "highest-value validation left" — but only one edge is done.
- **Two roles — hold both in mind:**
  - **(a) Validation (Bar-5, cheap first step):** does the learned **M** predict *where protein diverges from mRNA*
    (discordance)? This is the low-risk, Gibbs-compatible entry.
  - **(b) Channel (§I1, the frontier):** fold protein IN as a 2nd likelihood on the same latent M.
- **The link is the fork (axis L).** A plain-Gaussian protein link stays on Gibbs (cheap). The *honest* discordance
  model — protein carries a translational term the mRNA doesn't, possibly **saturating** — is **non-conjugate → it is
  the first thing that forces the engine from Gibbs (J1) to HMC (J2, numpyro/Stan)**. The design's own sequence: fit
  **protein-as-linear first**; if its residual structure vs mRNA is strong, graduate to J2.
- **Coverage is the constraint:** CPTAC-BRCA is n-limited (~120). This is a *small-n* channel — its value is
  **non-redundancy**, not power.
- **Gauge (axis E) must hold first.** The protein-channel β must land on the SAME scale as the mRNA β — the exact
  make-or-break single-edge equality that gated the CN channel. Run that check before any multi-edge claim.

## Intersections with the PROGRESSION / STATE session — coordinate at these seams
1. **Same fusion machinery.** Both fold a channel into ONE posterior via `spike_slab._gibbs_posterior(channels=)` +
   the `channel_cn.multi_family_terms` template. Don't both refactor the channel API differently — align on it.
2. **The inference engine (axis J) is shared and YOU are the one who moves it.** The state channel is
   Gaussian-conjugate (stays J1/Gibbs); your protein channel is the likely trigger for **J2/HMC**. If you go to HMC,
   the *whole* posterior (including their state channel) comes along — so decide the engine **jointly**. If you build a
   J2 sampler, they reuse it.
3. **Composition / deconvolution is shared.** The state session's **CIBERSORTx** run produces the composition fractions
   that feed the confounder block **C** (axis N) used by ALL channels — including yours. Protein discordance is often
   **composition-driven**, so you need the same composition adjustment. Their blocker is partly your blocker.
4. **State-conditional protein is where you two literally intersect.** Protein-vs-mRNA discordance may **differ by
   state** (tumour vs NAT). If both channels land, a **(protein × state)** fusion is the natural joint object — flag it,
   don't build it solo.
5. **Shared gauge (E).** Both sessions re-run the single-edge gauge for their own channel.

## What to read (in order)
0. **`HANDOFF_MODEL_INFRASTRUCTURE.md`** (model, units, code spine) **+ the WHOLE `PROGRAM_FORWARD_BOARD.md`** (all of
   §A–§F — the open-work landscape; "Top of the board" has you at #1). Read both end-to-end before the spec docs.
1. **`LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md`** — the frame + your spec:
   - §0 / thesis + the **status header** (current built state, the "immaterial-at-n≈1000" through-line).
   - **§7 channel catalog** — the **protein** row IS your spec (latent, link, conjugacy, coverage, null).
   - **§I1** (protein channel), **§L** (link form: L1 linear → L2 saturating / L3 threshold), **§J** (engine J1→J2),
     **§E** (gauge — *must-run first*), and the §7 per-channel notes (the linear-first-then-J2 sequence).
2. **`LEARNED_MODEL_DESIGN_RESPONSE.md`** — **§2.6** (protein-level & discordance, literature) + **Decision B**
   (observation target Y — protein vs mRNA; the doc calls it the most load-bearing decision).
3. **`eval/cptac_validation.py`** (existing infra) + registry **MH-33/34** (the miR-135b→GATA3 precedent).
4. **`PROGRAM_FORWARD_BOARD.md`** — §A③ + §E (your board entries) + "Top of the board" (you are #1).
5. **Memory** — `external-cohorts-and-apa-sources` (CPTAC availability/scope), `cptac-batch-is-enrolling-site`
   (the batch proxy = enrolling site).
6. **`HANDOFF_PROGRESSION_STATE.md`** — the parallel session, for the intersections above.

## Your job
Produce a plan: choose **validation-first (Bar-5, cheap, Gibbs)** vs **channel-fold-in (the frontier, may force HMC)**;
pick the **link form** (linear first; graduate only if discordance structure is real); and settle the **engine +
composition** coordination with the state session before either of you commits the shared machinery.
