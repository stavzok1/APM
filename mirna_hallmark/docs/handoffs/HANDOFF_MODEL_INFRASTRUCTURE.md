# Handoff — the model, its fundamental units, and the shared infrastructure

> **Read this first — both parallel sessions (CPTAC/protein and progression/state) build on exactly this.** It is the
> foundation the axis handoffs assume. Nothing here is your axis's plan; it is what the model *is* and what code you
> will plug into. Verified against the current tree (2026-07-12).

---

## 1. The model object — one posterior, four readouts
The converged model (`LEARNED_MODEL_DISCOVERY_SYNTHESIS §6b`) is **ONE Bayesian posterior over a latent repression
field M**, read four ways: **coupling** (π≡1 dense) · **attribution** (Shapley) · **identifiability** (parallel-view) ·
**discovery** (evidence-π). It currently fuses **two channels**: a curated-evidence **prior** (π/μ/τ) and the
tumour-**mRNA likelihood**. Everything else (CN, protein, state, binding) is either external validation or a channel
waiting to be folded in. The broader frame (`CHANNEL_FUSION_DESIGN §0`): *M is a latent causal repression parameter;
every data source is a noisy **channel** observing it through its own link, fused in ONE posterior.* **Your axis adds a
channel to this object** — you do not build a separate model.

## 2. The fundamental units (get these exactly right)
- **The latent β is PER SEED FAMILY** (`β_f`). Same-seed members repress identically (§8) ⇒ **there is no per-arm
  repression β.** "Which arm delivered the dose" is a *dose-delivery* readout (`within_family.delta_pooling`), a
  separate/orthogonal latent — NOT a coordinate of the repression posterior. Do not add an arm-level repression β.
- **The exposure/predictor is `X_fam`** — the family-pooled arm abundance `log2(1+Σ_m RPM_m)`, built by
  `families.collapse_by_family`. This is *the* choke point: ~22 modules inherit it. (It now has an opt-in `isomir=`
  seed-correction, default OFF — a coupling wash, don't turn it on.)
- **The fit unit is A1 = per-gene** (RESOLVED 2026-07-12, memory `a-axis-fit-unit-resolved`): a full arm×gene field
  does NOT beat per-gene at n≈1000. So you fit one posterior per gene, with cross-gene structure entering only through
  the gene-independent first-stage γ and the μ-priors — **not** a joint field.
- **The estimand is family→gene repression coupling**, in a repression gauge (positive = repression, on the
  z-scored/residualized `r = −resid(Y|C)` scale).
- **The recurring empirical floor:** every internal lever washes at n≈1000 (CN channel ~0.7% weight, δ-pooling,
  coding/host gate, isomiR, Student-t, cross-gene pooling). Your channel matters only if it is **exogenous and
  non-redundant** — design for that, not for squeezing the mRNA.

## 3. The code spine — the pipeline you plug into
```
data.assemble_gene(gene, he_only=False)         →  (Y, X, C, w_prior)      # per-gene data assembly (the entry point)
     Y : Series[participant]        target gene log2(TPM+1)
     X : DataFrame[participant×arm] regulator arm log2(RPM+1)
     C : DataFrame[participant×conf] confounder block  (CPE purity, mal_prolif, target_cn, deconv fractions…)
     w_prior : Series[arm]          adaptive-penalty prior weight (evidence)
families.family_of(arms)                         →  arm → seed family
families.collapse_by_family(X, w_prior, fam)     →  (X_fam[part×family], w_fam, members)   # THE predictor choke point
attribution_eb._prep(Y, X_fam, C)                →  (yr, Xz, cols)   # residualize on C + z-score → the sampler's inputs
spike_slab._gibbs_posterior(Xz, yr, pi, *, channels=, nu=, return_samples=)   →  (mean, sd, PIP) of β   # THE posterior
```
- **`_gibbs_posterior`** is the ONE sampler. Args that matter for you: **`channels=`** (the fusion API — a list of
  channel terms, see §4), **`nu=`** (Student-t robustness, built MH-92), **`return_samples=`** (posterior draws for
  Shapley/width). `pi≡1` gives the dense coupling readout; an evidence-`pi` gives the discovery readout.
- **The confounder block `C`** (axis N) is shared by every channel — it is where composition/purity/proliferation
  adjustment lives. Whatever you add to C affects *all* channels.

## 4. The channel template — how a channel is added (copy this)
The CN channel is the worked example: **`channel_cn.multi_family_terms(gene, cols, members)`** returns the `channels=`
list. **Each channel term is a Gaussian observation of a linear combination of β:**
```python
{ "members": [(col_index, gamma), …],   # which β's this observation constrains, and their loadings γ
  "pihat":   float,                       # the observed reduced-form value  (e.g. Z→Y)
  "s2":      float }                      # its variance  (sampling SE² + a bias/pleiotropy floor)
```
`_gibbs_posterior` adds the term `N(pihat; Σ γ·β, s2)` to the Gibbs sweep — so a channel that constrains a linear
combination of β's is a **cheap conjugate Gaussian term** and STAYS on Gibbs (J1). Your channel is: *what linear
function of β does my data observe, with what noise?*
- **Gauge (axis E) is make-or-break, run it FIRST.** `pihat` must be on the **same β scale** as the mRNA channel
  (the raw-`r` convention: z-score Z and X but NOT `r`; see `channel_cn` §E gauge fix). The single-edge check: on a
  clean strong edge, does the *channel-only* posterior-mean β equal the *mRNA-only* β? If not, the gauge is wrong and
  nothing downstream is valid.
- **Validation ladder (every channel passes it):** (i) **nests §6b** — channel off ⇒ bit-identical to the current
  posterior; (ii) **gauge** equality; (iii) **shuffled null collapses** the channel — the null primitive is built:
  `instrument.exclusion(shuffle_cn=True)` / `channel_cn.multi_family_terms(shuffle_cn=True)` (permute the source →
  dead channel; the channel's weight must collapse). Build the analogous shuffle for your source.

## 5. Shared cross-cutting infrastructure (both axes touch these — coordinate)
- **Confounder block C / composition.** The state session's **CIBERSORTx** deconvolution produces composition
  fractions that feed C — used by the protein channel too. One deconvolution, shared.
- **The gauge (E).** Both channels must land on the mRNA β scale; both re-run the single-edge check.
- **The inference engine (J).** Gaussian-conjugate channels stay on **Gibbs (J1)**. A non-conjugate link (the protein
  discordance/saturating model) forces **HMC (J2, numpyro/Stan)** — and moving to J2 moves the *whole* posterior,
  including the other axis's channel. Decide the engine jointly.
- **The evidence prior** (`priors.py` / `learned.evidence.ledger`) supplies π/μ/τ; a channel enters θ (magnitude) or π
  (inclusion/discovery) — know which readout yours feeds.
- **The null primitive** (shuffle the source) is reusable — the CN version just landed as `exclusion(shuffle_cn=)`.

## 6. Where each axis plugs in
- **CPTAC / protein** → a new likelihood on `β·(mRNA-repression)` observed at protein; link is the fork (linear→Gibbs,
  discordance/saturating→HMC); coverage ~120. Enters `_gibbs_posterior(channels=)` as a term (linear) or forces J2.
- **Progression / state** → `M_T = M_H + Δ`: healthy GTEx `M_H` as a Gaussian **prior mean** (Gibbs-cheap) or a 2nd
  mRNA-likelihood on healthy Y; the trajectory/parallel-view infra (`parallel_view.py`, `state.py`) is already BUILT —
  read it, don't rebuild; blocked on CIBERSORTx.

**Read next:** the **WHOLE `PROGRAM_FORWARD_BOARD.md`** (the consolidated open-work map, topics A–F — the master index of
everything standing and moving in parallel), then your axis handoff (`HANDOFF_CPTAC_PROTEIN.md` or
`HANDOFF_PROGRESSION_STATE.md`), then `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` §7 (channel catalog — your row) + the
sections it names.
