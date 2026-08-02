"""THE FOUR READOUTS of the ONE posterior, genome-wide — the canonical per-edge / per-gene table.

`LEARNED_MODEL_DISCOVERY_SYNTHESIS §6b`: the model is **ONE Bayesian posterior over the latent repression field
M, per gene**, read four ways. Those readouts were scattered across method-comparison diagnostics
(`attribution_eb.{full,identifiability_full,selection_full}` compare estimators; they are NOT the readouts) and
per-gene helpers. This module emits them ONCE, together, on ONE confounder block.

Per gene, THREE Gibbs runs (`spike_slab._gibbs_posterior`, the one sampler) — 2 on the core-C design, 1 on
the deconv design:
  * **π≡1 DENSE**  → coupling · attribution (magnitude + identity) · identifiability   [core AND deconv]
  * **evidence-π**  → discovery                                                        [core ONLY]
⚡ Discovery is defined on the CANONICAL core-C design, so the evidence-π chain runs for core only
(`discovery=False` on the deconv call). It used to run for BOTH and `dec["pip_disc"]` was never read — 1 of
4 chains per gene, ~25% of the module's runtime, computed and discarded. Gating it is OUTPUT-IDENTICAL
(verified: max|Δ| = 0.0 across all 23 numeric columns, 5,117 rows) and took the run 3.0 → 2.3 min.

Read four ways (all per SEED FAMILY — the latent is family-level, §8; arms are mapped back at the end):

  COUPLING            `beta` (posterior mean), `beta_sd`
  ATTRIB / MAGNITUDE  `beta` — the force  (the WHO-HOW-MUCH)
  ATTRIB / IDENTITY   `identity` — Shapley/LMG credit for **R²(Xz·m, yr)**, m = **bagged NNLS** (MH-138).
                      The DOCTRINE's identity: the value function v(S) = R²(X_S·m_S, y) is **NON-additive** —
                      a collinear pair's joint R² is far LESS than the sum of its marginals — so the
                      ordering-average genuinely **SPLITS the overlap**. Normalised to a share of explained
                      variance. ⚠ **SCOPE:** identity says WHO, magnitude says HOW MUCH, **neither says
                      WHETHER**. Per MH-126 the model **RANKS** the canonical family above chance
                      (0.370, p=0.018) but does **NOT NAME** it (argmax at chance at every k) — read the
                      ranking, not the argmax.
  ATTRIB / MAGNITUDE  `beta_frac` ± `beta_frac_sd` — β_f/Σβ, the gene-budget FRACTION.
                      ⚠ **RENAMED from `share` 2026-07-17 (MH-138), and DEMOTED from "identity".** It was
                      labelled "the Bayesian Shapley" but for the **ADDITIVE** aggregate ρ = Σ_f β_f·X_f the
                      Shapley value of f is **TRIVIALLY β_f** ⇒ normalising **splits NOTHING** under
                      collinearity (MH-122) — the one thing a Shapley exists to do. It is a β-fraction; the
                      name now says so. ⚠ Its denominator is majority-unidentified mass (**73.0% pooled**),
                      because the half-normal slab has a strictly positive mean and so **cannot zero** an
                      un-informed family (`beta == 0` in **0/5,117**). Read it via `beta_frac_reliable`.
  IDENTIFIABILITY     `z = beta/beta_sd` (|z|>2 = identified; MH-83: precision 0.86 / recall 0.89 vs the
                      conditioned-partial ground truth) and `beta_frac_sd` (WIDE = genuine non-identifiability
                      — MH-94's PTEN miR-141/200a 0.77±0.41 that a point Shapley hid)
  DISCOVERY           `pip_discovery` — inclusion posterior under the evidence-π prior

⚠⚠ **THE `deconv` BLOCK IS NOT A FREE "MORE CONTROL IS BETTER" KNOB — COMPOSITION MAY BE THE SIGNAL ITSELF.**
`LEARNED_MODEL_METHODS.md §1`: *"tumour attribution = fractions; **gate/coupling = core only unless deconv
requested**"* and *"**Cancer-Epithelial is deliberately excluded from the `deconv` block** — conditioning on the
compartment the target is expressed in **over-controls the signal**."* `DESIGN_RESPONSE` Decision B: the
deconvolution fraction is *"a composition covariate and a **robustness variant**"*. **Rationale:** a miRNA acting
through a cell-STATE shift (e.g. EMT) *produces* a composition change — so the composition is partly the MECHANISM,
not a nuisance. Conditioning it away destroys real signal (the same confound-vs-mechanism trap as MH-100/MH-101).

**So this module runs BOTH blocks and reports the RETENTION, rather than picking one:**
  * **core C**   = CPE + target_cn + mal_prolif                      → the CANONICAL coupling/identifiability
  * **+deconv**  = core + the 8 NON-malignant Wu-major lineages
                  (Cancer-Epithelial HELD OUT — see above)           → the CELL-INTRINSIC attribution
  * `retention = beta_deconv / beta_core` per family, classed as in `card.py`:
        **cell_intrinsic (≥0.7) · partial · composition_explained (<0.4)**
An edge whose coupling is `composition_explained` is FLAGGED, not deleted — the reader decides.
(This is the check `cptac_validation` never did, which is how MH-107's EMT/ZEB1 result was presented as
cell-intrinsic when it was compartment arithmetic.)

CLI:  .venv/bin/python3 -m mirna_hallmark.learned.readouts [--workers 8] [--limit N]
Out:  output/learned/readouts_edges.tsv   (one row per gene × family, with member arms)
      output/learned/readouts_genes.tsv   (one row per gene — the aggregate roll-up)
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT_EDGES = ROOT / "mirna_hallmark/output/learned/readouts_edges.tsv"
OUT_GENES = ROOT / "mirna_hallmark/output/learned/readouts_genes.tsv"

# ⚡ 2000/700 → 1000/350 (MH-141, 2026-07-17). LICENSED BY THE decoy_bench DESIGN: the arbiter is the
# sampler's OWN seed-to-seed jitter at the reference config — if a shorter chain moves the readouts LESS
# than reseeding does, the cut is free. MEASURED over 38 genes (p=2..64, incl. every hub):
#   * beta and beta_sd deviations both INSIDE the reseed floor; 500/200 and 200/80 are NOT (1.4–2.1× it).
#   * the deviation is NOISE, not BIAS — signed mean −9.7e-06 (beta, p=0.85) / +4.2e-05 (beta_sd, p=0.47),
#     i.e. 0.007σ / 0.029σ of the floor ⇒ burn-in at 350 fully forgets the init. (This was the real risk: a
#     short burn-in is SYSTEMATIC, and a max|Δ| comparison cannot tell a bias from noise.)
#   * genome-wide: `identified` flips 42/5117 = 0.82%, BELOW the reseed floor (1.6–1.9%); beta corr 0.9994,
#     beta_sd corr 0.9967; `identity` is BIT-IDENTICAL (corr 1.000000) — it is computed from the NNLS
#     weights, not the Gibbs mean (MH-138), so it does not depend on chain length at all.
# ⚠ My prediction was BACKWARDS: I expected beta to converge fast and beta_sd to be the constraint. At
# 200/80 beta is 2.1× the floor while beta_sd is only 1.3× — burn-in governs where the chain SITS more than
# how wide it is. Do not shorten further without re-running the floor test (`scratchpad/niter_sweep.py`).
N_ITER, BURN = 1000, 350
# Sampled-permutation count for the IDENTITY readout (`attribution.shapley_identity`). 400 is the value
# `identity_vs_magnitude` has always used. Cost measured 2026-07-17: 0.16 s at 4 families → 1.43 s at 64
# (mean 0.85 s/gene incl. the NNLS fit) ⇒ ~2.7 min genome-wide on 8 workers. Exact Shapley is 2^p and
# infeasible at p=64; the sampled estimator is the reason this column can exist at all.
IDENTITY_PERM = 400

_AUX: dict = {}


def add_reliability(d: pd.DataFrame) -> pd.DataFrame:
    """MH-119 — GATE THE RATIO READOUTS. `beta_frac` and `retention` are both X/Y with a Y that can vanish; a
    vanishing denominator makes the QUESTION ill-posed, not just the estimate noisy.

    **`beta_frac` (β_f/Σβ; called `share` and mislabelled "the Bayesian Shapley" until MH-138).** MEASURED
    mechanism: it explodes ONLY when a gene's βs CANCEL — repressive (β>0) and anti-repressive (β<0) units
    summing to ≈0 while Σ|β| stays large. Of the rows with
    |beta_frac|>2, **100%** sit in a gene with ≥1 negative β (only 8% of genes have any), their |Σβ|/Σ|β| is
    **0.172** vs **1.000** for normal rows, and their |Σβ| is 10x smaller. CDH11's βs cancel to −0.0018, so a
    β of +0.072 reads as "999% of the total" — arithmetically right, semantically meaningless. **The model
    already knew: `beta_frac_sd` (then `share_sd`) is 316 on those rows vs ~0.04 normally.** It was never read.

    **`retention` (β_deconv/β_core).** Same family, different mechanism: a SMALL |β_core| (not cancellation).
    Gate on the core β actually being identified.

    ⚠ **MH-124 CORRECTED THE CAUSE, and MH-138 the LABEL.** The negative βs behind the explosion were a
    SAMPLER BUG (`_rtnorm_pos` inverse-CDF saturating for mu/s < −7.0345 ⇒ deterministic NEGATIVE draws,
    impossible under a half-normal slab), not "anti-repressive biology" as first recorded. The gates below
    are correct regardless — but the `identity` column, not `beta_frac`, is now the identity readout.

    Emits (FLAG, don't delete — the raw values are kept):
      `net_pressure`     |Σβ|/Σ|β| ∈ [0,1] — the gene's SIGN COHERENCE. 1 = all units push the same way.
      `beta_frac_abs`      |β_f|/Σ|β_f| ∈ [0,1] — a BOUNDED companion to `beta_frac`, always well-defined.
      `beta_frac_reliable` the signed fraction is meaningful (net_pressure ≥ 0.5 AND beta_frac_sd ≤ 1).
      `retention_reliable`  the core β is identified (|z|>2), so the ratio has a real denominator.
    """
    g = d.groupby("gene")["beta"]
    sb = g.transform("sum")
    sa = g.transform(lambda s: s.abs().sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        d["net_pressure"] = (sb.abs() / sa.replace(0, np.nan)).clip(0, 1)
        d["beta_frac_abs"] = (d["beta"].abs() / sa.replace(0, np.nan)).clip(0, 1)
    d["beta_frac_reliable"] = (d["net_pressure"] >= 0.5) & (d["beta_frac_sd"].abs() <= 1.0)
    if "retention" in d.columns:
        d["retention_reliable"] = d.get("identified", pd.Series(False, index=d.index)).astype(bool)
    return d


def _detection() -> pd.Series:
    """Fraction of patients in which each arm is measured. MH-117: this is the ONE measurement proxy the
    arm-level β tracks within a family (+0.305, p=4e-3) — abundance (p=0.65) and variance (p=0.95) do not.
    Emitted so a reader can GATE on it; NOT used to drop arms (axiom 2a: flag, don't delete)."""
    if "det" not in _AUX:
        from mirna_hallmark.learned import data as LD
        X = LD._load()["X"]
        _AUX["det"] = pd.Series(np.isfinite(X.to_numpy(float)).mean(1), index=X.index)
    return _AUX["det"]


def _ago(arms) -> pd.Series:
    """AGO-loading (Manakov chimeric eCLIP) — the INDEPENDENT biological guide label."""
    if "ago" not in _AUX:
        from mirna_hallmark.learned import ago_loading as AG, data as LD
        _AUX["ago"] = AG.loading(list(LD._load()["X"].index))
    return _AUX["ago"]


def _posteriors(gene, deconv, n_iter, burn, seed, level="family", *, discovery: bool = True):
    """One design (given a C block) → the dense posterior (+draws) and the evidence-π discovery posterior.

    `level="arm"` runs the SAME estimator with the §8 SEED-FAMILY COLLAPSE REMOVED — the design is the raw
    per-ARM matrix. Nothing else changes, so any difference is attributable to the collapse alone.

    ⚠⚠ **SCOPE — READ BEFORE USING THIS FOR ATTRIBUTION (MH-122).** The arm level is a **COUPLING/PREDICTION**
    variant, NOT an attribution one:
      * ✅ **What it is for.** The arm design predicts held-out mRNA BETTER out-of-sample (gene-level, CPTAC,
        β TCGA-fitted: Δρ −0.011, p=2.8e-11) and its member choice is fold-STABLE (87.3% vs 43.2% chance).
        Σβ is conserved (family↔arm ρ=0.994); what improves is the per-sample FIELD.
      * ⛔ **What it is NOT for: naming WHICH member represses.** `Design §F` — *"same-seed arms are
        near-collinear ⇒ the identified estimand is family→gene (arm = NOMINATION)"* — **STANDS.** Same seed ⇒
        identical sites ⇒ identical per-molecule repression ⇒ **only DOSE DELIVERY can differ**, and dose is a
        LEVEL quantity, not a variation one (`CN_INSTRUMENT §7`). MEASURED against the chimeric occupancy gold
        (MH-122, n=207, permutation null): CN dose 63.3% / abundance 63.8% (p<1e-4) vs **this β 44.9% and the
        collapsed `resolution_report` driver 41.5% — BOTH AT CHANCE** (and they agree with each other 85%: the
        same estimator, the same wrong question). It misses the canonical deliverers (CCND1←miR-16,
        ZEB1←miR-200c) that the dose lens recovers.
      * ⇒ For WITHIN-FAMILY attribution use the doctrine's ladder: `chimeric_evidence` (measured occupancy),
        else abundance / `instrument.cn_copy`. NOT this β. (`ATTRIBUTION_PRIMITIVE §3/§7`.)
      * ⚠ An earlier claim that this level "recovers within-family resolution the collapse discards"
        (MH-118/120/121) is **RETRACTED** — the Hansen-J "causal heterogeneity" was PLEIOTROPY (|δ_s| 0.153 vs
        0.088, p=0.028), and the CN corroboration was DOSE-DELIVERY, which is what a within-family share IS.
    """
    from mirna_hallmark.learned import data as LD, families as FAM, attribution_eb as AE, spike_slab as SS
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=deconv)
    fam = FAM.family_of(pd.Index(X.columns))
    if level == "arm":
        Xf, wf = X, w
        members = {a: [a] for a in X.columns}          # each arm is its own unit
    else:
        Xf, wf, members = FAM.collapse_by_family(X, w, fam)
    if Xf.shape[1] < 1:
        return None
    yr, Xz, cols = AE._prep(Y, Xf, C)
    n, p = Xz.shape
    b, sd, pip_d, draws = SS._gibbs_posterior(Xz, yr, np.ones(p), n_iter=n_iter, burn=burn, seed=seed,
                                              return_samples=True)                  # π≡1 DENSE
    # ⚡ The evidence-π chain is the DISCOVERY readout, and discovery is defined on the CANONICAL (core-C)
    # design only — `gene_readouts` reads `core["pip_disc"]` and never `dec["pip_disc"]`. It used to run
    # unconditionally for BOTH blocks, so every gene paid for a full n_iter chain whose result was discarded:
    # 1 of the 4 Gibbs runs per gene = ~25% of the module's entire runtime, thrown away. (cProfile 2026-07-17:
    # `_gibbs_posterior` is 87% of `gene_readouts`.) Gated on `discovery=`; the deconv block passes False.
    if discovery:
        pi = np.clip(AE._evidence_pi(gene, cols), 0.0, 1.0)
        _, _, pip_disc = SS._gibbs_posterior(Xz, yr, pi, n_iter=n_iter, burn=burn, seed=seed)   # evidence-π
    else:
        pi = pip_disc = None
    # ⭐ FIXED WEIGHTS for the IDENTITY readout (MH-138). Identity = Shapley/LMG on R² of a FIXED-weight
    # aggregate, and the weights MUST be able to reach zero. The Gibbs mean CANNOT: the half-normal slab
    # N⁺(0,ν²) has a strictly positive mean, so an un-informed family relaxes TO THE PRIOR rather than to 0
    # (MEASURED: `beta == 0` in 0/5,117 edges, 100% positive; pooled 73.0% of β-mass sits on |z|≤2 edges).
    # Feed that to a Shapley and you credit families the model cannot distinguish from noise. Bagged NNLS
    # returns exact zeros — this is RATIONALE §2e's "NNLS-style point estimate + EB posterior sd", and the
    # form `attribution.bayes_shapley_identity` already uses. Fit on the SAME prepped design as the Gibbs
    # (`AE._prep` "matches `_bagged_nnls`"; for TCGA tumour x_floor("01") == the 0.1 floor _prep applies), so
    # m and Xz are aligned by construction — no cross-module design mismatch.
    m_nnls, _ = AE._bagged_nnls_meansd(Xz, yr, seed=seed)
    return dict(cols=cols, n=n, p=p, b=b, sd=sd, pip_d=pip_d, draws=draws, pip_disc=pip_disc,
                pi=pi, members=members, fam=fam, Xz=Xz, yr=yr, m_nnls=m_nnls)


def gene_readouts(gene: str, *, n_iter: int = N_ITER, burn: int = BURN, seed: int = 0,
                  level: str = "family") -> Optional[pd.DataFrame]:
    """The four readouts for one gene → one row per seed family (or per ARM if `level="arm"`), under BOTH C
    blocks + the retention tag."""
    core = _posteriors(gene, False, n_iter, burn, seed, level)   # CANONICAL coupling/identifiability (design)
    if core is None:
        return None
    # `discovery=False`: the DISCOVERY readout is defined on the canonical core-C design only, and
    # `dec["pip_disc"]` was never read — running it here cost a full Gibbs chain per gene for nothing.
    dec = _posteriors(gene, True, n_iter, burn, seed, level, discovery=False)  # CELL-INTRINSIC attribution

    cols, b, sd, draws = core["cols"], core["b"], core["sd"], core["draws"]

    # ---- ATTRIBUTION / MAGNITUDE, as a FRACTION of the gene's budget. ----
    # ⚠ RENAMED `share` → `beta_frac` (MH-138). This was labelled "ATTRIB / IDENTITY — the Bayesian Shapley"
    # and it is NOT identity: for the ADDITIVE value function Σβ_f the Shapley value of f is TRIVIALLY β_f, so
    # normalising splits NOTHING under collinearity — the exact thing a Shapley exists to do (MH-122). It is a
    # β-fraction; the name now says so. It is kept because the FRACTION is a real magnitude readout, but read
    # it with `beta_frac_reliable` — its denominator is majority-unidentified mass (73.0% pooled).
    def _frac(dr):
        tot = dr.sum(1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            sh = np.where(np.abs(tot) > 1e-9, dr / tot, np.nan)                     # per-draw β_f/Σβ
        return np.nanmean(sh, 0), np.nanstd(sh, 0)
    beta_frac, beta_frac_sd = _frac(draws)

    # ---- ⭐ ATTRIBUTION / IDENTITY — the DOCTRINE's identity, for the first time in this table (MH-138). ----
    # Shapley/LMG credit for R²(Xz·m, yr) under the NON-additive value function v(S) = R²(X_S·m_S, y) — a
    # collinear pair's joint R² is far LESS than the sum of its marginals, so the ordering-average genuinely
    # DIVIDES the overlap. Weights `m` = bagged NNLS (MH-138: the Gibbs mean cannot zero an un-informed
    # family, so it would credit noise). Normalised to a share of the explained variance.
    # ⚠ SCOPE: identity says WHO, magnitude says HOW MUCH, NEITHER says WHETHER. A high identity share is the
    # DISTRIBUTION of coupling that exists — never evidence a driver exists. And per MH-126, the model RANKS
    # the canonical family above chance but does NOT NAME it (argmax at chance at every k) — so read the
    # ranking, not the argmax.
    from mirna_hallmark.learned import attribution as AT
    phi = AT.shapley_identity(core["Xz"], core["m_nnls"], core["yr"], n_perm=IDENTITY_PERM, seed=seed)
    tot_phi = phi.sum()
    identity = phi / tot_phi if tot_phi > 1e-12 else np.full_like(phi, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 1e-12, b / sd, 0.0)

    d = pd.DataFrame({
        "gene": gene, "family": cols, "n": core["n"], "p_fam": core["p"],
        "beta": b, "beta_sd": sd,                        # COUPLING (core C — the design's canonical)
        "z": z, "identified": np.abs(z) > 2.0,           # IDENTIFIABILITY (MH-83: prec .86 / rec .89)
        "identity": identity,                            # ⭐ ATTRIBUTION/IDENTITY — Shapley/LMG on R² (MH-138)
        "m_nnls": core["m_nnls"],                        #    the fixed weights identity is computed from
        "beta_frac": beta_frac, "beta_frac_sd": beta_frac_sd,   # ATTRIBUTION/MAGNITUDE as a fraction (NOT identity)
        "pip_dense": core["pip_d"],
        "pip_discovery": core["pip_disc"],               # DISCOVERY (evidence-π)
        "prior_pi": core["pi"],
        "arms": [";".join(sorted(core["members"].get(c, []))) for c in cols],
        "n_arms": [len(core["members"].get(c, [])) for c in cols],
    })

    if level == "arm":
        # Arm-level rows carry what makes them INTERPRETABLE and AUDITABLE:
        #   seed_family  — which members are competing for the credit (the §8 unit we just un-collapsed)
        #   detection    — MH-117: the ONE artifact proxy the arm-level beta tracks (+0.305, p=4e-3). Abundance
        #                  and variance do NOT predict it. FLAG, don't delete (axiom 2a) — the reader gates.
        #   ago_loading  — the independent BIOLOGICAL guide label (Manakov chimeric eCLIP). |beta| tracks it
        #                  within-family at +0.209 (p=1.2e-03) ⇒ the member picked is the loaded guide.
        d = d.rename(columns={"family": "arm"})
        fam = core["fam"]
        d["seed_family"] = d.arm.map(fam)
        d["family_size"] = d.seed_family.map(d.groupby("seed_family").size())
        d["detection"] = d.arm.map(_detection())
        d["ago_loading"] = d.arm.map(_ago(list(d.arm)))
        d = d.drop(columns=["arms", "n_arms"])

    # ---- COMPOSITION RETENTION: how much of the coupling is CELL-INTRINSIC? (card.py convention) ----
    if dec is not None:
        m = pd.Series(dec["b"], index=dec["cols"])
        sdd = pd.Series(dec["sd"], index=dec["cols"])
        shd, _ = _frac(dec["draws"])
        # The CELL-INTRINSIC identity — the same Shapley/LMG on R², computed on the deconv design. The
        # docstring promises "+deconv → the CELL-INTRINSIC attribution"; before MH-138 that promise was met
        # only by `share_deconv` (= β_f/Σβ), which is not identity for the same additive-value-function
        # reason. Doubles the identity cost (~2.7 → ~5.4 min genome-wide, 8 workers) — cheap for a readout
        # the table has never had.
        phid = AT.shapley_identity(dec["Xz"], dec["m_nnls"], dec["yr"], n_perm=IDENTITY_PERM, seed=seed)
        tphid = phid.sum()
        idd = phid / tphid if tphid > 1e-12 else np.full_like(phid, np.nan)
        _key = "arm" if level == "arm" else "family"
        d["beta_deconv"] = d[_key].map(m)
        d["beta_sd_deconv"] = d[_key].map(sdd)
        d["identity_deconv"] = d[_key].map(pd.Series(idd, index=dec["cols"]))   # ⭐ cell-intrinsic IDENTITY
        d["beta_frac_deconv"] = d[_key].map(pd.Series(shd, index=dec["cols"]))  # cell-intrinsic β-FRACTION
        with np.errstate(divide="ignore", invalid="ignore"):
            d["retention"] = np.where(d.beta.abs() > 1e-9, d.beta_deconv / d.beta, np.nan)
        d["composition_class"] = pd.cut(d.retention, [-np.inf, 0.4, 0.7, np.inf],
                                        labels=["composition_explained", "partial", "cell_intrinsic"])
    return add_reliability(d)          # MH-119: gate the RATIO readouts (share / retention)


_LEVEL = "family"          # set by run(); module-level so forked workers inherit it


def _one(g):
    try:
        return gene_readouts(g, level=_LEVEL)
    except Exception:
        return None


def run(genes: Optional[Sequence[str]] = None, *, workers: int = 8, limit: Optional[int] = None,
        level: str = "family") -> pd.DataFrame:
    global _LEVEL
    _LEVEL = level                                       # forks inherit this
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")                    # parallelise over GENES, not BLAS
    from multiprocessing import Pool
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.evidence import ledger as LG

    LD._load(); LG.pooled_he_edges()                     # warm caches → forks inherit them (copy-on-write)
    if genes is None:
        genes = sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = list(genes)[:limit]
    _one(genes[0])                                       # force lazy paths before the fork

    with Pool(workers) as pool:
        parts = [d for d in pool.imap_unordered(_one, genes, chunksize=8) if d is not None]
    E = pd.concat(parts, ignore_index=True)

    # ---- GENE AGGREGATES (the roll-up) ----
    _u = "arm" if level == "arm" else "family"          # the readout UNIT

    def agg(d: pd.DataFrame) -> pd.Series:
        top_mag = d.loc[d.beta.idxmax()] if len(d) else None
        # ⭐ `top_family_identity` is now the argmax of the REAL identity (Shapley/LMG on R²), not of
        # `share` (MH-138). Before this it was `d.share.idxmax()`, and since share = β_f/Σβ with Σβ a
        # per-gene CONSTANT, argmax(share) ≡ argmax(β) — so the "identity" column was argmax β wearing an
        # identity label, and MEASURED IDENTICAL to `top_family_magnitude` in 99.35% of 1,549 genes (the
        # 0.65% that differed did so only because share is E[β_f/Σβ] per draw, not E[β_f]/E[Σβ] — Jensen,
        # not information). It duplicated magnitude AND argmax β is the one readout MH-124/126 measured
        # AT CHANCE for naming the canonical family. ⚠ These two columns will now genuinely DIFFER.
        top_id = d.loc[d.identity.idxmax()] if d.identity.notna().any() else None
        return pd.Series({
            "n_fam": len(d), "n_arms": int(d.n_arms.sum()) if "n_arms" in d else len(d),
            "total_pressure": float(d.beta.sum()),                       # gene AGGREGATE (Σ_f β_f)
            "n_identified": int(d.identified.sum()),                     # IDENTIFIABILITY
            "frac_identified": float(d.identified.mean()),
            "top_family_magnitude": None if top_mag is None else top_mag[_u],
            "top_beta": float(d.beta.max()),
            "top_family_identity": None if top_id is None else top_id[_u],
            "top_identity": float(d.identity.max()) if d.identity.notna().any() else np.nan,
            "identity_eq_magnitude": (None if (top_id is None or top_mag is None)
                                      else bool(top_id[_u] == top_mag[_u])),   # do WHO and HOW-MUCH agree?
            "top_beta_frac": float(d.beta_frac.max()) if d.beta_frac.notna().any() else np.nan,
            "max_beta_frac_sd": float(d.beta_frac_sd.max()) if d.beta_frac_sd.notna().any() else np.nan,
            "n_discovered": int((d.pip_discovery > 0.5).sum()),          # DISCOVERY
            "n_dense_included": int((d.pip_dense > 0.5).sum()),
            "concentration": float(d.beta.max() / d.beta.sum()) if d.beta.sum() > 0 else np.nan,
            # COMPOSITION RETENTION (the tag the design demands — never silently condition it away)
            "median_retention": float(d.retention.median()) if "retention" in d else np.nan,
            "n_composition_explained": int((d.get("composition_class") == "composition_explained").sum()),
            "n_cell_intrinsic": int((d.get("composition_class") == "cell_intrinsic").sum()),
        })

    G = E.groupby("gene").apply(agg).reset_index()
    # ARM level writes to its OWN files — the canonical §8 family tables are never clobbered.
    oe = OUT_EDGES if level == "family" else OUT_EDGES.with_name("readouts_arm_edges.tsv")
    og = OUT_GENES if level == "family" else OUT_GENES.with_name("readouts_arm_genes.tsv")
    oe.parent.mkdir(parents=True, exist_ok=True)
    E.to_csv(oe, sep="\t", index=False, float_format="%.5f")
    G.to_csv(og, sep="\t", index=False, float_format="%.5f")
    print(f"[readouts:{level}] {len(E)} rows x {E.gene.nunique()} genes -> {oe.name} / {og.name}")
    return E, G


def main() -> None:
    import argparse, time
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--level", choices=["family", "arm"], default="family",
                    help="family = one row per (gene, seed-family); arm = per (gene, arm) [§8 collapse removed]")
    a = ap.parse_args()
    t0 = time.time()
    E, G = run(workers=a.workers, limit=a.limit, level=a.level)
    grain = "family" if a.level == "family" else "arm"
    arm_slots = f" | {int(E.n_arms.sum())} arm-slots" if "n_arms" in E else ""
    print(f"\n=== THE FOUR READOUTS — genome-wide (core C canonical; +deconv retention reported) ===")
    print(f"  {len(G)} genes × {len(E)} (gene,{grain}) rows{arm_slots} | {(time.time()-t0)/60:.1f} min")
    print(f"\n  COUPLING       : median |beta| {E.beta.abs().median():.4f} | beta>0 {100*(E.beta>0).mean():.1f}%")
    print(f"  IDENTIFIABILITY: |z|>2 in {int(E.identified.sum())}/{len(E)} families ({100*E.identified.mean():.1f}%)")
    print(f"                   median beta_frac_sd {E.beta_frac_sd.median():.3f} (wide = genuine non-identifiability)")
    print(f"  ATTRIB/IDENTITY: median top_identity per gene {G.top_identity.median():.3f}"
          f"  | WHO == HOW-MUCH in {G.identity_eq_magnitude.mean():.1%} of genes")
    print(f"  DISCOVERY      : PIP>0.5 in {int((E.pip_discovery>0.5).sum())}/{len(E)} ({100*(E.pip_discovery>0.5).mean():.1f}%)")
    print(f"\n  GENE AGGREGATE : median total_pressure {G.total_pressure.median():.3f} | "
          f"median concentration {G.concentration.median():.3f}")
    if "composition_class" in E:
        vc = E.composition_class.value_counts()
        tot = int(vc.sum())
        print(f"\n  ⚠ COMPOSITION RETENTION (coupling core-C vs +deconv; the tag the design demands):")
        for k in ["cell_intrinsic", "partial", "composition_explained"]:
            n = int(vc.get(k, 0))
            print(f"     {k:22s} {n:5d}  ({100*n/max(tot,1):.1f}%)")
        print(f"     median retention {E.retention.median():.3f}  "
              f"(1.0 = fully cell-intrinsic; <0.4 = composition-explained)")
    print(f"\n  -> {OUT_EDGES}\n  -> {OUT_GENES}")


if __name__ == "__main__":
    main()
