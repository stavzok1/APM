"""THE FOUR READOUTS of the ONE posterior, genome-wide — the canonical per-edge / per-gene table.

`LEARNED_MODEL_DISCOVERY_SYNTHESIS §6b`: the model is **ONE Bayesian posterior over the latent repression field
M, per gene**, read four ways. Those readouts were scattered across method-comparison diagnostics
(`attribution_eb.{full,identifiability_full,selection_full}` compare estimators; they are NOT the readouts) and
per-gene helpers. This module emits them ONCE, together, on ONE confounder block.

Per gene, TWO Gibbs runs on the SAME design (`spike_slab._gibbs_posterior`, the one sampler):
  * **π≡1 DENSE**  → coupling · attribution (magnitude + identity) · identifiability
  * **evidence-π**  → discovery

Read four ways (all per SEED FAMILY — the latent is family-level, §8; arms are mapped back at the end):

  COUPLING            `beta` (posterior mean), `beta_sd`
  ATTRIB / MAGNITUDE  `beta` — the force  (the WHO-HOW-MUCH)
  ATTRIB / IDENTITY   `share` ± `share_sd` — the Bayesian Shapley (MH-94). For the linear aggregate
                      ρ = Σ_f β_f·X_f on z-scored X, the Shapley value of f is EXACTLY β_f, so the per-draw
                      share is `β_f / Σβ` → posterior mean ± sd. Abundance-removed ⇒ IDENTITY, not magnitude.
  IDENTIFIABILITY     `z = beta/beta_sd` (|z|>2 = identified; MH-83: precision 0.86 / recall 0.89 vs the
                      conditioned-partial ground truth) and `share_sd` (WIDE share = genuine non-identifiability
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

N_ITER, BURN = 2000, 700

_AUX: dict = {}


def add_reliability(d: pd.DataFrame) -> pd.DataFrame:
    """MH-119 — GATE THE RATIO READOUTS. `share` and `retention` are both X/Y with a Y that can vanish; a
    vanishing denominator makes the QUESTION ill-posed, not just the estimate noisy.

    **`share` (Bayesian Shapley, β_f/Σβ).** MEASURED mechanism: it explodes ONLY when a gene's βs CANCEL —
    repressive (β>0) and anti-repressive (β<0) units summing to ≈0 while Σ|β| stays large. Of the rows with
    |share|>2, **100%** sit in a gene with ≥1 negative β (only 8% of genes have any), their |Σβ|/Σ|β| is
    **0.172** vs **1.000** for normal rows, and their |Σβ| is 10x smaller. CDH11's βs cancel to −0.0018, so a
    β of +0.072 reads as "999% of the total" — arithmetically right, semantically meaningless. **The model
    already knew: `share_sd` is 316 on those rows vs ~0.04 normally.** It was simply never read.

    **`retention` (β_deconv/β_core).** Same family, different mechanism: a SMALL |β_core| (not cancellation).
    Gate on the core β actually being identified.

    Emits (FLAG, don't delete — the raw values are kept):
      `net_pressure`     |Σβ|/Σ|β| ∈ [0,1] — the gene's SIGN COHERENCE. 1 = all units push the same way.
      `share_abs`        |β_f|/Σ|β_f| ∈ [0,1] — a BOUNDED companion to `share`, always well-defined.
      `share_reliable`   the signed share is meaningful (net_pressure ≥ 0.5 AND share_sd ≤ 1).
      `retention_reliable`  the core β is identified (|z|>2), so the ratio has a real denominator.
    """
    g = d.groupby("gene")["beta"]
    sb = g.transform("sum")
    sa = g.transform(lambda s: s.abs().sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        d["net_pressure"] = (sb.abs() / sa.replace(0, np.nan)).clip(0, 1)
        d["share_abs"] = (d["beta"].abs() / sa.replace(0, np.nan)).clip(0, 1)
    d["share_reliable"] = (d["net_pressure"] >= 0.5) & (d["share_sd"].abs() <= 1.0)
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


def _posteriors(gene, deconv, n_iter, burn, seed, level="family"):
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
    pi = np.clip(AE._evidence_pi(gene, cols), 0.0, 1.0)
    _, _, pip_disc = SS._gibbs_posterior(Xz, yr, pi, n_iter=n_iter, burn=burn, seed=seed)   # evidence-π
    return dict(cols=cols, n=n, p=p, b=b, sd=sd, pip_d=pip_d, draws=draws, pip_disc=pip_disc,
                pi=pi, members=members, fam=fam)


def gene_readouts(gene: str, *, n_iter: int = N_ITER, burn: int = BURN, seed: int = 0,
                  level: str = "family") -> Optional[pd.DataFrame]:
    """The four readouts for one gene → one row per seed family (or per ARM if `level="arm"`), under BOTH C
    blocks + the retention tag."""
    core = _posteriors(gene, False, n_iter, burn, seed, level)   # CANONICAL coupling/identifiability (design)
    if core is None:
        return None
    dec = _posteriors(gene, True, n_iter, burn, seed, level)     # CELL-INTRINSIC attribution (composition removed)

    cols, b, sd, draws = core["cols"], core["b"], core["sd"], core["draws"]

    # ---- ATTRIBUTION / IDENTITY: Bayesian Shapley (MH-94). Linear aggregate ⇒ Shapley(f) = β_f exactly. ----
    def _shapley(dr):
        tot = dr.sum(1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            sh = np.where(np.abs(tot) > 1e-9, dr / tot, np.nan)                     # per-draw share β_f/Σβ
        return np.nanmean(sh, 0), np.nanstd(sh, 0)
    share, share_sd = _shapley(draws)

    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 1e-12, b / sd, 0.0)

    d = pd.DataFrame({
        "gene": gene, "family": cols, "n": core["n"], "p_fam": core["p"],
        "beta": b, "beta_sd": sd,                        # COUPLING (core C — the design's canonical)
        "z": z, "identified": np.abs(z) > 2.0,           # IDENTIFIABILITY (MH-83: prec .86 / rec .89)
        "share": share, "share_sd": share_sd,            # ATTRIBUTION/identity — Bayesian Shapley w/ width
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
        shd, _ = _shapley(dec["draws"])
        _key = "arm" if level == "arm" else "family"
        d["beta_deconv"] = d[_key].map(m)
        d["beta_sd_deconv"] = d[_key].map(sdd)
        d["share_deconv"] = d[_key].map(pd.Series(shd, index=dec["cols"]))      # cell-intrinsic IDENTITY
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
        top_id = d.loc[d.share.idxmax()] if d.share.notna().any() else None
        return pd.Series({
            "n_fam": len(d), "n_arms": int(d.n_arms.sum()) if "n_arms" in d else len(d),
            "total_pressure": float(d.beta.sum()),                       # gene AGGREGATE (Σ_f β_f)
            "n_identified": int(d.identified.sum()),                     # IDENTIFIABILITY
            "frac_identified": float(d.identified.mean()),
            "top_family_magnitude": None if top_mag is None else top_mag[_u],
            "top_beta": float(d.beta.max()),
            "top_family_identity": None if top_id is None else top_id[_u],
            "top_share": float(d.share.max()) if d.share.notna().any() else np.nan,
            "top_share_sd": float(top_id.share_sd) if top_id is not None else np.nan,
            "max_share_sd": float(d.share_sd.max()) if d.share_sd.notna().any() else np.nan,  # widest = least identifiable
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
    a = ap.parse_args()
    t0 = time.time()
    E, G = run(workers=a.workers, limit=a.limit)
    print(f"\n=== THE FOUR READOUTS — genome-wide (core C canonical; +deconv retention reported) ===")
    print(f"  {len(G)} genes × {len(E)} (gene,family) rows | {int(E.n_arms.sum())} arm-slots | {(time.time()-t0)/60:.1f} min")
    print(f"\n  COUPLING       : median |beta| {E.beta.abs().median():.4f} | beta>0 {100*(E.beta>0).mean():.1f}%")
    print(f"  IDENTIFIABILITY: |z|>2 in {int(E.identified.sum())}/{len(E)} families ({100*E.identified.mean():.1f}%)")
    print(f"                   median share_sd {E.share_sd.median():.3f} (wide = genuine non-identifiability)")
    print(f"  ATTRIB/IDENTITY: median top-share per gene {G.top_share.median():.3f} ± {G.top_share_sd.median():.3f}")
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
