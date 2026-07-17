"""EB posterior vs bagged-NNLS for ATTRIBUTION (2026-07-07). Does a per-gene learned-τ² Bayesian shrinkage
give the same stable per-edge weights as the canonical bagged NNLS — and hand us calibrated uncertainty
(posterior sd) for free?

Both estimators read the SAME object the attribution layer needs: the per-seed-family coefficient M_fam,
on **identical** inputs — C-residualised, z-scored family predictors with the `sd<0.1` variance floor
(exactly `states._bagged_nnls`). They differ only in HOW they stabilise the collinearity-unstable single
fit:
    bagged NNLS  — resample-and-average (variance ↓ ~1/B); sd = across-bootstrap spread.
    EB posterior — half-normal ridge with the slab variance ν² LEARNED (spike-and-slab Gibbs, π≡1);
                   mean = posterior mean, sd = posterior sd → the |mean/sd|>2 'identified edge' rule (§10).

Compared on the criteria that actually matter for attribution:
  (1) AGREEMENT   — do the two mean-weight vectors rank families the same? (Spearman per gene)
  (2) REPRODUCIBILITY — mean weights across two seeds (the property single-lasso fails at, corr 0.03)
  (3) DRIVER RECOVERY — is the top family the curated driver? (PTEN miR-21/17~92, ZEB1 miR-200/429, …)
  (4) UNCERTAINTY vs §9 — does EB |mean/sd| stay LOW on the family-only/unidentified genes (ESR1
      miR-221/222) and HIGH on single-driver genes? (the payoff NNLS-bagging only approximates)

CLI: `.venv/bin/python3 -m mirna_hallmark.learned.attribution_eb`
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import spike_slab as SS
from mirna_hallmark.learned import states as ST

HUB_PANEL = ["PTEN", "GATA3", "ESR1", "ZEB1", "CDKN1A"]
# curated driver families (seed-family label as FAM.family_of renders it) + §9 identifiability class
CURATED = {   # canonical driver families (as FAM renders the seed label) + expected §9 identifiability
    "PTEN": {"drivers": ["miR-21-5p/590-5p", "miR-103-3p/107", "miR-181-5p", "miR-182-5p"], "class": "multi-driver"},
    "ZEB1": {"drivers": ["miR-200bc-3p/429", "miR-141-3p/200a-3p"], "class": "co-drivers(§9 hard split)"},
    "ESR1": {"drivers": ["miR-18-5p"], "class": "miR-18 owner; 221/222 §9-unidentified"},
    "GATA3": {"drivers": ["miR-27-3p"], "class": "single"},
    "CDKN1A": {"drivers": ["miR-15-5p/16-5p/195-5p/424-5p/497-5p", "miR-17-5p/20-5p/93-5p/106-5p/519-3p"], "class": "multi"},
}


def _resid(v, Cm):
    return v - LinearRegression().fit(Cm, v).predict(Cm)


def _gene_family_data(gene):
    """(Y, Xf family predictors, C) exactly as canonical_M assembles them (tumour '01', deconv, HE support)."""
    he = D.high_evidence_edges()
    arms = pd.Index(sorted(set(he.loc[he["gene"] == gene, "miRNA"])))
    if arms.empty:
        return None
    d = ST._state_family_data(gene, "01", FAM.family_of(arms))
    return d                                                      # (Y, Xf, C) or None


def _prep(Y, Xf, C):
    """C-residualised target r and z-scored, variance-floored family predictors — matches `_bagged_nnls`."""
    Cm = C.to_numpy(float)
    yr = -_resid(Y.to_numpy(float), Cm)
    sd = Xf.std(ddof=0)
    Xz = ((Xf - Xf.mean()) / (sd + 1e-9)).fillna(0.0)
    Xz.loc[:, sd < 0.1] = 0.0
    return yr, Xz.to_numpy(float), list(Xf.columns)


def _bagged_nnls_meansd(Xz, yr, *, n_boot=40, seed=0):
    n = len(yr); rng = np.random.default_rng(seed)
    boots = np.zeros((n_boot, Xz.shape[1]))
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[b], _ = nnls(Xz[idx], yr[idx])
    return boots.mean(0), boots.std(0)                            # bootstrap mean, bootstrap sd


def _eb_posterior(Xz, yr, *, n_iter=2000, burn=700, seed=0):
    """Non-negative learned-τ² ridge posterior (half-normal slab, ν² sampled, π≡1). Returns (mean, sd) of β."""
    m, s, _ = SS._gibbs_posterior(Xz, yr, np.ones(Xz.shape[1]), n_iter=n_iter, burn=burn, seed=seed)
    return m, s


def _blasso_posterior(Xz, yr, *, n_iter=2000, burn=700, seed=0):
    """Non-negative Bayesian-lasso posterior (Laplace prior, per-edge τ_j² + learned λ²). Returns (mean, sd).
    Harder shrinkage-toward-0 than the ridge → should give a CLEANER point estimate on un-informed edges."""
    return SS._gibbs_blasso(Xz, yr, n_iter=n_iter, burn=burn, seed=seed)


def _ss_posterior(Xz, yr, pi, *, n_iter=2000, burn=700, seed=0):
    """SPIKE-AND-SLAB posterior with evidence-graded inclusion π (genuine selection). Returns (mean, sd, PIP).
    The inclusion indicator should ZERO un-informed edges (PIP→0) → clean raw mean, fixing the ridge bias."""
    return SS._gibbs_posterior(Xz, yr, np.clip(pi, 0.0, 1.0), n_iter=n_iter, burn=burn, seed=seed)


def _evidence_pi(gene, cols):
    """Evidence-graded inclusion prior π per family, aligned to `cols` (from priors.edge_priors)."""
    from mirna_hallmark.learned import priors as PR
    try:
        pr = PR.edge_priors(gene, by="family")               # no-evidence families → LOW π (0.05), not base rate,
        return pr["pi"].reindex(cols).fillna(0.05).to_numpy(float)  # else they leak into the raw mean
    except Exception:
        return np.full(len(cols), 0.05)


def compare_gene(gene: str) -> dict:
    d = _gene_family_data(gene)
    if d is None:
        return {"gene": gene, "error": "no HE support"}
    yr, Xz, cols = _prep(*d)
    pi = _evidence_pi(gene, cols)
    nn_m, nn_s = _bagged_nnls_meansd(Xz, yr, seed=0)
    eb_m, eb_s = _eb_posterior(Xz, yr, seed=0)
    sb_m, sb_s, sb_pip = SS._gibbs_ss_blasso(Xz, yr, pi, n_iter=2000, burn=700, seed=0)   # spike-slab-LASSO
    sb_m2, _, _ = SS._gibbs_ss_blasso(Xz, yr, pi, n_iter=2000, burn=700, seed=1)          # reproducibility
    tab = pd.DataFrame({"family": cols, "pi": pi, "nnls_mean": nn_m,
                        "eb_mean": eb_m, "eb_sd": eb_s, "sb_mean": sb_m, "sb_sd": sb_s, "sb_pip": sb_pip})
    tab["eb_z"] = tab["eb_mean"] / (tab["eb_sd"] + 1e-9)
    # KEY: does the spike-slab-LASSO (inclusion + Laplace slab) ZERO un-informed edges → clean raw mean +
    # PIP identifiability that beats the ridge |z|>2 on the hard co-driver genes?
    def _corr(a, b):
        return float(spearmanr(a, b).correlation) if np.std(a) > 0 and np.std(b) > 0 else np.nan
    cur = CURATED.get(gene, {})
    drivers = set(cur.get("drivers", []))
    ident_eb = set(tab.loc[tab["eb_z"] > 2, "family"])          # ridge identifiability = |z|>2
    ident_sb = set(tab.loc[tab["sb_pip"] > 0.5, "family"])      # spike-slab-lasso identifiability = PIP>0.5
    top_eb_rawmean = cols[int(np.argmax(eb_m))]
    top_sb_rawmean = cols[int(np.argmax(sb_m))]
    tab = tab.sort_values("sb_mean", ascending=False).reset_index(drop=True)
    return {"gene": gene, "n_fam": len(cols), "class": cur.get("class", ""),
            "agree_nnls_sb": round(_corr(nn_m, sb_m), 3), "repro_sb": round(_corr(sb_m, sb_m2), 3),
            "n_pip>0.5": int((sb_pip > 0.5).sum()),
            "eb_rawtop_driver": (top_eb_rawmean in drivers) if drivers else None,
            "sb_rawtop_driver": (top_sb_rawmean in drivers) if drivers else None,
            "rec_eb_z": f"{len(drivers & ident_eb)}/{len(drivers)}" if drivers else "",
            "rec_sb_pip": f"{len(drivers & ident_sb)}/{len(drivers)}" if drivers else "",
            "_table": tab}


def run(genes=None):
    genes = genes or HUB_PANEL
    rows = []
    for g in genes:
        r = compare_gene(g)
        rows.append(r)
        if "_table" in r:
            print(f"\n=== {g}  ({r['class']}, {r['n_fam']} families) — top-4 by spike-slab-lasso mean ===")
            print(r["_table"].head(4).to_string(index=False,
                  formatters={c: (lambda v: f"{v:.3f}") for c in
                              ["pi", "nnls_mean", "eb_mean", "eb_z", "sb_mean", "sb_sd", "sb_pip"]},
                  columns=["family", "pi", "nnls_mean", "eb_mean", "eb_z", "sb_mean", "sb_pip"]))
    df = pd.DataFrame([{k: v for k, v in r.items() if k != "_table"} for r in rows])
    print("\n=== SUMMARY: spike-slab-LASSO (inclusion+Laplace) vs EB-ridge vs bagged NNLS (attribution) ===")
    with pd.option_context("display.width", 190):
        print(df.to_string(index=False))
    print("\nKEY — does inclusion+Laplace fix the un-informed-edge RAW-mean bias?")
    print(f"  EB-ridge raw-top is a driver:        {sum(bool(r.get('eb_rawtop_driver')) for r in rows)}/"
          f"{sum(r.get('eb_rawtop_driver') is not None for r in rows)}")
    print(f"  spike-slab-lasso raw-top is a driver: {sum(bool(r.get('sb_rawtop_driver')) for r in rows)}/"
          f"{sum(r.get('sb_rawtop_driver') is not None for r in rows)}")
    print(f"  driver recovery — ridge |z|>2 {[r.get('rec_eb_z') for r in rows if r.get('rec_eb_z')]}"
          f"  vs  ss-lasso PIP>0.5 {[r.get('rec_sb_pip') for r in rows if r.get('rec_sb_pip')]}")
    return df


def full(*, out: str = "mirna_hallmark/output/learned/attribution_eb_full.tsv",
         limit: int | None = None, progress: int = 50) -> pd.DataFrame:
    """GENOME-WIDE EB-ridge-posterior vs bagged-NNLS attribution agreement + reproducibility (the scale check
    that justifies 'one model for attribution too'). Per multi-family gene: agreement ρ(NNLS mean, EB mean),
    reproducibility of each across two seeds, and the EB |z|>2 identified-set size vs the spike-slab PIP>0.5
    set (do the two identifiability signals concur without curated labels?). Lean: NNLS + EB-ridge + spike-π."""
    from pathlib import Path
    from mirna_hallmark.learned import priors as PR  # noqa
    from mirna_hallmark.learned.evidence import ledger as _LG
    genes = sorted(set(_LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = genes[:limit]
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[attribution_eb full] {i}/{len(genes)}", flush=True)
        try:
            d = _gene_family_data(g)
            if d is None or d[1].shape[1] < 2:
                continue
            yr, Xz, cols = _prep(*d)
            pi = _evidence_pi(g, cols)
            nn_m, _ = _bagged_nnls_meansd(Xz, yr, seed=0)
            nn_m2, _ = _bagged_nnls_meansd(Xz, yr, seed=1)
            eb_m, eb_s, _ = SS._gibbs_posterior(Xz, yr, np.ones(len(cols)), n_iter=1200, burn=400, seed=0)
            eb_m2, _, _ = SS._gibbs_posterior(Xz, yr, np.ones(len(cols)), n_iter=1200, burn=400, seed=1)
            _, _, pip = SS._gibbs_posterior(Xz, yr, np.clip(pi, 0, 1), n_iter=1200, burn=400, seed=0)
            ez = eb_m / (eb_s + 1e-9)

            def _c(a, b):
                return float(spearmanr(a, b).correlation) if np.std(a) > 0 and np.std(b) > 0 else np.nan
            id_z = set(np.where(ez > 2)[0]); id_p = set(np.where(pip > 0.5)[0])
            jac = len(id_z & id_p) / len(id_z | id_p) if (id_z | id_p) else np.nan
            rows.append({"gene": g, "n_fam": len(cols),
                         "agree_nnls_eb": _c(nn_m, eb_m), "repro_nnls": _c(nn_m, nn_m2), "repro_eb": _c(eb_m, eb_m2),
                         "n_ident_z": len(id_z), "n_ident_pip": len(id_p), "ident_jaccard": jac})
        except Exception:
            pass
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    d = df[df["agree_nnls_eb"].notna()]
    print(f"\n=== GENOME-WIDE ATTRIBUTION: EB-ridge vs bagged NNLS ({len(d)} multi-family genes) → {out} ===")
    print(f"agreement ρ(NNLS,EB): mean {d['agree_nnls_eb'].mean():.3f} median {d['agree_nnls_eb'].median():.3f} "
          f"| ≥0.7 on {100*(d['agree_nnls_eb']>=0.7).mean():.0f}% of genes")
    print(f"reproducibility: NNLS {d['repro_nnls'].mean():.3f}  EB {d['repro_eb'].mean():.3f} (single-lasso ≈0.03)")
    print(f"identifiability |z|>2 vs PIP>0.5: mean Jaccard {d['ident_jaccard'].mean():.3f} "
          f"| median idents z={int(d['n_ident_z'].median())} pip={int(d['n_ident_pip'].median())}")
    return df


def _partial_identified(Y, Xf, C, thr: float = -0.1):
    """Reference family-level identifiability (§7/§9 conditioned-partial): each family's anti-corr NET of the
    OTHER families + C. Identified = partial-Spearman < thr (the data resolves it). Returns boolean array."""
    Cm = C.to_numpy(float); yv = Y.to_numpy(float); Xv = Xf.to_numpy(float); p = Xv.shape[1]
    out = np.zeros(p, dtype=bool)
    for j in range(p):
        others = [k for k in range(p) if k != j]
        Z = np.column_stack([Cm, Xv[:, others]]) if others else Cm
        xr = _resid(Xv[:, j], Z)
        if np.std(xr) < 1e-6:
            continue                                                   # collinear → no independent variation → unidentified
        rho = spearmanr(xr, _resid(yv, Z)).correlation
        out[j] = (rho == rho) and rho < thr
    return out


def identifiability_full(*, out: str = "mirna_hallmark/output/learned/identifiability_vs_partial.tsv",
                         limit: int | None = None, progress: int = 50) -> pd.DataFrame:
    """GENOME-WIDE: do the Bayes identifiability signals (dense-ridge |z|>2 · spike-slab PIP>0.5) agree with
    the established conditioned-partial (§9-style) verdict? Per gene, precision/recall of each Bayes signal
    against the partial-identified set. The scale check the hub-only §9 concordance was missing."""
    from pathlib import Path
    from mirna_hallmark.learned.evidence import ledger as _LG
    genes = sorted(set(_LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = genes[:limit]
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[identifiability_full] {i}/{len(genes)}", flush=True)
        try:
            d = _gene_family_data(g)
            if d is None or d[1].shape[1] < 2:
                continue
            Y, Xf, C = d
            yr, Xz, cols = _prep(Y, Xf, C)
            pi = _evidence_pi(g, cols)
            eb_m, eb_s, _ = SS._gibbs_posterior(Xz, yr, np.ones(len(cols)), n_iter=1200, burn=400, seed=0)
            _, _, pip = SS._gibbs_posterior(Xz, yr, np.clip(pi, 0, 1), n_iter=1200, burn=400, seed=0)
            ref = _partial_identified(Y, Xf, C)                       # §9-style reference
            idz = (eb_m / (eb_s + 1e-9)) > 2
            idp = pip > 0.5
            nref = int(ref.sum())

            def _pr(sig):
                tp = int((sig & ref).sum()); ns = int(sig.sum())
                return (tp / ns if ns else np.nan, tp / nref if nref else np.nan)  # precision, recall
            pz, rz = _pr(idz); pp, rp = _pr(idp)
            rows.append({"gene": g, "n_fam": len(cols), "n_ref": nref,
                         "z_prec": pz, "z_rec": rz, "pip_prec": pp, "pip_rec": rp,
                         "n_z": int(idz.sum()), "n_pip": int(idp.sum())})
        except Exception:
            pass
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    d = df[df["n_ref"] > 0]
    print(f"\n=== IDENTIFIABILITY vs conditioned-partial (§9-style), {len(d)} genes with ≥1 partial-identified ===")
    print(f"  |z|>2 : precision {d['z_prec'].mean():.2f}  recall {d['z_rec'].mean():.2f}  "
          f"(does the DATA-signal match the data reference)")
    print(f"  PIP>.5: precision {d['pip_prec'].mean():.2f}  recall {d['pip_rec'].mean():.2f}  "
          f"(evidence-informed → recall extras the partial misses)")
    return df


def selection_full(*, out: str = "mirna_hallmark/output/learned/selection_lasso_vs_bayes.tsv",
                   limit: int | None = None, progress: int = 50) -> pd.DataFrame:
    """GENOME-WIDE cross-estimator: which FAMILIES survive in the LASSO (nonzero M, prior-guided selection)
    vs the BAYES (|z|>2 significance · PIP>0.5 inclusion), on the SAME family support. Quantitative overlap /
    lasso-only / bayes-only per gene → the concrete 'which survive and which not' the shuffle demo couldn't give."""
    from pathlib import Path
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import regression as LR
    from mirna_hallmark.learned.evidence import ledger as _LG
    genes = sorted(set(_LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = genes[:limit]
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[selection_full] {i}/{len(genes)}", flush=True)
        try:
            Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger", deconv=True)
            X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
            cols = list(X.columns)
            if len(cols) < 2:
                continue
            Ml = LR.fit_gene(Y, X, C, w)
            L = set(np.array(cols)[Ml.reindex(cols).fillna(0).to_numpy() > 0])
            yr, Xz, _ = _prep(Y, X, C)
            eb_m, eb_s, _ = SS._gibbs_posterior(Xz, yr, np.ones(len(cols)), n_iter=1200, burn=400, seed=0)
            _, _, pip = SS._gibbs_posterior(Xz, yr, np.clip(_evidence_pi(g, cols), 0, 1), n_iter=1200, burn=400, seed=0)
            Bz = set(np.array(cols)[(eb_m / (eb_s + 1e-9)) > 2])
            Bp = set(np.array(cols)[pip > 0.5])
            j = lambda a, b: len(a & b) / len(a | b) if (a | b) else np.nan
            rows.append({"gene": g, "n_fam": len(cols), "n_lasso": len(L), "n_bz": len(Bz), "n_pip": len(Bp),
                         "lasso∩bz": len(L & Bz), "lasso_only_vs_bz": len(L - Bz), "bz_only_vs_lasso": len(Bz - L),
                         "jac_lasso_bz": j(L, Bz), "jac_lasso_pip": j(L, Bp), "lasso∩pip": len(L & Bp)})
        except Exception as e:
            if i < 3:
                print(f"  [skip {g}] {e!r}", flush=True)
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    if "jac_lasso_bz" not in df:
        print("no rows collected"); return df
    d = df.dropna(subset=["jac_lasso_bz"])
    print(f"\n=== LASSO selection vs BAYES identification, {len(d)} multi-family genes ===")
    print(f"  per-gene sizes: lasso-nonzero {d['n_lasso'].mean():.1f} | Bayes |z|>2 {d['n_bz'].mean():.1f} | "
          f"Bayes PIP>.5 {d['n_pip'].mean():.1f} | candidates {d['n_fam'].mean():.1f}")
    print(f"  lasso vs |z|>2 : Jaccard {d['jac_lasso_bz'].mean():.2f} | lasso-only {d['lasso_only_vs_bz'].sum()} "
          f"| Bayes-only {d['bz_only_vs_lasso'].sum()} | overlap {d['lasso∩bz'].sum()}")
    print(f"  lasso vs PIP>.5: Jaccard {d['jac_lasso_pip'].mean():.2f} | overlap {d['lasso∩pip'].sum()}")
    print(f"  ⇒ the lasso SELECTS MORE (prior-permissive); Bayes |z|>2 is the significance-gated subset "
          f"({100*d['lasso∩bz'].sum()/max(d['n_bz'].sum(),1):.0f}% of Bayes calls are lasso-selected)")
    return df


if __name__ == "__main__":
    import sys as _sys
    if "--selection" in _sys.argv:
        lim = int(_sys.argv[_sys.argv.index("--limit") + 1]) if "--limit" in _sys.argv else None
        selection_full(limit=lim); _sys.exit(0)
    if "--identifiability" in _sys.argv:
        lim = int(_sys.argv[_sys.argv.index("--limit") + 1]) if "--limit" in _sys.argv else None
        identifiability_full(limit=lim)
    elif "--full" in _sys.argv:
        lim = int(_sys.argv[_sys.argv.index("--limit") + 1]) if "--limit" in _sys.argv else None
        full(limit=lim)
    else:
        run()
