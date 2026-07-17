"""Phase-1 MVP driver (Design §6 Phase 1 + §5 Bars 1–2).

For a hub-gene panel: fit the gene-focused NN adaptive-lasso and run the **out-of-fold gate** —
does the learned aggregate ρ = X·M anti-correlate with held-out target mRNA *better than raw arm
abundance* (partial-Spearman adjusting C)? Model repression ⇒ ρ↑ → Y↓, so a good rho is NEGATIVE and
`rho_model` should be MORE negative than `rho_abund`.

    .venv/bin/python3 -m mirna_hallmark.learned.mvp
    .venv/bin/python3 -m mirna_hallmark.learned.mvp --genes PTEN ESR1 --alpha 0.005 --folds 5

This is the cheapest diagnostic, not the full model. It nests the heuristic (Design §0): freeze the
prior and it reduces to curated fixed-M pressure. Gate to proceed to Phase 2 (Design §6): recovers
curated drivers AND beats raw-abundance coupling out-of-fold on the panel.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import regression as LR

HUB_PANEL = ["PTEN", "GATA3", "ESR1", "ZEB1", "CDKN1A"]


def _resid(v: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    return v - LinearRegression().fit(Cmat, v).predict(Cmat)


def _mean_jaccard(sets: list) -> float:
    """Mean pairwise Jaccard of the selected (nonzero) predictor sets across folds — selection stability."""
    js = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            js.append(len(sets[i] & sets[j]) / len(u) if u else 1.0)
    return float(np.mean(js)) if js else float("nan")


def oof_gate(gene: str, *, alpha: float = 0.005, folds: int = 5, seed: int = 0,
             w_prior_source: str = "evidence_score", family: bool = False, deconv: bool = False,
             n_latent: int = 0, occ: bool = False, kd: bool = False, n_tf: int = 0,
             tf_activity: bool = False, permute_prior: int | None = None, batch: bool = False,
             edges=None) -> dict:
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source=w_prior_source, deconv=deconv,
                                  n_latent=n_latent, n_tf=n_tf, tf_activity=tf_activity, edges=edges)
    if family:                                                 # collapse to seed-family predictors (Design §F)
        fam = FAM.family_of(X.columns)
        X, w, _members = FAM.collapse_by_family(X, w, fam)
    if permute_prior is not None:                              # Bar 3 null: shuffle prior→edge map (keep marginal)
        rng = np.random.default_rng(permute_prior)
        w = pd.Series(rng.permutation(w.to_numpy()), index=w.index)
    raw_abund = X.mean()                                       # pre-transform abundance (for the gauge metric)
    if occ:                                                    # occupancy link + threshold gauge (Design §C)
        from mirna_hallmark.learned import occupancy as OCC
        kappa = None
        if kd:                                                 # scanMiR affinity-aware per-edge κ (arm mode)
            from mirna_hallmark.learned import kd as KD
            kappa = KD.edge_kappa(X, gene)                     # family labels fall back to κ0 (naive)
        X = OCC.transform(X, kappa=kappa)
    n = len(Y)
    oof_model = np.full(n, np.nan)
    oof_abund = np.full(n, np.nan)
    sel_sets: list = []
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    for tr, te in kf.split(X):
        M = LR.fit_gene(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w, alpha=alpha)
        sel_sets.append(frozenset(M[M > 0].index))            # which predictors survive this fold
        oof_model[te] = LR.aggregate(X.iloc[te], M)
        oof_abund[te] = X.iloc[te].to_numpy(dtype=float).mean(axis=1)  # raw mean abundance baseline

    Cmat = C.to_numpy(dtype=float)
    if batch:                                                  # POST-FIT batch conditioning (plate_both) —
        from mirna_hallmark import tcga_batch as TB           # added to the coupling test ONLY, never the
        bd = TB.tcga_batch_dummies(list(Y.index), kind="plate_both").reindex(Y.index).fillna(0.0)
        bd = bd.loc[:, bd.sum(axis=0) > 0]                    # lasso fit (in-lasso strips heterogeneity)
        Cmat = np.column_stack([Cmat, bd.to_numpy(dtype=float)])
    yr = _resid(Y.to_numpy(dtype=float), Cmat)
    rho_model = spearmanr(_resid(oof_model, Cmat), yr).correlation      # expect < 0 (pressure ↓ target)
    rho_abund = spearmanr(_resid(oof_abund, Cmat), yr).correlation
    # Bar 2 (Design §0 ref 2): the CURATED fixed-M aggregate = X·w_prior with the weights FROZEN at the
    # prior (the nesting limit — no fitting), over ALL regulators. Does LEARNING the magnitudes beat, or
    # match with fewer arms, freezing them?
    w_vec = w.reindex(X.columns).fillna(0.0).to_numpy(dtype=float)
    rho_curated = spearmanr(_resid(X.to_numpy(dtype=float) @ w_vec, Cmat), yr).correlation

    M_full = LR.fit_gene(Y, X, C, w, alpha=alpha)
    top = M_full[M_full > 0].sort_values(ascending=False)
    # gauge metric (Design §C): share of the aggregate held by the single most-abundant arm. Linear a·w
    # lets the most-abundant arm dominate; the occupancy link saturates it → lower domination = D(m) retired.
    c = (M_full.clip(lower=0) * X.mean()).clip(lower=0)
    abund_dom = float(c.get(raw_abund.idxmax(), 0.0) / c.sum()) if c.sum() > 0 else float("nan")
    return {
        "gene": gene, "unit": "family" if family else "arm",
        "n": n, "n_pred": X.shape[1], "nonzero": int((M_full > 0).sum()), "abund_dom": abund_dom,
        "rho_abund": rho_abund, "rho_curated": rho_curated, "rho_model": rho_model,
        "vs_abund": bool(rho_model < rho_abund),
        "vs_curated": bool(rho_model <= rho_curated + 1e-6),   # learned ≥ curated coupling (sparser)
        "stability": _mean_jaccard(sel_sets),                 # cross-fold selection Jaccard (↑ = more stable)
        "top": ", ".join(f"{m}={v:.2f}" for m, v in top.head(4).items()),
    }


def run(genes=None, *, alpha: float = 0.005, folds: int = 5,
        w_prior_source: str = "evidence_score", family: bool = False, deconv: bool = False,
        n_latent: int = 0, occ: bool = False, kd: bool = False) -> pd.DataFrame:
    genes = genes or HUB_PANEL
    rows = []
    for g in genes:
        try:
            rows.append(oof_gate(g, alpha=alpha, folds=folds, w_prior_source=w_prior_source,
                                 family=family, deconv=deconv, n_latent=n_latent, occ=occ, kd=kd))
        except Exception as e:  # keep the panel going; report the failure
            rows.append({"gene": g, "error": repr(e)[:80]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 170, "display.max_colwidth", 62):
        print(df.to_string(index=False))
    if "vs_abund" in df:
        print(f"\ngate: learned beats raw-abundance {int(df['vs_abund'].sum())}/{df['vs_abund'].notna().sum()}"
              f" | matches/beats curated fixed-M {int(df['vs_curated'].sum())}/{df['vs_curated'].notna().sum()}"
              + (f" | mean stability {df['stability'].mean():.2f}" if "stability" in df else ""))
    return df


def gate_fdr(genes=None, *, alpha: float = 0.005, family: bool = True, deconv: bool = False,
             w_prior_source: str = "ledger", out: str = "mirna_hallmark/output/learned/gate_fdr.tsv",
             progress: int = 100) -> pd.DataFrame:
    """Genome-wide gate WITH multiplicity control (user 2026-07-05: full numbers + FDR, precursor-style).
    Runs the OOF gate on ALL HE genes; per gene, a one-sided analytic p that rho_model < 0 (real repression
    coupling; partial-t on df≈n−8). Then **BH** and **Benjamini-Yekutieli** q across genes — BY because the
    gene tests are DEPENDENT through shared regulators (the seed-family dependence the precursor's family-FDR
    is built for; `coupling_inference.benjamini_yekutieli`). Reports, among FDR-significant genes, how many the
    learned model BEATS raw abundance / matches curated. The honest full-scale replacement for the panel gate."""
    from pathlib import Path
    from scipy import stats as _st
    from mirna_hallmark import data_loaders as _D
    from mirna_hallmark import stats as _S
    from mirna_hallmark.coupling_inference import benjamini_yekutieli
    from mirna_hallmark.learned.evidence import ledger as _LG
    genes = genes or sorted(set(_LG.pooled_he_edges()["gene"].dropna().astype(str)))  # POOLED-HE gene universe (migration)
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[gate_fdr] {i}/{len(genes)}", flush=True)
        try:
            rows.append(oof_gate(g, alpha=alpha, family=family, deconv=deconv, w_prior_source=w_prior_source))
        except Exception:
            pass
    df = pd.DataFrame(rows)
    df = df[df["rho_model"].notna() & df["n"].notna()].copy()

    def _p_neg(rho, n):                                        # one-sided p for rho<0, partial-t, df≈n−8 (C+aggregate)
        dfree = int(n) - 8
        if not (rho == rho) or dfree <= 1:
            return np.nan
        if rho >= 0:
            return 1.0
        t = rho * np.sqrt(dfree / max(1.0 - rho * rho, 1e-9))
        return float(_st.t.cdf(t, dfree))                     # left tail for negative t

    df["p_coupling"] = [_p_neg(r.rho_model, r.n) for r in df.itertuples()]
    df["q_bh"] = _S.bh_fdr(df["p_coupling"].values)
    df["q_by"] = benjamini_yekutieli(df["p_coupling"].values)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    sig = df[df["q_by"] < 0.05]
    print(f"\n=== GENOME-WIDE GATE + FDR ({len(df)} HE genes fit) ===")
    print(f"FDR-significant coupling (q_BY<0.05 & rho_model<0): {len(sig)} ({100*len(sig)/max(len(df),1):.0f}%) "
          f"| q_BH<0.05: {int((df['q_bh']<0.05).sum())}")
    if len(sig):
        print(f"among FDR-significant genes: beats raw abundance {int((sig['rho_model']<sig['rho_abund']).sum())}/{len(sig)} "
              f"| ≥ curated {int((sig['rho_model']<=sig['rho_curated']+1e-6).sum())}/{len(sig)}")
    print(f"all fit: beats abundance {int(df['vs_abund'].sum())}/{len(df)} ({100*df['vs_abund'].mean():.0f}%) "
          f"| mean rho_model {df['rho_model'].mean():.3f} | median n {int(df['n'].median())}  → wrote {out}")
    return df


if __name__ == "__main__":
    import sys as _sys
    if "--fdr" in _sys.argv:
        gate_fdr(); _sys.exit(0)
    ap = argparse.ArgumentParser(description="Phase-1 learned-model MVP: OOF gate vs raw abundance")
    ap.add_argument("--genes", nargs="*", default=None, help="gene symbols (default: hub panel)")
    ap.add_argument("--alpha", type=float, default=0.005, help="lasso penalty")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--prior", choices=["evidence_score", "ledger", "ledger_mrna", "scanmir", "fused"],
                    default="evidence_score",
                    help="adaptive-penalty prior: hand-weighted evidence_score | deduped method-centric ledger "
                         "| ledger_mrna (transfection-calibrated) | scanmir (biochemical K_D affinity) "
                         "| fused (ledger existence + scanMiR magnitude)")
    ap.add_argument("--family", action="store_true", help="collapse regulators to seed-family predictors (Design §F)")
    ap.add_argument("--deconv", action="store_true", help="add CIBERSORTx non-malignant composition to C (Design §B)")
    ap.add_argument("--latent", type=int, default=0, help="add K transcription-state latent factors to C (Design §B)")
    ap.add_argument("--occ", action="store_true", help="occupancy link + threshold gauge on abundance (Design §C)")
    ap.add_argument("--kd", action="store_true", help="scanMiR affinity-aware per-edge κ for --occ (Design §C; arm mode)")
    a = ap.parse_args()
    run(a.genes, alpha=a.alpha, folds=a.folds, w_prior_source=a.prior, family=a.family,
        deconv=a.deconv, n_latent=a.latent, occ=a.occ, kd=a.kd)
