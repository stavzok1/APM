"""E — TarBase-v9 DIRECTION channel (Design §Decision E).

TarBase v9 records, per validated interaction, a `regulation` label (Negative = repressive / Positive =
activating). A count-based HE set (miRTarBase-derived) admits an edge on *evidence of interaction* without
its sign — so it can wrongly include an edge that is validated-**activating** (miRNA UP-regulates the target),
which contradicts the non-negative repression model. This module scores each edge's net direction and flags
the validated-activating HE edges for exclusion, then checks their TCGA coupling (an activating edge should
NOT show repressive coupling).

CLI: `python -m mirna_hallmark.learned.evidence.direction`
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import data_loaders as D

_TARBASE = Path("data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz")


def tarbase_direction(*, min_total: int = 1) -> pd.DataFrame:
    """Per (arm, gene): distinct-PMID counts of Negative vs Positive interactions + a net-direction class."""
    t = pd.read_csv(_TARBASE, sep="\t", low_memory=False, usecols=[
        "mirna_name", "gene_name", "regulation", "article_pubmed_id"])
    t = t[t["regulation"].isin(["Negative", "Positive"])].dropna(subset=["mirna_name", "gene_name"])
    # count DISTINCT PMIDs per (edge, direction) so one heavily-reported study doesn't dominate
    g = (t.drop_duplicates(["mirna_name", "gene_name", "regulation", "article_pubmed_id"])
           .groupby(["mirna_name", "gene_name", "regulation"]).size().unstack(fill_value=0))
    g = g.rename(columns={"Negative": "n_neg", "Positive": "n_pos"})
    for c in ("n_neg", "n_pos"):
        if c not in g:
            g[c] = 0
    g = g.reset_index().rename(columns={"mirna_name": "arm", "gene_name": "gene"})
    g = g[(g["n_neg"] + g["n_pos"]) >= min_total]
    g["direction"] = np.where(g["n_pos"] > g["n_neg"], "activating",
                              np.where(g["n_neg"] > g["n_pos"], "repressive", "mixed"))
    return g


def he_direction_audit() -> pd.DataFrame:
    """Map TarBase-v9 direction onto the HE edges → how many are validated-ACTIVATING (exclude candidates)?"""
    he = D.high_evidence_edges()[["miRNA", "gene"]].drop_duplicates()
    d = tarbase_direction()
    m = he.merge(d, left_on=["miRNA", "gene"], right_on=["arm", "gene"], how="left")
    cov = m[m["direction"].notna()]
    print(f"HE edges: {len(he)} | TarBase-v9 covered: {len(cov)} ({100*len(cov)/len(he):.0f}%)")
    print(cov["direction"].value_counts().to_string())
    act = cov[cov["direction"] == "activating"]
    print(f"\nvalidated-ACTIVATING HE edges (drop candidates): {len(act)}  "
          f"(net-positive; {int((act['n_neg']==0).sum())} with ZERO negative evidence)")
    with pd.option_context("display.width", 160):
        print(act.sort_values("n_pos", ascending=False).head(15)[["miRNA", "gene", "n_neg", "n_pos"]].to_string(index=False))
    return m


def coupling_by_direction(sample: int = 400) -> None:
    """Do the flagged-activating HE edges actually show non-repressive TCGA coupling? (validates the channel)."""
    from scipy.stats import spearmanr, mannwhitneyu
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import data as LD
    m = D.high_evidence_edges()[["miRNA", "gene"]].drop_duplicates().merge(
        tarbase_direction(), left_on=["miRNA", "gene"], right_on=["arm", "gene"], how="left")
    m = m[m["direction"].isin(["activating", "repressive"])]
    Xall, Yall = LD._load()["X"], LD._load()["Y"]
    Cbase = LD._load()["C"]
    conf = [c for c in LD.confounder_columns() if c in Cbase.columns]

    def _coup(arm, gene):
        if arm not in Xall.index or gene not in Yall.index:
            return np.nan
        parts = Cbase.index.intersection(Xall.columns).intersection(Yall.columns)
        x = Xall.loc[arm, parts].astype(float).fillna(0).to_numpy()
        y = Yall.loc[gene, parts].astype(float).to_numpy()
        C = Cbase.loc[parts, conf].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
        ok = np.isfinite(y)
        xr = x[ok] - LinearRegression().fit(C[ok], x[ok]).predict(C[ok])
        yr = y[ok] - LinearRegression().fit(C[ok], y[ok]).predict(C[ok])
        return float(spearmanr(xr, yr).correlation)

    rng = np.random.default_rng(0)
    sub = m.sample(min(sample, len(m)), random_state=0) if len(m) > sample else m
    sub = sub.copy(); sub["rho"] = [_coup(r.miRNA, r.gene) for r in sub.itertuples()]
    sub = sub.dropna(subset=["rho"])
    rep = sub.loc[sub["direction"] == "repressive", "rho"]
    act = sub.loc[sub["direction"] == "activating", "rho"]
    print(f"\nTCGA coupling by TarBase direction (n={len(sub)}):")
    print(f"  repressive edges: median rho {rep.median():+.3f} ({(rep<0).mean()*100:.0f}% negative), n={len(rep)}")
    print(f"  activating edges: median rho {act.median():+.3f} ({(act<0).mean()*100:.0f}% negative), n={len(act)}")
    if len(rep) > 5 and len(act) > 5:
        p = mannwhitneyu(rep, act, alternative="less").pvalue
        print(f"  MWU repressive < activating (rho more negative): p={p:.3f}  "
              f"→ {'direction is real (activating edges less repressive)' if p < 0.05 else 'not separable'}")


if __name__ == "__main__":
    he_direction_audit()
    coupling_by_direction()
