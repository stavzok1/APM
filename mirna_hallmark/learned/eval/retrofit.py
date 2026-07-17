"""Orphan RETROFIT, genome-wide + gridded over discovery sets. Fold each gene's discovered orphans into the
support, re-fit the dense ridge over [HE ∪ discovered], and measure: (a) do the inserted orphans EARN |z|>2
in the joint fit (per-set validation rate = a precision metric), (b) how much the HE regulators' attribution
SHIFTS (omitted-variable bias), (c) how many HE families FLIP |z|>2 status. Baseline HE fit computed ONCE
per gene and reused across sets.

    python -m mirna_hallmark.learned.eval.retrofit          # genome-wide, learned vs uniform∪densez
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import spike_slab as SS


def _resid(v, Z):
    return v - LinearRegression().fit(Z, v).predict(Z)


def _ridge_z(Y, Xf, C):
    Cm = C.to_numpy(float)
    yr = -_resid(Y.to_numpy(float), Cm)
    Xr = _resid(Xf.to_numpy(float), Cm)
    sd = Xr.std(0)
    Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    m, s, _ = SS._gibbs_posterior(Xz, yr, np.ones(Xz.shape[1]), n_iter=1000, burn=300, seed=0)
    return pd.Series(np.abs(m / (s + 1e-9)), index=Xf.columns)


def _load_set(fname):
    d = pd.read_csv(f"mirna_hallmark/output/learned/{fname}.tsv", sep="\t")
    d = d[d["robust"]] if "robust" in d else d
    return d.groupby("gene")["arm"].apply(list).to_dict()


def retrofit_grid(*, out: str = "mirna_hallmark/output/learned/retrofit_grid.tsv",
                  limit: int | None = None, progress: int = 50) -> pd.DataFrame:
    from pathlib import Path
    learned = _load_set("discovery_bayes_learned")
    uni = _load_set("discovery_bayes_uniform")
    dz = _load_set("discovery_bayes_densez")
    unidz = {g: sorted(set(uni.get(g, [])) | set(dz.get(g, []))) for g in set(uni) | set(dz)}
    SETS = {"learned": learned, "uniform∪densez": unidz}
    genes = sorted(set().union(*[set(s) for s in SETS.values()]))
    if limit:
        genes = genes[:limit]
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[retrofit] {i}/{len(genes)}", flush=True)
        try:
            Y, Xo, C, w = LD.assemble_gene(g, w_prior_source="ledger", orphans=True, deconv=True)
            _, Xhe, _, _ = LD.assemble_gene(g, w_prior_source="ledger", deconv=True)
            he_arms = [a for a in Xhe.columns if a in Xo.columns]
            if len(he_arms) < 2:
                continue
            Xbf, _, _ = FAM.collapse_by_family(Xo[he_arms], w.reindex(he_arms), FAM.family_of(he_arms))
            zb = _ridge_z(Y, Xbf, C)                                    # baseline HE fit (once)
            for sname, sset in SETS.items():
                darms = [a for a in sset.get(g, []) if a in Xo.columns and a not in he_arms]
                if not darms:
                    continue
                keep = he_arms + darms
                Xak, _, _ = FAM.collapse_by_family(Xo[keep], w.reindex(keep), FAM.family_of(keep))
                za = _ridge_z(Y, Xak, C)
                he_fams = [f for f in zb.index if f in za.index]
                disc_fams = [f for f in za.index if f not in zb.index]
                shifts = [za[f] - zb[f] for f in he_fams]
                flips = sum(1 for f in he_fams if (zb[f] > 2) != (za[f] > 2))
                val = [za[f] > 2 for f in disc_fams]
                rows.append({"gene": g, "set": sname, "n_disc": len(disc_fams), "n_he": len(he_fams),
                             "disc_valid": int(np.sum(val)), "he_shift_mean": float(np.mean(shifts)) if shifts else np.nan,
                             "he_flips": flips})
        except Exception:
            pass
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\n=== ORPHAN RETROFIT (genome-wide), gridded over discovery sets ===")
    for sname in SETS:
        d = df[df["set"] == sname]
        if not len(d):
            continue
        nd, nv = d["n_disc"].sum(), d["disc_valid"].sum()
        nhe, nfl = d["n_he"].sum(), d["he_flips"].sum()
        print(f"  {sname:14s}: genes {len(d)} | inserted orphans {nd} | **validated |z|>2 in joint fit "
              f"{nv}/{nd} ({100*nv/max(nd,1):.0f}%)** | HE-family mean Δ|z| {d['he_shift_mean'].mean():+.3f} | "
              f"HE flips {nfl}/{nhe} ({100*nfl/max(nhe,1):.0f}%)")
    return df


if __name__ == "__main__":
    retrofit_grid()
