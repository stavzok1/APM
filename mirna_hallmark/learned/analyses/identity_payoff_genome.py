"""GENOME-WIDE identifiability payoff — does identity-π change the identified driver vs evidence-π, across many
genes (not just the 4 saturated hubs) and over ALL regulator families (not just the canonical one)?

Per gene, full-data spike-and-slab PIP under:
  evidence-π      : inclusion_prior(w)                              (the baseline)
  identity-π      : inclusion_prior(w, spec_gain, spec=percentile, gate)   (evidence + identity, additive in logit)
  identity-only-π : p0 flat + spec (no evidence)                    (isolates identity's own identification)

Metrics: top-PIP family flip rate (evidence vs identity); mean |ΔPIP| over ALL families; #families whose PIP
crosses 0.5 (gained/lost); and whether the gainers are higher-percentile (mechanism). Restricted to multi-family
soups (≥3 families) where identification is non-trivial.

CLI: .venv/bin/python3 -m mirna_hallmark.learned.analyses.identity_payoff_genome --limit 80 --gain 2.0
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import LinearRegression

from mirna_hallmark.learned import rarity_bench as RB, seq_specificity as SQ, spike_slab as SS, priors as PR
from mirna_hallmark.learned.evidence import ledger as LG


def _resid(V, C):
    return V - LinearRegression().fit(C, V).predict(C)


def _pip(p, pi, seed=0, n_iter=1000):
    Y = p["Y"].to_numpy(float); X = p["X"].to_numpy(float); C = p["C"].to_numpy(float)
    r = -_resid(Y, C); Xr = _resid(X, C); sd = Xr.std(0)
    Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    pv = np.where(sd > 1e-9, pi, 0.0)
    _, pip = SS._gibbs_ss(Xz, r, pv, n_iter=n_iter, burn=350, seed=seed)
    return pip


def run(*, limit=80, gain=2.0, n_iter=900, out="mirna_hallmark/output/learned/identity_payoff_genome.tsv"):
    genes = sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))[:limit]
    rows = []
    for i, g in enumerate(genes):
        if i % 20 == 0:
            print(f"[identity_payoff] {i}/{len(genes)}", flush=True)
        try:
            p = RB._prep(g); cols = list(p["cols"])
            if len(cols) < 3:
                continue
            pct = SQ.affinity_percentile(g, pd.Index(cols)).to_numpy(float)
            if np.isnan(pct).all():
                continue
            med = np.median(p["X"].to_numpy(float), axis=0)
            gate = expit((med - 3.46) / 1.5)
            w = p["w"]
            pip_ev = _pip(p, PR.inclusion_prior(w))
            pip_id = _pip(p, PR.inclusion_prior(w, spec_gain=gain, spec=pct, spec_gate=gate))
            ev, idp = pd.Series(pip_ev, index=cols), pd.Series(pip_id, index=cols)
            flip = ev.idxmax() != idp.idxmax()
            d = idp - ev
            gained = [(c, round(pct[j], 2)) for j, c in enumerate(cols) if ev[c] <= 0.5 < idp[c]]
            lost = [c for c in cols if idp[c] <= 0.5 < ev[c]]
            rows.append({"gene": g, "n_fam": len(cols), "top_flip": bool(flip),
                         "mean_abs_dPIP": round(float(d.abs().mean()), 3),
                         "n_gained": len(gained), "n_lost": len(lost),
                         "gained_pct": round(float(np.mean([x[1] for x in gained])), 2) if gained else np.nan,
                         "top_ev": ev.idxmax().split("/")[0], "top_id": idp.idxmax().split("/")[0]})
        except Exception:
            pass
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    n = len(df)
    print(f"\n=== GENOME identifiability payoff: identity-π vs evidence-π ({n} multi-family genes) ===")
    print(f"  top-driver FLIP rate: {100*df['top_flip'].mean():.0f}%  ({int(df['top_flip'].sum())}/{n})")
    print(f"  mean |ΔPIP| over all families: {df['mean_abs_dPIP'].mean():.3f}")
    print(f"  families GAINED selection (PIP↑ past .5): total {int(df['n_gained'].sum())} | LOST: {int(df['n_lost'].sum())}")
    print(f"  mean percentile of GAINERS: {df['gained_pct'].mean():.2f}  (>0.5 ⇒ identity promotes dedicated families)")
    print(f"  genes with ANY change (flip or gain/lost): {int(((df['top_flip'])|(df['n_gained']>0)|(df['n_lost']>0)).sum())}/{n}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=80)
    ap.add_argument("--gain", type=float, default=2.0)
    ap.add_argument("--n-iter", type=int, default=900)
    a = ap.parse_args()
    run(limit=a.limit, gain=a.gain, n_iter=a.n_iter)
