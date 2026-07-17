"""Genome-wide coupling + FDR for the NON-INERT π model — the missing test. Does making π ≠ 1 carry the
specificity prior (seq-spec) change the coupling gate vs (a) the dense π≡1 winner and (b) the evidence-π inclusion
baseline? Same genes / folds / C / slab for all three, so the ONLY difference is the inclusion π.

  dense    : π ≡ 1                          (the §2d coupling winner; prior-inert)
  evid-π   : inclusion, evidence-only π      (the MH-83 inclusion baseline)
  seqspec-π: inclusion, evidence + gated seq-spec π   (the non-inert specificity model)

Per gene: 5-fold OOF ρ for each; then one-sided partial-t p(ρ<0), df≈n−8, BY-FDR across genes (same recipe as
spike_slab.gate_fdr). Reports mean ρ and FDR-significant fraction per method.

CLI: .venv/bin/python3 -m mirna_hallmark.learned.genome_seq_fdr --limit 120 --gain 1.2
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as _st

from mirna_hallmark.learned import rarity_bench as RB
from mirna_hallmark.learned import seq_specificity as SQ
from mirna_hallmark.learned import seq_spec_bench as SB
from mirna_hallmark.coupling_inference import benjamini_yekutieli


def _p_neg(rho, n):
    dfree = int(n) - 8
    if not (rho == rho) or dfree <= 1 or rho >= 0:
        return np.nan if rho != rho else 1.0
    t = rho * np.sqrt(dfree / max(1.0 - rho * rho, 1e-9))
    return float(_st.t.cdf(t, dfree))


def run(*, limit=120, gain=1.2, n_iter=800, out="mirna_hallmark/output/learned/genome_seq_fdr.tsv"):
    from mirna_hallmark.learned.evidence import ledger as LG
    genes = sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))[:limit]
    rows = []
    for i, g in enumerate(genes):
        if i % 25 == 0:
            print(f"[genome_seq_fdr] {i}/{len(genes)}", flush=True)
        try:
            p = RB._prep(g)
            n = len(p["Y"])
            if p["X"].shape[1] < 2:
                continue
            slab = RB._slab(p, locus="base", gain=0.0)
            rho_dense = RB._oof(p, dense=True, pi_full=None, slab=slab, n_iter=n_iter)[0]
            pi_ev = RB._pi(p, rarity_on=False, gain=0.0)                      # evidence-only π
            rho_ev = RB._oof(p, dense=False, pi_full=pi_ev, slab=slab, n_iter=n_iter)[0]
            seq = SQ.seq_spec(g, p["cols"], universe="HE").to_numpy(float)
            med = np.median(p["X"].to_numpy(float), axis=0)
            floor = 3.46
            pi_sq = SB._pi_channel(p, seq, med, gain=gain, gated=True, floor=floor, scale=1.5)
            rho_sq = RB._oof(p, dense=False, pi_full=pi_sq, slab=slab, n_iter=n_iter)[0]
            rows.append({"gene": g, "n": n, "n_pred": p["X"].shape[1],
                         "rho_dense": rho_dense, "rho_evpi": rho_ev, "rho_seqpi": rho_sq})
        except Exception:
            pass
    df = pd.DataFrame(rows)
    for tag in ["dense", "evpi", "seqpi"]:
        df[f"p_{tag}"] = [_p_neg(r, n) for r, n in zip(df[f"rho_{tag}"], df["n"])]
        df[f"q_{tag}"] = benjamini_yekutieli(df[f"p_{tag}"].values)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\n=== NON-INERT π: coupling + FDR ({len(df)} genes) → {out} ===")
    for tag, lbl in [("dense", "dense π≡1  "), ("evpi", "evidence-π "), ("seqpi", "seqspec-π  ")]:
        sig = int((df[f"q_{tag}"] < 0.05).sum())
        print(f"  {lbl}: mean ρ {df[f'rho_{tag}'].mean():+.4f} | FDR-sig(q_BY<0.05 & ρ<0) {sig}/{len(df)} "
              f"({100*sig/max(len(df),1):.0f}%)")
    print(f"  Δρ(seqspec−evidence) = {(df['rho_seqpi']-df['rho_evpi']).mean():+.4f} | "
          f"Δρ(seqspec−dense) = {(df['rho_seqpi']-df['rho_dense']).mean():+.4f}")
    print(f"  seqspec beats evidence on {int((df['rho_seqpi']<df['rho_evpi']-1e-6).sum())}/{len(df)} | "
          f"seqspec beats dense on {int((df['rho_seqpi']<df['rho_dense']-1e-6).sum())}/{len(df)}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=120)
    ap.add_argument("--gain", type=float, default=1.2)
    ap.add_argument("--n-iter", type=int, default=800)
    a = ap.parse_args()
    run(limit=a.limit, gain=a.gain, n_iter=a.n_iter)
