"""Which inclusion-prior knob owns COUPLING? A p0 × spec_gain grid for the spike-and-slab (π≠1), same
genes/folds/C/slab, vs the dense π≡1 reference. Tests the claim:

  * p0 (inclusion LEVEL) MOVES coupling — as p0→1 the spike-and-slab → dense (the coupling winner).
  * spec_gain (specificity ORDERING) is ~NULL to coupling — it re-nominates the specialist (identity), not the
    magnitude the abundant arms carry.

Expect: rho_ss improves monotonically with p0 (→ dense) and is FLAT across spec_gain at each p0.

Per config: 5-fold OOF ρ per gene; mean ρ + BY-FDR-significant fraction across genes (spike_slab.gate_fdr recipe).

CLI: .venv/bin/python3 -m mirna_hallmark.learned.genome_p0_spec_sweep --limit 80
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

P0S = [0.2, 0.4, 0.6, 0.8]
SPECS = [0.0, 2.0]
FLOOR, SCALE = 3.46, 1.5


def _p_neg(rho, n):
    dfree = int(n) - 8
    if not (rho == rho) or dfree <= 1:
        return np.nan
    if rho >= 0:
        return 1.0
    t = rho * np.sqrt(dfree / max(1.0 - rho * rho, 1e-9))
    return float(_st.t.cdf(t, dfree))


def run(*, limit=80, n_iter=700, out="mirna_hallmark/output/learned/genome_p0_spec_sweep.tsv"):
    from mirna_hallmark.learned.evidence import ledger as LG
    genes = sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))[:limit]
    configs = [("dense", None, None)] + [(f"p0={p0}_sg={sg}", p0, sg) for p0 in P0S for sg in SPECS]
    rows = []
    for i, g in enumerate(genes):
        if i % 20 == 0:
            print(f"[p0_spec_sweep] {i}/{len(genes)}", flush=True)
        try:
            p = RB._prep(g)
            n = len(p["Y"])
            if p["X"].shape[1] < 2:
                continue
            slab = RB._slab(p, locus="base", gain=0.0)
            seq = SQ.seq_spec(g, p["cols"], universe="HE").to_numpy(float)
            med = np.median(p["X"].to_numpy(float), axis=0)
            rec = {"gene": g, "n": n, "n_pred": p["X"].shape[1]}
            for name, p0, sg in configs:
                if name == "dense":
                    rho = RB._oof(p, dense=True, pi_full=None, slab=slab, n_iter=n_iter)[0]
                else:
                    pi = SB._pi_channel(p, seq, med, gain=sg, gated=True, floor=FLOOR, scale=SCALE, p0=p0)
                    rho = RB._oof(p, dense=False, pi_full=pi, slab=slab, n_iter=n_iter)[0]
                rec[f"rho_{name}"] = rho
            rows.append(rec)
        except Exception:
            pass
    df = pd.DataFrame(rows)
    for name, *_ in configs:
        col = f"rho_{name}"
        df[f"p_{name}"] = [_p_neg(r, n) for r, n in zip(df[col], df["n"])]
        df[f"q_{name}"] = benjamini_yekutieli(df[f"p_{name}"].values)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\n=== p0 × spec_gain COUPLING+FDR grid ({len(df)} genes) → {out} ===")
    print(f"{'config':16s} {'mean_rho':>9s} {'FDR-sig%':>9s}")
    for name, *_ in configs:
        sig = int((df[f"q_{name}"] < 0.05).sum())
        print(f"{name:16s} {df[f'rho_{name}'].mean():+9.4f} {100*sig/max(len(df),1):8.0f}%")
    print("\nMEAN ρ table (rows=p0, cols=spec_gain) — want: DOWN each column (p0 raises coupling), FLAT across a row:")
    tab = pd.DataFrame({sg: [df[f"rho_p0={p0}_sg={sg}"].mean() for p0 in P0S] for sg in SPECS}, index=[f"p0={p0}" for p0 in P0S])
    print(tab.round(4).to_string())
    print(f"dense reference: {df['rho_dense'].mean():+.4f}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=80)
    ap.add_argument("--n-iter", type=int, default=700)
    a = ap.parse_args()
    run(limit=a.limit, n_iter=a.n_iter)
