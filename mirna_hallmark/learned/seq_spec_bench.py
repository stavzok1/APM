"""Test the VALIDATED seq-specificity as an expression-gated inclusion-π channel — the retry of the goal that
seed-rarity failed. Does promoting high-affinity-CONCENTRATION edges in π lift the RIGHT specialists (expressed +
coupling), unlike rarity (which lifted unexpressed/non-coupling orphans like miR-105/718)?

π channel:  logit(π_{g,m}) = logit(p0) + κ·z_w + gain · gate_m · z(seq_spec(m,g))
  gate_m = expression gate ∈(0,1] = expit((median_log2_expr_m − floor)/scale) — no promotion of undetectable arms.
  seq_spec = affinity concentration on g (HE universe), z-scored over the gene's covered candidates.

Audit (the thing rarity failed): for each variant, the edges whose PIP crosses 0.5 — are they EXPRESSED and do they
COUPLE at the edge level (edge_coupling.tsv partial-Spearman)? Plus OOF SS coupling (must stay neutral) on soups
{PTEN,CCND1,VEGFA} vs clean {ESR1,CDH1}.

CLI: .venv/bin/python3 -m mirna_hallmark.learned.seq_spec_bench --gains 0.8 1.6
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.special import expit

from mirna_hallmark.learned import rarity_bench as RB
from mirna_hallmark.learned import seq_specificity as SQ
from mirna_hallmark.learned import spike_slab as SS

SOUP = ["PTEN", "CCND1", "VEGFA"]
CLEAN = ["ESR1", "CDH1"]
_EC = None


def _edge_coupling():
    global _EC
    if _EC is None:
        _EC = pd.read_csv("mirna_hallmark/output/learned/edge_coupling.tsv", sep="\t")
        _EC["a"] = _EC["arm"].str.replace("hsa-", "", regex=False)
    return _EC


def _fam_edge_rho(gene, fam):
    ec = _edge_coupling()
    best = np.nan
    for m in str(fam).split("/"):
        k = m if m.startswith(("miR-", "let-")) else "miR-" + m
        e = ec[(ec["gene"] == gene) & (ec["a"] == k)]
        if len(e):
            r = float(e["rho"].iloc[0])
            best = r if best != best else min(best, r)
    return best


def _expr_gate(med, *, floor: float, scale: float = 1.0):
    return expit((med - floor) / scale)


def _pi_channel(p, seq, med, *, gain, gated, floor, scale, p0=0.3, kappa=1.5):
    """Evidence π + gain·gate·z(seq_spec). gated=False ⇒ gate≡1."""
    w = np.log1p(np.clip(p["w"], 0, None))
    medn = np.median(w); mad = np.median(np.abs(w - medn)) * 1.4826
    zw = (w - medn) / mad if mad > 1e-9 else np.zeros_like(w)
    logit = np.log(p0 / (1 - p0)) + kappa * zw
    z = SQ.zscore(seq)
    gate = _expr_gate(med, floor=floor, scale=scale) if gated else np.ones_like(med)
    return np.clip(expit(logit + gain * gate * z), 0.02, 0.98)


def _full_pip(p, pi):
    """Full-data spike-and-slab PIP per family (for the audit of who gets selected)."""
    Y = p["Y"].to_numpy(float); X = p["X"].to_numpy(float); C = p["C"].to_numpy(float)
    r = -RB._resid(Y, C); Xr = RB._resid(X, C); sd = Xr.std(0)
    Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    pv = np.where(sd > 1e-9, pi, 0.0)
    _, pip = SS._gibbs_ss(Xz, r, pv, n_iter=1000, burn=350, seed=0)
    return pd.Series(pip, index=p["cols"])


def bench(genes=None, *, gains=(0.8, 1.6), floor_q=0.20):
    genes = genes or (SOUP + CLEAN)
    audit_rows = []
    for g in genes:
        p = RB._prep(g)
        seq = SQ.seq_spec(g, p["cols"], universe="HE").to_numpy(float)
        med = np.median(p["X"].to_numpy(float), axis=0)                     # per-family median log2 expr
        floor = float(np.quantile(med, floor_q))                            # gate floor = 20th pct of candidates
        klass = "soup" if g in SOUP else "clean"
        # baseline evidence-only π
        pi_base = _pi_channel(p, seq, med, gain=0.0, gated=False, floor=floor, scale=1.0)
        pip_base = _full_pip(p, pi_base)
        rho_base, *_ = RB._oof(p, dense=False, pi_full=pi_base, slab=RB._slab(p, locus="base", gain=0.0))
        for gain in gains:
            for gated in (True, False):
                pi = _pi_channel(p, seq, med, gain=gain, gated=gated, floor=floor, scale=1.0)
                pip = _full_pip(p, pi)
                rho, *_ = RB._oof(p, dense=False, pi_full=pi, slab=RB._slab(p, locus="base", gain=0.0))
                # PROMOTED edges: PIP crossed 0.5 from base
                promoted = [(c, round(float(pip[c]), 2)) for c in p["cols"]
                            if pip[c] > 0.5 and pip_base.get(c, 0) <= 0.5]
                # audit each promoted edge: seq-spec, expression, edge coupling
                det = []
                for c, pp in promoted:
                    i = list(p["cols"]).index(c)
                    det.append(f"{c.split('/')[0]}(pip{pp},seq{seq[i]:.2f},exp{med[i]:.1f},ρ{_fmt(_fam_edge_rho(g,c))})")
                audit_rows.append({"gene": g, "class": klass, "gain": gain,
                                   "gate": "on" if gated else "off",
                                   "rho_ss": rho, "d_rho": round(rho - rho_base, 3),
                                   "n_promoted": len(promoted), "promoted": "; ".join(det) or "-"})
        print(f"[seq_spec_bench] {g} done (floor={floor:.1f})", flush=True)
    df = pd.DataFrame(audit_rows)
    with pd.option_context("display.width", 240, "display.max_colwidth", 120, "display.max_rows", 100):
        print("\n=== promoted edges + audit (want: expressed exp≳floor, coupling ρ<0) ===")
        print(df.to_string(index=False))
    print("\nΔρ_ss vs base by class×gate (want ≈0 or negative = coupling-safe):")
    print(df.groupby(["class", "gate", "gain"])["d_rho"].mean().round(3).to_string())
    df.to_csv("mirna_hallmark/output/learned/seq_spec_bench.tsv", sep="\t", index=False)
    return df


def _fmt(v):
    return f"{v:+.2f}" if v == v else "NA"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--gains", nargs="*", type=float, default=[0.8, 1.6])
    a = ap.parse_args()
    bench(a.genes, gains=a.gains)
