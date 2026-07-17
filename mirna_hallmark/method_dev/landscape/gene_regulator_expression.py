"""How does expression distribute *across a gene's regulators*? For every HE-edge gene we take its regulator
miRNAs and their cohort-median abundance, and ask whether the gene "sees" many comparable regulators or is
dominated by one. Companion to `cohort_overview.py` (which is per-arm / per-gene-count); this one is about the
**abundance shape within a gene's regulator set**.

Panels:
  A  within-gene rank-abundance — each gene's regulators ranked & normalized to its top regulator (=1.0),
     median + IQR across genes. The #2 regulator is typically a small fraction of #1.
  B  dominance — distribution of the top regulator's share of the gene's total regulator abundance.
  C  effective vs nominal regulators — Simpson effective-N (1/Σpᵢ²) vs the raw count; effective stays ~1–2
     however many regulators a gene nominally has (most are abundance-negligible).

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.landscape.gene_regulator_expression`
Out: `figures/gene_regulator_expression.png`
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

import mirna_hallmark.data_loaders as D
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges

FIG = Path(__file__).parent.parent / "figures"
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))


def build() -> None:
    M = D.load_mirna_arms().astype(float)
    med = (2.0 ** M - 1.0).clip(lower=0.0).median(axis=1)
    idx = set(med.index)
    rmed = lambda a: float(med[a]) if a in idx else (float(med[_strip(a)]) if _strip(a) in idx else 0.0)

    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    e["rpm"] = e["miRNA"].map(rmed)

    profiles, top1, effN, nomN = [], [], [], []
    for _, sub in e.groupby("gene"):
        v = np.sort(sub["rpm"].to_numpy())[::-1]
        nomN_g = len(v)
        if nomN_g < 2 or v[0] <= 0:
            continue
        p = v / v.sum()
        profiles.append(v / v[0])
        top1.append(p.max()); effN.append(1.0 / np.sum(p ** 2)); nomN.append(nomN_g)
    top1 = np.array(top1); effN = np.array(effN); nomN = np.array(nomN)

    fig = plt.figure(figsize=(16, 5.2))
    gs = GridSpec(1, 3, wspace=0.27, top=0.82, bottom=0.16, left=0.06, right=0.98)

    # A — within-gene normalized rank-abundance (median + IQR)
    axA = fig.add_subplot(gs[0, 0])
    maxr = 12
    rr = np.arange(1, maxr + 1)
    stk = [[p[r - 1] for p in profiles if len(p) >= r] for r in rr]
    med_p = [np.median(c) for c in stk]
    q1 = [np.percentile(c, 25) for c in stk]; q3 = [np.percentile(c, 75) for c in stk]
    axA.fill_between(rr, np.clip(q1, 1e-3, None), np.clip(q3, 1e-3, None), color="#9ecae1", alpha=0.45, label="IQR across genes")
    axA.plot(rr, np.clip(med_p, 1e-3, None), "o-", color="#08519c", lw=1.8, ms=5, label="median")
    axA.set_yscale("log"); axA.set_xticks(rr)
    axA.set_xlabel("regulator rank within a gene (1 = most abundant)", fontsize=9)
    axA.set_ylabel("abundance ÷ gene's top regulator", fontsize=9); axA.legend(fontsize=8)
    axA.set_title(f"A · within-gene rank-abundance ({len(profiles)} genes)\nthe #2 regulator is typically "
                  f"~{med_p[1]:.0%} of the #1, #3 ~{med_p[2]:.0%}", fontsize=10, loc="left")

    # B — dominance: top-1 regulator share
    axB = fig.add_subplot(gs[0, 1])
    axB.hist(top1, bins=np.linspace(0, 1, 26), color="#fc8d59", edgecolor="0.4", linewidth=0.3)
    axB.axvline(np.median(top1), color="navy", ls="--", lw=1.2)
    axB.text(np.median(top1), axB.get_ylim()[1] * 0.95, f" median {np.median(top1):.0%}", color="navy", fontsize=8, va="top")
    axB.set_xlabel("top regulator's share of the gene's total regulator abundance", fontsize=9)
    axB.set_ylabel("# genes", fontsize=9)
    axB.set_title(f"B · dominance — one regulator usually owns the gene\n{(top1>0.5).mean():.0%} of genes "
                  f">50%-dominated, {(top1>0.8).mean():.0%} >80%", fontsize=10, loc="left")

    # C — effective vs nominal regulators
    axC = fig.add_subplot(gs[0, 2])
    axC.scatter(nomN + np.random.default_rng(0).uniform(-0.2, 0.2, len(nomN)), effN, s=8, color="#74c476",
                alpha=0.35, edgecolors="none")
    bins = [2, 3, 4, 6, 9, 14, 22, 40, 95]
    cx, cy = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (nomN >= lo) & (nomN < hi)
        if m.sum() >= 5:
            cx.append(np.median(nomN[m])); cy.append(np.median(effN[m]))
    axC.plot(cx, cy, "o-", color="#00441b", lw=2, ms=6, label="median effective-N per bin")
    lim = nomN.max()
    axC.plot([0, lim], [0, lim], "--", color="0.6", lw=1, label="y = x (all equal)")
    axC.set_xscale("log"); axC.set_ylim(0.8, np.ceil(np.percentile(effN, 99)) + 1)
    axC.set_xlabel("nominal # regulators (in-degree)", fontsize=9)
    axC.set_ylabel("effective # regulators  (Simpson 1/Σpᵢ²)", fontsize=9); axC.legend(fontsize=8, loc="upper left")
    axC.set_title(f"C · effective vs nominal regulators\neffective-N stays ~{np.median(effN):.1f} "
                  f"(median) however large the in-degree", fontsize=10, loc="left")

    fig.suptitle("How expression distributes across a gene's regulators — abundance is concentrated in ~1 regulator per gene",
                 fontsize=13, fontweight="bold", y=0.965)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "gene_regulator_expression.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[gene-reg-expr] {len(profiles)} genes | top-1 share median {np.median(top1):.0%} | "
          f"effective-N median {np.median(effN):.1f} vs nominal {np.median(nomN):.0f}; wrote figures/gene_regulator_expression.png")


if __name__ == "__main__":
    build()
