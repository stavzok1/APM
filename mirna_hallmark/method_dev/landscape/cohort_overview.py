"""Cohort overview — "what do we deal with": the global shape of miRNA expression across the TCGA-BRCA
cohort, and the gene-centric regulation statistics of the HE edge set. One figure, two rows.

Band 1 — miRNA expression (all 2,236 measured arms):
  A  rank-abundance — cohort-median RPM vs rank (log-log), all measured vs HE-edge arms.
  B  tier composition (robust/conditional/silent) for all arms vs HE arms vs HE EDGES (silent arms low-degree).
  C  per-arm cohort-median abundance histogram, with the RPM-10 detectability floor.
Band 2 — gene regulation (HE edges):
  D  # regulators (miRNAs) per gene — in-degree.
  E  # EXPRESSED regulators per gene (silent removed), with the top convergence-hub genes labelled.
  F  # distinct seeds per gene — seed-diversity convergence.
  G  total functional 7mer+ sites per gene — site load (Σ over regulators).
Band 3 — across-sample spread:
  H  boxplot + per-tumor dots for the top 50 most-abundant arms (RPM, log).

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.landscape.cohort_overview`
Out: `figures/cohort_overview.png`
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
from mirna_hallmark import config as C
from mirna_hallmark.arm_expression import arm_expression_tiers, filter_expressed_edges
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges
from mirna_hallmark.method_dev.landscape.he_edge_arm_landscape import _seed_map

FIG = Path(__file__).parent.parent / "figures"
_SEEDSITES = Path(__file__).parent.parent / "site_ladder" / "utr_seed_sites.tsv.gz"
_LADDER = Path(__file__).parent.parent / "site_ladder" / "utr_site_ladder_genomic.tsv.gz"
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))
_TIER_C = {"robust": "#1b7837", "conditional": "#f1a340", "silent": "#bdbdbd"}


def build() -> None:
    # ---- miRNA expression ---------------------------------------------------
    M = D.load_mirna_arms().astype(float)
    lin = (2.0 ** M - 1.0).clip(lower=0.0)
    med = lin.median(axis=1).sort_values(ascending=False)
    t = arm_expression_tiers()
    thr = float(np.log2(C.ARM_EXPRESSED_MIN_RPM + 1.0))

    # ---- gene regulation (HE edges) ----------------------------------------
    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    seed = _seed_map()
    e["seed"] = e["miRNA"].map(lambda a: seed.get(a, seed.get(_strip(a))))
    indeg = e.groupby("gene")["miRNA"].nunique()
    seeddiv = e.groupby("gene")["seed"].nunique()
    indeg_e = filter_expressed_edges(e).groupby("gene")["miRNA"].nunique().reindex(indeg.index).fillna(0).astype(int)
    sites_pg = (pd.read_csv(_SEEDSITES, sep="\t").groupby("gene")["n_7mer_plus"].sum()      # total functional 7mer+ sites/gene
                .reindex(indeg.index).fillna(0).astype(int))
    he_arms = sorted(e["miRNA"].unique())
    he_tier = pd.Series([t["tier"].get(a, t["tier"].get(_strip(a), "silent")) for a in he_arms]).value_counts()
    e["tier"] = e["miRNA"].map(lambda a: t["tier"].get(a, t["tier"].get(_strip(a), "silent")))
    edge_tier = e["tier"].value_counts()                                          # edges carried by each tier
    he_set, he_strip = set(he_arms), {_strip(a) for a in he_arms}
    med_he = med[med.index.map(lambda a: a in he_set or _strip(a) in he_strip)]   # arms we have edges for (in matrix)

    fig = plt.figure(figsize=(18, 19))
    gs_top = fig.add_gridspec(1, 3, top=0.915, bottom=0.79, left=0.05, right=0.99, wspace=0.26)
    gs_mid = fig.add_gridspec(1, 4, top=0.745, bottom=0.625, left=0.05, right=0.99, wspace=0.34)
    gs_meta = fig.add_gridspec(1, 1, top=0.565, bottom=0.49, left=0.05, right=0.99)
    gs_h = fig.add_gridspec(1, 1, top=0.43, bottom=0.04, left=0.05, right=0.99)

    # A — rank-abundance (log-log) + cumulative pool share
    axA = fig.add_subplot(gs_top[0, 0])
    rank = np.arange(1, len(med) + 1)
    axA.fill_between(rank, med.values.clip(min=1e-2), 1e-2, color="#b8c9e0", alpha=0.85, step="mid",
                     label=f"all measured ({len(med):,})")
    mh = np.sort(med_he.values)[::-1].clip(min=1e-2)
    axA.plot(np.arange(1, len(mh) + 1), mh, color="#b2182b", lw=1.7, label=f"HE-edge arms ({len(mh)})")
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.axhline(C.ARM_EXPRESSED_MIN_RPM, color="crimson", ls="--", lw=1)
    tot = med.sum(); cs = med.cumsum() / tot
    for k in (10, 50):
        axA.axvline(k, color="0.5", lw=0.7, ls=":")
        axA.text(k, axA.get_ylim()[1], f" top {k}\n ={cs.iloc[k-1]:.0%} pool", fontsize=7, va="top", color="0.3")
    axA.set_xlabel("miRNA arm rank (by cohort-median RPM)", fontsize=9)
    axA.set_ylabel("cohort-median RPM", fontsize=9); axA.legend(fontsize=7, loc="lower left")
    axA.set_title(f"A · miRNA rank-abundance — all {len(med):,} measured vs the {len(mh)} HE-edge arms\n"
                  f"top 10 arms = 74% of pool; the arms we model span the whole range", fontsize=10, loc="left")

    # B — tier COMPOSITION (100%-stacked): arms vs edges — silent arms are low-degree
    axB = fig.add_subplot(gs_top[0, 1])
    groups = {f"All arms\n({len(t):,})": t["tier"].value_counts(),
              f"HE arms\n({len(he_arms)})": he_tier,
              f"HE EDGES\n({len(e):,})": edge_tier}
    order = ["robust", "conditional", "silent"]
    for i, (lab, vc) in enumerate(groups.items()):
        total = vc.sum(); bottom = 0
        for tier in order:
            n = int(vc.get(tier, 0)); frac = n / total
            axB.bar(i, frac, bottom=bottom, color=_TIER_C[tier], width=0.62)
            if frac > 0.05:
                axB.text(i, bottom + frac / 2, f"{frac:.0%}\n({n:,})", ha="center", va="center", fontsize=7.5,
                         color="white" if tier != "conditional" else "0.2", fontweight="bold")
            bottom += frac
    axB.set_xticks(range(len(groups))); axB.set_xticklabels(groups.keys(), fontsize=8)
    axB.set_ylabel("composition"); axB.set_ylim(0, 1)
    axB.set_title("B · tier composition — arms vs edges\nsilent = 39% of HE arms but only 11% of edges (silent arms are low-degree)",
                  fontsize=10, loc="left")
    axB.legend([plt.Rectangle((0, 0), 1, 1, color=_TIER_C[k]) for k in order], order, fontsize=7, loc="lower center", ncol=3)

    # C — per-arm median abundance histogram
    axC = fig.add_subplot(gs_top[0, 2])
    axC.hist(t["median_log2rpm"], bins=60, color="#7fbf7b", edgecolor="0.4", linewidth=0.3)
    axC.axvline(thr, color="crimson", ls="--", lw=1, label=f"floor (RPM {C.ARM_EXPRESSED_MIN_RPM:g})")
    axC.set_xlabel("cohort-median expression  log₂(RPM+1)", fontsize=9)
    axC.set_ylabel("# arms", fontsize=9); axC.legend(fontsize=7)
    axC.set_title("C · abundance distribution\nmost arms sit near zero (silent); a thin expressed tail",
                  fontsize=10, loc="left")

    # D/E/F — gene-centric distributions
    def _hist(ax, s, color, xlab, title, top_label=False):
        cap = int(np.percentile(s, 99)); n_last = int((s >= cap).sum())
        ax.hist(s.clip(upper=cap), bins=range(0, cap + 2), color=color, edgecolor="0.4", linewidth=0.3, align="left")
        ax.axvline(s.median(), color="navy", ls="--", lw=1)
        ax.text(s.median(), ax.get_ylim()[1] * 0.92, f" median {s.median():g}", fontsize=7.5, color="navy")
        ax.set_yscale("log"); ax.set_ylabel("# genes (log)", fontsize=9)
        ax.set_xlabel(xlab + (f"   (last bar pools ALL ≥{cap}: {n_last} genes)" if (s > cap).any() else ""), fontsize=9)
        ax.set_title(title, fontsize=10, loc="left")
        if top_label:
            top = s.sort_values(ascending=False).head(6)
            txt = "  ·  ".join(f"{g} {int(v)}" for g, v in top.items())
            ax.text(0.98, 0.78, "most-regulated:\n" + txt.replace("  ·  ", "\n"), transform=ax.transAxes,
                    ha="right", va="top", fontsize=6.8, color="0.25",
                    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    _hist(fig.add_subplot(gs_mid[0, 0]), indeg, "#9ecae1",
          "# regulators (miRNAs) per gene", f"D · regulators per gene (in-degree)\n{e['gene'].nunique():,} genes; mean {indeg.mean():.1f}, max {indeg.max()}")
    _hist(fig.add_subplot(gs_mid[0, 1]), indeg_e, "#fdae6b",
          "# EXPRESSED regulators per gene", f"E · expressed regulators per gene\n(silent arms removed); mean {indeg_e.mean():.1f}", top_label=True)
    _hist(fig.add_subplot(gs_mid[0, 2]), seeddiv, "#bcbddc",
          "# distinct seeds per gene", f"F · seed DIVERSITY per gene (# distinct seeds)\namong its HE-edge regulators, all tiers; max {seeddiv.max()}")
    _hist(fig.add_subplot(gs_mid[0, 3]), sites_pg, "#a1d99b",
          "# 3′UTR sites per gene (Σ over regulators)", f"G · site LOAD per gene (total functional 7mer+ sites)\nHE edges; counts sites not seeds (cf. F); mean {sites_pg.mean():.1f}, max {sites_pg.max()}")

    # H — across-sample distribution of the top 50 most-abundant arms (boxplot + a dot per tumor)
    axH = fig.add_subplot(gs_h[0, 0])
    top50 = med.head(50).index.tolist()
    data = [np.clip(lin.loc[a].to_numpy(), 0.1, None) for a in top50]
    rng = np.random.default_rng(0)
    for i, y in enumerate(data):
        axH.scatter(i + rng.uniform(-0.28, 0.28, len(y)), y, s=1.2, alpha=0.05, color="#08306b",
                    edgecolors="none", zorder=1)
    axH.boxplot(data, positions=range(len(top50)), widths=0.62, showfliers=False, patch_artist=True, zorder=3,
                medianprops=dict(color="crimson", lw=1.3), boxprops=dict(facecolor="none", edgecolor="0.15", lw=0.8),
                whiskerprops=dict(color="0.2", lw=0.7), capprops=dict(color="0.2", lw=0.7))
    axH.axhline(C.ARM_EXPRESSED_MIN_RPM, color="crimson", ls="--", lw=0.8, alpha=0.6)
    axH.set_yscale("log"); axH.set_xlim(-0.7, len(top50) - 0.3)
    axH.set_xticks(range(len(top50))); axH.set_xticklabels([a.replace("hsa-", "") for a in top50], rotation=90, fontsize=6.5)
    axH.set_ylabel("RPM  (log; one dot = one tumor)", fontsize=9)
    axH.set_title(f"H · across-sample expression of the top 50 arms — boxplot (IQR + red median) + a dot per tumor "
                  f"(n={lin.shape[1]:,}); even the dominant arms swing ~10–100× across tumors", fontsize=10, loc="left")

    # ---- I: metagene seed-site position map across all genes' 3'UTRs -------- #
    axM = fig.add_subplot(gs_meta[0, 0])
    lad = pd.read_csv(_LADDER, sep="\t")
    lad["seed"] = lad["miRNA"].map(lambda a: seed.get(a, seed.get(_strip(a))))
    uds = lad.dropna(subset=["seed"]).drop_duplicates(["gene", "utr_pos", "seed"]).copy()   # family-collapsed seed-site instances
    uds["rel"] = (uds["utr_pos"] / uds["utr_len"]).clip(0, 1)
    nb = 40; uds["b"] = (uds["rel"] * nb).clip(0, nb - 1).astype(int)
    ngm = uds["gene"].nunique()
    prof = uds.groupby("b").size().reindex(range(nb), fill_value=0) / ngm                     # mean instances per gene per bin
    xc = (np.arange(nb) + 0.5) * 100 / nb
    axM.fill_between(xc, prof.values, color="#6a51a3", alpha=0.5, step="mid")
    axM.plot(xc, prof.values, color="#4a1486", lw=1.3)
    axM.set_xlim(0, 100); axM.set_ylim(0, prof.max() * 1.15)
    axM.set_xlabel("relative position along 3′UTR (%; stop codon 0 → poly-A 100)", fontsize=9)
    axM.set_ylabel("family-collapsed\nseed sites / gene", fontsize=8)
    axM.set_title(f"I · seed-site position map across all {ngm:,} sited genes' 3′UTRs — family-collapsed instances (same-seed arms once), "
                  f"UTR-length-normalised; sites concentrate near the stop codon (5′ end) and dip mid-UTR", fontsize=10, loc="left")

    fig.suptitle("Cohort overview — miRNA expression (A–C) · gene regulation (D–G) · seed-site 3′UTR position map (I) · top-arm spread (H), TCGA-BRCA",
                 fontsize=14, fontweight="bold", y=0.965)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "cohort_overview.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[overview] {len(med):,} arms ({int(t['expressed'].sum())} expressed) · {e['gene'].nunique():,} genes; "
          f"in-degree med {indeg.median():.0f}/max {indeg.max()}; wrote figures/cohort_overview.png")


if __name__ == "__main__":
    build()
