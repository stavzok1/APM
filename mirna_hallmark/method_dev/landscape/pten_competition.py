"""PTEN miRNA-competition figure — layered overview of the competition mechanism for one gene.

Panels:
  A  expression hierarchy of all PTEN regulators (cohort-median RPM + cross-sample IQR, coloured by
     expression tier) — how many miRNAs compete at each abundance, and how few dominate the RISC pool.
  B  seed structure bars: (B1) #3'UTR sites per miRNA, (B2) #miRNAs per seed (redundancy), (B3) site-type
     and seed counts.
  C  PTEN 3'UTR competition map: the EXPRESSED + sited regulators placed on the 6,461-nt 3'UTR at their
     seed-site positions, coloured by seed; same-seed rows = direct (seed) competition, vertical site
     clusters (≤40 nt) = positional competition (shaded). Only expressed arms (those that can actually
     occupy RISC) are shown.

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.landscape.pten_competition [GENE]`
Out: `figures/<GENE>_competition.png` (+ regenerates `figures/hub_<GENE>.png`).
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

import mirna_hallmark.data_loaders as D
from mirna_hallmark.arm_expression import arm_expression_tiers
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges
from mirna_hallmark.method_dev.landscape.he_edge_arm_landscape import _seed_map

OUT = Path(__file__).parent.parent / "figures"                              # shared figures dir at method_dev root
_LADDER = Path(__file__).parent.parent / "site_ladder" / "utr_site_ladder_genomic.tsv.gz"   # ladder data in site_ladder/
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))
_TIER_C = {"robust": "#1b7837", "conditional": "#f1a340", "silent": "#bdbdbd"}
_COMP_WIN = 40        # nt: S1 upper bound — sites 8–40 nt apart = disjoint-proximal (cooperative AGO)
_OVERLAP_WIN = 8      # nt: S2 — seed footprints overlap = mutually exclusive (true competition)


def _data(gene: str):
    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    arms = sorted(e.loc[e["gene"] == gene, "miRNA"].unique())
    t = arm_expression_tiers(); seed = _seed_map()
    lin = (2.0 ** D.load_mirna_arms().astype(float) - 1.0).clip(lower=0.0)   # linear RPM
    idx = set(lin.index)

    def row(a):
        k = a if a in idx else (_strip(a) if _strip(a) in idx else None)
        v = lin.loc[k] if k else None
        return dict(
            arm=a, seed=seed.get(a, seed.get(_strip(a))),
            tier=t["tier"].get(a, t["tier"].get(_strip(a), "silent")),
            med=float(v.median()) if v is not None else 0.0,
            q25=float(v.quantile(0.25)) if v is not None else 0.0,
            q75=float(v.quantile(0.75)) if v is not None else 0.0)

    df = pd.DataFrame([row(a) for a in arms])
    lad = pd.read_csv(_LADDER, sep="\t")
    sites = lad[lad["gene"] == gene][["miRNA", "utr_pos", "type", "utr_len"]].copy()
    df["n_sites"] = df["arm"].map(sites.groupby("miRNA").size()).fillna(0).astype(int)
    return df, sites


def build(gene: str = "PTEN") -> None:
    df, sites = _data(gene)
    utr_len = int(sites["utr_len"].iloc[0]) if len(sites) else 0
    fig = plt.figure(figsize=(15, 22))
    gs = GridSpec(5, 3, height_ratios=[0.9, 0.6, 0.42, 3.0, 0.78], hspace=0.42, wspace=0.28,
                  top=0.97, bottom=0.03, left=0.055, right=0.985)

    # ---- A: expression hierarchy ------------------------------------------- #
    axA = fig.add_subplot(gs[0, :])
    d = df.sort_values("med", ascending=False).reset_index(drop=True)
    x = np.arange(len(d))
    axA.bar(x, d["med"].clip(lower=0.05), color=[_TIER_C[t] for t in d["tier"]], width=0.9, zorder=3)
    axA.errorbar(x, d["med"].clip(lower=0.05), yerr=[d["med"] - d["q25"], d["q75"] - d["med"]],
                 fmt="none", ecolor="0.35", elinewidth=0.5, capsize=0, zorder=4, alpha=0.6)
    axA.axhline(10, color="crimson", lw=1, ls="--", zorder=2, label="expression floor (RPM 10)")
    axA.set_yscale("log"); axA.set_ylabel("cohort-median RPM\n(bar; whisker = IQR across tumors)", fontsize=9)
    axA.set_xlabel(f"{gene} regulator miRNAs, ranked by abundance (n={len(d)})", fontsize=9)
    tot = d["med"].sum(); top50 = int((d["med"].cumsum() / tot <= 0.5).sum()) + 1 if tot else 0
    axA.set_xlim(-1, len(d)); axA.set_title(f"A · {gene} regulators — expression hierarchy "
        f"({(d['tier']!='silent').sum()} expressed / {(d['tier']=='silent').sum()} silent; "
        f"top {top50} miRNAs = 50% of the RISC-competing pool)", fontsize=11, loc="left")
    h = [plt.Rectangle((0, 0), 1, 1, color=_TIER_C[k]) for k in ("robust", "conditional", "silent")]
    axA.legend(h + [plt.Line2D([], [], color="crimson", ls="--")],
               ["robust (median≥10)", "conditional (spikes≥10)", "silent (never≥10)", "RPM 10 floor"],
               fontsize=7, ncol=4, loc="upper right")

    # ---- B: seed-structure bars -------------------------------------------- #
    sit = df[df["n_sites"] > 0]
    axB1 = fig.add_subplot(gs[1, 0])
    vc = sit["n_sites"].value_counts().sort_index()
    axB1.bar(vc.index, vc.values, color="#4575b4"); axB1.set_xlabel("# 3′UTR sites for a miRNA", fontsize=8)
    axB1.set_ylabel("# miRNAs", fontsize=8); axB1.set_title("B1 · sites per miRNA", fontsize=9, loc="left")

    axB2 = fig.add_subplot(gs[1, 1])
    fam = sit.dropna(subset=["seed"]).groupby("seed")["arm"].nunique()
    fc = fam.value_counts().sort_index()
    axB2.bar(fc.index, fc.values, color="#91bfdb"); axB2.set_xlabel("# miRNAs sharing a seed", fontsize=8)
    axB2.set_ylabel("# seeds", fontsize=8); axB2.set_title(f"B2 · miRNAs per seed ({fam.size} seeds)", fontsize=9, loc="left")

    axB3 = fig.add_subplot(gs[1, 2])
    tc = sites["type"].value_counts().reindex(["8mer", "7mer-m8", "7mer-A1"]).fillna(0)
    axB3.bar(tc.index, tc.values, color=["#d73027", "#fc8d59", "#fee090"])
    axB3.set_ylabel("# sites", fontsize=8); axB3.set_title(f"B3 · site types ({len(sites)} sites)", fontsize=9, loc="left")
    axB3.tick_params(axis="x", labelsize=8)

    # ---- C0: distinct-seed density along the 3'UTR (aligned with C below) --- #
    axSd = fig.add_subplot(gs[2, :])
    _sm = _seed_map()
    sd = sites.assign(seed=sites["miRNA"].map(lambda a: _sm.get(a, _sm.get(_strip(a)))))
    uds = sd.dropna(subset=["seed"]).drop_duplicates(["utr_pos", "seed"])       # one row per (position, seed): family arms collapsed, but each position of a seed kept
    _W, _st = 150, 20
    _xs = np.arange(0, utr_len + 1, _st)
    _dens = np.array([int(uds["utr_pos"].between(x - _W / 2, x + _W / 2).sum()) for x in _xs])   # COUNT seed-site instances in window
    axSd.fill_between(_xs, _dens, color="#6a51a3", alpha=0.55, step="mid")
    axSd.set_xlim(-utr_len * 0.14, utr_len * 1.02); axSd.set_ylim(0, (_dens.max() + 1) if len(_dens) else 1)
    axSd.set_xticks([]); axSd.set_ylabel(f"seed sites\nper {_W} nt", fontsize=8)
    axSd.set_title(f"C₀ · seed-site density along the {gene} 3′UTR — family-collapsed **instances** (each seed counted at every position; same-seed arms once); "
                   f"sliding {_W} nt window ({len(uds)} seed-site instances / {uds['seed'].nunique()} distinct seeds; peaks = binding-site hotspots)", fontsize=10, loc="left")

    # ---- C: 3'UTR competition map (ALL expressed + sited regulators) ------- #
    axC = fig.add_subplot(gs[3, :])
    exp = df[(df["tier"] != "silent") & (df["n_sites"] > 0)].copy()
    fam_e = exp.dropna(subset=["seed"]).groupby("seed")["arm"].nunique()
    multi = set(fam_e[fam_e >= 2].index)                 # seeds with ≥2 expressed members
    exp["is_multi"] = exp["seed"].isin(multi)
    n_single = int((~exp["is_multi"]).sum())
    smap = exp.set_index("arm")["seed"].to_dict()
    is_multi = exp.set_index("arm")["is_multi"].to_dict()
    es = sites[sites["miRNA"].isin(exp["arm"]) | sites["miRNA"].map(_strip).isin(exp["arm"].map(_strip))].copy()
    es["seed"] = es["miRNA"].map(lambda a: smap.get(a, smap.get(_strip(a))))
    # order: multi-member families first (family size desc, seed, abundance), then single-seed (abundance)
    exp["_fam"] = exp["seed"].map(fam_e).where(exp["is_multi"], 0)
    multi_order = exp[exp["is_multi"]].sort_values(["_fam", "seed", "med"], ascending=[False, True, False])["arm"].tolist()
    single_order = exp[~exp["is_multi"]].sort_values("med", ascending=False)["arm"].tolist()
    order = multi_order + single_order
    n_multi = len(multi_order)
    ypos = {a: i for i, a in enumerate(order)}
    seeds_sorted = list(dict.fromkeys(exp.set_index("arm").loc[multi_order, "seed"]))
    scol = {s: plt.cm.tab20(i % 20) for i, s in enumerate(seeds_sorted)}
    def acol(a):                                          # multi → seed colour; single-seed → grey
        if is_multi.get(a, is_multi.get(_strip(a), False)):
            return scol.get(smap.get(a, smap.get(_strip(a))), "0.62")
        return "0.62"
    # site-pair regimes between DIFFERENT miRNAs: S2 overlap (≤8 nt, mutual exclusion) vs S1 proximal (8–40 nt)
    ps = es.sort_values("utr_pos").reset_index(drop=True)
    for i in range(len(ps) - 1):
        if ps.at[i, "miRNA"] == ps.at[i + 1, "miRNA"]:
            continue
        x0, x1 = ps.at[i, "utr_pos"], ps.at[i + 1, "utr_pos"]
        gap = x1 - x0
        if gap <= _OVERLAP_WIN:                                       # S2 — overlapping → competition
            axC.axvspan(x0 - 5, x1 + 5, color="#e8a39c", alpha=0.55, zorder=0)
        elif gap <= _COMP_WIN:                                        # S1 — proximal → cooperative
            axC.axvspan(x0 - 2, x1 + 2, color="#bcd2ec", alpha=0.5, zorder=0)
    for r in es.itertuples():
        a = r.miRNA if r.miRNA in ypos else _strip(r.miRNA)
        y = ypos.get(a, ypos.get(_strip(r.miRNA)))
        if y is None:
            continue
        ms = {"8mer": 60, "7mer-m8": 38, "7mer-A1": 22}.get(r.type, 25)
        c = acol(a)
        axC.scatter(r.utr_pos, y, s=ms, color=c, edgecolors="0.3", linewidths=0.3, zorder=3)
        axC.plot([r.utr_pos, r.utr_pos], [-0.6, y], color=c, lw=0.3, alpha=0.4, zorder=1)
    axC.axhline(-0.6, color="0.4", lw=4, solid_capstyle="round", zorder=2)   # the 3'UTR
    if n_single and n_multi:                                                 # divider: families (coloured) ↓ | single-seed (grey) ↑
        axC.axhline(n_multi - 0.5, color="0.7", lw=0.8, ls=":", zorder=2)
        axC.text(utr_len * 1.0, n_multi - 0.5, "single-seed (grey) ↑\nseed families ↓", fontsize=6.5,
                 color="0.5", ha="right", va="center", style="italic", zorder=7)
    for a in order:
        axC.text(-utr_len * 0.012, ypos[a], a.replace("hsa-", ""), ha="right", va="center", fontsize=6, color=acol(a))
    axC.set_xlim(-utr_len * 0.14, utr_len * 1.02); axC.set_ylim(-1.5, len(order))
    axC.set_xlabel(f"{gene} 3′UTR position (nt; total {utr_len:,})", fontsize=9)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    leg = [Line2D([], [], marker="o", color="w", markerfacecolor="0.45", markeredgecolor="0.3",
                  markersize=np.sqrt(s), label=lab) for s, lab in [(60, "8mer"), (38, "7mer-m8"), (22, "7mer-A1")]]
    leg += [Patch(facecolor="#e8a39c", alpha=0.7, label=f"overlap ≤{_OVERLAP_WIN} nt — S2 mutual exclusion"),
            Patch(facecolor="#bcd2ec", alpha=0.7, label=f"proximal {_OVERLAP_WIN}–{_COMP_WIN} nt — S1 cooperative")]
    axC.legend(handles=leg, loc="lower right", fontsize=7, ncol=2, framealpha=1.0,
               title="dot size = site type   ·   band = site-pair regime", title_fontsize=7.5)
    axC.set_yticks([]); axC.set_title(f"C · {gene} 3′UTR competition map — all {len(order)} expressed+sited regulators: "
        f"{n_multi} in {len(multi)} multi-member seed families (coloured = direct seed competition) + "
        f"{n_single} single-seed (grey = positional competition only); dot=site (size∝type)", fontsize=11, loc="left")

    # ---- D: position-ordered zoom on the densest competition hotspot ------- #
    axD = fig.add_subplot(gs[4, :])
    W, fw = 90, {"8mer": 8, "7mer-m8": 7, "7mer-A1": 7}
    starts = sorted(es["utr_pos"].unique())
    best = max(starts, key=lambda p: es[(es["utr_pos"] >= p) & (es["utr_pos"] <= p + W)]["seed"].nunique())
    zoom = es[(es["utr_pos"] >= best - 3) & (es["utr_pos"] <= best + W + 3)].sort_values(
        ["utr_pos", "seed", "miRNA"]).reset_index(drop=True)
    cov = {}
    for r in zoom.itertuples():
        for xnt in range(int(r.utr_pos), int(r.utr_pos) + fw.get(r.type, 7)):
            cov[xnt] = cov.get(xnt, 0) + 1
    over = sorted(x for x, c in cov.items() if c >= 2)
    if over:                                              # shade merged overlap runs (≥2 footprints share nt)
        runs, s0, p0 = [], over[0], over[0]
        for xnt in over[1:]:
            if xnt == p0 + 1:
                p0 = xnt
            else:
                runs.append((s0, p0)); s0, p0 = xnt, xnt
        runs.append((s0, p0))
        for a, b in runs:
            axD.axvspan(a, b + 1, color="#f4c7c3", alpha=0.55, zorder=0)
    for i, r in enumerate(zoom.itertuples()):
        w = fw.get(r.type, 7)
        axD.barh(i, w, left=r.utr_pos, height=0.66, color=scol.get(r.seed, "0.5"), edgecolor="0.25",
                 linewidth=0.4, zorder=3)
        axD.text(r.utr_pos - 1.5, i, r.miRNA.replace("hsa-", ""), ha="right", va="center", fontsize=7.5,
                 color=scol.get(r.seed, "0.3"))
        axD.text(r.utr_pos + w + 1.5, i, r.seed, ha="left", va="center", fontsize=6.5, family="monospace", color="0.4")
    axD.set_xlim(best - 27, best + W + 16); axD.set_ylim(-1, len(zoom)); axD.set_yticks([])
    axD.set_xlabel(f"{gene} 3′UTR position (nt) — bar = seed-match footprint (~7–8 nt); pink = ≥2 footprints overlap (S2)", fontsize=9)
    axD.set_title(f"D · competition-hotspot zoom @ {best}–{best + W} nt — {zoom['seed'].nunique()} distinct seeds / "
                  f"{zoom['miRNA'].nunique()} miRNAs overlapping one window (scattered across families in C)", fontsize=11, loc="left")
    from matplotlib.patches import Rectangle
    axC.add_patch(Rectangle((best - 4, -0.7), W + 8, len(order) + 0.1, fill=False, edgecolor="#c0392b",
                            lw=1.4, ls="--", zorder=6))
    axC.annotate("zoom → D", (best + W / 2, len(order) - 0.3), fontsize=8.5, color="#c0392b", ha="center",
                 va="top", fontweight="bold")

    fig.suptitle(f"{gene} — miRNA competition for a shared 3′UTR", fontsize=14, fontweight="bold", y=0.995)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{gene}_competition.png", dpi=145, bbox_inches="tight"); plt.close(fig)
    print(f"[competition] {gene}: A/B/C/D panels — {len(df)} regulators, {len(exp)} expressed+sited "
          f"({n_multi} in seed families + {n_single} single-seed), {fam_e.size} seeds; wrote figures/{gene}_competition.png")
    from mirna_hallmark.method_dev.landscape.he_edge_arm_landscape import build_gene_hub
    build_gene_hub(gene)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "PTEN")
