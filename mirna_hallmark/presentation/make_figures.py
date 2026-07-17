#!/usr/bin/env python3
"""Generate / collect all figures for the miRNA-Hallmark WIP talk.

Pure pandas + matplotlib (no seaborn/networkx). Each figure is independent and
wrapped so one failure never blocks the rest. Run from repo root:

    .venv/bin/python3 -m mirna_hallmark.presentation.make_figures
    # or: .venv/bin/python3 mirna_hallmark/presentation/make_figures.py

Writes PNGs to mirna_hallmark/presentation/figures/ and copies the pre-rendered
pipeline figures (grant heatmap, CNV boxplot, NMF heatmap, edge-transition,
Buffa) into the same folder so the deck references one directory.
"""
from __future__ import annotations
import shutil
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "mirna_hallmark" / "output"
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(parents=True, exist_ok=True)

NEG = "#c0392b"   # repression / damage
POS = "#2471a3"   # positive / pro-tumor
GREY = "#7f8c8d"
ACC = "#16a085"
plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 140, "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False})


def _save(fig, name):
    p = FIG / name
    fig.tight_layout()
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {name}")


def _tsv(rel, **kw):
    return pd.read_csv(OUT / rel, sep="\t", **kw)


def figure(fn):
    """decorator: run, catch, report."""
    def wrap():
        try:
            print(f"[{fn.__name__}]")
            fn()
        except Exception as e:  # noqa
            print(f"  !! {fn.__name__} FAILED: {e}")
            traceback.print_exc()
    wrap.__name__ = fn.__name__
    return wrap


# ----------------------------------------------------------------------------- copies
@figure
def copies():
    pairs = {
        "F12a_basal_heatmap.png": "figures/grant_figure_basal_heatmap.png",
        "F13_cnv_boxplot.png": "mirna_locus_cnv/figures/mirna_cnv_expr_concordance_rho_boxplot.png",
        "F17_edge_transition.png": "tissue_reference/edge_transition_pressure_panels/figures/edge_transition_pressure_panel_page01.png",
        "F20_nmf_basal.png": "tissue_reference/mirna_comovement/nmf_sample_signatures/figures/gene_pressure_acquired_within_Basal.png",
        "F27_buffa.png": "tissue_reference/buffa_validation/figures/buffa_vs_tcga_replication.png",
    }
    # buffa path variant
    alt_buffa = OUT / "buffa_validation/figures/buffa_vs_tcga_replication.png"
    if alt_buffa.exists():
        pairs["F27_buffa.png"] = "buffa_validation/figures/buffa_vs_tcga_replication.png"
    for dst, src in pairs.items():
        s = OUT / src
        if s.exists():
            shutil.copy(s, FIG / dst)
            print(f"  copied {dst}")
        else:
            print(f"  (skip, missing) {src}")


# ----------------------------------------------------------------------------- F4 edges
@figure
def f4_edges_funnel():
    e = _tsv("edges/mirna_hallmark_edges.tsv.gz")
    pair = e.drop_duplicates(["miRNA", "gene"])
    he = pair[pair["high_evidence"] == True]  # noqa
    cats = ["(miRNA,gene)\nrows", "unique\n(miRNA,gene)", "high-evidence\npairs"]
    vals = [len(e), len(pair), len(he)]
    fig, ax = plt.subplots(figsize=(6, 3.6))
    bars = ax.bar(cats, vals, color=[GREY, POS, NEG])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("count")
    ax.set_title(f"Hallmark-scoped miRTarBase edge map  ·  {pair['gene'].nunique():,} genes, "
                 f"{pair['miRNA'].nunique():,} arms")
    _save(fig, "F4_edges_funnel.png")


# ----------------------------------------------------------------------------- F5 pressure overview
@figure
def f5_pressure_overview():
    pmat = _tsv("hallmark_interaction/hallmark_pressure_per_sample.tsv.gz").set_index(
        _tsv("hallmark_interaction/hallmark_pressure_per_sample.tsv.gz").columns[0])
    vals = pmat.to_numpy(dtype=float).ravel()
    vals = vals[np.isfinite(vals)]
    e = _tsv("edges/mirna_hallmark_edges.tsv.gz").drop_duplicates(["miRNA", "gene"])
    reg_per_gene = e.groupby("gene")["miRNA"].nunique()
    wt = _tsv("robustness/aim2_circularity/mirna_basal_load_by_weighting.tsv")
    load = wt["load_evidence"].sort_values(ascending=False).to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    axs[0].hist(vals, bins=60, color=POS, alpha=.85)
    axs[0].set_title("(a) hallmark pressure values"); axs[0].set_xlabel("gated pressure"); axs[0].set_ylabel("count")
    axs[1].hist(reg_per_gene, bins=40, color=ACC, alpha=.85)
    axs[1].set_title("(b) regulators per gene"); axs[1].set_xlabel("# miRNA arms targeting gene")
    axs[1].axvline(reg_per_gene.median(), color=NEG, ls="--", lw=1, label=f"median {reg_per_gene.median():.0f}")
    axs[1].legend(fontsize=8)
    csum = np.cumsum(load) / load.sum()
    axs[2].plot(np.arange(1, len(load) + 1), csum, color=NEG)
    axs[2].set_title("(c) miRNA load concentration"); axs[2].set_xlabel("miRNA arms (ranked)")
    axs[2].set_ylabel("cumulative basal load")
    k = int((csum < .8).sum())
    axs[2].axhline(.8, color=GREY, ls=":", lw=1)
    axs[2].text(len(load) * .5, .5, f"{k} arms → 80% of load", fontsize=8, color=NEG)
    _save(fig, "F5_pressure_overview.png")


# ----------------------------------------------------------------------------- F6 AGO gate
@figure
def f6_ago_gate():
    a = _tsv("ago_gate/per_sample_ago_capacity.tsv.gz").sort_values("ago_capacity_z")
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.scatter(a["ago_capacity_z"], a["ago_gate"], s=8, color=POS, alpha=.5)
    ax.set_xlabel("RISC capacity  (mean AGO z-score)")
    ax.set_ylabel("gate multiplier")
    ax.set_title("AGO/RISC saturating gate  ∈ [0.5, 1]")
    ax.axhline(0.75, color=GREY, ls=":", lw=1)
    ax.axvline(0, color=GREY, ls=":", lw=1)
    ax.text(a["ago_capacity_z"].min(), 0.76, "0.75 at cohort mean", fontsize=8, color=GREY)
    _save(fig, "F6_ago_gate.png")


# ----------------------------------------------------------------------------- F9 enrichment lollipop
@figure
def f9_enrichment():
    d = _tsv("hallmark_interaction/hallmark_highev_target_enrichment.tsv").copy()
    d = d.sort_values("fold_enrichment")
    d["lab"] = d["hallmark_set"].str.replace("HALLMARK_", "", regex=False).str.replace("_", " ").str.title()
    show = pd.concat([d.head(8), d.tail(10)])
    col = np.where(show["hypergeom_q"] < 0.05, NEG, GREY)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.hlines(show["lab"], 1, show["fold_enrichment"], color=col, lw=2, alpha=.6)
    ax.scatter(show["fold_enrichment"], show["lab"], color=col, s=40, zorder=3)
    ax.axvline(1, color="k", lw=.8)
    ax.set_xlabel("fold enrichment for high-evidence targets")
    ax.set_title("Where do miRNAs aim?  (hypergeometric, top & bottom programs)")
    _save(fig, "F9_enrichment.png")


# ----------------------------------------------------------------------------- F10b miR pressure fan
@figure
def f10b_mir_fan():
    e = _tsv("edges/mirna_hallmark_edges.tsv.gz").drop_duplicates(["miRNA", "gene"])
    arm = "hsa-miR-17-5p"
    sub = e[e["miRNA"] == arm].copy()
    if sub.empty:
        cand = e["miRNA"].value_counts().index[0]
        arm, sub = cand, e[e["miRNA"] == cand].copy()
    sub = sub.sort_values("evidence_score", ascending=False).head(22)
    n = len(sub)
    ang = np.linspace(0.2, np.pi - 0.2, n)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    sc = sub["evidence_score"].to_numpy(float)
    w = 0.8 + 3.5 * sc / sc.max()
    for a, g, ww in zip(ang, sub["gene"], w):
        x, y = np.cos(a), np.sin(a)
        ax.plot([0, x], [0, y], color=NEG, lw=ww, alpha=.5, solid_capstyle="round")
        ax.text(1.04 * x, 1.04 * y, g, ha="left" if x >= 0 else "right",
                va="center", rotation=np.degrees(a) - 90, fontsize=7)
    ax.scatter([0], [0], s=420, color=POS, zorder=5)
    ax.text(0, -0.08, arm, ha="center", va="top", fontsize=10, weight="bold")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-0.2, 1.35); ax.axis("off")
    ax.set_title(f"Pressure fan: {arm} → hallmark targets (edge width = evidence)")
    _save(fig, "F10b_mir_fan.png")


# ----------------------------------------------------------------------------- F11 coupling volcano
@figure
def f11_coupling_volcano():
    d = _tsv("hallmark_interaction/hallmark_anticorrelation.tsv")
    d = d[d["view"] == "gated"].copy()
    d["nlq"] = -np.log10(d["partial_q"].clip(lower=1e-300))
    col = np.where((d["partial_rho"] < 0) & (d["partial_q"] < .05), NEG,
                   np.where((d["partial_rho"] > 0) & (d["partial_q"] < .05), POS, GREY))
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.scatter(d["partial_rho"], d["nlq"], c=col, s=28, alpha=.8)
    ax.axvline(0, color="k", lw=.7); ax.axhline(-np.log10(.05), color=GREY, ls=":", lw=1)
    nneg = int(((d["partial_rho"] < 0) & (d["partial_q"] < .05)).sum())
    ax.set_xlabel("partial Spearman ρ  (pressure vs target expr, CPE+HRD-adj)")
    ax.set_ylabel("−log10 q")
    ax.set_title(f"Does pressure repress?  {nneg}/50 hallmarks neg & FDR-sig")
    _save(fig, "F11_coupling_volcano.png")


# ----------------------------------------------------------------------------- F12b PTEN map
@figure
def f12b_pten_map():
    d = _tsv("pam50_gene_resolution/gene_mirna_drivers_by_pam50.tsv")
    p = d[(d["gene"] == "PTEN") & (d["subtype"] == "Basal")].copy()
    if p.empty:
        p = d[d["gene"] == "PTEN"].copy()
    p = p.sort_values("structural_weight", ascending=False).head(12)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.barh(p["miRNA"][::-1], p["structural_weight"][::-1], color=NEG, alpha=.8)
    ax.set_xlabel("structural pressure weight on PTEN (Basal)")
    ax.set_title("Gene pressure map: PTEN regulators")
    _save(fig, "F12b_pten_map.png")


# ----------------------------------------------------------------------------- F15 acquired volcano
@figure
def f15_acquired_volcano():
    s = _tsv("tissue_reference/mirna_state_class/mirna_state_class.tsv")
    s = s.dropna(subset=["log2fc_tumor_vs_healthy", "dHT"])
    cmap = {"rank+magnitude": ACC, "magnitude_only": NEG, "rank_only": POS, "none": GREY}
    col = s["acquired_axis"].map(cmap).fillna(GREY)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.scatter(s["log2fc_tumor_vs_healthy"], s["dHT"], c=col, s=18, alpha=.7)
    ax.axhline(0, color="k", lw=.6); ax.axvline(0, color="k", lw=.6)
    ax.set_xlabel("log2FC tumor vs healthy  (magnitude axis)")
    ax.set_ylabel("dHT  (rank axis, saturates)")
    ng = int(s["acquired_gainer"].sum()) if "acquired_gainer" in s else 0
    ax.set_title(f"The acquired program  ·  {ng}/{len(s)} arms gain on ≥1 axis")
    for lab, c in cmap.items():
        ax.scatter([], [], color=c, label=lab, s=18)
    ax.legend(fontsize=7, loc="lower right")
    # label a few oncomiRs
    for m in ["hsa-miR-21-5p", "hsa-miR-182-5p", "hsa-miR-183-5p", "hsa-miR-375"]:
        r = s[s["miRNA"] == m]
        if not r.empty:
            ax.annotate(m.replace("hsa-", ""), (r["log2fc_tumor_vs_healthy"].iloc[0], r["dHT"].iloc[0]),
                        fontsize=7, color="k")
    _save(fig, "F15_acquired_volcano.png")


# ----------------------------------------------------------------------------- F18 decoupling ladder
@figure
def f18_decoupling():
    d = _tsv("tissue_reference/decoupling_validation/decoupling_validation.tsv")
    rows = [
        ("copy number", d["n_cn"].notna().sum(), int(d["cnv_masked"].sum())),
        ("promoter meth", d["n_meth"].notna().sum(), int(d["meth_masked"].sum())),
        ("ATAC (direct)", d["n_atac_direct"].notna().sum(), int(d["atac_direct_masked"].sum())),
    ]
    labels = [r[0] for r in rows]; tested = [r[1] for r in rows]; masked = [r[2] for r in rows]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.barh(y, tested, color="#d5dbdb", label="tested")
    ax.barh(y, masked, color=NEG, label="masked (confound)")
    for yi, t, m in zip(y, tested, masked):
        ax.text(t + 1, yi, f"{m}/{t}", va="center", fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("lost/decoupled edges")
    ax.set_title("Is de-coupling a confound?  →  no (CN/meth) , tiny ATAC minority")
    ax.legend(fontsize=8, loc="lower right")
    _save(fig, "F18_decoupling.png")


# ----------------------------------------------------------------------------- F21 architecture
@figure
def f21_architecture():
    d = _tsv("tissue_reference/geneset_architecture/architecture_all_set_summary.tsv").copy()
    d["lab"] = d["hallmark"].str.replace("HALLMARK_", "", regex=False).str.replace("_", " ").str.title()
    d = d.dropna(subset=["total_pro_tumor_pressure"])
    d = d[d["total_pro_tumor_pressure"] != 0].sort_values("total_pro_tumor_pressure")
    show = pd.concat([d.head(8), d.tail(10)])
    col = np.where(show["total_pro_tumor_pressure"] > 0, POS, NEG)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(show["lab"], show["total_pro_tumor_pressure"], color=col, alpha=.85)
    ax.axvline(0, color="k", lw=.8)
    ax.set_xlabel("net pro-tumor pressure   (+ brake-release   /   − engine-damage)")
    ax.set_title("Gene-set architecture: the brake-release thesis")
    _save(fig, "F21_architecture.png")


# ----------------------------------------------------------------------------- F23 hub ladder
@figure
def f23_hub_ladder():
    d = _tsv("robustness/aim1_proliferation/hub_route_partial_corr.tsv")
    scope = "Basal" if (d["scope"] == "Basal").any() else d["scope"].iloc[0]
    d = d[d["scope"] == scope]
    steps = ["raw_rho", "rho_CPE_HRD", "rho_CPE_HRD_CN", "rho_e2f_g2m_CN"]
    xlab = ["raw", "+CPE/HRD", "+target CN", "+prolif"]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    palette = [NEG, POS, ACC, "#8e44ad"]
    for i, tgt in enumerate(["PTEN", "CDKN1A", "TGFBR2", "VIM"]):
        sub = d[d["target"] == tgt]
        if sub.empty:
            continue
        row = sub.loc[sub["rho_CPE_HRD_CN"].idxmin()]
        ys = [row[s] for s in steps]
        ax.plot(xlab, ys, "-o", color=palette[i % 4], label=f"{tgt} ← {row['miRNA'].replace('hsa-','')}")
    ax.axhline(0, color="k", lw=.6)
    ax.set_ylabel("partial Spearman ρ")
    ax.set_title(f"Hub routes survive the confound ladder ({scope})")
    ax.legend(fontsize=8)
    _save(fig, "F23_hub_ladder.png")


# ----------------------------------------------------------------------------- F24 robustness 2-panel
@figure
def f24_robustness():
    fam = _tsv("robustness/aim2_circularity/family_load_share.tsv").sort_values("share_evidence_pct", ascending=False).head(8)
    wt = _tsv("robustness/aim2_circularity/mirna_basal_load_by_weighting.tsv")
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    axs[0].barh(fam["family"][::-1], fam["share_evidence_pct"][::-1], color=ACC, alpha=.85)
    axs[0].set_xlabel("% of Basal load"); axs[0].set_title("(a) family load share")
    axs[1].scatter(wt["rank_evidence"], wt["rank_targetscan"], s=12, color=POS, alpha=.5)
    rho = wt[["rank_evidence", "rank_targetscan"]].corr(method="spearman").iloc[0, 1]
    axs[1].set_xlabel("rank (evidence weighting)"); axs[1].set_ylabel("rank (TargetScan weighting)")
    axs[1].set_title(f"(b) rank concordance  ρ≈{rho:.2f}")
    _save(fig, "F24_robustness.png")


# ----------------------------------------------------------------------------- F25 CPTAC protein
@figure
def f25_cptac():
    """MH-114. The bulk-protein miRNA signal is ~57% COMPARTMENT ARITHMETIC.

    The old version of this figure was a volcano of the UNADJUSTED partial ρ (purity+CIN only) titled
    "Independent protein validation: 53 genes FDR-neg", annotating ZEB1/PTEN/**CALD1**/**MMP2** — and CALD1
    (caldesmon) and MMP2 are textbook CAF/myofibroblast markers. Bulk protein and bulk miRNA are both
    compartment-weighted averages, so a miRNA merely EXPRESSED IN A DIFFERENT CELL TYPE than the target
    anti-correlates with it across patients with no cell-autonomous repression. Proof (`eval/
    protein_compartment_null`): degree-preserving SHUFFLED non-edges — pairs that CANNOT repress — show the
    SAME compartment gradient as curated edges (−0.129 vs −0.127). So the honest figure SHOWS that, rather
    than quietly reporting a smaller number.
    """
    s = _tsv("learned/protein_compartment_null_summary.tsv")
    d = _tsv("cptac_validation/prospective/gene_level_associations.tsv.gz")
    d = d[(d["variant"] == "highev|gated") & (d["layer"] == "protein_anticorr")].copy()

    fig, axs = plt.subplots(1, 3, figsize=(13.2, 4.2))
    # (a,b) the shuffled null, before and after adjusting for cell composition
    for ax, col, ttl in [(axs[0], "rho_base", "(a) bulk protein, NOT adjusted"),
                         (axs[1], "rho_adj", "(b) adjusted for cell composition")]:
        w, x = 0.36, np.arange(2)
        for k, (grp, c) in enumerate([("REAL_HE", NEG), ("SHUFFLED", GREY)]):
            v = [float(s[(s.group == grp) & (s.orientation == o)][col].iloc[0]) for o in ("OPPOSITE", "SAME")]
            ax.bar(x + (k - .5) * w, v, w, color=c, label="curated edges" if k == 0 else "SHUFFLED non-edges")
        ax.axhline(0, color="k", lw=.7)
        ax.set_xticks(x); ax.set_xticklabels(["opposite\ncompartment", "same\ncompartment"])
        ax.set_ylim(-.09, .09); ax.set_ylabel("mean ρ (miRNA vs protein)"); ax.set_title(ttl, fontsize=10)
    axs[0].legend(fontsize=8, frameon=False, loc="lower right")
    axs[0].annotate("non-edges CANNOT repress —\nyet they look identical", (0, -.074), fontsize=7.5,
                    color=GREY, ha="center", va="top")

    # (c) the honest volcano: composition-adjusted
    d["nlq"] = -np.log10(d["partial_q"].clip(lower=1e-300))
    col = np.where((d["partial_rho"] < 0) & (d["partial_q"] < .05), NEG, GREY)
    axs[2].scatter(d["partial_rho"], d["nlq"], c=col, s=14, alpha=.7)
    axs[2].axvline(0, color="k", lw=.6); axs[2].axhline(-np.log10(.05), color=GREY, ls=":", lw=1)
    nneg = int(((d["partial_rho"] < 0) & (d["partial_q"] < .05)).sum())
    for g in ["ZEB1", "PTEN"]:
        r = d[d["gene_name"] == g]
        if not r.empty:
            axs[2].annotate(g, (r["partial_rho"].iloc[0], r["nlq"].iloc[0]), fontsize=8, weight="bold")
    axs[2].set_xlabel("partial ρ (composition-adjusted)"); axs[2].set_ylabel("−log10 q")
    axs[2].set_title(f"(c) what survives: {nneg} genes FDR-neg", fontsize=10)

    fig.suptitle("CPTAC-2 protein: ~57% of the apparent repression is COMPARTMENT ARITHMETIC "
                 "(real effect survives as a set, p=2e-3)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, .93))
    _save(fig, "F25_cptac.png")


# ----------------------------------------------------------------------------- schematics
def _box(ax, x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.03",
                                fc=fc, ec="#34495e", lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)


def _arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                 color="#34495e", lw=1.4))


@figure
def f2_entity_flow():
    fig, ax = plt.subplots(figsize=(11, 2.6)); ax.axis("off")
    ax.set_xlim(0, 11); ax.set_ylim(0, 2.6)
    nodes = [(0.1, "locus\n(CN)", "#fdebd0"), (2.0, "arm\nexpression", "#fdebd0"),
             (4.0, "pressure\nper gene", "#d6eaf8"), (6.0, "× AGO gate", "#d5f5e3"),
             (8.0, "hallmark\nroll-up", "#d6eaf8"), (9.9, "target\nexpression", "#fadbd8")]
    for x, t, c in nodes:
        _box(ax, x, 1.0, 1.6, 0.9, t, c)
    for i in range(len(nodes) - 1):
        _arrow(ax, nodes[i][0] + 1.6, 1.45, nodes[i + 1][0], 1.45)
    _arrow(ax, 8.8, 1.0, 9.9, 1.0)  # roll-up vs target (coupling)
    ax.text(9.0, 0.6, "coupling (↕, expect −)", fontsize=8, color=NEG)
    ax.set_title("From locus to coupling: the pressure pipeline", fontsize=11)
    _save(fig, "F2_entity_flow.png")


@figure
def f8_pipeline_dag():
    fig, ax = plt.subplots(figsize=(11, 3.4)); ax.axis("off")
    ax.set_xlim(0, 11); ax.set_ylim(0, 3.4)
    core = [(0.1, "build_edges"), (2.1, "ago_gate"), (4.1, "hallmark_\ninteraction"),
            (6.1, "pam50_gene_\nresolution"), (8.1, "robustness_\nchecks")]
    for x, t in core:
        _box(ax, x, 2.2, 1.7, 0.8, t, "#d6eaf8")
    for i in range(len(core) - 1):
        _arrow(ax, core[i][0] + 1.7, 2.6, core[i + 1][0], 2.6)
    land = [(2.1, "cross_state_\nlandscape"), (4.1, "mirna_state_\nclass (R/P/0)"),
            (6.1, "decoupling_\nvalidation"), (8.1, "geneset_\narchitecture")]
    for x, t in land:
        _box(ax, x, 0.9, 1.7, 0.8, t, "#d5f5e3")
    for i in range(len(land) - 1):
        _arrow(ax, land[i][0] + 1.7, 1.3, land[i + 1][0], 1.3)
    _arrow(ax, 4.95, 2.2, 4.95, 1.7)
    val = [(0.1, "CPTAC\nprotein"), (0.1, "")]
    _box(ax, 0.1, 0.9, 1.7, 0.8, "CPTAC / Buffa\nvalidation", "#fadbd8")
    ax.text(0.1, 3.05, "core spine", color=POS, fontsize=9)
    ax.text(2.1, 1.75, "landscape (tumor / NAT / GTEx)", color=ACC, fontsize=9)
    ax.set_title("Pipeline at a glance", fontsize=11)
    _save(fig, "F8_pipeline_dag.png")


@figure
def f16_rp0():
    fig, ax = plt.subplots(figsize=(9, 3.2)); ax.axis("off")
    ax.set_xlim(0, 9); ax.set_ylim(0, 3.2)
    for i, st in enumerate(["GTEx\n(healthy)", "NAT", "tumor"]):
        _box(ax, 0.4 + i * 1.4, 2.2, 1.2, 0.7, st, "#eaf2f8")
    ax.text(0.2, 1.7, "each state →  R (sig repressive ρ<−0.1, q<.05)  ·  P (sig positive)  ·  0 (neither)",
            fontsize=9)
    examples = [("0 0 R", "acquired_realized — switched on only in tumor", NEG),
                ("R 0 0", "released brake — held in healthy, gone in tumor", POS),
                ("0 R R", "field_established — set in NAT, kept in tumor", ACC)]
    for j, (code, txt, c) in enumerate(examples):
        ax.text(0.4, 1.25 - j * 0.42, code, fontsize=11, weight="bold", color=c, family="monospace")
        ax.text(1.7, 1.25 - j * 0.42, txt, fontsize=9, color="#2c3e50")
    ax.set_title("The R / P / 0 coupling-trajectory code  (letters = GTEx / NAT / tumor)", fontsize=11)
    _save(fig, "F16_rp0.png")


@figure
def f19_nmf():
    fig, ax = plt.subplots(figsize=(9, 3.0)); ax.axis("off")
    ax.set_xlim(0, 9); ax.set_ylim(0, 3.0)
    _box(ax, 0.3, 1.0, 1.8, 1.2, "M\npressure\n(genes×samples)", "#fdebd0")
    ax.text(2.4, 1.6, "≈", fontsize=20)
    _box(ax, 3.0, 1.0, 1.4, 1.2, "W\n(programs)", "#d6eaf8")
    ax.text(4.6, 1.6, "×", fontsize=18)
    _box(ax, 5.2, 1.0, 1.4, 1.2, "H\n(activity\n×sample)", "#d5f5e3")
    ax.text(0.3, 0.6, "menu: aggregate · within-tumor · per-subtype · gene-pressure · signed gain/loss · share",
            fontsize=8.5, color="#2c3e50")
    ax.set_title("NMF: factor pressure into co-varying programs", fontsize=11)
    _save(fig, "F19_nmf.png")


def main():
    print(f"figures → {FIG}")
    for fn in [copies, f4_edges_funnel, f5_pressure_overview, f6_ago_gate, f9_enrichment,
               f10b_mir_fan, f11_coupling_volcano, f12b_pten_map, f15_acquired_volcano,
               f18_decoupling, f21_architecture, f23_hub_ladder, f24_robustness, f25_cptac,
               f2_entity_flow, f8_pipeline_dag, f16_rp0, f19_nmf]:
        fn()
    print("done.")


if __name__ == "__main__":
    main()
