"""Seed-grant figures for the miRNA x Hallmark subproject.

Two panels, both reading only existing robustness outputs (safe to re-run):

  Fig. 1 (spider) -- the 8 key Hallmark programs, one trace per PAM50 subtype.
    The radial axis is the **anti-correlation** (-rho) of miRNA Hallmark pressure
    vs target expression (partial Spearman, purity/HRD/proliferation/CN adjusted,
    within-subtype). Reverted from raw rho so that the **strongest repression sits
    furthest from the centre** -- Basal is the outer envelope on nearly every axis,
    making subtype heterogeneity in miRNA pressure the immediate read.

  Fig. 2 (heatmap) -- the miR-17~92 / 106b~25 cluster's pressure->repression
    coupling for the same 8 programs across subtypes. One miRNA family applies
    monotonically deeper pressure LumA -> LumB -> Her2 -> Basal (e.g. p53 pathway
    -0.22 -> -0.31 -> -0.44 -> -0.58): the LumA/LumB deep-dive difference in a
    single image.

    .venv/bin/python3 -m mirna_hallmark.make_seed_grant_figures
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from mirna_hallmark import config as C

ROBUST = C.OUTPUT_ROOT / "robustness" / "aim1_proliferation"
# Fig. 1 uses the high-confidence (evidence + TargetScan) pressure universe
# (mirna_hallmark.build_highconf_pressure); falls back to the full-aggregate table
# if the high-confidence one has not been built.
BY_SUB = ROBUST / "hallmark_coupling_by_pam50_highconf.tsv"
if not BY_SUB.exists():
    BY_SUB = ROBUST / "hallmark_coupling_by_pam50_prolif.tsv"
FAM_SUB = ROBUST / "family_coupling_by_pam50.tsv"
OUT_DIR = C.OUTPUT_ROOT / "figures"

SUBTYPES = ["LumA", "LumB", "Her2", "Basal"]
RHO_COL = "rho_prolif_cn_wsd_adj"  # purity+HRD+proliferation+CN+within-subtype-z
CLUSTER = "miR-17~92/106b~25/106a~363"

SUB_COLOR = {"LumA": "#4575b4", "LumB": "#74add1", "Her2": "#fdae61", "Basal": "#d73027"}

NAME = {
    "HALLMARK_APOPTOSIS": "Apoptosis",
    "HALLMARK_P53_PATHWAY": "p53 pathway",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING": "PI3K\u2013AKT\u2013mTOR",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
    "HALLMARK_TGF_BETA_SIGNALING": "TGF-\u03b2",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF\u03b1/NF-\u03baB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-\u03b3 response",
    "HALLMARK_NOTCH_SIGNALING": "Notch",
}
# Display order around the spider (tumour-suppressor / cell-state grouped).
ORDER = [
    "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_NOTCH_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]


def _key_matrix(path: Path, *, family: str | None = None) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t")
    d = d[d["key_hallmark"] == True]  # noqa: E712
    if family is not None:
        d = d[d["family"] == family]
    mat = d.drop_duplicates(["hallmark_set", "subtype"]).pivot(
        index="hallmark_set", columns="subtype", values=RHO_COL
    )
    return mat.reindex(index=[h for h in ORDER if h in mat.index], columns=SUBTYPES)


def _spider() -> Path:
    """8 key Hallmarks; radial = -rho (anti-correlation); outer = stronger repression."""
    mat = _key_matrix(BY_SUB)
    anti = -mat  # revert gradient: strongest repression -> largest radius
    labels = [NAME.get(h, h) for h in anti.index]
    n = len(labels)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    ring = angles + angles[:1]

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(5.4, 5.4), subplot_kw={"polar": True})
    rmax = float(np.nanmax(anti.values)) * 1.05
    for st in SUBTYPES:
        vals = anti[st].astype(float).tolist()
        vals += vals[:1]
        lw = 2.4 if st == "Basal" else 1.5
        ax.plot(ring, vals, label=st, linewidth=lw, color=SUB_COLOR[st])
        ax.fill(ring, vals, alpha=0.07, color=SUB_COLOR[st])
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, rmax)
    ax.set_rlabel_position(90)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.set_title(
        "High-confidence miRNA pressure \u2192 repression by PAM50 subtype\n"
        "(8 key Hallmark programs; evidence+TargetScan edges; radius = anti-\n"
        "correlation \u2212\u03c1, further from centre = stronger repression)",
        fontsize=9, y=1.10,
    )
    ax.legend(loc="upper right", bbox_to_anchor=(1.30, 1.12), fontsize=8, title="PAM50")
    fig.tight_layout()
    pdf = OUT_DIR / "seed_spider_keyhallmarks.pdf"
    fig.savefig(pdf)
    fig.savefig(OUT_DIR / "seed_spider_keyhallmarks.png", dpi=600)
    plt.close(fig)
    print(f"[seed figs] wrote {pdf}")
    return pdf


def _family_heatmap() -> Path:
    """miR-17~92 cluster coupling, 8 key Hallmarks x subtype (LumA/LumB deep-dive)."""
    mat = _key_matrix(FAM_SUB, family=CLUSTER)
    rows = [NAME.get(h, h) for h in mat.index]
    grid = mat.values.astype(float)
    norm = TwoSlopeNorm(vmin=-0.6, vcenter=0.0, vmax=0.05)
    cmap = plt.get_cmap("RdBu_r")

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    ax.imshow(grid, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(SUBTYPES)))
    ax.set_xticklabels(SUBTYPES, fontsize=8.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=8)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid[i, j]
            if pd.isna(v):
                continue
            txt = f"{v:+.2f}".replace("+0.", ".").replace("-0.", "\u2212.")
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.4,
                    color="white" if v < -0.30 else "black")
    ax.set_xticks(np.arange(-0.5, len(SUBTYPES), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.9)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(
        "miR-17~92 / 106b~25 cluster pressure \u2192 repression\n"
        "(partial \u03c1; \u2212 = repression; deepens LumA\u2192LumB\u2192Her2\u2192Basal)",
        fontsize=8.5,
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03, ticks=[-0.6, -0.4, -0.2, 0.0])
    cb.set_label("partial-Spearman \u03c1", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    pdf = OUT_DIR / "seed_family_subtype_heatmap.pdf"
    fig.savefig(pdf)
    fig.savefig(OUT_DIR / "seed_family_subtype_heatmap.png", dpi=600)
    plt.close(fig)
    print(f"[seed figs] wrote {pdf}")
    return pdf


def run() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _spider()
    _family_heatmap()


if __name__ == "__main__":
    run()
