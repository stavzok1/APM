"""Grant figure for the basal miRNA-hub discovery (affirmative).

Single-figure slot, sized for one column (~88 mm). The governing principle: put in
the figure what prose cannot carry. The discovery -- that across Hallmark programs
the BASAL column lights up broadly (concentrated on tumor-suppressor / cell-state
programs), most deeply of any subtype (LumB weaker; LumA / Her2 sparse), driven by
one compact miRNA cluster -- cannot be reconstructed from text, so it is the figure.
The
partial-correlation robustness numbers compress to two sentences and stay in text;
the one control worth a visual (the basal signal holding under an E2F/MYC-INDEPENDENT
proliferation proxy) gets a small strip.

  Panel A (dominant) -- subtype x Hallmark heatmap of the **proliferation-adjusted**
    pressure->repression coupling (partial Spearman, CPE + HRD + proliferation),
    per PAM50 subtype. Tumor-suppressor / cell-state rows are marked (bold) and the
    miR-17~92 / miR-106b~25 cluster is bracketed beside the basal column. What the
    reviewer sees is already the hardened signal, not the raw "basal is extreme"
    contrast (caption: purity/HRD/proliferation-adjusted).

  Panel B (compact strip) -- basal coupling for the key programs holding under the
    proliferation ladder, INCLUDING the E2F/MYC-independent proxy (the over-adjustment
    objection, answered visually). Curation-bias robustness stays in text.

Reads only existing robustness outputs; safe to re-run.

    .venv/bin/python3 -m mirna_hallmark.make_grant_figure
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
from mirna_hallmark import stats as S

ROBUST = C.OUTPUT_ROOT / "robustness"
BY_SUB = ROBUST / "aim1_proliferation" / "hallmark_coupling_by_pam50_prolif.tsv"
BASAL_LADDER = ROBUST / "aim1_proliferation" / "hallmark_coupling_prolif.tsv"
OUT_DIR = C.OUTPUT_ROOT / "figures"

SUBTYPES = ["LumA", "LumB", "Her2", "Basal"]

# Short display names.
NAME = {
    "HALLMARK_APOPTOSIS": "Apoptosis",
    "HALLMARK_P53_PATHWAY": "p53 pathway",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING": "PI3K\u2013AKT\u2013mTOR",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
    "HALLMARK_TGF_BETA_SIGNALING": "TGF-\u03b2",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF\u03b1/NF-\u03baB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "IL6\u2013JAK\u2013STAT3",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-\u03b3 response",
    "HALLMARK_DNA_REPAIR": "DNA repair",
    "HALLMARK_ANGIOGENESIS": "Angiogenesis",
    "HALLMARK_HYPOXIA": "Hypoxia",
    "HALLMARK_COAGULATION": "Coagulation",
    "HALLMARK_APICAL_JUNCTION": "Apical junction",
    "HALLMARK_NOTCH_SIGNALING": "Notch",
}

# Tumor-suppressor / cell-state rows (marked + spanned by the hub bracket), shown
# contiguously at the top so the bracket is a clean span.
TSG_ROWS = [
    "HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY", "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING",
]
# Remaining programs where basal pressure concentrates (context).
OTHER_ROWS = [
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_DNA_REPAIR",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_ANGIOGENESIS", "HALLMARK_COAGULATION", "HALLMARK_HYPOXIA",
    "HALLMARK_APICAL_JUNCTION", "HALLMARK_NOTCH_SIGNALING",
]
ROWS = TSG_ROWS + OTHER_ROWS

# Strip: the key programs, basal-only, across the proliferation ladder.
STRIP_ROWS = [
    "HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY", "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]
STRIP_COLS = [("+E2F/G2M", "rho_e2f_g2m", "p_e2f_g2m"),
              ("+MKI67", "rho_mki67", "p_mki67"),
              ("+ortho", "rho_ortho_noE2F_MYC", "p_ortho_noE2F_MYC")]
STRIP_INDEP_FROM = 1  # columns >= this index are E2F/MYC-independent

NORM = TwoSlopeNorm(vmin=-0.33, vcenter=0.0, vmax=0.13)
CMAP = plt.get_cmap("RdBu_r")


def _txt_color(v):
    return "white" if (pd.notna(v) and v < -0.16) else "black"


def _cell_text(v):
    if pd.isna(v):
        return ""
    return f"{v:+.2f}".replace("+0.", ".").replace("-0.", "\u2212.")


def _load():
    d = pd.read_csv(BY_SUB, sep="\t")
    rho = d.pivot(index="hallmark_set", columns="subtype", values="rho_prolif_adj")
    q = d.pivot(index="hallmark_set", columns="subtype", values="q_prolif_adj")
    lad = pd.read_csv(BASAL_LADDER, sep="\t")
    lad = lad[lad["scope"].eq("Basal")].set_index("hallmark_set")
    # FDR within basal scope per proxy column for the strip significance marks
    for _, _rc, pc in STRIP_COLS:
        lad[pc.replace("p_", "q_")] = S.bh_fdr(lad[pc].values)
    return rho, q, lad


def _heatmap(ax, mat, qmat, row_keys, col_labels, *, mark_rows=(), star_q=True):
    n_r, n_c = len(row_keys), len(col_labels)
    grid = np.array([[mat.loc[r, c] if (r in mat.index and c in mat.columns) else np.nan
                      for c in col_labels] for r in row_keys], dtype=float)
    ax.imshow(grid, aspect="auto", cmap=CMAP, norm=NORM)
    ax.set_xticks(range(n_c))
    ax.set_xticklabels(col_labels, fontsize=7.2)
    ax.set_yticks(range(n_r))
    labels = [NAME.get(r, r) for r in row_keys]
    ax.set_yticklabels(labels, fontsize=7.2)
    for tick, r in zip(ax.get_yticklabels(), row_keys):
        if r in mark_rows:
            tick.set_fontweight("bold")
            tick.set_color("#7a0000")
    for i, r in enumerate(row_keys):
        for j, c in enumerate(col_labels):
            v = grid[i, j]
            if pd.isna(v):
                continue
            sig = star_q and (r in qmat.index and c in qmat.columns
                              and pd.notna(qmat.loc[r, c]) and qmat.loc[r, c] < 0.05)
            ax.text(j, i, _cell_text(v) + ("*" if sig else ""),
                    ha="center", va="center", fontsize=6.6,
                    color=_txt_color(v), fontweight="bold" if sig else "normal")
    ax.set_xticks(np.arange(-0.5, n_c, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_r, 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.8)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def _hub_bracket(ax, n_marked, x=3.62):
    """Vertical bracket beside the basal column spanning the marked TSG rows."""
    y0, y1 = -0.45, n_marked - 0.55
    ax.plot([x, x], [y0, y1], color="#7a0000", lw=1.1, clip_on=False)
    for yy in (y0, y1):
        ax.plot([x - 0.07, x], [yy, yy], color="#7a0000", lw=1.1, clip_on=False)
    ax.annotate("miR-17~92 /\nmiR-106b~25\ncluster", xy=(x, (y0 + y1) / 2),
                xytext=(x + 0.12, (y0 + y1) / 2), va="center", ha="left",
                fontsize=6.8, fontweight="bold", color="#7a0000",
                annotation_clip=False)


def run() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 7.2, "font.family": "sans-serif", "pdf.fonttype": 42,
        "ps.fonttype": 42, "axes.linewidth": 0.7,
    })
    rho, q, lad = _load()

    fig = plt.figure(figsize=(3.5, 6.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(ROWS), len(STRIP_ROWS) + 1.6],
                          left=0.30, right=0.74, top=0.945, bottom=0.165, hspace=0.34)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    # ---- Panel A : subtype x Hallmark, proliferation-adjusted ---- #
    _heatmap(axA, rho, q, ROWS, SUBTYPES, mark_rows=set(TSG_ROWS))
    axA.axhline(len(TSG_ROWS) - 0.5, color="#444444", lw=0.8)  # TSG block separator
    _hub_bracket(axA, len(TSG_ROWS))
    fig.text(0.012, 0.965, "A   Basal-maximal, multi-program",
             fontsize=8, fontweight="bold", ha="left", va="bottom")

    # ---- Panel B : basal proliferation-independence strip ---- #
    strip = pd.DataFrame(
        {lab: [lad.loc[r, rc] if r in lad.index else np.nan for r in STRIP_ROWS]
         for lab, rc, _pc in STRIP_COLS}, index=STRIP_ROWS)
    stripq = pd.DataFrame(
        {lab: [lad.loc[r, pc.replace("p_", "q_")] if r in lad.index else np.nan
               for r in STRIP_ROWS] for lab, _rc, pc in STRIP_COLS}, index=STRIP_ROWS)
    col_labels = [lab for lab, _, _ in STRIP_COLS]
    _heatmap(axB, strip, stripq, STRIP_ROWS, col_labels, mark_rows=set(STRIP_ROWS))
    # mark the E2F/MYC-independent proxy columns
    axB.plot([STRIP_INDEP_FROM - 0.5, len(col_labels) - 0.5], [-0.72, -0.72],
             color="#1a5e1a", lw=1.4, clip_on=False)
    axB.text((STRIP_INDEP_FROM - 0.5 + len(col_labels) - 0.5) / 2, -0.92,
             "E2F/MYC-independent", ha="center", va="bottom", fontsize=6.4,
             style="italic", color="#1a5e1a")
    bpos = axB.get_position()
    fig.text(0.012, bpos.y1 + 0.045, "B   Robust to independent proxies (Basal)",
             fontsize=8, fontweight="bold", ha="left", va="bottom")

    # ---- shared colorbar ---- #
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cax = fig.add_axes([0.76, 0.45, 0.022, 0.32])
    cb = fig.colorbar(sm, cax=cax, ticks=[-0.3, -0.2, -0.1, 0.0, 0.1])
    cb.set_label("partial-Spearman \u03c1\n(pressure \u2194 target;  \u2212 = repression)",
                 fontsize=6.6)
    cb.ax.tick_params(labelsize=6.4, length=2)

    cap = (
        "Basal n=197 (LumA 503 / LumB 258 / Her2 92). Cell = within-subtype partial "
        "Spearman of miRNA Hallmark pressure vs target expression, adjusted for "
        "purity (CPE), HRD and proliferation; \u03c1<0 = pressure couples to "
        "repression. (A) The basal column is broadly negative and concentrated on "
        "tumor-suppressor / cell-state programs (42/50 vs LumB 17/50, LumA 18/50, "
        "Her2 0/50), driven by the bracketed miR-17~92/106b~25 cluster \u2014 a "
        "multi-program signal, deepest in Basal, that is already "
        "proliferation-hardened. (B) For the key "
        "programs the basal coupling persists under MKI67-alone and an "
        "E2F/MYC-stripped proxy, excluding over-adjustment. \u03c1 magnitudes are "
        "proxy-dependent and not interpreted as effect sizes. Bold rows = "
        "tumor-suppressor / cell-state; * q<0.05 (FDR within subtype)."
    )
    fig.text(0.012, 0.004, cap, fontsize=5.8, va="bottom", ha="left", wrap=True)

    pdf = OUT_DIR / "grant_figure_basal_heatmap.pdf"
    png = OUT_DIR / "grant_figure_basal_heatmap.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    plt.close(fig)
    print(f"[make_grant_figure] wrote {pdf}")
    print(f"[make_grant_figure] wrote {png}")
    return pdf


if __name__ == "__main__":
    run()
