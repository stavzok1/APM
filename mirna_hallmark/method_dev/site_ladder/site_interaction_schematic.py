"""Conceptual schematic — the three site-pair regimes on a target 3'UTR, ordered by increasing
site overlap (Axis S). Two miRNA sites s_i, s_j with their loaded AGO/RISC complexes:

  S0  disjoint distal   (sep ≳40 nt)   independent       → log-additive repression
  S1  disjoint proximal (sep 8–40 nt)  cooperating AGO   → super-additive (synergy)
  S2  overlapping        (s_i ∩ s_j≠∅) mutually exclusive → site competition

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.site_interaction_schematic`
Out: `figures/site_interaction_regimes.png`
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, FancyArrowPatch, Rectangle

OUT = Path(__file__).parent.parent / "figures"     # shared figures dir at method_dev root
UTR0, UTR1 = 8, 86           # 3'UTR span (abstract nt)
C_I, C_J = "#E8743B", "#1E88A8"   # site i / site j colours
C_AGO = "#9aa6b2"


def _panel(ax, y, si, sj, sep_txt, regime, consequence, formula, overlap=False):
    # 3'UTR ribbon
    ax.add_patch(Rectangle((UTR0, y - 0.18), UTR1 - UTR0, 0.36, fc="#dfe7ee", ec="#b6c2cd", lw=1, zorder=1))
    ax.text(UTR0 - 1.5, y, "5′", ha="right", va="center", fontsize=8, color="0.4")
    ax.text(UTR1 + 1.5, y, "3′", ha="left", va="center", fontsize=8, color="0.4")
    ax.text((UTR0 + UTR1) / 2, y - 0.7, "target 3′UTR", ha="center", va="top", fontsize=7.5, color="0.5", style="italic")

    def _site(span, col, label):
        x0, x1 = span
        ax.add_patch(Rectangle((x0, y - 0.18), x1 - x0, 0.36, fc=col, ec="none", zorder=3))
        cx = (x0 + x1) / 2
        # miRNA guide (short bar) + AGO complex above the seed match
        ax.plot([cx, cx], [y + 0.18, y + 0.95], color=col, lw=1.2, zorder=2)
        ax.add_patch(Ellipse((cx, y + 1.45), 11, 1.5, fc=C_AGO, ec="0.55", lw=1, alpha=0.95, zorder=4))
        ax.add_patch(Ellipse((cx, y + 1.5), 4.2, 0.7, fc="white", ec="0.6", lw=0.7, zorder=5))
        ax.text(cx, y + 1.5, "AGO", ha="center", va="center", fontsize=6.5, color="0.35", zorder=6)
        ax.text(cx, y + 2.45, label, ha="center", va="center", fontsize=7.5, color=col, fontweight="bold")
        return cx

    ci = _site(si, C_I, r"$s_i$")
    cj = _site(sj, C_J, r"$s_j$")

    # separation / overlap annotation between the sites
    if overlap:
        ox0, ox1 = max(si[0], sj[0]), min(si[1], sj[1])
        ax.add_patch(Rectangle((ox0, y - 0.18), ox1 - ox0, 0.36, fc="#b3261e", ec="none", alpha=0.55, zorder=4))
        ax.annotate("", (ox0, y - 0.55), (ox1, y - 0.55), arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1))
        ax.text((ox0 + ox1) / 2, y - 0.92, sep_txt, ha="center", va="top", fontsize=7, color="#b3261e", fontweight="bold")
    else:
        ax.annotate("", (si[1], y + 0.0), (sj[0], y + 0.0), arrowprops=dict(arrowstyle="<->", color="0.45", lw=1))
        ax.text((si[1] + sj[0]) / 2, y + 0.30, sep_txt, ha="center", va="bottom", fontsize=7, color="0.4")

    # left regime tag + right consequence
    ax.text(-2.5, y + 0.6, regime, ha="left", va="center", fontsize=15, fontweight="bold", color="0.2")
    ax.text(93, y + 0.95, consequence, ha="left", va="center", fontsize=8.5, fontweight="bold", color="0.15")
    ax.text(93, y + 0.05, formula, ha="left", va="center", fontsize=8, color="0.3")


def build() -> None:
    fig, ax = plt.subplots(figsize=(13.5, 8.2))
    _panel(ax, 11, (20, 27), (70, 77), r"sep $\gtrsim 40$ nt", "S0",
           "disjoint · distal → independent",
           r"$F = F_i \cdot F_j$   (log-additive repression)")
    _panel(ax, 6.5, (39, 46), (57, 64), r"sep $\sim 8$–$40$ nt", "S1",
           "disjoint · proximal → cooperating AGO",
           r"$F < F_i \cdot F_j$   (super-additive · synergy)")
    _panel(ax, 2, (43, 50), (47, 54), r"$s_i \cap s_j \neq \varnothing$", "S2",
           "overlapping → mutually exclusive",
           r"$F \approx \min(F_i, F_j)$   (site competition; ≤1 bound)", overlap=True)

    # Axis S (increasing overlap) on the far left
    ax.annotate("", (-4.2, 12.5), (-4.2, 1.0), arrowprops=dict(arrowstyle="-|>", color="0.45", lw=1.6))
    ax.text(-5.2, 6.7, "Axis S  —  increasing site overlap", rotation=90, ha="center", va="center",
            fontsize=10, color="0.4", fontweight="bold")

    ax.text(40, 14.0, "Two miRNA sites on one 3′UTR: regimes by site overlap",
            ha="center", fontsize=13, fontweight="bold")
    ax.text(93, 13.2, "$F$ = fraction mRNA remaining\n(lower = stronger repression)", fontsize=7.5, color="0.4", va="top")
    ax.set_xlim(-7, 116); ax.set_ylim(-0.5, 14.6); ax.axis("off")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "site_interaction_regimes.png", dpi=160)
    print(f"[schematic] wrote {OUT/'site_interaction_regimes.png'}")


if __name__ == "__main__":
    build()
