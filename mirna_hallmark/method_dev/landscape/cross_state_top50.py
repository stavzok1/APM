"""Cross-state top-50 arm expression — tumour (by PAM50) / NAT / GTEx, stacked for comparison. 4 bands.

Band 1  TCGA primary tumour: the 50 most-abundant arms (cohort-median RPM), boxplot + a dot per tumour
        coloured by PAM50 subtype (LumA / LumB / Her2 / Basal; Normal-like excluded).
Band 2  TCGA tumour-adjacent NAT: the SAME 50 arms (tumour order), dots by matched-patient subtype. Labels
        carry **rank-Δ and log2 fold-change vs tumour** (NAT is the same RPM platform, so FC is meaningful).
Band 3  GTEx ∥ tumour-50: the SAME 50 tumour arms shown in healthy GTEx, **raw log2(TPM+1)** (own platform
        scale — compare by RANK). Labels carry rank-Δ vs tumour → which tumour-abundant arms are healthy-low
        (cancer-elevated, green) vs constitutively high.
Band 4  GTEx OWN top-50: GTEx's most-abundant arms, raw TPM. Labels green=more abundant in cancer / red=lowered.

GTEx scale = raw TPM (QN-onto-TCGA was rejected: it imposes TCGA's distribution *shape* on GTEx). The green/red
up-down calls and all Δ's use **platform-invariant ranks** (GTEx vs tumour over the 783 shared arms; NAT vs
tumour over the full 2,236-arm TCGA platform).

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.landscape.cross_state_top50`
Out: `figures/cross_state_top50.png`
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.data_loaders import load_clinical_strata

FIG = Path(__file__).parent.parent / "figures"
_SUB_C = {"LumA": "#2c7fb8", "LumB": "#41b6c4", "Her2": "#fe9929", "Basal": "#d7301f"}
_GRN, _RED = "#1b7837", "#b2182b"
_FLOOR = 10.0


def _strip_box(ax, arms, mat, colmap, alpha, ylab="RPM (log)", floor=True):
    rng = np.random.default_rng(0)
    data, pos = [], []
    for i, a in enumerate(arms):
        if a not in mat.index:                                   # arm absent in this state -> gap
            continue
        y = np.clip(mat.loc[a].to_numpy(), 0.1, None)
        ax.scatter(i + rng.uniform(-0.28, 0.28, len(y)), y, s=2, alpha=alpha,
                   c=[colmap(c) for c in mat.columns], edgecolors="none", zorder=1)
        data.append(y); pos.append(i)
    ax.boxplot(data, positions=pos, widths=0.62, showfliers=False, patch_artist=True, zorder=3,
               medianprops=dict(color="black", lw=1.2), boxprops=dict(facecolor="none", edgecolor="0.15", lw=0.8),
               whiskerprops=dict(color="0.3", lw=0.6), capprops=dict(color="0.3", lw=0.6))
    if floor:
        ax.axhline(_FLOOR, color="crimson", ls="--", lw=0.7, alpha=0.6)
    ax.set_yscale("log"); ax.set_xlim(-0.7, len(arms) - 0.3); ax.set_ylabel(ylab, fontsize=9)


def _color_ticks(ax, arms, val):
    for tick, a in zip(ax.get_xticklabels(), arms):
        v = val.get(a)
        tick.set_color("0.5" if v is None else (_GRN if v > 0 else _RED))


def build():
    st = _state_matrices()
    lin = lambda m: (2.0 ** m.astype(float).fillna(0.0) - 1.0).clip(lower=0.0)     # undetected (NaN) -> 0
    tl, nl, gl = lin(st["tumor"]), lin(st["nat"]), lin(st["gtex"])

    clin = load_clinical_strata().set_index("participant")["PAM50_final"]
    scol = lambda c: _SUB_C.get(clin.get(str(c)[:12], None), "0.6")

    t_med, n_med = tl.median(axis=1), nl.median(axis=1)
    top_t = t_med.sort_values(ascending=False).head(50).index.tolist()
    rt_full, rn_full = t_med.rank(ascending=False), n_med.rank(ascending=False)    # tumour & NAT: same 2,236-arm platform
    nat_dr = {a: int(rn_full[a] - rt_full[a]) for a in top_t}
    nat_fc = {a: float(np.log2(t_med[a] + 1) - np.log2(n_med[a] + 1)) for a in top_t}   # + = up in cancer (tumour>NAT)

    common = sorted(set(tl.index) & set(gl.index))                                 # GTEx ranks over shared arms
    rt_c = tl.loc[common].median(axis=1).rank(ascending=False)
    rg_c = gl.loc[common].median(axis=1).rank(ascending=False)
    gtc_dr = {a: (int(rg_c[a] - rt_c[a]) if a in rg_c.index else None) for a in top_t}   # + = cancer-elevated
    top_g = gl.loc[common].median(axis=1).sort_values(ascending=False).head(50).index.tolist()
    gown_dr = {a: int(rg_c[a] - rt_c[a]) for a in top_g}

    fig = plt.figure(figsize=(19, 24))
    gs = fig.add_gridspec(4, 1, hspace=0.85, top=0.955, bottom=0.04, left=0.05, right=0.99)
    axT, axN, axGc, axG = (fig.add_subplot(gs[i]) for i in range(4))

    # band 1 — tumour by subtype
    _strip_box(axT, top_t, tl, scol, alpha=0.14)
    axT.set_xticks(range(50)); axT.set_xticklabels([a.replace("hsa-", "") for a in top_t], rotation=90, fontsize=6.5)
    axT.set_title(f"Tumour — TCGA primary (n={tl.shape[1]:,}) · top-50 arms · dot per tumour coloured by PAM50 subtype",
                  fontsize=11, loc="left")
    axT.legend([plt.Line2D([], [], marker="o", color="w", markerfacecolor=_SUB_C[k], ms=7, label=k) for k in _SUB_C],
               list(_SUB_C), fontsize=8, ncol=4, loc="upper right", framealpha=0.95)

    # band 2 — NAT (same 50): rank-Δ + log2FC vs tumour
    _strip_box(axN, top_t, nl, scol, alpha=0.35)
    axN.set_xticks(range(50))
    axN.set_xticklabels([f"{a.replace('hsa-','')}  Δr{nat_dr[a]:+d} fc{nat_fc[a]:+.1f}" for a in top_t], rotation=90, fontsize=6)
    _color_ticks(axN, top_t, nat_fc)
    nun = sum(nat_fc[a] > 0 for a in top_t)
    axN.set_title(f"NAT — tumour-adjacent normal (n={nl.shape[1]}) · SAME 50 arms (tumour order) · dot by matched patient subtype · "
                  f"label = rank-Δ & log2FC vs tumour (green=up in cancer {nun} / red=down {50-nun})", fontsize=11, loc="left")

    # band 3 — GTEx for the SAME tumour-top-50 (raw TPM): rank-Δ vs tumour
    _strip_box(axGc, top_t, gl, lambda c: "#1b9e77", alpha=0.18, ylab="TPM (log)", floor=False)
    axGc.set_xticks(range(50))
    axGc.set_xticklabels([f"{a.replace('hsa-','')}  {('Δr%+d' % gtc_dr[a]) if gtc_dr[a] is not None else 'n/a'}" for a in top_t],
                         rotation=90, fontsize=6)
    _color_ticks(axGc, top_t, gtc_dr)
    nugc = sum(1 for a in top_t if gtc_dr[a] and gtc_dr[a] > 0)
    axGc.set_title(f"GTEx ∥ tumour's top-50 (n={gl.shape[1]}) · SAME 50 arms (tumour order) · raw TPM (own scale, compare by rank) · "
                   f"rank-Δ vs tumour (green=cancer-elevated/healthy-low {nugc} / red=healthy-higher)", fontsize=11, loc="left")

    # band 4 — GTEx OWN top-50 (raw TPM)
    _strip_box(axG, top_g, gl, lambda c: "#1b9e77", alpha=0.18, ylab="TPM (log)", floor=False)
    axG.set_xticks(range(50))
    axG.set_xticklabels([f"{a.replace('hsa-','')}  {'↑' if gown_dr[a] > 0 else '↓'}Δ{gown_dr[a]:+d}" for a in top_g],
                        rotation=90, fontsize=6)
    _color_ticks(axG, top_g, gown_dr)
    nug = sum(d > 0 for d in gown_dr.values())
    axG.set_title(f"Healthy GTEx breast (n={gl.shape[1]}) · GTEx's OWN top-50 · raw TPM (compare by rank) · "
                  f"label green=more abundant in cancer ({nug}) / red=lowered ({50-nug}); rank-Δ over 783 shared arms",
                  fontsize=11, loc="left")

    fig.suptitle("Top-50 miRNA-arm expression across states — tumour (by subtype) · NAT · GTEx∥tumour-50 · GTEx-own-50, TCGA-BRCA",
                 fontsize=14, fontweight="bold", y=0.978)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "cross_state_top50.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[cross-state] tumour {tl.shape[1]} / NAT {nl.shape[1]} / GTEx {gl.shape[1]}; "
          f"NAT up-in-cancer {nun}/50; GTEx∥tum cancer-elevated {nugc}/50; GTEx-own up {nug}/50; wrote figures/cross_state_top50.png")


if __name__ == "__main__":
    build()
