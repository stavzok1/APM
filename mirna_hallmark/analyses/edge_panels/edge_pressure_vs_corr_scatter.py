"""Scatter: edge pressure metrics vs partial-ρ with target expression.

Two X-axis views, one Y-axis, enriched with n_regulators tier and
evidence-flat share for comparison.

  Y  — partial Spearman ρ: per-sample edge pressure vs target gene expression,
       adjusted for CPE + HRD + E2F/G2M proliferation proxy.

  X₁ — **gene-centric share** (0–1, evidence-weighted):
           share(m, g) = mean_s[ c(m,g,s) / Σ_{m'} c(m',g,s) ]
       The pressure c() already encodes log1p(evidence_score) in both the
       edge weight and the normalization denominator — so share reflects both
       expression and literature support.

  X₂ — **gene-centric share, evidence-flat** (0–1):
       Same formula but all evidence_score set to 1.0, so the share is driven
       purely by relative expression across competing miRNA arms — no database
       weighting.  Difference between X₁ and X₂ isolates the literature-weight
       contribution.

  X₃ — **mean raw pressure** (log₂):
       mean_s[ c(m,g,s) ] — absolute un-normalised pressure, log₂-scaled.

Additional stratification:
  n_regulators = distinct miRNA arms targeting this gene in the high-evidence
  edge set.  Tier: 1 (singleton) | 2–4 | 5+.  Singletons trivially have
  share = 1.0 — not informative dominance.

Run:
  .venv/bin/python3 -m mirna_hallmark.edge_pressure_vs_corr_scatter
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.analyses.edge_panels.edge_pressure_panels import TARGET_NORM, ATTR_MODE, _load_edges
from mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels import _state_pressure_map

OUT_DIR = C.TISSUE_REFERENCE_DIR / "edge_pressure_vs_corr_scatter"
CORR_TABLE = (
    C.TISSUE_REFERENCE_DIR
    / "edge_partial_corr_panels"
    / "edge_partial_corr_panel_corr_table.tsv"
)


# --------------------------------------------------------------------------- #
# Share / pressure computation
# --------------------------------------------------------------------------- #

def _gene_centric_shares(
    tumor_map: Dict[tuple, pd.Series],
) -> Dict[tuple, float]:
    """Per-edge mean gene-centric share across tumor samples.

    share(m, g, s) = c(m, g, s) / Σ_{m'} c(m', g, s)
    Returns mean over samples.  Bounded in (0, 1).
    """
    gene_to_keys: dict = defaultdict(list)
    for (arm, gene) in tumor_map:
        gene_to_keys[gene].append((arm, gene))

    out: Dict[tuple, float] = {}
    for gene, keys in gene_to_keys.items():
        # align Series on common sample index
        series_list = [
            pd.to_numeric(tumor_map[k], errors="coerce").rename(k)
            for k in keys
        ]
        mat = pd.concat(series_list, axis=1)   # samples × arms
        total = mat.sum(axis=1).replace(0, np.nan)
        for k in keys:
            share_s = mat[k] / total
            out[k] = float(share_s.mean(skipna=True))
    return out


N_REG_TIERS = {
    "singleton (n=1)": (1, 1),
    "2–4 regulators":  (2, 4),
    "5+ regulators":   (5, 9999),
}
N_REG_TIER_COLORS = {
    "singleton (n=1)": "#E15759",
    "2–4 regulators":  "#4E79A7",
    "5+ regulators":   "#59A14F",
}


def _assign_tier(n: int) -> str:
    for label, (lo, hi) in N_REG_TIERS.items():
        if lo <= n <= hi:
            return label
    return "5+ regulators"


def build_scatter_table(
    *,
    high_evidence_only: bool = True,
    hallmark: Optional[str] = None,
    corr_table_path: Path = CORR_TABLE,
) -> pd.DataFrame:
    """One row per edge with evidence-weighted share, evidence-flat share,
    raw pressure, partial-ρ, and n_regulators tier."""
    edges = _load_edges(high_evidence_only=high_evidence_only, hallmark=hallmark)
    genes = sorted(edges["gene"].unique())
    print(f"[scatter] edges={len(edges)}  genes={len(genes)}")

    # n_regulators per gene
    n_reg = edges.groupby("gene")["miRNA"].nunique().rename("n_regulators")

    tissue = _state_matrices()
    tumor_mat = tissue["tumor"]
    print(f"[scatter] computing evidence-weighted pressure map (n_tumor={tumor_mat.shape[1]}) …")
    tumor_map = _state_pressure_map(edges, tumor_mat, genes)

    # evidence-flat pressure map: all scores = 1 → log1p(1) cancels in edge_w / norm
    edges_flat = edges.copy()
    edges_flat["evidence_score"] = 1.0
    print("[scatter] computing evidence-flat pressure map …")
    tumor_map_flat = _state_pressure_map(edges_flat, tumor_mat, genes)

    print("[scatter] computing gene-centric shares (weighted + flat) …")
    gc_shares      = _gene_centric_shares(tumor_map)
    gc_shares_flat = _gene_centric_shares(tumor_map_flat)

    rows = []
    for (arm, gene), s in tumor_map.items():
        vals = pd.to_numeric(s, errors="coerce").dropna()
        if len(vals) == 0:
            continue
        n_r = int(n_reg.get(gene, 1))
        rows.append({
            "miRNA": arm, "gene": gene,
            "n_regulators": n_r,
            "n_reg_tier": _assign_tier(n_r),
            "gene_centric_share": gc_shares.get((arm, gene), np.nan),
            "gene_centric_share_flat": gc_shares_flat.get((arm, gene), np.nan),
            "mean_pressure": float(vals.mean()),
            "n_tumor": len(vals),
        })
    pressure_df = pd.DataFrame(rows)
    print(f"[scatter] metrics for {len(pressure_df)} edges")

    corr_df = pd.read_csv(corr_table_path, sep="\t")
    corr_df = corr_df.dropna(subset=["rho_adj"])
    print(f"[scatter] partial-ρ available for {len(corr_df)} edges")

    merged = pressure_df.merge(
        corr_df[["miRNA", "gene", "evidence_score",
                 "rho_adj", "p_adj", "q_adj", "n_adj"]],
        on=["miRNA", "gene"], how="inner",
    )
    merged["log2_mean_pressure"] = np.log2(merged["mean_pressure"].clip(lower=1e-12))
    merged["sig"] = merged["q_adj"].lt(0.05)
    merged["label"] = merged["miRNA"].str.replace("hsa-", "") + "→" + merged["gene"]

    print(f"[scatter] final table: {len(merged)} edges")
    for tier in N_REG_TIERS:
        sub = merged[merged["n_reg_tier"].eq(tier)]
        if len(sub):
            print(f"  {tier}: n={len(sub)}"
                  f"  share=[{sub['gene_centric_share'].min():.3f},"
                  f"{sub['gene_centric_share'].max():.3f}]")
    return merged


# --------------------------------------------------------------------------- #
# Shared plot helpers
# --------------------------------------------------------------------------- #

def _sig_dir_colors(df: pd.DataFrame) -> list:
    """Color by significance + direction (used when n_reg tier is not the color)."""
    def _c(row):
        if not row["sig"]:
            return "#aaaaaa"
        return "#4E79A7" if row["rho_adj"] < 0 else "#E15759"
    return list(df.apply(_c, axis=1))


def _tier_colors(df: pd.DataFrame) -> list:
    return [N_REG_TIER_COLORS[t] for t in df["n_reg_tier"]]


def _scatter_sizes(df: pd.DataFrame, base: float = 8, scale: float = 40) -> pd.Series:
    ev = df["evidence_score"].fillna(1).clip(lower=1)
    return base + (ev - ev.min()) / (ev.max() - ev.min() + 1e-9) * scale


def _sig_legend(df: pd.DataFrame):
    from matplotlib.patches import Patch
    n_rep = int((df["sig"] & df["rho_adj"].lt(0)).sum())
    n_pos = int((df["sig"] & df["rho_adj"].gt(0)).sum())
    n_ns  = int((~df["sig"]).sum())
    return [
        Patch(facecolor="#4E79A7", alpha=0.7,
              label=f"q<0.05, ρ<0  repressive (n={n_rep})"),
        Patch(facecolor="#E15759", alpha=0.7,
              label=f"q<0.05, ρ>0  co-activated (n={n_pos})"),
        Patch(facecolor="#aaaaaa", alpha=0.7,
              label=f"q≥0.05  (n={n_ns})"),
    ]


def _tier_legend():
    from matplotlib.patches import Patch
    return [
        Patch(facecolor=N_REG_TIER_COLORS[t], alpha=0.7, label=t)
        for t in N_REG_TIERS
    ]


def _annotate_top(ax, df: pd.DataFrame, x_col: str, n: int = 20) -> None:
    score = df["sig"].astype(int) * df["rho_adj"].abs()
    top = df.assign(_score=score).nlargest(n, "_score")
    texts = []
    for _, row in top.iterrows():
        texts.append(ax.text(
            row[x_col], row["rho_adj"], row["label"],
            fontsize=5.5, color="#222222",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.5, ec="none"),
        ))
    try:
        from adjustText import adjust_text  # type: ignore
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    expand_points=(1.2, 1.4), force_text=(0.2, 0.4))
    except Exception:
        pass


def _base_scatter(ax, df, x_col, colors, sizes):
    ax.scatter(df[x_col], df["rho_adj"],
               c=colors, s=sizes, alpha=0.5, linewidths=0, zorder=2)
    ax.axhline(0, color="#555555", lw=0.8, linestyle="--", zorder=1)
    ax.axvline(float(df[x_col].median()), color="#888888", lw=0.5,
               linestyle=":", zorder=1)
    ax.grid(alpha=0.16, linestyle="--", zorder=0)


# --------------------------------------------------------------------------- #
# Figure 1: share vs ρ — 3 facets by n_regulators tier
# --------------------------------------------------------------------------- #

def plot_share_by_tier(
    df: pd.DataFrame,
    out_path: Path,
    *,
    x_col: str = "gene_centric_share",
    x_label: str = "Gene-centric pressure share (evidence-weighted, 0–1)",
    title_suffix: str = "",
    dpi: int = 150,
    n_label: int = 15,
) -> None:
    """3-panel figure: one facet per n_regulators tier, sig-direction coloring."""
    import matplotlib.pyplot as plt

    tiers = list(N_REG_TIERS.keys())
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.5), sharey=True)

    for ax, tier in zip(axes, tiers):
        sub = df[df["n_reg_tier"].eq(tier)]
        if sub.empty:
            ax.axis("off")
            continue
        colors = _sig_dir_colors(sub)
        sizes  = _scatter_sizes(sub)
        _base_scatter(ax, sub, x_col, colors, sizes)
        _annotate_top(ax, sub, x_col, n=n_label)

        n_sig = int(sub["sig"].sum())
        ax.set_title(
            f"{tier}  (n={len(sub)}, q<0.05={n_sig})\n"
            f"share range [{sub[x_col].min():.2f}, {sub[x_col].max():.2f}]",
            fontsize=9,
        )
        ax.set_xlim(-0.02, 1.02)
        ax.set_xlabel(x_label, fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel("Partial Spearman ρ", fontsize=9)

    for ax in axes:
        ax.legend(handles=_sig_legend(df), fontsize=7, frameon=False,
                  loc="upper right")

    n_sig_total = int(df["sig"].sum())
    fig.suptitle(
        f"Gene-centric share vs partial ρ — stratified by n_regulators{title_suffix}\n"
        f"n={len(df)} edges total  ·  q<0.05={n_sig_total}  ·  "
        "softmax_logrpm  ·  adj. CPE+HRD+E2F/G2M",
        fontsize=10,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[scatter] tier facets → {out_path}")


# --------------------------------------------------------------------------- #
# Figure 2: evidence-weighted vs evidence-flat share — scatter coloured by tier
# --------------------------------------------------------------------------- #

def plot_weighted_vs_flat(
    df: pd.DataFrame,
    out_path: Path,
    *,
    dpi: int = 150,
    n_label: int = 20,
) -> None:
    """Evidence-weighted share (X) vs evidence-flat share (Y).
    Diagonal = no effect of evidence weighting.
    Points above diagonal = evidence weight inflates share.
    Colour = n_reg_tier; opacity encodes |ρ_adj|.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 8))

    colors = _tier_colors(df)
    rho_abs = df["rho_adj"].abs().fillna(0)
    alphas = 0.25 + 0.55 * (rho_abs / rho_abs.max())
    sizes  = _scatter_sizes(df)

    for idx, (_, row) in enumerate(df.iterrows()):
        ax.scatter(row["gene_centric_share"], row["gene_centric_share_flat"],
                   c=colors[idx], s=sizes.iloc[idx],
                   alpha=float(alphas.iloc[idx]), linewidths=0, zorder=2)

    # identity line
    ax.plot([0, 1], [0, 1], color="#888888", lw=0.9, linestyle="--",
            label="no evidence-weight effect", zorder=1)

    # annotate top edges where weighted ≠ flat by |ρ|
    diff_score = (df["gene_centric_share"] - df["gene_centric_share_flat"]).abs() * rho_abs
    top = df.assign(_ds=diff_score).nlargest(n_label, "_ds")
    texts = []
    for _, row in top.iterrows():
        texts.append(ax.text(
            row["gene_centric_share"], row["gene_centric_share_flat"],
            row["label"], fontsize=5.5, color="#333333",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.5, ec="none"),
        ))
    try:
        from adjustText import adjust_text  # type: ignore
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.4))
    except Exception:
        pass

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(
        "Evidence-weighted gene-centric share\n"
        "[log1p(evidence_score) in edge weight + normalization]",
        fontsize=10,
    )
    ax.set_ylabel(
        "Evidence-flat gene-centric share\n"
        "[uniform evidence weight — pure expression-driven]",
        fontsize=10,
    )
    ax.set_title(
        "Where does evidence weighting move share relative to expression-only?\n"
        "Above diagonal: evidence inflates share  ·  "
        "Below diagonal: evidence suppresses share  ·  "
        "opacity ∝ |ρ_adj|",
        fontsize=9.5,
    )

    # tier + rho legend
    handles = _tier_legend()
    from matplotlib.lines import Line2D
    handles += [
        Line2D([0], [0], color="#888888", lw=0.9, linestyle="--",
               label="weighted = flat"),
    ]
    ax.legend(handles=handles, fontsize=8, frameon=False)
    ax.grid(alpha=0.15, linestyle="--", zorder=0)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[scatter] weighted vs flat → {out_path}")


# --------------------------------------------------------------------------- #
# Figure 3: raw pressure vs ρ, coloured by tier
# --------------------------------------------------------------------------- #

def plot_pressure_scatter(
    df: pd.DataFrame,
    out_path: Path,
    *,
    dpi: int = 150,
    n_label: int = 20,
) -> None:
    import matplotlib.pyplot as plt

    x_col = "log2_mean_pressure"
    fig, ax = plt.subplots(figsize=(10, 7.5))

    colors = _tier_colors(df)
    sizes  = _scatter_sizes(df)
    _base_scatter(ax, df, x_col, colors, sizes)
    _annotate_top(ax, df, x_col, n=n_label)

    ax.set_xlabel("Mean tumor edge pressure  [log₂ c(m,g,s)]", fontsize=10)
    ax.set_ylabel("Partial Spearman ρ  (edge pressure vs target expression)", fontsize=10)
    n_sig = int(df["sig"].sum())
    ax.set_title(
        "Raw edge pressure vs partial ρ  (coloured by n_regulators tier)\n"
        f"n={len(df)}  ·  q<0.05={n_sig}  ·  softmax_logrpm  ·  adj. CPE+HRD+E2F/G2M",
        fontsize=9.5,
    )

    handles = _tier_legend()
    ax.legend(handles=handles, fontsize=8, frameon=False)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[scatter] pressure scatter → {out_path}")


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(
    *,
    out_dir: Path = OUT_DIR,
    high_evidence_only: bool = True,
    hallmark: Optional[str] = None,
    dpi: int = 150,
    n_label: int = 25,
    corr_table_path: Path = CORR_TABLE,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    df = build_scatter_table(
        high_evidence_only=high_evidence_only,
        hallmark=hallmark,
        corr_table_path=corr_table_path,
    )
    df.to_csv(out_dir / "scatter_table.tsv", sep="\t", index=False)

    fig_dir = out_dir / "figures"

    # 1. share vs ρ — 3 facets by n_reg tier (evidence-weighted)
    plot_share_by_tier(
        df, fig_dir / "share_vs_rho_by_nreg_tier.png",
        dpi=dpi, n_label=n_label,
    )
    # 2. same facets but with evidence-flat share
    plot_share_by_tier(
        df, fig_dir / "share_flat_vs_rho_by_nreg_tier.png",
        x_col="gene_centric_share_flat",
        x_label="Gene-centric pressure share (evidence-FLAT, 0–1)",
        title_suffix=" [flat evidence]",
        dpi=dpi, n_label=n_label,
    )
    # 3. evidence-weighted vs flat share — where does evidence weighting move points?
    plot_weighted_vs_flat(df, fig_dir / "weighted_vs_flat_share.png",
                          dpi=dpi, n_label=n_label)
    # 4. raw pressure vs ρ coloured by tier
    plot_pressure_scatter(df, fig_dir / "pressure_vs_rho_by_nreg_tier.png",
                          dpi=dpi, n_label=n_label)

    n_per_tier = df.groupby("n_reg_tier").size().to_dict()
    manifest = {
        "module": "mirna_hallmark.edge_pressure_vs_corr_scatter",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_edges": len(df),
        "n_sig_q05": int(df["sig"].sum()),
        "n_repressive_sig": int((df["sig"] & df["rho_adj"].lt(0)).sum()),
        "n_positive_sig": int((df["sig"] & df["rho_adj"].gt(0)).sum()),
        "n_per_tier": n_per_tier,
        "figures": {
            "share_by_tier":       "figures/share_vs_rho_by_nreg_tier.png",
            "share_flat_by_tier":  "figures/share_flat_vs_rho_by_nreg_tier.png",
            "weighted_vs_flat":    "figures/weighted_vs_flat_share.png",
            "pressure_by_tier":    "figures/pressure_vs_rho_by_nreg_tier.png",
        },
        "x_cols": {
            "gene_centric_share":      "evidence-weighted gene-centric share ∈ (0,1)",
            "gene_centric_share_flat": "evidence-flat gene-centric share ∈ (0,1)",
            "mean_pressure":            "mean tumor pressure (linear)",
        },
        "y_axis": "partial Spearman rho adj. CPE+HRD+E2F/G2M",
        "n_reg_tier_def": {t: f"{lo}–{hi}" for t, (lo, hi) in N_REG_TIERS.items()},
        "attr_mode": ATTR_MODE,
        "target_norm": TARGET_NORM,
        "high_evidence_only": high_evidence_only,
        "hallmark_filter": hallmark,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[scatter] done  n={len(df)}  q<0.05={manifest['n_sig_q05']}")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--high-evidence-only", action="store_true", default=True)
    ap.add_argument("--no-high-evidence-filter", dest="high_evidence_only",
                    action="store_false")
    ap.add_argument("--hallmark", type=str, default=None)
    ap.add_argument("--n-label", type=int, default=25)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--corr-table", type=Path, default=CORR_TABLE)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(
        out_dir=args.out_dir,
        high_evidence_only=args.high_evidence_only,
        hallmark=args.hallmark,
        dpi=args.dpi,
        n_label=args.n_label,
        corr_table_path=args.corr_table,
    )


if __name__ == "__main__":
    main()
