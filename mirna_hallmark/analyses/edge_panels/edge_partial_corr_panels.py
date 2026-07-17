"""Per-edge scatter panels: tumor miRNA pressure vs target gene expression.

Each subplot shows one miRNA→gene pair in the primary-tumor cohort:

  - X: per-sample edge pressure c(m,g,s)  (softmax_logrpm / evidence_mass)
  - Y: target gene log2(TPM+1)
  - Dots colored by PAM50 subtype
  - Regression line across all samples
  - Annotated with partial Spearman ρ adjusted for CPE + HRD + E2F/G2M proliferation
    proxy — the most comprehensive confounder set used in this project's
    robustness spine (``robustness_checks.py``, Aim 1).

Negative ρ means pressure and target expression anti-correlate (repression
signal).  Edges are ranked by |ρ_adj| descending (strongest coupling first).

Run:
  .venv/bin/python3 -m mirna_hallmark.edge_partial_corr_panels
  .venv/bin/python3 -m mirna_hallmark.edge_partial_corr_panels --max-edges 300
  .venv/bin/python3 -m mirna_hallmark.edge_partial_corr_panels --hallmark HALLMARK_APOPTOSIS
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cross_state.cross_state_expression_panels import PAM50_COLORS, PAM50_ORDER
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.analyses.edge_panels.edge_pressure_panels import (
    OUT_DIR,
    TARGET_NORM,
    ATTR_MODE,
    _load_edges,
    _paginate,
    build_panel_context,
)
from mirna_hallmark.family_normal_reference import _participant
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_engine import compute_edge_pressure_map
from mirna_hallmark.stats import bh_fdr

OUT_DIR_CORR = C.TISSUE_REFERENCE_DIR / "edge_partial_corr_panels"


# --------------------------------------------------------------------------- #
# Covariate builders
# --------------------------------------------------------------------------- #

def _prolif_proxy(rna: pd.DataFrame, hs: HallmarkSets) -> pd.Series:
    """Mean z-score of HALLMARK_E2F_TARGETS ∪ HALLMARK_G2M_CHECKPOINT per sample."""
    panel = sorted(
        (set(hs.sets.get("HALLMARK_E2F_TARGETS", [])) | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", [])))
        & set(rna.index)
    )
    if not panel:
        return pd.Series(dtype=float, name="prolif_e2f_g2m")
    sub = rna.loc[panel].astype(float)
    mu = sub.mean(axis=1)
    sigma = sub.std(axis=1).replace(0, np.nan)
    z = sub.sub(mu, axis=0).div(sigma, axis=0)
    return z.mean(axis=0).rename("prolif_e2f_g2m")


def _build_cov(
    clin: pd.DataFrame,
    participants: pd.Index,
    prolif: pd.Series,
) -> pd.DataFrame:
    """CPE + HRD + prolif_e2f_g2m, aligned to 12-char *participants*."""
    clin_idx = clin.set_index("participant")
    return pd.concat([
        pd.to_numeric(clin_idx["CPE"].reindex(participants), errors="coerce").rename("CPE"),
        pd.to_numeric(clin_idx["thornsson_hrd_score"].reindex(participants), errors="coerce").rename("HRD"),
        prolif.reindex(participants).rename("prolif"),
    ], axis=1)


# --------------------------------------------------------------------------- #
# Per-edge partial correlation
# --------------------------------------------------------------------------- #

def _vial_to_participant(s: pd.Series) -> pd.Series:
    """Remap a pressure Series from TCGA vial IDs → 12-char participant IDs.

    Multi-vial participants are collapsed by mean (matches the parent repo convention).
    """
    out = pd.Series(s.values, index=[_participant(str(c)) for c in s.index], dtype=float)
    if out.index.duplicated().any():
        out = out.groupby(level=0).mean()
    return out


def _edge_partial_rho(
    pressure: pd.Series,
    rna_gene: pd.Series,
    cov: pd.DataFrame,
) -> Tuple[float, float, int]:
    from analysis.utils.common.loaders import partial_spearman

    shared = pressure.index.intersection(rna_gene.index).intersection(cov.index)
    if len(shared) < 20:
        return np.nan, np.nan, len(shared)
    return partial_spearman(rna_gene.loc[shared], pressure.loc[shared], cov.loc[shared])


def compute_edge_corr_table(
    edges: pd.DataFrame,
    tumor_mat: pd.DataFrame,
    rna: pd.DataFrame,
    clin: pd.DataFrame,
    prolif: pd.Series,
    genes: Sequence[str],
) -> pd.DataFrame:
    """Per (miRNA, gene): partial ρ + p under CPE+HRD+prolif, sorted by |ρ| desc."""
    print("[edge_partial_corr] computing pressure map for tumor …")
    pmap = compute_edge_pressure_map(
        edges, tumor_mat, genes=genes, expr_mode=ATTR_MODE, target_norm=TARGET_NORM,
    )
    # All downstream joins use 12-char participant IDs (deduplicated — multi-vial → mean in pressure)
    participants = pd.Index(sorted(set(_participant(str(c)) for c in tumor_mat.columns)))
    cov = _build_cov(clin, participants, prolif)

    from scipy.stats import spearmanr

    rows: List[dict] = []
    for (arm, gene), pressure_vial in pmap.items():
        if gene not in rna.index:
            continue
        # remap vial IDs → 12-char participant IDs so RNA and clinical join works
        pressure_s = _vial_to_participant(pressure_vial)
        rna_row = rna.loc[gene]
        if isinstance(rna_row, pd.DataFrame):
            rna_s = pd.to_numeric(rna_row.iloc[0], errors="coerce")
        else:
            rna_s = pd.to_numeric(rna_row, errors="coerce")
        shared = pressure_s.index.intersection(rna_s.index).intersection(cov.index)
        if len(shared) < 20:
            continue

        raw_rho, raw_p = spearmanr(
            pressure_s.loc[shared].values, rna_s.loc[shared].values,
        )
        rho_adj, p_adj, n_adj = _edge_partial_rho(pressure_s, rna_s, cov)

        ev_row = edges.loc[edges["miRNA"].eq(arm) & edges["gene"].eq(gene)]
        ev_score = float(ev_row["evidence_score"].iloc[0]) if not ev_row.empty else np.nan

        rows.append({
            "miRNA": arm,
            "gene": gene,
            "evidence_score": ev_score,
            "raw_rho": round(float(raw_rho), 4),
            "raw_p": float(raw_p),
            "rho_adj": round(float(rho_adj), 4) if pd.notna(rho_adj) else np.nan,
            "p_adj": float(p_adj) if pd.notna(p_adj) else np.nan,
            "n_adj": int(n_adj),
        })
        if len(rows) % 200 == 0:
            print(f"  … {len(rows)} edges processed")

    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    valid = df["p_adj"].notna()
    df.loc[valid, "q_adj"] = bh_fdr(df.loc[valid, "p_adj"].values)
    df["abs_rho_adj"] = df["rho_adj"].abs()
    df = df.sort_values("abs_rho_adj", ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    return df


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def _pam50_lookup(clin: pd.DataFrame, tumor_mat_cols: pd.Index) -> pd.Series:
    """PAM50 Series indexed by unique 12-char participant IDs derived from tumor vial cols."""
    participants = sorted(set(_participant(str(c)) for c in tumor_mat_cols))
    return clin.set_index("participant")["PAM50_final"].reindex(participants)


def _scatter_edge(
    ax,
    arm: str,
    gene: str,
    pressure_s: pd.Series,
    rna_s: pd.Series,
    pam50: pd.Series,
    *,
    rho_adj: float,
    p_adj: float,
    q_adj: float,
    n_adj: int,
    ev_score: float,
) -> None:
    from scipy.stats import linregress

    shared = pressure_s.index.intersection(rna_s.index)
    x_all = pd.to_numeric(pressure_s.loc[shared], errors="coerce")
    y_all = pd.to_numeric(rna_s.loc[shared], errors="coerce")

    for sub in PAM50_ORDER:
        mask = pam50.reindex(shared).eq(sub)
        xs = x_all.loc[mask]
        ys = y_all.loc[mask]
        ok = np.isfinite(xs.values) & np.isfinite(ys.values)
        if ok.sum() < 3:
            continue
        ax.scatter(
            xs.values[ok], ys.values[ok],
            s=5, alpha=0.30,
            color=PAM50_COLORS.get(sub, "#888888"),
            edgecolors="none", linewidths=0,
            label=sub, zorder=2,
        )

    # Regression line across all samples
    ok_all = np.isfinite(x_all.values) & np.isfinite(y_all.values)
    if ok_all.sum() > 20:
        slope, intercept, *_ = linregress(x_all.values[ok_all], y_all.values[ok_all])
        xr = np.array([x_all.values[ok_all].min(), x_all.values[ok_all].max()])
        ax.plot(xr, slope * xr + intercept, "k-", lw=0.9, alpha=0.55, zorder=3)

    # Annotation
    sig = ""
    if pd.notna(q_adj):
        if q_adj < 0.001:
            sig = "***"
        elif q_adj < 0.01:
            sig = "**"
        elif q_adj < 0.05:
            sig = "*"
    if pd.notna(rho_adj):
        rho_str = f"ρ={rho_adj:+.2f}{sig} (n={n_adj})"
    else:
        rho_str = "ρ=n.a."
    ax.annotate(
        rho_str,
        xy=(0.97, 0.03), xycoords="axes fraction",
        ha="right", va="bottom", fontsize=6.8,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"),
    )


def plot_partial_corr_panels(
    ranked_subset: Sequence[tuple],
    corr_table: pd.DataFrame,
    pmap: Dict[tuple, pd.Series],
    rna: pd.DataFrame,
    pam50: pd.Series,
    out_path: Path,
    *,
    cols: int = 4,
    row_height: float = 4.0,
    col_width: float = 5.0,
    dpi: int = 150,
    title: Optional[str] = None,
    page_label: Optional[str] = None,
) -> None:
    import matplotlib.pyplot as plt

    ranked_subset = list(ranked_subset)
    if not ranked_subset:
        return

    ncol = max(1, min(cols, len(ranked_subset)))
    nrow = int(np.ceil(len(ranked_subset) / ncol))

    fig_w = col_width * ncol
    fig_h = row_height * nrow + 1.2
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, fig_h), squeeze=False)

    suptitle = title or "Edge pressure → target expression  (partial Spearman)"
    if page_label:
        suptitle = f"{suptitle}  ({page_label})"
    fig.suptitle(
        f"{suptitle}\n"
        "X = tumor edge pressure (softmax_logrpm)  ·  Y = target log2(TPM+1)  ·  "
        "ρ adjusted: CPE + HRD + E2F/G2M  ·  * q<0.05  ** q<0.01  *** q<0.001",
        fontsize=9, y=1.01,
    )

    # Build lookup from corr_table: (miRNA, gene) → plain dict for safe access
    ct_idx: Dict[tuple, dict] = {}
    if not corr_table.empty:
        for _, row in corr_table.iterrows():
            ct_idx[(str(row["miRNA"]), str(row["gene"]))] = row.to_dict()

    for i, key in enumerate(ranked_subset):
        ax = axes[i // ncol][i % ncol]
        arm, gene = key
        pressure_raw = pmap.get(key, pd.Series(dtype=float))
        pressure_s = _vial_to_participant(pressure_raw) if len(pressure_raw) else pressure_raw
        if gene not in rna.index or len(pressure_s) == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8, color="#aaaaaa")
            ax.axis("off")
            continue

        rna_row = rna.loc[gene]
        if isinstance(rna_row, pd.DataFrame):
            rna_s = pd.to_numeric(rna_row.iloc[0], errors="coerce")
        else:
            rna_s = pd.to_numeric(rna_row, errors="coerce")
        row = ct_idx.get(key, {})
        rho_adj = float(row.get("rho_adj", np.nan)) if row else np.nan
        p_adj = float(row.get("p_adj", np.nan)) if row else np.nan
        q_adj = float(row.get("q_adj", np.nan)) if row else np.nan
        n_adj = int(row.get("n_adj", 0)) if row else 0
        ev_score = float(row.get("evidence_score", np.nan)) if row else np.nan

        _scatter_edge(
            ax, arm, gene, pressure_s, rna_s, pam50,
            rho_adj=rho_adj, p_adj=p_adj, q_adj=q_adj,
            n_adj=n_adj, ev_score=ev_score,
        )

        short_arm = arm.replace("hsa-", "")
        ev_str = f"  ev={ev_score:.0f}" if pd.notna(ev_score) else ""
        ax.set_title(f"{short_arm} → {gene}{ev_str}", fontsize=7.5, pad=4)
        if i % ncol == 0:
            ax.set_ylabel("target log2(TPM+1)", fontsize=7)
        if i >= (nrow - 1) * ncol:
            ax.set_xlabel("pressure (softmax_logrpm)", fontsize=7)
        ax.tick_params(labelsize=6.5)
        ax.grid(alpha=0.2, linestyle="--")
        ax.set_axisbelow(True)

    for j in range(len(ranked_subset), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    # PAM50 legend on the first page only
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=PAM50_COLORS[s],
               markersize=6, label=s)
        for s in PAM50_ORDER if s in PAM50_COLORS
    ]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.015))

    fig.subplots_adjust(top=0.92, hspace=0.52, wspace=0.35, bottom=0.06)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(
    *,
    out_dir: Path = OUT_DIR_CORR,
    out_name: str = "edge_partial_corr_panel",
    cols: int = 4,
    per_page: Optional[int] = None,
    row_height: float = 4.0,
    col_width: float = 5.0,
    max_edges: Optional[int] = 500,
    high_evidence_only: bool = True,
    hallmark: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    edges = _load_edges(high_evidence_only=high_evidence_only, hallmark=hallmark)
    if edges.empty:
        print("[edge_partial_corr_panels] no edges after filtering — exiting")
        return {}

    genes = sorted(edges["gene"].unique())
    print(f"[edge_partial_corr] edges={len(edges)}  genes={len(genes)}"
          f"  high_ev={high_evidence_only}  hallmark={hallmark}")

    hs = HallmarkSets.load()
    print("[edge_partial_corr] loading tumor miRNA + RNA + clinical …")
    states, _ = build_panel_context()
    tumor_mat = states["tumor"]
    rna = D.load_rna()
    clin = D.load_clinical_strata()
    prolif = _prolif_proxy(rna, hs)

    corr_table = compute_edge_corr_table(
        edges, tumor_mat, rna, clin, prolif, genes,
    )
    if corr_table.empty:
        print("[edge_partial_corr] no edges with sufficient n — exiting")
        return {}

    corr_table.to_csv(out_dir / f"{out_name}_corr_table.tsv", sep="\t", index=False)
    print(f"[edge_partial_corr] {len(corr_table)} edges with ρ; "
          f"q<0.05: {(corr_table['q_adj'].fillna(1) < 0.05).sum()}; "
          f"negative (repression): {(corr_table['rho_adj'].fillna(0) < 0).sum()}")

    # Re-compute tumor pressure map once more for the selected edges
    # (corr_table already filtered to edges with sufficient data, ranked by |ρ|)
    ranked_keys = [(r["miRNA"], r["gene"]) for _, r in corr_table.iterrows()]
    if max_edges is not None:
        ranked_keys = ranked_keys[:max_edges]

    print(f"[edge_partial_corr] computing pressure map for plotting ({len(ranked_keys)} edges) …")
    subset_edges = edges.loc[
        edges.apply(lambda r: (r["miRNA"], r["gene"]) in set(ranked_keys), axis=1)
    ]
    pmap = compute_edge_pressure_map(
        subset_edges, tumor_mat, genes=genes, expr_mode=ATTR_MODE, target_norm=TARGET_NORM,
    )

    pam50 = _pam50_lookup(clin, tumor_mat.columns)

    if per_page is None:
        per_page = cols * 10
    pages = _paginate(ranked_keys, per_page)
    n_pages = len(pages)

    figure_paths: List[str] = []
    for page_i, chunk in enumerate(pages, start=1):
        stem = out_name if n_pages == 1 else f"{out_name}_page{page_i:02d}"
        png = out_dir / "figures" / f"{stem}.png"
        page_label = (
            f"page {page_i}/{n_pages}, "
            f"edges {(page_i - 1) * per_page + 1}–{(page_i - 1) * per_page + len(chunk)}"
        )
        plot_partial_corr_panels(
            chunk,
            corr_table,
            pmap,
            rna,
            pam50,
            png,
            cols=cols,
            row_height=row_height,
            col_width=col_width,
            title=title,
            page_label=page_label,
        )
        figure_paths.append(str(png.relative_to(C.OUTPUT_ROOT)))
        print(f"[edge_partial_corr] page {page_i}/{n_pages} ({len(chunk)} edges) -> {png}")

    manifest = {
        "module": "mirna_hallmark.edge_partial_corr_panels",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_edges_total": len(edges),
        "n_edges_with_rho": len(corr_table),
        "n_edges_plotted": len(ranked_keys),
        "n_pages": n_pages,
        "per_page": per_page,
        "cols": cols,
        "attr_mode": ATTR_MODE,
        "target_norm": TARGET_NORM,
        "covariates": "CPE + HRD + E2F/G2M_prolif",
        "high_evidence_only": high_evidence_only,
        "hallmark_filter": hallmark,
        "max_edges": max_edges,
        "ranking": "|rho_adj| descending",
        "figures": figure_paths,
        "corr_table": f"{out_name}_corr_table.tsv",
    }
    (out_dir / f"{out_name}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"[edge_partial_corr] done — {len(ranked_keys)} edges, {n_pages} pages")
    return {"corr_table": corr_table, "pmap": pmap, "manifest": manifest}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--max-edges", type=int, default=500,
        help="cap on edges plotted, ranked by |ρ_adj| (default 500)"
    )
    ap.add_argument("--all-edges", action="store_true")
    ap.add_argument(
        "--high-evidence-only", action="store_true", default=True,
        help="restrict to high-evidence edges (default: on)"
    )
    ap.add_argument(
        "--no-high-evidence-filter", dest="high_evidence_only", action="store_false"
    )
    ap.add_argument("--hallmark", type=str, default=None)
    ap.add_argument("--cols", type=int, default=4)
    ap.add_argument("--per-page", type=int, default=None)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_CORR)
    ap.add_argument("--out-name", type=str, default="edge_partial_corr_panel")
    ap.add_argument("--row-height", type=float, default=4.0)
    ap.add_argument("--col-width", type=float, default=5.0)
    ap.add_argument("--title", type=str, default=None)
    args = ap.parse_args()

    C.ensure_output_dirs()
    run(
        out_dir=args.out_dir,
        out_name=args.out_name,
        cols=args.cols,
        per_page=args.per_page,
        row_height=args.row_height,
        col_width=args.col_width,
        max_edges=None if args.all_edges else args.max_edges,
        high_evidence_only=args.high_evidence_only,
        hallmark=args.hallmark,
        title=args.title,
    )


if __name__ == "__main__":
    main()
