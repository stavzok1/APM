"""Grid boxplot panels: per-edge (miRNA → gene) pressure across tissue states + PAM50.

Each subplot shows one miRNA→gene pair's realized pressure contribution c(m,g,s):

  - GTEx v10 healthy breast
  - TCGA-BRCA adjacent normal (NAT, 11)
  - TCGA-BRCA primary tumor (01, all or NAT-matched)
  - TCGA-BRCA primary tumor by PAM50 (LumA / LumB / Her2 / Basal)

Attribution mode: ``softmax_logrpm`` (no-z, cross-state comparable) with
``evidence_mass`` target normalisation — the same spine used in
``cross_state_landscape``.

Edges are ranked by |median_tumor_pressure - median_nat_pressure| (biggest
same-platform movers first), with mean tumor activity as a tiebreaker.  By
default only high-evidence edges are included.

Cross-platform note: GTEx uses TPM/URS vs TCGA RPM/MIMAT — GTEx pressure
reflects a different scale; use it for direction and rank, not absolute
comparison against TCGA values.  Tumor-vs-NAT is the fair same-platform
contrast.

Each subplot also shows the miRNA's median log2(RPM+1) expression per group as a
third tick-label line (↓med) so pressure magnitude and miRNA abundance can be read
together.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_pressure_panels
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_pressure_panels --max-edges 300 --cols 4
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_pressure_panels --hallmark HALLMARK_APOPTOSIS
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_pressure_panels --all-edges

Matched-sample run (NAT-paired tumor only, separate subfolder):
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_pressure_panels --matched
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Literal, NamedTuple, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_expression_panels import (
    PAM50_COLORS,
    PAM50_ORDER,
    STATE_COLORS,
    TumorScope,
    _filter_tumor_matrix,
    _overlay_dots,
    _paginate,
    _pam50_tumor_matrices,
)
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.pressure_engine import compute_edge_pressure_map

OUT_DIR = C.TISSUE_REFERENCE_DIR / "edge_pressure_panels"
ATTR_MODE = C.PRESSURE_ATTRIBUTION_EXPR_MODE  # softmax_logrpm — no-z, cross-state
TARGET_NORM = C.PRESSURE_TARGET_NORM           # evidence_mass


class PlotGroup(NamedTuple):
    key: str
    label: str
    color: str
    kind: str  # tissue | pam50


def _load_edges(*, high_evidence_only: bool = True, hallmark: Optional[str] = None) -> pd.DataFrame:
    """Load and (optionally) filter the Hallmark edge table.

    Deduplicates on (miRNA, gene), keeping the row with the highest
    evidence_score.  This prevents double-counting when a (miRNA, gene) pair
    spans multiple hallmark sets.
    """
    path = C.OUTPUT_ROOT / "edges" / "mirna_hallmark_edges.tsv.gz"
    e = pd.read_csv(path, sep="\t", low_memory=False)
    if hallmark:
        e = e.loc[e["hallmark_set"].eq(hallmark)]
    if high_evidence_only and "high_evidence" in e.columns:
        e = e.loc[e["high_evidence"]]
    # Deduplicate: multiple hallmark-set rows for the same (miRNA, gene) → keep best.
    e = (
        e.sort_values("evidence_score", ascending=False)
        .drop_duplicates(["miRNA", "gene"])
        .reset_index(drop=True)
    )
    return e


def _state_pressure_map(
    edges: pd.DataFrame,
    mat: pd.DataFrame,
    genes: Sequence[str],
) -> Dict[tuple, pd.Series]:
    """Compute per-sample pressure for every (miRNA, gene) in *edges* using *mat*."""
    if edges.empty or mat.empty:
        return {}
    return compute_edge_pressure_map(
        edges, mat, genes=genes, expr_mode=ATTR_MODE, target_norm=TARGET_NORM,
    )


def _score_and_rank(
    tumor_map: Dict[tuple, pd.Series],
    nat_map: Dict[tuple, pd.Series],
) -> List[tuple]:
    """Sort (miRNA, gene) keys by |tumor_median − nat_median|, then tumor activity."""
    all_keys = sorted(set(tumor_map.keys()) | set(nat_map.keys()))
    rows: List[dict] = []
    for k in all_keys:
        t = tumor_map.get(k, pd.Series(dtype=float))
        n = nat_map.get(k, pd.Series(dtype=float))
        tm = float(t.median()) if len(t) else 0.0
        nm = float(n.median()) if len(n) else 0.0
        ta = float(t.abs().mean()) if len(t) else 0.0
        rows.append({"edge": k, "delta": abs(tm - nm), "tumor_abs": ta})
    return (
        pd.DataFrame(rows)
        .sort_values(["delta", "tumor_abs"], ascending=False)
        ["edge"]
        .tolist()
    )


def build_panel_context(
    *, tumor_scope: TumorScope = "all",
) -> Tuple[Dict[str, pd.DataFrame], List[PlotGroup]]:
    tissue = _state_matrices()
    tumor = _filter_tumor_matrix(tissue["tumor"], tissue["nat"], tumor_scope=tumor_scope)
    tissue = {**tissue, "tumor": tumor}
    pam50 = _pam50_tumor_matrices(tumor)
    states: Dict[str, pd.DataFrame] = {**tissue, **pam50}

    tumor_label = "Tumor (matched)" if tumor_scope == "matched" else "Tumor"
    groups: List[PlotGroup] = [
        PlotGroup("gtex", "GTEx†", STATE_COLORS["gtex"], "tissue"),
        PlotGroup("nat", "NAT", STATE_COLORS["nat"], "tissue"),
        PlotGroup("tumor", tumor_label, STATE_COLORS["tumor"], "tissue"),
    ]
    for sub in PAM50_ORDER:
        if sub in states:
            groups.append(PlotGroup(sub, sub, PAM50_COLORS[sub], "pam50"))
    return states, groups


def _global_banner(states: Dict[str, pd.DataFrame], groups: Sequence[PlotGroup]) -> str:
    parts: List[str] = []
    for g in groups:
        mat = states.get(g.key)
        if mat is not None and not mat.empty:
            parts.append(f"{g.label} n={mat.shape[1]}")
    return "  |  ".join(parts)


def _edge_vals(
    key: tuple,
    state_maps: Dict[str, Dict[tuple, pd.Series]],
    group_key: str,
) -> np.ndarray:
    m = state_maps.get(group_key, {})
    s = m.get(key)
    if s is None or len(s) == 0:
        return np.array([], dtype=float)
    return pd.to_numeric(s, errors="coerce").fillna(0.0).to_numpy(dtype=float)


def _arm_expr_median(arm: str, states: Dict[str, pd.DataFrame], group_key: str) -> Optional[float]:
    """Median log2(RPM+1 or TPM+1) of *arm* in the given state matrix."""
    mat = states.get(group_key)
    if mat is None or mat.empty or arm not in mat.index:
        return None
    vals = pd.to_numeric(mat.loc[arm], errors="coerce").dropna()
    return float(vals.median()) if len(vals) else None


def plot_edge_pressure_panels(
    edges_subset: Sequence[tuple],
    edge_key_meta: pd.DataFrame,
    state_maps: Dict[str, Dict[tuple, pd.Series]],
    groups: Sequence[PlotGroup],
    states: Dict[str, pd.DataFrame],
    out_path: Path,
    *,
    cols: int = 4,
    row_height: float = 4.5,
    col_width: float = 5.8,
    dpi: int = 150,
    title: Optional[str] = None,
    page_label: Optional[str] = None,
) -> None:
    import matplotlib.pyplot as plt

    edges_subset = list(edges_subset)
    if not edges_subset:
        return

    active_groups = [g for g in groups if g.key in states and not states[g.key].empty]
    ncol = max(1, min(cols, len(edges_subset)))
    nrow = int(np.ceil(len(edges_subset) / ncol))

    fig_w = col_width * ncol
    fig_h = row_height * nrow + 1.4
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, fig_h), squeeze=False)

    banner = _global_banner(states, active_groups)
    suptitle = title or "Edge pressure across tissue states + PAM50"
    if page_label:
        suptitle = f"{suptitle}  ({page_label})"
    fig.suptitle(
        f"{suptitle}\n{banner}  ·  softmax_logrpm pressure  ·  ↓ = miRNA median log2(RPM+1)  ·  †GTEx TPM≠RPM",
        fontsize=10,
        y=1.01,
    )

    rng = np.random.default_rng(0)
    tissue_n = sum(1 for g in active_groups if g.kind == "tissue")

    # Build a lookup from key → display evidence score
    ev_lookup: dict = {}
    if not edge_key_meta.empty:
        for _, row in edge_key_meta.iterrows():
            k = (str(row["miRNA"]), str(row["gene"]))
            ev_lookup[k] = row.get("evidence_score", np.nan)

    for i, key in enumerate(edges_subset):
        ax = axes[i // ncol][i % ncol]
        arm, gene = key
        tick_labels: List[str] = []
        plot_groups: List[np.ndarray] = []
        colors: List[str] = []
        positions: List[float] = []

        pos = 1.0
        for g in active_groups:
            vals = _edge_vals(key, state_maps, g.key)
            mat = states.get(g.key)
            n_cohort = int(mat.shape[1]) if mat is not None and not mat.empty else 0
            has = len(vals) > 0
            n_label = n_cohort if has else 0
            expr_med = _arm_expr_median(arm, states, g.key)
            expr_str = f"\n↓{expr_med:.1f}" if expr_med is not None else "\n↓—"
            tick_labels.append(f"{g.label}\n(n={n_label}){expr_str}")
            plot_groups.append(vals)
            colors.append(g.color)
            positions.append(pos)
            pos += 1.0
            if g.kind == "tissue" and g.key == "tumor" and tissue_n == 3:
                pos += 0.45

        bp = ax.boxplot(
            plot_groups,
            positions=positions,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.1},
            zorder=2,
        )
        for patch, color, vals in zip(bp["boxes"], colors, plot_groups):
            if len(vals) == 0:
                patch.set_visible(False)
            else:
                patch.set_facecolor(color)
                patch.set_alpha(0.55)
                patch.set_edgecolor("#333333")

        for p, vals, color in zip(positions, plot_groups, colors):
            if len(vals) == 0:
                ax.text(p, 0.0, "—", ha="center", va="bottom", fontsize=10, color="#999999")
            else:
                _overlay_dots(ax, p, vals, color, rng=rng)

        if tissue_n == 3 and len(positions) >= 4:
            sep = (positions[2] + positions[3]) / 2.0
            ax.axvline(sep, color="#bbbbbb", linestyle=":", linewidth=0.9, zorder=1)

        ax.axhline(0, color="#aaaaaa", linewidth=0.6, linestyle="--", zorder=1)
        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels, fontsize=7)

        # Compact miRNA name for readability
        short_arm = arm.replace("hsa-", "")
        ev = ev_lookup.get(key, np.nan)
        ev_str = f"  ev={ev:.0f}" if pd.notna(ev) else ""
        ax.set_title(f"{short_arm} → {gene}{ev_str}", fontsize=7.5, pad=4)

        if i % ncol == 0:
            ax.set_ylabel("pressure (softmax_logrpm)")
        ax.grid(axis="y", alpha=0.25, linestyle="--")
        ax.set_axisbelow(True)
        ax.margins(y=0.08)

    for j in range(len(edges_subset), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    fig.subplots_adjust(top=0.93, hspace=0.55, wspace=0.32)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


OUT_DIR_MATCHED = C.TISSUE_REFERENCE_DIR / "edge_pressure_panels_matched"


def run(
    *,
    out_dir: Path = OUT_DIR,
    out_name: str = "edge_pressure_panel",
    cols: int = 4,
    per_page: Optional[int] = None,
    row_height: float = 4.5,
    col_width: float = 5.8,
    max_edges: Optional[int] = 500,
    high_evidence_only: bool = True,
    hallmark: Optional[str] = None,
    tumor_scope: TumorScope = "all",
    title: Optional[str] = None,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    edges = _load_edges(high_evidence_only=high_evidence_only, hallmark=hallmark)
    if edges.empty:
        print("[edge_pressure_panels] no edges found after filtering — exiting")
        return {}

    genes = sorted(edges["gene"].unique())
    print(f"[edge_pressure_panels] edges={len(edges)}  genes={len(genes)}"
          f"  high_ev_only={high_evidence_only}  hallmark={hallmark}")

    states, groups = build_panel_context(tumor_scope=tumor_scope)
    print(f"[edge_pressure_panels] states: { {k: v.shape[1] for k,v in states.items()} }")

    # Compute pressure maps per state (this is the expensive step)
    print("[edge_pressure_panels] computing pressure maps per state …")
    state_maps: Dict[str, Dict[tuple, pd.Series]] = {}
    for g in groups:
        mat = states.get(g.key)
        if mat is None or mat.empty:
            state_maps[g.key] = {}
            continue
        state_maps[g.key] = _state_pressure_map(edges, mat, genes)
        n = len(state_maps[g.key])
        print(f"  {g.key}: {n} (miRNA, gene) pairs with pressure")

    # Rank edges by interestingness (|tumor_median − nat_median| then tumor activity)
    ranked = _score_and_rank(state_maps.get("tumor", {}), state_maps.get("nat", {}))
    if max_edges is not None:
        ranked = ranked[:max_edges]

    print(f"[edge_pressure_panels] ranked edges to plot: {len(ranked)}")

    if per_page is None:
        per_page = cols * 12
    pages = _paginate(ranked, per_page)
    n_pages = len(pages)

    figure_paths: List[str] = []
    for page_i, chunk in enumerate(pages, start=1):
        if n_pages == 1:
            stem = out_name
        else:
            stem = f"{out_name}_page{page_i:02d}"
        png = fig_dir / f"{stem}.png"
        page_label = (
            f"page {page_i}/{n_pages}, "
            f"edges {(page_i - 1) * per_page + 1}–{(page_i - 1) * per_page + len(chunk)}"
        )
        plot_edge_pressure_panels(
            chunk,
            edges,
            state_maps,
            groups,
            states,
            png,
            cols=cols,
            row_height=row_height,
            col_width=col_width,
            title=title,
            page_label=page_label,
        )
        figure_paths.append(str(png.relative_to(C.OUTPUT_ROOT)))
        print(f"[edge_pressure_panels] page {page_i}/{n_pages} ({len(chunk)} edges) -> {png}")

    # Summary TSV: one row per ranked edge with key stats
    summary_rows: List[dict] = []
    for rank_i, k in enumerate(ranked, start=1):
        arm, gene = k
        ev_row = edges.loc[edges["miRNA"].eq(arm) & edges["gene"].eq(gene)]
        ev_score = float(ev_row["evidence_score"].iloc[0]) if not ev_row.empty else np.nan
        row: dict = {"rank": rank_i, "miRNA": arm, "gene": gene, "evidence_score": ev_score}
        for g in groups:
            vals = _edge_vals(k, state_maps, g.key)
            row[f"median_{g.key}"] = float(np.median(vals)) if len(vals) else np.nan
            row[f"n_{g.key}"] = len(vals)
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / f"{out_name}_summary.tsv", sep="\t", index=False)

    state_counts = {g.key: int(states[g.key].shape[1]) for g in groups if g.key in states}
    manifest = {
        "module": "mirna_hallmark.analyses.edge_panels.edge_pressure_panels",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_edges_total": len(edges),
        "n_edges_plotted": len(ranked),
        "n_pages": n_pages,
        "per_page": per_page,
        "cols": cols,
        "attr_mode": ATTR_MODE,
        "target_norm": TARGET_NORM,
        "high_evidence_only": high_evidence_only,
        "hallmark_filter": hallmark,
        "max_edges": max_edges,
        "tumor_scope": tumor_scope,
        "group_counts": state_counts,
        "ranking": "|median_tumor_pressure - median_nat_pressure|, then mean_abs_tumor",
        "cross_platform_note": "GTEx TPM/URS vs TCGA RPM/MIMAT — shape/rank only for GTEx",
        "figures": figure_paths,
        "summary": f"{out_name}_summary.tsv",
    }
    (out_dir / f"{out_name}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"[edge_pressure_panels] done — {len(ranked)} edges, {n_pages} pages")
    return {"edges": edges, "ranked": ranked, "summary": summary, "manifest": manifest}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--max-edges", type=int, default=500,
        help="cap on edges plotted (ranked by interestingness, default 500)"
    )
    ap.add_argument(
        "--all-edges", action="store_true",
        help="plot all filtered edges without a max cap"
    )
    ap.add_argument(
        "--matched", action="store_true",
        help="shortcut: matched tumor-scope + separate output subfolder (edge_pressure_panels_matched/)"
    )
    ap.add_argument(
        "--high-evidence-only", action="store_true", default=True,
        help="restrict to high-evidence edges only (default: on)"
    )
    ap.add_argument(
        "--no-high-evidence-filter", dest="high_evidence_only", action="store_false",
        help="include all evidence levels"
    )
    ap.add_argument(
        "--hallmark", type=str, default=None,
        help="restrict to a single Hallmark set (e.g. HALLMARK_APOPTOSIS)"
    )
    ap.add_argument("--cols", type=int, default=4, help="grid columns per page")
    ap.add_argument("--per-page", type=int, default=None, help="edges per figure page")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--out-name", type=str, default="edge_pressure_panel")
    ap.add_argument("--row-height", type=float, default=4.0)
    ap.add_argument("--col-width", type=float, default=5.8)
    ap.add_argument(
        "--tumor-scope", choices=["all", "matched"], default="all",
        help="primary tumor cohort: all TCGA 01 samples or NAT-matched only"
    )
    ap.add_argument("--title", type=str, default=None)
    args = ap.parse_args()

    C.ensure_output_dirs()
    out_dir = args.out_dir
    out_name = args.out_name
    tumor_scope = args.tumor_scope
    if args.matched:
        out_dir = OUT_DIR_MATCHED
        out_name = "edge_pressure_panel_matched"
        tumor_scope = "matched"
    run(
        out_dir=out_dir,
        out_name=out_name,
        cols=args.cols,
        per_page=args.per_page,
        row_height=args.row_height,
        col_width=args.col_width,
        max_edges=None if args.all_edges else args.max_edges,
        high_evidence_only=args.high_evidence_only,
        hallmark=args.hallmark,
        tumor_scope=tumor_scope,
        title=args.title,
    )


if __name__ == "__main__":
    main()
