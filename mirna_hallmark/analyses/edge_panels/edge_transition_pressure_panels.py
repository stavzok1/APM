"""Per-edge transition-pressure panels: sample vs reference-layer median.

For each edge (miRNA → gene) and each tissue-state transition, every individual
sample's pressure is compared against the **population median** of the previous
layer:

  H→NAT   :  Δ(s) = pressure_nat(s)    − median(pressure_healthy)
  NAT→T   :  Δ(s) = pressure_tumor(s)  − median(pressure_nat)
  H→T     :  Δ(s) = pressure_tumor(s)  − median(pressure_healthy)

Each subplot shows three boxplot groups (one per transition) so the pressure
cascade healthy → NAT → tumor can be read in a single panel.

Two variants are produced:
  * all samples   — all GTEx, all NAT, all primary tumors
  * matched       — restrict to the ~103 participants with both NAT and tumor;
                    reference medians are still computed from the full NAT and
                    GTEx populations (population-level baseline, not matched-only)

Attribution mode: ``softmax_logrpm`` (no-z, cross-state comparable).
Edges sorted by NAT→T median Δ descending (most acquired in tumor).

Run:
  # all-samples version
  .venv/bin/python3 -m mirna_hallmark.edge_transition_pressure_panels
  # matched-samples version
  .venv/bin/python3 -m mirna_hallmark.edge_transition_pressure_panels --matched
  # single hallmark
  .venv/bin/python3 -m mirna_hallmark.edge_transition_pressure_panels --hallmark HALLMARK_APOPTOSIS
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_expression_panels import _overlay_dots, _paginate
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.analyses.edge_panels.edge_pressure_panels import (
    TARGET_NORM,
    ATTR_MODE,
    _load_edges,
)
from mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels import (
    _matched_pair_matrices,
    _state_pressure_map,
    _wilcoxon_p,
)
from mirna_hallmark.stats import bh_fdr

OUT_DIR_ALL = C.TISSUE_REFERENCE_DIR / "edge_transition_pressure_panels"
OUT_DIR_MATCHED = C.TISSUE_REFERENCE_DIR / "edge_transition_pressure_panels_matched"

# Colours per transition
TRANS_COLORS = {
    "H→NAT":  "#EDC948",   # NAT yellow
    "NAT→T":  "#E15759",   # tumor red
    "H→T":    "#7B4F9E",   # deep purple (combined)
}


# --------------------------------------------------------------------------- #
# Transition delta computation
# --------------------------------------------------------------------------- #

class TransitionMaps(NamedTuple):
    h_to_nat:   Dict[tuple, pd.Series]   # nat - healthy_median
    nat_to_t:   Dict[tuple, pd.Series]   # tumor - nat_median
    h_to_t:     Dict[tuple, pd.Series]   # tumor - healthy_median
    n_healthy:  int
    n_nat:      int
    n_tumor:    int


def _build_transitions(
    gtex_map: Dict[tuple, pd.Series],
    nat_map: Dict[tuple, pd.Series],
    tumor_map: Dict[tuple, pd.Series],
) -> TransitionMaps:
    """Compute per-edge transition delta Series.

    Reference medians are derived from the full *nat_map* and *gtex_map*
    populations (the ``median_of_the_previous_layer``), regardless of whether
    *tumor_map* is restricted to matched samples.
    """
    all_keys = sorted(
        set(gtex_map.keys()) | set(nat_map.keys()) | set(tumor_map.keys())
    )
    h_to_nat: Dict[tuple, pd.Series] = {}
    nat_to_t: Dict[tuple, pd.Series] = {}
    h_to_t: Dict[tuple, pd.Series] = {}

    for k in all_keys:
        gtex_s = pd.to_numeric(gtex_map.get(k, pd.Series(dtype=float)), errors="coerce")
        nat_s  = pd.to_numeric(nat_map.get(k, pd.Series(dtype=float)),  errors="coerce")
        tum_s  = pd.to_numeric(tumor_map.get(k, pd.Series(dtype=float)), errors="coerce")

        healthy_med = float(gtex_s.median()) if len(gtex_s) else np.nan
        nat_med     = float(nat_s.median())  if len(nat_s)  else np.nan

        if len(nat_s) >= 5 and pd.notna(healthy_med):
            h_to_nat[k] = nat_s - healthy_med

        if len(tum_s) >= 5 and pd.notna(nat_med):
            nat_to_t[k] = tum_s - nat_med

        if len(tum_s) >= 5 and pd.notna(healthy_med):
            h_to_t[k] = tum_s - healthy_med

    # infer sample sizes from a representative key
    def _n(m: Dict[tuple, pd.Series]) -> int:
        for s in m.values():
            return len(s)
        return 0

    return TransitionMaps(
        h_to_nat=h_to_nat, nat_to_t=nat_to_t, h_to_t=h_to_t,
        n_healthy=_n(gtex_map), n_nat=_n(nat_map), n_tumor=_n(tumor_map),
    )


def _score_and_rank(tm: TransitionMaps) -> List[tuple]:
    """Sort by NAT→T median Δ descending (most tumor-acquired first)."""
    keys = sorted(set(tm.nat_to_t.keys()))
    rows = []
    for k in keys:
        med = float(tm.nat_to_t[k].median()) if len(tm.nat_to_t.get(k, [])) else 0.0
        rows.append({"edge": k, "median_nat_to_t": med})
    return (
        pd.DataFrame(rows)
        .sort_values("median_nat_to_t", ascending=False)["edge"]
        .tolist()
    )


# --------------------------------------------------------------------------- #
# Summary / stats
# --------------------------------------------------------------------------- #

def _build_summary(
    ranked: Sequence[tuple],
    tm: TransitionMaps,
    edges: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[dict] = []
    p_vals_h2n: List[float] = []
    p_vals_n2t: List[float] = []
    p_vals_h2t: List[float] = []

    for k in ranked:
        arm, gene = k
        ev_row = edges.loc[edges["miRNA"].eq(arm) & edges["gene"].eq(gene)]
        ev = float(ev_row["evidence_score"].iloc[0]) if not ev_row.empty else np.nan

        def _stats(m: Dict[tuple, pd.Series]) -> Tuple[float, float, float]:
            s = m.get(k, pd.Series(dtype=float))
            if len(s) == 0:
                return np.nan, np.nan, np.nan
            return float(s.median()), float(s.quantile(0.25)), float(s.quantile(0.75))

        m_h2n, q25_h2n, q75_h2n = _stats(tm.h_to_nat)
        m_n2t, q25_n2t, q75_n2t = _stats(tm.nat_to_t)
        m_h2t, q25_h2t, q75_h2t = _stats(tm.h_to_t)

        p_h2n = _wilcoxon_p(tm.h_to_nat.get(k, pd.Series(dtype=float)))
        p_n2t = _wilcoxon_p(tm.nat_to_t.get(k, pd.Series(dtype=float)))
        p_h2t = _wilcoxon_p(tm.h_to_t.get(k, pd.Series(dtype=float)))
        p_vals_h2n.append(p_h2n)
        p_vals_n2t.append(p_n2t)
        p_vals_h2t.append(p_h2t)

        rows.append({
            "miRNA": arm, "gene": gene, "evidence_score": ev,
            "n_h_to_nat": len(tm.h_to_nat.get(k, [])),
            "median_h_to_nat": m_h2n, "q25_h_to_nat": q25_h2n, "q75_h_to_nat": q75_h2n,
            "p_h_to_nat": p_h2n,
            "n_nat_to_t": len(tm.nat_to_t.get(k, [])),
            "median_nat_to_t": m_n2t, "q25_nat_to_t": q25_n2t, "q75_nat_to_t": q75_n2t,
            "p_nat_to_t": p_n2t,
            "n_h_to_t": len(tm.h_to_t.get(k, [])),
            "median_h_to_t": m_h2t, "q25_h_to_t": q25_h2t, "q75_h_to_t": q75_h2t,
            "p_h_to_t": p_h2t,
        })

    df = pd.DataFrame(rows)
    for col, pvals in [("q_h_to_nat", p_vals_h2n),
                       ("q_nat_to_t", p_vals_n2t),
                       ("q_h_to_t",   p_vals_h2t)]:
        valid = pd.Series(pvals).notna()
        q = pd.Series([np.nan] * len(pvals))
        if valid.any():
            q[valid] = bh_fdr(pd.Series(pvals)[valid].values)
        df[col] = q.values
    return df


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def _sig_str(q: float) -> str:
    if pd.isna(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def _draw_transition_box(
    ax,
    vals: np.ndarray,
    pos: float,
    color: str,
    rng: np.random.Generator,
    *,
    q: float = np.nan,
) -> None:
    if len(vals) == 0:
        ax.text(pos, 0, "—", ha="center", va="bottom", fontsize=9, color="#aaaaaa")
        return
    bp = ax.boxplot(
        [vals], positions=[pos], widths=0.48,
        patch_artist=True, showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2}, zorder=2,
    )
    bp["boxes"][0].set_facecolor(color)
    bp["boxes"][0].set_alpha(0.55)
    bp["boxes"][0].set_edgecolor("#333333")
    _overlay_dots(ax, pos, vals, color, rng=rng)
    sig = _sig_str(q)
    if sig:
        ymax = float(np.percentile(vals, 75)) + 1.5 * float(
            np.percentile(vals, 75) - np.percentile(vals, 25)
        )
        ax.text(pos, ymax, sig, ha="center", va="bottom", fontsize=9, color="#333333")


def plot_transition_panels(
    edges_subset: Sequence[tuple],
    tm: TransitionMaps,
    summary: pd.DataFrame,
    out_path: Path,
    *,
    cols: int = 4,
    row_height: float = 4.3,
    col_width: float = 5.5,
    dpi: int = 150,
    title: Optional[str] = None,
    page_label: Optional[str] = None,
    matched: bool = False,
) -> None:
    import matplotlib.pyplot as plt

    edges_subset = list(edges_subset)
    if not edges_subset:
        return

    ncol = max(1, min(cols, len(edges_subset)))
    nrow = int(np.ceil(len(edges_subset) / ncol))
    fig_w = col_width * ncol
    fig_h = row_height * nrow + 1.4
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, fig_h), squeeze=False)

    scope = "matched-pairs" if matched else "all samples"
    suptitle = title or f"Edge pressure cascade: healthy → NAT → tumor  ({scope})"
    if page_label:
        suptitle = f"{suptitle}  ({page_label})"
    fig.suptitle(
        f"{suptitle}\n"
        f"GTEx n={tm.n_healthy}  NAT n={tm.n_nat}  Tumor n={tm.n_tumor}  ·  "
        "Δ = sample pressure − reference-layer median  ·  "
        "softmax_logrpm  ·  * q<0.05 Wilcoxon",
        fontsize=9.5, y=1.01,
    )

    # pre-build lookup
    meta_idx: dict = {}
    if not summary.empty:
        for _, row in summary.iterrows():
            meta_idx[(str(row["miRNA"]), str(row["gene"]))] = row.to_dict()

    rng = np.random.default_rng(0)

    TRANSITIONS = [
        ("H→NAT",  tm.h_to_nat,  TRANS_COLORS["H→NAT"],  "q_h_to_nat"),
        ("NAT→T",  tm.nat_to_t,  TRANS_COLORS["NAT→T"],  "q_nat_to_t"),
        ("H→T",    tm.h_to_t,    TRANS_COLORS["H→T"],    "q_h_to_t"),
    ]

    for i, key in enumerate(edges_subset):
        ax = axes[i // ncol][i % ncol]
        arm, gene = key
        meta = meta_idx.get(key, {})

        positions = [1.0, 2.1, 3.2]
        tick_labels = []

        for pos, (label, dmap, color, q_col) in zip(positions, TRANSITIONS):
            vals = pd.to_numeric(dmap.get(key, pd.Series(dtype=float)),
                                 errors="coerce").dropna().values
            q = meta.get(q_col, np.nan)
            _draw_transition_box(ax, vals, pos, color, rng, q=q)
            med = float(np.median(vals)) if len(vals) else np.nan
            n = len(vals)
            med_str = f"\nΔ={med:+.1f}" if pd.notna(med) else ""
            tick_labels.append(f"{label}\n(n={n}){med_str}")

        # separator lines between transitions
        ax.axvline(1.55, color="#cccccc", lw=0.8, linestyle=":", zorder=1)
        ax.axvline(2.65, color="#cccccc", lw=0.8, linestyle=":", zorder=1)
        ax.axhline(0, color="#888888", lw=0.7, linestyle="--", zorder=1)

        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels, fontsize=6.8)
        ax.set_xlim(0.4, 3.8)

        ev = meta.get("evidence_score", np.nan)
        ev_str = f"  ev={ev:.0f}" if pd.notna(ev) else ""
        short_arm = arm.replace("hsa-", "")
        ax.set_title(f"{short_arm} → {gene}{ev_str}", fontsize=7.5, pad=4)
        if i % ncol == 0:
            ax.set_ylabel("Δ pressure (vs layer median)")
        ax.grid(axis="y", alpha=0.22, linestyle="--")
        ax.set_axisbelow(True)
        ax.margins(y=0.12)

    for j in range(len(edges_subset), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    # legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=TRANS_COLORS[t], alpha=0.6, label=t)
               for t in ("H→NAT", "NAT→T", "H→T")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))

    fig.subplots_adjust(top=0.93, hspace=0.65, wspace=0.35, bottom=0.05)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(
    *,
    out_dir: Path = OUT_DIR_ALL,
    out_name: str = "edge_transition_pressure_panel",
    cols: int = 4,
    per_page: Optional[int] = None,
    row_height: float = 4.3,
    col_width: float = 5.5,
    max_edges: Optional[int] = 500,
    high_evidence_only: bool = True,
    hallmark: Optional[str] = None,
    matched: bool = False,
    title: Optional[str] = None,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    edges = _load_edges(high_evidence_only=high_evidence_only, hallmark=hallmark)
    if edges.empty:
        print("[edge_transition] no edges after filtering — exiting")
        return {}
    genes = sorted(edges["gene"].unique())
    scope = "matched" if matched else "all"
    print(f"[edge_transition] edges={len(edges)}  genes={len(genes)}"
          f"  scope={scope}  hallmark={hallmark}")

    tissue = _state_matrices()
    gtex_mat  = tissue.get("gtex",  pd.DataFrame())
    nat_mat   = tissue["nat"]
    tumor_mat = tissue["tumor"]

    if matched:
        matched_tumor_mat, _, matched_participants = _matched_pair_matrices()
        # Use matched tumor columns; reference medians from full NAT + GTEx
        tumor_mat_for_map = matched_tumor_mat
    else:
        tumor_mat_for_map = tumor_mat

    print(f"[edge_transition] computing pressure maps "
          f"(gtex={gtex_mat.shape[1] if not gtex_mat.empty else 0}, "
          f"nat={nat_mat.shape[1]}, tumor={tumor_mat_for_map.shape[1]}) …")

    gtex_map  = _state_pressure_map(edges, gtex_mat, genes) if not gtex_mat.empty else {}
    nat_map   = _state_pressure_map(edges, nat_mat, genes)
    tumor_map = _state_pressure_map(edges, tumor_mat_for_map, genes)
    print(f"  gtex={len(gtex_map)}  nat={len(nat_map)}  tumor={len(tumor_map)} edge-keys")

    tm = _build_transitions(gtex_map, nat_map, tumor_map)
    ranked = _score_and_rank(tm)
    if max_edges is not None:
        ranked = ranked[:max_edges]
    print(f"[edge_transition] ranked edges: {len(ranked)}")

    summary = _build_summary(ranked, tm, edges)
    summary.to_csv(out_dir / f"{out_name}_summary.tsv", sep="\t", index=False)
    for col, label in [("q_h_to_nat", "H→NAT"), ("q_nat_to_t", "NAT→T"), ("q_h_to_t", "H→T")]:
        if col in summary:
            n_sig = int((summary[col].fillna(1) < 0.05).sum())
            print(f"  q<0.05 {label}: {n_sig}")

    if per_page is None:
        per_page = cols * 11
    pages = _paginate(ranked, per_page)
    n_pages = len(pages)

    figure_paths: List[str] = []
    for page_i, chunk in enumerate(pages, start=1):
        stem = out_name if n_pages == 1 else f"{out_name}_page{page_i:02d}"
        png = out_dir / "figures" / f"{stem}.png"
        page_label = (
            f"page {page_i}/{n_pages}, "
            f"edges {(page_i-1)*per_page+1}–{(page_i-1)*per_page+len(chunk)}"
        )
        plot_transition_panels(
            chunk, tm, summary, png,
            cols=cols, row_height=row_height, col_width=col_width,
            title=title, page_label=page_label, matched=matched,
        )
        figure_paths.append(str(png.relative_to(C.OUTPUT_ROOT)))
        print(f"[edge_transition] page {page_i}/{n_pages} ({len(chunk)} edges) -> {png}")

    manifest = {
        "module": "mirna_hallmark.edge_transition_pressure_panels",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": scope,
        "n_gtex": gtex_mat.shape[1] if not gtex_mat.empty else 0,
        "n_nat": nat_mat.shape[1],
        "n_tumor": tumor_mat_for_map.shape[1],
        "transitions": {
            "H→NAT":  "pressure_nat(s) − median(pressure_healthy)",
            "NAT→T":  "pressure_tumor(s) − median(pressure_nat)",
            "H→T":    "pressure_tumor(s) − median(pressure_healthy)",
        },
        "ranking": "median NAT→T delta descending",
        "attr_mode": ATTR_MODE,
        "target_norm": TARGET_NORM,
        "high_evidence_only": high_evidence_only,
        "hallmark_filter": hallmark,
        "max_edges": max_edges,
        "figures": figure_paths,
        "summary": f"{out_name}_summary.tsv",
    }
    (out_dir / f"{out_name}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"[edge_transition] done — {len(ranked)} edges, {n_pages} pages")
    return {"tm": tm, "summary": summary, "manifest": manifest}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-edges", type=int, default=500)
    ap.add_argument("--all-edges", action="store_true")
    ap.add_argument("--high-evidence-only", action="store_true", default=True)
    ap.add_argument("--no-high-evidence-filter", dest="high_evidence_only",
                    action="store_false")
    ap.add_argument("--hallmark", type=str, default=None)
    ap.add_argument("--matched", action="store_true",
                    help="restrict tumor samples to matched (NAT-paired) subset; "
                         "reference medians remain from full NAT + GTEx populations")
    ap.add_argument("--cols", type=int, default=4)
    ap.add_argument("--per-page", type=int, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--out-name", type=str, default=None)
    ap.add_argument("--row-height", type=float, default=4.3)
    ap.add_argument("--col-width", type=float, default=5.5)
    ap.add_argument("--title", type=str, default=None)
    args = ap.parse_args()

    out_dir = args.out_dir or (OUT_DIR_MATCHED if args.matched else OUT_DIR_ALL)
    out_name = args.out_name or (
        "edge_transition_pressure_panel_matched" if args.matched
        else "edge_transition_pressure_panel"
    )

    C.ensure_output_dirs()
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
        matched=args.matched,
        title=args.title,
    )


if __name__ == "__main__":
    main()
