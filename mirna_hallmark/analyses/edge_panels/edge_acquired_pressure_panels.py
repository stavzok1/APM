"""Per-edge acquired-pressure panels: paired (tumor − NAT) pressure delta.

For each participant who has both a primary tumor (01) and an adjacent normal
(NAT, 11) miRNA profile, we compute the realized edge pressure c(m,g,s) in
each tissue type and take the **within-participant difference**:

    Δ_acquired(m, g, p) = pressure_tumor(m, g, p) − pressure_nat(m, g, p)

This is the pressure the tumor *acquired* over the paired normal tissue.  It
controls for participant-level variation (germline, general miRNA landscape) in
a way that comparing two independent boxplots cannot.

Each subplot shows:
  - A single boxplot of Δ_acquired across all matched participants
  - Optionally split by PAM50 subtype (within matched subset)
  - Y = 0 reference line: above = acquired; below = lost relative to NAT
  - Annotated with median Δ and Wilcoxon signed-rank p (one-sample, H₀: Δ=0)

Attribution mode: ``softmax_logrpm`` (no-z, cross-state comparable).
Edges sorted by median Δ descending (most acquired first).

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels --max-edges 300
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels --hallmark HALLMARK_APOPTOSIS
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
from mirna_hallmark.analyses.cross_state.cross_state_expression_panels import (
    PAM50_COLORS,
    PAM50_ORDER,
    _overlay_dots,
    _paginate,
    _pam50_tumor_matrices,
)
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.analyses.edge_panels.edge_pressure_panels import (
    TARGET_NORM,
    ATTR_MODE,
    _load_edges,
)
from mirna_hallmark.family_normal_reference import _load_full_mirna, _participant, _split_types
from mirna_hallmark.pressure_engine import compute_edge_pressure_map
from mirna_hallmark.stats import bh_fdr

OUT_DIR = C.TISSUE_REFERENCE_DIR / "edge_acquired_pressure_panels"


# --------------------------------------------------------------------------- #
# Matched-pair miRNA matrices
# --------------------------------------------------------------------------- #

def _matched_pair_matrices() -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """(tumor_mat, nat_mat, matched_participants) for participants with both tissue types.

    Returns matrices restricted to the same column order so column i in tumor_mat
    and column i in nat_mat correspond to the SAME participant.
    """
    full = _load_full_mirna()
    split = _split_types(full)
    tumor_all = split["tumor"]
    nat_all = split["nat"]

    # Map participant → columns in each matrix
    t_by_p: Dict[str, List[str]] = {}
    for c in tumor_all.columns:
        t_by_p.setdefault(_participant(str(c)), []).append(str(c))
    n_by_p: Dict[str, List[str]] = {}
    for c in nat_all.columns:
        n_by_p.setdefault(_participant(str(c)), []).append(str(c))

    matched = sorted(set(t_by_p.keys()) & set(n_by_p.keys()))
    if not matched:
        raise ValueError("No participants found with both tumor and NAT samples")

    # One vial per participant: first sorted column
    t_cols = [sorted(t_by_p[p])[0] for p in matched]
    n_cols = [sorted(n_by_p[p])[0] for p in matched]

    tumor_mat = tumor_all[t_cols].copy()
    tumor_mat.columns = pd.Index(matched)  # reindex to participant IDs
    nat_mat = nat_all[n_cols].copy()
    nat_mat.columns = pd.Index(matched)    # same order → matched columns

    return tumor_mat, nat_mat, matched


# --------------------------------------------------------------------------- #
# Paired delta computation
# --------------------------------------------------------------------------- #

def _state_pressure_map(
    edges: pd.DataFrame,
    mat: pd.DataFrame,
    genes: Sequence[str],
) -> Dict[tuple, pd.Series]:
    if edges.empty or mat.empty:
        return {}
    return compute_edge_pressure_map(
        edges, mat, genes=genes, expr_mode=ATTR_MODE, target_norm=TARGET_NORM,
    )


def _paired_delta(
    tumor_map: Dict[tuple, pd.Series],
    nat_map: Dict[tuple, pd.Series],
) -> Dict[tuple, pd.Series]:
    """For each edge key, compute per-participant Δ = tumor_pressure − nat_pressure.

    Only includes participants present in both maps' Series (inner join on index).
    """
    all_keys = sorted(set(tumor_map.keys()) & set(nat_map.keys()))
    out: Dict[tuple, pd.Series] = {}
    for k in all_keys:
        t = pd.to_numeric(tumor_map[k], errors="coerce")
        n = pd.to_numeric(nat_map[k], errors="coerce")
        shared = t.index.intersection(n.index)
        if len(shared) < 5:
            continue
        delta = t.loc[shared] - n.loc[shared]
        out[k] = delta
    return out


def _score_and_rank(delta_map: Dict[tuple, pd.Series]) -> List[tuple]:
    """Sort edges by median Δ_acquired descending."""
    rows = []
    for k, d in delta_map.items():
        med = float(d.median()) if len(d) else 0.0
        abs_med = abs(med)
        rows.append({"edge": k, "median_delta": med, "abs_median": abs_med})
    return (
        pd.DataFrame(rows)
        .sort_values("median_delta", ascending=False)
        ["edge"]
        .tolist()
    )


def _wilcoxon_p(delta: pd.Series) -> float:
    """One-sample Wilcoxon signed-rank test: H₀ = median Δ = 0."""
    from scipy.stats import wilcoxon
    vals = pd.to_numeric(delta, errors="coerce").dropna().values
    if len(vals) < 10 or np.all(vals == 0):
        return np.nan
    try:
        _, p = wilcoxon(vals, zero_method="wilcox", alternative="two-sided")
        return float(p)
    except ValueError:
        return np.nan


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def _pam50_series(matched_participants: Sequence[str]) -> pd.Series:
    clin = D.load_clinical_strata()
    return clin.set_index("participant")["PAM50_final"].reindex(matched_participants)


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


def plot_acquired_pressure_panels(
    edges_subset: Sequence[tuple],
    delta_map: Dict[tuple, pd.Series],
    nat_vs_healthy_map: Dict[tuple, pd.Series],
    pam50: pd.Series,
    ranked_meta: pd.DataFrame,
    out_path: Path,
    *,
    cols: int = 4,
    row_height: float = 4.5,
    col_width: float = 6.0,
    dpi: int = 150,
    title: Optional[str] = None,
    page_label: Optional[str] = None,
    show_pam50: bool = True,
) -> None:
    """Two-section subplot: [H→NAT | NAT→T paired (+ PAM50)].

    Left section  — H→NAT: each NAT sample's pressure minus GTEx-population median.
    Right section — NAT→T: within-person (tumor_i − nat_i) for matched pairs,
                           with PAM50 breakdown if show_pam50.
    Both share the same y-axis (pressure delta units).
    """
    import matplotlib.pyplot as plt

    NAT_COLOR    = "#EDC948"  # yellow — healthy→NAT
    PAIRED_COLOR = "#555555"  # dark grey — cohort-wide paired delta

    edges_subset = list(edges_subset)
    if not edges_subset:
        return

    ncol = max(1, min(cols, len(edges_subset)))
    nrow = int(np.ceil(len(edges_subset) / ncol))
    fig_w = col_width * ncol
    fig_h = row_height * nrow + 1.4
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, fig_h), squeeze=False)

    n_matched = len(pam50)
    n_nat = max((len(s) for s in nat_vs_healthy_map.values()), default=0)
    suptitle = title or (
        "Acquired pressure cascade: H→NAT (vs GTEx median)  |  NAT→T (within-person paired)"
    )
    if page_label:
        suptitle = f"{suptitle}  ({page_label})"
    fig.suptitle(
        f"{suptitle}\n"
        f"NAT n={n_nat}  ·  matched pairs n={n_matched}  ·  softmax_logrpm  ·  "
        "Δ>0 = pressure acquired above reference  ·  * q<0.05 Wilcoxon BH",
        fontsize=9.5, y=1.01,
    )

    # pre-build meta lookup
    meta_idx: dict = {}
    if not ranked_meta.empty:
        for _, row in ranked_meta.iterrows():
            meta_idx[(str(row["miRNA"]), str(row["gene"]))] = row.to_dict()

    rng = np.random.default_rng(0)

    for i, key in enumerate(edges_subset):
        ax = axes[i // ncol][i % ncol]
        arm, gene = key
        delta = delta_map.get(key, pd.Series(dtype=float))
        h2n   = nat_vs_healthy_map.get(key, pd.Series(dtype=float))

        if len(delta) == 0 and len(h2n) == 0:
            ax.axis("off")
            continue

        meta = meta_idx.get(key, {})

        # ------------------------------------------------------------------ #
        # Layout: pos 1 = H→NAT, then gap, pos 2.4+ = NAT→T groups
        # ------------------------------------------------------------------ #
        pos_h2n = 1.0
        gap = 1.5  # gap between sections

        # Build NAT→T groups (right section)
        nat_to_t_groups: List[Tuple[str, np.ndarray, str]] = []
        d_all = pd.to_numeric(delta, errors="coerce").dropna().values
        nat_to_t_groups.append(("All", d_all, PAIRED_COLOR))
        if show_pam50:
            for sub in PAM50_ORDER:
                sub_p = pam50.index[pam50.eq(sub)]
                sub_d = pd.to_numeric(delta.reindex(sub_p), errors="coerce").dropna()
                if len(sub_d) >= 5:
                    nat_to_t_groups.append((sub, sub_d.values, PAM50_COLORS[sub]))

        n_right = len(nat_to_t_groups)
        pos_right_start = pos_h2n + gap
        positions_right = [pos_right_start + j * 1.1 for j in range(n_right)]

        # ---- draw H→NAT box ---- #
        h2n_vals = pd.to_numeric(h2n, errors="coerce").dropna().values
        if len(h2n_vals) > 0:
            bp = ax.boxplot([h2n_vals], positions=[pos_h2n], widths=0.52,
                            patch_artist=True, showfliers=False,
                            medianprops={"color": "black", "linewidth": 1.2}, zorder=2)
            bp["boxes"][0].set_facecolor(NAT_COLOR)
            bp["boxes"][0].set_alpha(0.6)
            bp["boxes"][0].set_edgecolor("#333333")
            _overlay_dots(ax, pos_h2n, h2n_vals, NAT_COLOR, rng=rng)

        # ---- draw NAT→T boxes ---- #
        bp2 = ax.boxplot(
            [g[1] for g in nat_to_t_groups],
            positions=positions_right,
            widths=0.52, patch_artist=True, showfliers=False,
            medianprops={"color": "black", "linewidth": 1.2}, zorder=2,
        )
        for patch, (_, vals, color) in zip(bp2["boxes"], nat_to_t_groups):
            if len(vals) == 0:
                patch.set_visible(False)
            else:
                patch.set_facecolor(color)
                patch.set_alpha(0.55)
                patch.set_edgecolor("#333333")
        for pos, (_, vals, color) in zip(positions_right, nat_to_t_groups):
            _overlay_dots(ax, pos, vals, color, rng=rng)

        # ---- decorations ---- #
        sep_x = (pos_h2n + pos_right_start) / 2
        ax.axvline(sep_x, color="#cccccc", lw=0.9, linestyle=":", zorder=1)
        ax.axhline(0, color="#888888", linewidth=0.75, linestyle="--", zorder=1)

        # section labels
        ax.text(pos_h2n, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.1,
                "H→NAT", ha="center", va="bottom", fontsize=6.5,
                color="#9A7D0A", fontstyle="italic")
        ax.text(np.mean(positions_right),
                ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.1,
                "NAT→T (paired)", ha="center", va="bottom", fontsize=6.5,
                color="#555555", fontstyle="italic")

        # significance annotations
        q_h2n = meta.get("q_h_to_nat", np.nan)
        q_n2t = meta.get("q_wilcoxon", np.nan)
        med_h2n = float(np.median(h2n_vals)) if len(h2n_vals) else np.nan
        med_n2t = float(np.median(d_all)) if len(d_all) else np.nan
        ann_parts = []
        if pd.notna(med_h2n):
            ann_parts.append(f"H→N Δ={med_h2n:+.2f}{_sig_str(q_h2n)}")
        if pd.notna(med_n2t):
            ann_parts.append(f"N→T Δ={med_n2t:+.2f}{_sig_str(q_n2t)}")
        if ann_parts:
            ax.annotate(
                "  ".join(ann_parts),
                xy=(0.97, 0.02), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=6.5,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"),
            )

        # x ticks
        all_positions = [pos_h2n] + positions_right
        h2n_lab = f"H→NAT\n(n={len(h2n_vals)})"
        right_labs = [
            f"{lab}\n(n={len(v)})" for lab, v, _ in nat_to_t_groups
        ]
        ax.set_xticks(all_positions)
        ax.set_xticklabels([h2n_lab] + right_labs, fontsize=6.5)
        ax.set_xlim(pos_h2n - 0.7, positions_right[-1] + 0.7)

        ev = meta.get("evidence_score", np.nan)
        ev_str = f"  ev={ev:.0f}" if pd.notna(ev) else ""
        short_arm = arm.replace("hsa-", "")
        ax.set_title(f"{short_arm} → {gene}{ev_str}", fontsize=7.5, pad=4)

        if i % ncol == 0:
            ax.set_ylabel("Δ pressure (vs reference)")
        ax.grid(axis="y", alpha=0.22, linestyle="--")
        ax.set_axisbelow(True)
        ax.margins(y=0.12)

    for j in range(len(edges_subset), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    # legend
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=NAT_COLOR,    alpha=0.6, label="H→NAT  (vs GTEx pop. median)"),
        Patch(facecolor=PAIRED_COLOR, alpha=0.55, label="NAT→T  (within-person paired)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))

    fig.subplots_adjust(top=0.93, hspace=0.65, wspace=0.35, bottom=0.04)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(
    *,
    out_dir: Path = OUT_DIR,
    out_name: str = "edge_acquired_pressure_panel",
    cols: int = 4,
    per_page: Optional[int] = None,
    row_height: float = 4.2,
    col_width: float = 5.5,
    max_edges: Optional[int] = 500,
    high_evidence_only: bool = True,
    hallmark: Optional[str] = None,
    show_pam50: bool = True,
    title: Optional[str] = None,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    edges = _load_edges(high_evidence_only=high_evidence_only, hallmark=hallmark)
    if edges.empty:
        print("[edge_acquired_pressure] no edges after filtering — exiting")
        return {}
    genes = sorted(edges["gene"].unique())
    print(f"[edge_acquired_pressure] edges={len(edges)}  genes={len(genes)}"
          f"  high_ev={high_evidence_only}  hallmark={hallmark}")

    print("[edge_acquired_pressure] loading matched-pair miRNA matrices …")
    tumor_mat, nat_mat, matched = _matched_pair_matrices()
    print(f"[edge_acquired_pressure] matched participants: {len(matched)}"
          f"  tumor_mat={tumor_mat.shape}  nat_mat={nat_mat.shape}")

    print("[edge_acquired_pressure] loading full NAT + GTEx matrices …")
    tissue = _state_matrices()
    nat_full_mat  = tissue["nat"]
    gtex_mat      = tissue.get("gtex", pd.DataFrame())

    print("[edge_acquired_pressure] computing pressure maps …")
    tumor_map    = _state_pressure_map(edges, tumor_mat, genes)
    nat_mat_map  = _state_pressure_map(edges, nat_mat, genes)  # matched NAT (paired)
    nat_full_map = _state_pressure_map(edges, nat_full_mat, genes)  # all NAT (for H→NAT)
    gtex_map     = _state_pressure_map(edges, gtex_mat, genes) if not gtex_mat.empty else {}
    print(f"  tumor={len(tumor_map)}  nat(matched)={len(nat_mat_map)}"
          f"  nat(full)={len(nat_full_map)}  gtex={len(gtex_map)}")

    # Paired within-person delta: Δ = tumor_i − nat_i (matched only)
    delta_map = _paired_delta(tumor_map, nat_mat_map)
    print(f"[edge_acquired_pressure] paired Δ edges: {len(delta_map)}")

    # H→NAT: each NAT sample − GTEx population median
    nat_vs_healthy_map: Dict[tuple, pd.Series] = {}
    for k, nat_s in nat_full_map.items():
        gtex_s = pd.to_numeric(gtex_map.get(k, pd.Series(dtype=float)), errors="coerce")
        nat_vals = pd.to_numeric(nat_s, errors="coerce")
        if len(nat_vals) >= 5 and len(gtex_s) >= 5:
            healthy_med = float(gtex_s.median())
            nat_vs_healthy_map[k] = nat_vals - healthy_med

    ranked = _score_and_rank(delta_map)
    if max_edges is not None:
        ranked = ranked[:max_edges]

    # Build summary with both Wilcoxon tests (NAT→T paired and H→NAT)
    pam50 = _pam50_series(matched)
    summary_rows: List[dict] = []
    p_vals_n2t: List[float] = []
    p_vals_h2n: List[float] = []
    for k in ranked:
        arm, gene = k
        d = delta_map[k]
        h = nat_vs_healthy_map.get(k, pd.Series(dtype=float))
        ev_row = edges.loc[edges["miRNA"].eq(arm) & edges["gene"].eq(gene)]
        ev = float(ev_row["evidence_score"].iloc[0]) if not ev_row.empty else np.nan
        p_n2t = _wilcoxon_p(d)
        p_h2n = _wilcoxon_p(h)
        p_vals_n2t.append(p_n2t)
        p_vals_h2n.append(p_h2n)
        summary_rows.append({
            "miRNA": arm, "gene": gene, "evidence_score": ev,
            "n_pairs": len(d),
            "median_nat_to_t": float(d.median()) if len(d) else np.nan,
            "q25_nat_to_t": float(d.quantile(0.25)) if len(d) else np.nan,
            "q75_nat_to_t": float(d.quantile(0.75)) if len(d) else np.nan,
            "frac_positive_nat_to_t": float((d > 0).mean()) if len(d) else np.nan,
            "p_wilcoxon": p_n2t,
            "n_h_to_nat": len(h),
            "median_h_to_nat": float(h.median()) if len(h) else np.nan,
            "p_h_to_nat": p_h2n,
        })
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        v_n2t = summary["p_wilcoxon"].notna()
        if v_n2t.any():
            summary.loc[v_n2t, "q_wilcoxon"] = bh_fdr(summary.loc[v_n2t, "p_wilcoxon"].values)
        v_h2n = summary["p_h_to_nat"].notna()
        if v_h2n.any():
            summary.loc[v_h2n, "q_h_to_nat"] = bh_fdr(summary.loc[v_h2n, "p_h_to_nat"].values)
    summary.to_csv(out_dir / f"{out_name}_summary.tsv", sep="\t", index=False)
    for col, label in [("q_wilcoxon", "NAT→T"), ("q_h_to_nat", "H→NAT")]:
        if col in summary:
            n_sig = int((summary[col].fillna(1) < 0.05).sum())
            print(f"[edge_acquired_pressure] q<0.05 {label}: {n_sig}")

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
        plot_acquired_pressure_panels(
            chunk, delta_map, nat_vs_healthy_map, pam50, summary,
            png,
            cols=cols, row_height=row_height, col_width=col_width,
            title=title, page_label=page_label, show_pam50=show_pam50,
        )
        figure_paths.append(str(png.relative_to(C.OUTPUT_ROOT)))
        print(f"[edge_acquired_pressure] page {page_i}/{n_pages} ({len(chunk)} edges) -> {png}")

    manifest = {
        "module": "mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_matched_participants": len(matched),
        "n_nat_full": nat_full_mat.shape[1],
        "n_gtex": gtex_mat.shape[1] if not gtex_mat.empty else 0,
        "n_edges_total": len(edges),
        "n_edges_with_paired_delta": len(delta_map),
        "n_edges_with_nat_vs_healthy": len(nat_vs_healthy_map),
        "n_edges_plotted": len(ranked),
        "n_pages": n_pages,
        "per_page": per_page,
        "attr_mode": ATTR_MODE,
        "target_norm": TARGET_NORM,
        "high_evidence_only": high_evidence_only,
        "hallmark_filter": hallmark,
        "sections": {
            "H→NAT": "pressure_nat(s) − median(pressure_healthy)  [all NAT, population median]",
            "NAT→T": "pressure_tumor(p) − pressure_nat(p)  [matched within-person]",
        },
        "ranking": "median NAT→T paired Δ descending",
        "stat_test": "Wilcoxon signed-rank (one-sample H0: Δ=0), BH FDR",
        "figures": figure_paths,
        "summary": f"{out_name}_summary.tsv",
    }
    (out_dir / f"{out_name}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"[edge_acquired_pressure] done — {len(ranked)} edges, {n_pages} pages")
    return {"delta_map": delta_map, "summary": summary, "manifest": manifest}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-edges", type=int, default=500)
    ap.add_argument("--all-edges", action="store_true")
    ap.add_argument("--high-evidence-only", action="store_true", default=True)
    ap.add_argument("--no-high-evidence-filter", dest="high_evidence_only",
                    action="store_false")
    ap.add_argument("--hallmark", type=str, default=None)
    ap.add_argument("--cols", type=int, default=4)
    ap.add_argument("--per-page", type=int, default=None)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--out-name", type=str, default="edge_acquired_pressure_panel")
    ap.add_argument("--row-height", type=float, default=4.2)
    ap.add_argument("--col-width", type=float, default=5.5)
    ap.add_argument("--no-pam50", dest="show_pam50", action="store_false",
                    help="show cohort-only boxplot without PAM50 split")
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
        show_pam50=args.show_pam50,
        title=args.title,
    )


if __name__ == "__main__":
    main()
