"""Grid boxplot panels: per-miRNA arm expression across tissue states + PAM50.

Each miRNA subplot shows boxplots with overlaid sample dots for:

  - GTEx v10 healthy breast
  - TCGA tumor-adjacent normal (NAT, 11)
  - TCGA primary tumor (01) — all tumors (default) or **NAT-matched only**
  - TCGA primary tumor by PAM50 (LumA / LumB / Her2 / Basal; Normal-like excluded)

All samples in the expression vector are plotted. Explicit zeros and missing
values (NaN in the GDC matrix, treated as undetected → 0 on log2(PM+1)) are
included so cohort n is always the full group size when the arm row exists.
Empty only when the arm is absent from that state's matrix.

Large arm lists are split across multiple figure pages (``--per-page``).

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_expression_panels --preset all
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_expression_panels --preset all --tumor-scope matched --out-name all_arms_matched_panel
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_expression_panels --preset all_hsa --cols 4 --per-page 48
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
from mirna_hallmark.analyses.cross_state.cross_state_coupling import FOCUS_ARMS
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
from mirna_hallmark.data_loaders import load_clinical_strata
from mirna_hallmark.family_normal_reference import _load_full_mirna, _participant, _split_types
from mirna_hallmark.analyses.mir301.mir301_family_depth import FAMILY_ARMS

OUT_DIR = C.TISSUE_REFERENCE_DIR / "cross_state_expression_panels"

TumorScope = Literal["all", "matched"]

PAM50_ORDER = list(C.SUBTYPE_SINGLETONS)  # LumA, LumB, Her2, Basal
PAM50_COLORS = {
    "LumA": "#4C72B0",
    "LumB": "#8172B2",
    "Her2": "#C44E52",
    "Basal": "#55A868",
}
STATE_COLORS = {"gtex": "#59A14F", "nat": "#EDC948", "tumor": "#E15759"}

PRESETS: Dict[str, Sequence[str]] = {
    "coupling_focus": FOCUS_ARMS,
    "family": FAMILY_ARMS,
}


class PlotGroup(NamedTuple):
    key: str
    label: str
    color: str
    kind: str  # tissue | pam50


def _tumor_by_participant(tumor: pd.DataFrame) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for c in tumor.columns:
        out.setdefault(_participant(str(c)), []).append(str(c))
    return out


def _matched_tumor_columns(nat: pd.DataFrame, tumor: pd.DataFrame) -> List[str]:
    """One primary tumor per NAT participant (same rule as ``nat_tumor_umap`` matched)."""
    nat_by_p = {_participant(str(c)): str(c) for c in nat.columns}
    tumor_by_p = _tumor_by_participant(tumor)
    cols: List[str] = []
    for p in sorted(nat_by_p):
        if p not in tumor_by_p:
            continue
        cols.append(sorted(tumor_by_p[p])[0])
    return cols


def _filter_tumor_matrix(
    tumor: pd.DataFrame,
    nat: pd.DataFrame,
    *,
    tumor_scope: TumorScope,
) -> pd.DataFrame:
    if tumor_scope == "all":
        return tumor
    cols = _matched_tumor_columns(nat, tumor)
    if not cols:
        raise ValueError("No NAT-matched primary tumors found in the miRNA matrix")
    return tumor[cols]


def _all_cohort_arms(*, hsa_only: bool = False) -> List[str]:
    """Every arm row in the primary-tumor matrix, sorted by full-cohort median (desc)."""
    full = _load_full_mirna()
    tumor = _split_types(full)["tumor"]
    med = pd.to_numeric(tumor.median(axis=1), errors="coerce")
    order = med.sort_values(ascending=False).index.astype(str).tolist()
    if hsa_only:
        order = [a for a in order if a.startswith("hsa-")]
    return order


def _load_arm_order_from_summary(path: Path) -> List[str]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    col = "arm" if "arm" in df.columns else ("miRNA" if "miRNA" in df.columns else df.columns[0])
    seen: set[str] = set()
    out: List[str] = []
    for a in df[col].astype(str):
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def _apply_arm_order(arms: Sequence[str], order_file: Path) -> List[str]:
    """Reorder *arms* to match a prior panel summary (page-aligned comparison)."""
    ref = _load_arm_order_from_summary(order_file)
    arm_set = set(arms)
    ordered = [a for a in ref if a in arm_set]
    trailing = [a for a in arms if a not in set(ref)]
    return ordered + trailing


def _default_arm_order_file(out_dir: Path, preset: Optional[str]) -> Optional[Path]:
    """Map preset → whole-cohort summary TSV used to pin paging for matched reruns."""
    names = {
        "all": "all_arms_panel_summary.tsv",
        "all_hsa": "all_hsa_panel_summary.tsv",
        "landscape_extended": "landscape_extended_panel_summary.tsv",
        "coupling_focus": "coupling_focus_panel_summary.tsv",
        "family": "family_panel_summary.tsv",
    }
    name = names.get(preset or "")
    if not name:
        return None
    cand = out_dir / name
    return cand if cand.is_file() else None


def _landscape_extended_arms(*, n_induced: int = 30, n_movers: int = 25) -> List[str]:
    """Union of coupling focus, top tumor-induced DE arms, and pressure movers."""
    arms: List[str] = list(FOCUS_ARMS)
    de_path = C.TISSUE_REFERENCE_DIR / "cross_state_landscape" / "mirna_arm_de_landscape.tsv"
    if de_path.is_file():
        de = pd.read_csv(de_path, sep="\t", low_memory=False)
        sub = de.loc[de["de_class"].eq("tumor_induced") & de["arm"].astype(str).str.startswith("hsa-")]
        arms.extend(sub.head(n_induced)["arm"].astype(str).tolist())
    share_path = C.TISSUE_REFERENCE_DIR / "cross_state_landscape" / "mirna_pressure_share_shift.tsv"
    if share_path.is_file():
        sh = pd.read_csv(share_path, sep="\t", low_memory=False)
        col = "miRNA" if "miRNA" in sh.columns else sh.columns[0]
        arms.extend(
            sh.loc[sh[col].astype(str).str.startswith("hsa-")].head(n_movers)[col].astype(str).tolist()
        )
    seen: set[str] = set()
    out: List[str] = []
    for a in arms:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def _parse_arms(
    arms: Optional[str],
    arms_file: Optional[Path],
    preset: Optional[str],
    *,
    max_arms: Optional[int] = None,
    arm_order_file: Optional[Path] = None,
) -> List[str]:
    out: List[str] = []
    if preset == "landscape_extended":
        out.extend(_landscape_extended_arms())
    elif preset == "all":
        out.extend(_all_cohort_arms(hsa_only=False))
    elif preset == "all_hsa":
        out.extend(_all_cohort_arms(hsa_only=True))
    elif preset:
        if preset not in PRESETS:
            raise ValueError(
                f"Unknown preset {preset!r}; choose from {sorted(PRESETS) + ['landscape_extended', 'all', 'all_hsa']}"
            )
        out.extend(PRESETS[preset])
    if arms_file:
        df = pd.read_csv(arms_file, sep="\t", low_memory=False)
        col = "arm" if "arm" in df.columns else ("miRNA" if "miRNA" in df.columns else df.columns[0])
        out.extend(df[col].astype(str).tolist())
    if arms:
        out.extend([a.strip() for a in arms.split(",") if a.strip()])
    if not out:
        out = _landscape_extended_arms()
    seen: set[str] = set()
    deduped: List[str] = []
    for a in out:
        if a not in seen:
            seen.add(a)
            deduped.append(a)
    if max_arms is not None:
        deduped = deduped[: max_arms]
    if arm_order_file is not None:
        deduped = _apply_arm_order(deduped, arm_order_file)
    return deduped


def _arm_series(mat: Optional[pd.DataFrame], arm: str) -> Tuple[np.ndarray, int, int]:
    """Per-sample values for an arm when present in the matrix.

    Returns (values, n_imputed_zero, n_cohort). NaN → 0 on log2(PM+1) scale so
    undetected / missing measurements stay in the cohort count and plot at y=0.
    """
    if mat is None or mat.empty or arm not in mat.index:
        return np.array([], dtype=float), 0, 0
    raw = pd.to_numeric(mat.loc[arm], errors="coerce")
    n_cohort = int(mat.shape[1])
    n_imputed = int(raw.isna().sum())
    vals = raw.fillna(0.0).to_numpy(dtype=float)
    return vals, n_imputed, n_cohort


def _arm_values(mat: Optional[pd.DataFrame], arm: str) -> np.ndarray:
    vals, _, _ = _arm_series(mat, arm)
    return vals


def _pam50_tumor_matrices(tumor: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Split primary-tumor columns by PAM50_final (Normal-like excluded)."""
    clin = load_clinical_strata()
    pam = clin.set_index("participant")["PAM50_final"]
    cols_by_sub: Dict[str, List[str]] = {s: [] for s in PAM50_ORDER}
    for col in tumor.columns:
        sub = pam.get(_participant(str(col)))
        if sub in cols_by_sub:
            cols_by_sub[sub].append(col)
    return {s: tumor[cols] for s, cols in cols_by_sub.items() if cols}


def build_panel_context(*, tumor_scope: TumorScope = "all") -> Tuple[Dict[str, pd.DataFrame], List[PlotGroup]]:
    tissue = _state_matrices()
    tumor = _filter_tumor_matrix(tissue["tumor"], tissue["nat"], tumor_scope=tumor_scope)
    tissue = {**tissue, "tumor": tumor}
    pam50 = _pam50_tumor_matrices(tumor)
    states: Dict[str, pd.DataFrame] = {**tissue, **pam50}

    tumor_label = "Tumor (matched)" if tumor_scope == "matched" else "Tumor"
    groups: List[PlotGroup] = [
        PlotGroup("gtex", "GTEx", STATE_COLORS["gtex"], "tissue"),
        PlotGroup("nat", "NAT", STATE_COLORS["nat"], "tissue"),
        PlotGroup("tumor", tumor_label, STATE_COLORS["tumor"], "tissue"),
    ]
    for sub in PAM50_ORDER:
        if sub in states:
            groups.append(PlotGroup(sub, sub, PAM50_COLORS[sub], "pam50"))
    return states, groups


def build_expression_summary(
    arms: Sequence[str],
    states: Dict[str, pd.DataFrame],
    groups: Sequence[PlotGroup],
) -> pd.DataFrame:
    rows: List[dict] = []
    for arm in arms:
        for g in groups:
            mat = states.get(g.key)
            vals, n_imputed, n_cohort = _arm_series(mat, arm)
            in_vector = mat is not None and not mat.empty and arm in mat.index
            rows.append(
                {
                    "arm": arm,
                    "group": g.key,
                    "group_kind": g.kind,
                    "in_expression_vector": bool(in_vector),
                    "n_cohort": n_cohort,
                    "n_plotted": int(len(vals)),
                    "n_imputed_zero": n_imputed,
                    "median_log2pm": float(np.median(vals)) if len(vals) else np.nan,
                    "has_data": bool(in_vector and len(vals) > 0),
                }
            )
    return pd.DataFrame(rows)


def _global_banner(states: Dict[str, pd.DataFrame], groups: Sequence[PlotGroup]) -> str:
    parts: List[str] = []
    for g in groups:
        mat = states.get(g.key)
        if mat is not None and not mat.empty:
            parts.append(f"{g.label} n={mat.shape[1]}")
    return "  |  ".join(parts)


def _overlay_dots(
    ax,
    position: float,
    vals: np.ndarray,
    color: str,
    *,
    rng: np.random.Generator,
    max_dots: int = 1200,
) -> None:
    if len(vals) == 0:
        return
    plot_vals = vals
    if len(vals) > max_dots:
        idx = rng.choice(len(vals), size=max_dots, replace=False)
        plot_vals = vals[idx]
    jitter = rng.uniform(-0.14, 0.14, size=len(plot_vals))
    ax.scatter(
        position + jitter,
        plot_vals,
        s=5,
        alpha=0.22,
        color=color,
        edgecolors="none",
        linewidths=0,
        zorder=3,
    )


def _paginate(arms: Sequence[str], per_page: int) -> List[List[str]]:
    per_page = max(1, per_page)
    return [list(arms[i : i + per_page]) for i in range(0, len(arms), per_page)]


def plot_expression_panels(
    arms: Sequence[str],
    states: Dict[str, pd.DataFrame],
    groups: Sequence[PlotGroup],
    out_path: Path,
    *,
    cols: int = 4,
    row_height: float = 4.0,
    col_width: float = 5.8,
    dpi: int = 150,
    title: Optional[str] = None,
    page_label: Optional[str] = None,
) -> None:
    import matplotlib.pyplot as plt

    arms = list(arms)
    if not arms:
        return

    active_groups = [g for g in groups if g.key in states and not states[g.key].empty]
    ncol = max(1, min(cols, len(arms)))
    nrow = int(np.ceil(len(arms) / ncol))

    fig_w = col_width * ncol
    fig_h = row_height * nrow + 1.4
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, fig_h), squeeze=False)

    banner = _global_banner(states, active_groups)
    suptitle = title or "miRNA arm expression: tissue reference + PAM50 subtypes"
    if page_label:
        suptitle = f"{suptitle}  ({page_label})"
    fig.suptitle(
        f"{suptitle}\n{banner}  ·  log2(PM+1); dots = samples; GTEx TPM vs TCGA RPM",
        fontsize=10,
        y=1.01,
    )

    rng = np.random.default_rng(0)
    tissue_n = sum(1 for g in active_groups if g.kind == "tissue")

    for i, arm in enumerate(arms):
        ax = axes[i // ncol][i % ncol]
        tick_labels: List[str] = []
        plot_groups: List[np.ndarray] = []
        colors: List[str] = []
        positions: List[float] = []

        pos = 1.0
        for g in active_groups:
            vals, _, n_cohort = _arm_series(states[g.key], arm)
            in_vector = arm in states[g.key].index
            n_label = n_cohort if in_vector else 0
            tick_labels.append(f"{g.label}\n(n={n_label})")
            plot_groups.append(vals)
            colors.append(g.color)
            positions.append(pos)
            pos += 1.0
            if g.kind == "tissue" and g.key == "tumor" and tissue_n == 3:
                pos += 0.45  # gap before PAM50 block

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

        for pos, vals, color in zip(positions, plot_groups, colors):
            if len(vals) == 0:
                ax.text(pos, 0.02, "—", ha="center", va="bottom", fontsize=10, color="#999999")
            else:
                _overlay_dots(ax, pos, vals, color, rng=rng)

        if tissue_n == 3 and len(positions) >= 4:
            sep = (positions[2] + positions[3]) / 2.0
            ax.axvline(sep, color="#bbbbbb", linestyle=":", linewidth=0.9, zorder=1)

        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels, fontsize=7)
        ax.set_title(arm, fontsize=8, pad=5)
        if i % ncol == 0:
            ax.set_ylabel("log2(PM+1)")
        ax.grid(axis="y", alpha=0.25, linestyle="--")
        ax.set_axisbelow(True)
        ax.margins(y=0.08)

    for j in range(len(arms), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    fig.subplots_adjust(top=0.93, hspace=0.55, wspace=0.32)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def run(
    *,
    arms: Sequence[str],
    out_dir: Path = OUT_DIR,
    out_name: str = "expression_panel",
    cols: int = 4,
    per_page: Optional[int] = None,
    row_height: float = 4.0,
    col_width: float = 5.8,
    title: Optional[str] = None,
    tumor_scope: TumorScope = "all",
    arm_order_file: Optional[Path] = None,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    arms = list(arms)
    if per_page is None:
        per_page = cols * 12
    pages = _paginate(arms, per_page)
    n_pages = len(pages)

    states, groups = build_panel_context(tumor_scope=tumor_scope)
    summary = build_expression_summary(arms, states, groups)
    summary.to_csv(out_dir / f"{out_name}_summary.tsv", sep="\t", index=False)

    figure_paths: List[str] = []
    for page_i, chunk in enumerate(pages, start=1):
        if n_pages == 1:
            stem = out_name
        else:
            stem = f"{out_name}_page{page_i:02d}"
        png = fig_dir / f"{stem}.png"
        page_label = f"page {page_i}/{n_pages}, arms {(page_i - 1) * per_page + 1}–{(page_i - 1) * per_page + len(chunk)}"
        plot_expression_panels(
            chunk,
            states,
            groups,
            png,
            cols=cols,
            row_height=row_height,
            col_width=col_width,
            title=title,
            page_label=page_label,
        )
        figure_paths.append(str(png.relative_to(C.OUTPUT_ROOT)))
        print(f"[expr_panels] page {page_i}/{n_pages} ({len(chunk)} arms) -> {png}")

    state_counts = {g.key: int(states[g.key].shape[1]) for g in groups if g.key in states}
    manifest = {
        "module": "mirna_hallmark.analyses.cross_state.cross_state_expression_panels",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_arms": len(arms),
        "n_pages": n_pages,
        "per_page": per_page,
        "cols": cols,
        "group_counts": state_counts,
        "pam50_subtypes": PAM50_ORDER,
        "pam50_excluded": list(C.SUBTYPE_EXCLUDE_FROM_CONTRAST),
        "tumor_scope": tumor_scope,
        "arm_order_file": str(arm_order_file) if arm_order_file else None,
        "figures": figure_paths,
        "summary": f"{out_name}_summary.tsv",
        "expression_rule": "arm in matrix → all cohort samples plotted; NaN imputed as 0 on log2(PM+1); empty only if arm absent from matrix",
        "cross_platform_note": "GTEx TPM/URS vs TCGA RPM/MIMAT — shape/rank only for GTEx vs TCGA",
    }
    (out_dir / f"{out_name}_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(f"[expr_panels] groups: {state_counts}; arms={len(arms)}; pages={n_pages}; grid={cols} cols")
    missing = summary.loc[~summary["has_data"] & summary["in_expression_vector"], ["arm", "group"]]
    if not missing.empty:
        print(f"[expr_panels] in vector but all-NaN: {len(missing)} arm×group cells")
    return summary, state_counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", type=str, default=None, help="comma-separated arm names")
    ap.add_argument("--arms-file", type=Path, default=None, help="TSV with arm/miRNA column")
    ap.add_argument(
        "--preset",
        choices=sorted(PRESETS) + ["landscape_extended", "all", "all_hsa"],
        default="landscape_extended",
        help="arm list preset; all = every tumor-matrix row; all_hsa = hsa- arms only",
    )
    ap.add_argument("--max-arms", type=int, default=None, help="cap number of arms plotted")
    ap.add_argument("--per-page", type=int, default=None, help="arms per figure page (default: cols×12)")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--out-name", type=str, default="expression_panel")
    ap.add_argument("--cols", type=int, default=4, help="grid columns")
    ap.add_argument("--row-height", type=float, default=4.0, help="inches per subplot row")
    ap.add_argument("--col-width", type=float, default=5.8, help="inches per subplot column")
    ap.add_argument(
        "--tumor-scope",
        choices=["all", "matched"],
        default="all",
        help="primary tumor cohort: all TCGA 01 samples (default) or NAT-matched only (~103)",
    )
    ap.add_argument(
        "--arm-order-file",
        type=Path,
        default=None,
        help="pin arm list + paging to a prior panel summary TSV (default for --tumor-scope matched: preset's whole-cohort summary in --out-dir)",
    )
    ap.add_argument("--title", type=str, default=None)
    args = ap.parse_args()

    C.ensure_output_dirs()
    arm_order_file = args.arm_order_file
    if arm_order_file is None and args.tumor_scope == "matched":
        arm_order_file = _default_arm_order_file(args.out_dir, args.preset)
    arm_list = _parse_arms(
        args.arms,
        args.arms_file,
        args.preset,
        max_arms=args.max_arms,
        arm_order_file=arm_order_file,
    )
    run(
        arms=arm_list,
        out_dir=args.out_dir,
        out_name=args.out_name,
        cols=args.cols,
        per_page=args.per_page,
        row_height=args.row_height,
        col_width=args.col_width,
        title=args.title,
        tumor_scope=args.tumor_scope,
        arm_order_file=arm_order_file,
    )


if __name__ == "__main__":
    main()
