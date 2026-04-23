#!/usr/bin/env python3
"""
Stratify omics module coverage by PAM50 subtype and by collapsed pathologic stage,
using `BRCA_clinical_immune_unified.tsv` and a `participant_presence_omics.tsv` from
`sample_module_coverage.py`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.sample_coverage.coverage_module_lists import OMICS_MODULES  # noqa: E402
from pipeline.config import PATHS  # noqa: E402


def _strat_table(merged: pd.DataFrame, group_col: str, mods: List[str]) -> pd.DataFrame:
    rows = []
    for g, sub in merged.groupby(group_col, dropna=False):
        gkey = str(g).strip() if g is not None and not (isinstance(g, float) and np.isnan(g)) else ""
        if not gkey or gkey.lower() in ("nan", "none", ""):
            continue
        n = len(sub)
        row: dict = {group_col: gkey, "n_participants": n}
        for m in mods:
            if m not in sub.columns:
                row[f"pct_{m}"] = 0.0
            else:
                row[f"pct_{m}"] = round(100.0 * float(sub[m].astype(bool).sum()) / n, 2) if n else 0.0
        hit = sub[mods].astype(bool).sum(axis=1) if all(c in sub.columns for c in mods) else pd.Series([0] * n)
        row["mean_n_omics_modules"] = round(float(hit.mean()), 3) if n else 0.0
        rows.append(row)
    out = pd.DataFrame(rows).sort_values("n_participants", ascending=False)
    return out


def _plot_prevalence_heatmap(
    table: pd.DataFrame,
    group_col: str,
    mods: List[str],
    out_png: Path,
    title: str,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"Skipping figure {out_png.name}: {exc}")
        return

    t = table.set_index(group_col)[["pct_" + m for m in mods if "pct_" + m in table.columns]]
    if t.empty:
        return
    t.columns = [c.replace("pct_", "") for c in t.columns]
    arr = t.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * t.shape[1]), max(4, 0.35 * t.shape[0])))
    im = ax.imshow(arr, aspect="auto", vmin=0, vmax=100, cmap="YlOrRd")
    ax.set_xticks(range(t.shape[1]))
    ax.set_xticklabels(list(t.columns), rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(t.shape[0]))
    ax.set_yticklabels(list(t.index), fontsize=7)
    ncell = max(arr.shape[0], arr.shape[1])
    fs = max(4.0, min(8.0, 180.0 / max(ncell, 1)))
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = float(arr[i, j])
            s = f"{v:.0f}" if abs(v - round(v)) < 0.05 else f"{v:.1f}"
            ax.text(j, i, s, ha="center", va="center", color="white" if v > 55 else "black", fontsize=fs)
    ax.set_title(title)
    ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="% participants with module")
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def run_stratification(
    *,
    presence_path: Path,
    out_tables_dir: Path,
    out_figures_dir: Path,
    clinical_path: Path | None = None,
) -> None:
    """
    Programmatic entrypoint used by `analysis/sample_module_coverage.py`.
    Writes stratification TSVs under out_tables_dir and heatmap figures under out_figures_dir.
    """
    clin_path = clinical_path or PATHS.brca_clinical_immune_unified
    clin = pd.read_csv(clin_path, sep="\t", low_memory=False)
    if "participant" not in clin.columns:
        raise ValueError(f"Missing participant column in {clin_path}")

    stage_col = "pathologic_stage_collapsed"
    if stage_col not in clin.columns:
        raise ValueError(
            f"Missing {stage_col}; run scripts/annotations/update_brca_clinical_unified_stage.py"
        )

    clin = clin.drop_duplicates(subset=["participant"], keep="first")

    pres = pd.read_csv(presence_path, sep="\t", low_memory=False)
    idc = pres.columns[0]
    pres = pres.rename(columns={idc: "participant"})

    merged = clin.merge(pres, on="participant", how="inner")
    print(f"Merged clinical + omics presence: {len(merged)} participants")

    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    pam_col = "PAM50_final"
    if pam_col not in merged.columns:
        raise ValueError(f"Missing {pam_col} in unified clinical")

    t_pam = _strat_table(merged, pam_col, OMICS_MODULES)
    t_pam.to_csv(out_tables_dir / "stratification_by_PAM50.tsv", sep="\t", index=False)

    stage_sub = merged.loc[merged[stage_col].astype(str).str.len() > 0].copy()
    t_stg = _strat_table(stage_sub, stage_col, OMICS_MODULES)
    t_stg.to_csv(
        out_tables_dir / "stratification_by_pathologic_stage_collapsed.tsv",
        sep="\t",
        index=False,
    )

    _plot_prevalence_heatmap(
        t_pam,
        pam_col,
        OMICS_MODULES,
        out_figures_dir / "fig_omics_prevalence_by_PAM50.png",
        "Omics module prevalence by PAM50 (participants in clinical ∩ presence)",
    )
    if not t_stg.empty:
        _plot_prevalence_heatmap(
            t_stg,
            stage_col,
            OMICS_MODULES,
            out_figures_dir / "fig_omics_prevalence_by_stage_collapsed.png",
            "Omics module prevalence by collapsed stage",
        )

    print(f"Wrote stratification tables under {out_tables_dir}")
    print(f"Wrote stratification figures under {out_figures_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description="PAM50 / stage stratification of omics coverage.")
    ap.add_argument(
        "--presence",
        type=Path,
        required=True,
        help="participant_presence_omics.tsv from a coverage run",
    )
    ap.add_argument(
        "--clinical",
        type=Path,
        default=None,
        help="Unified clinical TSV (default: PATHS.brca_clinical_immune_unified)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: sibling folder stratification_<stamp> next to presence file)",
    )
    args = ap.parse_args()

    out_dir = args.out
    if out_dir is None:
        out_dir = args.presence.parent / f"stratification_{pd.Timestamp.utcnow().strftime('%Y%m%d_%H%M%S')}"
    run_stratification(
        presence_path=args.presence,
        out_tables_dir=out_dir,
        out_figures_dir=out_dir,
        clinical_path=args.clinical,
    )
    print(f"Wrote stratification tables and figures under {out_dir}")


if __name__ == "__main__":
    main()
