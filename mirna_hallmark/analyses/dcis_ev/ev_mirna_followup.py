"""Follow-up EV analyses: GSE255660 cancer vs benign; GSE241785 stage × TCGA pressure.

Outputs under ``mirna_hallmark/output/plasma_ev/multi_cohort/``:
  ``gse255660_cancer_vs_benign_arm_contrasts.tsv``
  ``gse255660_cancer_vs_benign_focus.tsv``
  ``gse241785_stage_ev_medians.tsv``
  ``gse241785_stage_vs_tcga_pressure.tsv``
  ``figures/gse255660_cancer_vs_benign_focus.png``
  ``figures/gse241785_stage_ev_vs_tcga_pressure.png``

Run (after multi_cohort ingest):
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.ev_mirna_followup
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.ev_mirna_multi_cohort import (
    _read_log2_wide,
    run_cohort_contrast,
)
from mirna_hallmark.analyses.dcis_ev.ev_mirna_replication import REPLICATION_FOCUS_ARMS, _anticorr_summary_for_arms
from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import MIR29_ARMS
from pipeline.config import PATHS

STAGE_ORDER = ("I", "II", "III", "IV")
TCGA_STAGE_BRIDGE = C.PLASMA_EV_DIR / "tcga_pressure_bridge" / "ev_candidate_pressure_by_stage_tn.tsv"


def _prep_255660_cancer_benign(samples: pd.DataFrame) -> pd.DataFrame:
    s = samples.loc[~samples["is_dilution_replicate"].astype(bool)].copy()
    s = s[s["bc_case"].astype(bool) | s["bc_benign"].astype(bool)]
    s["cancer_vs_benign"] = s["bc_case"].astype(bool)
    return s


def run_gse255660_cancer_vs_benign(
    *,
    out_dir: Path,
    screening_delta: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Cancer vs benign on abnormal-mammogram serum cohort (excl. healthy, dilution reps)."""
    from mirna_hallmark.analyses.dcis_ev.ev_mirna_multi_cohort import CohortSpec

    spec = CohortSpec(
        "GSE255660_cancer_benign",
        "breast_cancer vs breast_benign (serum; excl. healthy + dilution reps)",
        PATHS.mirna_ev_gse255660_log2_expr,
        PATHS.mirna_ev_gse255660_samples,
        _prep_255660_cancer_benign,
        "cancer_vs_benign",
        "cancer",
        "benign",
        "serum",
    )

    harmonized, samples, _ = run_cohort_contrast(spec)
    harmonized.to_csv(out_dir / "gse255660_cancer_vs_benign_arm_contrasts.tsv", sep="\t", index=False)

    arms = [a for a in REPLICATION_FOCUS_ARMS if a in set(harmonized["arm"])]
    focus = harmonized[harmonized["arm"].isin(arms)].copy()
    focus = focus.rename(
        columns={
            "delta": "gse255660_cancer_benign_delta",
            "q": "gse255660_cancer_benign_q",
            "auc": "gse255660_cancer_benign_auc",
            "sig_fdr05": "gse255660_cancer_benign_sig",
        }
    )

    if screening_delta is None:
        p270 = C.PLASMA_EV_DIR / "gse270497_screening_arm_contrasts.tsv"
        if p270.is_file():
            screening_delta = pd.read_csv(p270, sep="\t")
        else:
            screening_delta = pd.DataFrame()

    if not screening_delta.empty:
        c270 = screening_delta.set_index("arm")
        dcol = "delta_median_log2tpm_pos_minus_neg"
        if dcol not in c270.columns:
            dcol = "delta_median_log2_case_minus_control"
        focus["gse270497_screening_delta"] = focus["arm"].map(c270[dcol])
        focus["same_direction_as_gse270497_screening"] = np.sign(focus["gse270497_screening_delta"]) == np.sign(
            focus["gse255660_cancer_benign_delta"]
        )

    anticorr = _anticorr_summary_for_arms(arms)
    if not anticorr.empty:
        focus = focus.merge(
            anticorr[
                [
                    "arm",
                    "n_anticorr_gene_x_pam50",
                    "n_ecm_anticorr_hits",
                    "median_abs_partial_rho",
                ]
            ],
            on="arm",
            how="left",
        )

    focus.to_csv(out_dir / "gse255660_cancer_vs_benign_focus.tsv", sep="\t", index=False)

    n_case = int(samples["cancer_vs_benign"].astype(bool).sum())
    n_benign = int((~samples["cancer_vs_benign"].astype(bool)).sum())
    n_sig = int(focus["gse255660_cancer_benign_sig"].fillna(False).sum()) if "gse255660_cancer_benign_sig" in focus.columns else 0
    print(
        f"[ev_followup] GSE255660 cancer vs benign: {n_case} cancer / {n_benign} benign; "
        f"{n_sig}/{len(focus)} focus arms FDR<0.05"
    )
    return harmonized, focus


def _load_tcga_stage_pressure(path: Path = TCGA_STAGE_BRIDGE) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}; run tcga pressure bridge first")
    df = pd.read_csv(path, sep="\t")
    return df[df["stratum_type"] == "pathologic_stage_collapsed"].copy()


def run_gse241785_stage_pressure(*, out_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per-stage vesicular EV medians × TCGA stage-stratified pressure proxies."""
    expr = _read_log2_wide(PATHS.mirna_ev_gse241785_log2_norm)
    samples = pd.read_csv(PATHS.mirna_ev_gse241785_samples, sep="\t")
    samples = samples[samples["clinical_stage"].isin(STAGE_ORDER)].copy()

    tcga_stage = _load_tcga_stage_pressure()
    arms = [a for a in REPLICATION_FOCUS_ARMS if a in expr.index]

    ev_rows: List[dict] = []
    for arm in arms:
        for stage in STAGE_ORDER:
            ids = samples.loc[samples["clinical_stage"] == stage, "sample_id"].astype(str).tolist()
            ids = [i for i in ids if i in expr.columns]
            vals = pd.to_numeric(expr.loc[arm, ids], errors="coerce").dropna()
            ev_rows.append(
                {
                    "arm": arm,
                    "clinical_stage": stage,
                    "n_ev_samples": len(vals),
                    "ev_median_log2_norm": float(vals.median()) if len(vals) else np.nan,
                    "ev_mean_log2_norm": float(vals.mean()) if len(vals) else np.nan,
                }
            )
    ev_stage = pd.DataFrame(ev_rows)
    ev_stage.to_csv(out_dir / "gse241785_stage_ev_medians.tsv", sep="\t", index=False)

    tcga_sub = tcga_stage[tcga_stage["arm"].isin(arms)].copy()
    tcga_sub = tcga_sub.rename(
        columns={
            "stratum": "clinical_stage",
            "n_pressure_anticorr_genes_p_lt_0_05": "tcga_n_anticorr_stage",
            "strongest_partial_rho": "tcga_strongest_partial_rho",
            "n_samples": "tcga_n_tumor_samples",
        }
    )

    merged = ev_stage.merge(
        tcga_sub[
            [
                "arm",
                "clinical_stage",
                "tcga_n_anticorr_stage",
                "tcga_strongest_partial_rho",
                "tcga_n_tumor_samples",
                "strongest_gene",
            ]
        ],
        on=["arm", "clinical_stage"],
        how="left",
    )

    # Per-arm Spearman: EV median vs TCGA stage anticorr count (I–III where both defined)
    summary_rows: List[dict] = []
    for arm in arms:
        sub = merged[(merged["arm"] == arm) & merged["clinical_stage"].isin(["I", "II", "III"])].copy()
        sub = sub.dropna(subset=["ev_median_log2_norm", "tcga_n_anticorr_stage"])
        rho_n = p_n = rho_rho = np.nan
        if len(sub) >= 3:
            rho_n, p_n = spearmanr(sub["ev_median_log2_norm"], sub["tcga_n_anticorr_stage"])
            rho_rho, _ = spearmanr(sub["ev_median_log2_norm"], sub["tcga_strongest_partial_rho"].abs())
        row: dict = {"arm": arm}
        for stage in STAGE_ORDER:
            st = merged[(merged["arm"] == arm) & (merged["clinical_stage"] == stage)]
            if len(st):
                row[f"n_ev_stage_{stage}"] = int(st["n_ev_samples"].iloc[0])
                row[f"ev_median_stage_{stage}"] = float(st["ev_median_log2_norm"].iloc[0])
                row[f"tcga_n_anticorr_stage_{stage}"] = float(st["tcga_n_anticorr_stage"].iloc[0])
                if pd.notna(st["tcga_strongest_partial_rho"].iloc[0]):
                    row[f"tcga_strongest_rho_stage_{stage}"] = float(st["tcga_strongest_partial_rho"].iloc[0])
        row["spearman_ev_vs_tcga_n_anticorr_I_II_III"] = float(rho_n) if np.isfinite(rho_n) else np.nan
        row["p_ev_vs_tcga_n_anticorr"] = float(p_n) if np.isfinite(p_n) else np.nan
        row["spearman_ev_vs_tcga_abs_rho_I_II_III"] = float(rho_rho) if np.isfinite(rho_rho) else np.nan
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)

    anticorr = _anticorr_summary_for_arms(arms)
    if not anticorr.empty:
        summary = summary.merge(
            anticorr[["arm", "n_anticorr_gene_x_pam50", "n_ecm_anticorr_hits"]],
            on="arm",
            how="left",
        )

    summary.to_csv(out_dir / "gse241785_stage_vs_tcga_pressure.tsv", sep="\t", index=False)
    merged.to_csv(out_dir / "gse241785_stage_ev_tcga_long.tsv", sep="\t", index=False)

    print(
        f"[ev_followup] GSE241785 stage strata: "
        f"{samples['clinical_stage'].value_counts().sort_index().to_dict()}"
    )
    return merged, summary


def plot_cancer_benign_focus(focus: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    if focus.empty or "gse255660_cancer_benign_delta" not in focus.columns:
        return
    sub = focus.sort_values("gse255660_cancer_benign_delta", ascending=False)
    fig, ax = plt.subplots(figsize=(max(7, len(sub) * 0.55), 4.5))
    colors = ["#C44E52" if v > 0 else "#4C72B0" for v in sub["gse255660_cancer_benign_delta"]]
    ax.bar(range(len(sub)), sub["gse255660_cancer_benign_delta"], color=colors, alpha=0.85)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub["arm"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Δ median log2 (cancer − benign)")
    ax.set_title("GSE255660 serum miRNA — cancer vs benign (focus arms)")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_stage_ev_vs_tcga(summary: pd.DataFrame, merged: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    arms = [a for a in list(MIR29_ARMS) + list(REPLICATION_FOCUS_ARMS) if a in set(summary["arm"])]
    arms = list(dict.fromkeys(a for a in arms if a in set(summary["arm"])))[:8]
    if not arms:
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: EV median by stage (focus subset)
    mat = merged[merged["arm"].isin(arms)].pivot_table(
        index="arm", columns="clinical_stage", values="ev_median_log2_norm"
    )
    mat = mat.reindex(columns=[c for c in STAGE_ORDER if c in mat.columns])
    vmax = np.nanpercentile(np.abs(mat.values.astype(float)), 90) if mat.size else 1
    vmax = max(vmax, 0.5)
    im0 = axes[0].imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[0].set_xticks(range(mat.shape[1]))
    axes[0].set_xticklabels(mat.columns)
    axes[0].set_yticks(range(mat.shape[0]))
    axes[0].set_yticklabels(mat.index, fontsize=8)
    axes[0].set_title("GSE241785 vesicular EV (median log2 norm)")
    # Right: TCGA stage anticorr count by stage
    tmat = merged[merged["arm"].isin(arms)].pivot_table(
        index="arm", columns="clinical_stage", values="tcga_n_anticorr_stage"
    )
    tmat = tmat.reindex(columns=[c for c in ("I", "II", "III") if c in tmat.columns])
    im1 = axes[1].imshow(tmat.values.astype(float), aspect="auto", cmap="YlOrRd")
    axes[1].set_xticks(range(tmat.shape[1]))
    axes[1].set_xticklabels([f"TCGA stage {c}" for c in tmat.columns], fontsize=8)
    axes[1].set_yticks(range(tmat.shape[0]))
    axes[1].set_yticklabels(tmat.index, fontsize=8)
    axes[1].set_title("TCGA pressure proxy (n anticorr genes, p<0.05)")
    for i in range(tmat.shape[0]):
        for j in range(tmat.shape[1]):
            val = tmat.values[i, j]
            if np.isfinite(val):
                axes[1].text(j, i, f"{int(val)}", ha="center", va="center", fontsize=7)
    fig.colorbar(im0, ax=axes[0], shrink=0.8, label="EV log2 norm")
    fig.colorbar(im1, ax=axes[1], shrink=0.8, label="n genes")
    fig.suptitle("GSE241785 stage stratification vs TCGA stage pressure", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(*, out_dir: Path = C.PLASMA_EV_DIR) -> Dict[str, object]:
    multi_dir = out_dir / "multi_cohort"
    multi_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    _, focus_cb = run_gse255660_cancer_vs_benign(out_dir=multi_dir)
    plot_cancer_benign_focus(focus_cb, fig_dir / "gse255660_cancer_vs_benign_focus.png")

    merged_stage, summary_stage = run_gse241785_stage_pressure(out_dir=multi_dir)
    plot_stage_ev_vs_tcga(summary_stage, merged_stage, fig_dir / "gse241785_stage_ev_vs_tcga_pressure.png")

    # Extend multi-cohort wide table if present
    wide_path = multi_dir / "ev_multi_cohort_pressure_summary.tsv"
    if wide_path.is_file() and not focus_cb.empty:
        wide = pd.read_csv(wide_path, sep="\t")
        add = focus_cb[
            ["arm", "gse255660_cancer_benign_delta", "gse255660_cancer_benign_q", "gse255660_cancer_benign_auc"]
        ]
        wide = wide.merge(add, on="arm", how="left")
        wide.to_csv(multi_dir / "ev_multi_cohort_pressure_summary_extended.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.analyses.dcis_ev.ev_mirna_followup",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analyses": [
            "GSE255660_cancer_vs_benign",
            "GSE241785_clinical_stage_vs_tcga_stage_pressure",
        ],
        "tcga_stage_bridge": str(TCGA_STAGE_BRIDGE),
    }
    (multi_dir / "followup_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\n[ev_followup] GSE255660 cancer vs benign (focus arms):\n")
    show_cb = [
        c
        for c in [
            "arm",
            "gse255660_cancer_benign_delta",
            "gse255660_cancer_benign_q",
            "gse255660_cancer_benign_auc",
            "gse270497_screening_delta",
            "same_direction_as_gse270497_screening",
            "n_anticorr_gene_x_pam50",
        ]
        if c in focus_cb.columns
    ]
    print(focus_cb[show_cb].to_string(index=False))

    print("\n[ev_followup] GSE241785 stage × TCGA pressure (focus arms):\n")
    show_st = [
        c
        for c in summary_stage.columns
        if c
        in (
            "arm",
            "ev_median_stage_I",
            "ev_median_stage_II",
            "ev_median_stage_III",
            "tcga_n_anticorr_stage_I",
            "tcga_n_anticorr_stage_II",
            "tcga_n_anticorr_stage_III",
            "spearman_ev_vs_tcga_n_anticorr_I_II_III",
            "n_ecm_anticorr_hits",
        )
    ]
    print(summary_stage[show_st].head(12).to_string(index=False))

    return {"cancer_benign_focus": focus_cb, "stage_summary": summary_stage}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.PLASMA_EV_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
