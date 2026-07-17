"""Six-cohort EV/plasma miRNA comparison with TCGA pressure spine.

Contrasts (primary binary design per cohort):
  GSE270497  screening+ vs screening− (plasma EV-seq)
  GSE197020  cancer vs normal (plasma exosome-seq)
  GSE255660  cancer vs healthy serum (excl. dilution replicates)
  GSE301416  BRCA1/2 carrier vs wild-type serum qPCR
  GSE241784  Basal vs LumA vesicular NanoString (small n)
  GSE241785  Basal vs LumA vesicular NanoString (young; stage metadata)

Outputs under ``mirna_hallmark/output/plasma_ev/``:
  ``multi_cohort/ev_multi_cohort_arm_long.tsv``
  ``multi_cohort/ev_multi_cohort_focus_wide.tsv``
  ``multi_cohort/ev_multi_cohort_pressure_summary.tsv``
  ``figures/ev_multi_cohort_focus_heatmap.png``

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.ev_mirna_multi_cohort
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.ev_mirna_replication import (
    ECM_ANTICORR_GENES,
    REPLICATION_FOCUS_ARMS,
    _anticorr_summary_for_arms,
    _arm_auc,
)
from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import (
    FEATURE_COL,
    MIR29_ARMS,
    case_control_arm_contrasts,
    load_ev_log2_matrix,
    load_ev_samples,
)

from mirna_hallmark.analyses.dcis_ev.ev_mirna_replication import load_gse197020_matrix, load_gse197020_samples


@dataclass(frozen=True)
class CohortSpec:
    geo: str
    design: str
    matrix_path: Path
    samples_path: Path
    prepare: Callable[[pd.DataFrame], pd.DataFrame]
    label_col: str
    case_label: str
    control_label: str
    compartment: str


def _read_log2_wide(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path, sep="\t")
    feat = df.pop(FEATURE_COL)
    df.index = feat.astype(str).str.strip()
    return df


def _prep_270497(samples: pd.DataFrame) -> pd.DataFrame:
    return samples.copy()


def _prep_197020(samples: pd.DataFrame) -> pd.DataFrame:
    return samples.copy()


def _prep_255660(samples: pd.DataFrame) -> pd.DataFrame:
    s = samples.loc[~samples["is_dilution_replicate"].astype(bool)].copy()
    s = s[s["bc_case"].astype(bool) | s["healthy_control"].astype(bool)]
    s["cancer_vs_healthy"] = s["bc_case"].astype(bool)
    return s


def _prep_301416(samples: pd.DataFrame) -> pd.DataFrame:
    s = samples[
        samples["brca_status"].isin(["BRCA1 mutation", "BRCA2 mutation", "BRCA wild-type"])
    ].copy()
    s["brca_carrier"] = s["brca_status"].astype(str).str.contains("mutation", case=False)
    return s


def _prep_basal_vs_nonbasal(samples: pd.DataFrame) -> pd.DataFrame:
    s = samples.copy()
    s["basal_vs_nonbasal"] = s["pam50_proxy"].astype(str) == "Basal"
    return s


def _prep_basal_vs_luma(samples: pd.DataFrame) -> pd.DataFrame:
    s = samples[samples["pam50_proxy"].isin(["Basal", "LumA"])].copy()
    s["basal_vs_luma"] = s["pam50_proxy"].astype(str) == "Basal"
    return s


def cohort_registry() -> List[CohortSpec]:
    from pipeline.config import PATHS

    return [
        CohortSpec(
            "GSE270497",
            "screening_positive vs screening_negative",
            PATHS.mirna_ev_gse270497_log2_tpm,
            PATHS.mirna_ev_gse270497_samples,
            _prep_270497,
            "screening_positive",
            "screening_positive",
            "screening_negative",
            "plasma_exosome",
        ),
        CohortSpec(
            "GSE197020",
            "breast_cancer vs normal",
            PATHS.mirna_ev_gse197020_log2_cpm,
            PATHS.mirna_ev_gse197020_samples,
            _prep_197020,
            "bc_case",
            "cancer",
            "normal",
            "plasma_exosome",
        ),
        CohortSpec(
            "GSE255660",
            "breast_cancer vs healthy (serum; excl. dilution reps)",
            PATHS.mirna_ev_gse255660_log2_expr,
            PATHS.mirna_ev_gse255660_samples,
            _prep_255660,
            "cancer_vs_healthy",
            "cancer",
            "healthy",
            "serum",
        ),
        CohortSpec(
            "GSE301416",
            "BRCA1/2 carrier vs wild-type (serum qPCR)",
            PATHS.mirna_ev_gse301416_log2_dcq,
            PATHS.mirna_ev_gse301416_samples,
            _prep_301416,
            "brca_carrier",
            "carrier",
            "wild_type",
            "serum",
        ),
        CohortSpec(
            "GSE241784",
            "Basal vs non-Basal (vesicular NanoString; n=12)",
            PATHS.mirna_ev_gse241784_log2_norm,
            PATHS.mirna_ev_gse241784_samples,
            _prep_basal_vs_nonbasal,
            "basal_vs_nonbasal",
            "basal",
            "non_basal",
            "plasma_vesicular",
        ),
        CohortSpec(
            "GSE241785",
            "Basal vs LumA young BC (vesicular; stage in metadata)",
            PATHS.mirna_ev_gse241785_log2_norm,
            PATHS.mirna_ev_gse241785_samples,
            _prep_basal_vs_luma,
            "basal_vs_luma",
            "basal",
            "lumA",
            "plasma_vesicular",
        ),
    ]


def _harmonize_contrasts(contrasts: pd.DataFrame, spec: CohortSpec) -> pd.DataFrame:
    out = contrasts.copy()
    delta_col = "delta_median_log2_case_minus_control"
    if delta_col not in out.columns:
        for c in out.columns:
            if c.startswith("delta_median_log2"):
                out = out.rename(columns={c: "delta"})
                break
    else:
        out = out.rename(columns={delta_col: "delta"})
    out["cohort"] = spec.geo
    out["design"] = spec.design
    out["compartment"] = spec.compartment
    qcol = "mannwhitney_q" if "mannwhitney_q" in out.columns else "mannwhitney_p"
    out["q"] = out[qcol]
    out["sig_fdr05"] = out["q"] < C.FDR_ALPHA
    # signed effect: positive = higher in case group
    out["direction"] = np.where(out["delta"] > 0, "case_higher", np.where(out["delta"] < 0, "control_higher", "flat"))
    return out


def run_cohort_contrast(spec: CohortSpec) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    samples = pd.read_csv(spec.samples_path, sep="\t")
    samples = spec.prepare(samples)
    expr = _read_log2_wide(spec.matrix_path)
    contrasts = case_control_arm_contrasts(
        expr,
        samples,
        label_col=spec.label_col,
        case_label=spec.case_label,
        control_label=spec.control_label,
        min_total_samples=6 if spec.geo.startswith("GSE241") else 10,
    )
    auc = _arm_auc(expr, samples, spec.label_col)
    harmonized = _harmonize_contrasts(contrasts, spec)
    harmonized["n_case"] = harmonized.get(f"n_{spec.case_label}", np.nan)
    harmonized["n_control"] = harmonized.get(f"n_{spec.control_label}", np.nan)
    harmonized["auc"] = harmonized["arm"].map(auc)
    return harmonized, samples, auc


def build_focus_wide(long_df: pd.DataFrame, arms: Sequence[str]) -> pd.DataFrame:
    ref = "GSE270497"
    rows: List[dict] = []
    for arm in arms:
        r: dict = {"arm": arm}
        sub = long_df[long_df["arm"] == arm]
        for _, row in sub.iterrows():
            c = row["cohort"]
            r[f"{c}_delta"] = row["delta"]
            r[f"{c}_q"] = row["q"]
            r[f"{c}_auc"] = row["auc"]
            r[f"{c}_sig"] = bool(row["sig_fdr05"]) if pd.notna(row["sig_fdr05"]) else False
        ref_row = sub[sub["cohort"] == ref]
        ref_sign = np.sign(ref_row["delta"].iloc[0]) if len(ref_row) else np.nan
        same_dir = 0
        comp = 0
        for _, row in sub.iterrows():
            if row["cohort"] == ref or pd.isna(row["delta"]):
                continue
            if pd.isna(ref_sign) or ref_sign == 0 or pd.isna(row["delta"]) or row["delta"] == 0:
                continue
            comp += 1
            if np.sign(row["delta"]) == ref_sign:
                same_dir += 1
        r["n_cohorts_vs_ref"] = comp
        r["n_same_direction_as_gse270497"] = same_dir
        r["frac_same_direction_as_gse270497"] = same_dir / comp if comp else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def _pressure_tier(row: pd.Series) -> str:
    n = row.get("n_anticorr_gene_x_pam50", np.nan)
    if pd.isna(n) or n == 0:
        return "none"
    if n >= 20:
        return "strong"
    if n >= 5:
        return "moderate"
    return "weak"


def build_pressure_summary(wide: pd.DataFrame, anticorr: pd.DataFrame) -> pd.DataFrame:
    if anticorr.empty:
        return wide.copy()
    out = wide.merge(anticorr, left_on="arm", right_on="arm", how="left")
    out["tcga_pressure_tier"] = out.apply(_pressure_tier, axis=1)
    return out


def plot_focus_heatmap(wide: pd.DataFrame, cohorts: Sequence[str], out_path: Path) -> None:
    import matplotlib.pyplot as plt

    delta_cols = [f"{c}_delta" for c in cohorts if f"{c}_delta" in wide.columns]
    if not delta_cols:
        return
    mat = wide.set_index("arm")[delta_cols].astype(float)
    mat.columns = [c.replace("_delta", "") for c in mat.columns]
    mat = mat.loc[[a for a in REPLICATION_FOCUS_ARMS if a in mat.index]]
    if mat.empty:
        return

    vmax = np.nanpercentile(np.abs(mat.values), 95)
    vmax = max(vmax, 0.25)
    fig, ax = plt.subplots(figsize=(max(7, len(mat.columns) * 1.1), max(4, len(mat) * 0.35)))
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_title("EV/plasma miRNA Δ (case − control) across cohorts — focus arms")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color="0.2")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Δ median log2 abundance")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(*, out_dir: Path = C.PLASMA_EV_DIR) -> Tuple[pd.DataFrame, pd.DataFrame]:
    multi_dir = out_dir / "multi_cohort"
    multi_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    specs = cohort_registry()
    long_parts: List[pd.DataFrame] = []
    cohort_meta: Dict[str, dict] = {}

    for spec in specs:
        try:
            harmonized, samples, _auc = run_cohort_contrast(spec)
        except FileNotFoundError as exc:
            print(f"[ev_multi_cohort] SKIP {spec.geo}: {exc}")
            continue
        harmonized.to_csv(multi_dir / f"{spec.geo.lower()}_arm_contrasts.tsv", sep="\t", index=False)
        long_parts.append(
            harmonized[
                [
                    "cohort",
                    "design",
                    "compartment",
                    "arm",
                    "delta",
                    "q",
                    "sig_fdr05",
                    "auc",
                    "n_case",
                    "n_control",
                    "direction",
                ]
            ]
        )
        n_case = int(samples[spec.label_col].astype(bool).sum())
        n_ctrl = int((~samples[spec.label_col].astype(bool)).sum())
        cohort_meta[spec.geo] = {
            "design": spec.design,
            "compartment": spec.compartment,
            "n_case": n_case,
            "n_control": n_ctrl,
            "n_arms_tested": int(len(harmonized)),
            "n_sig_fdr05": int(harmonized["sig_fdr05"].fillna(False).sum()),
        }
        print(
            f"[ev_multi_cohort] {spec.geo}: {n_case} case / {n_ctrl} control; "
            f"{cohort_meta[spec.geo]['n_sig_fdr05']} arms FDR<0.05"
        )

    if not long_parts:
        raise RuntimeError("No cohorts ingested; run download/ingest scripts first")

    long_df = pd.concat(long_parts, ignore_index=True)
    long_df.to_csv(multi_dir / "ev_multi_cohort_arm_long.tsv", sep="\t", index=False)

    arms = list(
        dict.fromkeys(
            list(REPLICATION_FOCUS_ARMS)
            + [a for a in long_df["arm"] if a in set(REPLICATION_FOCUS_ARMS)]
        )
    )
    arms = [a for a in arms if a in set(long_df["arm"])]

    wide = build_focus_wide(long_df, arms)
    anticorr = _anticorr_summary_for_arms(arms)
    summary = build_pressure_summary(wide, anticorr)
    summary.to_csv(multi_dir / "ev_multi_cohort_focus_wide.tsv", sep="\t", index=False)
    summary.to_csv(multi_dir / "ev_multi_cohort_pressure_summary.tsv", sep="\t", index=False)

    cohort_order = [s.geo for s in specs if s.geo in cohort_meta]
    plot_focus_heatmap(summary, cohort_order, out_dir / "figures" / "ev_multi_cohort_focus_heatmap.png")

    # miR-29 family slice
    m29 = summary[summary["arm"].isin(MIR29_ARMS)].copy()
    m29.to_csv(multi_dir / "ev_multi_cohort_mir29_wide.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.analyses.dcis_ev.ev_mirna_multi_cohort",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohorts": cohort_meta,
        "reference_cohort": "GSE270497",
        "focus_arms": list(REPLICATION_FOCUS_ARMS),
        "note": (
            "Deltas are per-cohort log2 scale (TPM/CPM/expr/norm/dCq); compare direction and "
            "AUC across cohorts, not absolute magnitude. TCGA pressure from target_combined_anticorr."
        ),
    }
    (multi_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # Console summary for focus arms
    show_cols = ["arm", "tcga_pressure_tier", "n_anticorr_gene_x_pam50", "frac_same_direction_as_gse270497"]
    show_cols += [f"{g}_delta" for g in cohort_order if f"{g}_delta" in summary.columns]
    print("\n[ev_multi_cohort] Focus arms × cohorts (Δ case−control):\n")
    print(summary[show_cols].to_string(index=False))

    from mirna_hallmark.analyses.dcis_ev.ev_mirna_followup import run as run_followup

    run_followup(out_dir=out_dir)

    return long_df, summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.PLASMA_EV_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
