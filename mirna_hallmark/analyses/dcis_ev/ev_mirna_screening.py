"""Plasma exosome miRNA: GSE270497 screening+ vs screening- contrasts.

Translational layer orthogonal to TCGA tumor tissue miRNA / CN dosage:
compares circulating EV miRNA abundance between breast cancer screening-positive
and screening-negative women (GEO GSE270497, n=180).

Inputs (after ingest):
  ``data/miRNA/ev/GSE270497/mirna_log2_tpm_plus1.tsv.gz``
  ``data/miRNA/ev/GSE270497/samples.tsv``

Outputs under ``mirna_hallmark/output/plasma_ev/``:
  ``gse270497_screening_arm_contrasts.tsv`` — all arms, Mann-Whitney + BH-FDR
  ``gse270497_focus_panel_contrasts.tsv`` — ``PANEL_MIRNA_TIER_ARM_IDS`` subset
  ``gse270497_tcga_tissue_rank_compare.tsv`` — cohort-level EV delta vs TCGA median
  ``figures/gse270497_focus_panel_log2tpm.png`` — focus-arm boxplot
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
from mirna_hallmark import stats as S
from pipeline.config import PANEL_MIRNA_TIER_ARM_IDS, PATHS
from pipeline.miRNA_exp.loaders import map_panel_arms_to_mimat_ids

FEATURE_COL = PATHS.mirna_ev_gene_col

# miR-29 family arms for cross-layer EV / tissue summaries.
MIR29_ARMS: Tuple[str, ...] = (
    "hsa-miR-29a-3p",
    "hsa-miR-29b-3p",
    "hsa-miR-29c-3p",
)


def load_ev_log2_matrix(path: Optional[Path] = None) -> pd.DataFrame:
    """Wide log2 abundance matrix: index = miRNA_arm, columns = sample_id."""
    p = Path(path or C.EV_MIRNA_GSE270497_LOG2_TPM)
    if not p.is_file():
        raise FileNotFoundError(
            f"Missing {p}. Run:\n"
            "  bash scripts/mirna/download_geo_ev_gse270497.sh\n"
            "  .venv/bin/python3 scripts/mirna/ingest_geo_ev_gse270497.py"
        )
    df = pd.read_csv(p, sep="\t")
    feat = df.pop(FEATURE_COL)
    df.index = feat.astype(str).str.strip()
    return df


def load_ev_samples(path: Optional[Path] = None) -> pd.DataFrame:
    p = Path(path or C.EV_MIRNA_GSE270497_SAMPLES)
    if not p.is_file():
        raise FileNotFoundError(f"Missing {p}; run ingest_geo_ev_gse270497.py first")
    return pd.read_csv(p, sep="\t")


def _mann_whitney_delta(
    values: pd.Series,
    labels: pd.Series,
) -> Tuple[float, float, float, float, int, int]:
    """Return (delta_median_pos_minus_neg, u_stat, p, rank_biserial, n_pos, n_neg)."""
    from scipy.stats import mannwhitneyu

    v = pd.to_numeric(values, errors="coerce")
    g = labels.astype(bool)
    pos = v.loc[g].dropna()
    neg = v.loc[~g].dropna()
    n_pos, n_neg = len(pos), len(neg)
    if n_pos < 3 or n_neg < 3:
        return np.nan, np.nan, np.nan, np.nan, n_pos, n_neg
    delta = float(pos.median() - neg.median())
    try:
        res = mannwhitneyu(pos, neg, alternative="two-sided")
        u = float(res.statistic)
        p = float(res.pvalue)
    except ValueError:
        return delta, np.nan, np.nan, np.nan, n_pos, n_neg
    # rank-biserial: 1 - 2U/(n1*n2)
    rbc = 1.0 - (2.0 * u / (n_pos * n_neg))
    return delta, u, p, float(rbc), n_pos, n_neg


def case_control_arm_contrasts(
    expr: pd.DataFrame,
    samples: pd.DataFrame,
    *,
    label_col: str,
    case_label: str = "case",
    control_label: str = "control",
    min_total_samples: int = 10,
) -> pd.DataFrame:
    """Per-arm Mann-Whitney on log2 abundance; ``label_col`` must be boolean-like."""
    sample_ids = [c for c in expr.columns if c in set(samples["sample_id"])]
    if len(sample_ids) < min_total_samples:
        raise ValueError(
            f"Too few overlapping samples between expression matrix and samples.tsv "
            f"({len(sample_ids)} < {min_total_samples})"
        )

    samp = samples.set_index("sample_id").loc[sample_ids]
    labels = samp[label_col].astype(bool)

    rows: List[Dict[str, object]] = []
    for arm, row in expr[sample_ids].iterrows():
        delta, u, p, rbc, n_case, n_ctrl = _mann_whitney_delta(row, labels)
        rows.append(
            {
                "arm": arm,
                "delta_median_log2_case_minus_control": delta,
                "mannwhitney_u": u,
                "mannwhitney_p": p,
                "rank_biserial": rbc,
                f"n_{case_label}": n_case,
                f"n_{control_label}": n_ctrl,
                f"median_log2_{case_label}": float(pd.to_numeric(row.loc[labels], errors="coerce").median()),
                f"median_log2_{control_label}": float(pd.to_numeric(row.loc[~labels], errors="coerce").median()),
            }
        )

    out = pd.DataFrame(rows)
    m = out["mannwhitney_p"].notna()
    out.loc[m, "mannwhitney_q"] = S.bh_fdr(out.loc[m, "mannwhitney_p"].values)
    return out.sort_values("mannwhitney_p", na_position="last")


def screening_arm_contrasts(
    expr: pd.DataFrame,
    samples: pd.DataFrame,
) -> pd.DataFrame:
    """Per-arm screening+ vs screening- on log2(TPM+1)."""
    out = case_control_arm_contrasts(
        expr,
        samples,
        label_col="screening_positive",
        case_label="screening_positive",
        control_label="screening_negative",
    )
    return out.rename(
        columns={
            "delta_median_log2_case_minus_control": "delta_median_log2tpm_pos_minus_neg",
            "median_log2_screening_positive": "median_log2tpm_positive",
            "median_log2_screening_negative": "median_log2tpm_negative",
        }
    )


def tcga_tissue_median_by_arm() -> pd.DataFrame:
    """Cohort median log2(RPM+1) per panel arm in TCGA primary tumors."""
    from pipeline.sample_ids import tcga_sample_type_two_digit

    tcga_path = PATHS.mirna_expression
    if not tcga_path.is_file():
        return pd.DataFrame()

    header = pd.read_csv(tcga_path, sep="\t", nrows=0)
    cols = [c for c in header.columns if c != PATHS.mirna_gene_col]
    primary_cols = [c for c in cols if tcga_sample_type_two_digit(c) == "01"]
    if not primary_cols:
        return pd.DataFrame()

    usecols = [PATHS.mirna_gene_col] + primary_cols
    tcga = pd.read_csv(tcga_path, sep="\t", usecols=usecols)
    mimat_ids = tcga[PATHS.mirna_gene_col].astype(str).tolist()
    arm_map = map_panel_arms_to_mimat_ids(PANEL_MIRNA_TIER_ARM_IDS, mimat_ids)

    rows = []
    tcga_idx = tcga.set_index(PATHS.mirna_gene_col)
    for arm, mimats in arm_map.items():
        if not mimats:
            continue
        sub = tcga_idx.loc[[m for m in mimats if m in tcga_idx.index]]
        if sub.empty:
            continue
        med = sub[primary_cols].astype(float).median(axis=1).median()
        rows.append({"arm": arm, "tcga_primary_median_log2rpm": float(med), "mimat": mimats[0]})
    return pd.DataFrame(rows)


def compare_ev_delta_to_tcga_tissue(
    contrasts: pd.DataFrame,
    tcga_medians: pd.DataFrame,
    *,
    delta_col: str = "delta_median_log2tpm_pos_minus_neg",
) -> pd.DataFrame:
    """Rank correlation: EV case-control delta vs TCGA tissue median (panel arms)."""
    if tcga_medians.empty or delta_col not in contrasts.columns:
        return pd.DataFrame()

    merged = contrasts.merge(tcga_medians, on="arm", how="inner")
    if len(merged) < 3:
        return merged

    from scipy.stats import spearmanr

    x = merged[delta_col].astype(float)
    y = merged["tcga_primary_median_log2rpm"].astype(float)
    rho, p = spearmanr(x, y)
    merged.attrs["spearman_rho"] = float(rho)
    merged.attrs["spearman_p"] = float(p)
    return merged.sort_values(delta_col, ascending=False)


def plot_focus_panel_boxplot(
    expr: pd.DataFrame,
    samples: pd.DataFrame,
    focus_arms: Sequence[str],
    out_path: Path,
) -> None:
    import matplotlib.pyplot as plt

    sample_ids = [c for c in expr.columns if c in set(samples["sample_id"])]
    samp = samples.set_index("sample_id").loc[sample_ids]
    labels = samp["screening_positive"].map({True: "screening+", False: "screening-"})

    arms = [a for a in focus_arms if a in expr.index]
    if not arms:
        return

    n = len(arms)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.0 * nrow), squeeze=False)
    status_order = ["screening-", "screening+"]
    palette = {"screening-": "#4C72B0", "screening+": "#C44E52"}

    for i, arm in enumerate(arms):
        ax = axes[i // ncol][i % ncol]
        vals = pd.to_numeric(expr.loc[arm, sample_ids], errors="coerce")
        df = pd.DataFrame({"log2(TPM+1)": vals.values, "group": labels.values})
        groups = [df.loc[df["group"] == g, "log2(TPM+1)"].dropna().values for g in status_order]
        bp = ax.boxplot(groups, tick_labels=status_order, patch_artist=True)
        for patch, g in zip(bp["boxes"], status_order):
            patch.set_facecolor(palette[g])
            patch.set_alpha(0.75)
        ax.set_title(arm, fontsize=9)
        ax.set_ylabel("log2(TPM+1)" if i % ncol == 0 else "")

    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    fig.suptitle("GSE270497 plasma exosome miRNA — focus panel", fontsize=11, y=1.02)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(
    *,
    out_dir: Path = C.PLASMA_EV_DIR,
    expr_path: Optional[Path] = None,
    samples_path: Optional[Path] = None,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    expr = load_ev_log2_matrix(expr_path)
    samples = load_ev_samples(samples_path)

    contrasts = screening_arm_contrasts(expr, samples)
    contrasts.to_csv(out_dir / "gse270497_screening_arm_contrasts.tsv", sep="\t", index=False)

    focus = contrasts.loc[contrasts["arm"].isin(PANEL_MIRNA_TIER_ARM_IDS)].copy()
    focus.to_csv(out_dir / "gse270497_focus_panel_contrasts.tsv", sep="\t", index=False)

    tcga_med = tcga_tissue_median_by_arm()
    compare = compare_ev_delta_to_tcga_tissue(focus, tcga_med)
    compare.to_csv(out_dir / "gse270497_tcga_tissue_rank_compare.tsv", sep="\t", index=False)

    plot_focus_panel_boxplot(
        expr,
        samples,
        PANEL_MIRNA_TIER_ARM_IDS,
        out_dir / "figures" / "gse270497_focus_panel_log2tpm.png",
    )

    n_sig = int((contrasts["mannwhitney_q"] < C.FDR_ALPHA).sum())
    print(
        f"[ev_mirna] GSE270497: {len(contrasts)} arms tested, "
        f"{n_sig} FDR<{C.FDR_ALPHA}; focus panel n={len(focus)}"
    )
    if not focus.empty:
        top = focus.nsmallest(3, "mannwhitney_p")
        tops = ", ".join(
            f"{r.arm} Δ={r.delta_median_log2tpm_pos_minus_neg:.2f}"
            for r in top.itertuples()
        )
        print(f"[ev_mirna] focus panel smallest-p: {tops}")

    manifest = {
        "module": "mirna_hallmark.analyses.dcis_ev.ev_mirna_screening",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort": "GSE270497",
        "compartment": "plasma_exosome",
        "n_samples": int(len(samples)),
        "n_arms_tested": int(len(contrasts)),
        "n_fdr_significant": n_sig,
        "fdr_alpha": C.FDR_ALPHA,
        "focus_panel": list(PANEL_MIRNA_TIER_ARM_IDS),
        "tcga_compare_spearman_rho": compare.attrs.get("spearman_rho"),
        "tcga_compare_spearman_p": compare.attrs.get("spearman_p"),
        "note": "Unpaired screening cohort; not joinable to TCGA participants.",
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[ev_mirna] wrote outputs under {out_dir}")
    return contrasts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.PLASMA_EV_DIR)
    ap.add_argument("--expr", type=Path, default=None, help="Override log2 TPM matrix")
    ap.add_argument("--samples", type=Path, default=None, help="Override samples.tsv")
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, expr_path=args.expr, samples_path=args.samples)


if __name__ == "__main__":
    main()
