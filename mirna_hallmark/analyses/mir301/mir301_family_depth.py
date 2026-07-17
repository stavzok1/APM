"""miR-301/130/454 family depth around published and TS-orphan targets.

This module fills two gaps left by the single-arm miR-301a screen:

1. per-(gene, miRNA) realized pressure shares under the default and zless
   pressure variants; and
2. single-arm plus family-aggregate partial correlations for the shared
   TargetScan seed family (301a/301b/130a/130b/454 3p arms).

Outputs under ``mirna_hallmark/output/tissue_reference/mir301_family_depth/``.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.mir301.mir301_family_depth
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.mir301.mir301_focus_genes import PRIOR_ORPHAN_HITS, PUBLISHED_BRCA
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import compute_gene_pressure_contributions
from mirna_hallmark.robustness_checks import _pam50_scope_iter, _partial_ladder, _proliferation_proxies
from mirna_hallmark.robustness_checks import _load_targetscan_weights
from mirna_hallmark.stats import bh_fdr, zscore_rows


OUT_DIR = C.TISSUE_REFERENCE_DIR / "mir301_family_depth"
MIR301A_OUT_DIR = C.TISSUE_REFERENCE_DIR / "mir301a_depth"
MIN_SCOPE_N = 20
FDR_ALPHA = 0.05

# Variability decision rule (log2(RPM+1) spread across the cohort). An arm below
# BOTH thresholds is flagged low-variability: its single-arm anti-correlation is
# under-powered (little dynamic range to track target expression) and the z-term
# in pressure near-zeroes it. Thresholds are deliberately conservative and are
# reported alongside data-driven tertiles so the cut is transparent, not magic.
LOW_VAR_SD = 0.5
LOW_VAR_IQR = 0.75

# Same-seed 3p family used for target attribution. 5p arms are reported from EV
# contrasts separately when present, but not treated as the same TargetScan family.
FAMILY_ARMS: tuple[str, ...] = (
    "hsa-miR-301a-3p",
    "hsa-miR-301b-3p",
    "hsa-miR-130a-3p",
    "hsa-miR-130b-3p",
    "hsa-miR-454-3p",
)

EV_RELATED_ARMS: tuple[str, ...] = (
    *FAMILY_ARMS,
    "hsa-miR-301a-5p",
    "hsa-miR-301b-5p",
    "hsa-miR-130a-5p",
    "hsa-miR-130b-5p",
    "hsa-miR-454-5p",
)

EXTRA_FOCUS_GENES: tuple[str, ...] = (
    "FRZB",
    "IGF1",
    "HOXA3",
    "NR3C2",
    "PPARG",
)

PRESSURE_VARIANTS: tuple[str, ...] = (
    "softmax_z_logrpm",
    "softmax_logrpm",
    "softmax_absratio",
    "softmax_z_blend",
    "softmax",
)


def _arm_variability(mirna: pd.DataFrame) -> pd.DataFrame:
    """Cohort-wide per-arm expression-variability table + low-variability flag.

    Variability is judged on log2(RPM+1) across all cohort samples (the same scale
    the pressure z-term and single-arm partials see). We report sd / IQR, their
    cohort tertiles, and a boolean ``low_variability`` (sd < LOW_VAR_SD AND
    IQR < LOW_VAR_IQR).
    """
    m = mirna.astype(float)
    sd = m.std(axis=1)
    iqr = m.quantile(0.75, axis=1) - m.quantile(0.25, axis=1)
    med = m.median(axis=1)
    out = pd.DataFrame(
        {
            "arm_median_log2rpm": med,
            "arm_sd_log2rpm": sd,
            "arm_iqr_log2rpm": iqr,
        }
    )
    # Data-driven tertiles among expressed arms (median log2rpm > 0) for context.
    expressed = out.loc[out["arm_median_log2rpm"] > 0, "arm_sd_log2rpm"]
    if len(expressed) >= 3:
        q1, q2 = expressed.quantile([1 / 3, 2 / 3]).tolist()
        out["arm_sd_tertile"] = pd.cut(
            out["arm_sd_log2rpm"], bins=[-np.inf, q1, q2, np.inf],
            labels=["low", "mid", "high"],
        ).astype(str)
    else:
        out["arm_sd_tertile"] = "na"
    out["low_variability"] = (out["arm_sd_log2rpm"] < LOW_VAR_SD) & (out["arm_iqr_log2rpm"] < LOW_VAR_IQR)
    return out


def _add_fdr_and_flags(partials: pd.DataFrame) -> pd.DataFrame:
    """BH-FDR within each (predictor_type, predictor, scope) test family + flags.

    Adds CN-adjusted and proliferation+CN-adjusted q-values and significance
    flags. ``neg_sig`` (raw p<0.05) is preserved for backward compatibility but
    should not be used for headline counts.
    """
    if partials.empty:
        return partials
    df = partials.copy()
    rho_cn = "rho_CPE_HRD_CN" if "rho_CPE_HRD_CN" in df.columns else "rho_CPE_HRD"
    p_cn = "p_CPE_HRD_CN" if "p_CPE_HRD_CN" in df.columns else "p_CPE_HRD"
    has_prolif = "rho_e2f_g2m_CN" in df.columns and "p_e2f_g2m_CN" in df.columns

    df["q_CPE_HRD_CN"] = np.nan
    if has_prolif:
        df["q_e2f_g2m_CN"] = np.nan

    grp_cols = ["predictor_type", "predictor", "scope"]
    for _, idx in df.groupby(grp_cols).groups.items():
        sub = df.loc[idx]
        df.loc[idx, "q_CPE_HRD_CN"] = bh_fdr(sub[p_cn].fillna(1.0).values)
        if has_prolif:
            df.loc[idx, "q_e2f_g2m_CN"] = bh_fdr(sub["p_e2f_g2m_CN"].fillna(1.0).values)

    df["neg_sig_cn_fdr"] = (df[rho_cn] < 0) & (df["q_CPE_HRD_CN"] < FDR_ALPHA)
    if has_prolif:
        df["neg_sig_prolif_cn_fdr"] = (df["rho_e2f_g2m_CN"] < 0) & (df["q_e2f_g2m_CN"] < FDR_ALPHA)
    else:
        df["neg_sig_prolif_cn_fdr"] = False
    return df


def _focus_genes(top_orphans: int = 25) -> list[str]:
    genes = set(PUBLISHED_BRCA) | set(PRIOR_ORPHAN_HITS) | set(EXTRA_FOCUS_GENES)
    cohort_path = MIR301A_OUT_DIR / "mir301a_full_ts_coupling_cohort.tsv"
    if cohort_path.is_file():
        c = pd.read_csv(cohort_path, sep="\t")
        rho_col = "rho_CPE_HRD_CN" if "rho_CPE_HRD_CN" in c.columns else "rho_CPE_HRD"
        orphan = c.loc[c.get("evidence_class", "").eq("ts_orphan")].copy()
        if not orphan.empty and rho_col in orphan.columns:
            top = orphan.sort_values(rho_col).head(top_orphans)["gene"].astype(str)
            genes.update(top.tolist())
    return sorted(genes)


def _targetscan_family_scores(genes: Sequence[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    ts = _load_targetscan_weights(list(genes))
    if ts.empty:
        return pd.DataFrame(), pd.DataFrame()
    long = ts.loc[ts["miRNA"].isin(FAMILY_ARMS) & ts["gene"].isin(set(genes))].copy()
    long = long.sort_values(["gene", "ts_weight"], ascending=[True, False])
    if long.empty:
        return long, pd.DataFrame()
    wide = long.pivot_table(index="gene", columns="miRNA", values="ts_weight", aggfunc="max")
    wide = wide.reindex(columns=list(FAMILY_ARMS))
    wide["family_ts_max"] = wide.max(axis=1)
    wide["family_ts_sum"] = wide.sum(axis=1)
    wide["n_family_arms_predicted"] = wide[list(FAMILY_ARMS)].notna().sum(axis=1)
    wide["best_family_arm"] = wide[list(FAMILY_ARMS)].idxmax(axis=1)
    return long, wide.reset_index().sort_values("family_ts_max", ascending=False)


def _single_arm_partials(genes: Sequence[str], ts_long: pd.DataFrame) -> pd.DataFrame:
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    clinical = D.load_clinical_strata()
    clin = clinical.set_index("participant")
    hs = HallmarkSets.load()
    proxies = _proliferation_proxies(rna, hs)
    cnv = D.load_cnv_target_genes(genes)
    ts_lookup = {
        (r["miRNA"], r["gene"]): r["ts_weight"]
        for _, r in ts_long.iterrows()
    } if not ts_long.empty else {}

    rows: List[dict] = []
    for arm in FAMILY_ARMS:
        if arm not in mirna.index:
            continue
        x_all = mirna.loc[arm]
        for gene in genes:
            if gene not in rna.index:
                continue
            y_all = rna.loc[gene]
            cn_row = cnv.loc[gene] if gene in cnv.index else None
            for scope, samples in _pam50_scope_iter(clinical):
                cols = x_all.index.intersection(y_all.index).intersection(clin.index)
                if samples is not None:
                    cols = cols.intersection(samples)
                cols = pd.Index(sorted(cols))
                if len(cols) < MIN_SCOPE_N:
                    continue
                ladder = _partial_ladder(y_all, x_all, clin, cols, proxies, target_cn=cn_row)
                rho = ladder.get("rho_CPE_HRD_CN", ladder.get("rho_CPE_HRD"))
                p = ladder.get("p_CPE_HRD_CN", ladder.get("p_CPE_HRD"))
                rows.append(
                    {
                        "predictor": arm,
                        "predictor_type": "single_arm",
                        "gene": gene,
                        "scope": scope,
                        "n": len(cols),
                        "ts_weight": ts_lookup.get((arm, gene), np.nan),
                        "neg_sig": bool(pd.notna(rho) and rho < 0 and pd.notna(p) and p < 0.05),
                        **ladder,
                    }
                )
    return pd.DataFrame(rows).sort_values(["gene", "scope", "predictor"]).reset_index(drop=True)


def _family_predictors(mirna: pd.DataFrame) -> Dict[str, pd.Series]:
    present = [a for a in FAMILY_ARMS if a in mirna.index]
    if not present:
        return {}
    m = mirna.loc[present].astype(float)
    z = zscore_rows(m).fillna(0.0)
    linear_sum = np.expm1(np.log(2.0) * m).clip(lower=0.0).sum(axis=0)
    return {
        "family_mean_log2rpm": m.mean(axis=0),
        "family_max_log2rpm": m.max(axis=0),
        "family_sum_log2rpm": np.log2(linear_sum + 1.0),
        "family_mean_z": z.mean(axis=0),
    }


def _family_partials(genes: Sequence[str]) -> pd.DataFrame:
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    clinical = D.load_clinical_strata()
    clin = clinical.set_index("participant")
    hs = HallmarkSets.load()
    proxies = _proliferation_proxies(rna, hs)
    cnv = D.load_cnv_target_genes(genes)
    predictors = _family_predictors(mirna)

    rows: List[dict] = []
    for pred_name, x_all in predictors.items():
        for gene in genes:
            if gene not in rna.index:
                continue
            y_all = rna.loc[gene]
            cn_row = cnv.loc[gene] if gene in cnv.index else None
            for scope, samples in _pam50_scope_iter(clinical):
                cols = x_all.index.intersection(y_all.index).intersection(clin.index)
                if samples is not None:
                    cols = cols.intersection(samples)
                cols = pd.Index(sorted(cols))
                if len(cols) < MIN_SCOPE_N:
                    continue
                ladder = _partial_ladder(y_all, x_all, clin, cols, proxies, target_cn=cn_row)
                rho = ladder.get("rho_CPE_HRD_CN", ladder.get("rho_CPE_HRD"))
                p = ladder.get("p_CPE_HRD_CN", ladder.get("p_CPE_HRD"))
                rows.append(
                    {
                        "predictor": pred_name,
                        "predictor_type": "family_aggregate",
                        "gene": gene,
                        "scope": scope,
                        "n": len(cols),
                        "neg_sig": bool(pd.notna(rho) and rho < 0 and pd.notna(p) and p < 0.05),
                        **ladder,
                    }
                )
    return pd.DataFrame(rows).sort_values(["gene", "scope", "predictor"]).reset_index(drop=True)


def _pressure_contributions(genes: Sequence[str]) -> pd.DataFrame:
    parts = []
    for expr_mode in PRESSURE_VARIANTS:
        tmp = compute_gene_pressure_contributions(
            genes,
            expr_mode=expr_mode,  # type: ignore[arg-type]
            target_norm=C.PRESSURE_TARGET_NORM,  # type: ignore[arg-type]
        )
        if not tmp.empty:
            tmp["is_family_arm"] = tmp["miRNA"].isin(FAMILY_ARMS)
            parts.append(tmp)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _ev_family_contrasts() -> pd.DataFrame:
    path = C.PLASMA_EV_DIR / "multi_cohort" / "ev_multi_cohort_arm_long.tsv"
    if not path.is_file():
        return pd.DataFrame()
    ev = pd.read_csv(path, sep="\t")
    return ev.loc[ev["arm"].isin(EV_RELATED_ARMS)].copy()


def _summarize_partials(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return rows
    rho_col = "rho_CPE_HRD_CN" if "rho_CPE_HRD_CN" in rows.columns else "rho_CPE_HRD"
    out = []
    for (predictor_type, predictor), sub in rows.groupby(["predictor_type", "predictor"]):
        cohort = sub.loc[sub["scope"].eq("cohort")]
        neg_raw = cohort.loc[cohort["neg_sig"]]
        neg_fdr = cohort.loc[cohort.get("neg_sig_cn_fdr", False)]
        neg_prolif = cohort.loc[cohort.get("neg_sig_prolif_cn_fdr", False)]
        out.append(
            {
                "predictor_type": predictor_type,
                "predictor": predictor,
                "n_genes_cohort": int(cohort["gene"].nunique()),
                "n_neg_raw_p05": int(neg_raw["gene"].nunique()),
                "n_neg_cn_fdr": int(neg_fdr["gene"].nunique()),
                "n_neg_prolif_cn_fdr": int(neg_prolif["gene"].nunique()),
                "median_neg_rho_fdr": float(neg_fdr[rho_col].median()) if not neg_fdr.empty else np.nan,
                "top_neg_genes_fdr": ";".join(
                    neg_fdr.sort_values(rho_col)["gene"].astype(str).head(10).tolist()
                ) if not neg_fdr.empty else "",
            }
        )
    return pd.DataFrame(out).sort_values(["predictor_type", "n_neg_cn_fdr"], ascending=[True, False])


def _share_comparison_table(contrib: pd.DataFrame, focus: Sequence[str]) -> pd.DataFrame:
    """Compact per-(gene, family arm) share comparison across modes + metrics."""
    if contrib.empty:
        return contrib
    sub = contrib.loc[
        contrib["gene"].isin(set(focus)) & contrib["miRNA"].isin(set(FAMILY_ARMS))
    ].copy()
    keep = [
        "gene", "miRNA", "expr_mode", "median_log2rpm", "sd_log2rpm",
        "mean_abs_share", "mean_pos_share", "global_abs_share", "global_signed_share",
    ]
    keep = [c for c in keep if c in sub.columns]
    return sub[keep].sort_values(["gene", "miRNA", "expr_mode"]).reset_index(drop=True)


def run(*, out_dir: Path = OUT_DIR, top_orphans: int = 25) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = _focus_genes(top_orphans=top_orphans)
    ts_long, ts_wide = _targetscan_family_scores(genes)

    var_tbl = _arm_variability(D.load_mirna_arms())

    var_reset = var_tbl.rename_axis("predictor").reset_index()
    arm_partials = _single_arm_partials(genes, ts_long)
    if not arm_partials.empty:
        arm_partials = arm_partials.merge(var_reset, on="predictor", how="left")
    fam_partials = _family_partials(genes)
    all_partials = pd.concat([arm_partials, fam_partials], ignore_index=True)
    all_partials = _add_fdr_and_flags(all_partials)
    # propagate FDR flags back to the single-arm split for its own export
    arm_partials = all_partials.loc[all_partials["predictor_type"].eq("single_arm")].copy()
    fam_partials = all_partials.loc[all_partials["predictor_type"].eq("family_aggregate")].copy()

    contrib = _pressure_contributions(genes)
    ev = _ev_family_contrasts()
    summary = _summarize_partials(all_partials)
    share_cmp = _share_comparison_table(contrib, ["PTEN", "MEOX2", "TGFBR2", "RAI2"])

    pd.Series(genes, name="gene").to_csv(out_dir / "focus_genes.tsv", sep="\t", index=False)
    ts_long.to_csv(out_dir / "family_targetscan_scores_long.tsv", sep="\t", index=False)
    ts_wide.to_csv(out_dir / "family_targetscan_scores_wide.tsv", sep="\t", index=False)
    var_tbl.rename_axis("arm").reset_index().to_csv(
        out_dir / "arm_variability.tsv", sep="\t", index=False
    )
    arm_partials.to_csv(out_dir / "family_single_arm_partials.tsv", sep="\t", index=False)
    fam_partials.to_csv(out_dir / "family_aggregate_partials.tsv", sep="\t", index=False)
    all_partials.to_csv(out_dir / "family_all_partials.tsv", sep="\t", index=False)
    contrib.to_csv(out_dir / "gene_mirna_pressure_contributions_by_variant.tsv", sep="\t", index=False)
    share_cmp.to_csv(out_dir / "family_share_metric_comparison.tsv", sep="\t", index=False)
    ev.to_csv(out_dir / "ev_family_arm_contrasts.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "family_partial_summary.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.analyses.mir301.mir301_family_depth",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "family_arms": list(FAMILY_ARMS),
        "ev_related_arms": list(EV_RELATED_ARMS),
        "pressure_variants": list(PRESSURE_VARIANTS),
        "share_metrics": ["mean_abs_share", "mean_pos_share", "global_abs_share", "global_signed_share"],
        "target_norm": C.PRESSURE_TARGET_NORM,
        "significance": {
            "fdr_alpha": FDR_ALPHA,
            "fdr_family": "BH within (predictor_type, predictor, scope) over genes",
            "headline_flag": "neg_sig_cn_fdr",
            "prolif_flag": "neg_sig_prolif_cn_fdr (E2F/G2M + CN adjusted)",
        },
        "variability_rule": {
            "scale": "log2(RPM+1) across cohort",
            "low_variability": f"sd<{LOW_VAR_SD} AND IQR<{LOW_VAR_IQR}",
        },
        "n_focus_genes": len(genes),
        "top_orphans_from_mir301a_depth": top_orphans,
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[mir301_family_depth] focus genes={len(genes)}")
    print("[mir301_family_depth] Partial summary (FDR-based):\n", summary.to_string(index=False))
    fam_var = var_tbl.loc[var_tbl.index.isin(FAMILY_ARMS)]
    print("\n[mir301_family_depth] Family-arm variability:\n", fam_var.to_string())
    return {
        "ts_long": ts_long,
        "ts_wide": ts_wide,
        "arm_partials": arm_partials,
        "family_partials": fam_partials,
        "contributions": contrib,
        "share_comparison": share_cmp,
        "variability": var_tbl,
        "ev": ev,
        "summary": summary,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--top-orphans", type=int, default=25)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, top_orphans=args.top_orphans)


if __name__ == "__main__":
    main()
