"""Combined miRNA-target anti-correlation at target-gene resolution.

For each target gene, combine all targeting miRNA MIMAT arms into the existing
evidence-weighted miRNA pressure score, then test whether that combined pressure
anti-correlates with the target gene's RNA expression.

This deliberately reuses ``pressure_build.compute_gene_pressure`` (``softmax_z_logrpm``
+ ``evidence_mass`` on M0 miRTarBase edges) and the same AGO-gating / partial
Spearman machinery as ``hallmark_interaction``, instead of reimplementing edge weighting.
Marginal and partial Pearson r are reported alongside Spearman as sensitivity columns.

Outputs under ``output/target_combined_anticorr/``:
  - ``target_combined_anticorr.tsv``
  - ``target_combined_anticorr_by_pam50.tsv`` (raw + CPE/HRD + target-CN partials)
  - ``target_combined_top_mirnas.tsv``
  - ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.hallmark_interaction import compute_gene_pressure
from mirna_hallmark.pressure_build import load_mirtar_edges, method_blurb
from mirna_hallmark.hallmark_sets import HallmarkSets


def _edge_table(genes: Sequence[str]) -> pd.DataFrame:
    return load_mirtar_edges(genes, resolve_arms=True)


def _edge_summaries(edges: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-gene edge counts and top miRNA contributors."""
    e = edges.copy()
    e["edge_weight"] = np.log1p(pd.to_numeric(e["evidence_score"], errors="coerce").fillna(0.0))
    counts = (
        e.groupby("gene")
        .agg(
            n_targeting_mirnas=("miRNA", "nunique"),
            n_edges=("miRNA", "size"),
            total_edge_weight=("edge_weight", "sum"),
        )
        .reset_index()
    )
    top = (
        e.sort_values(["gene", "edge_weight"], ascending=[True, False])
        .groupby("gene")
        .head(10)
        .loc[:, ["gene", "miRNA", "evidence_score", "edge_weight"]]
        .copy()
    )
    return counts, top


def _matrix_row(mat: pd.DataFrame, row_id: str) -> pd.Series:
    """Return one row as a Series, median-collapsing duplicated row labels."""
    row = mat.loc[row_id]
    if isinstance(row, pd.DataFrame):
        row = row.apply(pd.to_numeric, errors="coerce").median(axis=0)
    return pd.to_numeric(row, errors="coerce")


def _gene_target_cn(genes: Sequence[str]) -> pd.DataFrame:
    """gene x participant target CN matrix; empty if unavailable."""
    try:
        return D.load_cnv_target_genes(list(genes))
    except Exception as exc:
        print(f"[target_combined_anticorr] target CN unavailable; skipping CN-adjusted partials: {exc}")
        return pd.DataFrame()


def combined_target_anticorrelation(
    *,
    genes: Sequence[str] | None = None,
    include_gated: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (cohort rows, PAM50 rows, top miRNA contributors)."""
    hs = HallmarkSets.load()
    gene_list = list(genes or hs.universe)
    edges = _edge_table(gene_list)
    edge_counts, top_mirnas = _edge_summaries(edges)
    edge_counts = edge_counts.set_index("gene")

    print(f"[target_combined_anticorr] computing combined pressure for {len(gene_list):,} genes ...")
    pressure = compute_gene_pressure(gene_list)
    rna = D.load_rna()
    clinical = D.load_clinical_strata().set_index("participant")
    conf_cols = [c for c in C.CONFOUNDER_NUMERIC if c in clinical.columns]
    target_cn = _gene_target_cn(pressure.index)

    pressure_views: Dict[str, pd.DataFrame] = {"ungated": pressure}
    if include_gated:
        gate = compute_ago_gate(rna)
        shared = pressure.columns.intersection(gate.index)
        pressure_views["gated"] = pressure[shared].mul(gate.reindex(shared)["ago_gate"], axis=1)

    cohort_rows: List[dict] = []
    by_pam50_rows: List[dict] = []

    for view, mat in pressure_views.items():
        samples = mat.columns.intersection(rna.columns)
        for gene in mat.index.intersection(rna.index):
            p = mat.loc[gene, samples].rename("combined_pressure")
            e = _matrix_row(rna, gene).reindex(samples).rename("target_rna")
            clin_cov = D.augment_tcga_batch(clinical.reindex(samples)[conf_cols]) if conf_cols else None
            st = S.correlation_pair(e, p, clin_cov)
            if st["n"] < 20:
                continue

            prho_cn, pp_cn = np.nan, np.nan
            pr_cn_r, pr_cn_p = np.nan, np.nan
            if not target_cn.empty and gene in target_cn.index:
                cov = clin_cov.copy() if clin_cov is not None else pd.DataFrame(index=samples)
                cov["target_cn"] = _matrix_row(target_cn, gene).reindex(samples)
                st_cn = S.correlation_pair(e, p, cov)
                prho_cn = st_cn["partial_rho"]
                pp_cn = st_cn["partial_p"]
                pr_cn_r = st_cn["partial_pearson_r"]
                pr_cn_p = st_cn["partial_pearson_p"]

            ec = edge_counts.loc[gene] if gene in edge_counts.index else pd.Series(dtype=float)
            cohort_rows.append(
                {
                    "gene": gene,
                    "view": view,
                    "n": st["n"],
                    "n_targeting_mirnas": int(ec.get("n_targeting_mirnas", 0)),
                    "n_edges": int(ec.get("n_edges", 0)),
                    "total_edge_weight": round(float(ec.get("total_edge_weight", np.nan)), 4)
                    if pd.notna(ec.get("total_edge_weight", np.nan))
                    else np.nan,
                    "spearman_rho": round(float(st["spearman_rho"]), 4)
                    if np.isfinite(st["spearman_rho"])
                    else np.nan,
                    "spearman_p": st["spearman_p"],
                    "pearson_r": round(float(st["pearson_r"]), 4)
                    if np.isfinite(st["pearson_r"])
                    else np.nan,
                    "pearson_p": st["pearson_p"],
                    "partial_rho": round(float(st["partial_rho"]), 4)
                    if np.isfinite(st["partial_rho"])
                    else np.nan,
                    "partial_p": st["partial_p"],
                    "partial_pearson_r": round(float(st["partial_pearson_r"]), 4)
                    if np.isfinite(st["partial_pearson_r"])
                    else np.nan,
                    "partial_pearson_p": st["partial_pearson_p"],
                    "partial_rho_plus_target_cn": round(float(prho_cn), 4)
                    if np.isfinite(prho_cn)
                    else np.nan,
                    "partial_p_plus_target_cn": pp_cn,
                    "partial_pearson_r_plus_target_cn": round(float(pr_cn_r), 4)
                    if np.isfinite(pr_cn_r)
                    else np.nan,
                    "partial_pearson_p_plus_target_cn": pr_cn_p,
                    "mean_combined_pressure": round(float(p.mean()), 4),
                    "mean_target_rna": round(float(e.mean()), 4),
                }
            )

            merged = pd.concat(
                [p, e, clinical[["PAM50_final"] + conf_cols]],
                axis=1,
                join="inner",
            ).dropna(subset=["combined_pressure", "target_rna", "PAM50_final"])
            gene_cn = (
                _matrix_row(target_cn, gene).rename("target_cn")
                if not target_cn.empty and gene in target_cn.index
                else None
            )
            for pam50, sub in merged.groupby("PAM50_final"):
                if len(sub) < 20:
                    continue
                sub_cov = D.augment_tcga_batch(sub[conf_cols]) if conf_cols else None
                st_sub = S.correlation_pair(
                    sub["target_rna"], sub["combined_pressure"], sub_cov
                )
                prho_cn, pp_cn = np.nan, np.nan
                pr_cn_r, pr_cn_p = np.nan, np.nan
                if gene_cn is not None:
                    cov = (sub_cov.copy() if sub_cov is not None else pd.DataFrame(index=sub.index))
                    cov["target_cn"] = gene_cn.reindex(sub.index)
                    st_cn = S.correlation_pair(sub["target_rna"], sub["combined_pressure"], cov)
                    prho_cn = st_cn["partial_rho"]
                    pp_cn = st_cn["partial_p"]
                    pr_cn_r = st_cn["partial_pearson_r"]
                    pr_cn_p = st_cn["partial_pearson_p"]
                by_pam50_rows.append(
                    {
                        "gene": gene,
                        "view": view,
                        "pam50": str(pam50),
                        "n": int(st_sub["n"]),
                        "spearman_rho": round(float(st_sub["spearman_rho"]), 4),
                        "spearman_p": st_sub["spearman_p"],
                        "pearson_r": round(float(st_sub["pearson_r"]), 4),
                        "pearson_p": st_sub["pearson_p"],
                        "partial_rho": round(float(st_sub["partial_rho"]), 4)
                        if np.isfinite(st_sub["partial_rho"])
                        else np.nan,
                        "partial_p": st_sub["partial_p"],
                        "partial_pearson_r": round(float(st_sub["partial_pearson_r"]), 4)
                        if np.isfinite(st_sub["partial_pearson_r"])
                        else np.nan,
                        "partial_pearson_p": st_sub["partial_pearson_p"],
                        "partial_rho_plus_target_cn": round(float(prho_cn), 4)
                        if np.isfinite(prho_cn)
                        else np.nan,
                        "partial_p_plus_target_cn": pp_cn,
                        "partial_pearson_r_plus_target_cn": round(float(pr_cn_r), 4)
                        if np.isfinite(pr_cn_r)
                        else np.nan,
                        "partial_pearson_p_plus_target_cn": pr_cn_p,
                    }
                )

    cohort = pd.DataFrame(cohort_rows)
    if not cohort.empty:
        for view, idx in cohort.groupby("view").groups.items():
            cohort.loc[idx, "spearman_q"] = S.bh_fdr(cohort.loc[idx, "spearman_p"].values)
            cohort.loc[idx, "pearson_q"] = S.bh_fdr(cohort.loc[idx, "pearson_p"].values)
            cohort.loc[idx, "partial_q"] = S.bh_fdr(
                cohort.loc[idx, "partial_p"].fillna(1.0).values
            )
            cohort.loc[idx, "partial_pearson_q"] = S.bh_fdr(
                cohort.loc[idx, "partial_pearson_p"].fillna(1.0).values
            )
            cohort.loc[idx, "partial_q_plus_target_cn"] = S.bh_fdr(
                cohort.loc[idx, "partial_p_plus_target_cn"].fillna(1.0).values
            )
            cohort.loc[idx, "partial_pearson_q_plus_target_cn"] = S.bh_fdr(
                cohort.loc[idx, "partial_pearson_p_plus_target_cn"].fillna(1.0).values
            )
        cohort = cohort.sort_values(["view", "spearman_rho"])

    by_pam50 = pd.DataFrame(by_pam50_rows)
    if not by_pam50.empty:
        for key, idx in by_pam50.groupby(["view", "pam50"]).groups.items():
            by_pam50.loc[idx, "spearman_q"] = S.bh_fdr(by_pam50.loc[idx, "spearman_p"].values)
            by_pam50.loc[idx, "pearson_q"] = S.bh_fdr(by_pam50.loc[idx, "pearson_p"].values)
            by_pam50.loc[idx, "partial_q"] = S.bh_fdr(
                by_pam50.loc[idx, "partial_p"].fillna(1.0).values
            )
            by_pam50.loc[idx, "partial_pearson_q"] = S.bh_fdr(
                by_pam50.loc[idx, "partial_pearson_p"].fillna(1.0).values
            )
            by_pam50.loc[idx, "partial_q_plus_target_cn"] = S.bh_fdr(
                by_pam50.loc[idx, "partial_p_plus_target_cn"].fillna(1.0).values
            )
            by_pam50.loc[idx, "partial_pearson_q_plus_target_cn"] = S.bh_fdr(
                by_pam50.loc[idx, "partial_pearson_p_plus_target_cn"].fillna(1.0).values
            )
        by_pam50 = by_pam50.sort_values(["view", "pam50", "spearman_rho"])

    return cohort, by_pam50, top_mirnas


def run(
    *,
    out_dir: Path = C.TARGET_COMBINED_ANTICORR_DIR,
    include_gated: bool = True,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort, by_pam50, top_mirnas = combined_target_anticorrelation(include_gated=include_gated)

    cohort.to_csv(out_dir / "target_combined_anticorr.tsv", sep="\t", index=False)
    by_pam50.to_csv(out_dir / "target_combined_anticorr_by_pam50.tsv", sep="\t", index=False)
    top_mirnas.to_csv(out_dir / "target_combined_top_mirnas.tsv", sep="\t", index=False)

    sig = cohort.loc[
        (cohort["view"] == "gated")
        & (cohort["spearman_rho"] < 0)
        & (cohort["spearman_q"] < C.FDR_ALPHA)
    ] if not cohort.empty else pd.DataFrame()
    sig_cn = cohort.loc[
        (cohort["view"] == "gated")
        & (cohort["partial_rho_plus_target_cn"] < 0)
        & (cohort["partial_q_plus_target_cn"] < C.FDR_ALPHA)
    ] if not cohort.empty else pd.DataFrame()
    by_pam50_gated = by_pam50.loc[by_pam50["view"] == "gated"] if not by_pam50.empty else pd.DataFrame()
    sig_by_pam50 = (
        by_pam50_gated.loc[
            (by_pam50_gated["partial_rho_plus_target_cn"] < 0)
            & (by_pam50_gated["partial_q_plus_target_cn"] < C.FDR_ALPHA)
        ]
        if not by_pam50_gated.empty
        else pd.DataFrame()
    )
    sig_by_pam50_counts = (
        sig_by_pam50.groupby("pam50")["gene"].nunique().to_dict() if not sig_by_pam50.empty else {}
    )

    manifest = {
        "module": "mirna_hallmark.target_combined_anticorr",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pressure": method_blurb(),
        "edge_source": str(C.MIRTAR_HALLMARK_SUMMARY),
        "min_evidence": C.PRESSURE_MIN_EVIDENCE,
        "weight_mode": C.PRESSURE_WEIGHT_MODE,
        "include_gated": include_gated,
        "within_subtype_partials": {
            "confounders": list(C.CONFOUNDER_NUMERIC),
            "plus_target_cn": True,
            "note": "PAM50 excluded from within-subtype partials (stratum-fixed)",
        },
        "n_rows": int(len(cohort)),
        "n_genes": int(cohort["gene"].nunique()) if not cohort.empty else 0,
        "n_gated_negative_fdr": int(len(sig)),
        "n_gated_negative_fdr_plus_target_cn": int(len(sig_cn)),
        "n_gated_negative_fdr_plus_target_cn_by_pam50": sig_by_pam50_counts,
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(
        f"[target_combined_anticorr] wrote {len(cohort):,} cohort rows; "
        f"gated negative q<{C.FDR_ALPHA}: {len(sig):,}; "
        f"CN-adjusted: {len(sig_cn):,}; "
        f"within-PAM50 CN-adjusted: {sig_by_pam50_counts}"
    )
    return cohort


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.TARGET_COMBINED_ANTICORR_DIR)
    ap.add_argument("--no-gate", action="store_true", help="Only compute ungated pressure")
    args = ap.parse_args()
    run(out_dir=args.out_dir, include_gated=not args.no_gate)


if __name__ == "__main__":
    main()
