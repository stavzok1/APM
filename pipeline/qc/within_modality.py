from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from pipeline.qc.base import Finding, Level


def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _try_load_rna_subset(*, sample_ids: List[str], genes: List[str]) -> pd.DataFrame:
    try:
        from pipeline.config import PATHS
        from pipeline.RNA_exp.signatures import read_tpm_wide_subset

        return read_tpm_wide_subset(PATHS.rna_expression, genes=genes, sample_cols=sample_ids, gene_col=PATHS.rna_gene_col)
    except Exception:
        return pd.DataFrame()


def compute_within_modality_sanity(
    *,
    ctx,
    sample_ids: List[str],
    out_dir: Path | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """
    Tier-2 within-modality sanity checks (warning-oriented).

    Returns:
    - per_sample: one row per sample with summary stats
    - per_gene: empty for now (reserved)
    - summary: dict with aggregate stats + file presence
    """
    per_gene = pd.DataFrame()
    rows: List[Dict[str, object]] = []
    summary: Dict[str, object] = {}

    # --- Methylation promoter beta distributions (gene-level aggregation matrix) ---
    gene_meth_path = Path(ctx.methylation_working_dir) / "cohort" / "gene_meth_matrix.csv"
    summary["methylation_gene_matrix_exists"] = gene_meth_path.exists()
    if gene_meth_path.exists():
        df = pd.read_csv(gene_meth_path, low_memory=False)
        if "gene_name" in df.columns:
            df = df.set_index("gene_name")
        if "meta" in df.columns:
            df = df.drop(columns=["meta"])
        keep = [s for s in sample_ids if s in df.columns]
        # per-sample mean/median beta across genes
        for sid in keep:
            v = _num(df[sid])
            rows.append(
                {
                    "sample_id": sid,
                    "meth_prom_beta_mean": float(v.mean()) if v.notna().any() else np.nan,
                    "meth_prom_beta_median": float(v.median()) if v.notna().any() else np.nan,
                    "meth_prom_beta_p05": float(v.quantile(0.05)) if v.notna().any() else np.nan,
                    "meth_prom_beta_p95": float(v.quantile(0.95)) if v.notna().any() else np.nan,
                }
            )
        summary["methylation_gene_matrix_n_genes"] = int(len(df))
        summary["methylation_gene_matrix_n_samples"] = int(len(keep))

    per_sample = pd.DataFrame(rows).drop_duplicates(subset=["sample_id"])

    # --- RNA sanity: log2(TPM+1) summary stats on a small, stable gene set ---
    summary["rna_subset_stats"] = []
    try:
        from pipeline.RNA_exp.signatures import GeneSets

        # Use APM class I gene set (small but always present in this project).
        genes = list(GeneSets().apm_class_i)
        expr = _try_load_rna_subset(sample_ids=sample_ids, genes=genes)
        if not expr.empty:
            # Per-sample mean/std across the gene set
            for sid in [s for s in sample_ids if s in expr.columns]:
                v = _num(expr[sid])
                summary["rna_subset_stats"].append(
                    {
                        "sample_id": sid,
                        "n_genes": int(v.notna().sum()),
                        "mean_log2tpm1": float(v.mean()) if v.notna().any() else np.nan,
                        "std_log2tpm1": float(v.std(ddof=0)) if v.notna().any() else np.nan,
                    }
                )
    except Exception:
        pass

    # --- CNV sanity: presence of per-sample gene calls + summary scale ---
    # (We avoid loading full segment tables here; keep it cheap and cohort-slice friendly.)
    cnv_stats = []
    for sid in sample_ids:
        p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls_ascat3.csv"
        if not p.exists():
            p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls.csv"
        if not p.exists():
            continue
        try:
            gdf = pd.read_csv(p, low_memory=False)
            num = None
            for cand in ("cnv_log2", "segment_mean", "log2_copy_ratio", "log2_ratio"):
                if cand in gdf.columns:
                    num = cand
                    break
            if num is None and "copy_number" in gdf.columns:
                cn = pd.to_numeric(gdf["copy_number"], errors="coerce")
                with np.errstate(divide="ignore", invalid="ignore"):
                    v = np.log2(cn / 2.0)
            else:
                v = pd.to_numeric(gdf[num], errors="coerce") if num is not None else pd.Series([], dtype=float)
            v = pd.Series(v).dropna()
            cnv_stats.append(
                {
                    "sample_id": sid,
                    "n_genes": int(len(v)),
                    "mean": float(v.mean()) if len(v) else np.nan,
                    "std": float(v.std(ddof=0)) if len(v) else np.nan,
                }
            )
        except Exception:
            continue
    summary["cnv_gene_calls_stats"] = cnv_stats

    # --- RPPA sanity: panel scores should be roughly z-scored (mean~0, std~1) ---
    p_panel = Path(ctx.rppa_output_dir) / "panel_scores.csv"
    summary["rppa_panel_scores_exists"] = p_panel.exists()
    if p_panel.exists():
        ps = pd.read_csv(p_panel)
        if "sample_id" in ps.columns:
            ps = ps.set_index("sample_id")
        # evaluate a few columns if present
        cols = [c for c in ps.columns if c.startswith("IFN_") or c.startswith("DDR_") or c.startswith("cGAS_")]
        if not cols:
            cols = list(ps.columns[: min(8, ps.shape[1])])
        stats = []
        for c in cols:
            v = _num(ps[c])
            stats.append(
                {
                    "rppa_col": c,
                    "mean": float(v.mean()) if v.notna().any() else np.nan,
                    "std": float(v.std(ddof=0)) if v.notna().any() else np.nan,
                    "nn_rate": float(v.notna().mean()),
                }
            )
        summary["rppa_panel_score_stats"] = stats

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        per_sample.to_csv(out_dir / "within_modality__per_sample.csv", index=False)
        pd.DataFrame(summary.get("rppa_panel_score_stats", [])).to_csv(
            out_dir / "within_modality__rppa_panel_score_stats.csv", index=False
        )
        pd.DataFrame(summary.get("rna_subset_stats", [])).to_csv(
            out_dir / "within_modality__rna_subset_stats.csv", index=False
        )
        pd.DataFrame(summary.get("cnv_gene_calls_stats", [])).to_csv(
            out_dir / "within_modality__cnv_gene_calls_stats.csv", index=False
        )

    return per_sample, per_gene, summary


def assert_within_modality_sanity(summary: Dict[str, object]) -> List[Finding]:
    """
    Warning-oriented assertions. We avoid hard failures here because outliers can be real.
    """
    findings: List[Finding] = []

    if not summary.get("methylation_gene_matrix_exists", False):
        findings.append(
            Finding(
                check="within_modality.methylation_gene_matrix_present",
                level=Level.WARN,
                message="Missing methylation cohort gene matrix; methylation sanity checks skipped.",
                context={},
                hint="Expected at data/Methylation/cohort/gene_meth_matrix.csv (via scratch context).",
            )
        )

    if not summary.get("rppa_panel_scores_exists", False):
        findings.append(
            Finding(
                check="within_modality.rppa_panel_scores_present",
                level=Level.WARN,
                message="Missing RPPA panel_scores.csv; RPPA sanity checks skipped.",
                context={},
                hint="Expected under data/rppa/processed/ or ctx.rppa_output_dir.",
            )
        )
        return findings

    # RNA subset sanity: only catch obvious breakage (e.g., all-NaN or near-constant).
    rna_stats = summary.get("rna_subset_stats") or []
    if not rna_stats:
        findings.append(
            Finding(
                check="within_modality.rna_subset_present",
                level=Level.WARN,
                message="Missing RNA subset stats; RNA within-modality checks skipped.",
                context={},
                hint="If this is unexpected, check PATHS.rna_expression and scratch sample IDs.",
            )
        )
    else:
        df = pd.DataFrame(rna_stats)
        if "std_log2tpm1" in df.columns:
            med_std = float(pd.to_numeric(df["std_log2tpm1"], errors="coerce").median())
            if med_std != med_std or med_std < 0.05:
                findings.append(
                    Finding(
                        check="within_modality.rna_subset_variance",
                        level=Level.WARN,
                        message=f"RNA subset per-sample std dev is very low (median std={med_std:.3f}).",
                        context={"median_std": med_std},
                        hint="Could indicate wrong transform, wrong join, or a constant matrix slice.",
                    )
                )

    # CNV gene calls sanity: warn if missing for most samples.
    cnv_stats = summary.get("cnv_gene_calls_stats") or []
    if not cnv_stats:
        findings.append(
            Finding(
                check="within_modality.cnv_gene_calls_present",
                level=Level.WARN,
                message="Missing CNV gene calls for cohort; CNV within-modality checks skipped.",
                context={},
                hint="Expected per-sample *_cnv_gene_calls*.csv under ctx.cnv_gene_tables_dir.",
            )
        )

    stats = summary.get("rppa_panel_score_stats") or []
    for s in stats:
        try:
            std = float(s.get("std"))
            mean = float(s.get("mean"))
            col = str(s.get("rppa_col"))
        except Exception:
            continue
        # Panel scores are typically *means of per-target z-scores* (see RPPA panels),
        # so their std dev is not expected to be ~1.0. We only want to catch gross
        # degeneracy (nearly constant) or wildly inflated scaling.
        if not (0.15 <= std <= 2.5):
            findings.append(
                Finding(
                    check="within_modality.rppa_panel_score_scale",
                    level=Level.WARN,
                    message=f"RPPA panel score {col} std dev is {std:.2f} (expected ~1).",
                    context={"col": col, "std": std, "mean": mean},
                    hint=(
                        "Panel scores are means of z-scored targets; std can be <1. "
                        "Very low std may indicate a homogeneous cohort slice or upstream scaling issues."
                    ),
                )
            )
        if abs(mean) > 0.3:
            findings.append(
                Finding(
                    check="within_modality.rppa_panel_score_centering",
                    level=Level.WARN,
                    message=f"RPPA panel score {col} mean is {mean:.2f} (expected ~0).",
                    context={"col": col, "std": std, "mean": mean},
                )
            )

    return findings

