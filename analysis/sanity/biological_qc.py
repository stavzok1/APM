from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS
from pipeline.sample_ids import normalize_tcga_id
from pipeline.RNA_exp.signatures import (
    GeneSets,
    compute_gene_set_scores_from_tpm,
    read_tpm_wide_subset,
)

from analysis.sanity._io import (
    ensure_out_dir,
    load_scratch_context,
    load_sample_ids,
    utc_run_id,
)


def _load_gene_meth_matrix(path: Path, *, genes: Sequence[str], sample_ids: Sequence[str]) -> pd.DataFrame:
    """
    gene_meth_matrix.csv is genes x samples + 'meta' column.
    """
    df = pd.read_csv(path, low_memory=False)
    if "gene_name" in df.columns:
        df = df.set_index("gene_name")
    # drop meta if present
    if "meta" in df.columns:
        df = df.drop(columns=["meta"])
    # subset
    keep = [s for s in sample_ids if s in df.columns]
    out = df.reindex(index=[g for g in genes if g in df.index], columns=keep)
    for c in keep:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def _load_cnv_gene_calls(path: Path, *, genes: Sequence[str]) -> pd.DataFrame:
    """
    Expect a per-sample gene table with `gene_name` and copy-number fields.

    Returns a small DF indexed by gene with columns:
    - copy_number: float
    - cn_state: int (0/1/2/3 where 3 means 3+)
    - cnv_log2: float (log2(copy_number/2))
    """
    df = pd.read_csv(path, low_memory=False)
    # heuristics: find gene column + best numeric column
    gcol = "gene_name" if "gene_name" in df.columns else ("gene" if "gene" in df.columns else None)
    if gcol is None:
        return pd.DataFrame()
    df[gcol] = df[gcol].astype(str)
    sub = df[df[gcol].isin(list(genes))].copy()

    # ASCAT3 gene-level tables: prefer copy_number and convert to log2 ratio vs diploid (2.0)
    if "copy_number" in sub.columns:
        cn = pd.to_numeric(sub["copy_number"], errors="coerce")
        # avoid divide-by-zero and log2 of 0
        with np.errstate(divide="ignore", invalid="ignore"):
            log2r = np.log2(cn / 2.0)
        cn_state = pd.Series(pd.NA, index=sub.index)
        ok = cn.notna()
        if ok.any():
            # state: round to nearest integer, clip to [0, 3] where 3 means 3+
            st = np.rint(cn[ok].astype(float)).astype(int)
            st = np.clip(st, 0, 3)
            cn_state.loc[ok] = st
        out = sub[[gcol]].copy()
        out["copy_number"] = cn
        out["cn_state"] = cn_state
        out["cnv_log2"] = log2r
        return out.drop_duplicates(subset=[gcol]).set_index(gcol)[["copy_number", "cn_state", "cnv_log2"]]

    # try common fields
    num_col = None
    for cand in ("segment_mean", "log2_copy_ratio", "log2_ratio", "cnv_log2", "segment_mean_log2"):
        if cand in sub.columns:
            num_col = cand
            break
    if num_col is None:
        # fallback: first numeric-ish
        for c in sub.columns:
            if c == gcol:
                continue
            if pd.api.types.is_numeric_dtype(sub[c]) or c.endswith("_log2") or c.endswith("_mean"):
                num_col = c
                break
    if num_col is None:
        return pd.DataFrame()
    sub[num_col] = pd.to_numeric(sub[num_col], errors="coerce")
    out = sub[[gcol, num_col]].drop_duplicates(subset=[gcol]).set_index(gcol)
    out.columns = ["cnv_log2"]
    # no state info available in this fallback mode
    out["copy_number"] = np.nan
    out["cn_state"] = pd.NA
    return out[["copy_number", "cn_state", "cnv_log2"]]


def _spearman_r(x: pd.Series, y: pd.Series) -> float:
    xr = x.rank(method="average")
    yr = y.rank(method="average")
    m = xr.notna() & yr.notna()
    if int(m.sum()) < 8:
        return float("nan")
    return float(np.corrcoef(xr[m].astype(float), yr[m].astype(float))[0, 1])


def _to_tcga_sample_id(x: str) -> str:
    """
    Normalize any TCGA-ish id to a TCGA sample id (e.g. TCGA-XX-YYYY-01).
    """
    tid = normalize_tcga_id(str(x))
    return str(tid.sample or "").strip()


def _group_summary(
    df: pd.DataFrame,
    *,
    group_col: str,
    value_cols: List[str],
    min_n: int = 3,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for g, sub in df.groupby(group_col, dropna=False):
        rec: Dict[str, object] = {group_col: g, "n": int(len(sub))}
        if int(len(sub)) < int(min_n):
            rows.append(rec)
            continue
        for c in value_cols:
            if c not in sub.columns:
                continue
            v = pd.to_numeric(sub[c], errors="coerce")
            rec[f"{c}__mean"] = float(v.mean()) if v.notna().any() else np.nan
            rec[f"{c}__median"] = float(v.median()) if v.notna().any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["n", group_col], ascending=[False, True])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch-json", type=str, required=True)
    args = ap.parse_args()

    ctx = load_scratch_context(Path(args.scratch_json))
    run_id = utc_run_id("biological_qc")
    out_dir = ensure_out_dir(run_id)

    sample_ids = load_sample_ids(ctx)
    if not sample_ids:
        raise SystemExit(
            "No sample_ids found. Expected either an existing methylation subset manifest "
            f"({ctx.methylation_subset_manifest}) or a prepare log next to {ctx.scratch_json}."
        )
    (out_dir / "sample_ids.txt").write_text("\n".join(sample_ids) + "\n", encoding="utf-8")

    gs = GeneSets()
    key_genes = list(gs.apm_class_i)

    # --- RNA signatures + key genes ---
    sig = compute_gene_set_scores_from_tpm(
        PATHS.rna_expression,
        sample_ids=sample_ids,
        gene_col=PATHS.rna_gene_col,
    )
    # signatures: returned as samples x scores with index.name="sample_id"
    sig = sig.reset_index()
    sig.to_csv(out_dir / "rna_signatures.csv", index=False)

    expr_key = read_tpm_wide_subset(
        PATHS.rna_expression,
        genes=key_genes,
        sample_cols=sample_ids,
        gene_col=PATHS.rna_gene_col,
    )
    expr_key.to_csv(out_dir / "rna_apm_key_genes_tpm.csv")

    # --- Methylation promoter beta ---
    gene_meth_path = ctx.methylation_working_dir / "cohort" / "gene_meth_matrix.csv"
    if gene_meth_path.exists():
        meth_key = _load_gene_meth_matrix(gene_meth_path, genes=key_genes, sample_ids=sample_ids)
        meth_key.to_csv(out_dir / "methylation_promoter_beta_key_genes.csv")
    else:
        meth_key = pd.DataFrame()

    # --- CNV (gene-level copy number) for key genes ---
    cnv_rows: List[Dict[str, object]] = []
    cnv_mat: Dict[str, pd.Series] = {}
    cn_state_mat: Dict[str, pd.Series] = {}
    cn_copy_mat: Dict[str, pd.Series] = {}
    for sid in sample_ids:
        p = ctx.cnv_gene_tables_dir / f"{sid}_cnv_gene_calls_ascat3.csv"
        if not p.exists():
            p = ctx.cnv_gene_tables_dir / f"{sid}_cnv_gene_calls.csv"
        if not p.exists():
            cnv_rows.append({"sample_id": sid, "cnv_path": str(p), "exists": False})
            continue
        cnv_rows.append({"sample_id": sid, "cnv_path": str(p), "exists": True})
        gdf = _load_cnv_gene_calls(p, genes=key_genes)
        if not gdf.empty:
            if "cnv_log2" in gdf.columns:
                cnv_mat[sid] = gdf["cnv_log2"]
            if "cn_state" in gdf.columns:
                cn_state_mat[sid] = gdf["cn_state"]
            if "copy_number" in gdf.columns:
                cn_copy_mat[sid] = gdf["copy_number"]
    pd.DataFrame(cnv_rows).to_csv(out_dir / "cnv_gene_calls_paths.csv", index=False)
    if cnv_mat:
        cnv_df = pd.DataFrame(cnv_mat).reindex(index=key_genes)
        cnv_df.to_csv(out_dir / "cnv_log2_key_genes.csv")
    else:
        cnv_df = pd.DataFrame()
    if cn_state_mat:
        cn_state_df = pd.DataFrame(cn_state_mat).reindex(index=key_genes)
        cn_state_df.to_csv(out_dir / "cnv_state_key_genes.csv")
    else:
        cn_state_df = pd.DataFrame()
    if cn_copy_mat:
        cn_copy_df = pd.DataFrame(cn_copy_mat).reindex(index=key_genes)
        cn_copy_df.to_csv(out_dir / "cnv_copy_number_key_genes.csv")
    else:
        cn_copy_df = pd.DataFrame()

    # --- Simple sanity associations (no plots; write correlations) ---
    # 1) CNV gain vs expression (per gene, across samples)
    if not cnv_df.empty and not expr_key.empty:
        corr_rows = []
        for g in key_genes:
            if g not in cnv_df.index or g not in expr_key.index:
                continue
            x = pd.to_numeric(cnv_df.loc[g], errors="coerce")
            y = pd.to_numeric(expr_key.loc[g], errors="coerce")
            m = x.notna() & y.notna()
            if int(m.sum()) < 8:
                continue
            r = float(np.corrcoef(x[m].astype(float), y[m].astype(float))[0, 1])
            corr_rows.append({"gene": g, "n": int(m.sum()), "pearson_r_cnv_vs_expr": r})
        pd.DataFrame(corr_rows).to_csv(out_dir / "corr_cnv_vs_expression_key_genes.csv", index=False)

    # 1b) CNV state vs expression (0/1/2/3+ bins)
    if not cn_state_df.empty and not expr_key.empty:
        rows: List[Dict[str, object]] = []
        low_tpm_thresh = 0.1
        for g in key_genes:
            if g not in cn_state_df.index or g not in expr_key.index:
                continue
            st = pd.to_numeric(cn_state_df.loc[g], errors="coerce")
            tpm = pd.to_numeric(expr_key.loc[g], errors="coerce")
            # per-state summaries
            rec: Dict[str, object] = {"gene": g}
            for s in (0, 1, 2, 3):
                m = (st == s) & tpm.notna()
                rec[f"n_cn{s}"] = int(m.sum())
                rec[f"median_tpm_cn{s}"] = float(tpm[m].median()) if int(m.sum()) else np.nan
                rec[f"mean_tpm_cn{s}"] = float(tpm[m].mean()) if int(m.sum()) else np.nan
            # deletion check: CN=0 should have very low expression
            m0 = (st == 0) & tpm.notna()
            rec["frac_low_tpm_given_cn0"] = (
                float((tpm[m0] <= low_tpm_thresh).mean()) if int(m0.sum()) else np.nan
            )
            # monotonic-ish association measure: Spearman(copy_number, TPM) if available
            if not cn_copy_df.empty and g in cn_copy_df.index:
                cn = pd.to_numeric(cn_copy_df.loc[g], errors="coerce")
                rec["spearman_r_copy_number_vs_tpm"] = _spearman_r(cn, tpm)
            else:
                rec["spearman_r_copy_number_vs_tpm"] = np.nan
            rows.append(rec)
        pd.DataFrame(rows).to_csv(out_dir / "cnv_state_vs_expression_key_genes.csv", index=False)

    # 2) Promoter methylation vs expression (per gene)
    if not meth_key.empty and not expr_key.empty:
        corr_rows = []
        for g in key_genes:
            if g not in meth_key.index or g not in expr_key.index:
                continue
            x = pd.to_numeric(meth_key.loc[g], errors="coerce")  # beta
            y = pd.to_numeric(expr_key.loc[g], errors="coerce")  # tpm
            m = x.notna() & y.notna()
            if int(m.sum()) < 8:
                continue
            r = float(np.corrcoef(x[m].astype(float), y[m].astype(float))[0, 1])
            corr_rows.append({"gene": g, "n": int(m.sum()), "pearson_r_beta_vs_expr": r})
        pd.DataFrame(corr_rows).to_csv(out_dir / "corr_promoter_beta_vs_expression_key_genes.csv", index=False)

    # --- Subtype comparisons: pipeline "advanced" table + original Thorsson table ---
    # Normalize to TCGA sample ids (…-01) so `…-01A` matches `…-01`.
    sig2 = sig.copy()
    if "sample_id" not in sig2.columns:
        raise SystemExit("rna_signatures missing sample_id")
    sig2["tcga_sample_id"] = sig2["sample_id"].map(_to_tcga_sample_id)

    sig_cols = [c for c in ("APM_classI_mean", "CD8_mean", "NK_mean", "IFNG_mean", "CYT") if c in sig2.columns]

    # 1) Pipeline-produced BRCA immune advanced table (PAM50 + thornsson_subtype)
    adv_path = Path("annotations") / "BRCA_immune_subtypes_advanced.tsv"
    if adv_path.exists():
        adv = pd.read_csv(adv_path, sep="\t", low_memory=False)
        if "sample_id" in adv.columns:
            adv["tcga_sample_id"] = adv["sample_id"].map(_to_tcga_sample_id)
        else:
            adv["tcga_sample_id"] = ""
        keep_cols = [
            c
            for c in (
                "PAM50_final",
                "PAM50_recomputed",
                "PAM50Call_RNAseq",
                "thornsson_subtype",
                "CPE",
            )
            if c in adv.columns
        ]
        adv_small = adv[["tcga_sample_id", *keep_cols]].drop_duplicates(subset=["tcga_sample_id"])
        joined_adv = sig2.merge(adv_small, on="tcga_sample_id", how="left")
        joined_adv.to_csv(out_dir / "signatures_joined_brca_immune_advanced.csv", index=False)

        if "thornsson_subtype" in joined_adv.columns:
            _group_summary(joined_adv, group_col="thornsson_subtype", value_cols=sig_cols).to_csv(
                out_dir / "signatures_by_thornsson_subtype_pipeline_table.csv", index=False
            )
        if "PAM50_final" in joined_adv.columns:
            _group_summary(joined_adv, group_col="PAM50_final", value_cols=sig_cols).to_csv(
                out_dir / "signatures_by_pam50_pipeline_table.csv", index=False
            )

    # 2) Original Thorsson table: immune subtype + pathway/program scores + cell fractions
    th_path = Path("annotations") / "Thornsson_immune_table.tsv"
    if th_path.exists():
        th = pd.read_csv(th_path, sep="\t", low_memory=False)
        if "TCGA Participant Barcode" in th.columns:
            th["tcga_sample_id"] = th["TCGA Participant Barcode"].map(_to_tcga_sample_id)
        else:
            th["tcga_sample_id"] = ""

        th_keep = [
            c
            for c in (
                "Immune Subtype",
                "TCGA Subtype",
                "PAM50_final",
                "IFN-gamma Response",
                "TGF-beta Response",
                "Wound Healing",
                "Proliferation",
                "Macrophage Regulation",
                "CTA Score",
                "Leukocyte Fraction",
                "Lymphocyte Infiltration Signature Score",
                "T Cells CD8",
                "T Cells CD4 Memory Activated",
                "NK Cells Activated",
            )
            if c in th.columns
        ]
        th_small = th[["tcga_sample_id", *th_keep]].drop_duplicates(subset=["tcga_sample_id"])
        joined_th = sig2.merge(th_small, on="tcga_sample_id", how="left")
        joined_th.to_csv(out_dir / "signatures_joined_thornsson_original.csv", index=False)

        extra_cols = [c for c in th_keep if c not in ("Immune Subtype", "TCGA Subtype", "PAM50_final")]
        val_cols = sig_cols + extra_cols

        if "Immune Subtype" in joined_th.columns:
            _group_summary(joined_th, group_col="Immune Subtype", value_cols=val_cols).to_csv(
                out_dir / "signatures_and_thornsson_scores_by_immune_subtype_original.csv", index=False
            )
        if "PAM50_final" in joined_th.columns:
            _group_summary(joined_th, group_col="PAM50_final", value_cols=val_cols).to_csv(
                out_dir / "signatures_and_thornsson_scores_by_pam50_original.csv", index=False
            )

    print(f"[OK] Wrote biological QC to: {out_dir}")


if __name__ == "__main__":
    main()

