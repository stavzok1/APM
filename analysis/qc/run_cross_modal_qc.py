#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.qc._context import load_context, load_sample_ids, to_tcga_sample_id, utc_run_id  # noqa: E402
from analysis.qc._report_utils import ensure_out_dir, md_table_preview, write_csv  # noqa: E402
from pipeline.config import PATHS  # noqa: E402
from pipeline.RNA_exp.signatures import GeneSets, compute_gene_set_scores_from_tpm, read_tpm_wide_subset  # noqa: E402
from pipeline.regulatory_elements import load_regulatory_element_focus  # noqa: E402


def _spearman_r(x: pd.Series, y: pd.Series) -> float:
    xr = x.rank(method="average")
    yr = y.rank(method="average")
    m = xr.notna() & yr.notna()
    if int(m.sum()) < 8:
        return float("nan")
    return float(np.corrcoef(xr[m].astype(float), yr[m].astype(float))[0, 1])


def _load_gene_meth_matrix(path: Path, *, genes: List[str], sample_ids: List[str]) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if "gene_name" in df.columns:
        df = df.set_index("gene_name")
    if "meta" in df.columns:
        df = df.drop(columns=["meta"])
    keep = [s for s in sample_ids if s in df.columns]
    out = df.reindex(index=[g for g in genes if g in df.index], columns=keep)
    for c in keep:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def _load_cnv_gene_calls(path: Path, *, genes: List[str]) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    gcol = "gene_name" if "gene_name" in df.columns else ("gene" if "gene" in df.columns else None)
    if gcol is None:
        return pd.DataFrame()
    df[gcol] = df[gcol].astype(str)
    sub = df[df[gcol].isin(list(genes))].copy()
    if sub.empty:
        return pd.DataFrame()
    if "copy_number" in sub.columns:
        cn = pd.to_numeric(sub["copy_number"], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            log2r = np.log2(cn / 2.0)
        out = sub[[gcol]].copy()
        out["cnv_log2"] = log2r
        return out.drop_duplicates(subset=[gcol]).set_index(gcol)[["cnv_log2"]]
    # fallback numeric
    num_col = None
    for cand in ("segment_mean", "log2_copy_ratio", "log2_ratio", "cnv_log2", "segment_mean_log2"):
        if cand in sub.columns:
            num_col = cand
            break
    if num_col is None:
        return pd.DataFrame()
    sub[num_col] = pd.to_numeric(sub[num_col], errors="coerce")
    out = sub[[gcol, num_col]].drop_duplicates(subset=[gcol]).set_index(gcol)
    out.columns = ["cnv_log2"]
    return out


def _qc_scanning_columns_integrity(out_dir: Path, *, n_rows: int = 200) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """
    Tier-1 structural check: derived `scan_*` columns match their nested sources.

    Reads only a small subset of rows/columns from the element focus parquet and recomputes.
    """
    # Load a small slice with nested decoded (pandas object) plus scan columns.
    import pyarrow.parquet as pq

    want_cols = [
        "cCRE_id",
        "gene_links",
        "screen_exp",
        "ABC_enhancers",
        "hichip",
        "chip_hits",
        "scan_gene_links_n_genes",
        "scan_screen_exp_MCF7_max_score",
        "scan_ABC_max_score_any_celltype",
        "scan_hichip_n_loops_MCF7",
        "scan_chip_hits_has_STAT1",
    ]
    pf = pq.ParquetFile(Path(PATHS.regulatory_elements_table_with_evidence_parquet))
    available = set(pf.schema.names)
    cols = [c for c in want_cols if c in available]
    missing = [c for c in want_cols if c not in available]

    elem = load_regulatory_element_focus(
        path=PATHS.regulatory_elements_table_with_evidence_parquet,
        columns=cols,
        decode_nested=True,
    )
    if len(elem) > n_rows:
        elem = elem.head(n_rows).copy()

    from pipeline.scanning_columns import derive_elem_focus_scanning_columns

    base_cols = [c for c in ("cCRE_id", "gene_links", "screen_exp", "ABC_enhancers", "hichip", "chip_hits") if c in elem.columns]
    recomputed = derive_elem_focus_scanning_columns(elem[base_cols])

    scan_cols = [c for c in cols if c.startswith("scan_")]
    joined = elem[["cCRE_id"] + scan_cols].merge(
        recomputed[["cCRE_id"] + [c for c in recomputed.columns if c.startswith("scan_")]],
        on="cCRE_id",
        how="left",
        suffixes=("", "__recomputed"),
    )

    checks = []
    def _eq(a: pd.Series, b: pd.Series) -> pd.Series:
        # compare with NaNs allowed
        return (a == b) | (a.isna() & b.isna())

    for c in [c for c in cols if c.startswith("scan_")]:
        if c + "__recomputed" not in joined.columns:
            continue
        ok = _eq(joined[c], joined[c + "__recomputed"])
        checks.append({"check": c, "n": int(len(ok)), "n_ok": int(ok.sum()), "ok_rate": float(ok.mean())})

    per_check = pd.DataFrame(checks)
    if not per_check.empty and "ok_rate" in per_check.columns:
        per_check = per_check.sort_values("ok_rate")
    write_csv(per_check, out_dir / "scanning_columns_integrity.csv")
    summary = {
        "n_rows_checked": int(len(elem)),
        "n_scan_cols_checked": int(len(per_check)),
        "missing_requested_columns": missing,
        "min_ok_rate": float(per_check["ok_rate"].min()) if not per_check.empty else 1.0,
    }
    return per_check, summary


def _qc_methylation_vs_rna(
    out_dir: Path,
    *,
    ctx,
    sample_ids: List[str],
    genes: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """
    Promoter methylation (beta) vs RNA expression (log2TPM+1) for a key gene set.
    """
    gene_meth_path = Path(ctx.methylation_working_dir) / "cohort" / "gene_meth_matrix.csv"
    meth = _load_gene_meth_matrix(gene_meth_path, genes=genes, sample_ids=sample_ids) if gene_meth_path.exists() else pd.DataFrame()
    expr = read_tpm_wide_subset(PATHS.rna_expression, genes=genes, sample_cols=sample_ids, gene_col=PATHS.rna_gene_col)

    # per-gene correlation
    per_gene_rows = []
    for g in genes:
        if meth.empty or g not in meth.index or g not in expr.index:
            continue
        x = pd.to_numeric(meth.loc[g], errors="coerce")
        y = pd.to_numeric(expr.loc[g], errors="coerce")
        m = x.notna() & y.notna()
        if int(m.sum()) < 8:
            continue
        per_gene_rows.append(
            {
                "gene": g,
                "n": int(m.sum()),
                "spearman_r_promoter_beta_vs_expr": _spearman_r(x, y),
                "pearson_r_promoter_beta_vs_expr": float(np.corrcoef(x[m].astype(float), y[m].astype(float))[0, 1]),
            }
        )
    per_gene = pd.DataFrame(per_gene_rows).sort_values("spearman_r_promoter_beta_vs_expr")
    write_csv(per_gene, out_dir / "methylation_vs_rna__per_gene.csv")

    # per-sample discordance counts (high beta + high expr) / (low beta + low expr)
    per_sample_rows = []
    if not meth.empty and not expr.empty:
        # percentile thresholds computed per gene across cohort
        expr_q25 = expr.apply(lambda r: pd.to_numeric(r, errors="coerce").quantile(0.25), axis=1)
        expr_q75 = expr.apply(lambda r: pd.to_numeric(r, errors="coerce").quantile(0.75), axis=1)
        for sid in sample_ids:
            rec = {"sample_id": sid, "tcga_sample_id": to_tcga_sample_id(sid)}
            n_hi_hi = 0
            n_lo_lo = 0
            n_total = 0
            for g in genes:
                if g not in meth.index or g not in expr.index or sid not in meth.columns or sid not in expr.columns:
                    continue
                b = _safe_float(meth.at[g, sid])
                e = _safe_float(expr.at[g, sid])
                if b is None or e is None:
                    continue
                n_total += 1
                if b > 0.7 and e > float(expr_q75.get(g, np.nan)):
                    n_hi_hi += 1
                if b < 0.2 and e < float(expr_q25.get(g, np.nan)):
                    n_lo_lo += 1
            rec["meth_expr_n_total"] = n_total
            rec["meth_expr_n_beta_hi_expr_hi"] = n_hi_hi
            rec["meth_expr_n_beta_lo_expr_lo"] = n_lo_lo
            per_sample_rows.append(rec)
    per_sample = pd.DataFrame(per_sample_rows)
    write_csv(per_sample, out_dir / "methylation_vs_rna__per_sample.csv")

    summary = {
        "n_genes_with_corr": int(len(per_gene)),
        "median_spearman_r": float(per_gene["spearman_r_promoter_beta_vs_expr"].median()) if not per_gene.empty else np.nan,
    }
    return per_sample, per_gene, summary


def _safe_float(x: object) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, float) and x != x:
            return None
        return float(x)
    except Exception:
        return None


def _qc_cnv_vs_rna(
    out_dir: Path,
    *,
    ctx,
    sample_ids: List[str],
    genes: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    expr = read_tpm_wide_subset(PATHS.rna_expression, genes=genes, sample_cols=sample_ids, gene_col=PATHS.rna_gene_col)

    cnv_mat: Dict[str, pd.Series] = {}
    for sid in sample_ids:
        p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls_ascat3.csv"
        if not p.exists():
            p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls.csv"
        if not p.exists():
            continue
        gdf = _load_cnv_gene_calls(p, genes=genes)
        if not gdf.empty and "cnv_log2" in gdf.columns:
            cnv_mat[sid] = gdf["cnv_log2"]
    cnv = pd.DataFrame(cnv_mat).reindex(index=genes)
    cnv.to_csv(out_dir / "cnv_log2_matrix_genes.csv")

    per_gene_rows = []
    for g in genes:
        if g not in cnv.index or g not in expr.index:
            continue
        x = pd.to_numeric(cnv.loc[g], errors="coerce")
        y = pd.to_numeric(expr.loc[g], errors="coerce")
        m = x.notna() & y.notna()
        if int(m.sum()) < 8:
            continue
        per_gene_rows.append(
            {
                "gene": g,
                "n": int(m.sum()),
                "pearson_r_cnv_vs_expr": float(np.corrcoef(x[m].astype(float), y[m].astype(float))[0, 1]),
                "spearman_r_cnv_vs_expr": _spearman_r(x, y),
            }
        )
    per_gene = pd.DataFrame(per_gene_rows).sort_values("pearson_r_cnv_vs_expr")
    write_csv(per_gene, out_dir / "cnv_vs_rna__per_gene.csv")

    # per-sample: fraction of genes with concordant sign (cnv_log2 and expr deviation from median)
    per_sample_rows = []
    if not cnv.empty and not expr.empty:
        expr_med = expr.apply(lambda r: pd.to_numeric(r, errors="coerce").median(), axis=1)
        for sid in sample_ids:
            rec = {"sample_id": sid, "tcga_sample_id": to_tcga_sample_id(sid)}
            ok = 0
            tot = 0
            for g in genes:
                if g not in cnv.index or g not in expr.index or sid not in cnv.columns or sid not in expr.columns:
                    continue
                c = _safe_float(cnv.at[g, sid])
                e = _safe_float(expr.at[g, sid])
                if c is None or e is None:
                    continue
                med = _safe_float(expr_med.get(g))
                if med is None:
                    continue
                tot += 1
                # dosage concordance heuristic: cnv gain -> expr above median; loss -> below
                if (c > 0 and e >= med) or (c < 0 and e <= med) or (abs(c) < 1e-6):
                    ok += 1
            rec["cnv_expr_n_total"] = tot
            rec["cnv_expr_concordant_frac"] = float(ok / tot) if tot else np.nan
            per_sample_rows.append(rec)
    per_sample = pd.DataFrame(per_sample_rows)
    write_csv(per_sample, out_dir / "cnv_vs_rna__per_sample.csv")

    summary = {
        "n_genes_with_corr": int(len(per_gene)),
        "median_pearson_r": float(per_gene["pearson_r_cnv_vs_expr"].median()) if not per_gene.empty else np.nan,
    }
    return per_sample, per_gene, summary


def _load_thornsson() -> pd.DataFrame:
    p = Path("annotations") / "Thornsson_immune_table.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", low_memory=False)
    if "TCGA Participant Barcode" in df.columns:
        df["tcga_sample_id"] = df["TCGA Participant Barcode"].map(to_tcga_sample_id)
    return df


def _load_atac_case_level_subset(
    *,
    peak_ids: List[str],
    sample_cols: List[str],
) -> pd.DataFrame:
    """
    Load a small subset of the wide ATAC case-level matrix for just *peak_ids* × *sample_cols*.

    The matrix is large, so this streams in chunks and filters rows by peak id.
    """
    path = Path(PATHS.atac_case_level_matrix)
    if not path.exists():
        return pd.DataFrame()

    hdr = pd.read_csv(path, nrows=0, low_memory=False)
    cols = list(hdr.columns)

    if "peak_id" in cols:
        id_cols = ["peak_id"]
        construct_peak_id = False
    else:
        cand = [c for c in ("chrom", "chr", "Chromosome", "seqnames") if c in cols]
        if not cand:
            return pd.DataFrame()
        chrom_col = cand[0]
        start_col = "start" if "start" in cols else ("Start" if "Start" in cols else None)
        end_col = "end" if "end" in cols else ("End" if "End" in cols else None)
        if start_col is None or end_col is None:
            return pd.DataFrame()
        id_cols = [chrom_col, start_col, end_col]
        construct_peak_id = True

    # Sample columns can be full aliquot barcodes (e.g. TCGA-XX-YYYY-01A-..),
    # while cohort ids are usually vials (TCGA-XX-YYYY-01A) or samples (TCGA-XX-YYYY-01).
    # Accept exact matches and prefix matches.
    prefixes = [str(x) for x in sample_cols if str(x)]
    want_samples = []
    for c in cols:
        if c in id_cols:
            continue
        if c in prefixes:
            want_samples.append(c)
            continue
        for p in prefixes:
            if c.startswith(p):
                want_samples.append(c)
                break
    if not want_samples:
        return pd.DataFrame()

    usecols = id_cols + want_samples
    wanted = set(peak_ids)
    out_rows: List[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=50_000, low_memory=False):
        if construct_peak_id:
            pid = (
                chunk[id_cols[0]].astype(str)
                + ":"
                + pd.to_numeric(chunk[id_cols[1]], errors="coerce").astype("Int64").astype(str)
                + "-"
                + pd.to_numeric(chunk[id_cols[2]], errors="coerce").astype("Int64").astype(str)
            )
            chunk = chunk.copy()
            chunk["peak_id"] = pid
        else:
            pid = chunk["peak_id"].astype(str)

        m = pid.isin(wanted)
        if bool(m.any()):
            sub = chunk.loc[m, ["peak_id"] + want_samples].copy()
            out_rows.append(sub)

    if not out_rows:
        return pd.DataFrame()

    df = pd.concat(out_rows, ignore_index=True)
    df = df.drop_duplicates(subset=["peak_id"], keep="first").set_index("peak_id")
    # Normalize sample column names to TCGA sample ids when possible (helps alignment).
    df = df.rename(columns={c: (to_tcga_sample_id(c) or c) for c in df.columns})
    # If multiple ATAC aliquots map to the same TCGA sample id, average them.
    if len(set(df.columns)) < len(df.columns):
        df = df.T.groupby(level=0).mean(numeric_only=True).T
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _qc_atac_vs_rna_cis(
    out_dir: Path,
    *,
    sample_ids: List[str],
    genes: List[str],
    max_pairs: int = 600,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """
    ATAC accessibility ↔ RNA expression in cis, for peaks linked via element-focus `atac_peak_links`.
    """
    # ATAC case-level matrices can use varying TCGA barcode granularities. Provide a broad
    # candidate column list and later normalize both ATAC+RNA to TCGA sample ids.
    tcga_samples = []
    for s in sample_ids:
        tcga_samples.append(str(s))
        tid = to_tcga_sample_id(s)
        if tid:
            tcga_samples.append(tid)
    tcga_samples = [s for s in tcga_samples if s]

    elem = load_regulatory_element_focus(
        path=PATHS.regulatory_elements_table_with_evidence_parquet,
        columns=["cCRE_id", "gene_links", "atac_peak_links"],
        decode_nested=True,
    )
    if elem.empty:
        return pd.DataFrame(), pd.DataFrame(), {"n_pairs": 0, "atac_loaded": False}

    def _ccre_gene_evidence_flags(gd: object) -> Dict[str, object]:
        """
        Extract lightweight evidence flags from an element-focus `gene_links[gene]` bundle.

        We *do not* hard-filter pairs by evidence (user prefers including distance-tier links),
        but we keep tags so downstream summaries can stratify by "supported vs distance-only".
        """
        if not isinstance(gd, dict):
            return {
                "has_screen_strong": False,
                "has_hichip": False,
                "has_abc": False,
                "dist_to_tss": None,
                "tier": None,
            }

        has_screen_strong = False
        se = gd.get("screen_exp")
        if isinstance(se, dict):
            pb = se.get("per_biosample")
            if isinstance(pb, dict):
                for _, assays in pb.items():
                    if not isinstance(assays, dict):
                        continue
                    for _, rec in assays.items():
                        if isinstance(rec, dict) and str(rec.get("strength", "")).lower() == "strong":
                            has_screen_strong = True
                            break
                    if has_screen_strong:
                        break

        hh = gd.get("hichip")
        has_hichip = isinstance(hh, dict) and len(hh) > 0

        # Keep ABC as an informational tag only (do not require it).
        abc = gd.get("ABC_enhancers")
        has_abc = isinstance(abc, list) and len(abc) > 0

        dist = gd.get("dist_to_tss")
        tier = gd.get("tier")
        try:
            dist_to_tss = int(dist) if dist is not None and dist == dist else None
        except Exception:
            dist_to_tss = None
        tier_s = str(tier) if tier is not None else None

        return {
            "has_screen_strong": bool(has_screen_strong),
            "has_hichip": bool(has_hichip),
            "has_abc": bool(has_abc),
            "dist_to_tss": dist_to_tss,
            "tier": tier_s,
        }

    pairs = []
    for _, r in elem.iterrows():
        gl = r.get("gene_links")
        apl = r.get("atac_peak_links")
        if not isinstance(gl, dict) or not isinstance(apl, list) or not apl:
            continue
        for g, gd in gl.items():
            if g not in genes:
                continue
            ev = _ccre_gene_evidence_flags(gd)
            ev_support = bool(ev.get("has_screen_strong")) or bool(ev.get("has_hichip"))
            for p in apl:
                if not isinstance(p, dict):
                    continue
                pid = p.get("peak_id")
                if not pid:
                    continue
                pairs.append(
                    {
                        "gene": str(g),
                        "peak_id": str(pid),
                        "abs_dist": p.get("abs_dist"),
                        "overlaps": p.get("overlaps"),
                        "ccre_gene_tier": ev.get("tier"),
                        "ccre_gene_dist_to_tss": ev.get("dist_to_tss"),
                        "ccre_has_screen_strong": ev.get("has_screen_strong"),
                        "ccre_has_hichip": ev.get("has_hichip"),
                        "ccre_has_abc": ev.get("has_abc"),
                        "ccre_evidence_supported": ev_support,
                    }
                )
    if not pairs:
        return pd.DataFrame(), pd.DataFrame(), {"n_pairs": 0, "atac_loaded": False}

    pairs_df = pd.DataFrame(pairs).drop_duplicates(subset=["gene", "peak_id"])
    if len(pairs_df) > max_pairs:
        pairs_df = pairs_df.head(max_pairs).copy()

    peak_ids = pairs_df["peak_id"].astype(str).unique().tolist()
    atac = _load_atac_case_level_subset(peak_ids=peak_ids, sample_cols=tcga_samples)
    if atac.empty:
        return pd.DataFrame(), pd.DataFrame(), {"n_pairs": int(len(pairs_df)), "atac_loaded": False}

    expr = read_tpm_wide_subset(PATHS.rna_expression, genes=genes, sample_cols=sample_ids, gene_col=PATHS.rna_gene_col)
    expr = expr.rename(columns={c: (to_tcga_sample_id(c) or c) for c in expr.columns})

    per_pair_rows = []
    for _, pr in pairs_df.iterrows():
        g = str(pr["gene"])
        pid = str(pr["peak_id"])
        if pid not in atac.index or g not in expr.index:
            continue
        x = pd.to_numeric(atac.loc[pid], errors="coerce")
        y = pd.to_numeric(expr.loc[g], errors="coerce")
        m = x.notna() & y.notna()
        if int(m.sum()) < 8:
            continue
        per_pair_rows.append(
            {
                "gene": g,
                "peak_id": pid,
                "n": int(m.sum()),
                "spearman_r_atac_vs_expr": _spearman_r(x, y),
                "pearson_r_atac_vs_expr": float(np.corrcoef(x[m].astype(float), y[m].astype(float))[0, 1]),
                "abs_dist": pr.get("abs_dist"),
                "overlaps": pr.get("overlaps"),
                "ccre_gene_tier": pr.get("ccre_gene_tier"),
                "ccre_gene_dist_to_tss": pr.get("ccre_gene_dist_to_tss"),
                "ccre_evidence_supported": pr.get("ccre_evidence_supported"),
            }
        )
    per_peak = pd.DataFrame(per_pair_rows)
    if not per_peak.empty:
        per_peak = per_peak.sort_values("spearman_r_atac_vs_expr")
    write_csv(per_peak, out_dir / "atac_vs_rna_cis__per_peak.csv")

    per_gene = pd.DataFrame()
    if not per_peak.empty:
        per_gene = (
            per_peak.groupby("gene", as_index=False)
            .agg(
                n_peaks=("peak_id", "nunique"),
                median_spearman=("spearman_r_atac_vs_expr", "median"),
                median_pearson=("pearson_r_atac_vs_expr", "median"),
            )
            .sort_values("median_spearman")
        )
    write_csv(per_gene, out_dir / "atac_vs_rna_cis__per_gene.csv")

    # Gene-level aggregation: per sample, take max ATAC across linked peaks for the gene.
    per_gene_agg = pd.DataFrame()
    if not per_peak.empty:
        gene_rows = []
        for g, sub in pairs_df.groupby("gene", sort=False):
            pids = [p for p in sub["peak_id"].astype(str).unique().tolist() if p in atac.index]
            if not pids or g not in expr.index:
                continue
            # Evidence support rate for this gene's linked cCREs (informational)
            ev_rate = float(pd.to_numeric(sub.get("ccre_evidence_supported"), errors="coerce").fillna(False).mean())
            atac_gene = atac.loc[pids]
            # collapse peak x sample -> sample
            atac_score = atac_gene.max(axis=0)
            rna = pd.to_numeric(expr.loc[g], errors="coerce")
            m = atac_score.notna() & rna.notna()
            if int(m.sum()) < 8:
                continue
            gene_rows.append(
                {
                    "gene": str(g),
                    "n_peaks": int(len(pids)),
                    "n": int(m.sum()),
                    "spearman_r_gene_atacmax_vs_expr": _spearman_r(atac_score, rna),
                    "pearson_r_gene_atacmax_vs_expr": float(np.corrcoef(atac_score[m].astype(float), rna[m].astype(float))[0, 1]),
                    "ccre_evidence_supported_rate": ev_rate,
                }
            )
        per_gene_agg = pd.DataFrame(gene_rows).sort_values("spearman_r_gene_atacmax_vs_expr")
    write_csv(per_gene_agg, out_dir / "atac_vs_rna_cis__per_gene_agg.csv")

    summary = {
        "n_pairs": int(len(pairs_df)),
        "n_peaks_loaded": int(len(atac)),
        "n_peaks_with_corr": int(per_peak["peak_id"].nunique()) if not per_peak.empty else 0,
        "median_spearman_r": float(per_peak["spearman_r_atac_vs_expr"].median()) if not per_peak.empty else np.nan,
        "n_genes_with_gene_agg_corr": int(len(per_gene_agg)) if not per_gene_agg.empty else 0,
        "median_spearman_r_gene_agg": float(per_gene_agg["spearman_r_gene_atacmax_vs_expr"].median()) if not per_gene_agg.empty else np.nan,
        "atac_loaded": True,
    }
    return per_peak, per_gene, summary


def _qc_ifn_axis(out_dir: Path, *, sample_ids: List[str]) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """
    IFN-gamma Response (Thorsson) vs APM response (RNA signatures).
    """
    sig = compute_gene_set_scores_from_tpm(PATHS.rna_expression, sample_ids=sample_ids, gene_col=PATHS.rna_gene_col).reset_index()
    sig["tcga_sample_id"] = sig["sample_id"].map(to_tcga_sample_id)
    th = _load_thornsson()
    if th.empty or "IFN-gamma Response" not in th.columns:
        out = sig[["sample_id", "tcga_sample_id"]].copy()
        out["ifn_gamma_response"] = np.nan
        out["apm_response_to_ifn_score"] = np.nan
        write_csv(out, out_dir / "ifn_axis__per_sample.csv")
        return out, {"n_joined": 0}

    th_small = th[["tcga_sample_id", "IFN-gamma Response"]].drop_duplicates(subset=["tcga_sample_id"])
    joined = sig.merge(th_small, on="tcga_sample_id", how="left")
    joined = joined.rename(columns={"IFN-gamma Response": "ifn_gamma_response"})

    # Define APM response as APM_classI_mean (already computed in GeneSets), normalized by cohort median.
    if "APM_classI_mean" in joined.columns:
        med = float(pd.to_numeric(joined["APM_classI_mean"], errors="coerce").median())
        joined["apm_response_to_ifn_score"] = pd.to_numeric(joined["APM_classI_mean"], errors="coerce") - med
    else:
        joined["apm_response_to_ifn_score"] = np.nan

    # Flag IFN-high but APM-low
    q75 = float(pd.to_numeric(joined["ifn_gamma_response"], errors="coerce").quantile(0.75))
    q25_apm = float(pd.to_numeric(joined["apm_response_to_ifn_score"], errors="coerce").quantile(0.25))
    joined["ifn_high_apm_low_flag"] = (
        (pd.to_numeric(joined["ifn_gamma_response"], errors="coerce") >= q75)
        & (pd.to_numeric(joined["apm_response_to_ifn_score"], errors="coerce") <= q25_apm)
    )
    write_csv(joined, out_dir / "ifn_axis__per_sample.csv")

    summary = {
        "n_joined": int(joined["ifn_gamma_response"].notna().sum()),
        "n_ifn_high_apm_low": int(joined["ifn_high_apm_low_flag"].sum()),
    }
    return joined, summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch-json", type=str, required=True)
    ap.add_argument("--n-rows-scan-check", type=int, default=200)
    args = ap.parse_args()

    ctx = load_context(Path(args.scratch_json))
    run_id = utc_run_id("cross_modal_qc")
    out_dir = ensure_out_dir(run_id)

    sample_ids = load_sample_ids(ctx)
    if not sample_ids:
        raise SystemExit("No sample_ids found for cohort.")

    gs = GeneSets()
    key_genes = list(gs.apm_class_i)

    # --- Tier 1: scanning columns integrity (hard-fail if ok_rate too low) ---
    scan_checks, scan_summary = _qc_scanning_columns_integrity(out_dir, n_rows=int(args.n_rows_scan_check))
    # conservative hard-fail threshold (structural bug)
    if int(scan_summary.get("n_scan_cols_checked", 0)) > 0 and scan_summary.get("min_ok_rate", 1.0) < 0.999:
        raise SystemExit(
            "Scanning-column integrity check failed (min_ok_rate < 0.999). "
            f"See: {out_dir / 'scanning_columns_integrity.csv'}"
        )

    # --- Cross-modal checks ---
    meth_samp, meth_gene, meth_sum = _qc_methylation_vs_rna(out_dir, ctx=ctx, sample_ids=sample_ids, genes=key_genes)
    cnv_samp, cnv_gene, cnv_sum = _qc_cnv_vs_rna(out_dir, ctx=ctx, sample_ids=sample_ids, genes=key_genes)
    ifn_samp, ifn_sum = _qc_ifn_axis(out_dir, sample_ids=sample_ids)
    atac_peak, atac_gene, atac_sum = _qc_atac_vs_rna_cis(out_dir, sample_ids=sample_ids, genes=key_genes)

    # --- Assemble per-sample QC card (pass/warn/fail left to downstream filtering for now) ---
    card = pd.DataFrame({"sample_id": sample_ids})
    card["tcga_sample_id"] = card["sample_id"].map(to_tcga_sample_id)
    if not meth_samp.empty:
        card = card.merge(meth_samp[["sample_id", "meth_expr_n_total", "meth_expr_n_beta_hi_expr_hi", "meth_expr_n_beta_lo_expr_lo"]], on="sample_id", how="left")
    if not cnv_samp.empty:
        card = card.merge(cnv_samp[["sample_id", "cnv_expr_n_total", "cnv_expr_concordant_frac"]], on="sample_id", how="left")
    if not ifn_samp.empty and "ifn_high_apm_low_flag" in ifn_samp.columns:
        card = card.merge(ifn_samp[["sample_id", "ifn_gamma_response", "apm_response_to_ifn_score", "ifn_high_apm_low_flag"]], on="sample_id", how="left")
    write_csv(card, out_dir / "per_sample_qc_card.csv")

    # per-gene card (limited set for now)
    gene_card = meth_gene.merge(cnv_gene, on="gene", how="outer", suffixes=("", "__cnv"))
    write_csv(gene_card, out_dir / "per_gene_qc_card.csv")

    # summary json + report
    summary = {
        "run_id": run_id,
        "n_samples": int(len(sample_ids)),
        "scratch_json": str(Path(args.scratch_json)),
        "tier1_scanning_columns_integrity": scan_summary,
        "methylation_vs_rna": meth_sum,
        "cnv_vs_rna": cnv_sum,
        "ifn_axis": ifn_sum,
        "atac_vs_rna_cis": atac_sum,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines: List[str] = []
    lines.append("## Cross-modal QC report\n\n")
    lines.append(f"- **run_id**: `{run_id}`\n")
    lines.append(f"- **n_samples**: {len(sample_ids)}\n")
    lines.append(f"- **scratch_json**: `{Path(args.scratch_json)}`\n\n")

    lines.append("### Tier 1 — scanning columns integrity (structural)\n\n")
    lines.append(f"- min_ok_rate: `{scan_summary.get('min_ok_rate')}`\n\n")
    if not scan_checks.empty:
        lines.append(md_table_preview(scan_checks.head(10)))
        lines.append("\n")

    lines.append("### Promoter methylation ↔ RNA (key APM genes)\n\n")
    if not meth_gene.empty:
        lines.append(md_table_preview(meth_gene.head(12)))
        lines.append("\n")

    lines.append("### CNV ↔ RNA dosage (key APM genes)\n\n")
    if not cnv_gene.empty:
        lines.append(md_table_preview(cnv_gene.head(12)))
        lines.append("\n")

    lines.append("### IFN axis ↔ APM response\n\n")
    if not ifn_samp.empty:
        cols = [c for c in ("sample_id", "ifn_gamma_response", "apm_response_to_ifn_score", "ifn_high_apm_low_flag") if c in ifn_samp.columns]
        lines.append(md_table_preview(ifn_samp[cols].head(12)))
        lines.append("\n")

    lines.append("### ATAC ↔ RNA (cis; linked peaks via element-focus `atac_peak_links`)\n\n")
    if not atac_peak.empty:
        lines.append(md_table_preview(atac_peak.head(12)))
        lines.append("\n")
    p_gene_agg = out_dir / "atac_vs_rna_cis__per_gene_agg.csv"
    if p_gene_agg.exists():
        try:
            gd = pd.read_csv(p_gene_agg)
            if not gd.empty:
                lines.append("Gene-level aggregation: per-sample `max(ATAC across linked peaks)` vs RNA.\n\n")
                lines.append(md_table_preview(gd.head(12)))
                lines.append("\n")
        except Exception:
            pass

    (out_dir / "report.md").write_text("".join(lines), encoding="utf-8")
    print(f"[OK] wrote: {out_dir / 'report.md'}")


if __name__ == "__main__":
    main()

