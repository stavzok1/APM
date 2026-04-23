from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline.config import PATHS
from pipeline.qc.base import Finding, Level
from pipeline.qc.context import to_tcga_sample_id
from pipeline.qc.hla_snv_vep import build_tumor_vial_to_snv_per_sample_csv, load_hla_snv_summary_for_vial
from pipeline.qc.hla_loh_ascat_segments import (
    build_tumor_vial_to_ascat_seg_path,
    hla_loh_from_segments_for_genes,
    load_ascat_allelic_seg,
    load_hla_gene_intervals_from_gene_table,
    load_hla_types_by_aliquot,
    match_hla_types_row,
)

HLA_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M"]


def _load_thornsson() -> pd.DataFrame:
    p = Path("annotations") / "Thornsson_immune_table.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", low_memory=False)
    if "TCGA Participant Barcode" in df.columns:
        df["tcga_sample_id"] = df["TCGA Participant Barcode"].map(to_tcga_sample_id)
    return df


def compute_hla_loh_vs_neoantigen_and_infiltration(
    *,
    ctx,
    sample_ids: List[str],
    out_dir: Path | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """
    Tier-3 check:
    - HLA LOH from ASCAT3 **allelic segments** (preferred when available), else gene-table ``loh_flag``.
    - Gene-table flags are kept because some installs may lack segment files.
    - Optional HLA typing columns from ``PATHS.hla_types_tsv`` (matched by aliquot prefix).
    - Compare neoantigen counts and infiltration proxies from Thorsson table.

    This is a *sanity* check; it won't prove causality.
    """
    per_gene = pd.DataFrame()
    th = _load_thornsson()
    summary: Dict[str, object] = {"thornsson_exists": not th.empty}

    vial_to_seg = build_tumor_vial_to_ascat_seg_path(
        cnv_seg_dir=PATHS.cnv_dir,
        ann_path=PATHS.cnv_annotations_path,
    )
    summary["n_ascat_seg_files_indexed"] = int(len(vial_to_seg))

    vial_to_snv = build_tumor_vial_to_snv_per_sample_csv(
        snv_output_dir=Path(ctx.snv_output),
        snv_manifest_path=PATHS.annotations_dir / "SNV" / "samples.tsv",
    )
    summary["n_snv_per_sample_files_indexed"] = int(len(vial_to_snv))

    hla_types_df = load_hla_types_by_aliquot(PATHS.hla_types_tsv)

    rows: List[Dict[str, object]] = []
    n_seg_ok = 0
    n_seg_mismatch_gene = 0
    for sid in sample_ids:
        p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls_ascat3.csv"
        if not p.exists():
            p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls.csv"
        tcga = to_tcga_sample_id(sid)
        base: Dict[str, object] = {"sample_id": sid, "tcga_sample_id": tcga}

        if not p.exists():
            base["hla_loh_gene_table_any"] = np.nan
            base["hla_loh_seg_any"] = np.nan
            base["hla_loh_any"] = np.nan
            rows.append(base)
            continue

        df = pd.read_csv(p, usecols=["gene_name", "loh_flag"], low_memory=False)
        df["gene_name"] = df["gene_name"].astype(str)
        sub = df[df["gene_name"].isin(HLA_GENES)].copy()
        gene_loh = (
            bool(pd.to_numeric(sub["loh_flag"], errors="coerce").fillna(0).astype(int).max() > 0)
            if not sub.empty
            else False
        )
        base["hla_loh_gene_table_any"] = bool(gene_loh)

        seg_loh: Optional[bool] = None
        seg_path = vial_to_seg.get(str(sid).strip())
        if seg_path is not None and seg_path.is_file():
            try:
                intervals = load_hla_gene_intervals_from_gene_table(p)
                if intervals:
                    segs = load_ascat_allelic_seg(seg_path)
                    seg_rec = hla_loh_from_segments_for_genes(segs, intervals)
                    for k, v in seg_rec.items():
                        base[k] = v
                    seg_loh = bool(seg_rec.get("hla_loh_seg_any", False))
                    n_seg_ok += 1
                    if seg_loh != gene_loh:
                        n_seg_mismatch_gene += 1
            except Exception as exc:
                base["hla_loh_seg_error"] = str(exc)[:500]
                seg_loh = None

        if seg_loh is None:
            base["hla_loh_seg_any"] = np.nan
            base["hla_loh_any"] = bool(gene_loh)
        else:
            base["hla_loh_seg_any"] = bool(seg_loh)
            # Prefer segment-based LOH when available (captures minor==0 including copy-neutral LOH).
            base["hla_loh_any"] = bool(seg_loh)

        hrow = match_hla_types_row(hla_types_df, str(sid).strip())
        if hrow is not None:
            base["hla_typing_aliquot_id"] = str(hrow.get("aliquot_id", ""))
            for col in ("A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"):
                if col in hrow.index:
                    base[f"hla_type_{col}"] = hrow.get(col)

        snv_h = load_hla_snv_summary_for_vial(tumor_vial=str(sid).strip(), vial_to_snv_csv=vial_to_snv)
        for k, v in snv_h.items():
            base[k] = v

        # Combined presentation-restriction proxy (CN/segment LOH **or** SNV damage on HLA loci with
        # population + pathogenicity gates; see ``pipeline/qc/hla_snv_vep.py`` constants).
        def _b(x: object) -> bool:
            if isinstance(x, bool):
                return x
            try:
                if x is None or (isinstance(x, float) and np.isnan(x)):
                    return False
            except Exception:
                return False
            return bool(x)

        base["hla_damage_cnv_or_snv_any"] = bool(
            _b(base.get("hla_loh_any"))
            or _b(snv_h.get("hla_snv_lof_rare_any"))
            or _b(snv_h.get("hla_snv_high_missense_damage_any"))
            or _b(snv_h.get("hla_snv_clin_path_damage_any"))
        )

        rows.append(base)

    per_sample = pd.DataFrame(rows)

    summary["n_seg_samples_with_hla_intervals"] = int(n_seg_ok)
    summary["n_seg_vs_gene_mismatch"] = int(n_seg_mismatch_gene)

    if "hla_loh_any" in per_sample.columns and per_sample["hla_loh_any"].notna().any():
        m = per_sample["hla_loh_any"].notna()
        summary["n_with_loh_call"] = int(m.sum())
        summary["loh_rate"] = float(per_sample.loc[m, "hla_loh_any"].mean())

    if "hla_loh_seg_any" in per_sample.columns and per_sample["hla_loh_seg_any"].notna().any():
        m = per_sample["hla_loh_seg_any"].notna()
        summary["loh_rate_seg"] = float(per_sample.loc[m, "hla_loh_seg_any"].astype(bool).mean())

    if "hla_loh_gene_table_any" in per_sample.columns and per_sample["hla_loh_gene_table_any"].notna().any():
        m = per_sample["hla_loh_gene_table_any"].notna()
        summary["loh_rate_gene_table"] = float(per_sample.loc[m, "hla_loh_gene_table_any"].astype(bool).mean())

    # Join Thorsson neoantigen / infiltration columns if present
    if not th.empty:
        keep = [
            c
            for c in [
                "SNV Neoantigens",
                "Indel Neoantigens",
                "Leukocyte Fraction",
                "Lymphocyte Infiltration Signature Score",
                "T Cells CD8",
                "IFN-gamma Response",
                "Immune Subtype",
            ]
            if c in th.columns
        ]
        th_small = th[["tcga_sample_id", *keep]].drop_duplicates(subset=["tcga_sample_id"])
        per_sample = per_sample.merge(th_small, on="tcga_sample_id", how="left")

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        per_sample.to_csv(out_dir / "hla_loh_vs_neoantigen__per_sample.csv", index=False)

    return per_sample, per_gene, summary


def assert_hla_loh_vs_neoantigen(summary: Dict[str, object]) -> List[Finding]:
    findings: List[Finding] = []
    if not summary.get("thornsson_exists", False):
        findings.append(
            Finding(
                check="hla_loh_vs_neoantigen.thornsson_present",
                level=Level.WARN,
                message="Thorsson immune table missing; LOH↔neoantigen QC can’t run.",
                context={},
            )
        )
    if int(summary.get("n_ascat_seg_files_indexed", 0) or 0) == 0:
        findings.append(
            Finding(
                check="hla_loh.ascat_seg_index",
                level=Level.WARN,
                message="No ASCAT3 allelic segment files indexed from CNV manifest; HLA LOH uses gene-table fallback only.",
                context={"cnv_dir": str(PATHS.cnv_dir)},
            )
        )
    if int(summary.get("n_snv_per_sample_files_indexed", 0) or 0) == 0:
        findings.append(
            Finding(
                check="hla_snv.per_sample_index",
                level=Level.WARN,
                message="No per-sample SNV tables indexed from SNV manifest; HLA SNV/VEP damage flags are skipped.",
                context={"snv_manifest": str(PATHS.annotations_dir / "SNV" / "samples.tsv")},
            )
        )
    return findings
