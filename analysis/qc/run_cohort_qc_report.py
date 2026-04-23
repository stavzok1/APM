#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.RNA_exp.signatures import (  # noqa: E402
    GeneSets,
    compute_gene_set_scores_from_tpm,
    read_tpm_wide_subset,
)

from analysis.qc._context import load_context, load_sample_ids, to_tcga_sample_id, utc_run_id  # noqa: E402
from analysis.qc._report_utils import ensure_out_dir, group_summary, md_table_preview, write_csv  # noqa: E402


def _load_brca_immune_advanced() -> pd.DataFrame:
    p = Path("annotations") / "BRCA_immune_subtypes_advanced.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", low_memory=False)
    if "sample_id" in df.columns:
        df["tcga_sample_id"] = df["sample_id"].map(to_tcga_sample_id)
    return df


def _load_thornsson_original() -> pd.DataFrame:
    p = Path("annotations") / "Thornsson_immune_table.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", low_memory=False)
    if "TCGA Participant Barcode" in df.columns:
        df["tcga_sample_id"] = df["TCGA Participant Barcode"].map(to_tcga_sample_id)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch-json", type=str, required=True)
    args = ap.parse_args()

    ctx = load_context(Path(args.scratch_json))
    run_id = utc_run_id("cohort_qc_report")
    out_dir = ensure_out_dir(run_id)

    sample_ids = load_sample_ids(ctx)
    if not sample_ids:
        raise SystemExit("No sample_ids found for cohort.")
    tcga_samples = [to_tcga_sample_id(x) for x in sample_ids]

    # --- RNA signatures ---
    sig = compute_gene_set_scores_from_tpm(
        PATHS.rna_expression,
        sample_ids=sample_ids,
        gene_col=PATHS.rna_gene_col,
    ).reset_index()
    sig["tcga_sample_id"] = sig["sample_id"].map(to_tcga_sample_id)
    write_csv(sig, out_dir / "rna_signatures.csv")

    # --- Immune annotation tables ---
    adv = _load_brca_immune_advanced()
    th = _load_thornsson_original()

    # keep only cohort samples in joined views
    sig_cohort = sig[sig["tcga_sample_id"].isin(set(tcga_samples))].copy()

    joined_adv = pd.DataFrame()
    if not adv.empty and "tcga_sample_id" in adv.columns:
        adv_small_cols = [c for c in ("PAM50_final", "thornsson_subtype", "CPE") if c in adv.columns]
        adv_small = adv[["tcga_sample_id", *adv_small_cols]].drop_duplicates(subset=["tcga_sample_id"])
        joined_adv = sig_cohort.merge(adv_small, on="tcga_sample_id", how="left")
        write_csv(joined_adv, out_dir / "signatures_joined_brca_immune_advanced.csv")

    joined_th = pd.DataFrame()
    if not th.empty and "tcga_sample_id" in th.columns:
        th_keep = [
            c
            for c in (
                "Immune Subtype",
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
        joined_th = sig_cohort.merge(th_small, on="tcga_sample_id", how="left")
        write_csv(joined_th, out_dir / "signatures_joined_thornsson_original.csv")

    sig_cols = [c for c in ("APM_classI_mean", "CD8_mean", "NK_mean", "IFNG_mean", "CYT") if c in sig.columns]

    # summaries
    by_files: List[Path] = []
    if not joined_adv.empty:
        if "PAM50_final" in joined_adv.columns:
            df = group_summary(joined_adv, group_col="PAM50_final", value_cols=sig_cols)
            p = out_dir / "signatures_by_pam50_pipeline_table.csv"
            write_csv(df, p)
            by_files.append(p)
        if "thornsson_subtype" in joined_adv.columns:
            df = group_summary(joined_adv, group_col="thornsson_subtype", value_cols=sig_cols)
            p = out_dir / "signatures_by_thornsson_subtype_pipeline_table.csv"
            write_csv(df, p)
            by_files.append(p)

    if not joined_th.empty:
        extra_cols = [
            c
            for c in (
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
            if c in joined_th.columns
        ]
        val_cols = sig_cols + extra_cols
        if "Immune Subtype" in joined_th.columns:
            df = group_summary(joined_th, group_col="Immune Subtype", value_cols=val_cols)
            p = out_dir / "signatures_and_thornsson_scores_by_immune_subtype_original.csv"
            write_csv(df, p)
            by_files.append(p)
        if "PAM50_final" in joined_th.columns:
            df = group_summary(joined_th, group_col="PAM50_final", value_cols=val_cols)
            p = out_dir / "signatures_and_thornsson_scores_by_pam50_original.csv"
            write_csv(df, p)
            by_files.append(p)

    # markdown report
    lines: List[str] = []
    lines.append("## Cohort QC report\n\n")
    lines.append(f"- **run_id**: `{run_id}`\n")
    lines.append(f"- **n_samples**: {len(sample_ids)}\n")
    lines.append(f"- **scratch_json**: `{Path(args.scratch_json)}`\n\n")

    lines.append("### RNA signatures (computed from TPM)\n\n")
    lines.append(f"- Output: `rna_signatures.csv`\n\n")
    lines.append(md_table_preview(sig_cohort[sig_cols + ['sample_id']].head(8)))
    lines.append("\n")

    lines.append("### Immune / subtype stratification\n\n")
    if by_files:
        lines.append("- Outputs:\n")
        for p in by_files:
            lines.append(f"  - `{p.name}`\n")
        lines.append("\n")
    if (out_dir / "signatures_and_thornsson_scores_by_immune_subtype_original.csv").exists():
        df = pd.read_csv(out_dir / "signatures_and_thornsson_scores_by_immune_subtype_original.csv")
        lines.append("#### By Thorsson Immune Subtype (C1–C6; original table)\n\n")
        lines.append(md_table_preview(df, max_rows=10))
        lines.append("\n")
    if (out_dir / "signatures_by_pam50_pipeline_table.csv").exists():
        df = pd.read_csv(out_dir / "signatures_by_pam50_pipeline_table.csv")
        lines.append("#### By PAM50 (pipeline advanced table)\n\n")
        lines.append(md_table_preview(df, max_rows=10))
        lines.append("\n")

    (out_dir / "report.md").write_text("".join(lines), encoding="utf-8")
    print(f"[OK] wrote: {out_dir / 'report.md'}")


if __name__ == "__main__":
    main()

