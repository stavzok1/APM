#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.qc._report_utils import ensure_out_dir, md_table_preview, write_csv  # noqa: E402
from pipeline.qc.base import Finding, Level  # noqa: E402
from pipeline.qc.context import load_context, load_sample_ids, to_tcga_sample_id  # noqa: E402
from pipeline.qc.within_modality import compute_within_modality_sanity, assert_within_modality_sanity  # noqa: E402
from pipeline.qc.sv_vs_rna import compute_sv_disruption_vs_rna, assert_sv_disruption_vs_rna  # noqa: E402
from pipeline.qc.loh_vs_neoantigen import (  # noqa: E402
    compute_hla_loh_vs_neoantigen_and_infiltration,
    assert_hla_loh_vs_neoantigen,
)


def _utc_run_id(prefix: str) -> str:
    from analysis.qc._context import utc_run_id

    return utc_run_id(prefix)


def _findings_to_df(findings: List[Finding]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "check": f.check,
                "level": f.level.value,
                "message": f.message,
                "hint": f.hint or "",
                "context_json": json.dumps(f.context, sort_keys=True),
            }
            for f in findings
        ]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch-json", type=str, required=True)
    ap.add_argument(
        "--build-rppa-if-missing",
        action="store_true",
        help="If ctx.rppa_output_dir is empty, run the RPPA pipeline to populate it.",
    )
    args = ap.parse_args()

    ctx = load_context(Path(args.scratch_json))
    run_id = _utc_run_id("full_qc")
    out_dir = ensure_out_dir(run_id)

    sample_ids = load_sample_ids(ctx)
    if not sample_ids:
        raise SystemExit("No sample_ids found for cohort (manifest or prepare log).")

    if args.build_rppa_if_missing:
        out_dir = Path(ctx.rppa_output_dir)
        has_any = out_dir.exists() and (any(out_dir.glob("*.csv")) or any(out_dir.glob("*.parquet")))
        if not has_any:
            from pipeline.rppa.rppa_main import run_rppa_pipeline
            from pipeline.rppa.rppa_config import RPPAPathConfig
            from pipeline.config import PATHS

            out_dir.mkdir(parents=True, exist_ok=True)
            run_rppa_pipeline(
                sample_dir=RPPAPathConfig.sample_data_dir,
                annotation_path=(Path(PATHS.annotations_dir) / "rppa" / "TCGA_antibodies_descriptions.gencode.v36.tsv"),
                output_dir=out_dir,
                metadata_path=RPPAPathConfig.sample_metadata_tumor,
                rna_expr=None,
                save_outputs=True,
            )

    # --- Compute ---
    wm_samp, wm_gene, wm_sum = compute_within_modality_sanity(ctx=ctx, sample_ids=sample_ids, out_dir=out_dir)
    sv_samp, sv_gene, sv_sum = compute_sv_disruption_vs_rna(ctx=ctx, sample_ids=sample_ids, out_dir=out_dir)
    loh_samp, loh_gene, loh_sum = compute_hla_loh_vs_neoantigen_and_infiltration(ctx=ctx, sample_ids=sample_ids, out_dir=out_dir)

    # --- Assert (pass/warn/fail list) ---
    findings: List[Finding] = []
    findings += assert_within_modality_sanity(wm_sum)
    findings += assert_sv_disruption_vs_rna(sv_sum)
    findings += assert_hla_loh_vs_neoantigen(loh_sum)

    # Cohort-size / “smoke run” disclosure (this matters for interpretation).
    if int(len(sample_ids)) < 200:
        findings.append(
            Finding(
                check="qc.cohort_size_low_power",
                level=Level.WARN,
                message=f"Cohort size is n={len(sample_ids)}; treat this as smoke-QC (low statistical power).",
                context={"n_samples": int(len(sample_ids))},
                hint="Run QC on a larger cohort before interpreting correlation directions/effect sizes.",
            )
        )

    # We append some findings after computing SV summaries; keep the initial df local only.
    findings_df = _findings_to_df(findings)

    # --- Per-sample QC card (numeric + some booleans) ---
    card = pd.DataFrame({"sample_id": sample_ids})
    card["tcga_sample_id"] = card["sample_id"].map(to_tcga_sample_id)
    if not wm_samp.empty:
        card = card.merge(wm_samp, on="sample_id", how="left")
    if not sv_samp.empty:
        card = card.merge(sv_samp, on="sample_id", how="left", suffixes=("", "__sv"))
    if not loh_samp.empty:
        card = card.merge(loh_samp.drop(columns=["tcga_sample_id"], errors="ignore"), on="sample_id", how="left")
    write_csv(card, out_dir / "per_sample_qc_card.csv")

    # --- Per-gene card ---
    gene_card = pd.DataFrame()
    if not sv_gene.empty:
        gene_card = sv_gene.copy()
    write_csv(gene_card, out_dir / "per_gene_qc_card.csv")

    # --- SV↔RNA delta inspection (typical effect sizes across genes) ---
    sv_delta_summary = {}
    sv_delta_table = pd.DataFrame()
    if not sv_gene.empty:
        sv_delta_table = sv_gene.copy()
        if "strat_kind" in sv_delta_table.columns and "strat_value" in sv_delta_table.columns:
            sv_delta_table = sv_delta_table[
                (sv_delta_table["strat_kind"].astype(str) == "ANY")
                & (sv_delta_table["strat_value"].astype(str) == "ANY")
            ].copy()
        # Keep only genes with adequate support for adjusted comparisons
        m = (
            pd.to_numeric(sv_delta_table.get("n_sv_hit"), errors="coerce").fillna(0) >= 5
        ) & (
            pd.to_numeric(sv_delta_table.get("n_sv_no_hit"), errors="coerce").fillna(0) >= 10
        )
        sv_supported = sv_delta_table[m].copy()
        cols = [
            c
            for c in (
                "gene",
                "n_sv_hit",
                "n_sv_no_hit",
                "delta_median_unadjusted",
                "delta_median_cnv_neutral",
                "reg_beta_sv_hit",
                "reg_beta_cnv_log2",
                "reg_beta_ifn",
                "reg_r2",
                "n_cnv_neutral_sv_hit",
                "n_cnv_neutral_no_sv",
            )
            if c in sv_supported.columns
        ]
        sv_supported = sv_supported[cols]
        write_csv(sv_supported, out_dir / "sv_vs_rna__supported_genes.csv")

        # Stratified views (informational; can be sparse).
        if "strat_kind" in sv_gene.columns and "strat_value" in sv_gene.columns:
            sv_by = sv_gene[(sv_gene["strat_kind"].astype(str) != "ANY")].copy()
            if not sv_by.empty:
                cols_by = [
                    c
                    for c in (
                        "strat_kind",
                        "strat_value",
                        "gene",
                        "n_sv_hit",
                        "n_sv_no_hit",
                        "delta_median_unadjusted",
                        "delta_median_cnv_neutral",
                        "reg_beta_sv_hit",
                        "reg_beta_cnv_log2",
                        "reg_beta_ifn",
                        "reg_r2",
                    )
                    if c in sv_by.columns
                ]
                write_csv(sv_by[cols_by], out_dir / "sv_vs_rna__per_gene_stratified.csv")

        def _q(x: pd.Series) -> dict:
            x = pd.to_numeric(x, errors="coerce")
            x = x[x.notna()]
            if x.empty:
                return {}
            return {
                "n": int(len(x)),
                "min": float(x.min()),
                "p10": float(x.quantile(0.10)),
                "p25": float(x.quantile(0.25)),
                "median": float(x.median()),
                "p75": float(x.quantile(0.75)),
                "p90": float(x.quantile(0.90)),
                "max": float(x.max()),
            }

        for c in ("delta_median_unadjusted", "delta_median_cnv_neutral", "reg_beta_sv_hit"):
            if c in sv_supported.columns:
                sv_delta_summary[c] = _q(sv_supported[c])

        # Directionality rates (how often SV-hit looks lower)
        if "delta_median_unadjusted" in sv_supported.columns:
            d = pd.to_numeric(sv_supported["delta_median_unadjusted"], errors="coerce")
            sv_delta_summary["frac_negative_unadjusted"] = float((d < 0).mean()) if d.notna().any() else np.nan
        if "delta_median_cnv_neutral" in sv_supported.columns:
            d = pd.to_numeric(sv_supported["delta_median_cnv_neutral"], errors="coerce")
            sv_delta_summary["frac_negative_cnv_neutral"] = float((d < 0).mean()) if d.notna().any() else np.nan

        # Convert key cohort-level SV↔RNA interpretation points into explicit findings.
        try:
            n_supported = int(len(sv_supported))
            if n_supported < 15:
                findings.append(
                    Finding(
                        check="sv_vs_rna.support_thin",
                        level=Level.WARN,
                        message=f"Only {n_supported} genes meet SV↔RNA minimum-support thresholds (hit>=5, no-hit>=10).",
                        context={"n_supported_genes": n_supported},
                        hint="This cohort may be too small or too SV-sparse for gene-level SV↔RNA QC on most APM genes.",
                    )
                )
        except Exception:
            pass

        try:
            if "frac_negative_unadjusted" in sv_delta_summary:
                frac = float(sv_delta_summary.get("frac_negative_unadjusted"))
                if frac == frac and frac < 0.5:
                    findings.append(
                        Finding(
                            check="sv_vs_rna.directionality_unexpected",
                            level=Level.WARN,
                            message=f"SV↔RNA unadjusted direction is unexpected for many genes (frac_negative={frac:.2f}).",
                            context={"frac_negative_unadjusted": frac},
                            hint="Check SV class/overlap definitions; consider separating promoter/CDS vs intronic hits, and interpret with CNV/IFN/purity stratification.",
                        )
                    )
            q = sv_delta_summary.get("delta_median_unadjusted") or {}
            med = q.get("median")
            if med is not None:
                med = float(med)
                if med == med and med > 0:
                    findings.append(
                        Finding(
                            check="sv_vs_rna.delta_positive_unadjusted",
                            level=Level.WARN,
                            message=f"Median unadjusted SV↔RNA delta is positive (median={med:.3f}); expected negative for true disruptive SVs.",
                            context={"median_delta_median_unadjusted": med},
                            hint="Often indicates confounding (CNV/IFN/subtype) or overly-broad SV-hit definitions; regression-adjusted betas may be more informative.",
                        )
                    )
        except Exception:
            pass

    # Now that we've appended post-compute findings, materialize and save findings.csv.
    findings_df = _findings_to_df(findings)
    write_csv(findings_df, out_dir / "findings.csv")

    # --- Summary JSON + report ---
    summary = {
        "run_id": run_id,
        "n_samples": int(len(sample_ids)),
        "scratch_json": str(Path(args.scratch_json)),
        "within_modality": wm_sum,
        "sv_vs_rna": sv_sum,
        "hla_loh_vs_neoantigen": loh_sum,
        "sv_vs_rna_delta_summary": sv_delta_summary,
        "n_findings": int(len(findings)),
        "n_fail": int((findings_df["level"] == Level.FAIL.value).sum()) if not findings_df.empty else 0,
        "n_warn": int((findings_df["level"] == Level.WARN.value).sum()) if not findings_df.empty else 0,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines: List[str] = []
    lines.append("## Full QC report\n\n")
    lines.append(f"- **run_id**: `{run_id}`\n")
    lines.append(f"- **n_samples**: {len(sample_ids)}\n")
    lines.append(f"- **scratch_json**: `{Path(args.scratch_json)}`\n\n")

    lines.append("### Findings (pass/warn/fail)\n\n")
    if not findings_df.empty:
        lines.append(md_table_preview(findings_df.head(20)))
    else:
        lines.append("_(none)_\n")
    lines.append("\n")

    lines.append("### Per-sample QC card (preview)\n\n")
    lines.append(md_table_preview(card.head(12)))
    lines.append("\n")

    lines.append("### SV ↔ RNA: typical deltas across genes (supported genes only)\n\n")
    if sv_delta_summary:
        lines.append("Derived from `sv_vs_rna__supported_genes.csv`.\n\n")
        # small, readable table
        rows = []
        for k in ("delta_median_unadjusted", "delta_median_cnv_neutral", "reg_beta_sv_hit"):
            if k in sv_delta_summary:
                q = sv_delta_summary[k]
                rows.append(
                    {
                        "metric": k,
                        "n_genes": q.get("n"),
                        "median": q.get("median"),
                        "p25": q.get("p25"),
                        "p75": q.get("p75"),
                        "min": q.get("min"),
                        "max": q.get("max"),
                    }
                )
        if rows:
            lines.append(md_table_preview(pd.DataFrame(rows), max_rows=10, max_cols=10))
            lines.append("\n")
        if "frac_negative_unadjusted" in sv_delta_summary:
            lines.append(f"- frac_negative_unadjusted: `{sv_delta_summary.get('frac_negative_unadjusted')}`\n")
        if "frac_negative_cnv_neutral" in sv_delta_summary:
            lines.append(f"- frac_negative_cnv_neutral: `{sv_delta_summary.get('frac_negative_cnv_neutral')}`\n")
        lines.append("\n")

    p_by = out_dir / "sv_vs_rna__per_gene_stratified.csv"
    if p_by.exists():
        try:
            by = pd.read_csv(p_by)
            if not by.empty:
                lines.append("### SV ↔ RNA stratified (region + SVTYPE; informational; can be sparse)\n\n")
                lines.append(md_table_preview(by.head(20)))
                lines.append("\n")
        except Exception:
            pass

    (out_dir / "report.md").write_text("".join(lines), encoding="utf-8")
    print(f"[OK] wrote: {out_dir / 'report.md'}")


if __name__ == "__main__":
    main()

