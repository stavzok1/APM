from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.qc import summarize_suspicious_columns, summarize_suspicious_nested_cells
from pipeline.Methylation.methylation_table import load_existing_methylation_probe_reference

from analysis.sanity._io import (
    ensure_out_dir,
    load_scratch_context,
    load_sample_ids,
    utc_run_id,
)


def _safe_read_csv(path: Path, nrows: int = 50_000) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, nrows=nrows)


def _qc_one_df(df: pd.DataFrame, *, name: str, out_dir: Path) -> None:
    cols = summarize_suspicious_columns(df)
    cols.to_csv(out_dir / f"{name}__suspicious_columns.csv", index=False)

    nested = summarize_suspicious_nested_cells(df)
    nested.to_csv(out_dir / f"{name}__nested_cells.csv", index=False)


def _scan_csv_files(
    files: Iterable[Path],
    *,
    out_dir: Path,
    prefix: str,
    nrows: int = 50_000,
    max_files: Optional[int] = None,
) -> None:
    """
    Scan many CSVs and write:
    - per-file summaries
    - aggregated counts of "almost always empty" columns across files
    """
    per_file_rows: List[Dict[str, object]] = []
    col_bad_counts: Dict[str, int] = {}
    col_seen_counts: Dict[str, int] = {}

    picked = []
    for p in files:
        picked.append(Path(p))
        if max_files is not None and len(picked) >= max_files:
            break

    for p in picked:
        if not p.is_file():
            continue
        try:
            df = _safe_read_csv(p, nrows=nrows)
        except Exception as exc:
            per_file_rows.append(
                {"file": str(p), "ok": False, "error": str(exc)[:300], "nrows_preview": 0}
            )
            continue

        qc = summarize_suspicious_columns(df)
        bad = qc[(qc["nan_frac"] > 0.99) | (qc["emptyish_frac"] > 0.99)]
        per_file_rows.append(
            {
                "file": str(p),
                "ok": True,
                "error": "",
                "nrows_preview": int(len(df)),
                "n_cols": int(df.shape[1]),
                "n_suspicious_cols": int(len(bad)),
            }
        )

        for _, row in qc.iterrows():
            c = str(row["column"])
            col_seen_counts[c] = col_seen_counts.get(c, 0) + 1
            if float(row.get("nan_frac", 0.0) or 0.0) > 0.99 or float(row.get("emptyish_frac", 0.0) or 0.0) > 0.99:
                col_bad_counts[c] = col_bad_counts.get(c, 0) + 1

    pd.DataFrame(per_file_rows).to_csv(out_dir / f"{prefix}__per_file_qc.csv", index=False)
    agg = pd.DataFrame(
        [
            {
                "column": c,
                "files_seen": int(col_seen_counts.get(c, 0)),
                "files_bad": int(col_bad_counts.get(c, 0)),
                "bad_frac": (col_bad_counts.get(c, 0) / col_seen_counts.get(c, 1)),
            }
            for c in sorted(col_seen_counts.keys())
        ]
    ).sort_values(["bad_frac", "files_bad", "files_seen"], ascending=[False, False, False])
    agg.to_csv(out_dir / f"{prefix}__aggregated_suspicious_columns.csv", index=False)


def _filter_files_for_samples(files: Iterable[Path], sample_ids: List[str]) -> List[Path]:
    """
    Keep only files whose basename contains one of the cohort sample ids.
    This prevents scanning historical outputs for other cohorts in shared data dirs.
    """
    sids = [str(s).strip() for s in sample_ids if str(s).strip()]
    out: List[Path] = []
    for p in files:
        name = Path(p).name
        if any(sid in name for sid in sids):
            out.append(Path(p))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch-json", type=str, required=True)
    args = ap.parse_args()

    ctx = load_scratch_context(Path(args.scratch_json))
    run_id = utc_run_id("structural_qc")
    out_dir = ensure_out_dir(run_id)

    sample_ids = load_sample_ids(ctx)
    if not sample_ids:
        raise SystemExit(
            "No sample_ids found. Expected either an existing methylation subset manifest "
            f"({ctx.methylation_subset_manifest}) or a prepare log next to {ctx.scratch_json}."
        )
    (out_dir / "sample_ids.txt").write_text("\n".join(sample_ids) + "\n", encoding="utf-8")

    # 1) Probe reference load projection sanity
    ref_path = ctx.methylation_working_dir / "reference" / "probe_annotations.parquet"
    if ref_path.exists():
        probe_ref = load_existing_methylation_probe_reference(ref_path)
        _qc_one_df(probe_ref, name="probe_reference_subset", out_dir=out_dir)

    # 2) Methylation cohort matrices
    cohort_dir = ctx.methylation_working_dir / "cohort"
    for fname in ("gene_meth_matrix.csv", "ccre_meth_matrix.csv", "lncrna_meth_matrix.csv", "atac_meth_matrix.csv"):
        p = cohort_dir / fname
        if p.exists():
            df = _safe_read_csv(p, nrows=200_000)
            _qc_one_df(df, name=f"methylation_{p.stem}", out_dir=out_dir)

    # 3) Methylation per-sample table shape + empties
    per_sample = ctx.methylation_working_dir / "per_sample"
    per_sample_rows: List[Dict] = []
    for sid in sample_ids:
        sdir = per_sample / sid
        rec: Dict[str, object] = {"sample_id": sid, "exists": bool(sdir.is_dir())}
        for suf in ("gene_meth.csv", "ccre_meth.csv", "lncrna_meth.csv", "atac_meth.csv", "probes.csv"):
            f = sdir / f"{sid}_{suf}"
            rec[f"has_{suf}"] = bool(f.is_file())
            if f.is_file() and suf != "probes.csv":
                df = _safe_read_csv(f, nrows=200_000)
                qc = summarize_suspicious_columns(df)
                bad = qc[(qc["nan_frac"] > 0.99) | (qc["emptyish_frac"] > 0.99)]
                rec[f"n_suspicious_cols_{suf}"] = int(len(bad))
            if f.is_file() and suf == "probes.csv":
                # probes can be huge; just count and do minimal empties scan
                df = _safe_read_csv(f, nrows=50_000)
                qc = summarize_suspicious_columns(df)
                bad = qc[(qc["nan_frac"] > 0.99) | (qc["emptyish_frac"] > 0.99)]
                rec["probe_preview_rows"] = int(len(df))
                rec["probe_preview_suspicious_cols"] = int(len(bad))
        per_sample_rows.append(rec)
    pd.DataFrame(per_sample_rows).to_csv(out_dir / "methylation_per_sample_presence_and_qc.csv", index=False)

    # 4) CNV gene calls (subset) — presence only (structural QC is heavy; do subset preview)
    cnv_rows: List[Dict] = []
    for sid in sample_ids:
        p = ctx.cnv_gene_tables_dir / f"{sid}_cnv_gene_calls_ascat3.csv"
        if not p.exists():
            p = ctx.cnv_gene_tables_dir / f"{sid}_cnv_gene_calls.csv"
        cnv_rows.append({"sample_id": sid, "cnv_gene_calls_path": str(p), "exists": bool(p.exists())})
    pd.DataFrame(cnv_rows).to_csv(out_dir / "cnv_gene_calls_presence.csv", index=False)

    # 5) SNV / SV / CNV / RPPA module outputs (preview scan)
    # SNV: many per-sample CSVs in PATHS.snv_output_dir/per_sample plus combined file.
    snv_root = ctx.snv_output
    if snv_root.exists():
        snv_files = list((snv_root / "per_sample").glob("*.csv")) if (snv_root / "per_sample").is_dir() else []
        combined = snv_root / "combined_snv_variants.csv"
        if combined.exists():
            snv_files = [combined] + snv_files
        snv_files = _filter_files_for_samples(snv_files, sample_ids) + ([combined] if combined.exists() else [])
        _scan_csv_files(snv_files, out_dir=out_dir, prefix="snv", nrows=50_000, max_files=120)

    # SV: scan the "final" outputs first (07), and chip/neojunction enrichments if present.
    sv_root = ctx.sv_output_root
    if sv_root.exists():
        for sub, pref in (
            ("07_final_sv_with_fimo", "sv_final"),
            ("08_neojunction_enriched", "sv_neojunction"),
            ("09_chip_enriched", "sv_chip_enriched"),
            ("02_processed_sv_csv", "sv_processed"),
        ):
            d = sv_root / sub
            if d.is_dir():
                files = _filter_files_for_samples(sorted(d.glob("*.csv")), sample_ids)
                _scan_csv_files(files, out_dir=out_dir, prefix=pref, nrows=50_000, max_files=120)

    # CNV annotated segments (not the per-gene tables above)
    cnv_ann = ctx.cnv_annotated_dir
    if cnv_ann.is_dir():
        files = _filter_files_for_samples(sorted(cnv_ann.glob("*.csv")), sample_ids)
        _scan_csv_files(files, out_dir=out_dir, prefix="cnv_annotated", nrows=50_000, max_files=120)

    # RPPA processed tables
    # When scratch json is methylation-only, this dir may still exist but wasn't regenerated.
    rppa_dir = ctx.rppa_output_dir
    if rppa_dir.is_dir():
        files = _filter_files_for_samples(sorted(rppa_dir.glob("*.csv")), sample_ids)
        _scan_csv_files(
            files,
            out_dir=out_dir,
            prefix="rppa",
            nrows=50_000,
            max_files=200,
        )

    print(f"[OK] Wrote structural QC to: {out_dir}")


if __name__ == "__main__":
    main()

