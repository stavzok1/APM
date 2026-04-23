#!/usr/bin/env python3
"""
Filter a pan-cancer isoform expression table (rows = Ensembl transcripts, columns = samples)
to BRCA participants, then to analysis-relevant transcripts, and attach GENCODE transcript tags.

Default input matches the Xena-style RSEM isoform file shipped as
``data/RNAexp_TCGA/tcga_rsem_isoform_tpm`` (first column = transcript id, header row with
``sample`` + ``TCGA-…`` columns).

The matrix is read in **chunks** to limit RAM. Outputs:
  - Filtered table (parquet or TSV)
  - JSON QC summary (column coverage, row coverage, annotation match rates)

Examples (from repo root)::

    .venv/bin/python3 pipeline/RNA_exp/filter_isoform_expression_brca.py \\
        --max-rows 20000 --out-table data/RNAexp_TCGA/isoform_brca_panel_sample.parquet

    .venv/bin/python3 pipeline/RNA_exp/filter_isoform_expression_brca.py \\
        --out-table data/RNAexp_TCGA/isoform_brca_panel.parquet \\
        --qc-json data/RNAexp_TCGA/isoform_brca_panel_qc.json

    # Genes-only transcript universe (no lncRNA tables); assess retention with:
    .venv/bin/python3 pipeline/RNA_exp/report_isoform_brca_genes.py --out-json /tmp/iso_gene_report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PIPELINE_GENE_PANEL  # noqa: E402
from pipeline.sample_ids import normalize_tcga_id, tcga_sample_type_two_digit  # noqa: E402


def _stem_ensembl(s: str) -> str:
    s = str(s).strip()
    if not s:
        return ""
    return s.split(".", 1)[0]


def load_brca_isoform_column_whitelist(
    brca_clinical_tsv: Path,
    sample_types: tuple[str, ...] = ("01", "10", "11"),
) -> tuple[Set[str], Set[str]]:
    """
    From BRCA clinical ``sampleID`` rows, build (1) strings that may appear as Xena matrix
    column headers (primary tumor 01, blood normal 10, solid normal 11) and (2) participants.

    This avoids matching **all** pan-cancer columns: only specimen barcodes implied by the
    BRCA manifest (plus ``-01`` / ``-01A`` style aliases) are considered.
    """
    df = pd.read_csv(brca_clinical_tsv, sep="\t", usecols=["sampleID"], low_memory=False)
    matrix_ids: Set[str] = set()
    participants: Set[str] = set()
    for sid in df["sampleID"].dropna().astype(str).str.strip():
        if not sid.startswith("TCGA"):
            continue
        code = tcga_sample_type_two_digit(sid)
        if code not in sample_types:
            continue
        matrix_ids.add(sid)
        tid = normalize_tcga_id(sid)
        if tid.sample:
            matrix_ids.add(tid.sample)
        if tid.sample_vial:
            matrix_ids.add(tid.sample_vial)
        if tid.participant:
            participants.add(tid.participant)
    return matrix_ids, participants


def load_cohort_participants(tsv: Optional[Path]) -> Optional[Set[str]]:
    if tsv is None or not tsv.is_file():
        return None
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if "participant" not in df.columns:
        return None
    return set(df["participant"].dropna().astype(str).str.strip())


def _transcript_rows(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if "feature" not in df.columns or "transcript_id" not in df.columns:
        return pd.DataFrame()
    tr = df[df["feature"].astype(str) == "transcript"].copy()
    tr = tr[tr["transcript_id"].notna() & (tr["transcript_id"].astype(str).str.len() > 3)]
    tr["stem"] = tr["transcript_id"].astype(str).map(_stem_ensembl)
    return tr


def build_annotation_and_allowed_stems(
    panel: List[str],
    probemap_path: Path,
    primary_features_path: Path,
    lnc_features_path: Path,
    lnc_window_path: Path,
    *,
    genes_only: bool = False,
) -> Tuple[Set[str], pd.DataFrame]:
    """
    Build allowed Ensembl transcript stems and a one-row-per-stem annotation table.

    When ``genes_only=True``, lncRNA window / lncRNA GENCODE features are ignored: the
    gene universe is **panel + primary gene names** only (no lncRNA-derived symbols or
    ENSG stems from the lncRNA distance table).
    """
    panel_set = {str(g).strip() for g in panel if str(g).strip()}

    lnc_symbols: Set[str] = set()
    lnc_gene_stems: Set[str] = set()
    if not genes_only:
        lncw = pd.read_csv(lnc_window_path, low_memory=False)
        for c in ("gene_name", "gene_id"):
            if c not in lncw.columns:
                continue
            for v in lncw[c].dropna().astype(str):
                v = v.strip()
                if not v:
                    continue
                if v.startswith("ENSG"):
                    lnc_gene_stems.add(_stem_ensembl(v))
                else:
                    lnc_symbols.add(v)

    prim = pd.read_csv(primary_features_path, low_memory=False)
    prim_syms = (
        set(prim["gene_name"].dropna().astype(str).str.strip())
        if "gene_name" in prim.columns
        else set()
    )
    gene_union = panel_set | prim_syms | lnc_symbols

    pm = pd.read_csv(probemap_path, sep="\t", low_memory=False)
    if "id" not in pm.columns or "gene" not in pm.columns:
        raise ValueError(f"probemap needs id,gene columns: {list(pm.columns)}")
    pm = pm.copy()
    pm["stem"] = pm["id"].astype(str).map(_stem_ensembl)
    stems_prob = set(pm.loc[pm["gene"].astype(str).isin(gene_union), "stem"])

    tr_p = _transcript_rows(primary_features_path)
    tr_l = (
        pd.DataFrame()
        if genes_only
        else _transcript_rows(lnc_features_path)
    )
    tr_all = pd.concat([tr_p, tr_l], ignore_index=True)
    if tr_all.empty:
        stems_feat: Set[str] = set()
    else:
        gn = (
            tr_all["gene_name"].astype(str).isin(gene_union)
            if "gene_name" in tr_all.columns
            else pd.Series(False, index=tr_all.index)
        )
        gid_ok = pd.Series(False, index=tr_all.index)
        if "gene_id" in tr_all.columns:
            gid_ok = tr_all["gene_id"].astype(str).map(_stem_ensembl).isin(lnc_gene_stems)
        stems_feat = set(tr_all.loc[gn | gid_ok, "stem"])

    allowed = stems_prob | stems_feat

    ann = tr_all[tr_all["stem"].isin(allowed)].copy()
    if ann.empty:
        ann = pd.DataFrame(columns=["stem", "tag", "is_MANE"])
    else:
        if "tag" not in ann.columns:
            ann["tag"] = ""
        if "is_MANE" not in ann.columns:
            ann["is_MANE"] = False
        ann["is_MANE"] = ann["is_MANE"].astype(str).str.lower().isin(("true", "1", "yes"))
        ann = ann[["stem", "tag", "is_MANE"]].sort_values(
            ["stem", "is_MANE"], ascending=[True, False]
        )
        ann = ann.drop_duplicates(subset=["stem"], keep="first")

    return allowed, ann


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--isoform-tsv",
        type=Path,
        default=ROOT / "data/RNAexp_TCGA/tcga_rsem_isoform_tpm",
        help="Input matrix (tab-separated; row0 = header).",
    )
    ap.add_argument(
        "--out-table",
        type=Path,
        default=ROOT / "data/RNAexp_TCGA/isoform_tpm_brca_panel.parquet",
    )
    ap.add_argument(
        "--qc-json",
        type=Path,
        default=ROOT / "data/RNAexp_TCGA/isoform_tpm_brca_panel_qc.json",
    )
    ap.add_argument(
        "--brca-clinical",
        type=Path,
        default=ROOT / "annotations/BRCA_clinical",
        help="TSV with TCGA BRCA ``sampleID`` column.",
    )
    ap.add_argument(
        "--cohort-summary",
        type=Path,
        default=None,
        help="Optional cohort_summary.tsv: intersect participants with BRCA clinical.",
    )
    ap.add_argument(
        "--probemap",
        type=Path,
        default=ROOT / "annotations/RNA/probeMap_gencode.v23.annotation.transcript.probemap",
    )
    ap.add_argument(
        "--primary-features",
        type=Path,
        default=ROOT / "data/primary_genes_all_features.csv",
    )
    ap.add_argument(
        "--lnc-features",
        type=Path,
        default=ROOT / "data/lncRNAs_genes_all_features.csv",
    )
    ap.add_argument(
        "--lnc-window",
        type=Path,
        default=ROOT / "data/lncRNA_matching/lncRNAs_with_genes_1000000bp.csv",
    )
    ap.add_argument(
        "--genes-only",
        action="store_true",
        help="Restrict transcript universe to panel + primary protein-coding GENCODE genes "
        "(exclude lncRNA tables and lncRNA-derived gene symbols).",
    )
    ap.add_argument("--chunksize", type=int, default=5000, help="Rows per read chunk.")
    ap.add_argument(
        "--max-rows",
        type=int,
        default=0,
        help="Stop after this many **kept** transcript rows (0 = no limit). For dry runs.",
    )
    ap.add_argument(
        "--max-chunks",
        type=int,
        default=0,
        help="Process at most this many input chunks (0 = read entire file). For smoke tests.",
    )
    ap.add_argument(
        "--format",
        choices=("parquet", "tsv"),
        default="parquet",
    )
    args = ap.parse_args()

    matrix_ids, brca_participants = load_brca_isoform_column_whitelist(args.brca_clinical)
    cohort_p = load_cohort_participants(args.cohort_summary)

    with args.isoform_tsv.open("rb") as fh:
        header_line = fh.readline().decode("utf-8", errors="replace")
    cols = header_line.rstrip("\r\n").split("\t")
    id_col = cols[0]
    kept_sample_cols: List[str] = []
    dropped_tcga = 0
    for c in cols[1:]:
        c = str(c).strip()
        if not c.startswith("TCGA"):
            continue
        tid = normalize_tcga_id(c)
        if (
            c not in matrix_ids
            and (tid.sample or "") not in matrix_ids
            and (tid.sample_vial or "") not in matrix_ids
        ):
            dropped_tcga += 1
            continue
        part = tid.participant
        if part is None or part not in brca_participants:
            dropped_tcga += 1
            continue
        if cohort_p is not None and part not in cohort_p:
            dropped_tcga += 1
            continue
        kept_sample_cols.append(c)

    allowed_stems, ann = build_annotation_and_allowed_stems(
        PIPELINE_GENE_PANEL,
        args.probemap,
        args.primary_features,
        args.lnc_features,
        args.lnc_window,
        genes_only=args.genes_only,
    )

    usecols = [id_col] + kept_sample_cols
    qc: Dict[str, object] = {
        "input": str(args.isoform_tsv),
        "n_header_columns": len(cols),
        "n_brca_manifest_matrix_id_aliases": len(matrix_ids),
        "n_brca_cohort_sample_columns": len(kept_sample_cols),
        "n_tcga_columns_dropped": dropped_tcga,
        "genes_only": bool(args.genes_only),
        "n_allowed_transcript_stems": len(allowed_stems),
        "n_annotation_rows_unique_stem": int(len(ann)),
        "n_rows_read": 0,
        "n_rows_kept_after_transcript_filter": 0,
        "n_annotation_matched_on_output": 0,
    }

    if not kept_sample_cols:
        raise SystemExit(
            "No TCGA sample columns passed BRCA (+ optional cohort) filters. "
            "Check --brca-clinical sampleIDs vs matrix header barcodes."
        )

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    kept_total = 0
    ann_matched = 0
    pq_writer: Optional[object] = None
    first_csv = True

    for ci, chunk in enumerate(
        pd.read_csv(
            args.isoform_tsv,
            sep="\t",
            usecols=lambda c: c in set(usecols),
            chunksize=args.chunksize,
            low_memory=False,
        )
    ):
        if args.max_chunks and ci >= args.max_chunks:
            break
        chunk = chunk.rename(columns={id_col: "transcript_id"})
        chunk = chunk.assign(
            stem=chunk["transcript_id"].astype(str).map(_stem_ensembl)
        )
        sub = chunk[chunk["stem"].isin(allowed_stems)].copy()
        qc["n_rows_read"] = int(qc["n_rows_read"]) + len(chunk)
        if sub.empty:
            continue
        sub = sub.merge(ann, on="stem", how="left")
        ann_matched += int(
            (
                sub["tag"].notna()
                & (sub["tag"].astype(str).str.strip().str.len() > 0)
            ).sum()
        )
        sub = sub.drop(columns=["stem"])
        if args.max_rows and args.max_rows > 0:
            room = args.max_rows - kept_total
            if room <= 0:
                break
            if len(sub) > room:
                sub = sub.iloc[:room].copy()
        kept_total += len(sub)
        if args.format == "parquet":
            import pyarrow as pa
            import pyarrow.parquet as pq

            table = pa.Table.from_pandas(sub, preserve_index=False)
            if pq_writer is None:
                pq_writer = pq.ParquetWriter(
                    args.out_table, table.schema, compression="snappy"
                )
            pq_writer.write_table(table)
        else:
            sub.to_csv(
                args.out_table,
                sep="\t",
                index=False,
                mode="w" if first_csv else "a",
                header=first_csv,
            )
            first_csv = False
        if args.max_rows and kept_total >= args.max_rows:
            break

    if pq_writer is not None:
        pq_writer.close()

    qc["n_rows_kept_after_transcript_filter"] = kept_total
    qc["n_annotation_matched_on_output"] = ann_matched

    args.qc_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.qc_json, "w", encoding="utf-8") as fh:
        json.dump(qc, fh, indent=2)
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
