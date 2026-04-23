#!/usr/bin/env python3
"""
Summarize BRCA isoform matrix filtering when using **genes-only** transcript rules
(``--genes-only`` on ``filter_isoform_expression_brca.py``).

Streams the isoform TPM table in chunks (same column whitelist as the filter script) and reports:
  - rows read vs rows that pass the transcript stem filter
  - distinct transcripts (stems) retained
  - genes represented, transcripts per gene, MANE-select counts per gene
  - global tag histogram

Does **not** load lncRNA feature tables. For a quick smoke pass, use ``--max-chunks``.

Examples::

    .venv/bin/python3 pipeline/RNA_exp/report_isoform_brca_genes.py --out-json /tmp/iso_gene_report.json
    .venv/bin/python3 pipeline/RNA_exp/report_isoform_brca_genes.py --max-chunks 3
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Set

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PIPELINE_GENE_PANEL  # noqa: E402
from pipeline.sample_ids import normalize_tcga_id  # noqa: E402


def _load_filter_isoform_module():
    p = Path(__file__).resolve().parent / "filter_isoform_expression_brca.py"
    spec = importlib.util.spec_from_file_location("filter_isoform_expression_brca", p)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_fiso = _load_filter_isoform_module()
_stem_ensembl = _fiso._stem_ensembl
_transcript_rows = _fiso._transcript_rows
build_annotation_and_allowed_stems = _fiso.build_annotation_and_allowed_stems
load_brca_isoform_column_whitelist = _fiso.load_brca_isoform_column_whitelist
load_cohort_participants = _fiso.load_cohort_participants


def _stem_to_gene_for_allowed(
    allowed: Set[str],
    probemap_path: Path,
    primary_features_path: Path,
    gene_union: Set[str],
) -> Dict[str, str]:
    """Map Ensembl transcript stem -> HGNC-style gene symbol (probemap first, then GENCODE primary)."""
    pm = pd.read_csv(probemap_path, sep="\t", low_memory=False)
    if "id" not in pm.columns or "gene" not in pm.columns:
        raise ValueError(f"probemap needs id,gene: {list(pm.columns)}")
    pm = pm.copy()
    pm["stem"] = pm["id"].astype(str).map(_stem_ensembl)
    out: Dict[str, str] = {}
    for _, row in pm.iterrows():
        st = str(row["stem"]).strip()
        g = str(row.get("gene", "") or "").strip()
        if not st or not g or g not in gene_union:
            continue
        if st in allowed and st not in out:
            out[st] = g

    trp = _transcript_rows(primary_features_path)
    if not trp.empty and "gene_name" in trp.columns:
        for _, row in trp.iterrows():
            st = str(row["stem"]).strip()
            if st not in allowed or st in out:
                continue
            g = str(row.get("gene_name", "") or "").strip()
            if g and g in gene_union:
                out[st] = g
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--isoform-tsv", type=Path, default=ROOT / "data/RNAexp_TCGA/tcga_rsem_isoform_tpm")
    ap.add_argument(
        "--brca-clinical",
        type=Path,
        default=ROOT / "annotations/BRCA_clinical",
    )
    ap.add_argument(
        "--cohort-summary",
        type=Path,
        default=None,
        help="Optional cohort_summary.tsv (participant filter).",
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
    ap.add_argument("--chunksize", type=int, default=5000)
    ap.add_argument("--max-chunks", type=int, default=0, help="0 = scan entire matrix.")
    ap.add_argument("--out-json", type=Path, default=None, help="Write full report JSON here.")
    ap.add_argument(
        "--top-genes",
        type=int,
        default=40,
        help="Include top N genes by transcript count in the printed summary.",
    )
    args = ap.parse_args()

    matrix_ids, brca_participants = load_brca_isoform_column_whitelist(args.brca_clinical)
    cohort_p = load_cohort_participants(args.cohort_summary)

    with args.isoform_tsv.open("rb") as fh:
        header_line = fh.readline().decode("utf-8", errors="replace")
    cols = header_line.rstrip("\r\n").split("\t")
    id_col = cols[0]
    kept_sample_cols: List[str] = []
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
            continue
        part = tid.participant
        if part is None or part not in brca_participants:
            continue
        if cohort_p is not None and part not in cohort_p:
            continue
        kept_sample_cols.append(c)

    if not kept_sample_cols:
        raise SystemExit("No sample columns after BRCA / cohort filters.")

    allowed, ann = build_annotation_and_allowed_stems(
        PIPELINE_GENE_PANEL,
        args.probemap,
        args.primary_features,
        args.lnc_features,
        args.lnc_window,
        genes_only=True,
    )
    ann_by = ann.set_index("stem", drop=False)

    prim = pd.read_csv(args.primary_features, low_memory=False)
    prim_syms = (
        set(prim["gene_name"].dropna().astype(str).str.strip())
        if "gene_name" in prim.columns
        else set()
    )
    gene_union = {str(g).strip() for g in PIPELINE_GENE_PANEL if str(g).strip()} | prim_syms
    stem_gene = _stem_to_gene_for_allowed(allowed, args.probemap, args.primary_features, gene_union)

    usecols = [id_col] + kept_sample_cols
    usecols_set = set(usecols)

    n_read = 0
    n_kept_rows = 0
    retained_stems: Set[str] = set()

    for ci, chunk in enumerate(
        pd.read_csv(
            args.isoform_tsv,
            sep="\t",
            usecols=lambda c: c in usecols_set,
            chunksize=args.chunksize,
            low_memory=False,
        )
    ):
        if args.max_chunks and ci >= args.max_chunks:
            break
        n_read += len(chunk)
        chunk = chunk.rename(columns={id_col: "transcript_id"})
        chunk = chunk.assign(
            stem=chunk["transcript_id"].astype(str).map(_stem_ensembl),
        )
        sub = chunk.loc[chunk["stem"].isin(allowed)].copy()
        if sub.empty:
            continue
        n_kept_rows += len(sub)
        retained_stems.update(sub["stem"].astype(str).unique())

    tag_counts: Counter[str] = Counter()
    mane_by_gene: Dict[str, int] = defaultdict(int)
    stems_by_gene: Dict[str, Set[str]] = defaultdict(set)
    n_mane_select_transcripts = 0
    for st in retained_stems:
        row = ann_by.loc[st] if st in ann_by.index else None
        tag = ""
        is_mane = False
        if row is not None:
            tv = row.get("tag")
            if tv is None or (isinstance(tv, float) and pd.isna(tv)):
                tag = ""
            else:
                tag = str(tv).strip()
            im = row.get("is_MANE", False)
            if isinstance(im, str):
                is_mane = im.strip().lower() in ("true", "1", "yes")
            else:
                is_mane = bool(im) and not (isinstance(im, float) and pd.isna(im))
        tag_key = tag if tag else "(empty)"
        tag_counts[tag_key] += 1
        if is_mane:
            n_mane_select_transcripts += 1
        g = stem_gene.get(st, "") or "(unmapped)"
        stems_by_gene[g].add(st)
        if is_mane and g != "(unmapped)":
            mane_by_gene[g] += 1

    per_gene: List[Dict[str, Any]] = []
    for g, stems in sorted(stems_by_gene.items(), key=lambda x: (-len(x[1]), x[0])):
        if g == "(unmapped)":
            continue
        per_gene.append(
            {
                "gene": g,
                "n_transcripts": len(stems),
                "n_mane_select_transcripts": int(mane_by_gene[g]),
            }
        )

    tx_per_gene_hist: Dict[str, int] = {"1": 0, "2-5": 0, "6-15": 0, ">15": 0}
    for row in per_gene:
        n = int(row["n_transcripts"])
        if n == 1:
            tx_per_gene_hist["1"] += 1
        elif n <= 5:
            tx_per_gene_hist["2-5"] += 1
        elif n <= 15:
            tx_per_gene_hist["6-15"] += 1
        else:
            tx_per_gene_hist[">15"] += 1

    report: Dict[str, Any] = {
        "isoform_tsv": str(args.isoform_tsv),
        "genes_only": True,
        "n_header_columns": len(cols),
        "n_brca_cohort_sample_columns": len(kept_sample_cols),
        "n_allowed_transcript_stems_panel_genes_only": len(allowed),
        "chunksize": args.chunksize,
        "max_chunks": args.max_chunks or None,
        "n_matrix_rows_read": n_read,
        "n_matrix_rows_retained": n_kept_rows,
        "n_distinct_transcripts_retained": len(retained_stems),
        "n_genes_with_at_least_one_transcript": len(
            [g for g in stems_by_gene if g != "(unmapped)"]
        ),
        "n_transcripts_unmapped_gene_symbol": len(stems_by_gene.get("(unmapped)", set())),
        "n_mane_select_transcripts_total": n_mane_select_transcripts,
        "n_genes_with_at_least_one_mane_transcript": sum(
            1 for g in stems_by_gene if g != "(unmapped)" and mane_by_gene[g] > 0
        ),
        "tag_counts_on_distinct_transcripts": dict(tag_counts.most_common()),
        "genes_by_transcript_count_bucket": tx_per_gene_hist,
        "per_gene": per_gene,
    }

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_json, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
        print(f"Wrote {args.out_json}")

    topn = max(0, int(args.top_genes))
    print(json.dumps({k: report[k] for k in report if k != "per_gene"}, indent=2))
    if topn and per_gene:
        print(f"\nTop {topn} genes by transcript count (genes-only panel):")
        for row in per_gene[:topn]:
            print(
                f"  {row['gene']}: transcripts={row['n_transcripts']}, "
                f"mane_select_transcripts={row['n_mane_select_transcripts']}"
            )


if __name__ == "__main__":
    main()
