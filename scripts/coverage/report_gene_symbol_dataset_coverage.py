#!/usr/bin/env python3
"""
Cross-check panel gene symbols in ``pipeline/config.py`` against datasets on disk
and summarize alias coverage (HGNC/UCSC/legacy) for RNA-style gaps.

Writes a markdown report under ``analysis/output/`` by default.

Usage (from repo root):
  python scripts/coverage/report_gene_symbol_dataset_coverage.py
  python scripts/coverage/report_gene_symbol_dataset_coverage.py --out /tmp/coverage.md
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


def _gencode_gene_names(path: Path, max_chunks: int = 40) -> Set[str]:
    if not path.is_file():
        return set()
    names: Set[str] = set()
    if path.suffix.lower() in (".parquet", ".pq"):
        try:
            import pyarrow.parquet as pq

            pf = pq.ParquetFile(path)
            for i in range(min(max_chunks, pf.num_row_groups)):
                tbl = pf.read_row_group(i, columns=["gene_name"])
                df = tbl.to_pandas()
                names.update(df["gene_name"].dropna().astype(str).unique())
        except Exception:
            return set()
        return names

    reader = pd.read_csv(
        path,
        chunksize=100_000,
        low_memory=False,
        encoding="latin-1",
        encoding_errors="replace",
        on_bad_lines="skip",
    )
    for i, chunk in enumerate(reader):
        if "gene_name" not in chunk.columns:
            return set()
        names.update(chunk["gene_name"].dropna().astype(str).unique())
        if i + 1 >= max_chunks:
            break
    return names


def _rna_row_symbols(path: Path, gene_col: str, max_rows: int = 200_000) -> Set[str]:
    if not path.is_file():
        return set()
    try:
        hdr = pd.read_csv(
            path, sep="\t", nrows=0, encoding="latin-1", encoding_errors="replace"
        )
        use = gene_col if gene_col in hdr.columns else list(hdr.columns)[0]
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=[use],
            nrows=max_rows,
            low_memory=False,
            encoding="latin-1",
            encoding_errors="replace",
            on_bad_lines="skip",
        )
        return set(df[use].dropna().astype(str).unique())
    except Exception:
        return set()


def _methylation_probe_gene_tokens(path: Path, nrows: int = 80_000) -> Set[str]:
    if not path.is_file():
        return set()
    from pipeline.Methylation.probe_loader import _parse_gene_list

    df = pd.read_csv(
        path, sep="\t", nrows=nrows, low_memory=False, encoding="latin-1", encoding_errors="replace"
    )
    if "geneNames" not in df.columns:
        return set()
    out: Set[str] = set()
    for val in df["geneNames"].dropna():
        for g in _parse_gene_list(val):
            out.add(g)
    return out


def _rppa_genes(path: Path) -> Set[str]:
    if not path.is_file():
        return set()
    df = pd.read_csv(
        path, sep="\t", low_memory=False, encoding="latin-1", encoding_errors="replace"
    )
    if "gene_name" not in df.columns:
        return set()
    return set(df["gene_name"].dropna().astype(str).unique())


def _mirtar_target_genes(path: Path, nrows: int = 200_000) -> Set[str]:
    if not path.is_file():
        return set()
    df = pd.read_csv(
        path,
        nrows=nrows,
        low_memory=False,
        encoding="latin-1",
        encoding_errors="replace",
    )
    col = None
    for c in df.columns:
        if c.lower() in ("target gene", "genes", "gene", "targets"):
            col = c
            break
    if col is None:
        for c in df.columns:
            if "gene" in c.lower():
                col = c
                break
    if col is None:
        return set()
    return set(df[col].dropna().astype(str).str.strip().unique())


def _cnv_gene_file(path: Path) -> Set[str]:
    if not path.is_file():
        return set()
    df = pd.read_csv(
        path, nrows=50_000, low_memory=False, encoding="latin-1", encoding_errors="replace"
    )
    for c in ("gene_name", "Gene", "symbol", "SYMBOL"):
        if c in df.columns:
            return set(df[c].dropna().astype(str).unique())
    return set()


def _tier_label(gene: str, tiers: List[Tuple[str, List[str]]]) -> str:
    for label, genes in tiers:
        if gene in genes:
            return label
    return "?"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        type=Path,
        default=_REPO / "analysis" / "output" / "gene_symbol_dataset_coverage.md",
    )
    args = ap.parse_args()

    from pipeline.config import (
        PATHS,
        PRIMARY_GENES,
        EXTENDED_PRIMARY_GENES,
        TIER2_MEDIUM_GENES,
        TIER3_CNV_ONLY_GENES,
        TIER4_READOUT_GENES,
        PANEL_ALIAS_SEED_SYMBOLS,
    )
    from pipeline.genes.symbol_normalization import get_panel_symbol_mapping

    tiers: List[Tuple[str, List[str]]] = [
        ("PRIMARY", list(PRIMARY_GENES)),
        ("EXTENDED_T1", list(EXTENDED_PRIMARY_GENES)),
        ("TIER2", list(TIER2_MEDIUM_GENES)),
        ("TIER3", list(TIER3_CNV_ONLY_GENES)),
        ("TIER4", list(TIER4_READOUT_GENES)),
    ]

    gencode_pq = Path(PATHS.gencode_gtf_pq)
    gencode_csv = gencode_pq.parent / "gencode.v49.annotation.gtf.csv"
    gencode_path = gencode_pq if gencode_pq.is_file() else gencode_csv

    datasets = {
        "gencode": _gencode_gene_names(gencode_path),
        "rna_matrix": _rna_row_symbols(Path(PATHS.rna_expression), PATHS.rna_gene_col)
        | _rna_row_symbols(Path(PATHS.rna_expression_raw), PATHS.rna_gene_col),
        "meth_probes": _methylation_probe_gene_tokens(Path(PATHS.methylation_probe_reference)),
        "rppa_antibodies": _rppa_genes(Path(PATHS.rppa_antibody_annotation_csv)),
        "mirtar_targets": _mirtar_target_genes(Path(PATHS.mirtarbase_csv)),
        "cnv_gene_table": _cnv_gene_file(Path(PATHS.cnv_genes)),
    }

    mapping = get_panel_symbol_mapping(None)
    reverse = defaultdict(list)
    for old, new in mapping.items():
        reverse[new].append(old)

    all_panel = list(PANEL_ALIAS_SEED_SYMBOLS)

    lines: List[str] = []
    lines.append("# Gene symbol coverage vs datasets\n")
    lines.append(f"- GENCODE source: `{gencode_path}`\n")
    lines.append(f"- RNA: `{PATHS.rna_expression}` / raw\n")
    lines.append(f"- Methylation probe ref (sampled rows): `{PATHS.methylation_probe_reference}`\n")
    lines.append(f"- RPPA antibodies: `{PATHS.rppa_antibody_annotation_csv}`\n")
    lines.append(f"- miRTarBase (sampled rows): `{PATHS.mirtarbase_csv}`\n")
    lines.append(f"- CNV gene table: `{PATHS.cnv_genes}`\n")
    lines.append(f"- HGNC-derived alias entries (panel-seeded): **{len(mapping)}**\n")

    lines.append("\n## Summary counts (genes in dataset symbol sets)\n")
    lines.append("| Tier | n | in GENCODE | in RNA | in meth probes | in RPPA | in miRTar | in CNV table |\n")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|\n")

    for label, genes in tiers:
        gset = set(genes)
        lines.append(
            f"| {label} | {len(gset)} | "
            f"{sum(1 for g in gset if g in datasets['gencode'])} | "
            f"{sum(1 for g in gset if g in datasets['rna_matrix'])} | "
            f"{sum(1 for g in gset if g in datasets['meth_probes'])} | "
            f"{sum(1 for g in gset if g in datasets['rppa_antibodies'])} | "
            f"{sum(1 for g in gset if g in datasets['mirtar_targets'])} | "
            f"{sum(1 for g in gset if g in datasets['cnv_gene_table'])} |\n"
        )

    lines.append("\n## Per-gene detail (canonical symbol in config)\n")
    lines.append(
        "| Gene | Tier | GENCODE | RNA | Meth | RPPA | miRTar | CNVtbl | "
        "aliases→gene hitting RNA |\n|---|---|---|:---|:---|:---|:---|:---|:---|\n"
    )

    rna = datasets["rna_matrix"]
    for gene in sorted(set(all_panel)):
        tier = _tier_label(gene, tiers)
        in_g = "✓" if gene in datasets["gencode"] else ""
        in_r = "✓" if gene in rna else ""
        in_m = "✓" if gene in datasets["meth_probes"] else ""
        in_rp = "✓" if gene in datasets["rppa_antibodies"] else ""
        in_mt = "✓" if gene in datasets["mirtar_targets"] else ""
        in_cnv = "✓" if gene in datasets["cnv_gene_table"] else ""

        alias_hits = []
        if not in_r and gene in reverse:
            for alt in sorted(set(reverse[gene])):
                if alt in rna:
                    alias_hits.append(alt)
        alias_cell = ", ".join(alias_hits[:6]) + ("…" if len(alias_hits) > 6 else "")

        lines.append(
            f"| {gene} | {tier} | {in_g} | {in_r} | {in_m} | {in_rp} | {in_mt} | {in_cnv} | {alias_cell} |\n"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
