#!/usr/bin/env python3
"""
Panel coverage vs major **inputs**: how many core genes, pipeline genes, tier-1 lncRNAs,
and tier miRNAs appear; list symbols **missing** per source.

Outputs (default):
  ``analysis/output/panel_input_coverage_report.md``
  ``analysis/output/panel_input_coverage_missing_by_source.tsv``

Run from repo root (uses ``pipeline.config.PATHS``):

  python scripts/coverage/panel_input_coverage_report.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Set, Tuple

import pandas as pd

SourceKind = Literal["gene", "mirna"]

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

def _fmt_missing(xs: Iterable[str], max_list: int = 120) -> str:
    s = sorted(set(xs))
    if len(s) <= max_list:
        return ", ".join(s) if s else "—"
    return ", ".join(s[:max_list]) + f", … (+{len(s) - max_list} more)"


def _scan_set(
    wanted: Set[str],
    found: Set[str],
) -> Tuple[int, int, List[str]]:
    w = set(wanted)
    f = w & set(found)
    miss = sorted(w - f)
    return len(f), len(miss), miss


def _rna_symbols(path: Path, gene_col: str, nrows: Optional[int], flat: Dict[str, str]) -> Set[str]:
    from pipeline.genes.symbol_normalization import apply_symbol_mapping_series

    hdr = pd.read_csv(path, sep="\t", nrows=0, encoding="latin-1", encoding_errors="replace")
    col = gene_col if gene_col in hdr.columns else list(hdr.columns)[0]
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=[col],
        nrows=nrows,
        encoding="latin-1",
        encoding_errors="replace",
        on_bad_lines="skip",
    )
    s = apply_symbol_mapping_series(df[col].astype(str), flat)
    return set(s.unique())


def _mirna_arm_symbols(path: Path, nrows: Optional[int]) -> Set[str]:
    from pipeline.genes.panel_alias_registry import normalize_mirna_symbol

    mdf = pd.read_csv(path, sep="\t", nrows=nrows, encoding="latin-1", encoding_errors="replace")
    idx = mdf.iloc[:, 0].astype(str)
    return set(idx.map(normalize_mirna_symbol))


def _genes_only_symbols(path: Path, nrows: Optional[int], flat: Dict[str, str]) -> Set[str]:
    from pipeline.genes.symbol_normalization import apply_symbol_mapping_series

    df = pd.read_csv(path, nrows=nrows, low_memory=False)
    if "gene_name" not in df.columns:
        return set()
    g = apply_symbol_mapping_series(df["gene_name"].astype(str), flat)
    return set(g.unique())


def _methylation_gene_union(path: Path, flat: Dict[str, str]) -> Set[str]:
    from pipeline.Methylation.probe_loader import load_probe_reference

    sub = load_probe_reference(path)
    out: Set[str] = set()
    for lst in sub.get("gene_list", []):
        if not isinstance(lst, list):
            continue
        for t in lst:
            tok = str(t).strip()
            if tok:
                out.add(flat.get(tok, tok))
    return out


def _targetscan_genes(path: Path, max_rows: int, flat: Dict[str, str]) -> Set[str]:
    from pipeline.genes.symbol_normalization import apply_symbol_mapping_series

    df = pd.read_csv(
        path,
        sep="\t",
        nrows=max_rows,
        encoding="latin-1",
        encoding_errors="replace",
    )
    if "Gene Symbol" not in df.columns:
        return set()
    s = apply_symbol_mapping_series(df["Gene Symbol"].astype(str), flat)
    return set(s.unique())


def _mirtar_target_genes(path: Path, max_rows: int, flat: Dict[str, str]) -> Set[str]:
    from pipeline.genes.symbol_normalization import apply_symbol_mapping_series

    df = pd.read_csv(path, nrows=max_rows, encoding="latin-1", encoding_errors="replace")
    col = "Target Gene" if "Target Gene" in df.columns else None
    if not col:
        return set()
    s = apply_symbol_mapping_series(df[col].astype(str).str.strip(), flat)
    return set(s.unique())


def _rppa_genes(path: Path, flat: Dict[str, str]) -> Set[str]:
    from pipeline.genes.symbol_normalization import apply_symbol_mapping_series

    df = pd.read_csv(path, sep="\t", encoding="latin-1", encoding_errors="replace")
    if "gene_name" not in df.columns:
        return set()
    s = apply_symbol_mapping_series(df["gene_name"].astype(str), flat)
    return set(s.unique())


def _atac_linked_genes(path: Path, max_rows: int, flat: Dict[str, str]) -> Set[str]:
    """Parse ``linked_genes`` from ATAC peak table (list, JSON list, or comma-separated)."""
    if not path.is_file():
        return set()
    from pipeline.atac_peaks import load_atac_peaks_annotated

    if path.suffix.lower() == ".parquet":
        df = load_atac_peaks_annotated(path, columns=["linked_genes"])
    else:
        df = pd.read_csv(path, usecols=["linked_genes"], nrows=max_rows, low_memory=False)
    df = df.head(max_rows)
    out: Set[str] = set()
    for cell in df["linked_genes"]:
        genes: List[str] = []
        if isinstance(cell, list):
            genes = [str(x) for x in cell]
        elif isinstance(cell, str) and cell.strip().startswith("["):
            try:
                genes = [str(x) for x in json.loads(cell)]
            except Exception:
                genes = []
        elif isinstance(cell, str) and cell.strip() and cell.lower() not in ("nan", "none"):
            genes = [p.strip() for p in cell.replace(";", ",").split(",") if p.strip()]
        for tok in genes:
            out.add(flat.get(tok, tok))
    return out


def _count_data_lines(path: Path) -> int:
    n = 0
    with path.open("rb") as fh:
        for _ in fh:
            n += 1
    return max(0, n - 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-md",
        type=Path,
        default=_REPO / "analysis" / "output" / "panel_input_coverage_report.md",
    )
    ap.add_argument(
        "--out-tsv",
        type=Path,
        default=_REPO / "analysis" / "output" / "panel_input_coverage_missing_by_source.tsv",
    )
    ap.add_argument("--rna-nrows", type=int, default=60_000, help="RNA matrix rows to scan (None=all)")
    ap.add_argument("--mirna-nrows", type=int, default=20_000, help="miRNA matrix rows to scan")
    ap.add_argument("--genes-only-nrows", type=int, default=50_000)
    ap.add_argument("--targetscan-rows", type=int, default=400_000)
    ap.add_argument("--mirtar-rows", type=int, default=500_000)
    ap.add_argument("--atac-rows", type=int, default=200_000)
    args = ap.parse_args()

    from pipeline.config import (
        OUTPUT_SUBDIRS,
        PATHS,
        PRIMARY_GENES,
        PIPELINE_GENE_PANEL,
        PANEL_MIRNA_TIER_ARM_IDS,
        TIER1_LNCRNA_GENES,
    )
    from pipeline.genes.panel_alias_registry import get_gene_flat_alias_map

    flat = get_gene_flat_alias_map(None)

    primary = set(PRIMARY_GENES)
    panel = set(PIPELINE_GENE_PANEL)
    lncs = set(TIER1_LNCRNA_GENES)
    tier_mirs = set(PANEL_MIRNA_TIER_ARM_IDS)

    sources: List[Tuple[str, Set[str], str, SourceKind]] = []

    rna_p = Path(PATHS.rna_expression) if Path(PATHS.rna_expression).is_file() else Path(
        PATHS.rna_expression_raw
    )
    if rna_p.is_file():
        sources.append(
            (
                "RNA TPM matrix (row `gene_symbol`, canonicalized)",
                _rna_symbols(rna_p, PATHS.rna_gene_col, args.rna_nrows, flat),
                str(rna_p),
                "gene",
            )
        )

    mir_p = Path(PATHS.mirna_expression_tsv)
    if mir_p.is_file():
        sources.append(
            (
                "miRNA Xena matrix (row ids → mature arm via registry)",
                _mirna_arm_symbols(mir_p, args.mirna_nrows),
                str(mir_p),
                "mirna",
            )
        )

    go = Path(PATHS.genes_only)
    if go.is_file():
        sources.append(
            (
                "genes_only annotation (`gene_name`)",
                _genes_only_symbols(go, args.genes_only_nrows, flat),
                str(go),
                "gene",
            )
        )

    rppa = Path(PATHS.rppa_antibody_annotation_csv)
    if rppa.is_file():
        sources.append(("RPPA antibody table (`gene_name`)", _rppa_genes(rppa, flat), str(rppa), "gene"))

    tsp = Path(PATHS.targetscan_predictions)
    if tsp.is_file():
        sources.append(
            (
                f"TargetScan (`Gene Symbol`, first {args.targetscan_rows:,} rows)",
                _targetscan_genes(tsp, args.targetscan_rows, flat),
                str(tsp),
                "gene",
            )
        )

    mt = Path(PATHS.mirtarbase_csv)
    if mt.is_file():
        sources.append(
            (
                f"miRTarBase (`Target Gene`, first {args.mirtar_rows:,} rows)",
                _mirtar_target_genes(mt, args.mirtar_rows, flat),
                str(mt),
                "gene",
            )
        )

    pref = Path(PATHS.methylation_probe_reference)
    if pref.is_file():
        sources.append(
            (
                "Methylation probe reference (union of `gene_list` tokens)",
                _methylation_gene_union(pref, flat),
                str(pref),
                "gene",
            )
        )

    atac_dir = Path(PATHS.working_dir) / OUTPUT_SUBDIRS.get("atac_peaks", "atac_peaks")
    atac_tbl = atac_dir / "atac_peaks_annotated.parquet"
    if not atac_tbl.is_file():
        atac_tbl = atac_dir / "atac_peaks_annotated.csv"
    if atac_tbl.is_file():
        sources.append(
            (
                f"ATAC peak table `linked_genes` (first {args.atac_rows:,} rows)",
                _atac_linked_genes(atac_tbl, args.atac_rows, flat),
                str(atac_tbl),
                "gene",
            )
        )

    # CNV / SNV miRNA **locus** tables (interval overlap; not miRTarBase motifs).
    mir_mat = Path(PATHS.mirna_mature_loci_csv)
    mir_pre = Path(PATHS.mirna_precursor_loci_csv)
    cnv_note_lines: List[str] = []
    if mir_mat.is_file():
        n_m = _count_data_lines(mir_mat)
        cnv_note_lines.append(
            f"- **CNV/SNV overlap miRNA table** (default in code: `PATHS.mirna_mature_loci_csv`): "
            f"**{n_m}** data rows in `{mir_mat.name}` (full table is loaded for CNV/SNV; not reduced by symbol normalization)."
        )
    if mir_pre.is_file():
        n_p = _count_data_lines(mir_pre)
        cnv_note_lines.append(
            f"- Precursor hairpin table `cnv_miRNA.csv`: **{n_p}** rows (`PATHS.mirna_precursor_loci_csv`). "
            "CNV runner prefers **mature loci** when that path exists, then falls back to `mirna_path`."
        )

    md: List[str] = []
    md.append("# Panel vs input data — coverage summary\n\n")
    md.append("## Gene / lncRNA / miRNA sets\n\n")
    md.append(f"- **PRIMARY_GENES**: {len(primary)}\n")
    md.append(f"- **PIPELINE_GENE_PANEL**: {len(panel)}\n")
    md.append(f"- **Tier-1 lncRNAs** (extended panel): {len(lncs)}\n")
    md.append(f"- **PANEL_MIRNA_TIER_ARM_IDS**: {len(tier_mirs)}\n\n")

    if cnv_note_lines:
        md.append("## CNV / SNV miRNA locus coverage\n\n")
        md.extend(f"{ln}\n" for ln in cnv_note_lines)
        md.append(
            "\nSymbol normalization is **not** applied to miRNA **interval** tables (would corrupt MIMAT-style labels). "
            "Gene / lncRNA annotation tables are normalized.\n\n"
        )

    md.append("## Quantitative coverage — **gene-level** sources\n\n")
    md.append(
        "| Source | n PRIMARY hit | n PRIMARY miss | n PANEL hit | n PANEL miss | "
        "n lncRNA hit | n lncRNA miss |\n"
        "|---|---:|---:|---:|---:|---:|---:|\n"
    )

    tsv_rows: List[str] = []

    for label, found, path, kind in sources:
        if kind != "gene":
            continue
        p_hit, p_miss_n, p_miss = _scan_set(primary, found)
        g_hit, g_miss_n, g_miss = _scan_set(panel, found)
        l_hit, l_miss_n, l_miss = _scan_set(lncs, found)
        md.append(
            f"| {label} | {p_hit} | {p_miss_n} | {g_hit} | {g_miss_n} | {l_hit} | {l_miss_n} |\n"
        )
        tsv_rows.append(
            "\t".join([label[:120], str(path), "gene", "PRIMARY_MISS", _fmt_missing(p_miss, 200)])
        )
        tsv_rows.append("\t".join([label[:120], str(path), "gene", "PANEL_MISS", _fmt_missing(g_miss, 200)]))
        tsv_rows.append("\t".join([label[:120], str(path), "gene", "LNCRNA_MISS", _fmt_missing(l_miss, 200)]))

    md.append("\n## Quantitative coverage — **miRNA arm** sources\n\n")
    md.append("| Source | n tier-miRNA hit | n tier-miRNA miss |\n|---|---:|---:|\n")
    for label, found, path, kind in sources:
        if kind != "mirna":
            continue
        m_hit, m_miss_n, m_miss = _scan_set(tier_mirs, found)
        md.append(f"| {label} | {m_hit} | {m_miss_n} |\n")
        tsv_rows.append(
            "\t".join([label[:120], str(path), "mirna", "TIER_MIRNA_MISS", _fmt_missing(m_miss, 200)])
        )

    md.append("\n## Missing symbols — gene-level sources\n\n")
    for label, found, path, kind in sources:
        if kind != "gene":
            continue
        _, _, p_miss = _scan_set(primary, found)
        _, _, g_miss = _scan_set(panel, found)
        _, _, l_miss = _scan_set(lncs, found)
        md.append(f"### {label}\n\n")
        md.append(f"- Path: `{path}`\n")
        md.append(f"- **PRIMARY missing** ({len(p_miss)}): {_fmt_missing(p_miss)}\n")
        md.append(f"- **PANEL missing** ({len(g_miss)}): {_fmt_missing(g_miss)}\n")
        md.append(f"- **lncRNA missing** ({len(l_miss)}): {_fmt_missing(l_miss)}\n\n")

    md.append("## Missing symbols — miRNA arm sources\n\n")
    for label, found, path, kind in sources:
        if kind != "mirna":
            continue
        _, _, m_miss = _scan_set(tier_mirs, found)
        md.append(f"### {label}\n\n")
        md.append(f"- Path: `{path}`\n")
        md.append(f"- **Tier miRNA missing** ({len(m_miss)}): {_fmt_missing(m_miss)}\n\n")

    md.append("## SV FIMO motifs vs miRTarBase\n\n")
    md.append(
        "- **FIMO / MEME** in the SV pipeline uses ``PATHS.sv_meme_file`` (curated motif library), "
        "not motifs exported from miRTarBase.\n"
        "- **miRTarBase** is used for **gene–miRNA target evidence** (e.g. ``get_mirtarbase_targets`` / methylation / SNV miRNA **interval** overlap), "
        "not for building the SV MEME file.\n"
        "- **SNV** miRNA overlap uses ``mirna_mature_loci_csv`` (or fallback ``mirna_path``), same geometry idea as CNV.\n\n"
    )

    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text("".join(md), encoding="utf-8")

    header = "source_label\tpath\tkind\twhich\tmissing_symbols\n"
    args.out_tsv.write_text(header + "\n".join(tsv_rows) + "\n", encoding="utf-8")
    print(f"Wrote {args.out_md}")
    print(f"Wrote {args.out_tsv}")


if __name__ == "__main__":
    main()
