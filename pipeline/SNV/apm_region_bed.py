from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import pandas as pd


def _genes_for_apm_region_bed() -> List[str]:
    """
    Genes to include in the SNV pre-VEP bedtools-intersect region BED.

    This BED controls which variants end up in ``*.APM_1Mb(.vep).vcf`` outputs.
    We include:
    - the main SNV panel genes (PIPELINE_GENE_PANEL; Tier-1 extended by default)
    - plus SNV/CNV stratifiers (TP53/PTEN/PIK3CA, etc.)
    """
    from pipeline.config import PIPELINE_GENE_PANEL, TIER_SNV_CNV_STRATIFIER_GENES

    genes = list(PIPELINE_GENE_PANEL) + list(TIER_SNV_CNV_STRATIFIER_GENES)
    # Stable, unique, uppercase while preserving original order is not essential for BED.
    out: List[str] = []
    seen = set()
    for g in genes:
        gU = str(g).strip().upper()
        if not gU or gU in seen:
            continue
        seen.add(gU)
        out.append(gU)
    return out


def build_gene_flank_bed(
    *,
    genes: Iterable[str],
    gencode_parquet: Path,
    flank_bp: int,
) -> pd.DataFrame:
    """
    Build a BED (0-based start, half-open end) containing +/- flank_bp around gene bodies.

    Uses only ``feature == 'gene'`` rows from the parquet to keep it light.
    """
    gencode_parquet = Path(gencode_parquet)
    if not gencode_parquet.exists():
        raise FileNotFoundError(f"Missing GENCODE parquet: {gencode_parquet}")

    genes = [str(g).strip().upper() for g in genes if str(g).strip()]
    if not genes:
        raise ValueError("No genes provided")

    # Fast parquet predicate pushdown (pyarrow engine).
    # Parquet schema uses ``seqname`` (GTF-like) rather than ``chrom``.
    cols = ["seqname", "start", "end", "feature", "gene_name", "strand"]
    df = pd.read_parquet(
        gencode_parquet,
        columns=cols,
        filters=[("feature", "==", "gene"), ("gene_name", "in", genes)],
    )
    if df.empty:
        raise ValueError("No gene rows found in parquet for requested gene set")

    # Coordinates in our tables are 1-based inclusive (GTF-like). Convert to BED:
    # start0 = start1 - 1, end_excl = end1 (since end1 inclusive -> end_excl = end1).
    df = df.copy()
    df["gene_name"] = df["gene_name"].astype(str).str.upper()
    df["chrom"] = df["seqname"].astype(str)
    df["start0"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64") - 1
    df["end_excl"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start0", "end_excl"])

    df["start0"] = (df["start0"] - int(flank_bp)).clip(lower=0)
    df["end_excl"] = df["end_excl"] + int(flank_bp)

    bed = df[["chrom", "start0", "end_excl", "gene_name", "strand"]].copy()
    bed = bed.drop_duplicates(subset=["chrom", "start0", "end_excl", "gene_name"])
    bed = bed.sort_values(["chrom", "start0", "end_excl", "gene_name"], kind="mergesort")
    return bed


def write_bed(df: pd.DataFrame, out_path: Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # BED columns: chrom, start, end, name, strand (tab-separated; no header)
    df.to_csv(out_path, sep="\t", header=False, index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build +/- flank BED for SNV pre-VEP APM-region subsetting.")
    ap.add_argument("--out", type=str, default="", help="Output BED path.")
    ap.add_argument("--flank-bp", type=int, default=1_000_000, help="+/- flank around gene body.")
    ap.add_argument("--gencode-parquet", type=str, default="", help="Override PATHS.gencode_gtf_pq.")
    args = ap.parse_args()

    from pipeline.config import PATHS

    genes = _genes_for_apm_region_bed()
    gtf_pq = Path(args.gencode_parquet) if args.gencode_parquet else Path(PATHS.gencode_gtf_pq)
    out = Path(args.out) if args.out else (Path(PATHS.working_dir) / "SNV" / "apm_genes_1Mb.bed")

    bed = build_gene_flank_bed(genes=genes, gencode_parquet=gtf_pq, flank_bp=int(args.flank_bp))
    write_bed(bed, out)
    print(f"Wrote BED: {out}  (genes={len(set(bed['gene_name']))}, rows={len(bed)})")


if __name__ == "__main__":
    main()

