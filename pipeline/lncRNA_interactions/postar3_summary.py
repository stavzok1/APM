from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from pipeline.config import PATHS, PRIMARY_GENES, TIER1_LNCRNA_GENES
from pipeline.lncRNA_interactions.encori import build_encori_lncrna_target_list


@dataclass(frozen=True)
class Postar3SummaryPaths:
    postar3_parquet: Path = Path(PATHS.working_dir) / "RBP-RNA" / "POSTAR3.parquet"
    lncrna_pairs_csv: Path = Path(PATHS.working_dir) / "lncRNA_matching" / "genes_lncRNAs_1000000bp_distances.csv"
    out_overlap_parquet: Path = Path(PATHS.working_dir) / "RBP-RNA" / "POSTAR3.overlaps_selected_lncrnas.parquet"
    out_rbp_summary_parquet: Path = Path(PATHS.working_dir) / "RBP-RNA" / "POSTAR3.selected_lncrnas.rbp_summary.parquet"
    out_lncrna_summary_parquet: Path = Path(PATHS.working_dir) / "RBP-RNA" / "POSTAR3.selected_lncrnas.lncrna_summary.parquet"


Postar3LncRegionMode = Literal["gene", "exons", "introns", "promoter"]


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    iv = sorted(intervals, key=lambda x: (x[0], x[1]))
    out: list[tuple[int, int]] = []
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out


def _subtract_intervals(span: tuple[int, int], sub: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """
    Return span minus union(sub). Assumes sub is merged/sorted.
    """
    s0, e0 = span
    if not sub:
        return [(s0, e0)]
    out: list[tuple[int, int]] = []
    cur = s0
    for s, e in sub:
        if e < cur:
            continue
        if s > e0:
            break
        if s > cur:
            out.append((cur, min(e0, s - 1)))
        cur = max(cur, e + 1)
        if cur > e0:
            break
    if cur <= e0:
        out.append((cur, e0))
    return [(s, e) for s, e in out if s <= e]


def _load_lncrna_intervals_from_gencode(
    lncrnas: Iterable[str],
    gencode_parquet: Path,
) -> pd.DataFrame:
    """
    Load gene-level intervals for the requested lncRNA gene symbols.

    Uses the full gencode parquet (preferred) and falls back to slim if needed.
    """
    wanted = {str(x).strip() for x in lncrnas if str(x).strip()}
    if not wanted:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "gene_name", "gene_id", "gene_type"])

    d = ds.dataset(str(gencode_parquet), format="parquet")
    # Columns are either (chrom,start,end) or (seqname,start,end) depending on source; normalize to chrom.
    cols = set(d.schema.names)
    chrom_col = "chrom" if "chrom" in cols else "seqname" if "seqname" in cols else "chromosome" if "chromosome" in cols else None
    if chrom_col is None:
        raise ValueError(f"Could not find chrom column in {gencode_parquet}; columns={sorted(cols)[:30]}")

    # Pull only gene rows; gene_type should contain "lncRNA" in this repo's derived tables.
    base_cols = [chrom_col, "start", "end", "strand", "gene_id", "gene_type", "feature"]
    if "gene_name" in cols:
        base_cols.insert(4, "gene_name")
    elif "gene" in cols:
        base_cols.insert(4, "gene")
    else:
        raise ValueError(f"Missing gene symbol column (expected gene_name/gene) in {gencode_parquet}")

    table = d.to_table(columns=base_cols, filter=(ds.field("feature") == "gene"))
    df = table.to_pandas()
    # Some GENCODE parquet builds store HGNC symbol under `gene_name` OR `gene`.
    if "gene_name" in df.columns:
        sym = df["gene_name"].astype(str)
    elif "gene" in df.columns:
        sym = df["gene"].astype(str)
    else:
        raise ValueError(f"Missing gene symbol column (expected gene_name/gene) in {gencode_parquet}")
    df = df[sym.isin(wanted)].copy()
    df = df.rename(columns={chrom_col: "chrom"})
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    if "gene_name" not in df.columns and "gene" in df.columns:
        df["gene_name"] = df["gene"].astype(str)
    return df[["chrom", "start", "end", "strand", "gene_name", "gene_id", "gene_type"]].drop_duplicates()


def _load_lncrna_feature_intervals_from_gencode(
    lncrnas: Iterable[str],
    gencode_parquet: Path,
    *,
    mode: Postar3LncRegionMode,
    promoter_upstream_bp: int = 2000,
    promoter_downstream_bp: int = 500,
) -> pd.DataFrame:
    """
    Build intervals for feature-aware overlap:
    - gene: gene span
    - exons: merged exon intervals per lncRNA gene
    - introns: (gene span) minus (merged exon intervals)
    - promoter: TSS-centered window using strand
    """
    gene_iv = _load_lncrna_intervals_from_gencode(lncrnas, gencode_parquet)
    if gene_iv.empty:
        return gene_iv
    if mode == "gene":
        gene_iv["region_mode"] = "gene"
        return gene_iv

    d = ds.dataset(str(gencode_parquet), format="parquet")
    cols = set(d.schema.names)
    chrom_col = "chrom" if "chrom" in cols else "seqname" if "seqname" in cols else "chromosome" if "chromosome" in cols else None
    if chrom_col is None:
        raise ValueError(f"Could not find chrom column in {gencode_parquet}; columns={sorted(cols)[:30]}")

    wanted = set(gene_iv["gene_name"].astype(str).tolist())

    if mode in ("exons", "introns"):
        # Load exon rows for these gene symbols, then merge to per-gene exon unions.
        table = d.to_table(
            columns=[chrom_col, "start", "end", "strand", "gene_name", "gene_id", "gene_type", "feature"],
            filter=(ds.field("feature") == "exon"),
        )
        ex = table.to_pandas()
        ex = ex[ex["gene_name"].astype(str).isin(wanted)].copy()
        ex = ex.rename(columns={chrom_col: "chrom"})
        ex["start"] = pd.to_numeric(ex["start"], errors="coerce").astype("Int64")
        ex["end"] = pd.to_numeric(ex["end"], errors="coerce").astype("Int64")
        ex = ex.dropna(subset=["start", "end"])
        ex["start"] = ex["start"].astype(int)
        ex["end"] = ex["end"].astype(int)

        ex_by_gene: dict[str, dict[str, object]] = {}
        for g, sub in ex.groupby("gene_name", sort=False):
            chrom = str(sub["chrom"].iloc[0])
            strand = str(sub["strand"].iloc[0]) if "strand" in sub.columns else "+"
            gene_id = str(sub["gene_id"].iloc[0]) if "gene_id" in sub.columns else ""
            gene_type = str(sub["gene_type"].iloc[0]) if "gene_type" in sub.columns else "lncRNA"
            iv = _merge_intervals([(int(r.start), int(r.end)) for r in sub.itertuples(index=False)])
            ex_by_gene[str(g)] = {"chrom": chrom, "strand": strand, "gene_id": gene_id, "gene_type": gene_type, "exons": iv}

        rows: list[dict[str, object]] = []
        for r in gene_iv.itertuples(index=False):
            g = str(r.gene_name)
            chrom = str(r.chrom)
            strand = str(r.strand)
            gene_id = str(r.gene_id)
            gene_type = str(r.gene_type)
            span = (int(r.start), int(r.end))
            exons = ex_by_gene.get(g, {}).get("exons", [])
            if mode == "exons":
                for s, e in exons:
                    rows.append({"chrom": chrom, "start": s, "end": e, "strand": strand, "gene_name": g, "gene_id": gene_id, "gene_type": gene_type, "region_mode": "exons"})
            else:
                intr = _subtract_intervals(span, list(exons))
                for s, e in intr:
                    rows.append({"chrom": chrom, "start": s, "end": e, "strand": strand, "gene_name": g, "gene_id": gene_id, "gene_type": gene_type, "region_mode": "introns"})
        return pd.DataFrame(rows)

    if mode == "promoter":
        rows = []
        for r in gene_iv.itertuples(index=False):
            tss = int(r.start) if str(r.strand) == "+" else int(r.end)
            if str(r.strand) == "+":
                s = max(0, tss - int(promoter_upstream_bp))
                e = tss + int(promoter_downstream_bp)
            else:
                s = max(0, tss - int(promoter_downstream_bp))
                e = tss + int(promoter_upstream_bp)
            rows.append({"chrom": str(r.chrom), "start": int(s), "end": int(e), "strand": str(r.strand), "gene_name": str(r.gene_name), "gene_id": str(r.gene_id), "gene_type": str(r.gene_type), "region_mode": "promoter"})
        return pd.DataFrame(rows)

    raise ValueError(f"Unknown mode: {mode}")


def _overlap_mask_for_intervals(starts: np.ndarray, ends: np.ndarray, intervals: list[tuple[int, int]]) -> np.ndarray:
    """
    Boolean mask of rows overlapping ANY of the given intervals.
    """
    if not intervals:
        return np.zeros_like(starts, dtype=bool)
    m = np.zeros_like(starts, dtype=bool)
    for s, e in intervals:
        # overlap if start <= e and end >= s
        m |= (starts <= e) & (ends >= s)
    return m


def build_postar3_overlaps_for_selected_lncrnas(
    *,
    n_extra_close_lncrnas: int = 20,
    paths: Postar3SummaryPaths = Postar3SummaryPaths(),
    gencode_parquet: Path | None = None,
    chunk_rows: int = 1_500_000,
    region_mode: Postar3LncRegionMode = "gene",
    promoter_upstream_bp: int = 2000,
    promoter_downstream_bp: int = 500,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Stream scan POSTAR3.parquet and keep only peaks overlapping selected lncRNA gene loci.
    """
    gencode_parquet = gencode_parquet or Path(PATHS.gencode_gtf_full_pq)

    lncrnas = build_encori_lncrna_target_list(
        panel_lncrnas=TIER1_LNCRNA_GENES,
        lncrna_proximity_pairs_csv=paths.lncrna_pairs_csv,
        primary_genes=PRIMARY_GENES,
        n_extra_close_lncrnas=int(n_extra_close_lncrnas),
    )

    lncrna_iv = _load_lncrna_feature_intervals_from_gencode(
        lncrnas,
        gencode_parquet,
        mode=region_mode,
        promoter_upstream_bp=int(promoter_upstream_bp),
        promoter_downstream_bp=int(promoter_downstream_bp),
    )
    if lncrna_iv.empty:
        raise ValueError("No lncRNA intervals loaded; cannot overlap POSTAR3.")

    # group intervals by chromosome for fast masking
    iv_by_chr: dict[str, list[tuple[int, int]]] = {}
    for chrom, sub in lncrna_iv.groupby("chrom", sort=False):
        iv_by_chr[str(chrom)] = [(int(r.start), int(r.end)) for r in sub.itertuples(index=False)]

    d = ds.dataset(str(paths.postar3_parquet), format="parquet")
    keep_cols = ["chrom", "start", "end", "strand", "rbp", "assay", "cell_tissue", "source_accession", "score"]

    out_rows: list[pd.DataFrame] = []
    for chrom in sorted(iv_by_chr.keys()):
        # scan per chromosome so we can do in-memory overlap masks cheaply
        scanner = d.scanner(columns=keep_cols, filter=(ds.field("chrom") == chrom), batch_size=int(chunk_rows))
        intervals = iv_by_chr[chrom]
        for batch in scanner.to_batches():
            df = batch.to_pandas()
            starts = df["start"].to_numpy(dtype=np.int64, copy=False)
            ends = df["end"].to_numpy(dtype=np.int64, copy=False)
            m = _overlap_mask_for_intervals(starts, ends, intervals)
            if not m.any():
                continue
            hit = df.loc[m].copy()
            hit["__selected_set"] = f"tier1_plus_{int(n_extra_close_lncrnas)}close"
            hit["__region_mode"] = str(region_mode)
            out_rows.append(hit)

    overlaps = pd.concat(out_rows, ignore_index=True) if out_rows else pd.DataFrame(columns=keep_cols + ["__selected_set"])
    paths.out_overlap_parquet.parent.mkdir(parents=True, exist_ok=True)
    overlaps.to_parquet(paths.out_overlap_parquet, index=False)

    # Per-RBP summary (across selected lncRNA loci)
    if overlaps.empty:
        rbp_summary = pd.DataFrame(columns=["rbp", "n_peaks", "n_assays", "n_cell_tissue", "assays", "cell_tissue"])
        lncrna_summary = pd.DataFrame(columns=["gene_name", "chrom", "start", "end", "n_overlapping_peaks"])
    else:
        rbp_summary = (
            overlaps.groupby("rbp", sort=False)
            .agg(
                n_peaks=("rbp", "size"),
                n_assays=("assay", lambda s: int(pd.Series(s).nunique(dropna=True))),
                n_cell_tissue=("cell_tissue", lambda s: int(pd.Series(s).nunique(dropna=True))),
                assays=("assay", lambda s: sorted(pd.unique(pd.Series(s).dropna().astype(str)))[:50]),
                cell_tissue=("cell_tissue", lambda s: sorted(pd.unique(pd.Series(s).dropna().astype(str)))[:50]),
            )
            .reset_index()
            .sort_values(["n_peaks", "rbp"], ascending=[False, True])
        )

        # Per-lncRNA locus summary (how many POSTAR peaks fall in each lncRNA gene interval)
        # Here we do a simple interval count by chrom again; n is small.
        counts = []
        by_chr = overlaps.groupby("chrom", sort=False)
        for chrom, ivs in lncrna_iv.groupby("chrom", sort=False):
            chrom = str(chrom)
            ov = by_chr.get_group(chrom) if chrom in by_chr.groups else None
            if ov is None or ov.empty:
                for r in ivs.itertuples(index=False):
                    counts.append({"gene_name": r.gene_name, "chrom": chrom, "start": int(r.start), "end": int(r.end), "n_overlapping_peaks": 0})
                continue
            ov_s = ov["start"].to_numpy(dtype=np.int64, copy=False)
            ov_e = ov["end"].to_numpy(dtype=np.int64, copy=False)
            for r in ivs.itertuples(index=False):
                m = (ov_s <= int(r.end)) & (ov_e >= int(r.start))
                counts.append({"gene_name": r.gene_name, "chrom": chrom, "start": int(r.start), "end": int(r.end), "n_overlapping_peaks": int(m.sum())})
        lncrna_summary = pd.DataFrame(counts).sort_values(["n_overlapping_peaks", "gene_name"], ascending=[False, True])

    rbp_summary.to_parquet(paths.out_rbp_summary_parquet, index=False)
    lncrna_summary.to_parquet(paths.out_lncrna_summary_parquet, index=False)

    return overlaps, rbp_summary, lncrna_summary

