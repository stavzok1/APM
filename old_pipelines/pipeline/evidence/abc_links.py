"""
ABC model enhancer-gene link processing.

Complete pipeline: load → filter → aggregate → map to cCREs → build nested dicts
"""

from pathlib import Path
from typing import List, Set, Dict, Any, Optional, Iterable

import numpy as np
import pandas as pd

from ..utils import safe_float, safe_int
from ..schemas import empty_abc_celltype_entry, empty_abc_block


# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================

COLS_ABC_RAW = [
    "chr", "start", "end", "name", "class", "activity_base",
    "TargetGene", "TargetGeneTSS", "TargetGeneExpression",
    "TargetGenePromoterActivityQuantile", "TargetGeneIsExpressed",
    "distance", "isSelfPromomer", "hic_contact", "powerlaw_contact",
    "powerlaw_contact_reference", "hic_contact_pl_scaled", "hic_pseudocount",
    "hic_contact_pl_scaled_adj", "ABC.Score.Numerator", "ABC.Score",
    "powerlaw.Score.Numerator", "powerlaw.Score", "CellType",
]

RENAME_ABC = {
    "name": "enh_id",
    "chr": "chrom",
    "class": "element_class",
    "activity_base": "activity",
    "TargetGene": "gene_name",
    "TargetGeneExpression": "gene_expr",
    "TargetGenePromoterActivityQuantile": "promoter_activity_q",
    "TargetGeneIsExpressed": "gene_is_expressed",
    "isSelfPromoter": "is_self_promoter",
    "hic_contact_pl_scaled_adj": "hic_pl_scaled",
    "ABC.Score.Numerator": "ABC_num",
    "ABC.Score": "ABC_score",
    "powerlaw.Score": "powerlaw_score",
    "CellType": "cell_type",
}

ABC_NUMERIC_COLS = [
    "start", "end", "activity", "gene_expr", "promoter_activity_q",
    "distance", "hic_pl_scaled", "ABC_num", "ABC_score", "powerlaw_score"
]


# =============================================================================
# STREAMING & FILTERING
# =============================================================================

def _stream_and_filter_abc(
    file_path: Path,
    gene_set: Set[str],
    celltype_set: Optional[Set[str]],
    chunksize: int,
) -> pd.DataFrame:
    """Stream ABC file and retain rows matching genes and cell types."""
    chunks = []
    
    # Determine which columns we actually need
    keep_cols = [
        "chr", "start", "end", "name", "class", "activity_base", "TargetGene",
        "TargetGeneExpression", "TargetGenePromoterActivityQuantile", "TargetGeneIsExpressed",
        "distance", "isSelfPromoter", "hic_contact_pl_scaled_adj",
        "ABC.Score.Numerator", "ABC.Score", "powerlaw.Score", "CellType",
    ]
    
    # Read with header row
    reader = pd.read_csv(
        file_path, sep="\t", header=0,
        chunksize=chunksize, low_memory=False, compression="infer",
    )

    for i, chunk in enumerate(reader):
        # Standardize column names if needed
        if "chr" not in chunk.columns and "chrom" in chunk.columns:
            chunk = chunk.rename(columns={"chrom": "chr"})
        
        chunk["TargetGene"] = chunk["TargetGene"].astype(str).str.strip()
        
        if "CellType" in chunk.columns:
            chunk["CellType"] = chunk["CellType"].astype(str).str.strip()

        mask = chunk["TargetGene"].isin(gene_set)
        if celltype_set and "CellType" in chunk.columns:
            mask &= chunk["CellType"].isin(celltype_set)

        if mask.any():
            # Only keep columns we need
            cols_to_keep = [c for c in keep_cols if c in chunk.columns]
            chunks.append(chunk.loc[mask, cols_to_keep].copy())

        if (i + 1) % 20 == 0:
            print(f"  Processed {i + 1} chunks...")

    if not chunks:
        raise ValueError("No ABC rows found for selected genes/cell types.")

    df = pd.concat(chunks, ignore_index=True)
    print(f"  Total rows extracted: {len(df)}")
    return df


# =============================================================================
# AGGREGATION
# =============================================================================

def _aggregate_per_enhancer_gene_celltype(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate ABC data to one row per (enhancer, gene, cell_type)."""
    group_cols = ["chrom", "start", "end", "enh_id", "gene_name", "cell_type"]

    return df.groupby(group_cols).agg(
        element_class=("element_class", "first"),
        activity=("activity", "max"),
        distance=("distance", "median"),
        is_self_promoter=("is_self_promoter", "max"),
        hic_pl_scaled=("hic_pl_scaled", "max"),
        ABC_num=("ABC_num", "max"),
        ABC_score=("ABC_score", "max"),
        powerlaw_score=("powerlaw_score", "max"),
        gene_expr=("gene_expr", "first"),
        promoter_activity_q=("promoter_activity_q", "first"),
        gene_is_expressed=("gene_is_expressed", "max"),
    ).reset_index()


def _add_strength_flags(
    df_agg: pd.DataFrame,
    present_threshold: float,
    strong_threshold: float,
) -> pd.DataFrame:
    """Add presence/strength flags and rank within gene."""
    df_agg = df_agg.copy()
    
    # Presence: ABC score above threshold
    df_agg["is_present"] = df_agg["ABC_score"] >= present_threshold
    
    # Strong: ABC score above strong threshold AND not a promoter
    df_agg["is_strong"] = (
        (df_agg["ABC_score"] >= strong_threshold) & 
        (df_agg["element_class"] != "promoter")
    )
    
    # Rank within gene (normalized to max score per gene/celltype)
    max_score = df_agg.groupby(["gene_name", "cell_type"])["ABC_score"].transform("max")
    df_agg["rank_within_gene"] = np.where(max_score > 0, df_agg["ABC_score"] / max_score, 0.0)
    
    return df_agg


# =============================================================================
# DICT BUILDERS
# =============================================================================

def _build_abc_full_dict(
    df_agg: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse to ABC_full dict per (enhancer, gene)."""
    
    def build_celltype_dict(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        out = {}
        for _, r in df_sub.iterrows():
            out[r["cell_type"]] = {
                "ABC_score": safe_float(r["ABC_score"]),
                "ABC_num": safe_float(r["ABC_num"]),
                "activity": safe_float(r["activity"]),
                "distance": safe_int(r["distance"]),
                "element_class": r["element_class"],
                "is_self_promoter": bool(r["is_self_promoter"]),
                "hic_pl_scaled": safe_float(r["hic_pl_scaled"]),
                "powerlaw_score": safe_float(r["powerlaw_score"]),
                "gene_expr": safe_float(r["gene_expr"]),
                "promoter_activity_q": safe_float(r["promoter_activity_q"]),
                "gene_is_expressed": bool(r["gene_is_expressed"]),
                "rank_within_gene": safe_float(r["rank_within_gene"]),
                "is_present": bool(r["is_present"]),
                "is_strong": bool(r["is_strong"]),
            }
        return out

    return (
        df_agg.groupby(["chrom", "start", "end", "enh_id", "gene_name"], sort=False)
        .apply(build_celltype_dict)
        .rename("ABC_full")
        .reset_index()
    )


# =============================================================================
# PUBLIC API
# =============================================================================

def build_abc_links(
    abc_path: Path,
    gene_list: Iterable[str],
    celltypes: Optional[List[str]] = None,
    present_threshold: float = 0.015,
    strong_threshold: float = 0.05,
    chunksize: int = 500_000,
) -> pd.DataFrame:
    """
    Build enhancer-gene links from ABC model predictions.
    
    Returns DataFrame with columns:
    - chrom, start, end, enh_id, gene_name
    - ABC_full: dict of {celltype: {ABC_score, activity, ...}}
    """
    print("Building ABC enhancer links...")

    gene_set = set(gene_list)
    celltype_set = set(celltypes) if celltypes else None

    df = _stream_and_filter_abc(abc_path, gene_set, celltype_set, chunksize)
    
    # Rename columns
    rename_map = {k: v for k, v in RENAME_ABC.items() if k in df.columns}
    df = df.rename(columns=rename_map)

    # Convert numeric columns
    for col in ABC_NUMERIC_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Convert boolean columns
    if "gene_is_expressed" in df.columns:
        df["gene_is_expressed"] = df["gene_is_expressed"].astype(bool)
    if "is_self_promoter" in df.columns:
        df["is_self_promoter"] = df["is_self_promoter"].astype(bool)

    # Aggregate per enhancer-gene-celltype
    df_agg = _aggregate_per_enhancer_gene_celltype(df)
    
    # Add strength flags
    df_agg = _add_strength_flags(df_agg, present_threshold, strong_threshold)
    
    # Build ABC_full dict
    abc_full = _build_abc_full_dict(df_agg)

    print(f"  Built {len(abc_full)} enhancer-gene links")
    return abc_full


# =============================================================================
# MAP TO cCREs
# =============================================================================

def map_abc_to_ccres(
    df_abc: pd.DataFrame,
    ccres: pd.DataFrame,
) -> pd.DataFrame:
    """
    Map ABC enhancers to cCREs using center-in-interval rule.
    
    Args:
        df_abc: ABC links from build_abc_links
        ccres: cCRE DataFrame with chrom, start, end, ENCODE_id
    
    Returns:
        DataFrame with ENCODE_id added, matching ABC enhancers to cCREs
    """
    print("Mapping ABC enhancers to cCREs...")

    # Prepare cCREs
    reg = ccres[["chrom", "start", "end", "ENCODE_id"]].copy()
    reg["start"] = pd.to_numeric(reg["start"], errors="coerce")
    reg["end"] = pd.to_numeric(reg["end"], errors="coerce")
    reg = reg.dropna(subset=["start", "end"])
    reg = reg.rename(columns={"start": "reg_start", "end": "reg_end"})
    reg[["reg_start", "reg_end"]] = reg[["reg_start", "reg_end"]].astype("int64")

    # Prepare ABC enhancers
    abc = df_abc.copy()
    abc["start"] = pd.to_numeric(abc["start"], errors="coerce")
    abc["end"] = pd.to_numeric(abc["end"], errors="coerce")
    abc = abc.dropna(subset=["start", "end"])
    abc = abc.rename(columns={"start": "enh_start", "end": "enh_end"})
    abc[["enh_start", "enh_end"]] = abc[["enh_start", "enh_end"]].astype("int64")
    abc["enh_center"] = ((abc["enh_start"] + abc["enh_end"]) // 2).astype("int64")

    # Map per chromosome using merge_asof
    mapped_chunks = []
    for chrom, reg_chr in reg.groupby("chrom"):
        abc_chr = abc[abc["chrom"] == chrom]
        if abc_chr.empty:
            continue

        reg_chr = reg_chr.sort_values("reg_start").reset_index(drop=True)
        abc_chr = abc_chr.sort_values("enh_center").reset_index(drop=True)

        merged = pd.merge_asof(
            abc_chr, reg_chr,
            left_on="enh_center", right_on="reg_start",
            direction="backward",
        )

        # Keep only where center is inside the cCRE
        inside = (merged["enh_center"] >= merged["reg_start"]) & (merged["enh_center"] < merged["reg_end"])
        merged = merged[inside]

        if not merged.empty:
            merged["chrom"] = chrom
            keep_cols = ["ENCODE_id", "gene_name", "ABC_full", "chrom", "enh_start", "enh_end"]
            mapped_chunks.append(merged[[c for c in keep_cols if c in merged.columns]])

    if not mapped_chunks:
        print("  Warning: No ABC enhancers could be mapped to cCREs")
        return pd.DataFrame(columns=["ENCODE_id", "gene_name", "ABC_full"])

    result = pd.concat(mapped_chunks, ignore_index=True)
    result = result.drop_duplicates(
        subset=["ENCODE_id", "gene_name", "chrom", "enh_start", "enh_end"],
        keep="first",
    )

    print(f"  Mapped {len(result)} enhancer-cCRE pairs")
    return result


def collapse_abc_per_link(df_mapped: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse multiple ABC enhancers per (cCRE_id, gene_name) into a list.
    
    Returns DataFrame with:
    - ENCODE_id, gene_name
    - ABC_enhancers: list of {start, end, ABC_full}
    """
    def build_enhancer_list(df_sub: pd.DataFrame) -> List[Dict[str, Any]]:
        enhancers = []
        for _, r in df_sub.iterrows():
            enhancers.append({
                "start": int(r["enh_start"]) if pd.notna(r.get("enh_start")) else None,
                "end": int(r["enh_end"]) if pd.notna(r.get("enh_end")) else None,
                "ABC_full": r.get("ABC_full", {}),
            })
        return enhancers

    result = (
        df_mapped
        .groupby(["ENCODE_id", "gene_name"], as_index=False)
        .apply(lambda df_sub: pd.Series({"ABC_enhancers": build_enhancer_list(df_sub)}))
    )
    
    return result
