"""
Methylation probe reference loading and annotation.

Loads GDC methylation probe reference and annotates with:
- Promoter overlap for panel genes
- Gene body overlap
- cCRE overlap
- lncRNA promoter overlap
- ATAC peak overlap
- TAD domain context

Main functions:
    load_probe_reference: Load and parse GDC probe reference
    annotate_probes_with_promoters: Add gene promoter annotations
    annotate_probes_with_gene_bodies: Add gene body annotations
    annotate_probes_with_ccres: Add cCRE overlap annotations
    annotate_probes_with_lncrnas: Add lncRNA promoter annotations
    annotate_probes_with_atac: Add ATAC peak overlap annotations
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple

import pandas as pd
import numpy as np

from ..config import THRESHOLDS
from ..utils import harmonize_chrom_column
from .meth_schemas import METHYLATION_COLUMNS as MC, CGI_CONTEXTS


# =============================================================================
# PROBE REFERENCE LOADING
# =============================================================================

def load_probe_reference(
    path: Path,
    sep: str = "\t",
    gene_symbol_mapping: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """
    Load GDC methylation probe reference.
    
    Expected columns:
        CpG_chrm, CpG_beg, CpG_end, probe_strand, probeID,
        genesUniq, geneNames, transcriptTypes, transcriptIDs,
        distToTSS, CGI, CGIposition
    
    Args:
        path: Path to probe reference file
        sep: Column separator (default: tab)
        gene_symbol_mapping: Optional old_symbol -> canonical for ``gene_list`` tokens.
            When None and env ``APM_USE_GENE_SYMBOL_MAPPING`` is on (default), uses
            the pipeline-wide HGNC/UCSC/legacy map.

    Returns:
        DataFrame with standardized column names and parsed fields
    """
    print(f"Loading probe reference from {path}...")
    
    df = pd.read_csv(path, sep=sep, low_memory=False)
    
    # Standardize column names
    rename_map = {
        "CpG_chrm": "chrom",
        "CpG_beg": "start",
        "CpG_end": "end",
        "probe_strand": "strand",
    }
    df = df.rename(columns=rename_map)
    
    # Harmonize chromosomes (ensure chr prefix)
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].where(
        df["chrom"].str.startswith("chr"),
        "chr" + df["chrom"]
    )
    
    # Ensure numeric coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    
    # Parse CGI context
    df["in_CGI"] = df["CGI"].notna() & (df["CGI"] != "")
    df["CGI_context"] = df["CGIposition"].fillna("OpenSea").replace("", "OpenSea")
    
    # Parse gene lists (semicolon or comma separated in GDC format)
    df["gene_list"] = df["geneNames"].apply(_parse_gene_list)

    _map = gene_symbol_mapping
    if _map is None:
        if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() not in (
            "0",
            "false",
            "no",
        ):
            from ..genes.symbol_normalization import default_symbol_mapping

            _map = default_symbol_mapping()
    if _map:
        from ..genes.symbol_normalization import normalize_gene_name_list

        df["gene_list"] = df["gene_list"].apply(
            lambda xs: normalize_gene_name_list(list(xs), _map)
            if isinstance(xs, list)
            else normalize_gene_name_list([], _map)
        )

    # Parse distance to TSS (may have multiple values)
    df["distToTSS_parsed"] = df["distToTSS"].apply(_parse_dist_to_tss)
    
    # Compute probe center (for point-based operations)
    df["center"] = ((df["start"] + df["end"]) // 2).astype("Int64")
    
    print(f"  Loaded {len(df)} probes")
    print(f"  Chromosomes: {df['chrom'].nunique()}")
    print(f"  In CGI: {df['in_CGI'].sum()} ({100*df['in_CGI'].mean():.1f}%)")
    
    return df


def _parse_gene_list(val) -> List[str]:
    """Parse semicolon/comma separated gene names."""
    if pd.isna(val) or val == "":
        return []
    # Handle both ; and , separators
    genes = str(val).replace(",", ";").split(";")
    return [g.strip() for g in genes if g.strip()]


def _parse_dist_to_tss(val) -> List[int]:
    """Parse semicolon-separated distances to TSS."""
    if pd.isna(val) or val == "":
        return []
    try:
        parts = str(val).replace(",", ";").split(";")
        return [int(float(p.strip())) for p in parts if p.strip()]
    except (ValueError, TypeError):
        return []


# =============================================================================
# PROMOTER ANNOTATION
# =============================================================================

def annotate_probes_with_promoters(
    probes: pd.DataFrame,
    genes: pd.DataFrame,
    upstream_bp: int = THRESHOLDS.meth_promoter_upstream_bp,
    downstream_bp: int = THRESHOLDS.meth_promoter_downstream_bp,
    gene_panel: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Annotate probes with promoter overlap for panel genes.
    
    Adds columns:
        - in_promoter: bool - probe overlaps any panel gene promoter
        - promoter_genes: list[str] - gene names whose promoters overlap
        - promoter_gene_ids: list[str] - gene IDs whose promoters overlap
    
    Args:
        probes: Probe reference DataFrame
        genes: Gene DataFrame with chrom, start, end, strand, gene_name, gene_id
        upstream_bp: Promoter upstream of TSS
        downstream_bp: Promoter downstream of TSS
        gene_panel: List of gene names to consider (None = all)
    
    Returns:
        DataFrame with promoter annotations added
    """
    print("Annotating probes with gene promoters...")
    
    probes = probes.copy()
    
    # Filter to panel genes if specified
    panel_genes = genes.copy()
    if gene_panel is not None:
        panel_genes = panel_genes[panel_genes["gene_name"].isin(gene_panel)]
    
    # Filter to gene-level features
    if "feature" in panel_genes.columns:
        panel_genes = panel_genes[panel_genes["feature"] == "gene"]
    
    if panel_genes.empty:
        probes["in_promoter"] = False
        probes["promoter_genes"] = [[] for _ in range(len(probes))]
        probes["promoter_gene_ids"] = [[] for _ in range(len(probes))]
        return probes
    
    # Compute TSS if not present
    if "tss" not in panel_genes.columns:
        panel_genes["tss"] = np.where(
            panel_genes["strand"] == "+",
            panel_genes["start"],
            panel_genes["end"]
        )
    
    # Compute promoter regions
    panel_genes["prom_start"] = np.where(
        panel_genes["strand"] == "+",
        (panel_genes["tss"] - upstream_bp).clip(lower=0),
        (panel_genes["tss"] - downstream_bp).clip(lower=0),
    )
    panel_genes["prom_end"] = np.where(
        panel_genes["strand"] == "+",
        panel_genes["tss"] + downstream_bp,
        panel_genes["tss"] + upstream_bp,
    )
    
    # Build probe-to-promoter mapping
    promoter_map = _map_probes_to_intervals(
        probes,
        panel_genes,
        interval_start_col="prom_start",
        interval_end_col="prom_end",
        name_col="gene_name",
        id_col="gene_id" if "gene_id" in panel_genes.columns else None,
    )
    
    # Apply mapping
    probes["promoter_genes"] = probes["probeID"].map(
        lambda x: promoter_map.get(x, {}).get("names", [])
    )
    probes["promoter_gene_ids"] = probes["probeID"].map(
        lambda x: promoter_map.get(x, {}).get("ids", [])
    )
    probes["in_promoter"] = probes["promoter_genes"].apply(len) > 0
    
    n_in_promoter = probes["in_promoter"].sum()
    print(f"  Probes in promoters: {n_in_promoter} ({100*n_in_promoter/len(probes):.2f}%)")
    
    return probes


# =============================================================================
# GENE BODY ANNOTATION
# =============================================================================

def annotate_probes_with_gene_bodies(
    probes: pd.DataFrame,
    genes: pd.DataFrame,
    gene_panel: Optional[List[str]] = None,
    exclude_promoter: bool = True,
    upstream_bp: int = THRESHOLDS.meth_promoter_upstream_bp,
    downstream_bp: int = THRESHOLDS.meth_promoter_downstream_bp,
) -> pd.DataFrame:
    """
    Annotate probes with gene body overlap.
    
    Gene body is defined as the gene interval excluding the promoter region.
    
    Adds columns:
        - in_gene_body: bool
        - gene_body_genes: list[str]
    
    Args:
        probes: Probe reference DataFrame
        genes: Gene DataFrame
        gene_panel: List of gene names to consider
        exclude_promoter: If True, gene body excludes promoter region
        upstream_bp: Promoter upstream (for exclusion)
        downstream_bp: Promoter downstream (for exclusion)
    """
    print("Annotating probes with gene bodies...")
    
    probes = probes.copy()
    
    # Filter to panel genes
    panel_genes = genes.copy()
    if gene_panel is not None:
        panel_genes = panel_genes[panel_genes["gene_name"].isin(gene_panel)]
    
    if "feature" in panel_genes.columns:
        panel_genes = panel_genes[panel_genes["feature"] == "gene"]
    
    if panel_genes.empty:
        probes["in_gene_body"] = False
        probes["gene_body_genes"] = [[] for _ in range(len(probes))]
        return probes
    
    # For gene body, we want gene interval but possibly excluding promoter
    panel_genes = panel_genes.copy()
    
    if exclude_promoter:
        # Compute TSS
        panel_genes["tss"] = np.where(
            panel_genes["strand"] == "+",
            panel_genes["start"],
            panel_genes["end"]
        )
        
        # Adjust gene body to exclude promoter
        # For + strand: body starts at TSS + downstream_bp
        # For - strand: body ends at TSS - downstream_bp
        panel_genes["body_start"] = np.where(
            panel_genes["strand"] == "+",
            panel_genes["tss"] + downstream_bp,
            panel_genes["start"],
        )
        panel_genes["body_end"] = np.where(
            panel_genes["strand"] == "+",
            panel_genes["end"],
            panel_genes["tss"] - downstream_bp,
        )
        
        # Ensure valid intervals
        panel_genes = panel_genes[panel_genes["body_end"] > panel_genes["body_start"]]
    else:
        panel_genes["body_start"] = panel_genes["start"]
        panel_genes["body_end"] = panel_genes["end"]
    
    if panel_genes.empty:
        probes["in_gene_body"] = False
        probes["gene_body_genes"] = [[] for _ in range(len(probes))]
        return probes
    
    # Build mapping
    body_map = _map_probes_to_intervals(
        probes,
        panel_genes,
        interval_start_col="body_start",
        interval_end_col="body_end",
        name_col="gene_name",
    )
    
    probes["gene_body_genes"] = probes["probeID"].map(
        lambda x: body_map.get(x, {}).get("names", [])
    )
    probes["in_gene_body"] = probes["gene_body_genes"].apply(len) > 0
    
    n_in_body = probes["in_gene_body"].sum()
    print(f"  Probes in gene bodies: {n_in_body} ({100*n_in_body/len(probes):.2f}%)")
    
    return probes

# =============================================================================
# GENE CDS ANNOTATION
# =============================================================================
def annotate_probes_with_gene_cds(
    probes: pd.DataFrame,
    genes_cds: pd.DataFrame,
    gene_panel: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Annotate probes with gene CDS (coding sequence) overlap.
    
    Adds columns:
        - in_cds: bool
        - cds_genes: list[str]
    """
    print("Annotating probes with gene CDS...")
    
    probes = probes.copy()
    
    # Filter to CDS features from protein-coding transcripts
    cds_df = genes_cds[
        (genes_cds["feature"] == "CDS") &
        (genes_cds["transcript_type"] == "protein_coding")
    ].copy()
    
    if gene_panel is not None:
        cds_df = cds_df[cds_df["gene_name"].isin(gene_panel)]
    
    if cds_df.empty:
        probes["in_cds"] = False
        probes["cds_genes"] = [[] for _ in range(len(probes))]
        print("  No CDS features found after filtering")
        return probes
    
    # Ensure numeric coordinates
    cds_df["start"] = pd.to_numeric(cds_df["start"], errors="coerce")
    cds_df["end"] = pd.to_numeric(cds_df["end"], errors="coerce")
    cds_df = cds_df.dropna(subset=["start", "end"])
    
    # Build mapping
    cds_map = _map_probes_to_intervals(
        probes,
        cds_df,
        interval_start_col="start",
        interval_end_col="end",
        name_col="gene_name",
    )
    
    probes["cds_genes"] = probes["probeID"].map(
        lambda x: list(set(cds_map.get(x, {}).get("names", [])))
    )
    probes["in_cds"] = probes["cds_genes"].apply(len) > 0
    
    n_in_cds = probes["in_cds"].sum()
    print(f"  Probes in CDS regions: {n_in_cds} ({100*n_in_cds/len(probes):.2f}%)")
    
    return probes


# =============================================================================
# CCRE ANNOTATION
# =============================================================================

def _lower_col_lookup(df: pd.DataFrame) -> Dict[str, str]:
    return {str(c).lower(): c for c in df.columns}


def _ccre_id_column_nonempty(cc: pd.DataFrame, col: str) -> bool:
    if col not in cc.columns:
        return False
    s = cc[col].astype(str).str.strip()
    return bool(((s != "") & (s.str.lower() != "nan")).any())


def resolve_methylation_ccre_id_column(
    ccres: pd.DataFrame,
    *,
    ccre_id_col: str = "cCRE_id",
    encode_id_col: str = "ENCODE_id",
) -> Tuple[pd.DataFrame, str]:
    """
    Choose a stable per-element ID column for probe ↔ cCRE overlap.

    ``regulatory_element_focus.csv`` and raw SCREEN exports differ: some tables
    omit ``cCRE_id`` / ``ENCODE_id``. We try common names (case-insensitive),
    then fall back to ``chrom:start-end`` on a copy of the frame.
    """
    cc = ccres.copy()
    low = _lower_col_lookup(cc)
    if "chrom" not in cc.columns:
        for key in ("chromosome", "seqname", "seqnames"):
            if key in low:
                cc = cc.rename(columns={low[key]: "chrom"})
                low = _lower_col_lookup(cc)
                break

    def _resolve_exact(name: str) -> Optional[str]:
        if name in cc.columns:
            return name
        ln = name.lower()
        return low.get(ln)

    ordered_prefs: List[str] = []
    for pref in (ccre_id_col, encode_id_col):
        c = _resolve_exact(pref)
        if c is not None and c not in ordered_prefs:
            ordered_prefs.append(c)
    for n in (
        "element_id",
        "elementId",
        "screen_id",
        "ccre_accession",
        "accession",
        "encode_accession",
        "name",
        "id",
    ):
        c = _resolve_exact(n)
        if c is not None and c not in ordered_prefs:
            ordered_prefs.append(c)
    for col in cc.columns:
        s = str(col).lower()
        if "ccre" in s and "id" in s and col not in ordered_prefs:
            ordered_prefs.append(col)

    for c in ordered_prefs:
        if _ccre_id_column_nonempty(cc, c):
            return cc, c

    # Any declared ID column even if sparse (better than failing outright).
    for c in ordered_prefs:
        if c in cc.columns:
            return cc, c

    if "chrom" in cc.columns and "start" in cc.columns and "end" in cc.columns:
        ss = pd.to_numeric(cc["start"], errors="coerce")
        ee = pd.to_numeric(cc["end"], errors="coerce")
        mask = ss.notna() & ee.notna()
        sid = pd.Series("", index=cc.index, dtype=object)
        sid.loc[mask] = (
            cc.loc[mask, "chrom"].astype(str)
            + ":"
            + ss.loc[mask].astype(int).astype(str)
            + "-"
            + ee.loc[mask].astype(int).astype(str)
        )
        cc["_meth_ccre_id"] = sid
        return cc, "_meth_ccre_id"

    raise ValueError(
        "Regulatory element / cCRE table has no usable ID column "
        f"(tried {ccre_id_col!r}, {encode_id_col!r}, and common alternatives) "
        "and could not synthesize chrom:start-end (need chrom/chromosome + start + end). "
        f"Columns: {list(cc.columns)}"
    )


def annotate_probes_with_ccres(
    probes: pd.DataFrame,
    ccres: pd.DataFrame,
    ccre_id_col: str = "cCRE_id",
    encode_id_col: str = "ENCODE_id",
) -> pd.DataFrame:
    """
    Annotate probes with overlapping cCREs.
    
    Adds columns:
        - overlapping_ccres: list[str] - cCRE IDs that overlap
        - ccre_types: list[str] - types of overlapping cCREs
        - n_overlapping_ccres: int
    
    Args:
        probes: Probe reference DataFrame
        ccres: cCRE DataFrame with chrom, start, end, type
        ccre_id_col: Column name for cCRE ID
        encode_id_col: Column name for ENCODE ID
    """
    print("Annotating probes with cCREs...")
    
    probes = probes.copy()
    
    if ccres.empty:
        probes["overlapping_ccres"] = [[] for _ in range(len(probes))]
        probes["ccre_types"] = [[] for _ in range(len(probes))]
        probes["n_overlapping_ccres"] = 0
        return probes

    cc_work, id_col = resolve_methylation_ccre_id_column(
        ccres, ccre_id_col=ccre_id_col, encode_id_col=encode_id_col
    )
    cc_work, _ = harmonize_chrom_column(cc_work, chrom_col="chrom")
    print(f"  Using element ID column: {id_col!r}")

    # Build mapping
    ccre_map = _map_probes_to_ccres(probes, cc_work, id_col)
    
    probes["overlapping_ccres"] = probes["probeID"].map(
        lambda x: ccre_map.get(x, {}).get("ids", [])
    )
    probes["ccre_types"] = probes["probeID"].map(
        lambda x: ccre_map.get(x, {}).get("types", [])
    )
    probes["n_overlapping_ccres"] = probes["overlapping_ccres"].apply(len)
    
    n_with_ccre = (probes["n_overlapping_ccres"] > 0).sum()
    print(f"  Probes overlapping cCREs: {n_with_ccre} ({100*n_with_ccre/len(probes):.2f}%)")
    
    return probes


def _map_probes_to_ccres(
    probes: pd.DataFrame,
    ccres: pd.DataFrame,
    id_col: str,
) -> Dict[str, Dict[str, List]]:
    """
    Map probes to overlapping cCREs.
    
    Returns:
        dict: probeID -> {"ids": [...], "types": [...]}
    """
    result = {}
    
    # Process per chromosome
    for chrom in probes["chrom"].dropna().unique():
        probes_chr = probes[probes["chrom"] == chrom]
        ccres_chr = ccres[ccres["chrom"] == chrom]
        
        if ccres_chr.empty:
            continue
        
        # Sort cCREs for efficient searching
        ccres_chr = ccres_chr.sort_values("start").reset_index(drop=True)
        ccre_starts = ccres_chr["start"].values
        ccre_ends = ccres_chr["end"].values
        ccre_ids = ccres_chr[id_col].values
        type_col = None
        for tc in ("type", "raw_type", "primary_class", "element_class"):
            if tc in ccres_chr.columns:
                type_col = tc
                break
        if type_col is not None:
            ccre_types = ccres_chr[type_col].values
        else:
            ccre_types = np.array([""] * len(ccres_chr), dtype=object)
        
        for _, probe in probes_chr.iterrows():
            probe_pos = probe["start"]  # CpG is essentially a point
            
            # Find overlapping cCREs
            overlaps = (ccre_starts <= probe_pos) & (ccre_ends >= probe_pos)
            
            if overlaps.any():
                overlap_ids = ccre_ids[overlaps].tolist()
                overlap_types = [t for t in ccre_types[overlaps] if t]
                result[probe["probeID"]] = {
                    "ids": overlap_ids,
                    "types": overlap_types,
                }
    
    return result


# =============================================================================
# LNCRNA ANNOTATION
# =============================================================================

def annotate_probes_with_lncrnas(
    probes: pd.DataFrame,
    lncrnas: pd.DataFrame,
    upstream_bp: int = THRESHOLDS.meth_promoter_upstream_bp,
    downstream_bp: int = THRESHOLDS.meth_promoter_downstream_bp,
    lncrna_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Annotate probes with lncRNA promoter overlap.
    
    Adds columns:
        - in_lncrna_promoter: bool
        - lncrna_promoter_genes: list[str]
    
    Args:
        probes: Probe reference DataFrame
        lncrnas: lncRNA DataFrame with gene-like structure
        upstream_bp: Promoter upstream of TSS
        downstream_bp: Promoter downstream of TSS
        lncrna_names: Specific lncRNAs to consider
    """
    print("Annotating probes with lncRNA promoters...")
    
    probes = probes.copy()
    
    # Filter lncRNAs
    panel_lnc = lncrnas.copy()
    if lncrna_names is not None:
        panel_lnc = panel_lnc[panel_lnc["gene_name"].isin(lncrna_names)]
    
    # GENCODE-style tables include multiple feature types; gene-centric tables omit ``feature``.
    if "feature" in panel_lnc.columns:
        panel_lnc = panel_lnc[panel_lnc["feature"] == "gene"]
    
    if panel_lnc.empty:
        probes["in_lncrna_promoter"] = False
        probes["lncrna_promoter_genes"] = [[] for _ in range(len(probes))]
        return probes
    
    # Compute TSS and promoter
    panel_lnc["tss"] = np.where(
        panel_lnc["strand"] == "+",
        panel_lnc["start"],
        panel_lnc["end"]
    )
    panel_lnc["prom_start"] = np.where(
        panel_lnc["strand"] == "+",
        (panel_lnc["tss"] - upstream_bp).clip(lower=0),
        (panel_lnc["tss"] - downstream_bp).clip(lower=0),
    )
    panel_lnc["prom_end"] = np.where(
        panel_lnc["strand"] == "+",
        panel_lnc["tss"] + downstream_bp,
        panel_lnc["tss"] + upstream_bp,
    )
    
    # Build mapping
    lnc_map = _map_probes_to_intervals(
        probes,
        panel_lnc,
        interval_start_col="prom_start",
        interval_end_col="prom_end",
        name_col="gene_name",
    )
    
    probes["lncrna_promoter_genes"] = probes["probeID"].map(
        lambda x: lnc_map.get(x, {}).get("names", [])
    )
    probes["in_lncrna_promoter"] = probes["lncrna_promoter_genes"].apply(len) > 0
    
    n_in_lnc = probes["in_lncrna_promoter"].sum()
    print(f"  Probes in lncRNA promoters: {n_in_lnc} ({100*n_in_lnc/len(probes):.2f}%)")
    
    return probes


# =============================================================================
# ATAC PEAK ANNOTATION
# =============================================================================

def annotate_probes_with_atac(
    probes: pd.DataFrame,
    atac_peaks: pd.DataFrame,
    peak_id_col: str = "peak_id",
) -> pd.DataFrame:
    """
    Annotate probes with overlapping ATAC peaks.
    
    Adds columns:
        - overlapping_atac_peaks: list[str]
        - n_overlapping_atac: int
    
    Args:
        probes: Probe reference DataFrame
        atac_peaks: ATAC peaks DataFrame with chrom, start, end
        peak_id_col: Column name for peak ID
    """
    print("Annotating probes with ATAC peaks...")
    
    probes = probes.copy()
    
    if atac_peaks.empty:
        probes["overlapping_atac_peaks"] = [[] for _ in range(len(probes))]
        probes["n_overlapping_atac"] = 0
        return probes
    
    # Build mapping
    atac_map = _map_probes_to_intervals(
        probes,
        atac_peaks,
        interval_start_col="start",
        interval_end_col="end",
        name_col=peak_id_col,
    )
    
    probes["overlapping_atac_peaks"] = probes["probeID"].map(
        lambda x: atac_map.get(x, {}).get("names", [])
    )
    probes["n_overlapping_atac"] = probes["overlapping_atac_peaks"].apply(len)
    
    n_with_atac = (probes["n_overlapping_atac"] > 0).sum()
    print(f"  Probes overlapping ATAC peaks: {n_with_atac} ({100*n_with_atac/len(probes):.2f}%)")
    
    return probes


# =============================================================================
# GENERIC INTERVAL MAPPING HELPER
# =============================================================================

def _map_probes_to_intervals(
    probes: pd.DataFrame,
    intervals: pd.DataFrame,
    interval_start_col: str,
    interval_end_col: str,
    name_col: str,
    id_col: Optional[str] = None,
) -> Dict[str, Dict[str, List]]:
    """
    Generic helper to map probes (points) to intervals.
    
    Returns:
        dict: probeID -> {"names": [...], "ids": [...]}
    """
    result = {}
    
    # Process per chromosome
    for chrom in probes["chrom"].dropna().unique():
        probes_chr = probes[probes["chrom"] == chrom]
        intervals_chr = intervals[intervals["chrom"] == chrom]
        
        if intervals_chr.empty:
            continue
        
        # Convert to numpy for speed
        int_starts = intervals_chr[interval_start_col].values
        int_ends = intervals_chr[interval_end_col].values
        int_names = intervals_chr[name_col].values
        int_ids = intervals_chr[id_col].values if id_col and id_col in intervals_chr.columns else None
        
        for _, probe in probes_chr.iterrows():
            probe_pos = probe["start"]  # CpG is a point
            
            # Find overlapping intervals
            overlaps = (int_starts <= probe_pos) & (int_ends >= probe_pos)
            
            if overlaps.any():
                names = int_names[overlaps].tolist()
                ids = int_ids[overlaps].tolist() if int_ids is not None else []
                result[probe["probeID"]] = {"names": names, "ids": ids}
    
    return result


# =============================================================================
# TAD ANNOTATION (REUSES EXISTING INFRASTRUCTURE)
# =============================================================================

def annotate_probes_with_tads(
    probes: pd.DataFrame,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
    *,
    slim_tad_payload: Optional[bool] = None,
) -> pd.DataFrame:
    """
    Annotate probes with TAD domain context.
    
    Reuses the pipeline's TAD annotation infrastructure.
    Probes use ``kind="interval"`` (CpG coordinates). **Default:** *slim* TAD
    payload (compact dicts) so EPIC-scale probes × many biosamples stay bounded;
    set ``APM_METH_TAD_SLIM=0`` or ``slim_tad_payload=False`` for full gene-style
    payloads (``domains`` + full boundary objects).

    Args:
        probes: Probe reference DataFrame
        tad_domains: TAD domains DataFrame
        domain_flanks: Domain flanks DataFrame
        boundaries: Boundaries DataFrame
        biosample: Biosample name for the TAD data
        slim_tad_payload: If None, slim unless ``APM_METH_TAD_SLIM`` is ``0``/``false``/``no``.
    
    Returns:
        DataFrame with TAD_domains column updated
    """
    # Import from main pipeline - lazy import to avoid circular deps
    try:
        from ..tad_annotation import annotate_df_with_tads
        from ..tad_annotation.annotator import slim_tad_payload_for_probe
    except ImportError:
        print("  Warning: TAD annotation not available, skipping")
        return probes

    if slim_tad_payload is None:
        slim_tad_payload = os.environ.get("APM_METH_TAD_SLIM", "1").strip().lower() not in (
            "0",
            "false",
            "no",
        )

    print(f"Annotating probes with TAD context ({biosample})...")

    return annotate_df_with_tads(
        probes,
        kind="interval",
        tad_domains=tad_domains,
        domain_flanks=domain_flanks,
        boundaries=boundaries,
        biosample=biosample,
        out_col="TAD_domains",
        payload_postprocess=(slim_tad_payload_for_probe if slim_tad_payload else None),
        use_vectorized_interval=True,
    )


# =============================================================================
# COMBINED ANNOTATION
# =============================================================================

def annotate_probes_full(
    probes: pd.DataFrame,
    genes: pd.DataFrame,
    genes_cds: pd.DataFrame,
    ccres: pd.DataFrame,
    lncrnas: Optional[pd.DataFrame] = None,
    lncrnas_cds: Optional[pd.DataFrame] = None,
    atac_peaks: Optional[pd.DataFrame] = None,
    gene_panel: Optional[List[str]] = None,
    lncrna_names: Optional[List[str]] = None,
    upstream_bp: int = THRESHOLDS.meth_promoter_upstream_bp,
    downstream_bp: int = THRESHOLDS.meth_promoter_downstream_bp,
) -> pd.DataFrame:
    """
    Apply all annotations to probe reference in one call.
    
    This is the main entry point for building the annotated probe reference.
    """
    print("=" * 60)
    print("ANNOTATING PROBE REFERENCE")
    print("=" * 60)
    
    # Gene promoters
    probes = annotate_probes_with_promoters(
        probes, genes,
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
        gene_panel=gene_panel,
    )
    
    # Gene bodies
    probes = annotate_probes_with_gene_bodies(
        probes, genes,
        gene_panel=gene_panel,
        exclude_promoter=True,
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
    )
    
    # Gene CDS
    probes = annotate_probes_with_gene_cds(
        probes, genes_cds,
        gene_panel=gene_panel,
    )
    
    # cCREs
    probes = annotate_probes_with_ccres(probes, ccres)
    
  
    
    # lncRNAs
    if lncrnas is not None and not lncrnas.empty:
        probes = annotate_probes_with_lncrnas(
            probes, lncrnas,
            upstream_bp=upstream_bp,
            downstream_bp=downstream_bp,
            lncrna_names=lncrna_names,
        )
    else:
        probes["in_lncrna_promoter"] = False
        probes["lncrna_promoter_genes"] = [[] for _ in range(len(probes))]

 
    
    # ATAC peaks
    if atac_peaks is not None and not atac_peaks.empty:
        probes = annotate_probes_with_atac(probes, atac_peaks)
    else:
        probes["overlapping_atac_peaks"] = [[] for _ in range(len(probes))]
        probes["n_overlapping_atac"] = 0
    
    print("=" * 60)
    print("PROBE ANNOTATION COMPLETE")
    print(f"  Total probes: {len(probes)}")
    print(f"  In gene promoters: {probes['in_promoter'].sum()}")
    print(f"  In gene bodies: {probes['in_gene_body'].sum()}")
    print(f"  Overlapping cCREs: {(probes['n_overlapping_ccres'] > 0).sum()}")
    print(f"  In lncRNA promoters: {probes['in_lncrna_promoter'].sum()}")
    print(f"  Overlapping ATAC: {(probes['n_overlapping_atac'] > 0).sum()}")
    print("=" * 60)
    
    return probes
