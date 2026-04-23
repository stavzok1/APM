from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline.atac_peaks.peak_table import load_atac_peaks_annotated
from pipeline.config import OUTPUT_SUBDIRS, PATHS, PIPELINE_GENE_PANEL
from pipeline.genes.gene_loader import load_genes, filter_genes_by_names, add_promoter_columns
from pipeline.Methylation.probe_loader import load_probe_reference
from pipeline.regulatory_elements import load_ccres

from .loader import discover_tad_sources, load_tad_source
from .annotator import annotate_genes_with_tads
from .mirroring import mirror_genes_into_domains


def _overlap_hits_by_chrom(
    boundaries: pd.DataFrame,
    features: pd.DataFrame,
    *,
    feature_id_col: str,
    chrom_col_b: str = "chrom",
    start_col_b: str = "start",
    end_col_b: str = "end",
    chrom_col_f: str = "chrom",
    start_col_f: str = "start",
    end_col_f: str = "end",
) -> Dict[str, List[str]]:
    """
    Return mapping boundary_id -> list[feature_id] for interval overlaps.

    Efficient enough for (n_boundaries ~ few 10k, n_features ~ few 100k) by
    per-chromosome sorting + searchsorted on feature starts.
    """
    out: Dict[str, List[str]] = {bid: [] for bid in boundaries["boundary_id"].astype(str).tolist()}
    if boundaries.empty or features.empty:
        return out

    b = boundaries[[chrom_col_b, start_col_b, end_col_b, "boundary_id"]].copy()
    f = features[[chrom_col_f, start_col_f, end_col_f, feature_id_col]].copy()

    b[chrom_col_b] = b[chrom_col_b].astype(str)
    f[chrom_col_f] = f[chrom_col_f].astype(str)
    b[start_col_b] = pd.to_numeric(b[start_col_b], errors="coerce").astype("Int64")
    b[end_col_b] = pd.to_numeric(b[end_col_b], errors="coerce").astype("Int64")
    f[start_col_f] = pd.to_numeric(f[start_col_f], errors="coerce").astype("Int64")
    f[end_col_f] = pd.to_numeric(f[end_col_f], errors="coerce").astype("Int64")
    b = b.dropna(subset=[start_col_b, end_col_b])
    f = f.dropna(subset=[start_col_f, end_col_f])

    for chrom, bchr in b.groupby(chrom_col_b, sort=False):
        fchr = f[f[chrom_col_f] == chrom]
        if fchr.empty:
            continue
        fchr = fchr.sort_values(start_col_f, kind="mergesort")
        f_starts = fchr[start_col_f].to_numpy(dtype=np.int64, copy=False)
        f_ends = fchr[end_col_f].to_numpy(dtype=np.int64, copy=False)
        f_ids = fchr[feature_id_col].astype(str).to_numpy(copy=False)

        for _, row in bchr.iterrows():
            bs = int(row[start_col_b])
            be = int(row[end_col_b])
            bid = str(row["boundary_id"])
            # candidates with start <= be
            hi = int(np.searchsorted(f_starts, be, side="right"))
            if hi <= 0:
                continue
            cand_ids = f_ids[:hi]
            cand_ends = f_ends[:hi]
            keep = cand_ends >= bs
            if np.any(keep):
                out[bid].extend(cand_ids[keep].tolist())
    return out


def _adjacent_domains_for_boundary(
    domains: pd.DataFrame,
    boundary_pos: int,
    *,
    chrom: str,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Find left/right adjacent domains around a boundary anchor position.
    """
    d = domains[domains["chrom"].astype(str) == str(chrom)].copy()
    if d.empty:
        return None, None
    d["start"] = pd.to_numeric(d["start"], errors="coerce").astype("Int64")
    d["end"] = pd.to_numeric(d["end"], errors="coerce").astype("Int64")
    d = d.dropna(subset=["start", "end"]).sort_values(["start", "end"], kind="mergesort")
    # left: max end <= pos
    left = d[d["end"] <= boundary_pos].tail(1)
    # right: min start >= pos
    right = d[d["start"] >= boundary_pos].head(1)
    left_id = str(left["domain_id"].iloc[0]) if not left.empty else None
    right_id = str(right["domain_id"].iloc[0]) if not right.empty else None
    return left_id, right_id


@dataclass(frozen=True)
class BoundaryEnrichmentPaths:
    processed_dir: Path
    out_path: Path


def build_boundaries_enriched_table(
    *,
    processed_dir: Path,
    biosamples: Optional[Iterable[str]] = None,
    gene_panel: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Build a *single* table of boundaries enriched with:
    - overlapping cCREs (IDs)
    - overlapping ATAC peaks (interval IDs from atac_peaks_annotated)
    - overlapping methylation probes (probeID)
    - gene lists in left/right adjacent domains (biosample-specific)

    Output granularity: one row per boundary (across all biosamples), with a `biosample` column.
    """
    processed_dir = Path(processed_dir)
    avail = discover_tad_sources(processed_dir, required_files=True)
    if biosamples is None:
        biosamples = list(avail.keys())
    else:
        biosamples = [b for b in biosamples if b in avail]

    # Features (genome-wide; we’ll overlap per biosample boundary sets)
    ccres = load_ccres(PATHS.ccre_csv)[["chrom", "start", "end", "cCRE_id"]].copy()
    ccres.rename(columns={"cCRE_id": "ccre_id"}, inplace=True)

    atac_path = Path(PATHS.working_dir) / OUTPUT_SUBDIRS["atac_peaks"] / "atac_peaks_annotated.parquet"
    atac = load_atac_peaks_annotated(atac_path, columns=["chrom", "start", "end"]).copy()
    atac = atac.reset_index(drop=True)
    atac["atac_peak_id"] = atac.index.astype(str)

    probes = load_probe_reference(PATHS.methylation_probe_reference)[
        ["chrom", "start", "end", "probeID"]
    ].copy()
    probes.rename(columns={"probeID": "probe_id"}, inplace=True)

    gene_panel = gene_panel or list(PIPELINE_GENE_PANEL)
    genes_all = load_genes(PATHS.gencode_gtf_pq)
    genes = filter_genes_by_names(genes_all, gene_panel)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)

    rows: List[pd.DataFrame] = []
    for b in biosamples:
        paths = avail[b]
        domains, boundaries, flanks = load_tad_source(paths)

        # Normalize columns in boundaries/domains
        boundaries = boundaries.copy()
        domains = domains.copy()
        # some sources use 'name' for ids; loader should normalize, but be defensive
        if "boundary_id" not in boundaries.columns and "name" in boundaries.columns:
            boundaries["boundary_id"] = boundaries["name"]
        if "domain_id" not in domains.columns and "name" in domains.columns:
            domains["domain_id"] = domains["name"]

        boundaries["boundary_id"] = boundaries["boundary_id"].astype(str)
        boundaries["chrom"] = boundaries["chrom"].astype(str)
        boundaries["start"] = pd.to_numeric(boundaries["start"], errors="coerce").astype("Int64")
        boundaries["end"] = pd.to_numeric(boundaries["end"], errors="coerce").astype("Int64")
        if "pos" in boundaries.columns:
            boundaries["pos"] = pd.to_numeric(boundaries["pos"], errors="coerce").astype("Int64")
        else:
            boundaries["pos"] = ((boundaries["start"] + boundaries["end"]) // 2).astype("Int64")

        # Overlap features
        ccre_hits = _overlap_hits_by_chrom(boundaries, ccres, feature_id_col="ccre_id")
        atac_hits = _overlap_hits_by_chrom(boundaries, atac, feature_id_col="atac_peak_id")
        probe_hits = _overlap_hits_by_chrom(boundaries, probes, feature_id_col="probe_id")

        # Genes in adjacent domains: annotate genes with this biosample’s TADs, then mirror to domains.
        g_ann = annotate_genes_with_tads(genes.copy(), domains, flanks, boundaries, biosample=b)
        dom_with_genes = mirror_genes_into_domains(domains.copy(), g_ann, biosample=b, mode="primary")
        dom_gene_map: Dict[str, List[str]] = {}
        if "gene_hits" in dom_with_genes.columns:
            for _, dr in dom_with_genes.iterrows():
                did = str(dr.get("domain_id"))
                payload = dr.get("gene_hits")
                if isinstance(payload, dict):
                    bucket = payload.get(b, []) or []
                    dom_gene_map[did] = sorted({str(x.get("gene_name")) for x in bucket if x.get("gene_name")})

        # Build output frame for this biosample
        out = boundaries[["boundary_id", "chrom", "start", "end", "pos"]].copy()
        out["biosample"] = b
        out["ccre_ids"] = out["boundary_id"].map(lambda x: ccre_hits.get(str(x), []))
        out["atac_peak_ids"] = out["boundary_id"].map(lambda x: atac_hits.get(str(x), []))
        out["methyl_probe_ids"] = out["boundary_id"].map(lambda x: probe_hits.get(str(x), []))

        def _left_right_genes(r) -> Tuple[List[str], List[str], Optional[str], Optional[str]]:
            left_d, right_d = _adjacent_domains_for_boundary(domains, int(r["pos"]), chrom=str(r["chrom"]))
            left_genes = dom_gene_map.get(left_d, []) if left_d else []
            right_genes = dom_gene_map.get(right_d, []) if right_d else []
            return left_genes, right_genes, left_d, right_d

        lr = out.apply(_left_right_genes, axis=1, result_type="expand")
        out["left_domain_genes"] = lr[0]
        out["right_domain_genes"] = lr[1]
        out["left_domain_id"] = lr[2]
        out["right_domain_id"] = lr[3]

        rows.append(out)

    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

