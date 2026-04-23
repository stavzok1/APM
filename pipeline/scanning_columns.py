from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

from pipeline.config import PIPELINE_GENE_PANEL, TIER1_LNCRNA_GENES


def _as_dict(x: Any) -> Dict[str, Any]:
    return x if isinstance(x, dict) else {}


def _as_list(x: Any) -> List[Any]:
    return x if isinstance(x, list) else []


def _safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _max_or_none(xs: Iterable[Optional[float]]) -> Optional[float]:
    m: Optional[float] = None
    for x in xs:
        if x is None:
            continue
        if m is None or x > m:
            m = x
    return m


def _count_screen_biosamples_with_strength(screen_block: Dict[str, Any], *, want: str) -> int:
    """
    screen_block: { per_biosample: { bio: { assay: {strength,...}, ... } } }
    want: "any" or "strong"
    """
    per = _as_dict(screen_block.get("per_biosample"))
    n = 0
    for _bio, assays in per.items():
        assays = _as_dict(assays)
        strengths = []
        for _assay, entry in assays.items():
            entry = _as_dict(entry)
            strengths.append(str(entry.get("strength") or "none"))
        if want == "any":
            if any(s != "none" for s in strengths):
                n += 1
        elif want == "strong":
            if any(s == "strong" for s in strengths):
                n += 1
    return n


def _screen_biosample_max_score(screen_block: Dict[str, Any], biosample: str) -> Optional[float]:
    per = _as_dict(screen_block.get("per_biosample"))
    assays = _as_dict(per.get(biosample))
    scores: List[Optional[float]] = []
    for _assay, entry in assays.items():
        entry = _as_dict(entry)
        scores.append(_safe_float(entry.get("score")))
    return _max_or_none(scores)


def _screen_conservation_breast_n_strong(screen_block: Dict[str, Any]) -> int:
    cons = _as_dict(screen_block.get("conservation_breast"))
    total = 0
    for _assay, entry in cons.items():
        entry = _as_dict(entry)
        try:
            total += int(entry.get("n_strong") or 0)
        except Exception:
            continue
    return int(total)


def _abc_stats(abc_enhancers: Any, *, celltype: Optional[str] = None) -> Tuple[Optional[float], int, bool]:
    """
    Returns: (max_abc_score, n_celltypes_strong, any_self_promoter)
    If celltype is provided, max_abc_score is restricted to that celltype.
    """
    enhs = _as_list(abc_enhancers)
    max_score: Optional[float] = None
    strong_celltypes = set()
    any_self_promoter = False

    for e in enhs:
        e = _as_dict(e)
        abc_full = _as_dict(e.get("ABC_full"))
        for ct, ct_entry in abc_full.items():
            if celltype is not None and ct != celltype:
                continue
            ct_entry = _as_dict(ct_entry)
            sc = _safe_float(ct_entry.get("ABC_score"))
            if sc is not None and (max_score is None or sc > max_score):
                max_score = sc
            if bool(ct_entry.get("is_strong")):
                strong_celltypes.add(ct)
            if bool(ct_entry.get("is_self_promoter")):
                any_self_promoter = True

    return max_score, len(strong_celltypes), any_self_promoter


def _hichip_stats(hichip_block: Any) -> Tuple[int, Optional[float], Dict[str, int]]:
    """
    Returns: (n_celltypes_with_loops, max_score_any, n_loops_by_celltype)
    """
    hb = _as_dict(hichip_block)
    n_with = 0
    max_score: Optional[float] = None
    n_loops: Dict[str, int] = {}
    for ct, entry in hb.items():
        entry = _as_dict(entry)
        n = int(entry.get("n_loops") or 0)
        n_loops[str(ct)] = n
        if n > 0:
            n_with += 1
        sc = _safe_float(entry.get("max_score"))
        if sc is not None and (max_score is None or sc > max_score):
            max_score = sc
    return n_with, max_score, n_loops


def _chip_stats(chip_hits: Any, *, cell_type: Optional[str] = None) -> Tuple[int, bool, bool, int]:
    """
    chip_hits (elem_focus): dict[tf][cell_type][source] = summary
    Returns: (n_tfs, has_STAT1, has_CTCF, n_tfs_in_cell_type)
    """
    ch = _as_dict(chip_hits)
    tfs = list(ch.keys())
    has_stat1 = "STAT1" in ch
    has_ctcf = "CTCF" in ch
    n_tfs_in_cell = 0
    if cell_type is not None:
        for tf, per_cell in ch.items():
            per_cell = _as_dict(per_cell)
            if cell_type in per_cell:
                n_tfs_in_cell += 1
    return len(tfs), has_stat1, has_ctcf, n_tfs_in_cell


def derive_elem_focus_scanning_columns(elem_focus: pd.DataFrame) -> pd.DataFrame:
    """
    Add scan_* scalar columns to the regulatory element focus table.

    The function is intentionally defensive: missing nested columns yield null/0 outputs.
    """
    df = elem_focus.copy()

    panel = set(PIPELINE_GENE_PANEL)
    tier1_lncrna = set(TIER1_LNCRNA_GENES)

    gene_links = df["gene_links"] if "gene_links" in df.columns else pd.Series([{}] * len(df))
    df["scan_gene_links_n_genes"] = gene_links.apply(lambda x: len(x) if isinstance(x, dict) else 0)
    df["scan_gene_links_any_panel_gene"] = gene_links.apply(
        lambda x: any(g in panel for g in x.keys()) if isinstance(x, dict) else False
    )
    df["scan_gene_links_any_tier1_lncrna"] = gene_links.apply(
        lambda x: any(g in tier1_lncrna for g in x.keys()) if isinstance(x, dict) else False
    )

    # SCREEN blocks (row-level)
    if "screen_exp" in df.columns:
        df["scan_screen_exp_n_biosamples_any"] = df["screen_exp"].apply(
            lambda x: _count_screen_biosamples_with_strength(_as_dict(x), want="any")
        )
        df["scan_screen_exp_n_biosamples_strong"] = df["screen_exp"].apply(
            lambda x: _count_screen_biosamples_with_strength(_as_dict(x), want="strong")
        )
        for bio in ("MCF7", "MCF10A"):
            df[f"scan_screen_exp_{bio}_max_score"] = df["screen_exp"].apply(
                lambda x, b=bio: _screen_biosample_max_score(_as_dict(x), b)
            )
        df["scan_screen_exp_breast_n_strong"] = df["screen_exp"].apply(
            lambda x: _screen_conservation_breast_n_strong(_as_dict(x))
        )
    else:
        df["scan_screen_exp_n_biosamples_any"] = 0
        df["scan_screen_exp_n_biosamples_strong"] = 0
        df["scan_screen_exp_MCF7_max_score"] = None
        df["scan_screen_exp_MCF10A_max_score"] = None
        df["scan_screen_exp_breast_n_strong"] = 0

    # ABC enhancers (row-level list)
    if "ABC_enhancers" in df.columns:
        df["scan_ABC_max_score_any_celltype"] = df["ABC_enhancers"].apply(
            lambda x: _abc_stats(x)[0]
        )
        df["scan_ABC_n_celltypes_strong"] = df["ABC_enhancers"].apply(
            lambda x: _abc_stats(x)[1]
        )
        df["scan_ABC_self_promoter_any"] = df["ABC_enhancers"].apply(
            lambda x: _abc_stats(x)[2]
        )
        df["scan_ABC_max_score_MCF7_ENCODE"] = df["ABC_enhancers"].apply(
            lambda x: _abc_stats(x, celltype="MCF7_ENCODE")[0]
        )
    else:
        df["scan_ABC_max_score_any_celltype"] = None
        df["scan_ABC_n_celltypes_strong"] = 0
        df["scan_ABC_self_promoter_any"] = False
        df["scan_ABC_max_score_MCF7_ENCODE"] = None

    # HiChIP (row-level dict)
    if "hichip" in df.columns:
        df["scan_hichip_n_celltypes_with_loops"] = df["hichip"].apply(lambda x: _hichip_stats(x)[0])
        df["scan_hichip_max_score_any"] = df["hichip"].apply(lambda x: _hichip_stats(x)[1])
        df["scan_hichip_n_loops_MCF7"] = df["hichip"].apply(lambda x: _hichip_stats(x)[2].get("MCF7", 0))
    else:
        df["scan_hichip_n_celltypes_with_loops"] = 0
        df["scan_hichip_max_score_any"] = None
        df["scan_hichip_n_loops_MCF7"] = 0

    # ChIP hits (row-level dict)
    if "chip_hits" in df.columns:
        df["scan_chip_hits_n_TFs"] = df["chip_hits"].apply(lambda x: _chip_stats(x)[0])
        df["scan_chip_hits_has_STAT1"] = df["chip_hits"].apply(lambda x: _chip_stats(x)[1])
        df["scan_chip_hits_has_CTCF"] = df["chip_hits"].apply(lambda x: _chip_stats(x)[2])
        df["scan_chip_hits_MCF7_n_TFs"] = df["chip_hits"].apply(lambda x: _chip_stats(x, cell_type="MCF7")[3])
    else:
        df["scan_chip_hits_n_TFs"] = 0
        df["scan_chip_hits_has_STAT1"] = False
        df["scan_chip_hits_has_CTCF"] = False
        df["scan_chip_hits_MCF7_n_TFs"] = 0

    return df


def derive_atac_peak_scanning_columns(peak_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add scan_* scalar columns to the ATAC peak annotated table.
    """
    df = peak_df.copy()
    panel = set(PIPELINE_GENE_PANEL)
    tier1_lncrna = set(TIER1_LNCRNA_GENES)

    # gene_links is an ATAC-specific nested dict (when present)
    if "gene_links" in df.columns:
        df["scan_gene_links_n_genes"] = df["gene_links"].apply(lambda x: len(x) if isinstance(x, dict) else 0)
        df["scan_gene_links_any_panel_gene"] = df["gene_links"].apply(
            lambda x: any(g in panel for g in x.keys()) if isinstance(x, dict) else False
        )
    else:
        df["scan_gene_links_n_genes"] = df.get("n_genes_total", 0)
        df["scan_gene_links_any_panel_gene"] = False

    if "lncrna_links" in df.columns:
        df["scan_lncrna_links_n_lncrnas"] = df["lncrna_links"].apply(lambda x: len(x) if isinstance(x, dict) else 0)
        df["scan_lncrna_links_any_tier1_lncrna"] = df["lncrna_links"].apply(
            lambda x: any(g in tier1_lncrna for g in x.keys()) if isinstance(x, dict) else False
        )
    else:
        df["scan_lncrna_links_n_lncrnas"] = df.get("n_lncrnas_total", 0)
        df["scan_lncrna_links_any_tier1_lncrna"] = False

    if "ccre_links" in df.columns:
        df["scan_ccre_links_n_ccres"] = df["ccre_links"].apply(lambda x: len(x) if isinstance(x, list) else 0)
    else:
        df["scan_ccre_links_n_ccres"] = df.get("n_ccres_total", 0)

    if "TAD_domains" in df.columns:
        df["scan_TAD_domains_n_biosamples"] = df["TAD_domains"].apply(lambda x: len(x) if isinstance(x, dict) else 0)
    else:
        df["scan_TAD_domains_n_biosamples"] = 0

    if "TAD_boundary_overlaps" in df.columns:
        df["scan_TAD_boundary_overlaps_any"] = df["TAD_boundary_overlaps"].apply(
            lambda x: any(v.get("overlaps_boundary", False) for v in x.values()) if isinstance(x, dict) else False
        )
    else:
        df["scan_TAD_boundary_overlaps_any"] = False

    return df

