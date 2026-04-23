"""
Canonical gene / miRNA symbol registries for the APM panel.

- **Flat map** ``alias -> canonical`` merges: HGNC (streamed), hand-curated
  per-canonical alias lists, UCSC RNA fixes, legacy dataset symbols, and
  **Ensembl gene ids** (``ENSG…`` / version-stripped) from the slim GENCODE table
  for every symbol in ``PANEL_ALIAS_SEED_SYMBOLS``.
- **Inverse map** ``canonical -> frozenset(aliases)`` is derived for inspection,
  reports, and tests.

**miRNA levels** (by design):

- **Xena expression matrix rows** are **MIMAT** accessions. Resolve those first via
  ``data/miRNA/mirna_mature_loci.csv`` (``mature_accession`` → ``mirbase_mature_id``),
  yielding **mature-arm** ids (``hsa-miR-…-5p`` / ``-3p``).
- **TargetScan / miRTarBase / literature** often use short **arm** names (``miR-155-5p``)
  or **family**-style MIR symbols; those are folded into the same **arm-level**
  canonical ids listed in ``TIER_MIRNA_ARM_CANONICAL_TO_ALIASES`` and
  ``MANUAL_MIRNA_CANONICAL_TO_ALIASES``.

All modules should obtain the flat gene map via ``get_gene_flat_alias_map()`` rather
than rebuilding HGNC merges ad hoc.
"""

from __future__ import annotations

import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from collections.abc import Collection
from typing import Dict, FrozenSet, Mapping, Optional, Tuple

import pandas as pd

from .gene_aliases import (
    build_rna_expression_symbol_mapping,
    merge_symbol_maps,
    strip_ensembl_version,
)

# Hand-curated extras (HGNC file is still authoritative for most symbols).
# Keys are **canonical** GENCODE/HGNC symbols used in ``pipeline/config.py``.
MANUAL_CANONICAL_TO_ALIASES: Dict[str, Tuple[str, ...]] = {
    "NECTIN2": ("PVRL2", "PRR2", "CD112", "HVEB"),
    "NCR3LG1": ("DKFZp686O24166",),
    "STING1": ("TMEM173", "STING", "MITA", "MPYS", "ERIS"),
    "CGAS": ("MB21D1", "cGAS"),
}


def _manual_alias_to_canonical() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for canon, aliases in MANUAL_CANONICAL_TO_ALIASES.items():
        for a in aliases:
            if a and a != canon:
                out[str(a).strip()] = canon
    return out


def _invert_flat(flat: Mapping[str, str]) -> Dict[str, FrozenSet[str]]:
    inv: dict[str, set[str]] = defaultdict(set)
    for alias, canon in flat.items():
        inv[canon].add(alias)
        inv[canon].add(canon)
    return {c: frozenset(xs) for c, xs in inv.items()}


def _max_lines_key(max_lines: Optional[int]) -> int:
    return -1 if max_lines is None else int(max_lines)


def _ensembl_gene_id_flat_map(
    gencode_path: Path,
    seed_symbols: Collection[str],
) -> Dict[str, str]:
    """
    Map ``ENSG…`` (with or without version) → **GENCODE gene_name** for rows whose
    symbol is in ``seed_symbols``.

    HGNC streaming intentionally skips Ensembl ids in ``content``; this layer
    restores them so matrices keyed by Ensembl still join to canonical symbols.
    """
    p = Path(gencode_path)
    if not p.is_file() or not seed_symbols:
        return {}
    out: Dict[str, str] = {}
    try:
        if p.suffix.lower() == ".parquet":
            df = pd.read_parquet(p, columns=["gene_name", "gene_id"])
        else:
            df = pd.read_csv(p, usecols=lambda c: c in ("gene_name", "gene_id"))
    except (ValueError, KeyError, OSError):
        try:
            df = pd.read_parquet(p) if p.suffix.lower() == ".parquet" else pd.read_csv(p)
        except Exception:
            return {}
    if "gene_name" not in df.columns or "gene_id" not in df.columns:
        return {}
    sub = df.loc[df["gene_name"].astype(str).isin(seed_symbols)]
    for gn, gid in zip(sub["gene_name"].astype(str), sub["gene_id"].astype(str)):
        gn, gid = gn.strip(), gid.strip()
        if not gn or gn == "nan" or not gid or gid == "nan" or not gid.startswith("ENSG"):
            continue
        out[gid] = gn
        base = strip_ensembl_version(gid)
        if base and base != gid:
            out[base] = gn
    return out


@lru_cache(maxsize=8)
def _gene_maps_cached(max_lines_key: int) -> Tuple[Dict[str, str], Dict[str, FrozenSet[str]]]:
    max_lines: Optional[int] = None if max_lines_key < 0 else max_lines_key
    from ..config import (
        PATHS,
        PANEL_ALIAS_SEED_SYMBOLS,
        UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE,
        LEGACY_DATASET_SYMBOL_RENAMES,
    )

    hgnc_ucsc_legacy = build_rna_expression_symbol_mapping(
        PATHS.hgnc_alias_table,
        PANEL_ALIAS_SEED_SYMBOLS,
        UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE,
        LEGACY_DATASET_SYMBOL_RENAMES,
        max_lines=max_lines,
        verbose=False,
    )
    manual = _manual_alias_to_canonical()
    flat = merge_symbol_maps(hgnc_ucsc_legacy, manual)
    ens = _ensembl_gene_id_flat_map(Path(PATHS.gencode_gtf_pq), PANEL_ALIAS_SEED_SYMBOLS)
    flat = merge_symbol_maps(flat, ens)
    inv = _invert_flat(flat)
    return flat, inv


def get_gene_flat_alias_map(max_lines: Optional[int] = None) -> Dict[str, str]:
    """All non-canonical tokens that should fold into a panel gene symbol."""
    return dict(_gene_maps_cached(_max_lines_key(max_lines))[0])


def get_gene_canonical_to_aliases(max_lines: Optional[int] = None) -> Dict[str, FrozenSet[str]]:
    """Canonical symbol -> every known synonym (including self and Ensembl ids)."""
    return dict(_gene_maps_cached(_max_lines_key(max_lines))[1])


# --- mature miRNA arm ids (TargetScan / miRTarBase / literature spellings) ---

MANUAL_MIRNA_CANONICAL_TO_ALIASES: Dict[str, Tuple[str, ...]] = {
    "hsa-miR-155-5p": ("miR-155-5p", "MIR155", "hsa-mir-155-5p"),
    "hsa-miR-21-5p": ("miR-21-5p", "hsa-mir-21-5p"),
    "hsa-let-7a-5p": ("let-7a-5p", "LET-7A-5P", "hsa-let-7a", "hsa-mir-let-7a-5p"),
}

# Hypothesis-tier miRNAs from ``research_plan/01_gene_panel_extended.md`` — **arm-level**
# canonical ids (aligned with ``mirna_mature_loci.csv`` / miRBase) plus common aliases.
TIER_MIRNA_ARM_CANONICAL_TO_ALIASES: Dict[str, Tuple[str, ...]] = {
    "hsa-miR-148a-3p": ("miR-148a-3p", "MIR148A", "hsa-mir-148a-3p"),
    "hsa-miR-152-3p": ("miR-152-3p", "MIR152", "hsa-mir-152-3p"),
    "hsa-miR-9-5p": ("miR-9-5p", "MIR9", "MIR9-1", "MIR9-2", "MIR9-3", "hsa-mir-9-5p"),
    "hsa-miR-125a-5p": ("miR-125a-5p", "MIR125A", "hsa-mir-125a-5p"),
    "hsa-miR-125b-5p": ("miR-125b-5p", "MIR125B1", "MIR125B2", "hsa-mir-125b-5p"),
    "hsa-miR-27a-3p": ("miR-27a-3p", "MIR27A", "hsa-mir-27a-3p"),
    "hsa-miR-34a-5p": ("miR-34a-5p", "MIR34A", "hsa-mir-34a-5p"),
    "hsa-miR-346-5p": (
        "miR-346",
        "miR-346-5p",
        "MIR346",
        "hsa-mir-346",
        "hsa-miR-346",
        "hsa-mir-346-5p",
    ),
    "hsa-miR-200c-3p": ("miR-200c-3p", "MIR200C", "hsa-mir-200c-3p"),
    "hsa-miR-146a-5p": (
        "miR-146a",
        "miR-146a-5p",
        "MIR146A",
        "hsa-mir-146a",
        "hsa-miR-146a",
        "hsa-mir-146a-5p",
    ),
}


def _mirna_alias_token_flat() -> Dict[str, str]:
    m: Dict[str, str] = {}
    for src in (MANUAL_MIRNA_CANONICAL_TO_ALIASES, TIER_MIRNA_ARM_CANONICAL_TO_ALIASES):
        for canon, aliases in src.items():
            c = canon.strip()
            for a in aliases:
                if a and str(a).strip() != c:
                    m[str(a).strip()] = c
    return m


@lru_cache(maxsize=1)
def get_mirna_flat_alias_map() -> Dict[str, str]:
    """Alias / alternate spellings -> canonical mature-arm miRNA id (hsa-miR-…)."""
    return dict(_mirna_alias_token_flat())


@lru_cache(maxsize=1)
def get_mirna_mimat_to_arm_id_map() -> Dict[str, str]:
    """
    MIMAT accession → miRBase mature-arm id (``hsa-miR-…``), from ``mirna_mature_loci.csv``.

    **Level:** mature transcript / arm (not pre-miRNA hairpin locus).
    """
    from ..config import PATHS

    p = Path(getattr(PATHS, "mirna_mature_loci_csv", ""))
    if not p.is_file():
        return {}
    try:
        df = pd.read_csv(p, usecols=["mature_accession", "mirbase_mature_id"], low_memory=False)
    except (ValueError, KeyError, OSError):
        try:
            df = pd.read_csv(p, low_memory=False)
        except Exception:
            return {}
    if "mature_accession" not in df.columns or "mirbase_mature_id" not in df.columns:
        return {}
    out: Dict[str, str] = {}
    for raw_acc, raw_arm in zip(df["mature_accession"], df["mirbase_mature_id"]):
        if raw_acc is None or (isinstance(raw_acc, float) and pd.isna(raw_acc)):
            continue
        if raw_arm is None or (isinstance(raw_arm, float) and pd.isna(raw_arm)):
            continue
        acc, arm = str(raw_acc).strip(), str(raw_arm).strip()
        if acc.startswith("MIMAT") and arm and arm.lower() not in ("nan", "none", ""):
            out[acc] = _harmonize_hsa_mir_token(arm)
    return out


def _harmonize_hsa_mir_token(s: str) -> str:
    """Normalize common ``hsa-mir-`` → ``hsa-miR-`` spelling for mature ids."""
    if not s:
        return s
    return re.sub(r"^hsa-mir-", "hsa-miR-", s, flags=re.IGNORECASE)


def normalize_mirna_symbol(
    token: str,
    mapping: Optional[Mapping[str, str]] = None,
    *,
    mimat_map: Optional[Mapping[str, str]] = None,
) -> str:
    """
    Normalize one miRNA label to a **mature-arm** canonical id.

    1. If the token looks like a **MIMAT** accession (Xena row ids), map via
       ``mirna_mature_loci.csv``.
    2. Harmonize ``hsa-mir-`` / ``hsa-miR-``.
    3. Apply alias → canonical map (short names, MIR symbols, tier spellings).
    """
    if token is None:
        return ""
    s = str(token).strip()
    if not s:
        return ""
    mm = mimat_map if mimat_map is not None else get_mirna_mimat_to_arm_id_map()
    if s.startswith("MIMAT"):
        arm = mm.get(s)
        if arm:
            s = arm
    s = _harmonize_hsa_mir_token(s)
    m = mapping if mapping is not None else get_mirna_flat_alias_map()
    if s in m:
        return m[s]
    low = s.lower()
    for k, v in m.items():
        if k.lower() == low:
            return v
    return s


def clear_gene_alias_caches() -> None:
    """Test helper: drop lru_cache contents after env / path changes."""
    _gene_maps_cached.cache_clear()
    get_mirna_flat_alias_map.cache_clear()
    get_mirna_mimat_to_arm_id_map.cache_clear()
