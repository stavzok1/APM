"""
HGNC-style alias resolution for gene symbols in wide matrices (RNA, etc.).

The file ``annotations/Alias_v5.22.xls`` is a tab-header, comma-separated-body
table listing approved ``gene_symbol`` and alternate ``content`` strings.
We stream it once and build old_symbol -> approved_symbol for tokens that
look like gene symbols, restricted to a caller-provided approved set.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, MutableMapping, Optional

_TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{1,39}$")


def strip_ensembl_version(gid: Any) -> str:
    """
    Strip Ensembl version suffix from gene / transcript ids (``ENSG…``.12 → ``ENSG…``).

    Used for matching external tables that may omit or include the version number.
    """
    s = str(gid).strip()
    if not s or s == "nan":
        return ""
    if s.startswith("ENS") and "." in s:
        return s.split(".", 1)[0]
    return s

_ALIAS_TYPES = frozenset({"Alias", "PreviousIdentifier"})


def _is_alias_token(approved: str, content: str) -> bool:
    if not content or content == approved:
        return False
    if len(content) > 48:
        return False
    if any(ch.isspace() for ch in content):
        return False
    if content.startswith("ENSG") or content.startswith("ENST"):
        return False
    if content.isdigit():
        return False
    return bool(_TOKEN.match(content))


def load_hgnc_alias_renames(
    alias_path: Path,
    approved_symbols: frozenset[str],
    *,
    max_lines: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, str]:
    """
    Map alternate symbols -> approved HGNC/GENCODE gene_symbol for rows
    whose approved symbol appears in ``approved_symbols``.

    Later rows overwrite earlier ones for the same alias (NCBI/HGNC
    duplicates are common).
    """
    path = Path(alias_path)
    if not path.is_file():
        return {}

    out: Dict[str, str] = {}
    n = 0
    with path.open("r", encoding="latin-1", errors="replace", newline="") as handle:
        header = handle.readline()
        if not header:
            return {}
        # Header is tab-separated; body lines are comma-separated CSV rows.
        for line in handle:
            n += 1
            if max_lines is not None and n > max_lines:
                break
            line = line.rstrip("\r\n")
            if not line:
                continue
            try:
                row = next(csv.reader([line]))
            except StopIteration:
                continue
            if len(row) < 4:
                continue
            approved, content, _acc, typ = row[0], row[1], row[2], row[3]
            if typ not in _ALIAS_TYPES:
                continue
            if approved not in approved_symbols:
                continue
            if not _is_alias_token(approved, content):
                continue
            out[content] = approved

    if verbose:
        print(f"  HGNC alias map: {len(out)} entries (scanned {n} data lines from {path.name})")
    return out


def merge_symbol_maps(
    base: Mapping[str, str],
    *extra: Mapping[str, str],
) -> Dict[str, str]:
    """Later mappings override earlier keys."""
    merged: Dict[str, str] = dict(base)
    for part in extra:
        merged.update(part)
    return merged


def build_rna_expression_symbol_mapping(
    alias_path: Path,
    approved_symbols: FrozenSet[str],
    ucsc_overrides: Mapping[str, str],
    legacy_renames: Mapping[str, str],
    *,
    max_lines: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, str]:
    """
    Full replacement dict for ``normalize_expression_mat``:
    HGNC aliases -> approved, then UCSC manual fixes, then legacy symbols
    seen in some TCGA/Xena builds.
    """
    hgnc = load_hgnc_alias_renames(
        alias_path, approved_symbols, max_lines=max_lines, verbose=verbose
    )
    return merge_symbol_maps(hgnc, ucsc_overrides, legacy_renames)


def symbols_missing_from_index(
    wanted: Iterable[str],
    index_symbols: Iterable[str],
) -> list[str]:
    """Return wanted symbols absent from the given index (e.g. RNA matrix rows)."""
    idx = set(index_symbols)
    return [g for g in wanted if g not in idx]
