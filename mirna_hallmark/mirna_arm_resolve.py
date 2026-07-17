"""Map miRTarBase / TargetScan mature-arm names onto GDC expression matrix rows.

Default harmonization (Tier A + guarded Tier B):

**Tier A — always applied**
- exact matrix name
- isoform suffixes ``.1`` / ``.2`` / ``.3`` (same mature product; GDC MIMAT labels)
- star-suffix strip (``*``)
- loci-file map: mirBase mature ID / MIMAT accession → expressible row

**Tier B — guarded**
- precursor stem without ``-3p``/``-5p`` → append mature arm **only when exactly one**
  candidate exists in the expression matrix (e.g. ``544a`` → ``544a-3p``)

**Not applied by default** (different molecules / paralog risk):
- opposite mature arm (``-3p`` ↔ ``-5p``)
- letter slim (``151a`` → ``151``)
- precursor sibling collapse (map all arms from one pre-miRNA to one expressible row)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from mirna_hallmark import config as C

MatchKind = str

_PREC_STEM = re.compile(r"^hsa-miR-\d+[a-z]?$")


def mirna_arm_alias_candidates(arm: str) -> List[str]:
    """Tier-A ordered candidates: exact, isoform suffixes, star strip, isoform strip."""
    arm = str(arm).strip()
    if not arm or arm.lower() == "nan":
        return []

    out: List[str] = []

    def add(x: str) -> None:
        x = str(x).strip()
        if x and x not in out:
            out.append(x)

    add(arm)

    for sfx in (".1", ".2", ".3"):
        add(f"{arm}{sfx}")

    m = re.match(r"^(.*)-([35]p(\*)?)$", arm)
    if m:
        stem, tail = m.group(1), m.group(2)
        for sfx in (".1", ".2", ".3"):
            add(f"{stem}-{tail}{sfx}")

    if arm.endswith("*"):
        base = arm[:-1]
        add(base)
        for sfx in (".1", ".2", ".3"):
            add(f"{base}{sfx}")

    if re.search(r"\.[0-9]+$", arm):
        add(re.sub(r"\.[0-9]+$", "", arm))

    return out


def _precursor_append_hit(arm: str, expr_index: set[str]) -> Optional[str]:
    """Tier B: map stem → single mature arm when unambiguous."""
    if not _PREC_STEM.match(arm) or arm in expr_index:
        return None
    hits: List[str] = []
    for tail in ("-3p", "-5p", "-2-3p"):
        for sfx in ("", ".1", ".2", ".3"):
            cand = f"{arm}{tail}{sfx}" if sfx else f"{arm}{tail}"
            if cand in expr_index and cand not in hits:
                hits.append(cand)
    if len(hits) == 1:
        return hits[0]
    return None


def _resolve_one(query: str, expr_index: set[str]) -> Tuple[Optional[str], MatchKind]:
    for i, cand in enumerate(mirna_arm_alias_candidates(query)):
        if cand not in expr_index:
            continue
        if i == 0:
            return cand, "exact"
        if cand.endswith(".1") or cand.endswith(".2") or cand.endswith(".3"):
            return cand, "isoform_suffix"
        if query.endswith("*") and cand == query[:-1]:
            return cand, "star_strip"
        return cand, "alias"

    prec = _precursor_append_hit(query, expr_index)
    if prec:
        return prec, "precursor_append"
    return None, "missing"


def build_arm_resolver(
    expr_index: Sequence[str],
    *,
    loci_path: Optional[Path] = None,
) -> Dict[str, str]:
    """Pre-register mirbase_mature_id / MIMAT → expressible arm from loci file."""
    idx = {str(x) for x in expr_index}
    resolver: Dict[str, str] = {a: a for a in idx}

    loci_p = Path(loci_path or C.MIRNA_MATURE_LOCI)
    if not loci_p.exists():
        return resolver

    loci = pd.read_csv(loci_p, low_memory=False)
    for _, row in loci.iterrows():
        mb = str(row.get("mirbase_mature_id", "")).strip()
        if not mb or mb.lower() == "nan":
            continue
        hit, _ = _resolve_one(mb, idx)
        if not hit:
            continue
        resolver[mb] = hit
        if mb.endswith(".1"):
            base = mb[:-2]
            if base not in idx:
                resolver.setdefault(base, hit)
        base_no_iso = re.sub(r"\.[0-9]+$", "", mb)
        if base_no_iso != mb and base_no_iso not in idx:
            resolver.setdefault(base_no_iso, hit)
        mimat = str(row.get("mature_accession", "")).strip()
        if mimat and mimat.lower() != "nan":
            resolver.setdefault(mimat, hit)

    return resolver


def resolve_arm(
    arm: str,
    expr_index: Sequence[str],
    resolver: Optional[Dict[str, str]] = None,
    idx_set: Optional[set] = None,
) -> Tuple[str, MatchKind]:
    """Return (expressible_arm, match_kind); *missing* keeps the query name.

    Pass ``idx_set`` (a prebuilt ``{str(x) for x in expr_index}``) to avoid rebuilding
    the expression-index set on every call -- the hot path resolves one arm per edge.
    """
    arm = str(arm).strip()
    idx = idx_set if idx_set is not None else {str(x) for x in expr_index}
    if arm in idx:
        return arm, "exact"
    if resolver and arm in resolver and resolver[arm] in idx:
        kind = "exact" if resolver[arm] == arm else "loci_map"
        return resolver[arm], kind
    hit, kind = _resolve_one(arm, idx)
    if hit:
        return hit, kind
    return arm, "missing"


def resolve_edges_mirna(
    edges: pd.DataFrame,
    mirna_log2: pd.DataFrame,
    *,
    loci_path: Optional[Path] = None,
    weight_col: str = "evidence_score",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Rewrite ``miRNA`` to expressible arms; collapse duplicates; emit audit."""
    if edges.empty:
        return edges.copy(), pd.DataFrame(
            columns=["miRNA_query", "miRNA_resolved", "match", "n_edges", "n_genes"]
        )

    idx = mirna_log2.index
    resolver = build_arm_resolver(idx, loci_path=loci_path)
    e = edges.copy()
    e["miRNA_query"] = e["miRNA"].astype(str)
    # Resolve each UNIQUE query arm once and map back, instead of once per edge.
    # (Previously each edge rebuilt the ~2k expr-index set and re-resolved repeated
    # arm names -- ~913k redundant calls on the all-edge set.) Result is identical.
    idx_set = {str(x) for x in idx}
    res_map: Dict[str, str] = {}
    kind_map: Dict[str, str] = {}
    for q in e["miRNA_query"].unique():
        r, k = resolve_arm(q, idx, resolver, idx_set=idx_set)
        res_map[q] = r
        kind_map[q] = k
    e["miRNA"] = e["miRNA_query"].map(res_map)
    e["arm_match"] = e["miRNA_query"].map(kind_map)

    agg: Dict[str, str | type] = {
        "miRNA_query": "first",
        "arm_match": "first",
    }
    if weight_col in e.columns:
        agg[weight_col] = "max"
    for c in e.columns:
        if c not in ("miRNA", "gene", "miRNA_query", "arm_match", weight_col):
            agg[c] = "first"

    out = e.groupby(["miRNA", "gene"], as_index=False).agg(agg)

    audit = (
        e.groupby(["miRNA_query", "miRNA", "arm_match"], as_index=False)
        .agg(n_edges=("gene", "size"), n_genes=("gene", "nunique"))
        .rename(columns={"miRNA": "miRNA_resolved", "arm_match": "match"})
        .sort_values(["match", "n_edges"], ascending=[True, False])
    )
    return out, audit


def write_arm_resolve_audit(
    audit: pd.DataFrame,
    path: Optional[Path] = None,
) -> Path:
    """Persist edge arm resolution audit (missing + resolved aliases)."""
    out = Path(path or C.OUTPUT_ROOT / "logs" / "mirna_arm_resolve_audit.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out, sep="\t", index=False)
    return out
