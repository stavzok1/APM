"""Load and index the MSigDB Hallmark gene sets.

The annotation file ``annotations/GSEA/h.all.v2026.1.Hs.txt`` is a single-line
JSON object (despite the ``.txt`` extension): ``{set_name: {..., "geneSymbols":
[...]}}`` for the 50 ``H`` (Hallmark) collection sets.

Exposes:
- ``load_hallmark_sets()`` -> ``{set_name: [gene, ...]}``
- ``gene_to_hallmarks()`` -> ``{gene: [set_name, ...]}`` (genes may be in many sets)
- ``hallmark_universe()`` -> sorted unique gene list across all sets
- ``HallmarkSets`` convenience container
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from mirna_hallmark.config import HALLMARK_JSON


def _read_raw(path: Path) -> dict:
    text = path.read_text(encoding="utf-8").strip()
    return json.loads(text)


def load_hallmark_sets(path: Optional[Path] = None) -> Dict[str, List[str]]:
    """Return ``{HALLMARK_SET: [gene_symbol, ...]}`` sorted by set name."""
    raw = _read_raw(Path(path or HALLMARK_JSON))
    out: Dict[str, List[str]] = {}
    for set_name, payload in raw.items():
        genes = payload.get("geneSymbols") or []
        # de-dup, keep deterministic order
        out[set_name] = sorted(dict.fromkeys(str(g).strip() for g in genes if g))
    return dict(sorted(out.items()))


def gene_to_hallmarks(sets: Optional[Dict[str, List[str]]] = None) -> Dict[str, List[str]]:
    """Return ``{gene: [hallmark_set, ...]}`` (a gene can belong to many sets)."""
    sets = sets or load_hallmark_sets()
    mapping: Dict[str, List[str]] = {}
    for set_name, genes in sets.items():
        for g in genes:
            mapping.setdefault(g, []).append(set_name)
    return {g: sorted(v) for g, v in mapping.items()}


def hallmark_universe(sets: Optional[Dict[str, List[str]]] = None) -> List[str]:
    """Return the sorted union of all genes across the 50 Hallmark sets."""
    sets = sets or load_hallmark_sets()
    universe: set = set()
    for genes in sets.values():
        universe.update(genes)
    return sorted(universe)


@dataclass
class HallmarkSets:
    """Container bundling the three indexed views."""

    sets: Dict[str, List[str]]
    gene_to_sets: Dict[str, List[str]]
    universe: List[str]

    @classmethod
    def load(cls, path: Optional[Path] = None) -> "HallmarkSets":
        sets = load_hallmark_sets(path)
        return cls(
            sets=sets,
            gene_to_sets=gene_to_hallmarks(sets),
            universe=hallmark_universe(sets),
        )

    @property
    def n_sets(self) -> int:
        return len(self.sets)

    @property
    def n_genes(self) -> int:
        return len(self.universe)


if __name__ == "__main__":
    hs = HallmarkSets.load()
    print(f"[hallmark_sets] {hs.n_sets} sets, {hs.n_genes} unique genes")
    multi = sum(1 for g, s in hs.gene_to_sets.items() if len(s) > 1)
    print(f"[hallmark_sets] genes in >1 set: {multi}")
    sizes = {k: len(v) for k, v in hs.sets.items()}
    biggest = max(sizes, key=sizes.get)
    print(f"[hallmark_sets] largest set: {biggest} ({sizes[biggest]} genes)")
