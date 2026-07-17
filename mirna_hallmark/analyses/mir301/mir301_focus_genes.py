"""Shared miR-301a focus-gene panel and evidence-class helpers.

Used by tissue-reference modules (``mir301a_target_depth``, ``mir301_family_depth``,
``mir301_family_network``). Kept separate from plasma-EV screening so TCGA tissue
analyses do not import the EV machinery.
"""

from __future__ import annotations

from typing import Set, Tuple

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D

MIR301A_ARM = "hsa-miR-301a-3p"

# Luciferase / exosome-validated or strongly cited BRCA-relevant targets.
PUBLISHED_BRCA: Tuple[str, ...] = (
    "BTG1",
    "PTEN",
    "MEOX2",
    "SMAD4",
    "TGFBR2",
    "CDKN1A",
    "FOSL2",
    "TIMP2",
)

# Prior TS-orphan cohort hits (top tier from capped screen).
PRIOR_ORPHAN_HITS: Tuple[str, ...] = (
    "RAI2",
    "MAF",
    "ENPP5",
    "SYBU",
    "ACVR1",
    "KLF7",
    "TSHZ1",
    "RPS6KA5",
    "SLAIN1",
    "PAN3",
    "CHST1",
)


def he_pairs() -> Set[Tuple[str, str]]:
    edges = D.load_hallmark_edges()
    he = edges.loc[edges["high_evidence"]].drop_duplicates(["miRNA", "gene"])
    return set(zip(he["miRNA"], he["gene"]))


def mirtar_low_genes(arm: str) -> Set[str]:
    if not C.MIRTAR_HALLMARK_SUMMARY.is_file():
        return set()
    m = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY)
    sub = m[m["miRNA"].astype(str).str.contains("301a", case=False, na=False)]
    he = he_pairs()
    out: Set[str] = set()
    for _, r in sub.iterrows():
        g, mir = r["gene"], r["miRNA"]
        if (mir, g) not in he:
            out.add(g)
    return out


def evidence_class(gene: str, arm: str, he: Set[Tuple[str, str]]) -> str:
    if gene in PUBLISHED_BRCA:
        return "published_brca"
    if (arm, gene) in he:
        return "mirtar_he"
    if gene in mirtar_low_genes(arm):
        return "mirtar_low"
    return "ts_orphan"
