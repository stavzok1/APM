"""Pathway/program enrichment of the discovery output — the biological-validation layer (user 2026-07-05).
Reuses the precursor's `stats.hypergeom_enrichment` (one-sided hypergeometric) + `stats.bh_fdr`.

Two questions:
  1. discovery_enrichment  — are the (robust / fully-novel) discovery TARGET genes over-represented in any
     Hallmark program vs the HE-gene universe?
  2. mirna_target_enrichment — for a hub miRNA (family), are its robust-discovery target genes concentrated in
     a particular program? (turns the program-network hubs into an enrichment-tested claim)

CLI: `python -m mirna_hallmark.learned.analyses.enrichment`
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.learned import families as FAM

_DISC = "mirna_hallmark/output/learned/discoveries.tsv"


def _hallmark_sets() -> dict:
    from mirna_hallmark.hallmark_sets import HallmarkSets
    hs = HallmarkSets.load()
    return hs.sets if hasattr(hs, "sets") else {}


def _universe() -> set:
    return set(D.high_evidence_edges()["gene"].dropna().astype(str))


def _enrich_geneset(hits: set, universe: set, sets: dict, *, min_prog: int = 5, min_hit: int = 2) -> pd.DataFrame:
    hits = hits & universe
    rows = []
    for prog, genes in sets.items():
        genes = set(genes) & universe
        k = len(hits & genes)
        if k < min_hit or len(genes) < min_prog:
            continue
        fold, p = S.hypergeom_enrichment(k, len(genes), len(hits), len(universe))
        rows.append({"program": prog.replace("HALLMARK_", ""), "n_prog": len(genes),
                     "n_hit": k, "fold": round(fold, 2), "p": p})
    df = pd.DataFrame(rows)
    if len(df):
        df["q"] = S.bh_fdr(df["p"].values)
        df = df.sort_values("p")
    return df


def discovery_enrichment(*, robust_only: bool = True) -> pd.DataFrame:
    """Are the discovery target genes enriched in Hallmark programs? (robust + fully-novel sets)."""
    disc = pd.read_csv(_DISC, sep="\t")
    if robust_only and "robust" in disc:
        disc = disc[disc["robust"]]
    uni, sets = _universe(), _hallmark_sets()
    out = []
    for label, sub in [("robust", disc), ("novel", disc[disc["sub_he_evidence"] == 0])]:
        e = _enrich_geneset(set(sub["gene"].astype(str)), uni, sets)
        if len(e):
            e.insert(0, "set", label)
            e.insert(1, "n_targets", sub["gene"].nunique())
            out.append(e)
    df = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    if len(df):
        with pd.option_context("display.width", 160):
            print("=== discovery TARGET-gene enrichment across Hallmark programs (BH per set) ===")
            print(df[df["q"] < 0.1].to_string(index=False))
    return df


def mirna_target_enrichment(arms=("hsa-miR-106b-5p", "hsa-miR-29a-3p", "hsa-miR-17-5p", "hsa-let-7b-5p"),
                            *, robust_only: bool = True) -> pd.DataFrame:
    """For each hub arm's seed family: are its robust-discovery target genes concentrated in a program?"""
    disc = pd.read_csv(_DISC, sep="\t")
    if robust_only and "robust" in disc:
        disc = disc[disc["robust"]]
    fam_of = FAM.family_of(sorted(disc["arm"].unique()))
    disc["family"] = disc["arm"].map(fam_of)
    uni, sets = _universe(), _hallmark_sets()
    rows = []
    for arm in arms:
        fam = FAM.family_of(pd.Index([arm])).iloc[0]
        tg = set(disc.loc[disc["family"] == fam, "gene"].astype(str))
        if len(tg) < 3:
            continue
        e = _enrich_geneset(tg, uni, sets)
        if len(e):
            top = e.iloc[0]
            rows.append({"arm": arm.replace("hsa-", ""), "family": fam, "n_targets": len(tg),
                         "top_program": top["program"], "n_hit": int(top["n_hit"]),
                         "fold": top["fold"], "p": top["p"], "q": round(float(top["q"]), 4)})
    df = pd.DataFrame(rows)
    if len(df):
        with pd.option_context("display.width", 160):
            print("\n=== hub-miRNA (family) robust-discovery target enrichment — top program per family ===")
            print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    discovery_enrichment()
    mirna_target_enrichment()
