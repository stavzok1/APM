"""ENCORI collapsed pair table and M0′ edge builder for dual-spine comparison."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Set

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.evidence_scoring import (
    compute_enc_depth,
    compute_encori_edge_score,
    load_encori_mrna_table,
)
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna
from mirna_hallmark.pressure_engine import filter_edges_by_abundance
from mirna_hallmark.robustness_checks import HUB_ROUTES


@dataclass(frozen=True)
class EncoriM0Spec:
    clip_min: int = 2
    top_k_program: int = 25
    top_k_hub: int = 50
    abundance_floor: Optional[float] = None
    score_mode: str = "depth"  # depth | reliability | two_channel


def _top_k(gene: str, hub_genes: Set[str], spec: EncoriM0Spec) -> int:
    return spec.top_k_hub if gene in hub_genes else spec.top_k_program


def collapse_encori_pairs(
    enc: pd.DataFrame,
    *,
    genes: Optional[Sequence[str]] = None,
    score_mode: str = "depth",
) -> pd.DataFrame:
    """One row per (miRNA, gene) with max clipExpNum and computed enc score."""
    if enc.empty:
        return pd.DataFrame(columns=["miRNA", "gene", "enc_depth", "clipExpNum"])

    df = enc.copy()
    df["gene"] = df["geneName"].astype(str).str.strip()
    df["miRNA"] = df["miRNAname"].astype(str).str.strip()
    df = df.loc[df["miRNA"].str.startswith("hsa-", na=False)]
    if genes is not None:
        df = df.loc[df["gene"].isin(set(genes))]

    for c in ["clipExpNum", "degraExpNum", "pancancerNum", "TargetScan", "miRanda", "PITA"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    df = df.sort_values(
        ["miRNA", "gene", "clipExpNum", "degraExpNum", "pancancerNum"],
        ascending=[True, True, False, False, False],
    )
    out = df.drop_duplicates(["miRNA", "gene"], keep="first").copy()
    out["enc_depth"] = compute_encori_edge_score(out, score_mode)
    out["rank_score"] = compute_enc_depth(out) if score_mode == "two_channel" else out["enc_depth"]
    return out[["miRNA", "gene", "enc_depth", "rank_score", "clipExpNum"] + [c for c in ("degraExpNum", "pancancerNum") if c in out.columns]]


def load_collapsed_encori_pairs(
    *,
    parquet_path: Path | None = None,
    genes: Optional[Sequence[str]] = None,
    score_mode: str = "depth",
) -> pd.DataFrame:
    enc = load_encori_mrna_table(parquet_path)
    return collapse_encori_pairs(enc, genes=genes, score_mode=score_mode)


def enc_depth_lookup(
    *,
    parquet_path: Path | None = None,
    genes: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """(miRNA, gene) -> enc_depth for S1 joins."""
    pairs = load_collapsed_encori_pairs(parquet_path=parquet_path, genes=genes)
    if pairs.empty:
        return pairs
    return pairs[["miRNA", "gene", "enc_depth"]].drop_duplicates()


def build_encori_m0_edges(
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
    *,
    spec: EncoriM0Spec,
    parquet_path: Path | None = None,
    hub_genes: Optional[Set[str]] = None,
) -> pd.DataFrame:
    """ENCORI M0′: clip gate, abundance floor, per-gene top-K by rank_score."""
    pairs = load_collapsed_encori_pairs(
        parquet_path=parquet_path, genes=genes, score_mode=spec.score_mode,
    )
    if pairs.empty:
        return pd.DataFrame(columns=["miRNA", "gene", "evidence_score"])

    sub = pairs.loc[pd.to_numeric(pairs["clipExpNum"], errors="coerce").fillna(0) >= spec.clip_min].copy()
    if sub.empty:
        return pd.DataFrame(columns=["miRNA", "gene", "evidence_score"])

    floor = C.PRESSURE_ABUNDANCE_FLOOR if spec.abundance_floor is None else spec.abundance_floor
    sub = filter_edges_by_abundance(sub, mirna_df, floor)
    if sub.empty:
        return pd.DataFrame(columns=["miRNA", "gene", "evidence_score"])

    hubs = hub_genes if hub_genes is not None else set(HUB_ROUTES.keys())
    rank_col = "rank_score" if "rank_score" in sub.columns else "enc_depth"
    sub = sub.sort_values(["gene", rank_col], ascending=[True, False])
    capped = []
    for gene, grp in sub.groupby("gene", sort=False):
        capped.append(grp.head(_top_k(str(gene), hubs, spec)))
    sub = pd.concat(capped, ignore_index=True) if capped else sub.iloc[:0]

    edges = sub[["miRNA", "gene", "enc_depth"]].rename(columns={"enc_depth": "evidence_score"})
    resolved, _audit = resolve_edges_mirna(edges, mirna_df, weight_col="evidence_score")
    return resolved


def write_encori_pair_table(
    *,
    out_path: Path | None = None,
    genes: Optional[Sequence[str]] = None,
) -> Path:
    out = Path(out_path or C.OUTPUT_ROOT / "encori" / "encori_pair_table.tsv.gz")
    pairs = load_collapsed_encori_pairs(genes=genes)
    out.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(out, sep="\t", index=False, compression="gzip")
    return out
