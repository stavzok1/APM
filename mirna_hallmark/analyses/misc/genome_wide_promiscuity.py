"""Canonical GENOME-WIDE miRNA binding budget D(m) for promiscuity normalization.

`D(m)` = Σ over arm *m*'s **entire** targetome of the per-edge binding weight — an
**evidence-mass** denominator, NOT a target count. The budget-split weight
`W(m,g) = w_eff(m,g) / D(m)` is then the *fraction* of arm m's total binding budget
spent on gene g: a genome-wide promiscuity normalization. (Within-universe degree/mass
would massively overstate each gene's share — that was the bug in the first force build;
see `[[aggregate-force-vs-abundance-design]]`.)

Two currencies, each matched to its numerator `w_eff`:
  - ``validated``  : w = log1p(evidence_score)  — miRTarBase functional MTI, genome-wide
  - ``targetscan`` : w = |weighted context++|   — TargetScan hsa predictions, genome-wide

Keyed by **expressible arm** (resolved onto the miRNA matrix index, so it joins the
analysis edges). Cached to ``output/matrices/gw_budget_{source}.tsv``.

Run standalone to (re)build the caches:
  ``.venv/bin/python3 -m mirna_hallmark.analyses.misc.genome_wide_promiscuity``
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna

TARGETSCAN_CONTEXT = C.REPO_ROOT / "data" / "miRNA" / "Predicted_Targets_Context_Scores.default_predictions.txt"
_CACHE_DIR = C.OUTPUT_ROOT / "matrices"


def _resolve_arm_weights(frame: pd.DataFrame) -> pd.Series:
    """frame[miRNA,gene,w] (mature names) -> Series arm -> Σ w over the genome-wide targetome.

    Resolves raw mature names onto expressible arms first, so the budget keys match the
    analysis edges; weights summed per resolved arm."""
    mirna_df = D.load_mirna_arms()
    f = frame.dropna(subset=["miRNA", "gene", "w"]).copy()
    f = f.groupby(["miRNA", "gene"], as_index=False)["w"].sum()          # one weight per (mature,gene)
    f = f.rename(columns={"w": "evidence_score"})                        # resolver's default weight_col
    resolved, _audit = resolve_edges_mirna(f, mirna_df)
    return resolved.groupby("miRNA")["evidence_score"].sum().astype(float)


def _gw_validated_budget() -> pd.Series:
    """Genome-wide Σ log1p(evidence_score) per arm, with evidence_score on the **CANONICAL spine
    construction** — replicates `load_mirtar_edges` genome-wide: `confidence_logclass` scorer
    (`config.PRESSURE_EVIDENCE_SCORER`), score ≥ `PRESSURE_MIN_EVIDENCE`, + ENCORI α-boost
    (`PRESSURE_ENCORI_ALPHA`). So the budget currency == the edge w_eff currency."""
    from pipeline.genes.mirtarbase import load_mirtarbase
    from mirna_hallmark.build_edges import compute_interaction_summary_fast
    from mirna_hallmark.evidence_scoring import apply_scorer, apply_encori_boost_to_edges
    from mirna_hallmark.encori_edges import enc_depth_lookup

    raw = pd.read_csv(C.MIRTARBASE_CSV, usecols=["Target Gene", "Species (miRNA)"], low_memory=False)
    all_genes = raw.loc[raw["Species (miRNA)"] == "hsa", "Target Gene"].dropna().astype(str).unique().tolist()
    df = load_mirtarbase(C.MIRTARBASE_CSV, C.MIR_FAMILY_INFO, gene_panel=all_genes)
    summary = apply_scorer(compute_interaction_summary_fast(df), C.PRESSURE_EVIDENCE_SCORER)
    summary = summary[pd.to_numeric(summary["evidence_score"], errors="coerce").fillna(0) >= C.PRESSURE_MIN_EVIDENCE]
    edges = summary[["miRNA", "gene", "evidence_score"]].copy()
    alpha = float(getattr(C, "PRESSURE_ENCORI_ALPHA", 0.0))
    if alpha > 0:
        edges = apply_encori_boost_to_edges(edges, enc_depth_lookup(genes=all_genes), alpha=alpha)
    frame = pd.DataFrame({"miRNA": edges["miRNA"].astype(str), "gene": edges["gene"].astype(str),
                          "w": np.log1p(pd.to_numeric(edges["evidence_score"], errors="coerce").fillna(0.0))})
    return _resolve_arm_weights(frame)


def _gw_targetscan_budget() -> pd.Series:
    """Genome-wide Σ |weighted context++| per arm (TargetScan hsa predictions, all genes)."""
    parts = []
    for chunk in pd.read_csv(TARGETSCAN_CONTEXT, sep="\t",
                             usecols=lambda c: c in ("Gene Symbol", "miRNA", "weighted context++ score"),
                             chunksize=2_000_000):
        chunk = chunk.loc[chunk["miRNA"].astype(str).str.startswith("hsa-")].copy()
        chunk["w"] = pd.to_numeric(chunk["weighted context++ score"], errors="coerce").abs()
        chunk = chunk.dropna(subset=["w"])
        if not chunk.empty:
            parts.append(chunk[["miRNA", "Gene Symbol", "w"]].rename(columns={"Gene Symbol": "gene"}))
    if not parts:
        return pd.Series(dtype=float)
    return _resolve_arm_weights(pd.concat(parts, ignore_index=True))


@lru_cache(maxsize=4)
def genome_wide_budget(source: str = "validated", force: bool = False) -> pd.Series:
    """D(m) = genome-wide binding budget per expressible arm. ``source`` ∈ {validated, targetscan}."""
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = _CACHE_DIR / f"gw_budget_{source}.tsv"
    if cache.exists() and not force:
        s = pd.read_csv(cache, sep="\t", index_col=0).iloc[:, 0]
        s.index = s.index.astype(str)
        return s.astype(float)
    if source == "validated":
        s = _gw_validated_budget()
    elif source == "targetscan":
        s = _gw_targetscan_budget()
    else:
        raise ValueError(f"unknown source: {source!r}")
    s.name = f"D_{source}"
    s.to_csv(cache, sep="\t")
    print(f"[gw-budget] {source}: {len(s):,} arms, "
          f"median D={s.median():.2f} IQR[{s.quantile(.25):.2f},{s.quantile(.75):.2f}] -> {cache}")
    return s


def main() -> None:
    for src in ("validated", "targetscan"):
        s = genome_wide_budget(src, force=True)
        top = s.sort_values(ascending=False).head(5)
        print(f"  {src} top-budget arms: " + ", ".join(f"{a}={v:.1f}" for a, v in top.items()))


if __name__ == "__main__":
    main()
