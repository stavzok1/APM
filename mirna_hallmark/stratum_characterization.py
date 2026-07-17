"""Characterize each Hallmark set across strata by target-gene dosage and miRNA expression.

For every Hallmark set x stratum (PAM50 / TNBC / immune / stage) we summarize:

1. **Target-gene CNV burden** -- ASCAT3 copy number of the Hallmark's member
   genes (mean CN, % gain / % loss / % deep-deletion), answering "is this
   program's gene dosage skewed in this stratum?".

2. **Regulatory miRNA expression** -- expression (log2(RPM+1)) of the
   *high-evidence* miRNAs that target the Hallmark's genes, answering "how much
   repressive miRNA is present against this program in this stratum?".

Together these two views let the interaction module ask whether strata with high
regulatory-miRNA expression also show the expected target-program suppression.

Outputs (``output/stratum_characterization/``):
- ``hallmark_target_cnv_by_stratum.tsv``
- ``hallmark_mirna_expr_by_stratum.tsv``
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.hallmark_sets import HallmarkSets


def _gene_hallmark_membership(hs: HallmarkSets) -> pd.DataFrame:
    """Long (gene, hallmark_set) membership frame (a gene may appear in many sets)."""
    rows = [(g, s) for g, sets in hs.gene_to_sets.items() for s in sets]
    return pd.DataFrame(rows, columns=["gene", "hallmark_set"])


def target_cnv_by_stratum(
    hs: HallmarkSets,
    clinical: pd.DataFrame,
) -> pd.DataFrame:
    """Mean CN + alteration prevalence of each Hallmark's target genes per stratum."""
    cnv = D.load_cnv_target_genes(hs.universe)  # gene x participant (cached)
    long = (
        cnv.rename_axis("gene").reset_index().melt(
            id_vars="gene", var_name="participant", value_name="copy_number"
        )
    ).dropna(subset=["copy_number"])
    long["cn_state"] = D.classify_cn(long["copy_number"])
    membership = _gene_hallmark_membership(hs)
    long = long.merge(membership, on="gene", how="inner")
    long = long.merge(clinical, on="participant", how="left")

    parts = []
    for stratum_col, layer in C.STRATUM_SPECS:
        if stratum_col not in clinical.columns:
            continue
        sub = long.dropna(subset=[stratum_col])
        g = sub.groupby(["hallmark_set", stratum_col], observed=True)
        agg = g.agg(
            n_gene_sample=("copy_number", "size"),
            n_genes=("gene", "nunique"),
            n_participants=("participant", "nunique"),
            mean_copy_number=("copy_number", "mean"),
        ).reset_index()
        st = sub.assign(
            _gain=sub["cn_state"].isin(["gain", "amp"]),
            _loss=sub["cn_state"].eq("loss"),
            _deep=sub["cn_state"].eq("deep_del"),
        ).groupby(["hallmark_set", stratum_col], observed=True)[["_gain", "_loss", "_deep"]].mean().reset_index()
        agg = agg.merge(st, on=["hallmark_set", stratum_col])
        agg = agg.rename(columns={stratum_col: "stratum"})
        agg.insert(0, "stratification_layer", layer)
        agg["mean_copy_number"] = agg["mean_copy_number"].round(4)
        agg["pct_gain"] = (100 * agg.pop("_gain")).round(2)
        agg["pct_loss"] = (100 * agg.pop("_loss")).round(2)
        agg["pct_deep_del"] = (100 * agg.pop("_deep")).round(2)
        parts.append(agg)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def mirna_expr_by_stratum(
    edges: pd.DataFrame,
    clinical: pd.DataFrame,
    *,
    high_evidence_only: bool = True,
) -> pd.DataFrame:
    """Expression of Hallmark-targeting miRNAs per stratum."""
    e = edges.loc[edges["high_evidence"]] if high_evidence_only else edges
    he = e[["miRNA", "hallmark_set"]].drop_duplicates()

    expr = D.load_mirna_arms()  # arm x participant log2(RPM+1)
    he = he.loc[he["miRNA"].isin(expr.index)]
    expr_long = (
        expr.loc[he["miRNA"].unique()].rename_axis("miRNA").reset_index().melt(
            id_vars="miRNA", var_name="participant", value_name="log2rpm"
        )
    ).dropna(subset=["log2rpm"])
    long = he.merge(expr_long, on="miRNA", how="inner").merge(clinical, on="participant", how="left")

    parts = []
    for stratum_col, layer in C.STRATUM_SPECS:
        if stratum_col not in clinical.columns:
            continue
        sub = long.dropna(subset=[stratum_col])
        agg = sub.groupby(["hallmark_set", stratum_col], observed=True).agg(
            n_mirnas=("miRNA", "nunique"),
            n_participants=("participant", "nunique"),
            mean_log2rpm=("log2rpm", "mean"),
            median_log2rpm=("log2rpm", "median"),
        ).reset_index().rename(columns={stratum_col: "stratum"})
        agg.insert(0, "stratification_layer", layer)
        agg["mean_log2rpm"] = agg["mean_log2rpm"].round(4)
        agg["median_log2rpm"] = agg["median_log2rpm"].round(4)
        parts.append(agg)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def run(*, out_dir: Path = C.STRATUM_DIR, high_evidence_only: bool = True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    edges = D.load_hallmark_edges()

    print("[stratum] target-gene CNV burden per Hallmark x stratum ...")
    cnv_tab = target_cnv_by_stratum(hs, clinical)
    cnv_tab.to_csv(out_dir / "hallmark_target_cnv_by_stratum.tsv", sep="\t", index=False)
    print(f"[stratum]   rows: {len(cnv_tab):,}")

    print("[stratum] Hallmark-targeting miRNA expression per stratum ...")
    expr_tab = mirna_expr_by_stratum(edges, clinical, high_evidence_only=high_evidence_only)
    expr_tab.to_csv(out_dir / "hallmark_mirna_expr_by_stratum.tsv", sep="\t", index=False)
    print(f"[stratum]   rows: {len(expr_tab):,}")

    manifest = {
        "module": "mirna_hallmark.stratum_characterization",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_hallmark_sets": hs.n_sets,
        "n_universe_genes": hs.n_genes,
        "strata": [layer for _, layer in C.STRATUM_SPECS],
        "mirna_filter": "high_evidence only" if high_evidence_only else "all edges",
        "target_cnv": "ASCAT3 integer CN of Hallmark member genes; mean CN + % gain/loss/deep_del",
        "mirna_expr": "log2(RPM+1) of Hallmark-targeting miRNAs (arm-level)",
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[stratum] wrote outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all-edges", action="store_true", help="Use all edges (not just high-evidence)")
    ap.add_argument("--out-dir", type=Path, default=C.STRATUM_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, high_evidence_only=not args.all_edges)


if __name__ == "__main__":
    main()
