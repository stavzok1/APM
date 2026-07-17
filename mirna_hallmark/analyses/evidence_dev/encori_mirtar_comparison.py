"""ENCORI vs miRTarBase pair overlap statistics (Hallmark universe)."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.encori_edges import collapse_encori_pairs
from mirna_hallmark.evidence_scoring import apply_scorer, build_m0_edges, load_encori_mrna_table
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.robustness_checks import HUB_ROUTES

OUT_DIR = C.OUTPUT_ROOT / "encori" / "mirtar_comparison"
HUB_GENES = tuple(HUB_ROUTES.keys())


def _pair_set(df: pd.DataFrame) -> set[tuple[str, str]]:
    if df.empty:
        return set()
    return set(zip(df["miRNA"].astype(str), df["gene"].astype(str)))


def hub_gene_comparison(
    mirtar_summary: pd.DataFrame,
    m0_edges: pd.DataFrame,
    enc_pairs: pd.DataFrame,
    *,
    hub_genes: Sequence[str] = HUB_GENES,
) -> pd.DataFrame:
    """Per-hub pair counts with separate overlap_m0 vs overlap_mirtar_all."""
    rows = []
    for gene in hub_genes:
        mt_g = mirtar_summary.loc[
            (mirtar_summary["gene"] == gene)
            & mirtar_summary["miRNA"].astype(str).str.startswith("hsa-", na=False)
        ]
        m0_g = m0_edges.loc[m0_edges["gene"] == gene]
        enc_g = enc_pairs.loc[enc_pairs["gene"] == gene] if not enc_pairs.empty else enc_pairs

        mt_set = _pair_set(mt_g.rename(columns={"miRNA": "miRNA", "gene": "gene"}))
        m0_set = _pair_set(m0_g)
        enc_set = _pair_set(enc_g)

        overlap_mirtar_all = len(mt_set & enc_set)
        overlap_m0 = len(m0_set & enc_set)

        clip_med = float("nan")
        if not enc_g.empty and "clipExpNum" in enc_g.columns:
            clip_med = float(pd.to_numeric(enc_g["clipExpNum"], errors="coerce").median())

        rows.append({
            "gene": gene,
            "mirtar_pairs": len(mt_set),
            "m0_pairs": len(m0_set),
            "encori_pairs": len(enc_set),
            "overlap_m0": overlap_m0,
            "overlap_mirtar_all": overlap_mirtar_all,
            "encori_only": len(enc_set - m0_set),
            "m0_only": len(m0_set - enc_set),
            "median_clipExpNum": clip_med,
        })
    return pd.DataFrame(rows)


def run(*, out_dir: Path | None = None) -> Path:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)

    hs = HallmarkSets.load()
    mirna_df = D.load_mirna_arms()
    mirtar_summary = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, low_memory=False)
    mirtar_summary = mirtar_summary.loc[mirtar_summary["gene"].isin(set(hs.universe))]
    mirtar_summary = mirtar_summary.loc[
        mirtar_summary["miRNA"].astype(str).str.startswith("hsa-", na=False)
    ]

    m0_edges = build_m0_edges(hs.universe, mirna_df, scorer_id="tiered_permissive")
    enc_raw = load_encori_mrna_table()
    enc_pairs = collapse_encori_pairs(enc_raw, genes=hs.universe)

    mt_set = _pair_set(mirtar_summary[["miRNA", "gene"]].drop_duplicates())
    m0_set = _pair_set(m0_edges)
    enc_set = _pair_set(enc_pairs)

    overlap_mirtar = len(mt_set & enc_set)
    overlap_m0 = len(m0_set & enc_set)

    summary = {
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "hallmark_genes": len(hs.universe),
        "mirtar_pairs": len(mt_set),
        "m0_pairs": len(m0_set),
        "encori_pairs": len(enc_set),
        "overlap_mirtar_all": overlap_mirtar,
        "overlap_m0": overlap_m0,
        "encori_only_vs_m0": len(enc_set - m0_set),
        "m0_only_vs_encori": len(m0_set - enc_set),
        "jaccard_encori_mirtar_all": overlap_mirtar / len(mt_set | enc_set) if (mt_set | enc_set) else 0.0,
        "jaccard_encori_m0": overlap_m0 / len(m0_set | enc_set) if (m0_set | enc_set) else 0.0,
        "frac_m0_with_encori": overlap_m0 / len(m0_set) if m0_set else 0.0,
    }
    (out_root / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8",
    )

    hub_tbl = hub_gene_comparison(mirtar_summary, m0_edges, enc_pairs)
    hub_tbl.to_csv(out_root / "hub_gene_comparison.tsv", sep="\t", index=False)

    per_gene = []
    for gene in sorted(hs.universe):
        mt_g = mirtar_summary.loc[mirtar_summary["gene"] == gene]
        m0_g = m0_edges.loc[m0_edges["gene"] == gene]
        enc_g = enc_pairs.loc[enc_pairs["gene"] == gene] if not enc_pairs.empty else enc_pairs
        mt_set_g = _pair_set(mt_g[["miRNA", "gene"]].drop_duplicates())
        m0_set_g = _pair_set(m0_g)
        enc_set_g = _pair_set(enc_g)
        per_gene.append({
            "gene": gene,
            "mirtar_pairs": len(mt_set_g),
            "m0_pairs": len(m0_set_g),
            "encori_pairs": len(enc_set_g),
            "overlap_m0": len(m0_set_g & enc_set_g),
            "overlap_mirtar_all": len(mt_set_g & enc_set_g),
            "encori_only_pairs": len(enc_set_g - m0_set_g),
            "mirtar_only_pairs": len(mt_set_g - enc_set_g),
        })
    pd.DataFrame(per_gene).to_csv(out_root / "per_gene_pair_counts.tsv", sep="\t", index=False)

    overlap_rows = []
    for m, g in sorted(m0_set & enc_set):
        overlap_rows.append({"miRNA": m, "gene": g})
    pd.DataFrame(overlap_rows).to_csv(out_root / "overlap_m0_pairs.tsv", sep="\t", index=False)

    print(f"[encori_mirtar_comparison] M0={len(m0_set):,} ENCORI={len(enc_set):,} overlap_m0={overlap_m0:,}")
    return out_root


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
