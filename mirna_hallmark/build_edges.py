"""Regenerate miRTarBase edges for the Hallmark gene universe and annotate them.

The on-disk ``data/miRNA/mirtarbase/mirtar_interaction_summary.csv`` is filtered
to ``PIPELINE_GENE_PANEL``. This subproject targets the full ~4,400-gene Hallmark
universe, so we rebuild a Hallmark-scoped interaction summary and write it into
the subproject output tree.

Performance: the canonical ``pipeline.genes.mirtarbase.get_mirtarbase_targets``
builds five tables via per-group ``.apply`` (incl. large nested-JSON gene/miRNA/
family summaries) and is far too slow over ~4,400 genes. We instead reuse only
the fast, vectorized normalizer ``load_mirtarbase`` and the canonical constants,
then compute the ``(miRNA, gene)`` interaction summary with vectorized groupby
aggregations on boolean flag columns. The output schema matches the columns the
pressure engine (``analysis.expression.mirna_target_integration.load_mirtar_edges``)
reads, so it stays drop-in compatible. Use ``--use-canonical-builder`` to fall
back to the slow reference implementation for parity checks.

We then emit ``output/edges/mirna_hallmark_edges.tsv.gz`` with one row per
(miRNA, target gene, hallmark_set). A gene in multiple Hallmark sets yields one
row per set (its (miRNA, gene) evidence repeated in each Hallmark context);
cohort-level counts dedup on (miRNA, gene).

"High-evidence" (configurable in ``config.py``) means the edge has at least one
``Functional MTI`` study AND a low-throughput functional validation (luciferase
reporter or protein/western), captured by the cross-count columns
``n_reporter__functional_mti_studies`` / ``n_protein__functional_mti_studies``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets


def compute_interaction_summary_fast(df_normalized: pd.DataFrame) -> pd.DataFrame:
    """Vectorized (miRNA, gene) interaction summary from the normalized miRTarBase table.

    Reproduces the column schema of the canonical
    ``mirtar_interaction_summary.csv`` for the fields downstream code consumes:
    per-experiment-class and per-support-type study counts, the
    experiment x support cross-counts, ``n_studies`` and ``evidence_score``.

    No per-group ``.apply`` -- all aggregation is groupby(max)/groupby(sum) on
    boolean indicator columns, so it scales to the full Hallmark universe.
    """
    from pipeline.genes.mirtarbase import EXPERIMENT_CLASSES, SUPPORT_TYPES, _support_slug

    cols = ["miRNA_norm", "gene_norm", "miRNA_family", "support_type_norm",
            "experiment_classes", "study_ids"]
    s = df_normalized[cols].copy()

    # Explode to one row per (miRNA, gene, study_id) ...
    s = s.explode("study_ids").rename(columns={"study_ids": "study_id"})
    s["study_id"] = s["study_id"].astype(str).str.strip()
    s = s.loc[s["study_id"] != ""]
    # ... then per experiment class (NaN-safe: empty lists -> NaN exp_class -> no flag)
    s = s.explode("experiment_classes").rename(columns={"experiment_classes": "exp_class"})

    support_slugs = {st: _support_slug(st) for st in SUPPORT_TYPES}
    exp_flag_cols = [f"has_{ex}" for ex in EXPERIMENT_CLASSES]
    sup_flag_cols = [f"has_{support_slugs[st]}" for st in SUPPORT_TYPES]
    for ex in EXPERIMENT_CLASSES:
        s[f"has_{ex}"] = s["exp_class"].eq(ex)
    for st in SUPPORT_TYPES:
        s[f"has_{support_slugs[st]}"] = s["support_type_norm"].eq(st)

    # Collapse to (miRNA, gene, study): OR the flags across raw rows of a study.
    keys = ["miRNA_norm", "gene_norm", "study_id"]
    study = (
        s.groupby(keys, sort=False)[exp_flag_cols + sup_flag_cols].max().astype(bool)
    )
    # study-level experiment x support cross flags
    cross_cols = []
    for ex in EXPERIMENT_CLASSES:
        for st in SUPPORT_TYPES:
            col = f"has_{ex}__{support_slugs[st]}"
            study[col] = study[f"has_{ex}"] & study[f"has_{support_slugs[st]}"]
            cross_cols.append(col)

    study = study.reset_index()
    grp = study.groupby(["miRNA_norm", "gene_norm"], sort=False)
    count_src = exp_flag_cols + sup_flag_cols + cross_cols
    summary = grp[count_src].sum().astype(int)
    summary["n_studies"] = grp.size()

    # rename to canonical n_*_studies schema
    rename = {f"has_{ex}": f"n_{ex}_studies" for ex in EXPERIMENT_CLASSES}
    rename.update({f"has_{support_slugs[st]}": f"n_{support_slugs[st]}_studies" for st in SUPPORT_TYPES})
    rename.update({c: f"n_{c[len('has_'):]}_studies" for c in cross_cols})
    summary = summary.rename(columns=rename)

    summary["evidence_score"] = (
        3 * summary["n_reporter_studies"]
        + 3 * summary["n_binding_studies"]
        + 2 * summary["n_protein_studies"]
        + 1 * summary["n_rna_studies"]
        + 1 * summary["n_perturbation_studies"]
    )

    summary = summary.reset_index().rename(columns={"miRNA_norm": "miRNA", "gene_norm": "gene"})
    # miRNA family is per-miRNA; map it back from the normalized table
    fam_map = (
        df_normalized.dropna(subset=["miRNA_norm"])
        .drop_duplicates("miRNA_norm")
        .set_index("miRNA_norm")["miRNA_family"]
    )
    summary.insert(1, "miRNA_family", summary["miRNA"].map(fam_map))
    return summary


def regenerate_mirtarbase_summary(
    *,
    universe: list[str],
    out_dir: Path = C.MIRTARBASE_HALLMARK_DIR,
    force: bool = False,
    use_canonical_builder: bool = False,
) -> Path:
    """Build a Hallmark-scoped miRTarBase (miRNA, gene) summary; return its path."""
    summary_path = out_dir / "mirtar_interaction_summary.csv"
    if summary_path.exists() and not force:
        print(f"[build_edges] reusing existing summary: {summary_path}")
        return summary_path

    out_dir.mkdir(parents=True, exist_ok=True)

    if use_canonical_builder:
        from pipeline.genes.mirtarbase import get_mirtarbase_targets

        print(f"[build_edges] (canonical/slow) regenerating miRTarBase for {len(universe)} genes")
        get_mirtarbase_targets(
            mirtarbase_csv=C.MIRTARBASE_CSV,
            family_info_tsv=C.MIR_FAMILY_INFO,
            gene_panel=universe,
            output_dir=out_dir,
        )
    else:
        from pipeline.genes.mirtarbase import load_mirtarbase

        print(f"[build_edges] (vectorized) normalizing miRTarBase for {len(universe)} Hallmark genes")
        df = load_mirtarbase(
            mirtarbase_csv=C.MIRTARBASE_CSV,
            family_info_tsv=C.MIR_FAMILY_INFO,
            gene_panel=universe,
        )
        summary = compute_interaction_summary_fast(df)
        summary.to_csv(summary_path, index=False)
        print(f"[build_edges] wrote vectorized summary: {summary_path} ({len(summary):,} miRNA-gene pairs)")

    if not summary_path.exists():
        raise FileNotFoundError(f"miRTarBase summary not produced at {summary_path}")
    return summary_path


def _high_evidence_mask(df: pd.DataFrame) -> pd.Series:
    """Boolean Series: Functional MTI + low-throughput validation (or score fallback)."""
    n_func = pd.to_numeric(df.get("n_functional_mti_studies"), errors="coerce").fillna(0)

    lowthr_cols = ["n_reporter__functional_mti_studies", "n_protein__functional_mti_studies"]
    have_cross = all(c in df.columns for c in lowthr_cols)
    if have_cross and C.HIGH_EVIDENCE_REQUIRE_LOWTHROUGHPUT:
        lowthr = sum(
            pd.to_numeric(df[c], errors="coerce").fillna(0) for c in lowthr_cols
        )
        return (n_func >= C.HIGH_EVIDENCE_MIN_FUNCTIONAL_MTI) & (lowthr >= 1)
    # fallback: functional MTI count + evidence score floor
    score = pd.to_numeric(df.get("evidence_score"), errors="coerce").fillna(0)
    return (n_func >= C.HIGH_EVIDENCE_MIN_FUNCTIONAL_MTI) & (score >= C.HIGH_EVIDENCE_MIN_SCORE)


def build_hallmark_edges(
    *,
    summary_path: Optional[Path] = None,
    hs: Optional[HallmarkSets] = None,
    write: bool = True,
) -> pd.DataFrame:
    """Annotate the miRTarBase summary with Hallmark membership + high_evidence.

    Returns one row per (miRNA, gene, hallmark_set).
    """
    hs = hs or HallmarkSets.load()
    summary_path = Path(summary_path or C.MIRTAR_HALLMARK_SUMMARY)
    summary = pd.read_csv(summary_path, low_memory=False)
    summary = summary.loc[summary["miRNA"].astype(str).str.startswith("hsa-", na=False)].copy()

    keep_cols = [
        c
        for c in [
            "miRNA",
            "miRNA_family",
            "gene",
            "n_studies",
            "n_functional_mti_studies",
            "n_functional_mti_weak_studies",
            "n_reporter_studies",
            "n_protein_studies",
            "n_reporter__functional_mti_studies",
            "n_protein__functional_mti_studies",
            "evidence_score",
        ]
        if c in summary.columns
    ]
    summary = summary[keep_cols].copy()
    summary["high_evidence"] = _high_evidence_mask(summary)

    # explode to (miRNA, gene, hallmark_set)
    g2h = hs.gene_to_sets
    summary["hallmark_set"] = summary["gene"].map(g2h)
    edges = summary.explode("hallmark_set").dropna(subset=["hallmark_set"]).reset_index(drop=True)

    if write:
        C.EDGES_DIR.mkdir(parents=True, exist_ok=True)
        edges.to_csv(C.HALLMARK_EDGES_TSV, sep="\t", index=False, compression="gzip")
        _write_manifest(edges, summary_path)
        print(f"[build_edges] wrote {len(edges):,} (miRNA,gene,hallmark) rows -> {C.HALLMARK_EDGES_TSV}")
    return edges


def _write_manifest(edges: pd.DataFrame, summary_path: Path) -> None:
    manifest = {
        "module": "mirna_hallmark.build_edges",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "mirtarbase_raw": str(C.MIRTARBASE_CSV),
            "mirtar_hallmark_summary": str(summary_path),
            "hallmark_json": str(C.HALLMARK_JSON),
        },
        "high_evidence_definition": {
            "min_functional_mti_studies": C.HIGH_EVIDENCE_MIN_FUNCTIONAL_MTI,
            "require_lowthroughput": C.HIGH_EVIDENCE_REQUIRE_LOWTHROUGHPUT,
            "lowthroughput_columns": [
                "n_reporter__functional_mti_studies",
                "n_protein__functional_mti_studies",
            ],
            "fallback_min_evidence_score": C.HIGH_EVIDENCE_MIN_SCORE,
        },
        "counts": {
            "n_rows_miRNA_gene_hallmark": int(len(edges)),
            "n_unique_edges_miRNA_gene": int(edges[["miRNA", "gene"]].drop_duplicates().shape[0]),
            "n_unique_mirnas": int(edges["miRNA"].nunique()),
            "n_unique_target_genes": int(edges["gene"].nunique()),
            "n_high_evidence_rows": int(edges["high_evidence"].sum()),
            "n_hallmark_sets_with_edges": int(edges["hallmark_set"].nunique()),
        },
    }
    (C.EDGES_DIR / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def run(*, force: bool = False, use_canonical_builder: bool = False) -> pd.DataFrame:
    """Regenerate the Hallmark-scoped summary and write the annotated edge table."""
    C.ensure_output_dirs()
    hs = HallmarkSets.load()
    print(f"[build_edges] Hallmark universe: {hs.n_genes} genes across {hs.n_sets} sets")
    summary_path = regenerate_mirtarbase_summary(
        universe=hs.universe,
        force=force,
        use_canonical_builder=use_canonical_builder,
    )
    edges = build_hallmark_edges(summary_path=summary_path, hs=hs, write=True)
    n_he = int(edges["high_evidence"].sum())
    print(
        f"[build_edges] done: {edges['miRNA'].nunique():,} miRNAs -> "
        f"{edges['gene'].nunique():,} Hallmark genes; high-evidence rows: {n_he:,}"
    )
    return edges


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force-rebuild", action="store_true", help="Re-run miRTarBase even if summary exists")
    ap.add_argument(
        "--use-canonical-builder",
        action="store_true",
        help="Use the slow reference pipeline.genes.mirtarbase.get_mirtarbase_targets (parity check)",
    )
    args = ap.parse_args()
    run(force=args.force_rebuild, use_canonical_builder=args.use_canonical_builder)


if __name__ == "__main__":
    main()
