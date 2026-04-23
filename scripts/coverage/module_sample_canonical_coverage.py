#!/usr/bin/env python3
"""
One-shot coverage check: apply canonical gene/miRNA alias registries to a **sample**
of each major input type and count how many symbols are remapped.

Output: ``analysis/output/module_sample_canonical_coverage.md`` (default).

Run from repo root:
  python scripts/coverage/module_sample_canonical_coverage.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

import pandas as pd

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


def _lines(title: str, rows: List[Tuple[str, str]]) -> List[str]:
    out = [f"\n## {title}\n", "| Metric | Value |\n|---|---|\n"]
    for k, v in rows:
        out.append(f"| {k} | {v} |\n")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        type=Path,
        default=_REPO / "analysis" / "output" / "module_sample_canonical_coverage.md",
    )
    args = ap.parse_args()

    from pipeline.config import (
        PATHS,
        PRIMARY_GENES,
        PIPELINE_GENE_PANEL,
        PANEL_ALIAS_SEED_SYMBOLS,
        PANEL_MIRNA_TIER_ARM_IDS,
    )
    from pipeline.genes.panel_alias_registry import (
        get_gene_flat_alias_map,
        get_gene_canonical_to_aliases,
        get_mirna_flat_alias_map,
        normalize_mirna_symbol,
    )
    from pipeline.genes.symbol_normalization import apply_symbol_mapping_series

    flat = get_gene_flat_alias_map(None)
    inv = get_gene_canonical_to_aliases(None)
    mflat = get_mirna_flat_alias_map()

    lines: List[str] = []
    lines.append("# Module sample — canonical symbol coverage\n\n")
    lines.append(f"- Panel genes (seed count): **{len(PANEL_ALIAS_SEED_SYMBOLS)}**\n")
    n_ens_flat = sum(1 for k in flat if str(k).startswith("ENSG"))
    lines.append(f"- Flat gene alias entries: **{len(flat)}** (includes **{n_ens_flat}** Ensembl gene-id keys)\n")
    lines.append(f"- Canonicals with ≥1 synonym set: **{len(inv)}**\n")
    lines.append(f"- Manual + tier miRNA alias keys: **{len(mflat)}**\n")
    lines.append(
        "- miRNA tier arms (``config.PANEL_MIRNA_TIER_ARM_IDS``): **"
        f"{len(PANEL_MIRNA_TIER_ARM_IDS)}** mature-arm ids\n"
    )

    # --- RNA (TPM matrix) ---
    rna_path = Path(PATHS.rna_expression) if Path(PATHS.rna_expression).is_file() else Path(
        PATHS.rna_expression_raw
    )
    if rna_path.is_file():
        col = PATHS.rna_gene_col
        hdr = pd.read_csv(
            rna_path, sep="\t", nrows=0, encoding="latin-1", encoding_errors="replace"
        )
        use = col if col in hdr.columns else list(hdr.columns)[0]
        df = pd.read_csv(
            rna_path,
            sep="\t",
            usecols=[use],
            nrows=25_000,
            encoding="latin-1",
            encoding_errors="replace",
            on_bad_lines="skip",
        )
        raw = df[use].astype(str)
        new = apply_symbol_mapping_series(raw, flat)
        n_change = int((raw != new).sum())
        tier_primary_seen = sum(1 for g in PRIMARY_GENES if g in set(raw))
        lines += _lines(
            "RNA expression (first 25k rows)",
            [
                ("path", str(rna_path)),
                ("rows scanned", str(len(df))),
                ("symbols remapped to canonical", str(n_change)),
                ("primary-panel genes appearing as row labels", str(tier_primary_seen)),
            ],
        )
    else:
        lines.append("\n## RNA expression\n\n*(file missing)*\n")

    # --- miRNA matrix (XENA) ---
    mir_path = Path(PATHS.mirna_expression_tsv)
    if mir_path.is_file():
        mdf = pd.read_csv(
            mir_path, sep="\t", nrows=8000, encoding="latin-1", encoding_errors="replace"
        )
        idx = mdf.iloc[:, 0].astype(str)
        new_idx = idx.map(lambda s: normalize_mirna_symbol(s))
        n_m = int((idx != new_idx).sum())
        norm_set = set(new_idx.astype(str))
        tier_seen = sum(1 for a in PANEL_MIRNA_TIER_ARM_IDS if a in norm_set)
        lines += _lines(
            "miRNA expression matrix (first 8000 row labels)",
            [
                ("path", str(mir_path)),
                ("labels remapped (MIMAT→arm + spelling aliases)", str(n_m)),
                ("tier arms (``PANEL_MIRNA_TIER_ARM_IDS``) present after normalize", str(tier_seen)),
            ],
        )
    else:
        lines.append("\n## miRNA expression\n\n*(file missing)*\n")

    # --- Methylation probe reference ---
    pref = Path(PATHS.methylation_probe_reference)
    if pref.is_file():
        from pipeline.Methylation.probe_loader import load_probe_reference

        sub = load_probe_reference(pref)
        raw_lists = sub["gene_list"].head(15_000)
        n_tokens = 0
        n_changed = 0
        for lst in raw_lists:
            if not isinstance(lst, list):
                continue
            for t in lst:
                n_tokens += 1
                c = flat.get(str(t).strip(), str(t).strip())
                if c != str(t).strip():
                    n_changed += 1
        lines += _lines(
            "Methylation probe reference (first 15k rows, token-wise)",
            [
                ("path", str(pref)),
                ("gene tokens scanned", str(n_tokens)),
                ("tokens remapped", str(n_changed)),
            ],
        )
    else:
        lines.append("\n## Methylation probes\n\n*(file missing)*\n")

    # --- RPPA annotation ---
    rppa_ann = Path(PATHS.rppa_antibody_annotation_csv)
    if rppa_ann.is_file():
        df = pd.read_csv(rppa_ann, sep="\t", encoding="latin-1", encoding_errors="replace")
        if "gene_name" in df.columns:
            g = df["gene_name"].astype(str)
            g2 = apply_symbol_mapping_series(g, flat)
            lines += _lines(
                "RPPA antibody table",
                [
                    ("path", str(rppa_ann)),
                    ("rows", str(len(df))),
                    ("gene_name cells remapped", str(int((g != g2).sum()))),
                ],
            )
    else:
        lines.append("\n## RPPA annotation\n\n*(file missing)*\n")

    # --- TargetScan ---
    tsp = Path(PATHS.targetscan_predictions)
    if tsp.is_file():
        df = pd.read_csv(tsp, sep="\t", nrows=60_000, encoding="latin-1", encoding_errors="replace")
        if "Gene Symbol" in df.columns:
            g = df["Gene Symbol"].astype(str)
            g2 = apply_symbol_mapping_series(g, flat)
            lines += _lines(
                "TargetScan (first 60k rows)",
                [
                    ("path", str(tsp)),
                    ("Gene Symbol remapped", str(int((g != g2).sum()))),
                ],
            )
    else:
        lines.append("\n## TargetScan\n\n*(file missing)*\n")

    # --- miRTarBase ---
    mt = Path(PATHS.mirtarbase_csv)
    if mt.is_file():
        df = pd.read_csv(mt, nrows=120_000, encoding="latin-1", encoding_errors="replace")
        gcol = "Target Gene" if "Target Gene" in df.columns else None
        if gcol:
            g = df[gcol].astype(str).str.strip()
            g2 = apply_symbol_mapping_series(g, flat)
            lines += _lines(
                "miRTarBase (first 120k rows)",
                [
                    ("path", str(mt)),
                    ("target gene cells remapped", str(int((g != g2).sum()))),
                ],
            )
    else:
        lines.append("\n## miRTarBase\n\n*(file missing)*\n")

    # --- genes_only (annotation used by ATAC / methylation / SV) ---
    gonly = Path(PATHS.genes_only)
    if gonly.is_file():
        gdf = pd.read_csv(gonly, nrows=20_000, low_memory=False)
        if "gene_name" in gdf.columns:
            g = gdf["gene_name"].astype(str)
            g2 = apply_symbol_mapping_series(g, flat)
            lines += _lines(
                "GENCODE genes_only CSV (first 20k rows)",
                [
                    ("path", str(gonly)),
                    ("gene_name cells remapped", str(int((g != g2).sum()))),
                ],
            )
    else:
        lines.append("\n## genes_only annotation\n\n*(file missing)*\n")

    # --- Panel tier presence in flat map keys (spot) ---
    lines.append("\n## Pipeline gene panel vs remap keys\n")
    lines.append("| set | n genes | panel symbols that are keys of an alias remap |\n|---|---:|---:|\n")
    for label, genes in (
        ("PRIMARY", PRIMARY_GENES),
        ("PIPELINE_GENE_PANEL", PIPELINE_GENE_PANEL),
    ):
        gset = set(genes)
        as_key = sum(1 for g in gset if g in flat)
        lines.append(f"| {label} | {len(gset)} | {as_key} |\n")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
