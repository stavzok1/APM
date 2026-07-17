"""Roadmap Phase 3 (spatial) — resolution layer: SINGLE-CELL in-situ RNA (10x Xenium breast).

The finest spatial-RNA resolution: 10x Xenium Human Breast Cancer (Rep1; ~168k cells, 313-gene panel,
x/y centroids). Unlike Visium spots, cells are typed directly (no deconvolution) — so program
LOCALIZATION here is true single-cell, the Xenium advantage. The cost is the small targeted panel
(~313 genes), so only part of each Hallmark program is measurable and most released-pressure anchor
genes (SLC2A1, FAP, VIM, HIF1A) are OFF the panel — only gained-pressure anchors ERBB2/ESR1/CDH1 are on
it, so the brake-release arm is weak here (reported, not headlined).

Readouts:
  1. cell typing → 5 coarse compartments via panel marker metagenes (argmax).
  2. program LOCALIZATION: per-cell program z (panel subset) averaged by compartment.
  3. niche structure: Moran's I per program on centroids (single-cell spatial autocorrelation) +
     kNN neighbourhood compartment composition.
  4. anchor co-localization (panel anchors only): gene vs program / epithelial-cell-ness.

Run: ``.venv/bin/python3 -m mirna_hallmark.spatial_xenium``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.spatial import spatial_common as SP
from mirna_hallmark.analyses.cross_state.cross_state_coupling import EPI_MARKERS, IMMUNE_MARKERS, STROMA_MARKERS
from mirna_hallmark.stats import zscore_rows

XDIR = C.REPO_ROOT / "data" / "external" / "brca_xenium"
MTX = XDIR / "cell_feature_matrix"
CELLS = XDIR / "Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz"
OUT_DIR = C.OUTPUT_ROOT / "spatial_xenium"

# coarse single-cell compartments from panel marker metagenes (endothelial markers if on panel)
COMPARTMENT_MARKERS = {"Epithelial": EPI_MARKERS, "Immune": IMMUNE_MARKERS, "CAF": STROMA_MARKERS,
                       "Endothelial": ["PECAM1", "VWF", "CLDN5", "CD34"]}


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[xenium] loading single-cell matrix + centroids ...")
    counts = SP.read_10x_mtx(MTX)
    cells = pd.read_csv(CELLS)
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    counts.columns = counts.columns.astype(str)
    shared = counts.columns.intersection(cells.index)
    counts = counts[shared]
    coords = cells.reindex(shared)[["x_centroid", "y_centroid"]].rename(
        columns={"x_centroid": "x", "y_centroid": "y"})
    log = SP.cpm_log(counts)
    z = zscore_rows(log)
    print(f"[xenium] {log.shape[1]:,} cells × {log.shape[0]} panel genes")

    # 1. cell typing → compartment (argmax marker metagene)
    metas = {}
    for comp, mk in COMPARTMENT_MARKERS.items():
        present = [g for g in mk if g in z.index]
        metas[comp] = z.loc[present].mean(axis=0) if present else pd.Series(-np.inf, index=z.columns)
    meta_df = pd.DataFrame(metas).fillna(-np.inf)
    compartment = meta_df.idxmax(axis=1)
    compartment[meta_df.max(axis=1) < 0.1] = "Other"            # weakly-typed cells
    comp_counts = compartment.value_counts().to_dict()

    # 2. program localization (per-cell z averaged by compartment) — true single-cell
    pg = SP.program_gene_sets()
    prog = SP.score_programs(log, pg)                           # cells × programs (panel subset)
    prog_by_comp = prog.groupby(compartment).mean().round(3)
    prog_by_comp.to_csv(out_dir / "program_by_compartment.tsv", sep="\t")
    panel_cov = {p: int(sum(g in log.index for g in genes)) for p, genes in pg.items()}

    # 3. niche structure: Moran's I per program (single-cell centroids) + kNN compartment mixing
    moran = {p: round(SP.morans_i(prog[p], coords), 3) for p in prog.columns}
    # kNN neighbourhood compartment enrichment (are compartments spatially clustered?)
    xy = coords.to_numpy(float)
    tree = cKDTree(xy)
    _, nbr = tree.query(xy, k=16)
    comp_code = compartment.to_numpy()
    same = np.mean([(comp_code[nbr[i, 1:]] == comp_code[i]).mean() for i in range(0, len(xy), 25)])  # subsample
    niche_homophily = round(float(same), 3)

    # 4. anchor co-localization (panel anchors only)
    delta = SP.load_pressure_delta()
    epi_ind = (compartment == "Epithelial").astype(float)
    anchor = {}
    for g in SP.ANCHOR_GENES:
        if g not in log.index:
            continue
        gv = log.loc[g]
        prog_rho = {p: float(spearmanr(gv, prog[p])[0]) for p in prog.columns}
        best = max(prog_rho, key=lambda k: abs(prog_rho[k]))
        anchor[g] = {"pressure": "released" if delta.get(g, 0) < 0 else "gained",
                     "best_program": best, "rho_best_program": round(prog_rho[best], 3),
                     "rho_epithelial_cell": round(float(spearmanr(gv, epi_ind)[0]), 3),
                     "morans_i": round(SP.morans_i(gv, coords), 3)}

    summary = {
        "module": "mirna_hallmark.spatial_xenium",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "resolution": "single-cell in-situ RNA (Xenium), PRIMARY breast tumour",
        "route": "measured mRNA (panel-limited); miRNA side not measured (localization layer)",
        "reference": "10x Genomics Xenium Human Breast Cancer Rep1 (313-gene panel)",
        "n_cells": int(log.shape[1]), "compartment_counts": comp_counts,
        "program_panel_coverage": panel_cov,
        "program_by_compartment": prog_by_comp.to_dict("index"),
        "morans_i_per_program": moran, "niche_homophily_knn15": niche_homophily,
        "anchor_colocalization_panel_only": anchor,
        "caveats": ["313-gene panel → programs partly measured (proliferation/hypoxia/glycolysis ~6-8 genes); "
                    "EMT/immune well-covered (~20)",
                    "released-pressure anchors (SLC2A1/FAP/VIM/HIF1A) OFF panel → brake-release arm weak here; "
                    "only gained anchors ERBB2/ESR1/CDH1 present",
                    "single-cell counts sparse (~28 transcripts/cell) → per-cell program z noisy, denoised by "
                    "compartment averaging; no healthy axis (localization layer)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[xenium] compartments: {comp_counts}")
    print(f"[xenium] program panel coverage: {panel_cov}")
    print(f"[xenium] program by compartment:\n{prog_by_comp}")
    print(f"[xenium] Moran's I: {moran} | niche homophily(k15): {niche_homophily}")
    print(f"[xenium] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
