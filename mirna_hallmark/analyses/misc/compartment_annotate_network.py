"""Roadmap Phase 2a: compartment-annotate the framework's regulatory network.

For every framework coupling edge, label *which cell compartment expresses the target gene* using the
Wu primary-breast scRNA atlas (GSE176078, `celltype_major` pseudobulk). This tells us, per edge,
whether a miRNA→target coupling is **tumour-epithelial-intrinsic** or **stromal/immune/endothelial
(paracrine / composition)** — the gene-resolved companion to the bulk composition retest (MH-57).
Prediction it tests: the edges that *attenuated* under composition adjustment (miR-29→collagen/ECM)
should have **CAF-expressed** targets, while the robust edges should be **epithelial-expressed**.

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.misc.compartment_annotate_network``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.brca_deconvolution import build_signature

OUT_DIR = C.OUTPUT_ROOT / "compartment_annotate_network"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"

# coarse compartments over Wu celltype_major
COARSE = {
    "Cancer Epithelial": "Epithelial", "Normal Epithelial": "Epithelial",
    "CAFs": "CAF",
    "T-cells": "Immune", "Myeloid": "Immune", "B-cells": "Immune", "Plasmablasts": "Immune",
    "Endothelial": "Endothelial", "PVL": "PVL",
}


def _tau(row: np.ndarray) -> float:
    """Tissue-specificity tau in [0,1]: 1 = compartment-specific, 0 = ubiquitous."""
    m = row.max()
    if m <= 0:
        return np.nan
    return float((1.0 - row / m).sum() / (len(row) - 1))


def gene_compartments() -> pd.DataFrame:
    sig = build_signature()                                       # genes x 9 celltype CPM
    # Assign each gene to its PEAK fine cell type, then map to compartment — unbiased by how many
    # fine subtypes a compartment has (sum inflates multi-subtype, mean dilutes it).
    frac = sig.div(sig.sum(axis=1).replace(0, np.nan), axis=0)    # fraction per fine cell type
    dom_fine = sig.idxmax(axis=1)
    out = pd.DataFrame(index=sig.index)
    out["dominant_celltype"] = dom_fine
    out["dominant_compartment"] = dom_fine.map(COARSE)
    out["dominant_fraction"] = frac.max(axis=1).round(3)          # fraction carried by the peak cell type
    out["tau_specificity"] = sig.apply(lambda r: _tau(r.to_numpy()), axis=1).round(3)
    return out


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[compartment] building Wu scRNA gene-compartment map ...")
    gc = gene_compartments()
    gc.to_csv(out_dir / "gene_compartment_map.tsv", sep="\t")

    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e["arm"] = e["key"].str.split(r"\|\|").str[0]
    e["gene"] = e["key"].str.split(r"\|\|").str[1]
    neg = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()          # headline BY-neg couplings
    ann = neg.join(gc, on="gene")
    ann.to_csv(out_dir / "edges_compartment_annotated.tsv", sep="\t", index=False)

    covered = ann["dominant_compartment"].notna()
    dist = ann.loc[covered, "dominant_compartment"].value_counts(normalize=True).round(3).to_dict()
    # miR-29 family edges (the MH-57 attenuated, ECM-targeting class) — where do their targets live?
    mir29 = ann[ann["arm"].str.contains("miR-29", na=False) & covered]
    summary = {
        "module": "mirna_hallmark.analyses.misc.compartment_annotate_network",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reference": "GSE176078 Wu 2021 scRNA (celltype_major -> 5 coarse compartments)",
        "n_headline_edges": int(len(neg)), "n_with_target_in_scRNA": int(covered.sum()),
        "headline_target_compartment_distribution": dist,
        "median_target_tau": round(float(ann.loc[covered, "tau_specificity"].median()), 3),
        "mir29_family": {
            "n_edges": int(len(mir29)),
            "target_compartment_distribution": mir29["dominant_compartment"].value_counts(normalize=True).round(3).to_dict(),
            "example_CAF_targets": mir29[mir29["dominant_compartment"] == "CAF"]["gene"].head(10).tolist(),
        },
        "epithelial_specific_examples": ann[(ann["dominant_compartment"] == "Epithelial") & (ann["dominant_fraction"] > 0.6)]
            .sort_values("dominant_fraction", ascending=False)[["arm", "gene", "dominant_fraction"]].head(8)
            .apply(lambda r: f"{r['arm']}->{r['gene']} ({r['dominant_fraction']})", axis=1).tolist(),
        "caveats": ["compartment = where the TARGET is expressed (scRNA), not where the miRNA acts",
                    "Wu tumour scRNA; coarse 5-compartment; ubiquitous genes (low tau) are ambiguous"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[compartment] {covered.sum()}/{len(neg)} headline targets mapped; distribution: {dist}")
    print(f"[compartment] miR-29 family targets: {summary['mir29_family']['target_compartment_distribution']} "
          f"(CAF examples: {summary['mir29_family']['example_CAF_targets'][:5]})")
    print(f"[compartment] wrote {out_dir}")
    return ann


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
