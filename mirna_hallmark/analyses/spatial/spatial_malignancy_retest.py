"""Roadmap Phase 3 (spatial) — push #1: the FULL MH-56(d) retest on PRIMARY tumour Visium.

`spatial_visium_wu`/`_public` could only test the brake-release arm on primary tissue, because a single
tumour section has no healthy reference — so the concordant-repression arm rested on MIBI (DCIS). This
module recovers the missing axis: an **inferCNV-style per-spot MALIGNANCY score** (genome-wide CNV energy
relative to the section's immune/endothelial = diploid reference) is a within-section *malignant→normal
epithelium* gradient on the **mRNA side** — the spatial analogue of `delta_tumor_gtex`. Restricting to
epithelial-rich spots removes the stroma contrast (so the axis is malignant-epi vs normal-epi, NOT
epithelial-identity), then the standard `gene_level_pressure_retest` runs **both arms** on primary
invasive tissue.

Circularity guard: a gene on a recurrently-altered arm could rise with malignancy via its OWN CNV. Each
anchor is therefore also scored against a malignancy computed with **its own chromosome excluded**
(`exclude_chr`); the clean call is the own-chromosome-removed one.

Still no spatial miRNA (binding constraint) — the pressure side is the projected bulk delta; only the
readout axis is now primary-invasive instead of DCIS.

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.spatial.spatial_malignancy_retest``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.spatial import spatial_common as SP
from mirna_hallmark.analyses.cross_state.cross_state_coupling import EPI_MARKERS
from mirna_hallmark.analyses.spatial.spatial_visium_wu import MTX_DIR, SPATIAL_DIR, _samples

OUT_DIR = C.OUTPUT_ROOT / "spatial_malignancy_retest"
TEN_X = C.REPO_ROOT / "data" / "external" / "brca_visium_10x"


def _load_visium(sample: str):
    if sample == "10x_public":
        counts = SP.read_10x_h5(TEN_X / "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
        pos = SP.read_visium_positions(TEN_X / "spatial" / "tissue_positions_list.csv")
    else:
        counts = SP.read_10x_mtx(MTX_DIR / f"{sample}_filtered_count_matrix")
        pos = SP.read_visium_positions(SPATIAL_DIR / f"{sample}_spatial" / "tissue_positions_list.csv")
    spots = counts.columns.intersection(pos.index[pos.get("in_tissue", 1) == 1]) if "in_tissue" in pos \
        else counts.columns.intersection(pos.index)
    counts = counts[spots]
    log = SP.cpm_log(counts)
    lin = np.power(2.0, log) - 1.0
    return log, lin, pos.reindex(spots)


def analyse_section(sample: str, gene_pos, anchors: List[str]) -> dict:
    log, lin, pos = _load_visium(sample)
    # Marker metagenes (RELIABLE, unlike NNLS epithelial fraction which collapses to ~1.0 everywhere, MH-58):
    # epithelial-marker score gates the epithelial compartment and picks the diploid reference.
    z = SP.zscore_rows(log)
    epi_mk = z.loc[[g for g in EPI_MARKERS if g in z.index]].mean()
    ref_mask = epi_mk < epi_mk.quantile(0.33)                    # bottom-tertile epithelial = stroma/immune = diploid
    malig = SP.infercnv_malignancy(log, ref_mask, gene_pos, window=100)
    epi_spots = epi_mk.index[epi_mk > epi_mk.median()]           # within-epithelium malignant→normal gradient
    # sanity: within epithelium, malignancy should track proliferation (more malignant epi proliferates more)
    prolif = SP.score_programs(log, SP.program_gene_sets())["proliferation"]
    rho_prolif = float(spearmanr(malig.reindex(epi_spots), prolif.reindex(epi_spots), nan_policy="omit")[0])

    retest = SP.gene_level_pressure_retest(log[epi_spots], malig.reindex(epi_spots), min_units=40)
    # own-chromosome-removed malignancy for the anchors (circularity guard)
    anchor_clean = {}
    for g in anchors:
        if g not in log.index or g not in gene_pos.index:
            continue
        chrn = int(gene_pos.loc[g, "chrn"])
        m_excl = SP.infercnv_malignancy(log, ref_mask, gene_pos, exclude_chr=chrn).reindex(epi_spots)
        v = log.loc[g, epi_spots]
        anchor_clean[g] = round(float(spearmanr(v, m_excl, nan_policy="omit")[0]), 3)
    return {"sample": sample, "n_spots": int(log.shape[1]), "n_epi_spots": int(len(epi_spots)),
            "malignancy_vs_proliferation": round(rho_prolif, 3), "retest": retest,
            "anchor_rho_excl_own_chr": anchor_clean}


def run(*, out_dir: Path = OUT_DIR, samples: List[str] | None = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = samples or (_samples() + ["10x_public"])
    print(f"[malig] gene positions; {len(samples)} sections (marker-based reference, no deconvolution)")
    gene_pos = SP.load_gene_positions()
    anchors = SP.ANCHOR_GENES

    per = []
    for s in samples:
        try:
            r = analyse_section(s, gene_pos, anchors)
            per.append(r)
            n_cls = r["retest"]["classification"].value_counts().to_dict() if not r["retest"].empty else {}
            print(f"   {s}: {r['n_epi_spots']}/{r['n_spots']} epi spots | malig~prolif ρ={r['malignancy_vs_proliferation']} | {n_cls}")
        except Exception as e:
            print(f"   {s}: SKIPPED ({type(e).__name__}: {e})")

    # pool the retest across sections by gene (median rho + modal classification)
    allret = pd.concat([r["retest"].assign(sample=r["sample"]) for r in per if not r["retest"].empty], ignore_index=True)
    pooled = (allret.groupby("gene")
              .agg(pressure=("pressure", "first"), role=("role", "first"), is_anchor=("is_anchor", "first"),
                   n_sections=("sample", "nunique"),
                   median_rho=("spatial_tumour_rho", "median"),
                   modal_class=("classification", lambda s: s.mode().iloc[0]))
              .reset_index())
    # require presence in >=4 sections for a stable call
    pooled = pooled[pooled["n_sections"] >= 4].sort_values(["is_anchor", "median_rho"], ascending=[False, True])
    pooled.to_csv(out_dir / "malignancy_retest_pooled.tsv", sep="\t", index=False)
    allret.to_csv(out_dir / "malignancy_retest_per_section.tsv.gz", sep="\t", index=False)

    # anchor table with own-chromosome-removed clean rho
    arows = []
    for g in anchors:
        clean = [r["anchor_rho_excl_own_chr"].get(g) for r in per if g in r["anchor_rho_excl_own_chr"]]
        clean = [c for c in clean if c is not None and np.isfinite(c)]
        prow = pooled[pooled["gene"] == g]
        if prow.empty or not clean:
            continue
        arows.append({"gene": g, "pressure": prow["pressure"].iloc[0],
                      "median_rho_vs_malignancy": round(float(prow["median_rho"].iloc[0]), 3),
                      "modal_class": prow["modal_class"].iloc[0],
                      "median_rho_excl_own_chr": round(float(np.median(clean)), 3),
                      "n_sections": int(prow["n_sections"].iloc[0])})
    anchor_df = pd.DataFrame(arows)
    anchor_df.to_csv(out_dir / "anchor_malignancy_retest.tsv", sep="\t", index=False)

    cls_counts = pooled["modal_class"].value_counts().to_dict()
    summary = {
        "module": "mirna_hallmark.analyses.spatial.spatial_malignancy_retest",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "resolution": "Visium spot, PRIMARY tumours — within-section malignancy axis (inferCNV-style)",
        "route": "measured mRNA; malignancy = genome-wide CNV energy vs diploid reference; pressure = projected bulk delta",
        "n_sections": len(per), "sections": [r["sample"] for r in per],
        "pooled_classification_counts": cls_counts,
        "anchor_calls": anchor_df.set_index("gene").to_dict("index") if not anchor_df.empty else {},
        "median_malignancy_vs_proliferation_within_epi": round(float(np.median([r["malignancy_vs_proliferation"] for r in per])), 3),
        "caveats": [
            "malignancy axis is mRNA-derived (no spatial miRNA) — pressure side is still the projected bulk delta",
            "restricted to epithelial-rich spots → malignant-epi vs normal-epi (removes the epithelial-identity confound, "
            "unlike the raw epithelial-frac correlation)",
            "own-CNV circularity controlled by `*_excl_own_chr` (malignancy recomputed without the gene's chromosome)",
            "inferCNV reference = immune/endothelial spots; CNV-from-expression is approximate; CID sections TNBC-enriched",
        ],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[malig] pooled classification (epi spots): {cls_counts}")
    if not anchor_df.empty:
        print(f"[malig] anchor calls (incl. own-chr-removed clean rho):\n{anchor_df.to_string(index=False)}")
    print(f"[malig] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, samples=args.samples)


if __name__ == "__main__":
    main()
