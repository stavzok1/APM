"""Roadmap Phase 3 (spatial) — resolution layer: SPOT (~55µm, multi-cell), PRIMARY-tumour RNA.

Wu et al. 2021 breast Visium (Zenodo 4739739; 6 primary tumours: 1142243F, 1160920F, CID4290, CID4465,
CID44971, CID4535) — the *primary-tumour* spot-resolution layer (the MIBI anchor is DCIS). Same atlas we
deconvolve in Phase 0, so the Wu celltype_major signature applies directly.

Per section: score the 5 programs per spot, deconvolve spots → 5 coarse compartments, localize each
program to its compartment (program×compartment-fraction Spearman), measure spatial niche structure
(Moran's I), and run the BRAKE-RELEASE co-localization arm of the spatial decoupling test — do
released-pressure target genes (SLC2A1, FAP) spatially track the program niches they serve?

What a single tumour section CANNOT give: the healthy→tumour axis (delta_tumor_gtex is a tumour-vs-
healthy contrast). So the full MH-56d concordant-repression arm stays with the MIBI anchor (which has
normal→DCIS→IBC); here the decoupling readout is the *within-section* program↔projected-pressure
co-variation, and the per-gene tumour-content correlation is reported but flagged epithelial-identity-
confounded (epithelial markers are trivially high in epithelial-rich spots).

Caveats: spots are multi-cell → epithelial deconvolution ρ=0.39 unreliable (MH-58), lean on immune
(0.85)/CAF (0.57); route (i) projected pressure is composition-driven by construction.

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.spatial.spatial_visium_wu``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.spatial import spatial_common as SP

WU_DIR = C.REPO_ROOT / "data" / "external" / "brca_visium_wu"
MTX_DIR = WU_DIR / "filtered_count_matrices"
SPATIAL_DIR = WU_DIR / "spatial"
OUT_DIR = C.OUTPUT_ROOT / "spatial_visium_wu"


def _samples() -> List[str]:
    return sorted(p.name.replace("_filtered_count_matrix", "") for p in MTX_DIR.glob("*_filtered_count_matrix"))


def analyse_section(sample: str, sig: pd.DataFrame, program_genes: Dict[str, List[str]]) -> dict:
    counts = SP.read_10x_mtx(MTX_DIR / f"{sample}_filtered_count_matrix")
    pos = SP.read_visium_positions(SPATIAL_DIR / f"{sample}_spatial" / "tissue_positions_list.csv")
    spots = counts.columns.intersection(pos.index[pos["in_tissue"] == 1])
    counts = counts[spots]
    log = SP.cpm_log(counts)
    lin = np.power(2.0, log) - 1.0                              # back to ~CPM-linear for NNLS
    coords = pos.reindex(spots)[["x", "y"]]
    return {"sample": sample, **SP.analyse_spots(log, lin, coords, sig, program_genes)}


def run(*, out_dir: Path = OUT_DIR, samples: List[str] | None = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = samples or _samples()
    print(f"[visium_wu] building Wu signature once; {len(samples)} sections: {samples}")
    sig = SP.build_signature()
    pg = SP.program_gene_sets()

    per = []
    for s in samples:
        try:
            r = analyse_section(s, sig, pg)
            per.append(r)
            print(f"   {s}: {r['n_spots']} spots | Moran's I {r['morans_i']}")
        except Exception as e:  # a bad section shouldn't sink the suite
            print(f"   {s}: SKIPPED ({type(e).__name__}: {e})")

    # pool: median localization + Moran's I across sections
    progs = list(pg.keys())
    comps = sorted({c for r in per for c in r["localization"][progs[0]]}) if per else []
    loc_med = pd.DataFrame(
        {p: {c: round(float(np.nanmedian([r["localization"][p][c] for r in per if c in r["localization"][p]])), 3)
             for c in comps} for p in progs}).T
    loc_med.to_csv(out_dir / "program_compartment_localization.tsv", sep="\t")
    moran_med = {p: round(float(np.nanmedian([r["morans_i"][p] for r in per])), 3) for p in progs} if per else {}

    # pool anchor brake-release co-localization
    arows = []
    for g in SP.ANCHOR_GENES:
        vals = [r["anchor"][g] for r in per if g in r["anchor"]]
        if not vals:
            continue
        arows.append({"gene": g, "pressure": vals[0]["pressure"], "n_sections": len(vals),
                      "median_rho_best_program": round(float(np.nanmedian([v["rho_best_program"] for v in vals])), 3),
                      "modal_best_program": pd.Series([v["best_program"] for v in vals]).mode().iloc[0],
                      "median_rho_epithelial_frac": round(float(np.nanmedian([v["rho_epithelial_frac"] for v in vals])), 3),
                      "median_morans_i": round(float(np.nanmedian([v["morans_i"] for v in vals])), 3)})
    anchor_df = pd.DataFrame(arows)
    anchor_df.to_csv(out_dir / "anchor_brake_release_localization.tsv", sep="\t", index=False)
    pd.DataFrame([{**{"sample": r["sample"], "n_spots": r["n_spots"]}, **r["morans_i"]} for r in per]).to_csv(
        out_dir / "per_section_morans_i.tsv", sep="\t", index=False)

    summary = {
        "module": "mirna_hallmark.analyses.spatial.spatial_visium_wu",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "resolution": "Visium spot ~55µm (multi-cell), PRIMARY breast tumours",
        "route": "measured mRNA (programs+targets); miRNA side = route(i) composition projection only",
        "reference": "Wu et al 2021 Visium (Zenodo 4739739)",
        "n_sections": len(per), "sections": [r["sample"] for r in per],
        "median_morans_i_per_program": moran_med,
        "program_compartment_localization_median": loc_med.round(3).to_dict("index"),  # {program: {compartment: rho}}
        "anchor_brake_release": anchor_df.set_index("gene").to_dict("index") if not anchor_df.empty else {},
        "caveats": ["spots multi-cell; epithelial deconv ρ=0.39 unreliable (MH-58) — lean on immune/CAF",
                    "single tumour sections lack a healthy axis → concordant-repression arm stays w/ MIBI anchor",
                    "per-gene epithelial-frac corr is epithelial-IDENTITY-confounded (reported, not headlined)",
                    "route(i) projected pressure is composition-driven by construction"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[visium_wu] median Moran's I: {moran_med}")
    print(f"[visium_wu] program→compartment (median):\n{loc_med.round(2)}")
    if not anchor_df.empty:
        print(f"[visium_wu] anchor brake-release co-localization:\n{anchor_df.to_string(index=False)}")
    print(f"[visium_wu] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, samples=args.samples)


if __name__ == "__main__":
    main()
