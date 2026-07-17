"""Roadmap Phase 3 (spatial) — resolution layer: SPOT, INDEPENDENT replication (public 10x Visium).

10x Genomics public Human Breast Cancer (Block A, Section 1) Visium block — an independent spot-resolution
section to replicate the Wu-Visium program localization + brake-release co-localization (`spatial_visium_wu`).
Different patient, different lab/pipeline → cross-dataset robustness for the spot layer.

Same shared analysis as Wu (`spatial_common.analyse_spots`): program scores per spot, deconvolution →
compartment localization, Moran's I niche structure, anchor brake-release co-localization. Same caveats
(multi-cell spots, no in-section healthy axis, route-(i) composition projection only).

Run: ``.venv/bin/python3 -m mirna_hallmark.spatial_visium_public``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.spatial import spatial_common as SP

X_DIR = C.REPO_ROOT / "data" / "external" / "brca_visium_10x"
H5 = X_DIR / "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
POS = X_DIR / "spatial" / "tissue_positions_list.csv"
OUT_DIR = C.OUTPUT_ROOT / "spatial_visium_public"


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[visium_10x] building Wu signature; loading 10x breast block ...")
    sig = SP.build_signature()
    pg = SP.program_gene_sets()

    counts = SP.read_10x_h5(H5)
    pos = SP.read_visium_positions(POS)
    spots = counts.columns.intersection(pos.index)
    counts = counts[spots]
    log = SP.cpm_log(counts)
    lin = np.power(2.0, log) - 1.0
    coords = pos.reindex(spots)[["x", "y"]]
    r = SP.analyse_spots(log, lin, coords, sig, pg)

    loc = pd.DataFrame(r["localization"]).T
    loc.to_csv(out_dir / "program_compartment_localization.tsv", sep="\t")
    anchor_df = pd.DataFrame([{"gene": g, **v} for g, v in r["anchor"].items()])
    anchor_df.to_csv(out_dir / "anchor_brake_release_localization.tsv", sep="\t", index=False)

    summary = {
        "module": "mirna_hallmark.spatial_visium_public",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "resolution": "Visium spot ~55µm (multi-cell), INDEPENDENT 10x breast block",
        "route": "measured mRNA; miRNA side = route(i) composition projection only",
        "reference": "10x Genomics Human Breast Cancer Block A Section 1 (Space Ranger 1.1.0)",
        "n_spots": r["n_spots"], "morans_i_per_program": r["morans_i"],
        "program_compartment_localization": r["localization"],
        "anchor_brake_release": r["anchor"],
        "caveats": ["single public block (n=1 section, no outcome) — replication of the spot pattern only",
                    "spots multi-cell; epithelial deconv unreliable (MH-58); no healthy axis in-section",
                    "route(i) projected pressure composition-driven"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[visium_10x] {r['n_spots']} spots | Moran's I {r['morans_i']}")
    print(f"[visium_10x] program→compartment:\n{loc.round(2)}")
    if not anchor_df.empty:
        print(f"[visium_10x] anchor brake-release:\n{anchor_df.to_string(index=False)}")
    print(f"[visium_10x] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
