"""Roadmap Phase 3 (spatial) — resolution layer: SINGLE-CELL PROTEIN anchor (Risom DCIS MIBI).

The measured-protein, in-situ anchor for the spatial program-localization + MH-56(d) decoupling retest.
Built FIRST because it needs no download (Risom is already staged) and validates the
`spatial_common.gene_level_pressure_retest` engine end-to-end: it must reproduce the MH-56d in-situ
calls — **SLC2A1/GLUT1 brake-release** (pressure released, protein up in tumour) and **ERBB2/HER2,
ESR1/ER, CDH1/ECAD concordant-repression** (pressure gained, protein down in tumour) — before the
RNA-Visium/Xenium layers are trusted.

Two readouts:
  1. PROGRAM LOCALIZATION (markers, localization-only per the MH-56 caveat): the 5 programs' single
     protein markers across normal→DCIS→IBC in the TUMOUR compartment (reproduces MH-56d Ki67/HIF1a/
     GLUT1/IDO1 rise at invasion) and across compartments.
  2. SPATIAL MH-56d RETEST (rigorous, gene-by-gene): for the framework-target proteins, correlate the
     per-cell protein with tumour-lineage membership → brake-release / concordant-repression / discordant.

Caveats: MIBI is **DCIS, not primary-invasive**, and **protein not miRNA** — this is the in-situ bridge/
anchor, not the primary-tumour claim. Single markers for programs (localization only).

Run: ``.venv/bin/python3 -m mirna_hallmark.spatial_mibi_anchor``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.spatial import spatial_common as SP

CELLTABLE = C.REPO_ROOT / "data" / "external" / "risom2022_dcis_mibi" / "risom2022_single_cell.csv"
OUT_DIR = C.OUTPUT_ROOT / "spatial_mibi_anchor"

# MIBI protein channel → framework target gene (the genes whose MH-56d pressure call we re-test).
PROTEIN_TO_GENE = {"GLUT1": "SLC2A1", "FAP": "FAP", "HER2": "ERBB2",
                   "ER": "ESR1", "ECAD": "CDH1", "VIM": "VIM", "HIF1a": "HIF1A"}
# Single-marker proxy per program (LOCALIZATION ONLY — not the pressure correspondence; MH-56 lesson).
PROGRAM_MARKER = {"proliferation": "Ki67", "hypoxia": "HIF1a", "glycolysis": "GLUT1",
                  "EMT": "VIM", "immune": "CD45"}
STAGE_ORDER = ["normal", "DCIS", "IBC"]


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CELLTABLE)
    lin = "celllineage"
    df["is_tumour"] = (df[lin].astype(str).str.lower() == "tumor").astype(float)
    print(f"[mibi] {len(df):,} cells | tumour {int(df['is_tumour'].sum()):,} | "
          f"stages {df['Tissue_Type'].value_counts().to_dict()}")

    # ---- 1. PROGRAM LOCALIZATION: program markers in tumour cells across normal→DCIS→IBC ----
    tum = df[df["is_tumour"] == 1]
    loc_rows = []
    for prog, mk in PROGRAM_MARKER.items():
        if mk not in df.columns:
            continue
        by_stage = {s: round(float(tum.loc[tum["Tissue_Type"] == s, mk].mean()), 4) for s in STAGE_ORDER}
        # monotone normal→IBC trend (Spearman of marker vs ordinal stage, tumour cells)
        ts = tum["Tissue_Type"].map({s: i for i, s in enumerate(STAGE_ORDER)})
        ok = ts.notna() & tum[mk].notna()
        rho = float(spearmanr(ts[ok], tum.loc[ok, mk])[0]) if ok.sum() > 20 else np.nan
        loc_rows.append({"program": prog, "marker": mk, **{f"tumour_{s}": by_stage[s] for s in STAGE_ORDER},
                         "stage_trend_rho": round(rho, 3)})
    loc = pd.DataFrame(loc_rows)
    loc.to_csv(out_dir / "program_localization_by_stage.tsv", sep="\t", index=False)

    # program markers by compartment (where each program lives)
    comp_col = "compartment" if "compartment" in df.columns else lin
    comp_means = {}
    for prog, mk in PROGRAM_MARKER.items():
        if mk in df.columns:
            comp_means[prog] = df.groupby(comp_col)[mk].mean().round(4)
    comp = pd.DataFrame(comp_means)
    comp.to_csv(out_dir / "program_marker_by_compartment.tsv", sep="\t")

    # ---- 2. SPATIAL MH-56d RETEST on framework-target proteins (gene-by-gene) ----
    # The faithful MH-56d axis is the healthy→tumour PROGRESSION (delta_tumor_gtex is a tumour-vs-healthy
    # contrast), NOT tumour-vs-stroma membership (that just re-reads epithelial identity — ER/HER2/ECAD are
    # epithelial proteins). So restrict to the epithelial/tumour lineage (which spans normal→DCIS→IBC) and
    # score each cell by ordinal stage. SLC2A1/GLUT1 (released) rise at invasion → brake-release; ERBB2/ESR1/
    # CDH1 (gained) fall at invasion (TNBC subtype shift) → concordant-repression.
    epi = df[df[lin].astype(str).str.lower() == "tumor"].copy()
    stage_ord = epi["Tissue_Type"].map({s: i for i, s in enumerate(STAGE_ORDER)})
    epi = epi[stage_ord.notna()]
    progression = stage_ord.dropna()                         # 0 normal → 1 DCIS → 2 IBC (epithelial cells)
    chans = [p for p in PROTEIN_TO_GENE if p in df.columns]
    readout = epi[chans].T.copy()
    readout.index = [PROTEIN_TO_GENE[p] for p in chans]      # genes × epithelial-cells (proteins)
    readout.columns = epi.index
    retest = SP.gene_level_pressure_retest(readout, progression,
                                           target_genes=list(PROTEIN_TO_GENE.values()), min_units=200)
    retest.to_csv(out_dir / "mh56d_spatial_retest.tsv", sep="\t", index=False)

    cls = retest.set_index("gene")["classification"].to_dict() if not retest.empty else {}
    # validation gate: the engine must reproduce the MH-56d anchor calls
    validation = {
        "SLC2A1_brake_release": cls.get("SLC2A1") == "brake_release",
        "ERBB2_concordant_repression": cls.get("ERBB2") == "concordant_repression",
        "ESR1_concordant_repression": cls.get("ESR1") == "concordant_repression",
    }
    summary = {
        "module": "mirna_hallmark.spatial_mibi_anchor",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "resolution": "single-cell PROTEIN (Risom DCIS MIBI)", "route": "measured protein (not miRNA)",
        "n_cells": int(len(df)), "n_tumour_cells": int(df["is_tumour"].sum()),
        "program_localization_stage_trend": loc.set_index("program")["stage_trend_rho"].to_dict(),
        "mh56d_classification": cls,
        "engine_validation_vs_MH56d": validation,
        "caveats": ["MIBI protein, NOT miRNA; DCIS not primary-invasive (in-situ anchor/bridge)",
                    "program markers single-protein = LOCALIZATION ONLY (MH-56 lesson: not pressure correspondence)",
                    "retest axis = normal→DCIS→IBC progression within epithelial/tumour-lineage cells (the healthy→tumour analogue of delta_tumor_gtex), NOT tumour-vs-stroma (epithelial-identity confound)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[mibi] program stage-trend (tumour): {summary['program_localization_stage_trend']}")
    print(f"[mibi] MH-56d spatial retest: {cls}")
    print(f"[mibi] engine validation vs MH-56d: {validation}")
    print(f"[mibi] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
