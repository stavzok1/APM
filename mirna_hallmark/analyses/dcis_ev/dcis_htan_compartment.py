"""Compartment test of the invasion miRNA programs in Risom 2022 DCIS MIBI (effector level).

Tests, in single-cell spatial proteomics (MIBI-TOF; 69,672 cells, matched normal/DCIS/IBC),
the two microenvironmental programs MH-48/MH-51 nominate as the LATE (invasion-coupled) axis:
  - **myoepithelial loss** (the miR-145/-143 program): myoepithelial-cell fraction + myoep
    markers (CK5/P63) across normal→DCIS→IBC.
  - **CAF / ECM desmoplasia** (the miR-29c program) — resolves the MH-51 open question:
      (ii) FRACTION shift : does the fibroblast/CAF fraction rise normal→DCIS→IBC?
      (i)  REGULATORY     : within fibroblasts, does the ECM program (COLI/SMA/FAP/VIM/MMP9)
                            intensify with stage, and are **CAF > NORMFIBRO** (activation)?
    Both ⇒ desmoplastic CAF activation (per-CAF regulatory change); only (ii) ⇒ pure fraction shift.

**Effector/compartment readout (MIBI protein) — NOT miR-29c itself** (HTAN/MIBI carry no miRNA;
WHATS_NEXT §5.5/§5.8). It tests the miR-29 *targets* (collagen-I etc.) + the CAF/myoepithelial
compartments — exactly what splits scenario (i) vs (ii). Cross-check vs the authors' own
`Tissue Feature Data/Table_S4_ECM_genes_CAF_high_low.csv` and `Table_S6_Myoepithelial_LDA`.

DATA: Risom single-cell table → ``data/external/risom2022_dcis_mibi/risom2022_single_cell.csv``
(Mendeley ``d87vg86zd8``; the file API is OAuth-gated so it was placed manually; see the dir's
``README.md`` for provenance + the tissue-feature side tables).

Run:
  .venv/bin/python3 -m mirna_hallmark.dcis_htan_compartment
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.dcis_timing import _jonckheere

DATA_DIR = C.REPO_ROOT / "data" / "external" / "risom2022_dcis_mibi"
OUT_DIR = C.OUTPUT_ROOT / "dcis_htan_compartment"
STAGE_MAP = {"normal": 0, "dcis": 1, "ibc": 2}        # Tissue_Type (drop 'lobular')
STAGE_NAME = {0: "normal", 1: "DCIS", 2: "IBC"}
ECM_MARKERS = ("COLI", "SMA", "FAP", "VIM", "MMP9")   # miR-29 effector axis (present-checked)
MYOEP_MARKERS = ("CK5", "P63")                        # myoepithelial integrity
FIBRO_PHENOS = ("CAF", "NORMFIBRO", "MYOFIBRO", "FIBRO_VIMonly")
MIN_CELLS = 20


def _find() -> Optional[Path]:
    for pat in ("risom2022_single_cell.csv", "Single_Cell_Data.csv", "*Single_Cell*.csv",
                "*cell*table*.csv", "*.csv"):
        hits = sorted(DATA_DIR.glob(pat))
        if hits:
            return hits[0]
    return None


def _pick(df: pd.DataFrame, prefer) -> Optional[str]:
    return next((c for c in prefer if c in df.columns), None)


def load_risom() -> Tuple[pd.DataFrame, dict]:
    path = _find()
    if path is None:
        raise FileNotFoundError(
            f"Risom single-cell table not in {DATA_DIR}. Download from Mendeley d87vg86zd8 "
            "(https://data.mendeley.com/datasets/d87vg86zd8 — file API OAuth-gated, manual) and "
            "place Single_Cell_Data.csv there.")
    df = pd.read_csv(path)
    cols = {"stage": _pick(df, ("Tissue_Type", "Tissue", "Status")),
            "lineage": _pick(df, ("celllineage", "manual_gating_cellineage", "sublineage")),
            "pheno": _pick(df, ("phenotype", "manual_gating_phenotype")),
            "sample": _pick(df, ("Point_Num", "point", "fov", "PointNumber"))}
    df["_stage"] = df[cols["stage"]].astype(str).str.lower().map(STAGE_MAP)
    return df, cols


def fraction_shift(df: pd.DataFrame, cols: dict) -> pd.DataFrame:
    """(ii) per-sample compartment fractions + CAF-activation ratio, by stage."""
    lin, ph, samp = cols["lineage"], cols["pheno"], cols["sample"]
    rows = []
    for sid, sub in df.dropna(subset=["_stage"]).groupby(samp):
        n = len(sub)
        nfib = sub[ph].isin(FIBRO_PHENOS).sum()
        rows.append({
            "sample": sid, "stage": int(sub["_stage"].mode().iloc[0]), "n_cells": n,
            "frac_fibroblast": float((sub[lin].astype(str).str.lower() == "fibroblast").mean()),
            "frac_myoep": float((sub[lin].astype(str).str.lower() == "myoep").mean()),
            "frac_CAF": float((sub[ph] == "CAF").mean()),
            # CAF as a share of all fibroblasts = activation, fraction-independent
            "CAF_activation_ratio": float((sub[ph] == "CAF").sum() / nfib) if nfib else np.nan,
        })
    return pd.DataFrame(rows)


def _trend(frame: pd.DataFrame, col: str) -> dict:
    f = frame.dropna(subset=[col, "stage"])
    if f["stage"].nunique() < 2 or len(f) < 6:
        return {"by_stage": {}, "jonckheere_p": None}
    jp = _jonckheere(f[col].to_numpy(float), f["stage"].to_numpy(float))[1]
    return {"by_stage": {STAGE_NAME[int(k)]: round(v, 3) for k, v in f.groupby("stage")[col].mean().items()},
            "jonckheere_p": round(float(jp), 4) if np.isfinite(jp) else None}


def within_fibroblast_ecm(df: pd.DataFrame, cols: dict) -> pd.DataFrame:
    """(i) within fibroblasts: ECM marker by stage (Kruskal+Jonckheere) + CAF-vs-NORMFIBRO MWU."""
    lin, ph = cols["lineage"], cols["pheno"]
    fib = df[(df[lin].astype(str).str.lower() == "fibroblast") & df["_stage"].notna()]
    caf, norm = fib[fib[ph] == "CAF"], fib[fib[ph] == "NORMFIBRO"]
    rows = []
    for m in (x for x in ECM_MARKERS if x in df.columns):
        bys = fib.dropna(subset=[m]).groupby("_stage")[m]
        groups = [g.to_numpy() for _, g in bys if len(g) >= MIN_CELLS]
        kw = kruskal(*groups)[1] if len(groups) >= 2 else np.nan
        sub = fib.dropna(subset=[m])
        jt = _jonckheere(sub[m].to_numpy(float), sub["_stage"].to_numpy(float))[1] if len(sub) > 30 else np.nan
        caf_p = (mannwhitneyu(caf[m].dropna(), norm[m].dropna(), alternative="greater")[1]
                 if len(caf) >= MIN_CELLS and len(norm) >= MIN_CELLS else np.nan)
        rows.append({"marker": m,
                     "means_by_stage": {STAGE_NAME[int(k)]: round(v, 3) for k, v in bys.mean().items()},
                     "kruskal_p_across_stage": round(float(kw), 5) if np.isfinite(kw) else None,
                     "jonckheere_p": round(float(jt), 5) if np.isfinite(jt) else None,
                     "CAF_gt_NORMFIBRO_mwu_p": round(float(caf_p), 5) if np.isfinite(caf_p) else None,
                     "CAF_mean": round(float(caf[m].mean()), 3) if len(caf) else None,
                     "NORMFIBRO_mean": round(float(norm[m].mean()), 3) if len(norm) else None})
    return pd.DataFrame(rows)


def myoepithelial_loss(df: pd.DataFrame, cols: dict) -> dict:
    """Myoepithelial-marker expression within the myoep compartment across stage (the loss program)."""
    lin = cols["lineage"]
    myo = df[(df[lin].astype(str).str.lower() == "myoep") & df["_stage"].notna()]
    out = {}
    for m in (x for x in MYOEP_MARKERS if x in df.columns):
        bys = myo.dropna(subset=[m]).groupby("_stage")[m]
        sub = myo.dropna(subset=[m])
        jt = _jonckheere(sub[m].to_numpy(float), sub["_stage"].to_numpy(float))[1] if len(sub) > 30 else np.nan
        out[m] = {"by_stage": {STAGE_NAME[int(k)]: round(v, 3) for k, v in bys.mean().items()},
                  "jonckheere_p": round(float(jt), 4) if np.isfinite(jt) else None}
    return out


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        df, cols = load_risom()
    except FileNotFoundError as e:
        print(f"[dcis_htan_compartment] PENDING DATA — {e}")
        return
    print(f"[dcis_htan_compartment] {len(df):,} cells; cols={cols}; "
          f"stages={df['_stage'].map(STAGE_NAME).value_counts(dropna=True).to_dict()}")

    frac = fraction_shift(df, cols)
    frac.to_csv(out_dir / "compartment_fraction_by_stage.tsv", sep="\t", index=False)
    ecm = within_fibroblast_ecm(df, cols)
    ecm.to_csv(out_dir / "fibroblast_ecm_by_stage.tsv", sep="\t", index=False)
    myo = myoepithelial_loss(df, cols)

    summ = {"module": "mirna_hallmark.dcis_htan_compartment",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "cohort": "Risom 2022 Cell DCIS MIBI-TOF (single-cell spatial proteomics)",
            "n_cells": int(len(df)), "n_samples": int(frac["sample"].nunique()),
            "fraction_shift_ii": {
                "fibroblast": _trend(frac, "frac_fibroblast"),
                "CAF": _trend(frac, "frac_CAF"),
                "CAF_activation_ratio": _trend(frac, "CAF_activation_ratio"),
                "myoepithelial": _trend(frac, "frac_myoep")},
            "within_CAF_ecm_i": ecm.to_dict("records"),
            "myoepithelial_loss": myo,
            "caveat": ("effector/compartment readout (MIBI protein) — tests miR-29 targets + CAF/"
                       "myoepithelial compartments, NOT miR-29c itself")}
    (out_dir / "method_manifest.json").write_text(json.dumps(summ, indent=2, default=str), encoding="utf-8")

    f2 = summ["fraction_shift_ii"]
    print(f"[ii] fibroblast frac {f2['fibroblast']['by_stage']} (JT p={f2['fibroblast']['jonckheere_p']})")
    print(f"[ii] CAF frac {f2['CAF']['by_stage']} (JT p={f2['CAF']['jonckheere_p']}); "
          f"CAF-activation ratio {f2['CAF_activation_ratio']['by_stage']} (JT p={f2['CAF_activation_ratio']['jonckheere_p']})")
    print(f"[ii] myoep frac {f2['myoepithelial']['by_stage']} (JT p={f2['myoepithelial']['jonckheere_p']})")
    print("[i] within-fibroblast ECM by stage:")
    for r in summ["within_CAF_ecm_i"]:
        print(f"    {r['marker']:5s} {r['means_by_stage']} JT p={r['jonckheere_p']}  "
              f"CAF>NORMFIBRO p={r['CAF_gt_NORMFIBRO_mwu_p']} (CAF {r['CAF_mean']} vs norm {r['NORMFIBRO_mean']})")
    print(f"[myoep loss] {myo}")
    print(f"[dcis_htan_compartment] wrote outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
