"""Leverage Risom MIBI further: DCIS PROGRESSION OUTCOME + SPATIAL NICHE of the invasion programs.

Beyond MH-52 (compartment events across normal→DCIS→IBC), the Risom single-cell table also carries
**outcome** (Status: case = DCIS that later progressed to ipsilateral invasion; ctrl = DCIS that did
not) and **spatial** info (per-cell neighbour distances + tissue masks). Two new questions:

  1. OUTCOME — do the invasion-program metrics (myoepithelial fraction, CAF-activation ratio,
     within-fibroblast collagen) differ between progressor (case) and non-progressor (ctrl) DCIS?
     (Risom's headline: myoepithelial disruption was *more* advanced in non-progressors — i.e.
     counter-intuitively protective; we test our metrics against that.)
  2. SPATIAL NICHE — is myoepithelial loss spatially coupled to CAF activation / ECM? i.e. are
     activated (FAP⁺/CAF) fibroblasts in myoepithelium-poor neighbourhoods (larger Neighbor_dist_MYOEP),
     and is collagen higher in fibroblasts far from myoepithelium?

Effector/compartment readout (MIBI protein), not miRNA (apm-rigor-protocol move E).
Run: ``.venv/bin/python3 -m mirna_hallmark.dcis_risom_outcome_niche``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from mirna_hallmark import config as C

CELLTABLE = C.REPO_ROOT / "data" / "external" / "risom2022_dcis_mibi" / "risom2022_single_cell.csv"
OUT_DIR = C.OUTPUT_ROOT / "dcis_risom_outcome_niche"
FIBRO_PHENOS = ("CAF", "NORMFIBRO", "MYOFIBRO", "FIBRO_VIMonly")


def _mwu(a, b, alt="two-sided"):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 4 or len(b) < 4:
        return np.nan, np.nan, np.nan
    return float(np.median(a)), float(np.median(b)), float(mannwhitneyu(a, b, alternative=alt)[1])


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CELLTABLE)
    lin, ph = "celllineage", "phenotype"
    dcis = df[df["Tissue_Type"] == "DCIS"].copy()
    case = dcis[dcis["Status"] == "case"]; ctrl = dcis[dcis["Status"] == "ctrl"]
    print(f"[risom_outcome] DCIS: case(progressor) {len(case):,} cells / ctrl(non) {len(ctrl):,}")

    # ---- 1. OUTCOME: per-FOV compartment metrics, case vs ctrl ----
    def per_fov(sub):
        rows = []
        for fov, s in sub.groupby("Point_Num"):
            nfib = s[ph].isin(FIBRO_PHENOS).sum()
            fib = s[s[lin].astype(str).str.lower() == "fibroblast"]
            rows.append({"fov": fov, "status": s["Status"].iloc[0],
                         "frac_myoep": float((s[lin].astype(str).str.lower() == "myoep").mean()),
                         "caf_activation": float((s[ph] == "CAF").sum() / nfib) if nfib else np.nan,
                         "fib_COLI": float(fib["COLI"].mean()) if len(fib) else np.nan})
        return pd.DataFrame(rows)
    fov = pd.concat([per_fov(case), per_fov(ctrl)], ignore_index=True)
    fov.to_csv(out_dir / "dcis_outcome_per_fov.tsv", sep="\t", index=False)
    outcome = {}
    for metric in ("frac_myoep", "caf_activation", "fib_COLI"):
        c = fov[fov.status == "case"][metric]; n = fov[fov.status == "ctrl"][metric]
        mc, mn, p = _mwu(c, n)
        outcome[metric] = {"case_median": round(mc, 3) if mc == mc else None,
                           "ctrl_median": round(mn, 3) if mn == mn else None,
                           "mwu_p": round(p, 4) if p == p else None}

    # ---- 2. SPATIAL NICHE: myoep loss <-> CAF activation / ECM ----
    fibd = dcis[(dcis[lin].astype(str).str.lower() == "fibroblast") & dcis["Neighbor_dist_MYOEP"].notna()]
    niche = {}
    if len(fibd) > 50:
        caf = fibd[fibd[ph] == "CAF"]; norm = fibd[fibd[ph] == "NORMFIBRO"]
        # activated CAFs farther from myoepithelium? (myoep-poor niche)
        mc, mn, p = _mwu(caf["Neighbor_dist_MYOEP"], norm["Neighbor_dist_MYOEP"], alt="greater")
        niche["CAF_farther_from_myoep_than_NORMFIBRO"] = {"caf_med_dist": round(mc, 1) if mc == mc else None,
            "normfibro_med_dist": round(mn, 1) if mn == mn else None, "mwu_p_greater": round(p, 4) if p == p else None}
        # collagen higher in fibroblasts far from myoepithelium?
        rho, pr = spearmanr(fibd["Neighbor_dist_MYOEP"], fibd["COLI"])
        niche["fibroblast_COLI_vs_dist_to_myoep_spearman"] = {"rho": round(float(rho), 3), "p": float(f"{pr:.2e}")}
        far = fibd[fibd["Neighbor_dist_MYOEP"] > fibd["Neighbor_dist_MYOEP"].median()]
        near = fibd[fibd["Neighbor_dist_MYOEP"] <= fibd["Neighbor_dist_MYOEP"].median()]
        mc, mn, p = _mwu(far["COLI"], near["COLI"], alt="greater")
        niche["COLI_far_vs_near_myoep"] = {"far_med": round(mc, 3) if mc == mc else None,
            "near_med": round(mn, 3) if mn == mn else None, "mwu_p_greater": round(p, 4) if p == p else None}

    summary = {"module": "mirna_hallmark.dcis_risom_outcome_niche",
               "generated_utc": datetime.now(timezone.utc).isoformat(),
               "n_case_fov": int(fov[fov.status == "case"]["fov"].nunique()),
               "n_ctrl_fov": int(fov[fov.status == "ctrl"]["fov"].nunique()),
               "outcome_case_vs_ctrl": outcome,
               "spatial_niche": niche,
               "caveat": "MIBI protein/compartment (not miRNA); case/ctrl = DCIS progressor/non; per-FOV n modest"}
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[risom_outcome] case/ctrl FOVs: {summary['n_case_fov']}/{summary['n_ctrl_fov']}")
    for m, v in outcome.items():
        print(f"   {m:16s} case {v['case_median']} vs ctrl {v['ctrl_median']}  p={v['mwu_p']}")
    print("[risom_niche]")
    for k, v in niche.items():
        print(f"   {k}: {v}")
    print(f"[risom_outcome_niche] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
