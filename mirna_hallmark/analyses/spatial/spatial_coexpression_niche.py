"""Roadmap Phase 3 (spatial) — push #2: neighbourhood co-occurrence for the brake-release co-localization.

Finding 3 (`spatial_visium_wu`) showed released-pressure genes track their programs as a *global* across-
spot correlation (SLC2A1↔glycolysis, FAP↔EMT). The stronger spatial claim is that they occupy the SAME
PHYSICAL NEIGHBOURHOODS — high-gene spots sitting next to high-program spots. This module tests that with
**bivariate Moran's I** (kNN-weighted, gene vs program) and a location-permutation p-value, per primary
Visium section, pooled. A significant positive bivariate I means co-occurrence in space, not just
co-variation across the whole section.

Run: ``.venv/bin/python3 -m mirna_hallmark.spatial_coexpression_niche``
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
from mirna_hallmark.analyses.spatial.spatial_malignancy_retest import _load_visium
from mirna_hallmark.analyses.spatial.spatial_visium_wu import _samples

OUT_DIR = C.OUTPUT_ROOT / "spatial_coexpression_niche"
# released-pressure anchors and the program they serve (from Finding 3); + the gained discordant VIM/HIF1A
ANCHOR_PROGRAM = {"SLC2A1": "glycolysis", "FAP": "EMT", "VIM": "EMT", "HIF1A": "hypoxia"}


def analyse_section(sample: str, pg) -> dict:
    log, lin, pos = _load_visium(sample)
    coords = pos[["x", "y"]]
    prog = SP.score_programs(log, pg)
    out = {}
    for g, p in ANCHOR_PROGRAM.items():
        if g not in log.index or p not in prog.columns:
            continue
        gv, pv = log.loc[g], prog[p]
        uni = float(spearmanr(gv, pv)[0])                       # global across-spot correlation (Finding 3)
        biv_I, biv_p = SP.bivariate_morans_i(gv, pv, coords)    # neighbourhood co-occurrence
        out[g] = {"program": p, "univariate_rho": round(uni, 3), "bivariate_morans_I": biv_I, "perm_p": biv_p}
    return {"sample": sample, "n_spots": int(log.shape[1]), "pairs": out}


def run(*, out_dir: Path = OUT_DIR, samples: List[str] | None = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = samples or (_samples() + ["10x_public"])
    pg = SP.program_gene_sets()
    print(f"[niche] bivariate Moran's I over {len(samples)} sections")

    per = []
    for s in samples:
        try:
            r = analyse_section(s, pg)
            per.append(r)
        except Exception as e:
            print(f"   {s}: SKIPPED ({type(e).__name__}: {e})")

    rows = []
    for g, p in ANCHOR_PROGRAM.items():
        vals = [r["pairs"][g] for r in per if g in r["pairs"]]
        if not vals:
            continue
        Is = [v["bivariate_morans_I"] for v in vals if v["bivariate_morans_I"] == v["bivariate_morans_I"]]
        nsig = sum(v["perm_p"] is not None and v["perm_p"] < 0.05 for v in vals)
        rows.append({"gene": g, "program": p, "n_sections": len(vals),
                     "median_univariate_rho": round(float(np.median([v["univariate_rho"] for v in vals])), 3),
                     "median_bivariate_morans_I": round(float(np.median(Is)), 3) if Is else np.nan,
                     "n_sections_perm_p_lt05": int(nsig)})
    tab = pd.DataFrame(rows)
    tab.to_csv(out_dir / "neighbourhood_cooccurrence.tsv", sep="\t", index=False)
    pd.DataFrame([{"sample": r["sample"], "gene": g, **v}
                  for r in per for g, v in r["pairs"].items()]).to_csv(
        out_dir / "per_section_bivariate.tsv", sep="\t", index=False)

    summary = {
        "module": "mirna_hallmark.spatial_coexpression_niche",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "test": "bivariate Moran's I (kNN) + 199-perm location null, primary Visium",
        "n_sections": len(per),
        "neighbourhood_cooccurrence": tab.set_index("gene").to_dict("index") if not tab.empty else {},
        "interpretation": "positive bivariate I + perm_p<0.05 ⇒ the released gene and its program occupy the "
                          "same physical neighbourhoods (co-occurrence), strengthening Finding 3 beyond global correlation",
        "caveats": ["bivariate I is symmetric (co-occurrence, not direction); permutation shuffles program locations",
                    "measured mRNA only; programs are metagene scores"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[niche] neighbourhood co-occurrence:\n{tab.to_string(index=False)}")
    print(f"[niche] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, samples=args.samples)


if __name__ == "__main__":
    main()
