"""LCM-microdissected breast progression (GSE162670) — composition-reduced replication.

GSE162670 is laser-capture-microdissected (LCM) miRNA by tissue **morphology** (TNBC FFPE):
normal ducts/lobules -> DCIS -> invasive (many TNBC morphologies). Because the material is
**microdissected**, stromal/immune contamination is physically reduced relative to bulk
tissue, so this is the tie-breaker for the composition question raised by the partial-
confounder adjustment (`dcis_timing`/`dcis_mrna_coupling`):

  - In *bulk* human DCIS tissue the miR-29c grade/coupling signal was **composition-driven**
    (collapsed under epi/immune/stroma adjustment).
  - In the *pure* MCF10A cell line (GSE93740) miR-29c declined cell-intrinsically (ρ≈−0.88).
  - Here, in *LCM-purified human epithelium*, does miR-29c (and the loss-leaders) still
    decline normal->DCIS->invasive? If yes -> epithelial-intrinsic; if no -> composition.

NB the planned epithelial-vs-stroma *compartment* labels are NOT present (the series is
morphology-resolved, not compartment-resolved), so per the plan contingency this is used as
a composition-REDUCED progression-replication cohort, not a compartment split.

Run:
  .venv/bin/python3 -m mirna_hallmark.dcis_lcm_compartment_gse162670
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.dcis_timing import _jonckheere, acquired_direction
from mirna_hallmark.stats import bh_fdr

GSE162670 = C.REPO_ROOT / "data" / "external" / "dcis_geo" / "GSE162670_RMA_normalized_matrix.txt.gz"
DE_LANDSCAPE = (C.TISSUE_REFERENCE_DIR / "cross_state_landscape" / "mirna_arm_de_landscape.tsv")
OUT_DIR = C.OUTPUT_ROOT / "dcis_lcm_gse162670"
LOSS_LEADERS = ("miR-145", "miR-126", "miR-29c", "miR-140", "miR-497", "miR-143")
_STAGE_NAME = {0: "normal_LCM", 1: "DCIS_LCM", 2: "invasive_LCM"}


def _morph_stage(col: str) -> float:
    """Map LCM sample name (e.g. 1N / 3T2_cis+lymph / 1T1_apo) -> morphology ordinal."""
    c = col.strip().lower()
    if re.match(r"^\d+n$", c):
        return 0.0                       # normal ducts/lobules
    if "cis" in c:
        return 1.0                       # ductal carcinoma in situ
    if "t" in c:
        return 2.0                       # invasive TNBC morphology
    return np.nan


def load_lcm() -> Tuple[pd.DataFrame, pd.Series]:
    """(LCM miRNA RMA-log2 x sample, morphology stage 0/1/2)."""
    with gzip.open(GSE162670, "rt") as fh:
        lines = fh.readlines()
    h = next(i for i, l in enumerate(lines) if l.startswith("ID_REF\t"))
    df = pd.read_csv(io.StringIO("".join(lines[h:])), sep="\t", index_col=0)
    df = df.loc[df.index.astype(str).str.startswith("hsa-")].apply(pd.to_numeric, errors="coerce")
    stage = pd.Series({c: _morph_stage(c) for c in df.columns}).dropna()
    return df[stage.index], stage


def lcm_trend(df: pd.DataFrame, stage: pd.Series) -> pd.DataFrame:
    """Per miRNA: Spearman vs morphology ordinal + Jonckheere, in LCM-purified tissue."""
    s = stage.reindex(df.columns).to_numpy(float)
    rows = []
    for arm in df.index:
        y = df.loc[arm].to_numpy(float)
        m = ~np.isnan(y)
        if m.sum() < 8 or np.unique(y[m]).size < 2:
            continue
        rho, p = spearmanr(y[m], s[m])
        jz, jp = _jonckheere(y[m], s[m])
        rows.append({"arm": arm, "lcm_trend_rho": rho, "lcm_trend_p": p, "lcm_jt_p": jp,
                     "n": int(m.sum())})
    res = pd.DataFrame(rows)
    res["lcm_trend_q"] = bh_fdr(res["lcm_trend_p"].to_numpy())
    acq = acquired_direction(pd.read_csv(DE_LANDSCAPE, sep="\t"), None).rename(
        columns={"modern_arm": "arm"})[["arm", "acq_dir"]]
    res = res.merge(acq, on="arm", how="left")
    res["loss_leader"] = res["arm"].str.contains("|".join(LOSS_LEADERS))
    return res.sort_values("lcm_trend_rho").reset_index(drop=True)


def _summary(t: pd.DataFrame, stage: pd.Series) -> dict:
    m29 = t[t["arm"].str.contains("miR-29")]
    ll = t[t["loss_leader"]]
    # does the LCM progression agree with the TCGA acquired direction (composition-reduced)?
    ax = t[t["acq_dir"].abs() > 0]
    conc = float((np.sign(ax["lcm_trend_rho"]) == np.sign(ax["acq_dir"])).mean()) if len(ax) else np.nan
    return {
        "n_mirna": int(len(t)),
        "stage_counts": {_STAGE_NAME[k]: int((stage == k).sum()) for k in sorted(stage.unique())},
        "miR29_LCM": [
            {"arm": r.arm, "lcm_trend_rho": round(r.lcm_trend_rho, 3),
             "lcm_trend_q": round(r.lcm_trend_q, 4), "lcm_jt_p": round(r.lcm_jt_p, 4)}
            for r in m29.sort_values("lcm_trend_rho").itertuples()],
        "loss_leaders_LCM": [
            {"arm": r.arm, "lcm_trend_rho": round(r.lcm_trend_rho, 3), "lcm_trend_q": round(r.lcm_trend_q, 4)}
            for r in ll.sort_values("lcm_trend_rho").itertuples()],
        "n_loss_leaders_declining_q<0.1": int(((ll["lcm_trend_rho"] < 0) & (ll["lcm_trend_q"] < 0.1)).sum()),
        "lcm_vs_tcga_acquired_direction_concordance": round(conc, 3),
        "verdict_miR29c": ("epithelial-intrinsic decline (survives microdissection)"
                           if (m29[m29["arm"].str.contains("miR-29c")]["lcm_trend_rho"] < 0).any()
                           and (m29[m29["arm"].str.contains("miR-29c")]["lcm_trend_q"] < 0.1).any()
                           else "not significant in LCM -> consistent with composition-driven in bulk"),
        "caveats": [
            "LCM morphology-resolved (not epi/stroma compartment) -> composition REDUCED, not eliminated",
            "TNBC-only FFPE; Affymetrix RMA; cross-platform -> rank/direction only",
            "invasive = pooled TNBC morphologies; normal = tumor-adjacent ducts/lobules (field effect possible)",
        ],
    }


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[dcis_lcm] loading GSE162670 LCM morphology progression ...")
    df, stage = load_lcm()
    print(f"[dcis_lcm] {df.shape[0]} hsa miRNA x {df.shape[1]} LCM samples; stages "
          f"{ {_STAGE_NAME[k]: int((stage==k).sum()) for k in sorted(stage.unique())} }")
    t = lcm_trend(df, stage)
    t.to_csv(out_dir / "lcm_progression_trend.tsv", sep="\t", index=False)

    summary = _summary(t, stage)
    manifest = {"module": "mirna_hallmark.dcis_lcm_compartment_gse162670",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "cohort": "GSE162670 (LCM miRNA by morphology, TNBC; normal->DCIS->invasive)", **summary}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("[dcis_lcm] miR-29 LCM trend:")
    for r in summary["miR29_LCM"][:6]:
        print(f"   {r['arm']:18s} ρ={r['lcm_trend_rho']:+.2f} q={r['lcm_trend_q']:.3f} jt_p={r['lcm_jt_p']:.3f}")
    print(f"[dcis_lcm] miR-29c verdict: {summary['verdict_miR29c']}")
    print(f"[dcis_lcm] loss-leaders declining (q<0.1): {summary['n_loss_leaders_declining_q<0.1']}; "
          f"LCM-vs-TCGA acquired concordance {summary['lcm_vs_tcga_acquired_direction_concordance']}")
    print(f"[dcis_lcm] wrote outputs under {out_dir}")
    return t


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
