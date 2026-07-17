"""Cellular MCF10A isogenic progression (GSE93740) — miR-29c validation + DIRECT export test.

GSE93740 profiles the MCF10A isogenic series as **cellular** miRNA across
P(normal) -> AT1(hyperplasia) -> NeoT(neoplastic) -> DCIS -> CA1d/CA1h(invasive),
3 reps/stage. Two uses:

1. **Independent validation** of the MH-48/MH-49 miR-29c (and loss-leader) decline: does
   miR-29c fall along the *cellular* progression in this separate isogenic cohort?
2. **DIRECT export test (resolves the MH-50 open caveat).** GSE93740 is the *cellular*
   counterpart of GSE297448's *EVs* from the **same MCF10A lineage**, so we no longer have
   to borrow the cellular direction from TCGA: per miRNA we cross the cellular trend
   (GSE93740) against the compositional EV-share trend (GSE297448). cellular-DOWN + EV-UP =
   *directly* exported (the cell sheds the suppressor); cellular-UP + EV-UP = co-accumulated.

Isogenic -> the ordered trajectory is legitimate. Small n (3/stage) -> direction, not power.

Run:
  .venv/bin/python3 -m mirna_hallmark.dcis_mir29c_gse93740
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.dcis_ev_progression import _clr, load_ev, progression_trend
from mirna_hallmark.analyses.dcis_ev.dcis_timing import _jonckheere
from mirna_hallmark.stats import bh_fdr

GSE93740 = C.REPO_ROOT / "data" / "external" / "dcis_geo" / "GSE93740_bedrosian-mirna.exact4.txt.gz"
OUT_DIR = C.OUTPUT_ROOT / "plasma_ev" / "mcf10a_cellular"
LOSS_LEADERS = ("miR-145", "miR-126", "miR-29c", "miR-140", "miR-497", "miR-143")
# native cellular stage ordinal (finer than the EV series; NeoT is the extra step)
_STAGE_NAME = {0: "normal_P", 1: "hyperplasia_AT1", 2: "neoplastic_NeoT", 3: "DCIS", 4: "invasive_CA1"}


def _stage_of(col: str) -> float:
    c = col.lower()
    if "ca1" in c:
        return 4.0
    if "dcis" in c:
        return 3.0
    if "neot" in c:
        return 2.0
    if "at1" in c:
        return 1.0
    if "-p-" in c or "(p)" in c:
        return 0.0
    return np.nan


def load_cellular() -> Tuple[pd.DataFrame, pd.Series]:
    """(cellular miRNA log2 x sample, stage ordinal 0..4)."""
    df = pd.read_csv(GSE93740, sep="\t", index_col=0)
    stage = pd.Series({c: _stage_of(c) for c in df.columns}).dropna()
    df = df[stage.index]
    ev = np.log2(df.astype(float).clip(lower=0) + 1.0) if (df.values > 50).any() else df.astype(float)
    return ev, stage


def cellular_trend(ev: pd.DataFrame, stage: pd.Series) -> pd.DataFrame:
    """Per miRNA: cellular Spearman vs stage ordinal + Jonckheere."""
    s = stage.reindex(ev.columns).to_numpy(float)
    rows = []
    for arm in ev.index:
        y = ev.loc[arm].to_numpy(float)
        if np.unique(y[~np.isnan(y)]).size < 2:
            continue
        rho, p = spearmanr(y, s, nan_policy="omit")
        jz, jp = _jonckheere(y[~np.isnan(y)], s[~np.isnan(y)])
        rows.append({"arm": arm, "cell_trend_rho": rho, "cell_trend_p": p, "cell_jt_p": jp})
    res = pd.DataFrame(rows)
    res["cell_trend_q"] = bh_fdr(res["cell_trend_p"].to_numpy())
    return res


def direct_export(cell: pd.DataFrame) -> pd.DataFrame:
    """Cross the cellular trend (GSE93740) with the EV-share trend (GSE297448), SAME lineage."""
    ev, evstage = load_ev()
    evt = progression_trend(_clr(ev), evstage, raw=ev)[["arm", "ev_trend_rho", "ev_trend_q"]]
    m = cell.merge(evt, on="arm", how="inner")

    def _cls(r):
        cu, cd = r["cell_trend_rho"] > 0, r["cell_trend_rho"] < 0
        eu = r["ev_trend_rho"] > 0
        if cd and eu:
            return "DIRECTLY_EXPORTED (cell-down, EV-up)"
        if cu and eu:
            return "co-accumulated (cell-up, EV-up)"
        if cd and not eu:
            return "depleted (cell-down, EV-down)"
        return "retained (cell-up, EV-down)"

    m["direct_class"] = m.apply(_cls, axis=1)
    m["loss_leader"] = m["arm"].str.contains("|".join(LOSS_LEADERS))
    return m.sort_values("cell_trend_rho").reset_index(drop=True)


def _summary(cell: pd.DataFrame, dx: pd.DataFrame) -> dict:
    ll = cell[cell["arm"].str.contains("|".join(LOSS_LEADERS))]
    m29 = cell[cell["arm"].str.contains("miR-29")]
    # export anti-correlation: across miRNAs, does cellular trend anti-correlate with EV trend?
    rho_ce, p_ce = spearmanr(dx["cell_trend_rho"], dx["ev_trend_rho"])
    return {
        "n_mirna": int(len(cell)),
        "stages": _STAGE_NAME,
        "miR29_cellular": [
            {"arm": r.arm, "cell_trend_rho": round(r.cell_trend_rho, 3),
             "cell_trend_q": round(r.cell_trend_q, 4), "cell_jt_p": round(r.cell_jt_p, 4)}
            for r in m29.sort_values("cell_trend_rho").itertuples()],
        "loss_leaders_cellular": [
            {"arm": r.arm, "cell_trend_rho": round(r.cell_trend_rho, 3),
             "cell_trend_q": round(r.cell_trend_q, 4)}
            for r in ll.sort_values("cell_trend_rho").itertuples()],
        "n_loss_leaders_declining_q<0.1": int(((ll["cell_trend_rho"] < 0) & (ll["cell_trend_q"] < 0.1)).sum()),
        "DIRECT_export": {
            "cellular_vs_EV_trend_spearman": round(float(rho_ce), 3),
            "cellular_vs_EV_trend_p": float(f"{p_ce:.3g}"),
            "class_counts": dx["direct_class"].value_counts().to_dict(),
            "directly_exported_loss_leaders": [
                r.arm for r in dx[dx["loss_leader"]
                                  & dx["direct_class"].str.startswith("DIRECTLY")].itertuples()],
        },
        "caveats": [
            "isogenic MCF10A cellular EVs/cells, n=3/stage -> mechanistic direction, not population",
            "GSE93740 (cellular) + GSE297448 (EV) are the SAME lineage -> direct export (no TCGA borrow)",
            "EV side is compositional (CLR share); a negative cellular-vs-EV trend correlation = selective export",
            "sublines differ slightly (CA1d/CA1h vs CA1a) but share the progression",
        ],
    }


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[mcf10a_cellular] loading GSE93740 cellular MCF10A progression ...")
    ev, stage = load_cellular()
    print(f"[mcf10a_cellular] {ev.shape[0]} miRNA x {ev.shape[1]} samples; stages "
          f"{ {(_STAGE_NAME[s]): int((stage==s).sum()) for s in sorted(stage.unique())} }")
    cell = cellular_trend(ev, stage)
    dx = direct_export(cell)
    dx.to_csv(out_dir / "cellular_vs_ev_direct_export.tsv", sep="\t", index=False)
    cell.to_csv(out_dir / "cellular_progression_trend.tsv", sep="\t", index=False)

    summary = _summary(cell, dx)
    manifest = {"module": "mirna_hallmark.dcis_mir29c_gse93740",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "cohort": "GSE93740 cellular MCF10A (P->AT1->NeoT->DCIS->CA1) + GSE297448 EV", **summary}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("[mcf10a_cellular] miR-29 cellular trend:")
    for r in summary["miR29_cellular"]:
        print(f"   {r['arm']:18s} cell ρ={r['cell_trend_rho']:+.2f} q={r['cell_trend_q']:.3f} jt_p={r['cell_jt_p']:.3f}")
    de = summary["DIRECT_export"]
    print(f"[mcf10a_cellular] DIRECT export: cellular-vs-EV trend ρ={de['cellular_vs_EV_trend_spearman']} "
          f"(p={de['cellular_vs_EV_trend_p']}); classes {de['class_counts']}")
    print(f"[mcf10a_cellular] directly-exported loss-leaders: {de['directly_exported_loss_leaders']}")
    print(f"[mcf10a_cellular] wrote outputs under {out_dir}")
    return dx


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
