"""Roadmap Phase 3 outcome — famous miRNAs: head-to-head ABUNDANCE vs ANTI-CORR vs PRESSURE.

Answers "which readout prognoses best, and does anti-corr beat raw abundance?" For each famous miRNA,
three per-patient readouts over its OWN BY-neg target set:
  * abundance       = the miRNA's expression (input / potential)
  * anti_corr       = -mean z-expression of its targets (realized repression of its targets)
  * pressure        = mean framework gene-pressure on its targets (predicted repression)
Cox vs outcome (TCGA DFI+OS, Buffa DRFS); then a head-to-head: per (miRNA,endpoint) which readout has
the smallest p, tallied across miRNAs, + median -log10 p per readout.

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.outcome.outcome_famous_compare``
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome
from mirna_hallmark.analyses.outcome.outcome_famous_mirnas import FAMOUS, _edges_by_arm

OUT_DIR = C.OUTPUT_ROOT / "outcome_famous_compare"


def _readouts(coh, arm, targets):
    a = coh["arm_of"](arm)
    if a is None or a not in coh["mi"].index:
        return None
    tg = [g for g in targets if g in coh["rna"].index]
    if len(tg) < 5:
        return None
    sub = coh["rna"].loc[tg]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    out = {"abundance": coh["mi"].loc[a], "anti_corr": -z.mean(axis=0)}
    gtg = [g for g in tg if g in coh["gp"].index]
    if len(gtg) >= 5:
        out["pressure"] = coh["gp"].loc[gtg].mean(axis=0)
    return out


def _screen(coh) -> pd.DataFrame:
    eba = _edges_by_arm()
    rows = []
    for arm, role in FAMOUS.items():
        targets, _ = eba.get(arm, ([], None))
        ro = _readouts(coh, arm, targets)
        if ro is None:
            continue
        for metric, x in ro.items():
            for ep, (t, e) in coh["endpoints"].items():
                c, p = _cox_one(t, e, x.astype(float), coh["cov"])
                rows.append({"cohort": coh["name"], "miRNA": arm, "role": role, "metric": metric,
                             "endpoint": ep, "hr": np.exp(c), "p": p})
    r = pd.DataFrame(rows)
    for key, idx in r.groupby(["endpoint", "metric"]).groups.items():
        ok = r.loc[idx, "p"].notna(); r.loc[r.loc[idx].index[ok], "q"] = bh_fdr(r.loc[idx, "p"].dropna().to_numpy())
    return r


def _tcga():
    mi = D.load_mirna_arms(); mi = mi[~mi.index.duplicated(keep="first")]
    rna = D.load_rna(); rna = rna[~rna.index.duplicated(keep="first")]
    gp = pd.read_csv(C.OUTPUT_ROOT / "pressure_layer_comparison" / "L1_coupling_spine" / "gene_pressure_per_sample.tsv.gz", sep="\t", index_col=0)
    cl = load_tcga_outcome()
    return {"name": "TCGA", "mi": mi, "rna": rna, "gp": gp, "arm_of": lambda a: a if a in mi.index else None,
            "endpoints": {"DFI": (cl["dfi_t"], cl["dfi_e"]), "OS": (cl["os_t"], cl["os_e"])},
            "cov": pd.concat([cl["base"], cl["comp"]], axis=1)}


def _buffa():
    from mirna_hallmark.eval.buffa_validation import load_buffa
    from mirna_hallmark.analyses.outcome.outcome_buffa_pressure import _buffa_pressure, _clinical
    mi, rna = load_buffa()
    gpg, _, _ = _buffa_pressure()
    base = {re.sub(r"-(3p|5p)$", "", x): x for x in mi.index}
    clin = _clinical(mi.columns)
    return {"name": "Buffa", "mi": mi, "rna": rna, "gp": gpg,
            "arm_of": lambda a: base.get(re.sub(r"-(3p|5p)$", "", a)),
            "endpoints": {"DRFS": (clin["T"], clin["E"])}, "cov": _z(clin["er"]).to_frame("er")}


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    res = pd.concat([_screen(_tcga()), _screen(_buffa())], ignore_index=True)
    res.to_csv(out_dir / "famous_compare.tsv", sep="\t", index=False)

    # head-to-head: per (cohort,miRNA,endpoint) which readout has min p
    wide = res.pivot_table(index=["cohort", "miRNA", "endpoint"], columns="metric", values="p")
    metrics = [m for m in ("abundance", "anti_corr", "pressure") if m in wide.columns]
    win = wide[metrics].idxmin(axis=1).value_counts().to_dict()
    medlogp = {m: round(float(np.median(-np.log10(res[res.metric == m]["p"].dropna()))), 3) for m in metrics}
    nfdr = {m: int(((res.metric == m) & (res.q < 0.1)).sum()) for m in metrics}
    # paired: anti_corr vs abundance per (cohort,miRNA,endpoint) -log10p delta
    paired = wide.dropna(subset=["abundance", "anti_corr"]) if {"abundance", "anti_corr"} <= set(wide.columns) else pd.DataFrame()
    anti_better = int((paired["anti_corr"] < paired["abundance"]).sum()) if len(paired) else 0
    summary = {
        "module": "mirna_hallmark.analyses.outcome.outcome_famous_compare",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_FDR_by_metric": nfdr,
        "median_neglog10p_by_metric": medlogp,
        "head_to_head_wins(min p)": win,
        "anti_corr_beats_abundance_n/total": f"{anti_better}/{len(paired)}",
        "buffa_FDR_by_metric": {m: int(((res.cohort == "Buffa") & (res.metric == m) & (res.q < 0.1)).sum()) for m in metrics},
        "top_pressure_hits": res[(res.metric == "pressure") & (res.q < 0.2)].sort_values("p").head(8).apply(lambda r: f"{r['cohort']}/{r['miRNA']}/{r['endpoint']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "caveats": ["famous-miRNA pre-specified set; pressure = mean framework gene-pressure on the miRNA's targets",
                    "Buffa DRFS / TCGA DFI+OS; abundance vs anti_corr vs pressure head-to-head"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[compare] FDR hits by metric: {nfdr}")
    print(f"[compare] median -log10(p) by metric: {medlogp}  (higher = stronger)")
    print(f"[compare] head-to-head wins (lowest p): {win}")
    print(f"[compare] anti_corr beats abundance: {summary['anti_corr_beats_abundance_n/total']}")
    print(f"[compare] Buffa FDR by metric: {summary['buffa_FDR_by_metric']}")
    print(f"[compare] top pressure hits: {summary['top_pressure_hits']}")
    print(f"[compare] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
