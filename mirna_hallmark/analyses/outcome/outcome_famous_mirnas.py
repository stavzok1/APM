"""Roadmap Phase 3 outcome — focused, hypothesis-driven: do FAMOUS breast-cancer miRNAs' pressure /
target-coupling relate to outcome?

The genome-wide sweeps were null partly from multiple-testing dilution. Here we pre-specify a curated set
of well-known oncomiRs / tumour-suppressor-miRs and test, per miRNA, two readouts vs survival in TCGA
(DFI+OS) and Buffa (DRFS), with small-set FDR:
  * EXPRESSION (the arm level)
  * TARGET-REPRESSION (per patient: -mean z-expression of the miRNA's BY-neg target genes = anti-corr
    activity on its own targets)
plus a descriptive panel (n BY-neg targets, median coupling rho, top Hallmark programs hit).

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.outcome.outcome_famous_mirnas``
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
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome

OUT_DIR = C.OUTPUT_ROOT / "outcome_famous_mirnas"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"

FAMOUS = {  # well-known breast-cancer miRNAs (arm-level) + role
    "hsa-miR-21-5p": "oncomiR", "hsa-miR-155-5p": "oncomiR", "hsa-miR-10b-5p": "oncomiR-metastasis",
    "hsa-miR-9-5p": "oncomiR", "hsa-miR-221-3p": "oncomiR", "hsa-miR-222-3p": "oncomiR",
    "hsa-miR-373-3p": "oncomiR", "hsa-miR-17-5p": "oncomiR(17-92)", "hsa-miR-18a-5p": "oncomiR(17-92)",
    "hsa-miR-19b-3p": "oncomiR(17-92)", "hsa-miR-210-3p": "oncomiR-hypoxia", "hsa-miR-103a-3p": "oncomiR",
    "hsa-miR-181a-5p": "oncomiR", "hsa-miR-93-5p": "oncomiR",
    "hsa-miR-200a-3p": "TSmiR-EMT", "hsa-miR-200b-3p": "TSmiR-EMT", "hsa-miR-200c-3p": "TSmiR-EMT",
    "hsa-miR-141-3p": "TSmiR-EMT", "hsa-miR-429": "TSmiR-EMT", "hsa-let-7a-5p": "TSmiR",
    "hsa-let-7b-5p": "TSmiR", "hsa-miR-205-5p": "TSmiR", "hsa-miR-34a-5p": "TSmiR-p53",
    "hsa-miR-145-5p": "TSmiR", "hsa-miR-143-3p": "TSmiR", "hsa-miR-31-5p": "TSmiR",
    "hsa-miR-126-3p": "TSmiR", "hsa-miR-146a-5p": "TSmiR", "hsa-miR-125b-5p": "TSmiR",
    "hsa-miR-335-5p": "TSmiR", "hsa-miR-30a-5p": "TSmiR", "hsa-miR-29c-3p": "TSmiR-ECM",
}


def _edges_by_arm() -> dict:
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()
    e["arm"] = e["key"].str.split(r"\|\|").str[0]; e["gene"] = e["key"].str.split(r"\|\|").str[1]
    return {a: (g["gene"].tolist(), float(g["rho"].median())) for a, g in e.groupby("arm")}


def _screen(coh: dict, edges_by_arm: dict) -> pd.DataFrame:
    rows = []
    for arm, role in FAMOUS.items():
        a = coh["arm_of"](arm)
        if a is None or a not in coh["mi"].index:
            continue
        targets, med_rho = edges_by_arm.get(arm, ([], np.nan))
        tg = [g for g in targets if g in coh["rna"].index]
        expr = coh["mi"].loc[a]
        if len(tg) >= 5:                                                        # target-repression: -mean z-expr over targets
            sub = coh["rna"].loc[tg]
            z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
            rep = -z.mean(axis=0)
        else:
            rep = None
        for ep, (t, e) in coh["endpoints"].items():
            ce, pe = _cox_one(t, e, expr, coh["cov"])
            rows.append({"cohort": coh["name"], "miRNA": arm, "role": role, "metric": "expression",
                         "endpoint": ep, "n_targets": len(tg), "hr": np.exp(ce), "p": pe})
            if rep is not None:
                cr, pr = _cox_one(t, e, rep, coh["cov"])
                rows.append({"cohort": coh["name"], "miRNA": arm, "role": role, "metric": "target_repression",
                             "endpoint": ep, "n_targets": len(tg), "hr": np.exp(cr), "p": pr})
    r = pd.DataFrame(rows)
    for key, idx in r.groupby(["endpoint", "metric"]).groups.items():
        ok = r.loc[idx, "p"].notna(); r.loc[r.loc[idx].index[ok], "q"] = bh_fdr(r.loc[idx, "p"].dropna().to_numpy())
    return r


def _tcga() -> dict:
    mi = D.load_mirna_arms(); mi = mi[~mi.index.duplicated(keep="first")]
    rna = D.load_rna(); rna = rna[~rna.index.duplicated(keep="first")]
    cl = load_tcga_outcome()
    return {"name": "TCGA", "mi": mi, "rna": rna, "arm_of": lambda a: a if a in mi.index else None,
            "endpoints": {"DFI": (cl["dfi_t"], cl["dfi_e"]), "OS": (cl["os_t"], cl["os_e"])},
            "cov": pd.concat([cl["base"], cl["comp"]], axis=1)}


def _buffa() -> dict:
    from mirna_hallmark.eval.buffa_validation import load_buffa
    from mirna_hallmark.analyses.outcome.outcome_buffa_pressure import _clinical
    mi, rna = load_buffa()
    base = {re.sub(r"-(3p|5p)$", "", x): x for x in mi.index}
    clin = _clinical(mi.columns)
    return {"name": "Buffa", "mi": mi, "rna": rna,
            "arm_of": lambda a: base.get(re.sub(r"-(3p|5p)$", "", a)),
            "endpoints": {"DRFS": (clin["T"], clin["E"])}, "cov": _z(clin["er"]).to_frame("er")}


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    eba = _edges_by_arm()
    res = pd.concat([_screen(_tcga(), eba), _screen(_buffa(), eba)], ignore_index=True)
    res.to_csv(out_dir / "famous_mirna_survival.tsv", sep="\t", index=False)
    # descriptive: programs each famous miRNA's targets hit
    hs = HallmarkSets.load()
    desc = []
    for arm, role in FAMOUS.items():
        tg, med = eba.get(arm, ([], np.nan))
        progs = sorted({h.replace("HALLMARK_", "") for h, m in hs.sets.items() for g in tg if g in m}, key=lambda x: x)[:4] if tg else []
        desc.append({"miRNA": arm, "role": role, "n_BYneg_targets": len(tg), "median_rho": round(med, 3) if tg else None, "top_programs": progs})
    pd.DataFrame(desc).to_csv(out_dir / "famous_mirna_descriptive.tsv", sep="\t", index=False)

    sig = res[res["p"] < 0.05].sort_values("p")
    fdr = res[res["q"] < 0.1].sort_values("q")
    summary = {
        "module": "mirna_hallmark.analyses.outcome.outcome_famous_mirnas",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_famous_tested": int(res["miRNA"].nunique()),
        "FDR_hits": fdr.apply(lambda r: f"{r['cohort']}/{r['miRNA']}({r['role']})/{r['metric']}/{r['endpoint']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "nominal_hits": sig.head(20).apply(lambda r: f"{r['cohort']}/{r['miRNA']}({r['role']})/{r['metric']}/{r['endpoint']} HR={r['hr']:.2f} p={r['p']:.3f}", axis=1).tolist(),
        "cross_cohort_concordant": [],
        "caveats": ["pre-specified famous set (~30) -> small-set FDR; expression + target-repression; both cohorts",
                    "Buffa base-name match, DRFS only; TCGA DFI+OS"],
    }
    # cross-cohort: same miRNA+metric nominal in TCGA(any ep) and Buffa DRFS, same direction
    for mirna in res["miRNA"].unique():
        for metric in ("expression", "target_repression"):
            tc = res[(res.cohort == "TCGA") & (res.miRNA == mirna) & (res.metric == metric) & (res.p < 0.05)]
            bf = res[(res.cohort == "Buffa") & (res.miRNA == mirna) & (res.metric == metric) & (res.p < 0.05)]
            if len(tc) and len(bf) and np.sign(np.log(tc.hr.iloc[0])) == np.sign(np.log(bf.hr.iloc[0])):
                summary["cross_cohort_concordant"].append(f"{mirna}/{metric}: TCGA {tc.endpoint.iloc[0]} HR{tc.hr.iloc[0]:.2f} / Buffa DRFS HR{bf.hr.iloc[0]:.2f}")
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[famous] FDR hits: {summary['FDR_hits']}")
    print(f"[famous] nominal: {summary['nominal_hits'][:12]}")
    print(f"[famous] cross-cohort concordant: {summary['cross_cohort_concordant']}")
    print(f"[famous] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
