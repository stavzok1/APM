"""Roadmap Phase 3 outcome — miR-200/EMT deep-dive + the pressure↔expression anti-corr readout +
target-set modelling by gene role / cancer-relevance.

Three things at once:
  (A) miR-200 family (200a/b/c-3p, 141-3p, 429) deep-dive vs DRFS/OS, resolved to individual TARGETS
      (which drive it — ZEB1/2?).
  (B) the missing readout — **pe_coherence** = per-patient anti-corr of PRESSURE vs EXPRESSION over the
      miRNA's targets (does its pressure actually suppress its targets in this patient).
  (C) **target-set by role / cancer-relevance**: compute every readout not on ALL targets but on subsets
      {all, oncogene, TSG, TF, high-dependency(DepMap)} — biology-weighted, not flat.

Readouts: abundance · anti_corr(expr) · pressure · pe_coherence · escape, each × target-subset.
Cox vs TCGA DFI+OS + Buffa DRFS, FDR within (cohort,endpoint). Context miRNAs: let-7a, miR-30a, miR-210.

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_mir200_deepdive``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.gene_roles import load_gene_roles, load_gene_dependency
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z
from mirna_hallmark.analyses.outcome.outcome_famous_mirnas import _edges_by_arm
from mirna_hallmark.analyses.outcome.outcome_famous_compare import _tcga, _buffa

OUT_DIR = C.OUTPUT_ROOT / "outcome_mir200_deepdive"
FOCUS = {"hsa-miR-200a-3p": "miR-200", "hsa-miR-200b-3p": "miR-200", "hsa-miR-200c-3p": "miR-200",
         "hsa-miR-141-3p": "miR-200", "hsa-miR-429": "miR-200",
         "hsa-let-7a-5p": "TS", "hsa-miR-30a-5p": "TS", "hsa-miR-210-3p": "onco"}
TF_TSV = C.REPO_ROOT / "annotations" / "humantfs_lambert2018_tf.tsv"


def _zrows(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def _role_sets():
    gr = load_gene_roles().set_index("gene")
    tf = set(pd.read_csv(TF_TSV, sep="\t").iloc[:, 0].astype(str))
    dep = load_gene_dependency().set_index("gene")["dep_role_weight"]
    hidep = set(dep[dep.abs() > dep.abs().quantile(0.75)].index)
    return {"oncogene": set(gr[gr.malignancy_sign == 1].index), "TSG": set(gr[gr.malignancy_sign == -1].index),
            "TF": tf, "highdep": hidep}


def _readouts(coh, arm, tg):
    a = coh["arm_of"](arm)
    if a is None or a not in coh["mi"].index or len(tg) < 5:
        return {}
    ez = _zrows(coh["rna"].loc[tg])
    out = {"abundance": _z(coh["mi"].loc[a]), "anti_corr": -ez.mean(axis=0)}
    gtg = [g for g in tg if g in coh["gp"].index]
    if len(gtg) >= 5:
        common = coh["gp"].columns.intersection(ez.columns)            # pressure & RNA sample sets differ
        pz = _zrows(coh["gp"].loc[gtg, common]); ezg = ez.loc[gtg, common]
        out["pressure"] = coh["gp"].loc[gtg, common].mean(axis=0)
        A, B = pz.to_numpy(), ezg.to_numpy()
        Am, Bm = A - np.nanmean(A, 0), B - np.nanmean(B, 0)
        den = np.sqrt(np.nansum(Am ** 2, 0) * np.nansum(Bm ** 2, 0))
        out["pe_coherence"] = pd.Series(-np.nansum(Am * Bm, 0) / np.where(den == 0, np.nan, den), index=pz.columns)
    df = pd.concat([out["abundance"].rename("ab"), out["anti_corr"].rename("an")], axis=1).dropna()
    if len(df) >= 30:
        b = np.polyfit(df["ab"], df["an"], 1)
        out["escape"] = -(out["anti_corr"] - (b[0] * out["abundance"] + b[1]))
    return out


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    eba = _edges_by_arm(); roles = _role_sets()
    rows, target_rows = [], []
    for loader in (_tcga, _buffa):
        coh = loader()
        for arm, fam in FOCUS.items():
            targets, _ = eba.get(arm, ([], None))
            tg_all = [g for g in targets if g in coh["rna"].index]
            subsets = {"all": tg_all, **{r: [g for g in tg_all if g in gs] for r, gs in roles.items()}}
            for sname, tg in subsets.items():
                for metric, x in _readouts(coh, arm, tg).items():
                    for ep, (t, e) in coh["endpoints"].items():
                        c, p = _cox_one(t, e, x.astype(float), coh["cov"])
                        rows.append({"cohort": coh["name"], "miRNA": arm, "family": fam, "subset": sname,
                                     "n_targets": len(tg), "metric": metric, "endpoint": ep, "hr": np.exp(c), "p": p})
        # per-target: which individual miR-200 targets drive DRFS/OS (expression)
        m200_tg = sorted({g for a in FOCUS if FOCUS[a] == "miR-200" for g in eba.get(a, ([], None))[0] if g in coh["rna"].index})
        ez = _zrows(coh["rna"].loc[m200_tg])
        for ep, (t, e) in coh["endpoints"].items():
            for g in m200_tg:
                c, p = _cox_one(t, e, ez.loc[g], coh["cov"])
                target_rows.append({"cohort": coh["name"], "target": g, "endpoint": ep, "hr": np.exp(c), "p": p})
        print(f"[m200] {coh['name']} done ({len(m200_tg)} miR-200 targets)")
    res = pd.DataFrame(rows)
    for key, idx in res.groupby(["cohort", "endpoint"]).groups.items():
        ok = res.loc[idx, "p"].notna(); res.loc[res.loc[idx].index[ok], "q"] = bh_fdr(res.loc[idx, "p"].dropna().to_numpy())
    res.to_csv(out_dir / "mir200_readouts.tsv", sep="\t", index=False)
    tgt = pd.DataFrame(target_rows)
    for key, idx in tgt.groupby(["cohort", "endpoint"]).groups.items():
        ok = tgt.loc[idx, "p"].notna(); tgt.loc[tgt.loc[idx].index[ok], "q"] = bh_fdr(tgt.loc[idx, "p"].dropna().to_numpy())
    tgt.to_csv(out_dir / "mir200_per_target.tsv", sep="\t", index=False)

    m2 = res[res.family == "miR-200"]
    summary = {
        "module": "mirna_hallmark.outcome_mir200_deepdive",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pe_coherence_FDR": res[(res.metric == "pe_coherence") & (res.q < 0.1)].sort_values("p").apply(lambda r: f"{r['cohort']}/{r['miRNA']}/{r['subset']}/{r['endpoint']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "miR200_best_readout_subset": m2[m2.q < 0.1].sort_values("p").head(12).apply(lambda r: f"{r['cohort']}/{r['miRNA']}/{r['metric']}/{r['subset']}/{r['endpoint']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "subset_FDR_counts": {s: int(((res.subset == s) & (res.q < 0.1)).sum()) for s in ["all", "oncogene", "TSG", "TF", "highdep"]},
        "metric_FDR_counts": {m: int(((res.metric == m) & (res.q < 0.1)).sum()) for m in res.metric.unique()},
        "miR200_top_targets": tgt[tgt.q < 0.2].sort_values("p").head(12).apply(lambda r: f"{r['cohort']}/{r['target']}/{r['endpoint']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "caveats": ["miR-200 family deep-dive; readouts incl. pe_coherence; target-subsets by role/dependency",
                    "TCGA DFI+OS + Buffa DRFS; FDR within (cohort,endpoint)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[m200] pe_coherence FDR: {summary['pe_coherence_FDR']}")
    print(f"[m200] subset FDR counts: {summary['subset_FDR_counts']}")
    print(f"[m200] metric FDR counts: {summary['metric_FDR_counts']}")
    print(f"[m200] miR-200 best readout×subset: {summary['miR200_best_readout_subset'][:8]}")
    print(f"[m200] miR-200 top driving targets: {summary['miR200_top_targets']}")
    print(f"[m200] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
