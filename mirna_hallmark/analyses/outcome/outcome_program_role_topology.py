"""Roadmap Phase 3 outcome — pressure per PROGRAM × gene-ROLE, with/without a TOPOLOGY prior, vs survival.

The finest resolution we hadn't tried: instead of whole-program pressure or *global* role aggregates,
test each program's genes **split by role** (all / TF / oncogene / TSG), aggregated two ways —
**unweighted** vs **topology-weighted** (architecture `w_arch` = reverse-PageRank centrality × redundancy,
from the signed OmniPath network in `geneset_architecture`). Asks: is a program's prognostic pressure
carried by its TFs / oncogenes / TSGs specifically, and does up-weighting topologically-central genes help?
Cox vs TCGA DFI+OS, composition-adjusted, BH-FDR. (Honest: any hit needs Buffa replication.)

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_program_role_topology``
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
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, load_tcga_outcome

OUT_DIR = C.OUTPUT_ROOT / "outcome_program_role_topology"
GP = C.OUTPUT_ROOT / "pressure_layer_comparison" / "L1_coupling_spine" / "gene_pressure_per_sample.tsv.gz"
ARCH = C.OUTPUT_ROOT / "tissue_reference" / "geneset_architecture" / "architecture_all_gene_weights.tsv"


def _features(gp: pd.DataFrame, arch: pd.DataFrame) -> pd.DataFrame:
    arch = arch.drop_duplicates(["hallmark", "gene"])
    feats = {}
    for prog, g in arch.groupby("hallmark"):
        g = g[g["gene"].isin(gp.index)]
        roles = {"all": g, "TF": g[g["is_tf"] == True], "onco": g[g["malignancy_sign"] == 1],
                 "TSG": g[g["malignancy_sign"] == -1]}
        for rname, gr in roles.items():
            if len(gr) < 5:
                continue
            P = gp.loc[gr["gene"].values]                                   # genes x sample
            feats[f"{prog}|{rname}|unw"] = P.mean(axis=0)
            w = pd.to_numeric(gr.set_index("gene")["w_arch"], errors="coerce").reindex(gr["gene"].values).fillna(0).to_numpy()
            if w.sum() > 0:
                feats[f"{prog}|{rname}|topo"] = pd.Series((P.to_numpy().T @ w) / w.sum(), index=P.columns)
    return pd.DataFrame(feats).T


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cl = load_tcga_outcome()
    gp = pd.read_csv(GP, sep="\t", index_col=0)
    arch = pd.read_csv(ARCH, sep="\t")
    feat = _features(gp, arch)
    print(f"[prt] {feat.shape[0]} program×role×weighting features")
    adj = pd.concat([cl["base"], cl["comp"]], axis=1)
    rows = []
    for f in feat.index:
        x = feat.loc[f].astype(float)
        if x.std(skipna=True) == 0 or x.notna().sum() < 50:
            continue
        prog, role, wt = f.split("|")
        dca, dpa = _cox_one(cl["dfi_t"], cl["dfi_e"], x, adj)
        oca, opa = _cox_one(cl["os_t"], cl["os_e"], x, adj)
        rows.append({"program": prog, "role": role, "weighting": wt, "feature": f,
                     "dfi_hr": np.exp(dca), "dfi_p": dpa, "os_hr": np.exp(oca), "os_p": opa})
    res = pd.DataFrame(rows)
    for col in ("dfi_p", "os_p"):
        ok = res[col].notna(); res.loc[ok, col.replace("_p", "_q")] = bh_fdr(res.loc[ok, col].to_numpy())
    res.to_csv(out_dir / "program_role_topology_survival.tsv", sep="\t", index=False)

    fdr = res[(res.dfi_q < 0.1) | (res.os_q < 0.1)]
    # did topology weighting help? compare min q within (program,role) topo vs unw
    piv = res.assign(minq=res[["dfi_q", "os_q"]].min(axis=1)).pivot_table(index=["program", "role"], columns="weighting", values="minq")
    topo_better = int((piv.get("topo", pd.Series(dtype=float)) < piv.get("unw", pd.Series(dtype=float))).sum()) if {"topo", "unw"} <= set(piv.columns) else 0
    summary = {
        "module": "mirna_hallmark.outcome_program_role_topology",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_features": int(feat.shape[0]),
        "n_FDR_hits": int(len(fdr)),
        "FDR_hits": fdr.sort_values([c for c in ("os_q", "dfi_q")]).head(20)
            .apply(lambda r: f"{r['program']}|{r['role']}|{r['weighting']} DFI HR={r['dfi_hr']:.2f} q={r['dfi_q']:.3f} | OS HR={r['os_hr']:.2f} q={r['os_q']:.3f}", axis=1).tolist(),
        "topology_helped_n_program_roles": topo_better,
        "best_overall": res.reindex(res[["dfi_q", "os_q"]].min(axis=1).sort_values().index).head(8)
            .apply(lambda r: f"{r['program']}|{r['role']}|{r['weighting']} minq={min(r['dfi_q'], r['os_q']):.3f}", axis=1).tolist(),
        "caveats": ["TCGA-only, composition-adjusted; FDR over all features (large surface)",
                    "topology prior = architecture w_arch (revPageRank×redundancy); any hit needs Buffa replication"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[prt] FDR hits: {summary['n_FDR_hits']} | topology beat unweighted in {topo_better} program-roles")
    print(f"[prt] top: {summary['best_overall']}")
    print(f"[prt] FDR list: {summary['FDR_hits'][:10]}")
    print(f"[prt] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
