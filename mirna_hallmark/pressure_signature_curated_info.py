"""Design builds 1c/1d — the two remaining pieces of PRESSURE_PROGNOSTIC_DESIGN.md:

  (1c) CURATED driver-targeting panel — a small interpretable composite (mean realized-repression over
       famous TS-miRNAs / oncomiRs) vs full clinical. Does a biology-first hand-panel match the learned one?
  (1d) INFORMATION-CAPTURE — the unsupervised objective: PCA/NMF of the realized-pressure matrix → the
       MINIMAL set of miRNAs that reconstructs the overall pressure landscape, + how much variance a compact
       panel captures, + whether that compact info-panel retains the METABRIC prognostic gain.

TCGA (DFI+OS) + Buffa/METABRIC (DRFS, vs full clinical). Built on the realized readout (MH-73/74 winner).

Run: ``.venv/bin/python3 -m mirna_hallmark.pressure_signature_curated_info``
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.gene_roles import load_gene_dependency
from mirna_hallmark.analyses.outcome.outcome_survival import _z, load_tcga_outcome
from mirna_hallmark.analyses.outcome.outcome_famous_mirnas import FAMOUS
from mirna_hallmark.analyses.pressure_dev.pressure_prognostic_signature import (
    _cv_concordance, _log2tpm_tcga, _edges_w,
)
from mirna_hallmark.analyses.pressure_dev.pressure_prognostic_gene_centric import _buffa_full_clin

OUT_DIR = C.OUTPUT_ROOT / "pressure_signature_curated_info"


def _realized(arm_of, mi, l2tpm, edges, dep):
    """miRNA x patient realized-repression matrix (target-relevance-weighted, top-8)."""
    expr_mean = l2tpm.mean(axis=1); feats = {}
    for arm, sub in edges.groupby("arm"):
        a = arm_of(arm)
        if a is None or a not in mi.index:
            continue
        sub = sub[sub["gene"].isin(l2tpm.index)]
        w = sub["rho"].abs().values * (1 + dep.reindex(sub["gene"]).abs().fillna(0).values) * (expr_mean.reindex(sub["gene"]).fillna(0).values > 1)
        w = pd.Series(w, index=sub["gene"].values); w = w[w > 0].sort_values(ascending=False).head(8)
        if len(w) < 3:
            continue
        g = list(w.index); wv = w.values / w.sum()
        E = l2tpm.loc[g]; z = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1).replace(0, np.nan), axis=0)
        feats[arm] = -(z.mul(wv, axis=0)).sum(axis=0)
    return pd.DataFrame(feats).T


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    edges = _edges_w(); dep = load_gene_dependency().set_index("gene")["dep_role_weight"]
    mi = D.load_mirna_arms(); mi = mi[~mi.index.duplicated(keep="first")]
    l2 = _log2tpm_tcga(); cl = load_tcga_outcome()
    tR = _realized(lambda a: a if a in mi.index else None, mi, l2, edges, dep)
    cov = pd.concat([cl["base"], cl["comp"]], axis=1)
    from mirna_hallmark.eval.buffa_validation import load_buffa
    bmi, brna = load_buffa(); base = {re.sub(r"-(3p|5p)$", "", x): x for x in bmi.index}
    bR = _realized(lambda a: base.get(re.sub(r"-(3p|5p)$", "", a)), bmi, brna, edges, dep)
    bT, bE, bfull = _buffa_full_clin(bmi)
    results = {}

    # --- (1c) curated composites ---
    ts = [m for m in FAMOUS if FAMOUS[m].startswith("TS")]; onc = [m for m in FAMOUS if FAMOUS[m].startswith("onco")]
    def composites(R):
        return pd.DataFrame({"curated_TS_realized": R.loc[R.index.intersection(ts)].mean(axis=0),
                             "curated_onco_realized": R.loc[R.index.intersection(onc)].mean(axis=0)}).T
    for ep in ("dfi", "os"):
        r = _cv_concordance(composites(tR), cov, cl[f"{ep}_t"], cl[f"{ep}_e"])
        if r: results[f"TCGA/curated/{ep.upper()}"] = r
    r = _cv_concordance(composites(bR), bfull, bT, bE)
    if r: results["Buffa/curated/DRFS"] = r

    # --- (1d) information-capture: PCA of TCGA realized -> minimal reconstructing panel ---
    X = tR.dropna(axis=1, how="any").T                                  # patients x miRNA
    X = X.loc[:, X.std() > 0]
    Xc = (X - X.mean()) / X.std()
    pca = PCA(n_components=15, random_state=0).fit(Xc)
    cum = np.cumsum(pca.explained_variance_ratio_)
    load = np.abs(pca.components_).max(axis=0)                          # max |loading| per miRNA across 15 PCs
    panel = list(X.columns[np.argsort(-load)[:20]])                    # minimal 20-miRNA info-panel
    # variance captured by reconstructing from the 20-miRNA panel's own PCA
    Xp = Xc[panel]
    var_panel = float(np.cumsum(PCA(n_components=min(10, len(panel)), random_state=0).fit(Xp).explained_variance_ratio_)[-1])
    info = {"n_miRNA_total": int(X.shape[1]), "var_top10_PCs": round(float(cum[9]), 3),
            "var_top15_PCs": round(float(cum[-1]), 3), "info_panel_20": panel,
            "var_captured_within_panel": round(var_panel, 3)}
    # does the 20-miRNA info-panel RETAIN the METABRIC prognostic gain?
    bpanel = bR.loc[bR.index.intersection(panel)]
    rp = _cv_concordance(bpanel, bfull, bT, bE)
    if rp: results["Buffa/info_panel20/DRFS"] = rp
    tpanel = tR.loc[tR.index.intersection(panel)]
    for ep in ("dfi", "os"):
        r = _cv_concordance(tpanel, cov, cl[f"{ep}_t"], cl[f"{ep}_e"])
        if r: results[f"TCGA/info_panel20/{ep.upper()}"] = r

    summary = {
        "module": "mirna_hallmark.pressure_signature_curated_info",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "information_capture": info,
        "results": results,
        "beats_clinical": {k: v for k, v in results.items() if v.get("delta", -1) > 0.01},
        "caveats": ["curated = mean realized over famous TS / onco miRNAs; info-panel = top-20 by PCA loading",
                    "TCGA sparse (overfit); Buffa/METABRIC vs full clinical is the test; METABRIC-full = powered"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[curated-info] info-capture: {info['n_miRNA_total']} miRNAs; top-10 PCs capture {info['var_top10_PCs']:.0%}; "
          f"20-miRNA panel internal var {info['var_captured_within_panel']:.0%}")
    print(f"[curated-info] info-panel (20): {panel[:10]} ...")
    print("[curated-info] CV concordance vs clinical:")
    for k, v in results.items():
        print(f"   {k:28s} {v['cv_concordance_pressure+clin']} vs {v['cv_concordance_clin_only']} (delta {v['delta']:+.3f})")
    print(f"[curated-info] beats clinical: {list(summary['beats_clinical'])}")
    print(f"[curated-info] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
