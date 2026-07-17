"""Design build 1e — can we improve the functional prognostic? Three levers:

  (A) LOCKED-PANEL generalization — train the realized LASSO-Cox on METABRIC(Buffa), FREEZE it, and SCORE
      TCGA (no fitting in sparse TCGA). The honest generalization test — does the METABRIC signature transfer?
  (B) COMBINE axes — realized + decoupling together (level + dysfunction) vs realized alone, METABRIC CV.
  (C) NONLINEAR (Random Survival Forest) on the realized features, METABRIC CV — interactions the LASSO misses.

Run: ``.venv/bin/python3 -m mirna_hallmark.pressure_prognostic_improve``
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.gene_roles import load_gene_dependency
from mirna_hallmark.analyses.outcome.outcome_survival import _z, load_tcga_outcome
from mirna_hallmark.analyses.pressure_dev.pressure_prognostic_signature import _cv_concordance, _log2tpm_tcga, _edges_w
from mirna_hallmark.analyses.pressure_dev.pressure_prognostic_gene_centric import _buffa_full_clin
from mirna_hallmark.pressure_signature_curated_info import _realized

OUT_DIR = C.OUTPUT_ROOT / "pressure_prognostic_improve"


def _decoupling(R, mi, arm_of):
    out = {}
    for arm in R.index:
        a = arm_of(arm)
        if a is None or a not in mi.index:
            continue
        ab = _z(mi.loc[a]); df = pd.concat([ab.rename("a"), R.loc[arm].rename("r")], axis=1).dropna()
        if len(df) < 30:
            continue
        b = np.polyfit(df["a"], df["r"], 1)
        out[f"dec_{arm}"] = -(R.loc[arm] - (b[0] * ab + b[1]))
    return pd.DataFrame(out).T


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    edges = _edges_w(); dep = load_gene_dependency().set_index("gene")["dep_role_weight"]
    mi = D.load_mirna_arms(); mi = mi[~mi.index.duplicated(keep="first")]
    l2 = _log2tpm_tcga(); cl = load_tcga_outcome(); cov = pd.concat([cl["base"], cl["comp"]], axis=1)
    tR = _realized(lambda a: a if a in mi.index else None, mi, l2, edges, dep)
    from mirna_hallmark.eval.buffa_validation import load_buffa
    bmi, brna = load_buffa(); base = {re.sub(r"-(3p|5p)$", "", x): x for x in bmi.index}
    barm_of = lambda a: base.get(re.sub(r"-(3p|5p)$", "", a))
    bR = _realized(barm_of, bmi, brna, edges, dep)
    bT, bE, bfull = _buffa_full_clin(bmi)
    results = {}

    # ---- (A) locked panel: fit on Buffa, score TCGA ----
    shared = bR.index.intersection(tR.index)
    Xb = bR.loc[shared].T; dfb = Xb.join(bT.rename("T")).join(bE.rename("E")).dropna(); dfb = dfb[dfb["T"] > 0]
    dfb = dfb.loc[:, dfb.std(numeric_only=True) > 0]
    feats = [c for c in shared if c in dfb.columns]
    cph = CoxPHFitter(penalizer=1.0, l1_ratio=0.0).fit(dfb[feats + ["T", "E"]], "T", "E")   # RIDGE (dense — signal is distributed)
    panel = feats; coef = cph.params_[panel]
    print(f"[improve] locked Buffa RIDGE panel: {len(panel)} miRNAs (dense)")
    Xt = tR.loc[panel].T                                                  # TCGA realized for the locked panel
    lp = Xt.mul(coef, axis=1).sum(axis=1)                                 # frozen linear predictor on TCGA
    for ep in ("dfi", "os"):
        t, e = cl[f"{ep}_t"], cl[f"{ep}_e"]
        df = pd.concat([lp.rename("lp"), t.rename("T"), e.rename("E"), cov], axis=1).dropna(); df = df[df["T"] > 0]
        if df["E"].sum() < 20:
            continue
        c_lp = concordance_index(df["T"], -df["lp"], df["E"])
        ccols = [c for c in cov.columns if c in df.columns and df[c].std() > 0]
        try:
            mc = CoxPHFitter(penalizer=0.1).fit(df[ccols + ["T", "E"]], "T", "E")
            c_clin = concordance_index(df["T"], -mc.predict_partial_hazard(df[ccols]), df["E"])
            mcl = CoxPHFitter(penalizer=0.1).fit(df[ccols + ["lp", "T", "E"]], "T", "E")
            c_both = concordance_index(df["T"], -mcl.predict_partial_hazard(df[ccols + ["lp"]]), df["E"])
        except Exception:
            c_clin = c_both = None
        results[f"LOCKED_Buffa->TCGA/{ep.upper()}"] = {"panel_only_C": round(float(c_lp), 3),
            "clinical_C": round(float(c_clin), 3) if c_clin else None,
            "clinical+panel_C": round(float(c_both), 3) if c_both else None,
            "delta_vs_clinical": round(float(c_both - c_clin), 3) if c_both else None}

    # ---- (B) combine realized + decoupling (Buffa CV) ----
    bDec = _decoupling(bR, bmi, barm_of)
    combined = pd.concat([bR, bDec])
    results["Buffa/realized_only/DRFS"] = _cv_concordance(bR, bfull, bT, bE)
    results["Buffa/realized+decoupling/DRFS"] = _cv_concordance(combined, bfull, bT, bE)

    # ---- (C) RSF nonlinear on Buffa realized ----
    X = bR.T; df = X.join(bfull).join(bT.rename("T")).join(bE.rename("E")).dropna(); df = df[df["T"] > 0]
    df = df.loc[:, df.std(numeric_only=True) > 0]
    fcols = [c for c in X.columns if c in df.columns]; ccols = [c for c in bfull.columns if c in df.columns]
    rng = np.random.default_rng(0); idx = rng.permutation(df.index.to_numpy()); folds = np.array_split(idx, 5)
    rf, rc = [], []
    for i in range(5):
        te = df.loc[folds[i]]; tr = df.drop(index=folds[i])
        if tr["E"].sum() < 10 or te["E"].sum() < 3:
            continue
        ytr = Surv.from_arrays(tr["E"].astype(bool), tr["T"])
        for cols, store in ((fcols + ccols, rf), (ccols, rc)):
            try:
                m = RandomSurvivalForest(n_estimators=200, min_samples_leaf=20, max_features="sqrt", random_state=0, n_jobs=4).fit(tr[cols], ytr)
                store.append(concordance_index_censored(te["E"].astype(bool), te["T"], m.predict(te[cols]))[0])
            except Exception:
                pass
    results["Buffa/RSF_realized/DRFS"] = {"rsf_realized+clin": round(float(np.mean(rf)), 3) if rf else None,
                                          "rsf_clin_only": round(float(np.mean(rc)), 3) if rc else None,
                                          "delta": round(float(np.mean(rf) - np.mean(rc)), 3) if rf and rc else None}

    summary = {"module": "mirna_hallmark.pressure_prognostic_improve",
               "generated_utc": datetime.now(timezone.utc).isoformat(), "results": results,
               "caveats": ["locked panel = trained on Buffa, scored on TCGA (no TCGA fitting); the honest generalization test"]}
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    for k, v in results.items():
        print(f"[improve] {k}: {v}")
    print(f"[improve] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
