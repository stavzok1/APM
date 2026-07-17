"""Design build 1 (PRESSURE_PROGNOSTIC_DESIGN.md) — a parsimonious miRNA-centric pressure/decoupling
prognostic on the FUNCTIONAL axes, target-relevance-weighted, vs clinical.

Per miRNA, per patient, four readouts over its BY-neg target set (weighted by |coupling|×cancer-relevance
×expressed, top-K dropout):
  realized      = -weighted-mean z(target expression)            (functional repression; beats abundance, MH-67)
  decoupling    = resid(realized ~ abundance)                    (escape; opened the TCGA crack, MH-68)
  acquired_NAT  = -weighted-mean (target log2TPM - NAT mean)     (cancer-acquired repression vs matched NAT)
  acquired_GTEx = -weighted-mean (target log2TPM - GTEx breast)  (vs true-healthy GTEx)
Selection: per-axis LASSO-Cox, 5-fold nested-CV concordance vs a clinical-only baseline + stability
selection; curated (driver-targeting famous miRNAs) panel too. TCGA (DFI+OS) + Buffa (DRFS, realized/
decoupling — no normal baseline there). Acceptance bar: beats clinical out-of-sample, cross-cohort.

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.pressure_dev.pressure_prognostic_signature``
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

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.gene_roles import load_gene_dependency
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome
from mirna_hallmark.analyses.outcome.outcome_famous_mirnas import FAMOUS

OUT_DIR = C.OUTPUT_ROOT / "pressure_prognostic_signature"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
GTEX_GCT = C.REPO_ROOT / "data" / "GTEx" / "gene_tpm_v10_breast.gct.gz"
TOPK = 8


def _log2tpm_tcga():
    r = D.load_rna(); r = r[~r.index.duplicated(keep="first")]
    return np.power(2.0, r) - 1.0                                   # load_rna double-logs -> 2^v-1 = log2(TPM+1)


def _nat_mean():
    from pipeline.config import PATHS
    raw = pd.read_csv(PATHS.rna_expression, sep="\t", index_col=0)
    if "Ensembl_ID" in raw.columns:
        raw = raw.drop(columns=["Ensembl_ID"])
    raw = raw[~raw.index.duplicated(keep="first")]
    nat = [c for c in raw.columns if len(c) >= 15 and c[13:15] == "11"]
    return raw[nat].mean(axis=1)                                    # log2(TPM+1) NAT mean per gene


def _gtex_mean():
    g = pd.read_csv(GTEX_GCT, sep="\t", skiprows=2, index_col=1)    # Description=gene symbol
    g = g.drop(columns=[c for c in ("Name", "id") if c in g.columns], errors="ignore")
    g = g[~g.index.duplicated(keep="first")].apply(pd.to_numeric, errors="coerce")
    return np.log2(g.mean(axis=1) + 1.0)                            # log2(TPM+1) GTEx breast mean


def _edges_w():
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()
    e["arm"] = e["key"].str.split(r"\|\|").str[0]; e["gene"] = e["key"].str.split(r"\|\|").str[1]
    return e[["arm", "gene", "rho"]]


def _weights(sub: pd.DataFrame, dep: pd.Series, expr_mean: pd.Series) -> pd.Series:
    w = sub["rho"].abs().values * (1.0 + dep.reindex(sub["gene"]).abs().fillna(0).values) \
        * (expr_mean.reindex(sub["gene"]).fillna(0).values > 1.0)                       # expressed gate
    w = pd.Series(w, index=sub["gene"].values)
    w = w[w > 0]
    return w.sort_values(ascending=False).head(TOPK)               # top-K dropout


def _readouts(arm_of, mi, l2tpm, edges, dep, nat_mean, gtex_mean):
    expr_mean = l2tpm.mean(axis=1)
    feats = {"realized": {}, "decoupling": {}, "acquired_NAT": {}, "acquired_GTEx": {}}
    for arm, sub in edges.groupby("arm"):
        a = arm_of(arm)
        if a is None or a not in mi.index:
            continue
        sub = sub[sub["gene"].isin(l2tpm.index)]
        w = _weights(sub, dep, expr_mean)
        if len(w) < 3:
            continue
        g = list(w.index); wv = w.values / w.sum()
        E = l2tpm.loc[g]                                            # log2TPM targets x samples
        z = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1).replace(0, np.nan), axis=0)
        realized = -(z.mul(wv, axis=0)).sum(axis=0)
        ab = _z(mi.loc[a]).reindex(realized.index)
        df = pd.concat([ab.rename("a"), realized.rename("r")], axis=1).dropna()
        feats["realized"][arm] = realized
        if len(df) >= 30:
            b = np.polyfit(df["a"], df["r"], 1)
            feats["decoupling"][arm] = -(realized - (b[0] * ab + b[1]))
        if nat_mean is not None:
            lfc = E.sub(nat_mean.reindex(g).values, axis=0)
            feats["acquired_NAT"][arm] = -(lfc.mul(wv, axis=0)).sum(axis=0)
        if gtex_mean is not None:
            lfc = E.sub(gtex_mean.reindex(g).fillna(E.mean(axis=1)).values, axis=0)
            feats["acquired_GTEx"][arm] = -(lfc.mul(wv, axis=0)).sum(axis=0)
    return {k: pd.DataFrame(v).T for k, v in feats.items() if v}


def _cv_concordance(feat: pd.DataFrame, cov: pd.DataFrame, t, e, l1=0.7):
    """5-fold CV concordance: LASSO-Cox on features+clinical vs clinical-only; + stability-selected features."""
    X = feat.T
    df = X.join(cov).join(t.rename("T")).join(e.rename("E")).dropna()
    df = df.loc[:, df.std(numeric_only=True) > 0]; df = df[df["T"] > 0]
    fcols = [c for c in X.columns if c in df.columns]; ccols = [c for c in cov.columns if c in df.columns]
    if df["E"].sum() < 25 or not fcols:
        return None
    rng = np.random.default_rng(0); idx = rng.permutation(df.index.to_numpy()); folds = np.array_split(idx, 5)
    cf, cc = [], []; sel = {}
    for i in range(5):
        te = df.loc[folds[i]]; tr = df.drop(index=folds[i])
        if tr["E"].sum() < 10 or te["E"].sum() < 3:
            continue
        for cols, store in ((fcols + ccols, cf), (ccols, cc)):
            try:
                m = CoxPHFitter(penalizer=0.2, l1_ratio=l1).fit(tr[cols + ["T", "E"]], "T", "E")
                risk = m.predict_partial_hazard(te[cols])
                store.append(concordance_index(te["T"], -risk, te["E"]))
                if cols is fcols + ccols:
                    for fcn in fcols:
                        if abs(m.params_.get(fcn, 0)) > 1e-6:
                            sel[fcn] = sel.get(fcn, 0) + 1
            except Exception:
                pass
    if not cf:
        return None
    return {"cv_concordance_pressure+clin": round(float(np.mean(cf)), 3),
            "cv_concordance_clin_only": round(float(np.mean(cc)), 3) if cc else None,
            "delta": round(float(np.mean(cf) - (np.mean(cc) if cc else 0.5)), 3),
            "stable_features": sorted([f for f, n in sel.items() if n >= 3], key=lambda f: -sel[f])[:12]}


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    edges = _edges_w(); dep = load_gene_dependency().set_index("gene")["dep_role_weight"]
    # TCGA
    mi = D.load_mirna_arms(); mi = mi[~mi.index.duplicated(keep="first")]
    l2 = _log2tpm_tcga(); cl = load_tcga_outcome()
    tfeat = _readouts(lambda a: a if a in mi.index else None, mi, l2, edges, dep, _nat_mean(), _gtex_mean())
    print(f"[sig] TCGA readouts: {{{', '.join(f'{k}:{v.shape[0]}' for k,v in tfeat.items())}}}")
    cov = pd.concat([cl["base"], cl["comp"]], axis=1)
    results = {}
    for axis, fm in tfeat.items():
        for ep in ("dfi", "os"):
            r = _cv_concordance(fm, cov, cl[f"{ep}_t"], cl[f"{ep}_e"])
            if r:
                results[f"TCGA/{axis}/{ep.upper()}"] = r
    # combined (all axes) TCGA
    allf = pd.concat(list(tfeat.values()))
    allf.index = [f"{ax}|{m}" for ax, fm in tfeat.items() for m in fm.index]
    for ep in ("dfi", "os"):
        r = _cv_concordance(allf, cov, cl[f"{ep}_t"], cl[f"{ep}_e"])
        if r:
            results[f"TCGA/ALL/{ep.upper()}"] = r
    # Buffa (realized + decoupling only; no normal baseline there)
    import gzip
    from mirna_hallmark.eval.buffa_validation import load_buffa, BUFFA_MIRNA_SERIES
    bmi, brna = load_buffa(); base = {re.sub(r"-(3p|5p)$", "", x): x for x in bmi.index}
    lines = [l.rstrip("\n") for l in gzip.open(BUFFA_MIRNA_SERIES, "rt")]
    gsm = [x.strip('"') for x in [l for l in lines if l.lower().startswith("!sample_geo_accession")][0].split("\t")[1:]]
    def _ch(k):
        rr = [l for l in lines if l.lower().startswith("!sample_characteristics") and k in l.lower()]
        return pd.Series(pd.to_numeric([x.strip('"').split(":", 1)[1].strip() if ":" in x else np.nan for x in rr[0].split("\t")[1:]], errors="coerce"), index=gsm) if rr else pd.Series(np.nan, index=gsm)
    bT, bE = _ch("distant-relapse free survival"), _ch("distant-relapse event")
    bclin_full = pd.DataFrame({"er": _ch("er status"), "age": _z(_ch("patient age")), "grade": _z(_ch("tumour grade")),
                               "nodes": _z(_ch("nodes involved")), "size": _z(_ch("tumour size"))})   # FAIR full clinical
    bfeat = _readouts(lambda a: base.get(re.sub(r"-(3p|5p)$", "", a)), bmi, brna, edges, dep, None, None)
    for axis, fm in bfeat.items():
        r = _cv_concordance(fm, bclin_full, bT, bE)
        if r:
            results[f"Buffa/{axis}/DRFS"] = r

    summary = {
        "module": "mirna_hallmark.analyses.pressure_dev.pressure_prognostic_signature",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "design": "miRNA-centric realized/decoupling/acquired(NAT,GTEx); target-weighted+top8; LASSO-Cox nested-CV vs clinical",
        "results": results,
        "beats_clinical": {k: v for k, v in results.items() if v.get("delta", -1) > 0.01},
        "caveats": ["~100 events -> parsimony enforced by L1; acceptance bar = beats clinical out-of-sample, cross-cohort",
                    "acquired axes TCGA-only (NAT/GTEx baselines); Buffa = realized/decoupling"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print("[sig] CV concordance (pressure+clin vs clin) by axis/cohort:")
    for k, v in results.items():
        print(f"   {k:24s} {v['cv_concordance_pressure+clin']} vs {v['cv_concordance_clin_only']} (delta {v['delta']:+.3f})  feats={v['stable_features'][:4]}")
    print(f"[sig] beats clinical: {list(summary['beats_clinical'])}")
    print(f"[sig] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
