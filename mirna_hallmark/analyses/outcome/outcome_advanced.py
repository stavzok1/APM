"""Roadmap Phase 3 outcome — regulatory-STATE survival (the right questions), both cohorts.

Reframes outcome away from pressure *magnitude* (exhaustively null, MH-60..64) toward regulatory
*state / configuration*, run identically in TCGA (DFI+OS) and Buffa/METABRIC-miRNA (DRFS):

  #1 COUPLING COHERENCE / DECOUPLING (per patient) — the framework-native "is control intact?":
       coherence_p = -corr over BY-neg edges of (miRNA-z, target-z) for patient p;
       decoupling   = per-edge cohort regression target~miRNA -> per-patient mean signed residual
                      (net de-repression / leaky control) + mean |residual| (incoherence).
  #2 PRESSURE CONFIGURATION — NMF latent pressure-states (combinations the mean misses) -> Cox.
  #3 NONLINEAR ML — Random Survival Forest + gradient-boosted Cox, 5-fold CV concordance,
       pressure(program+role) vs clinical-only (do interactions/combinations matter at all?).

Discipline: both cohorts; FDR within; pressure must beat clinical out-of-sample; TCGA hits need Buffa.

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_advanced``
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome

OUT_DIR = C.OUTPUT_ROOT / "outcome_advanced"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
TF_TSV = C.REPO_ROOT / "annotations" / "humantfs_lambert2018_tf.tsv"


def _zrows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def _neg_edges() -> pd.DataFrame:
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()
    e["arm"] = e["key"].str.split(r"\|\|").str[0]; e["gene"] = e["key"].str.split(r"\|\|").str[1]
    return e[["arm", "gene"]]


def _role_features(gene_p: pd.DataFrame) -> pd.DataFrame:
    from mirna_hallmark.gene_roles import load_gene_roles
    gr = load_gene_roles(); tf = set(pd.read_csv(TF_TSV, sep="\t").iloc[:, 0].astype(str))
    sets = {"TF": tf, "oncogene": set(gr[gr.malignancy_sign == 1].gene), "TSG": set(gr[gr.malignancy_sign == -1].gene)}
    rows = {"role_all": gene_p.mean(axis=0)}
    for r, gs in sets.items():
        g = gene_p.index.intersection(gs)
        if len(g) >= 5:
            rows[f"role_{r}"] = gene_p.loc[g].mean(axis=0)
    return pd.DataFrame(rows).T


# ---------------- cohort loaders ----------------
def _tcga_cohort() -> dict:
    mi = D.load_mirna_arms(); mi = mi[~mi.index.duplicated(keep="first")]
    rna = D.load_rna(); rna = rna[~rna.index.duplicated(keep="first")]
    cl = load_tcga_outcome()
    prog = pd.read_csv(C.OUTPUT_ROOT / "hallmark_interaction" / "hallmark_pressure_per_sample.tsv.gz", sep="\t", index_col=0)
    gp = pd.read_csv(C.OUTPUT_ROOT / "pressure_layer_comparison" / "L1_coupling_spine" / "gene_pressure_per_sample.tsv.gz", sep="\t", index_col=0)
    e = _neg_edges(); e = e[e["arm"].isin(mi.index) & e["gene"].isin(rna.index)]
    return {"name": "TCGA", "mi": mi, "rna": rna, "edges": e, "prog": prog, "gp": gp,
            "endpoints": {"DFI": (cl["dfi_t"], cl["dfi_e"]), "OS": (cl["os_t"], cl["os_e"])},
            "cov": pd.concat([cl["base"], cl["comp"]], axis=1)}


def _buffa_cohort() -> dict:
    from mirna_hallmark.eval.buffa_validation import load_buffa, BUFFA_MIRNA_SERIES
    from mirna_hallmark.analyses.outcome.outcome_buffa_pressure import _buffa_pressure, _clinical
    mi, rna = load_buffa()
    gpg, prog, _ = _buffa_pressure()
    e = _neg_edges()
    base = {re.sub(r"-(3p|5p)$", "", a): a for a in mi.index}                  # base-name match
    e["arm_b"] = e["arm"].map(lambda a: base.get(re.sub(r"-(3p|5p)$", "", a)))
    e = e.dropna(subset=["arm_b"]); e = e[e["gene"].isin(rna.index)]
    e = e.rename(columns={"arm": "arm_tcga", "arm_b": "arm"})[["arm", "gene"]]
    clin = _clinical(prog.columns)
    return {"name": "Buffa", "mi": mi, "rna": rna, "edges": e, "prog": prog, "gp": gpg,
            "endpoints": {"DRFS": (clin["T"], clin["E"])},
            "cov": pd.concat([_z(clin["er"]).rename("er")], axis=1)}


# ---------------- analyses ----------------
def coupling_coherence(coh: dict) -> pd.DataFrame:
    e = coh["edges"]
    samples = coh["mi"].columns.intersection(coh["rna"].columns)             # miRNA & RNA sample sets differ
    miz = _zrows(coh["mi"][samples].loc[e["arm"].values]).to_numpy()         # E x N (aligned)
    tgz = _zrows(coh["rna"][samples].loc[e["gene"].values]).to_numpy()
    ok = np.isfinite(miz).all(1) & np.isfinite(tgz).all(1)
    A, B = miz[ok], tgz[ok]
    Am = A - A.mean(0); Bm = B - B.mean(0)
    coher = pd.Series(-(Am * Bm).sum(0) / np.sqrt((Am ** 2).sum(0) * (Bm ** 2).sum(0)), index=samples)
    Ar = A - A.mean(1, keepdims=True); Br = B - B.mean(1, keepdims=True)
    slope = (Ar * Br).sum(1) / np.where((Ar ** 2).sum(1) == 0, np.nan, (Ar ** 2).sum(1))
    resid = B - (A * slope[:, None] + (B.mean(1) - slope * A.mean(1))[:, None])
    feat = pd.DataFrame({"coherence": coher, "decoupling_signed": pd.Series(np.nanmean(resid, 0), index=samples),
                         "decoupling_abs": pd.Series(np.nanmean(np.abs(resid), 0), index=samples)}).T
    return _cox_feats(feat, coh, "coupling_state")


def _pp_anticorr(Pz: pd.DataFrame, Ez: pd.DataFrame) -> pd.Series:
    """Per-patient -corr across genes between pressure-z and expression-z (functional repression coherence)."""
    A, B = Pz.to_numpy(), Ez.to_numpy()
    Am, Bm = A - np.nanmean(A, 0), B - np.nanmean(B, 0)
    den = np.sqrt(np.nansum(Am ** 2, 0) * np.nansum(Bm ** 2, 0))
    return pd.Series(-np.nansum(Am * Bm, 0) / np.where(den == 0, np.nan, den), index=Pz.columns)


def pressure_expression_coherence(coh: dict) -> pd.DataFrame:
    """#1b: per-patient anti-corr of gene PRESSURE vs its EXPRESSION — global, per-program, per-role."""
    from mirna_hallmark.hallmark_sets import HallmarkSets
    from mirna_hallmark.gene_roles import load_gene_roles
    samples = coh["gp"].columns.intersection(coh["rna"].columns)
    genes = coh["gp"].index.intersection(coh["rna"].index)
    Pz = _zrows(coh["gp"].loc[genes, samples]); Ez = _zrows(coh["rna"].loc[genes, samples])
    feats = {"pe_coh_global": _pp_anticorr(Pz, Ez)}
    for prog, members in HallmarkSets.load().sets.items():
        g = genes.intersection(members)
        if len(g) >= 10:
            feats[f"pe_coh_{prog.replace('HALLMARK_','')}"] = _pp_anticorr(Pz.loc[g], Ez.loc[g])
    gr = load_gene_roles(); tf = set(pd.read_csv(TF_TSV, sep="\t").iloc[:, 0].astype(str))
    for r, gs in {"TF": tf, "oncogene": set(gr[gr.malignancy_sign == 1].gene), "TSG": set(gr[gr.malignancy_sign == -1].gene)}.items():
        g = genes.intersection(gs)
        if len(g) >= 10:
            feats[f"pe_coh_role_{r}"] = _pp_anticorr(Pz.loc[g], Ez.loc[g])
    return _cox_feats(pd.DataFrame(feats).T, coh, "pe_coherence")


def nmf_factors(coh: dict, k: int = 12) -> pd.DataFrame:
    X = coh["prog"].T.dropna()                                                # sample x program
    Xn = (X - X.min(axis=0)).clip(lower=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        W = NMF(n_components=min(k, Xn.shape[1] - 1), init="nndsvda", max_iter=600, random_state=0).fit_transform(Xn)
    feat = pd.DataFrame(W, index=X.index, columns=[f"NMF{i+1}" for i in range(W.shape[1])]).T
    return _cox_feats(feat, coh, "nmf_factor")


def _cox_feats(feat: pd.DataFrame, coh: dict, fam: str) -> pd.DataFrame:
    rows = []
    for ep, (t, e) in coh["endpoints"].items():
        for f in feat.index:
            c, p = _cox_one(t, e, feat.loc[f].astype(float), coh["cov"])
            rows.append({"cohort": coh["name"], "analysis": fam, "feature": f, "endpoint": ep,
                         "hr": np.exp(c), "p": p})
    r = pd.DataFrame(rows)
    for ep, idx in r.groupby("endpoint").groups.items():
        okm = r.loc[idx, "p"].notna()
        r.loc[r.loc[idx].index[okm], "q"] = bh_fdr(r.loc[idx, "p"].dropna().to_numpy())
    return r


def ml_survival(coh: dict) -> list:
    feat = pd.concat([coh["prog"], _role_features(coh["gp"])])                # program + role features
    X = feat.T
    out = []
    for ep, (t, e) in coh["endpoints"].items():
        df = X.join(coh["cov"]).join(t.rename("T")).join(e.rename("E")).dropna()
        df = df[df["T"] > 0]
        if df["E"].sum() < 20:
            continue
        feats = [c for c in X.columns if c in df.columns]; clin = [c for c in coh["cov"].columns if c in df.columns]
        rng = np.random.default_rng(0); idx = rng.permutation(df.index.to_numpy()); folds = np.array_split(idx, 5)
        cfull, cclin = [], []
        for i in range(5):
            te = df.loc[folds[i]]; tr = df.drop(index=folds[i])
            if tr["E"].sum() < 10 or te["E"].sum() < 3:
                continue
            ytr = Surv.from_arrays(tr["E"].astype(bool), tr["T"]);
            for cols, store in ((feats + clin, cfull), (clin, cclin)):
                try:
                    m = RandomSurvivalForest(n_estimators=100, min_samples_leaf=15, max_features="sqrt", random_state=0, n_jobs=4).fit(tr[cols], ytr)
                    pred = m.predict(te[cols])
                    store.append(concordance_index_censored(te["E"].astype(bool), te["T"], pred)[0])
                except Exception:
                    pass
        out.append({"cohort": coh["name"], "endpoint": ep, "n_events": int(df["E"].sum()),
                    "rsf_concordance_pressure+clin": round(float(np.mean(cfull)), 3) if cfull else None,
                    "rsf_concordance_clin_only": round(float(np.mean(cclin)), 3) if cclin else None})
    return out


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cox_all, ml_all = [], []
    for loader in (_tcga_cohort, _buffa_cohort):
        coh = loader()
        print(f"[adv] {coh['name']}: {coh['edges'].shape[0]} edges, {coh['prog'].shape[0]} programs, endpoints {list(coh['endpoints'])}")
        cox_all.append(coupling_coherence(coh)); print(f"[adv]   {coh['name']} coherence done")
        cox_all.append(pressure_expression_coherence(coh)); print(f"[adv]   {coh['name']} pe-coherence done")
        cox_all.append(nmf_factors(coh)); print(f"[adv]   {coh['name']} NMF done")
        ml_all += ml_survival(coh); print(f"[adv]   {coh['name']} RSF done")
    cox = pd.concat(cox_all, ignore_index=True)
    cox.to_csv(out_dir / "advanced_cox.tsv", sep="\t", index=False)
    pd.DataFrame(ml_all).to_csv(out_dir / "advanced_ml.tsv", sep="\t", index=False)
    summary = {
        "module": "mirna_hallmark.outcome_advanced",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "coupling_state": cox[cox.analysis == "coupling_state"].apply(lambda r: f"{r['cohort']}/{r['endpoint']}/{r['feature']} HR={r['hr']:.2f} p={r['p']:.3f} q={r.get('q', np.nan):.3f}", axis=1).tolist(),
        "pe_coherence_FDR_hits": cox[(cox.analysis == "pe_coherence") & (cox.q < 0.1)].sort_values("q").apply(lambda r: f"{r['cohort']}/{r['endpoint']}/{r['feature']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "nmf_FDR_hits": cox[(cox.analysis == "nmf_factor") & (cox.q < 0.1)].apply(lambda r: f"{r['cohort']}/{r['endpoint']}/{r['feature']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "ml_concordance": ml_all,
        "caveats": ["coherence/decoupling are per-patient regulatory-state metrics (new); RSF 5-fold CV",
                    "pressure must beat clinical out-of-sample; TCGA findings need Buffa concordance"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[adv] coupling-state: {summary['coupling_state']}")
    print(f"[adv] pe-coherence FDR hits: {summary['pe_coherence_FDR_hits']}")
    print(f"[adv] NMF FDR hits: {summary['nmf_FDR_hits']}")
    print(f"[adv] ML concordance (pressure+clin vs clin): {ml_all}")
    print(f"[adv] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
