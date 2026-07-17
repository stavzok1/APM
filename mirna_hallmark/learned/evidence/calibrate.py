"""Calibrate the evidence-class weights against an ORTHOGONAL functional readout (Design §Decision E;
BUILD_PLAN §1b — "does the fusion learn its hyperparameters non-circularly?").

`methods.CLASS_WEIGHT` are **hand-set** (reporter 3.0 > protein 2.5 > qpcr 1.5 > chimeric/degradome 1.0 >
ago_clip 0.5). This module asks whether an INDEPENDENT dataset — direct miRNA transfection + RNA-seq
(McGeary et al. 2019, GSE140217: 17 transfected miRNAs × HeLa, batch-normalized logTPM) — recovers that
*ordering*. If it does, the assay-directness prior is data-grounded, not circular hand-tuning.

Ground-truth functional effect (mRNA level), leave-one-out reference:
    repression(g,m) = mean_{m'≠m} logTPM(g,m')  −  logTPM(g,m)          (+ = repressed by m)
Joined per **seed family** (the transfected duplex's guide seed = the identified estimand, Design §F) ×
gene to the ledger's per-assay-class distinct-PMID counts. Then:

  (A) curated vs non-curated repression — does being in the ledger at all capture function?
  (B) learned class weights vs hand-set — NNLS of repression on per-class log1p(#PMID), **within-family
      centered** (removes the per-miRNA potency confound) and controlling total evidence (tests assay
      *directness*, not merely *amount*). Report Spearman(learned order, hand-set order).
  (C) does the fused `ledger_weight` track repression better than a raw PMID count (is the weighting worth it)?

Caveat, stated up front: the readout is **mRNA**, so protein-only classes (western/proteomics =
translational repression with little mRNA change) are *expected* to under-predict here. The learned
weights are therefore "weights for the mRNA task" — which is exactly the model's Y (TCGA mRNA). We report
them as such, not as a universal truth.

Bridges: UCSC refGene (NM_→symbol, data/external_cache/refseq), mirGeneDB-confirmed guide arms (17 short
names → guide arm → seed family via families.py). CLI: python -m mirna_hallmark.learned.evidence.calibrate
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned.evidence import ledger as L
from mirna_hallmark.learned.evidence.methods import CLASS_WEIGHT

TF = Path("data/external_cache/transfection/GSE140217_HeLa_transfection_logtpm_batchnorm.txt.gz")
REFGENE = Path("data/external_cache/refseq/refGene_hg38.txt.gz")
CACHE = Path("mirna_hallmark/output/learned/transfection_calibration.parquet")

# 17 transfected short names → GUIDE arm (miRBase spelling; guide 5p/3p confirmed against mirGeneDB
# hsa_guide_mat.fas). lsy-6 is C. elegans (no human ortholog) → dropped.
GUIDE_ARM: dict[str, str | None] = {
    "let7": "hsa-let-7a-5p", "lsy6": None, "mir1": "hsa-miR-1-3p", "mir124": "hsa-miR-124-3p",
    "mir137": "hsa-miR-137-3p", "mir139": "hsa-miR-139-5p", "mir143": "hsa-miR-143-3p",
    "mir144": "hsa-miR-144-3p", "mir153": "hsa-miR-153-3p", "mir155": "hsa-miR-155-5p",
    "mir182": "hsa-miR-182-5p", "mir199a": "hsa-miR-199a-3p", "mir204": "hsa-miR-204-5p",
    "mir205": "hsa-miR-205-5p", "mir216b": "hsa-miR-216b-5p", "mir223": "hsa-miR-223-3p",
    "mir7": "hsa-miR-7-5p",
}
# assay classes carried into the fused weight (drop 'other' = no weight). Order = hand-set descending.
CLASSES = ["reporter", "western", "proteomics", "qpcr_rna", "degradome", "chimeric", "ago_clip"]


def _nm2sym() -> dict:
    ref = pd.read_csv(REFGENE, sep="\t", header=None, usecols=[1, 12], names=["acc", "sym"], dtype=str)
    ref["acc"] = ref["acc"].str.split(".").str[0]
    return dict(zip(ref["acc"], ref["sym"]))


def _repression_long() -> pd.DataFrame:
    """Tidy (family, gene, repression, curated?) — the guide arm's seed family × measured gene."""
    nm2sym = _nm2sym()
    tf = pd.read_csv(TF, sep="\t", index_col=0)
    tf.index = pd.Series(tf.index).str.split(".").str[0].map(nm2sym).values
    tf = tf[tf.index.notna()]
    tf = tf[~tf.index.duplicated()]
    k = tf.shape[1]
    loo = (tf.sum(axis=1).to_numpy()[:, None] - tf.to_numpy()) / (k - 1)
    rep = pd.DataFrame(loo - tf.to_numpy(), index=tf.index, columns=tf.columns)   # + = repressed
    arms = {sn: a for sn, a in GUIDE_ARM.items() if a}
    fam = FAM.family_of(pd.Index(list(arms.values())))
    long = []
    for sn, arm in arms.items():
        if sn not in rep.columns:
            continue
        f = fam.get(arm)
        s = rep[sn].rename("repression").reset_index().rename(columns={"index": "gene"})
        s["family"] = f
        long.append(s)
    out = pd.concat(long, ignore_index=True)
    out.columns = ["gene", "repression", "family"]
    return out[["family", "gene", "repression"]]


def _class_counts() -> pd.DataFrame:
    """Per (family, gene) distinct-PMID count for each assay class, from the ledger, restricted to the
    16 transfected guide-seed families. Plus the fused `ledger_weight` and total distinct PMIDs."""
    arms = {sn: a for sn, a in GUIDE_ARM.items() if a}
    fam = FAM.family_of(pd.Index(list(arms.values())))
    led = L.build_ledger()
    fam_all = FAM.family_of(pd.Index(led["arm"].unique()))
    target_fams = set(fam.values)
    led = led.copy()
    led["family"] = led["arm"].map(fam_all)
    led = led[led["family"].isin(target_fams)]
    per = (led.groupby(["family", "gene", "assay_class"])["pmid"].nunique()
           .reset_index(name="n_pmid"))
    wide = per.pivot_table(index=["family", "gene"], columns="assay_class",
                           values="n_pmid", fill_value=0).reset_index()
    for c in CLASSES:
        if c not in wide.columns:
            wide[c] = 0
    logc = {c: np.log1p(wide[c].to_numpy()) for c in CLASSES}
    wide["ledger_weight"] = sum(CLASS_WEIGHT[c] * logc[c] for c in CLASSES)
    wide["total_pmid"] = wide[CLASSES].sum(axis=1)
    return wide


def _within_family_center(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Subtract the per-family mean from each column (= family fixed effects / potency de-confounding)."""
    out = df.copy()
    for c in cols:
        out[c] = df[c] - df.groupby("family")[c].transform("mean")
    return out


def calibrate() -> dict:
    rep = _repression_long()
    cc = _class_counts()
    j = rep.merge(cc, on=["family", "gene"], how="left")
    j["curated"] = j["ledger_weight"].notna()
    j[CLASSES + ["ledger_weight", "total_pmid"]] = j[CLASSES + ["ledger_weight", "total_pmid"]].fillna(0.0)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    j.to_parquet(CACHE)

    # (A) curated vs non-curated repression — does curation capture function? (within-family z, then pool)
    z = j.copy()
    z["rep_z"] = z["repression"] - z.groupby("family")["repression"].transform("mean")
    z["rep_z"] = z["rep_z"] / z.groupby("family")["repression"].transform("std").replace(0, np.nan)
    cur = z.loc[z["curated"], "rep_z"].dropna()
    non = z.loc[~z["curated"], "rep_z"].dropna()
    u, p_ab = mannwhitneyu(cur, non, alternative="greater")
    A = {"n_curated": int(len(cur)), "n_noncurated": int(len(non)),
         "mean_rep_z_curated": float(cur.mean()), "mean_rep_z_noncurated": float(non.mean()),
         "mannwhitney_p_greater": float(p_ab)}

    # (B) learned class weights — only the curated edges (where class counts exist); within-family center
    #     both target and features (FE estimator) so it's the *directness* signal net of miRNA potency.
    e = j[j["curated"]].copy()
    logcols = []
    for c in CLASSES:
        e[f"log_{c}"] = np.log1p(e[c])
        logcols.append(f"log_{c}")
    e["log_total"] = np.log1p(e["total_pmid"])
    ec = _within_family_center(e, ["repression"] + logcols + ["log_total"])
    Xb = ec[logcols].to_numpy(float)
    yb = ec["repression"].to_numpy(float)
    fit = LinearRegression(positive=True, fit_intercept=False).fit(Xb, yb)
    learned = pd.Series(fit.coef_, index=CLASSES)
    handset = pd.Series([CLASS_WEIGHT[c] for c in CLASSES], index=CLASSES)
    rho_order = spearmanr(learned.values, handset.values).correlation
    # amount-vs-directness: does class composition survive controlling total evidence? (add log_total)
    Xb2 = ec[logcols + ["log_total"]].to_numpy(float)
    fit2 = LinearRegression(positive=True, fit_intercept=False).fit(Xb2, yb)
    learned_ctrl = pd.Series(fit2.coef_[:len(CLASSES)], index=CLASSES)
    B = {"learned": learned, "learned_ctrl_total": learned_ctrl, "handset": handset,
         "spearman_order_vs_handset": float(rho_order), "n_edges": int(len(e))}

    # functional-vs-binding partial test (robust, low-assumption): within-family, does a "functional"
    # evidence score beat the "binding" (ago_clip) score at predicting repression?
    e["func_score"] = np.log1p(e[["reporter", "western", "proteomics", "qpcr_rna", "degradome"]].sum(axis=1))
    e["bind_score"] = np.log1p(e["chimeric"] + e["ago_clip"])
    ecf = _within_family_center(e, ["repression", "func_score", "bind_score"])
    rho_func = spearmanr(ecf["func_score"], ecf["repression"]).correlation
    rho_bind = spearmanr(ecf["bind_score"], ecf["repression"]).correlation
    B["rho_functional_evidence"] = float(rho_func)
    B["rho_binding_evidence"] = float(rho_bind)

    # (C) fused weight vs raw count — is weighting worth it? (within-family corr with repression)
    ec2 = _within_family_center(e, ["repression"])
    ec2["ledger_weight"] = e["ledger_weight"].values
    ec2["log_total"] = np.log1p(e["total_pmid"]).values
    ec2 = _within_family_center(ec2, ["ledger_weight", "log_total"])
    C = {"rho_ledger_weight": float(spearmanr(ec2["ledger_weight"], ec2["repression"]).correlation),
         "rho_raw_count": float(spearmanr(ec2["log_total"], ec2["repression"]).correlation)}
    return {"A": A, "B": B, "C": C, "join": j}


def _report() -> None:
    r = calibrate()
    A, B, C = r["A"], r["B"], r["C"]
    print("=" * 78)
    print("TRANSFECTION CALIBRATION of evidence-class weights (McGeary 2019 GSE140217, 16 seed families)")
    print("=" * 78)
    print(f"\n(A) Does curation capture function? curated vs non-curated repression (within-family z)")
    print(f"    curated   n={A['n_curated']:>6}  mean rep_z = {A['mean_rep_z_curated']:+.3f}")
    print(f"    non-cur.  n={A['n_noncurated']:>6}  mean rep_z = {A['mean_rep_z_noncurated']:+.3f}")
    print(f"    Mann-Whitney (curated > non-curated) p = {A['mannwhitney_p_greater']:.2e}")
    print(f"\n(B) Learned class weights vs hand-set (within-family, {B['n_edges']} curated edges)")
    tab = pd.DataFrame({"hand_set": B["handset"], "learned": B["learned"].round(4),
                        "learned_ctrl_total": B["learned_ctrl_total"].round(4)})
    tab["learned_norm"] = (tab["learned"] / tab["learned"].max()).round(2) if tab["learned"].max() > 0 else 0
    print(tab.to_string())
    print(f"    Spearman(learned order, hand-set order) = {B['spearman_order_vs_handset']:+.3f}")
    print(f"    functional-evidence ρ={B['rho_functional_evidence']:+.3f}  vs  "
          f"binding-evidence ρ={B['rho_binding_evidence']:+.3f}  (within-family; expect func > bind)")
    print(f"\n(C) Is the fused weighting worth it? within-family corr with repression")
    print(f"    fused ledger_weight ρ = {C['rho_ledger_weight']:+.3f}")
    print(f"    raw log(total PMID) ρ = {C['rho_raw_count']:+.3f}")
    print(f"\n[cache] {CACHE}")


if __name__ == "__main__":
    _report()
