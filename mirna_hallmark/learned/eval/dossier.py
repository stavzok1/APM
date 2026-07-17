"""Per-edge DISCOVERY DOSSIER over the combined TS∪K_D consensus (6,744 edges). Turn each discovered edge
from 'N methods agreed' into a fully-characterized candidate. Tiered:

  TIER 1 (all edges, cheap join): method membership, candidacy source, scanMiR duplex, sub-HE PMIDs.
  TIER 2 (all edges, moderate per-gene compute): the REALIZED coupling — partial-Spearman(arm, target |
          C + HE-aggregate), the 'does this edge actually repress' number; + arm detectability (median RPM).
  TIER 3 (gold/cross-candidacy only, expensive): CPTAC protein coupling (if gene in CPTAC).

    python -m mirna_hallmark.learned.eval.dossier
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import regression as LR

COMBINED = "mirna_hallmark/output/learned/discovery_consensus_combined.tsv"
OUT = "mirna_hallmark/output/learned/discovery_dossier.tsv"


def _resid(v, Z):
    return v - LinearRegression().fit(Z, v).predict(Z)


def _one_block(gene: str, arms: list, deconv: bool) -> dict:
    """Realized partial coupling per arm under ONE confounder block: Spearman(resid(arm|C+HE-agg),
    resid(target|C+HE-agg))."""
    out = {}
    try:                                                          # assemble both candidacy sources → all arms present
        Y, Xt, C, wt = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True,
                                        orphan_source="targetscan", deconv=deconv)
        _, Xk, _, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True,
                                       orphan_source="kd", deconv=deconv)
        _, Xhe, _, whe = LD.assemble_gene(gene, w_prior_source="ledger", deconv=deconv)
    except Exception:
        return {}
    X = pd.concat([Xt, Xk.loc[:, [c for c in Xk.columns if c not in Xt.columns]]], axis=1)
    Cm = C.to_numpy(float)
    he_agg = Xhe.to_numpy(float) @ LR.fit_gene(Y, Xhe, C, whe).reindex(Xhe.columns).fillna(0).to_numpy()
    Cext = np.column_stack([Cm, he_agg])
    yr = _resid(Y.to_numpy(float), Cext)
    for a in arms:
        if a not in X.columns:
            continue
        xv = X[a].to_numpy(float)
        xr = _resid(xv, Cext)
        rho = spearmanr(xr, yr).correlation if np.std(xr) > 1e-6 else np.nan
        out[a] = {"rho": float(rho) if rho == rho else np.nan, "det": float(np.median(xv))}
    return out


def _gene_edges_coupling(gene: str, arms: list) -> dict:
    """Realized partial coupling per discovered arm, under **BOTH** confounder blocks, + the RETENTION tag.

    ⚠⚠ **MH-119 — THIS LANE WAS COMPOSITION-BLIND.** It called `assemble_gene` with the DEFAULT C
    (`deconv=False`: CPE + target_cn + mal_prolif), i.e. purity but **NO cell-composition block**, and its
    `realized_coupling` drives `couples` and `dossier_pass` (2,522 edges). MH-114 showed bulk mRNA carries the
    SAME compartment artifact as bulk protein — degree-preserving SHUFFLED non-edges reproduce the real edges'
    compartment gradient EXACTLY (−0.1355 vs −0.1389) — and CPE does NOT remove it (purity is only ~1% of the
    loss; the stromal MIX is the rest). So an edge could earn `dossier_pass` on compartment arithmetic alone.

    THE FIX IS **NOT** TO CONDITION (axiom 2a: a miRNA acting via a cell-STATE shift PRODUCES a composition
    change, so blanket conditioning OVER-CONTROLS). It is to run **BOTH** blocks and report
    `coupling_retention = rho_deconv / rho_core`, classed as `card.py` — **FLAG, DON'T DELETE.**
    """
    core = _one_block(gene, arms, deconv=False)                   # the design's canonical block
    dec = _one_block(gene, arms, deconv=True)                     # + the 8 Wu-major non-malignant lineages
    out = {}
    for a in arms:
        c = core.get(a)
        if not c:
            out[a] = {}
            continue
        d = dec.get(a, {})
        rc, rd = c["rho"], d.get("rho", np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            ret = rd / rc if (rc == rc and abs(rc) > 1e-9 and rd == rd) else np.nan
        cls = np.nan
        if ret == ret:
            cls = ("composition_explained" if ret < 0.4 else
                   "partial" if ret < 0.7 else "cell_intrinsic")
        out[a] = {"realized_coupling": round(rc, 3) if rc == rc else np.nan,
                  "realized_coupling_deconv": round(rd, 3) if rd == rd else np.nan,
                  "coupling_retention": round(float(ret), 3) if ret == ret else np.nan,
                  "composition_class": cls,
                  "arm_median_rpm": round(c["det"], 2),
                  "arm_expressed": bool(c["det"] > np.log2(11))}
    return out


_ARMS: dict = {}          # gene -> arms; module-level so forked workers inherit it (the `readouts.run` pattern)


def _one_gene(g: str):
    """One gene's coupling, for the worker pool. Returns (gene, {arm: {...}})."""
    try:
        return g, _gene_edges_coupling(g, _ARMS.get(g, []))
    except Exception:
        return g, {}


def build(*, limit: int | None = None, progress: int = 100, workers: int = 8,
          out: str | None = None) -> pd.DataFrame:
    """⚡ PARALLEL over genes (2026-07-17). The loop was serial at ~1.03 s/gene × 1,344 genes ≈ 23 min, and
    it is embarrassingly parallel — each gene's `_gene_edges_coupling` is independent (it re-assembles its
    own design). Same worker-pool pattern as `readouts.run`.

    ⚠ **`limit=` NO LONGER WRITES TO THE PRODUCTION PATH.** It used to: a `build(limit=40)` sizing run
    silently overwrote `discovery_dossier.tsv` with 40 rows (it did, to me, 2026-07-17 — and `output/` is
    gitignored, so there was no diff to catch it). A truncated run is a TEST; it must not clobber a
    production artifact. `limit` now forces a `.limitN` suffix unless `out=` is given explicitly.
    """
    c = pd.read_csv(COMBINED, sep="\t")
    if limit:
        c = c.head(limit)
    dest = out or (OUT if not limit else OUT.replace(".tsv", f".limit{limit}.tsv"))
    # bring scanMiR + sub-HE from the source consensus tables
    ts = pd.read_csv("mirna_hallmark/output/learned/discovery_consensus_targetscan.tsv", sep="\t")
    ts["edge"] = list(zip(ts.gene, ts.arm))
    subd = ts.drop_duplicates("edge").set_index("edge")["sub_he_evidence"].to_dict()  # dict avoids .loc[tuple] ambiguity
    genes = sorted(c["gene"].unique())
    _ARMS.clear()
    _ARMS.update({g: list(c.loc[c.gene == g, "arm"]) for g in genes})
    cpl = {}
    if workers and workers > 1:
        import multiprocessing as mp
        with mp.Pool(workers) as pool:
            for i, (g, res) in enumerate(pool.imap_unordered(_one_gene, genes, chunksize=4)):
                if progress and i % progress == 0:
                    print(f"[dossier] {i}/{len(genes)} genes", flush=True)
                cpl.update({(g, a): v for a, v in res.items()})
    else:
        for i, g in enumerate(genes):
            if progress and i % progress == 0:
                print(f"[dossier] {i}/{len(genes)} genes", flush=True)
            cpl.update({(g, a): v for a, v in _gene_edges_coupling(g, _ARMS[g]).items()})
    # assemble the dossier rows
    rows = []
    for _, r in c.iterrows():
        e = (r.gene, r.arm)
        cp = cpl.get(e, {})
        she = subd.get(e, np.nan)
        rows.append({
            "gene": r.gene, "arm": r.arm,
            "in_targetscan": r.in_targetscan, "in_kd": r.in_kd, "xcandidacy": r.xcandidacy,
            "ts_nmethods": int(r.ts_nmethods), "kd_nbayes": int(r.kd_nbayes),
            "scanmir_rep": r.scanmir_rep,
            "sub_he_evidence": float(she) if pd.notna(she) else np.nan,
            "realized_coupling": cp.get("realized_coupling", np.nan),
            "realized_coupling_deconv": cp.get("realized_coupling_deconv", np.nan),
            "coupling_retention": cp.get("coupling_retention", np.nan),
            "composition_class": cp.get("composition_class", np.nan),
            "arm_median_rpm": cp.get("arm_median_rpm", np.nan),
            "arm_expressed": cp.get("arm_expressed", np.nan),
        })
    d = pd.DataFrame(rows)
    # a headline "confidence" flag: couples (ρ<-0.1) AND expressed AND scanMiR duplex
    d["couples"] = d["realized_coupling"] < -0.1
    # MH-119: the coupling must ALSO survive the cell-composition block. An edge whose bulk coupling is
    # `composition_explained` is compartment arithmetic (MH-112: 90.5% of them are OPPOSITE-compartment
    # miRNA/target pairs, Fisher OR=0.09, p=1.6e-87) and must NOT earn a pass. `composition_class` is
    # reported for EVERY edge either way — the gate only tightens `dossier_pass`, nothing is deleted.
    d["couples_cell_intrinsic"] = d["couples"] & (d["composition_class"] != "composition_explained")
    d["dossier_pass"] = (d["couples_cell_intrinsic"] & (d["arm_expressed"] == True)
                         & (d["scanmir_rep"] < 0)).fillna(False)
    d.to_csv(dest, sep="\t", index=False)
    print(f"\n=== DISCOVERY DOSSIER: {len(d)} edges → {dest} ===")
    print(f"  realized coupling ρ<−0.1: {int(d['couples'].sum())} ({100*d['couples'].mean():.0f}%) | "
          f"arm expressed: {int((d['arm_expressed']==True).sum())} | "
          f"**dossier_pass (couples+expressed+duplex): {int(d['dossier_pass'].sum())} ({100*d['dossier_pass'].mean():.0f}%)**")
    print(f"  mean realized coupling: {d['realized_coupling'].mean():.3f} | median {d['realized_coupling'].median():.3f}")
    cc = d["composition_class"].value_counts(dropna=True)
    n_ce = int(cc.get("composition_explained", 0))
    lost = int((d["couples"] & (d["composition_class"] == "composition_explained")).sum())
    print(f"  COMPOSITION (MH-119, both C blocks): cell_intrinsic {int(cc.get('cell_intrinsic',0))} | "
          f"partial {int(cc.get('partial',0))} | composition_explained {n_ce} "
          f"({100*n_ce/max(d['composition_class'].notna().sum(),1):.0f}%)")
    print(f"  ⇒ edges that COUPLE but are COMPARTMENT ARITHMETIC (now denied a pass): {lost} "
          f"({100*lost/max(int(d['couples'].sum()),1):.0f}% of couplers)")
    # by source
    for src, m in [("TargetScan", d.in_targetscan), ("K_D", d.in_kd), ("cross-candidacy", d.xcandidacy)]:
        sub_d = d[m]
        print(f"  {src:14s}: {len(sub_d)} edges | couples {100*sub_d['couples'].mean():.0f}% | "
              f"dossier_pass {100*sub_d['dossier_pass'].mean():.0f}% | mean ρ {sub_d['realized_coupling'].mean():.3f}")
    return d


def tier3_protein(*, cohort: str = "prospective",
                  out: str = "mirna_hallmark/output/learned/discovery_dossier_gold.tsv") -> pd.DataFrame:
    """TIER 3 (gold sets): per-edge CPTAC PROTEIN coupling — the independent layer. For each gold edge,
    partial Spearman(arm abundance, target protein) and (arm, protein-vs-mRNA discordance) on the INDEPENDENT
    prospective CPTAC cohort. Restricted to the gold/ultragold/cross-candidacy edges where it matters.

    ⚠⚠ **CELL-COMPOSITION-ADJUSTED since MH-107 (2026-07-12). This was a RAW Spearman with NO covariates**, while
    the TCGA model these edges come from already conditions on the 8 Wu-major lineages. Bulk protein is heavily
    compartment-loaded (EMT proteins ↔ CAF fraction r=+0.509; ZEB1 protein +0.768) and epithelial miRNAs (miR-200)
    anti-correlate with stroma — which MANUFACTURES protein "coupling" out of compartment arithmetic.
    **Measured on these very gold edges: of 353 called protein-coupled by the RAW test, 259 (73%) keep their sign
    under composition but 94 (27%) SIGN-FLIP, and the mean |ρ| falls 60% (−0.190 → −0.075).** So ~1/4 of the
    "268 TRIPLE-VALIDATED" edges (MH-84) were composition artifacts. `protein_coupling_RAW` keeps the historical
    value for provenance."""
    from mirna_hallmark.learned.eval import ood_protein as OP
    d = pd.read_csv(OUT, sep="\t")
    d["edge"] = list(zip(d.gene, d.arm))
    ts = pd.read_csv("mirna_hallmark/output/learned/discovery_consensus_targetscan.tsv", sep="\t")
    ts["edge"] = list(zip(ts.gene, ts.arm))
    gold = set(ts[ts.in_partial & ts.in_learned & ts.in_densez].edge) | set(ts[ts.n_methods == 5].edge)
    gold |= set(d.loc[d.xcandidacy, "edge"])
    cp = OP._cptac(cohort)
    arms, prot, resid = cp["arms"], cp["protein"], cp["resid"]
    s_all = [c for c in arms.columns if c in prot.columns]
    D = OP._cov(cohort, s_all)                                     # MH-107: composition block

    def _sp(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return round(float(spearmanr(a[m], b[m]).correlation), 3) if m.sum() > 15 else np.nan
    rows = []
    for e in gold:
        g, a = e
        if a not in arms.index or g not in prot.index:
            rows.append({"gene": g, "arm": a, "prot_n": 0, "protein_coupling": np.nan, "discord_coupling": np.nan})
            continue
        s = s_all
        xv = arms.loc[a, s].astype(float).to_numpy()
        pv = prot.loc[g, s].astype(float).to_numpy()
        rv = resid.loc[g, s].astype(float).to_numpy()
        xr, pr_, rr = OP._resid(xv, D), OP._resid(pv, D), OP._resid(rv, D)
        rows.append({"gene": g, "arm": a, "prot_n": len(s),
                     "protein_coupling": _sp(xr, pr_),             # composition-adjusted (the value to trust)
                     "discord_coupling": _sp(xr, rr),
                     "protein_coupling_RAW": _sp(xv, pv),          # historical/confounded, provenance only
                     "discord_coupling_RAW": _sp(xv, rv)})
    t3 = pd.DataFrame(rows)
    m = d.merge(t3, on=["gene", "arm"], how="right")
    m.to_csv(out, sep="\t", index=False)
    ok = m[m["protein_coupling"].notna()]
    print(f"\n=== TIER 3: CPTAC {cohort} protein coupling on {len(m)} gold edges ({len(ok)} testable) ===")
    print(f"  protein-coupled (ρ<0): {int((ok['protein_coupling']<0).sum())}/{len(ok)} ({100*(ok['protein_coupling']<0).mean():.0f}%) "
          f"| mean ρ_protein {ok['protein_coupling'].mean():+.3f}")
    print(f"  discordance-coupled (ρ<0): {int((ok['discord_coupling']<0).sum())}/{len(ok)} "
          f"({100*(ok['discord_coupling']<0).mean():.0f}%) — beyond-mRNA translational signal")
    # triple validation: mRNA coupling (dossier) + protein + literature
    tri = ok[(ok["realized_coupling"] < -0.1) & (ok["protein_coupling"] < 0) & (ok["sub_he_evidence"] > 0)]
    print(f"  ★ TRIPLE-VALIDATED (mRNA-couples + protein-couples + PMIDs): {len(tri)} edges")
    print(tri.sort_values("protein_coupling").head(10)[
        ["gene", "arm", "realized_coupling", "protein_coupling", "sub_he_evidence"]].to_string(index=False))
    return m


if __name__ == "__main__":
    import sys
    if "--tier3" in sys.argv:
        tier3_protein()
    else:
        w = 8
        for a in sys.argv[1:]:
            if a.startswith("--workers"):
                w = int(a.split("=")[1]) if "=" in a else 8
        build(workers=w)
