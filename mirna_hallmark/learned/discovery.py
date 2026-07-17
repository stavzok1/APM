"""Discovery lane (Design §Decision D/E soft-discovery) — where the prior actually pays off.

Bar 3 showed the prior's value is INCLUSION (which arms), not magnitude. Curation's inclusion is
incomplete: the heavily-studied edges are often buffered in bulk breast while less-studied breast-active
edges carry the cross-sample variance. So use the *sequence* inclusion prior (TargetScan-predicted sites,
`assemble_gene(orphans=True)`) to NOMINATE arms outside the curated HE set, and keep only those that add
**repressive coupling beyond the curated aggregate** with **biochemical support** (scanMiR K_D).

Criterion for a discovery (all three):
  1. orphan (TargetScan site, NOT in curated HE for this gene);
  2. **partial** repressive coupling — corr(arm, Y | C, HE-aggregate) < threshold  (adds signal curation missed);
  3. scanMiR predicted repression < 0 (a biochemical duplex, not a spurious correlation).
Weak sub-HE ledger evidence (if any) is reported as corroboration, not required.

CLI: `python -m mirna_hallmark.learned.discovery`  (scan a panel, rank cohort-wide discoveries)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import kd as KD
from mirna_hallmark.learned import regression as LR
from mirna_hallmark.learned.evidence import ledger as LG

PANEL = ["PTEN", "ESR1", "ZEB1", "CDKN1A", "GATA3", "BCL2", "MYC", "CDKN1B", "RB1", "SMAD4"]


def _resid(v, Cmat):
    return v - LinearRegression().fit(Cmat, v).predict(Cmat)


def discover_gene(gene: str, *, min_partial: float = -0.12, alpha: float = 0.01,
                  min_scanmir: float = 0.0, permute: int | None = None,
                  deconv_check: bool = True) -> pd.DataFrame:
    """Orphan edges that add repressive coupling to `gene` beyond the curated HE aggregate + scanMiR support.
    `permute` shuffles the target (breaks the miRNA→target relation) → the FDR null. `min_scanmir` requires
    scanMiR predicted repression below it (stringency)."""
    try:
        Y, Xo, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True)
        _, Xhe, _, whe = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return pd.DataFrame()
    if permute is not None:
        Y = pd.Series(np.random.default_rng(permute).permutation(Y.to_numpy()), index=Y.index)
    he = set(Xhe.columns)
    orphans = [a for a in Xo.columns if a not in he]
    if not orphans:
        return pd.DataFrame()
    Cm = C.to_numpy(float)
    # curated HE aggregate (fit on HE only), added to C so orphan coupling is PARTIAL on what curation captures
    M_he = LR.fit_gene(Y, Xhe, C, whe, alpha=alpha)
    he_agg = Xhe.to_numpy(float) @ M_he.reindex(Xhe.columns).fillna(0).to_numpy()
    Cext = np.column_stack([Cm, he_agg])
    yr = _resid(Y.to_numpy(float), Cext)
    aff = KD.genome_affinity()      # genome-wide (unbiased) biochemical support — was HE-restricted KD.affinity()
    affg = aff.loc[aff["gene"] == gene].set_index("arm")["repression"]
    lw = LG.edge_weights(); lwg = lw.loc[lw["gene"] == gene].set_index("arm")["ledger_weight"]
    rows = []
    for a in orphans:
        xr = _resid(Xo[a].to_numpy(float), Cext)
        rho = spearmanr(xr, yr).correlation
        rows.append({"gene": gene, "arm": a, "partial_coupling": round(float(rho), 3),
                     "scanmir_rep": round(float(affg.get(a, np.nan)), 2) if a in affg.index else np.nan,
                     "sub_he_evidence": round(float(lwg.get(a, 0.0)), 2)})
    df = pd.DataFrame(rows)
    keep = df[(df["partial_coupling"] < min_partial) & (df["scanmir_rep"] < min_scanmir)].copy()
    if keep.empty or not deconv_check:
        return keep
    # composition control: recompute partial coupling with CIBERSORTx non-malignant fractions in C.
    # A real edge survives; a composition artefact (e.g. luminal-fraction-driven) collapses.
    try:
        Yd, Xod, Cd, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True, deconv=True)
        if permute is not None:
            Yd = pd.Series(np.random.default_rng(permute).permutation(Yd.to_numpy()), index=Yd.index)
        _, Xhed, _, whed = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
        Cmd = Cd.to_numpy(float)
        aggd = Xhed.to_numpy(float) @ LR.fit_gene(Yd, Xhed, Cd, whed, alpha=alpha).reindex(Xhed.columns).fillna(0).to_numpy()
        Cextd = np.column_stack([Cmd, aggd])
        yrd = _resid(Yd.to_numpy(float), Cextd)
        keep["partial_deconv"] = [round(float(spearmanr(_resid(Xod[a].to_numpy(float), Cextd), yrd).correlation), 3)
                                  if a in Xod.columns else np.nan for a in keep["arm"]]
    except Exception:
        keep["partial_deconv"] = np.nan
    keep["robust"] = keep["partial_deconv"] < min_partial * 0.6           # retains >60% of signal under deconv
    return keep.sort_values("partial_coupling")


_BATCH: dict = {}


def _batch_np(parts) -> np.ndarray:
    """plate_both batch dummies (post-fit conditioning only) for the given participants, cached full-cohort."""
    if "bd" not in _BATCH:
        from mirna_hallmark import tcga_batch as TB
        _BATCH["bd"] = TB.tcga_batch_dummies(list(LD._load()["X"].columns), kind="plate_both")
    bd = _BATCH["bd"].reindex(list(parts)).fillna(0.0)
    return bd.loc[:, bd.sum(axis=0) > 0].to_numpy(float)


def _he_arms_map() -> dict:
    """gene → set of curated POOLED-HE regulator arms present in the arm matrix (one build, avoids a 2nd
    assemble). POOLED-HE (migration) = the 'known' baseline, so orphans are edges beyond it."""
    he = LG.pooled_he_edges()
    xidx = set(LD._load()["X"].index)
    return {str(g): set(sub["miRNA"]) & xidx for g, sub in he.groupby("gene")}


def scan_all(genes=None, *, permute: int | None = None, min_partial: float = -0.15,
             min_scanmir: float = -0.5, alpha: float = 0.01, progress: int = 200,
             batch: bool = True) -> pd.DataFrame:
    """ONE assemble per gene: orphan partial coupling (beyond curated HE aggregate) + scanMiR prefilter.
    No deconv (that's applied to survivors afterwards). `permute` shuffles the target → the FDR null."""
    genes = genes or sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))  # POOLED-HE gene universe (migration)
    hemap = _he_arms_map()
    aff = KD.genome_affinity(); lw = LG.edge_weights()      # genome-wide biochemical support (was HE-restricted)
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[scan{'·null' if permute else ''}] {i}/{len(genes)} genes, {len(rows)} hits", flush=True)
        try:
            Y, Xo, C, w = LD.assemble_gene(g, w_prior_source="ledger", orphans=True)
        except Exception:
            continue
        if permute is not None:
            Y = pd.Series(np.random.default_rng(permute + i).permutation(Y.to_numpy()), index=Y.index)
        he_arms = [a for a in hemap.get(g, set()) if a in Xo.columns]
        orphans = [a for a in Xo.columns if a not in set(he_arms)]
        if not orphans or len(he_arms) < 1:
            continue
        affg = aff.loc[aff["gene"] == g].set_index("arm")["repression"]
        cand = [a for a in orphans if a in affg.index and affg[a] < min_scanmir]   # scanMiR prefilter (cheap)
        if not cand:
            continue
        Cm = C.to_numpy(float)
        Xhe = Xo[he_arms]
        M_he = LR.fit_gene(Y, Xhe, C, w.reindex(he_arms), alpha=alpha)
        he_agg = Xhe.to_numpy(float) @ M_he.reindex(Xhe.columns).fillna(0).to_numpy()
        Cext = np.column_stack([Cm, he_agg, _batch_np(Y.index)]) if batch else np.column_stack([Cm, he_agg])
        yr = _resid(Y.to_numpy(float), Cext)
        lwg = lw.loc[lw["gene"] == g].set_index("arm")["ledger_weight"]
        for a in cand:
            rho = spearmanr(_resid(Xo[a].to_numpy(float), Cext), yr).correlation
            if rho < min_partial:
                rows.append({"gene": g, "arm": a, "partial_coupling": round(float(rho), 3),
                             "scanmir_rep": round(float(affg[a]), 2),
                             "sub_he_evidence": round(float(lwg.get(a, 0.0)), 2)})
    return pd.DataFrame(rows)


def deconv_validate(cand: pd.DataFrame, *, min_partial: float = -0.15, alpha: float = 0.01,
                    progress: int = 50, batch: bool = True) -> pd.DataFrame:
    """Composition-robustness for candidate (gene, arm) edges — **ONE deconv assemble per gene**, validating
    all of that gene's candidates in the same pass (not a re-`discover_gene` per gene). Adds `partial_deconv`
    + `robust` (retains >60% of the coupling under CIBERSORTx non-malignant fractions)."""
    out = []
    for i, (g, grp) in enumerate(cand.groupby("gene")):
        if progress and i % progress == 0:
            print(f"[deconv] {i}/{cand['gene'].nunique()} genes", flush=True)
        try:
            Yd, Xod, Cd, _ = LD.assemble_gene(g, w_prior_source="ledger", orphans=True, deconv=True)
            _, Xhed, _, whed = LD.assemble_gene(g, w_prior_source="ledger", deconv=True)
        except Exception:
            continue
        Cmd = Cd.to_numpy(float)
        aggd = Xhed.to_numpy(float) @ LR.fit_gene(Yd, Xhed, Cd, whed, alpha=alpha).reindex(Xhed.columns).fillna(0).to_numpy()
        Cextd = np.column_stack([Cmd, aggd, _batch_np(Yd.index)]) if batch else np.column_stack([Cmd, aggd])
        yrd = _resid(Yd.to_numpy(float), Cextd)
        for _, r in grp.iterrows():
            a = r["arm"]
            pdc = (spearmanr(_resid(Xod[a].to_numpy(float), Cextd), yrd).correlation
                   if a in Xod.columns else np.nan)
            row = r.to_dict()
            row["partial_deconv"] = round(float(pdc), 3) if pdc == pdc else np.nan
            ret = (pdc / r["partial_coupling"]) if r["partial_coupling"] else np.nan
            row["retention"] = round(float(ret), 2) if ret == ret else np.nan
            # composition-robust = retains ≥60% of coupling under deconv AND still couples (not a stroma edge)
            row["robust"] = bool(ret >= 0.6 and pdc < -0.1)
            out.append(row)
    return pd.DataFrame(out).sort_values("partial_coupling") if out else pd.DataFrame()


def run_all(*, min_partial: float = -0.15, min_scanmir: float = -0.5, top: int | None = None,
            out="mirna_hallmark/output/learned/discoveries.tsv") -> pd.DataFrame:
    """Genome-wide discovery: single-pass real + null scans → FDR, persist candidates, then single-pass
    deconv-validation (`top` limits validation to the strongest, else all)."""
    from pathlib import Path
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    real = scan_all(min_partial=min_partial, min_scanmir=min_scanmir)
    real.to_csv(Path(out).parent / "discovery_candidates.tsv", sep="\t", index=False)   # persist BEFORE deconv
    null = scan_all(min_partial=min_partial, min_scanmir=min_scanmir, permute=101)
    fdr = len(null) / max(len(real), 1)
    print(f"[run_all] real hits={len(real)}  null hits={len(null)}  edge-level FDR≈{fdr:.3f}  (candidates persisted)")
    cand = real.nsmallest(top, "partial_coupling") if top else real
    dv = deconv_validate(cand, min_partial=min_partial)                                  # single-pass
    if len(dv):
        dv.to_csv(out, sep="\t", index=False)
        rob = dv[dv["robust"]] if "robust" in dv else dv
        print(f"[run_all] {len(dv)} validated → {len(rob)} deconv-robust across {rob['gene'].nunique()} genes "
              f"| {int((rob['sub_he_evidence']==0).sum())} fully novel | wrote {out}")
        print(rob.sort_values('partial_deconv').head(15).to_string(index=False))
    return dv


def run(genes=None, *, min_partial: float = -0.12) -> pd.DataFrame:
    genes = genes or PANEL
    out = pd.concat([discover_gene(g, min_partial=min_partial) for g in genes], ignore_index=True) \
        if genes else pd.DataFrame()
    if out.empty:
        print("no discoveries at this threshold")
        return out
    out = out.sort_values("partial_coupling")
    with pd.option_context("display.width", 160):
        print(f"=== Discovered orphan edges (partial coupling < {min_partial} beyond curated HE + scanMiR-supported) ===")
        print(out.to_string(index=False))
    nr = int(out["robust"].sum()) if "robust" in out else 0
    print(f"\n{len(out)} candidates across {out['gene'].nunique()} genes | {nr} SURVIVE composition control "
          f"(deconv) | {int((out['sub_he_evidence'] == 0).sum())} fully novel (sequence+expression+K_D only)")
    if nr:
        print("robust discoveries:", ", ".join(f"{r.arm}→{r.gene}" for r in out[out['robust']].itertuples()))
    return out


if __name__ == "__main__":
    run()
