"""Is the confounder block's IMMUNE control COMPLETE — or does it only look complete?

The question (and why it is NOT a panel question)
-------------------------------------------------
`build_C` controls composition with the 8 non-malignant **Wu-9** fractions, of which four are immune
(`T-cells`, `B-cells`, `Myeloid`, `Plasmablasts`). Two facts sit in tension:

* β is **panel-INVARIANT** — ρ=0.94 across entirely different atlases, ρ=0.99 across reference
  subsamples. So C's *span* is stable, and swapping the deconvolution changes nothing.
* But Wu-9's lymphoid tracks Thorsson's lymphocyte-infiltration score at only **r=0.31**, where **LM22
  reaches 0.70** — and Wu's `T-cells` is the ONE column that is technically unreliable (reference-resample
  r=0.685/0.731 vs ≥0.89 for everything else).

Panel-invariance says C's span is *stable*. It does **not** say C's span is *sufficient*. Those are
different claims, and only the second is at issue here. This is the same shape as **MH-100**, where C's
`mal_prolif` proved materially incomplete (a residual-proliferation axis carried |ρ|≤0.83 *after* C) yet a
fuller global control **over-controlled** the majority of genes ⇒ the deliverable was a per-gene FLAG, not
a change to C.

The three-part test (mirrors `prolif_verdict.py`)
-------------------------------------------------
1. **RESIDUAL** — fit `R = resid(immune_signal | C)`. If R still carries immune information, C is
   incomplete. Uses two yardsticks C never sees: **LM22** (a dedicated 22-type immune reference, TCGA
   tumours) and **Thorsson**'s published immune signature scores.
2. **LEVER** — does adding an immune term to C actually move β? (An incomplete control that doesn't move
   β is harmless.)
3. **OVER-CONTROL GUARD** — many Hallmark targets *are* immune genes (INTERFERON_GAMMA_RESPONSE,
   INFLAMMATORY_RESPONSE, …). Conditioning on immune content there removes the signal itself — the
   `data._latent` failure mode. Verdicts are therefore reported SPLIT by whether the target is an immune
   gene; a global "fix" that helps non-immune genes and destroys immune genes is not a fix.

Any fix must be **cross-cohort computable** (the state and protein channels both need a matched C), so the
candidate is a lymphocyte marker METAGENE — LM22/Thorsson exist for TCGA only and can never be covariates.

CLI: `python -m mirna_hallmark.learned.analyses.immune_completeness [--n-genes 200]`
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import gauge as G
from mirna_hallmark.learned.attribution_eb import _prep, _bagged_nnls_meansd

_LM22 = Path("data/deconvolution/cibersortx/hires_out/CIBERSORTxGEP_Job1_Fractions.txt")
_CLIN = Path("annotations/BRCA_clinical_immune_unified.tsv")

# Structural lymphocyte/leukocyte markers (pan-immune presence, NOT a signalling programme — so the
# metagene stays non-circular for immune-signalling target genes).
LYMPHO_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "LCK", "MS4A1", "CD19", "IL7R", "GZMB"]
IMMUNE_HALLMARKS = ["HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_ALLOGRAFT_REJECTION",
                    "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_COMPLEMENT",
                    "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
_CACHE: dict = {}


def _p3(ix) -> pd.Index:
    return pd.Index(["-".join(str(i).split("-")[:3]) for i in ix])


def lm22() -> Optional[pd.DataFrame]:
    """LM22 immune fractions (22 types), participant-keyed. A reference C NEVER sees."""
    if "lm22" not in _CACHE:
        if not _LM22.exists():
            _CACHE["lm22"] = None
        else:
            f = pd.read_csv(_LM22, sep="\t")
            f.index = _p3(f["Mixture"])
            f = f[~f.index.duplicated(keep="first")]
            cols = [c for c in f.columns if c not in ("Mixture", "P-value", "Correlation", "RMSE")]
            _CACHE["lm22"] = f[cols].apply(pd.to_numeric, errors="coerce")
    return _CACHE["lm22"]


def thorsson() -> Optional[pd.DataFrame]:
    """Thorsson's published immune signature scores. Also never seen by C."""
    if "thor" not in _CACHE:
        c = pd.read_csv(_CLIN, sep="\t").drop_duplicates("participant").set_index("participant")
        cols = [x for x in c.columns if x.startswith("thornsson_") and
                any(k in x for k in ("lymphocyte", "macrophage", "ifng", "tgfb", "wound"))]
        _CACHE["thor"] = c[cols].apply(pd.to_numeric, errors="coerce") if cols else None
    return _CACHE["thor"]


def lympho_metagene(parts) -> Optional[pd.Series]:
    """The CANDIDATE FIX: a structural lymphocyte marker metagene. Unlike LM22/Thorsson it is computable
    IDENTICALLY in TCGA / NAT / GTEx / CPTAC, so it can enter a matched cross-cohort C."""
    E = CF.expression("tcga")
    mk = [g for g in LYMPHO_MARKERS if g in E.index]
    sub = E.loc[mk]
    z = sub.sub(sub.mean(1), axis=0).div(sub.std(1) + 1e-9, axis=0)
    return z.mean(0).reindex(parts)


def _resid(v: np.ndarray, Cm: np.ndarray) -> np.ndarray:
    return v - LinearRegression().fit(Cm, v).predict(Cm)


def residual_immune(parts) -> pd.DataFrame:
    """PART 1 — how much immune signal SURVIVES C? Residualise each external yardstick on C and ask how
    much of it is left (R²) and whether the residual still tracks the other yardsticks."""
    C = CF.build_C("tcga", parts)
    Cm = np.c_[np.ones(len(parts)), C.to_numpy(float)]
    lm, th = lm22(), thorsson()
    probes = {}
    if lm is not None:
        ly = [c for c in lm.columns if c.startswith(("T cells", "B cells", "Plasma", "NK cells"))]
        my = [c for c in lm.columns if c.startswith(("Mono", "Macro", "Dendritic", "Mast", "Eosino", "Neutro"))]
        probes["LM22_lymphoid"] = lm[ly].sum(1).reindex(parts)
        probes["LM22_myeloid"] = lm[my].sum(1).reindex(parts)
    if th is not None:
        for c in th.columns:
            probes[c.replace("thornsson_", "TH_")] = th[c].reindex(parts)
    probes["metagene_lympho"] = lympho_metagene(parts)

    rows = []
    for name, s in probes.items():
        s = pd.to_numeric(s, errors="coerce")
        m = s.notna().to_numpy()
        if m.sum() < 100:
            continue
        v = s[m].to_numpy(float)
        r = _resid(v, Cm[m])
        rows.append({"probe": name, "n": int(m.sum()),
                     "var_left_after_C": float(np.var(r) / max(np.var(v), 1e-12)),
                     "corr_resid_vs_LM22lymph": (
                         float(spearmanr(r, _resid(probes["LM22_lymphoid"][m].to_numpy(float), Cm[m])).correlation)
                         if "LM22_lymphoid" in probes else np.nan)})
    return pd.DataFrame(rows).set_index("probe")


def beta_shift(genes=None, *, n_genes: int = 200, n_boot: int = 30) -> pd.DataFrame:
    """PARTS 2+3 — does an immune term MOVE β, and does it OVER-CONTROL immune target genes?"""
    from mirna_hallmark.learned.evidence import ledger as LG
    from mirna_hallmark.hallmark_sets import HallmarkSets
    hs = HallmarkSets.load()
    immune_genes = set().union(*[set(hs.sets.get(h, [])) for h in IMMUNE_HALLMARKS])

    X, Y = G.cohort_matrices("tcga")
    ed = LG.pooled_he_edges()
    dec = CF.deconv("tcga")
    parts = [p for p in Y.columns if p in X.columns and p in dec.index]
    genes = list(genes) if genes is not None else [g for g in ed["gene"].unique() if g in Y.index][:n_genes]

    C0 = CF.build_C("tcga", parts)
    mg = lympho_metagene(parts)
    C1 = C0.join(mg.rename("lympho"))
    D0 = pd.DataFrame(np.c_[np.ones(len(parts)), C0.to_numpy(float)], index=parts)
    D1 = pd.DataFrame(np.c_[np.ones(len(parts)), C1.fillna(C1.median()).to_numpy(float)], index=parts)

    # SCALE-FREE floor (an absolute 0.2 dropped 34% of TCGA genes; same bug class as gauge.beta_table's).
    _sds = np.array([float(Y.loc[g, parts].astype(float).std()) for g in genes if g in Y.index])
    _sds = _sds[np.isfinite(_sds) & (_sds > 0)]
    _floor = float(0.25 * np.median(_sds)) if len(_sds) else 0.0
    rows = []
    for g in genes:
        y = Y.loc[g, parts].astype(float)
        if y.std() < _floor:
            continue
        regs = [m for m in ed.loc[ed["gene"] == g, "miRNA"].unique() if m in X.index]
        if len(regs) < 2:
            continue
        Xg = X.loc[regs, parts].T.astype(float).fillna(0.0)
        Xf, _, _ = FAM.collapse_by_family(Xg, pd.Series(1.0, index=regs), FAM.family_of(pd.Index(regs)))
        if Xf.shape[1] < 1:
            continue
        b0 = _bagged_nnls_meansd(*_prep(y, Xf, D0)[:2][::-1], n_boot=n_boot)[0]
        b1 = _bagged_nnls_meansd(*_prep(y, Xf, D1)[:2][::-1], n_boot=n_boot)[0]
        cols = list(Xf.columns)
        for c, v0, v1 in zip(cols, b0, b1):
            rows.append({"gene": g, "family": c, "beta_base": float(v0), "beta_immune": float(v1),
                         "immune_gene": g in immune_genes})
    return pd.DataFrame(rows)


def run(n_genes: int = 200) -> None:
    X, Y = G.cohort_matrices("tcga")
    dec = CF.deconv("tcga")
    parts = [p for p in Y.columns if p in X.columns and p in dec.index]
    print(f"[immune] TCGA tumours: {len(parts)}")

    print("\n" + "=" * 88)
    print("PART 1 — RESIDUAL: how much immune signal SURVIVES the current C?")
    print("  var_left_after_C ≈ 1.0 ⇒ C removes almost none of it ⇒ C is immune-INCOMPLETE")
    print("=" * 88)
    with pd.option_context("display.width", 130):
        print(residual_immune(parts).round(3).to_string())

    print("\n" + "=" * 88)
    print("PART 2+3 — LEVER (does an immune term move β?) and OVER-CONTROL guard")
    print("=" * 88)
    b = beta_shift(n_genes=n_genes)
    if b.empty:
        print("  no edges"); return
    for label, sub in [("ALL", b), ("immune target genes", b[b.immune_gene]),
                       ("non-immune target genes", b[~b.immune_gene])]:
        if len(sub) < 10:
            continue
        nz = (sub[["beta_base", "beta_immune"]] > 1e-8).all(1)
        rho = spearmanr(sub.beta_base, sub.beta_immune).correlation
        rho_nz = spearmanr(sub.loc[nz, "beta_base"], sub.loc[nz, "beta_immune"]).correlation if nz.sum() > 10 else np.nan
        d = (sub.beta_immune - sub.beta_base)
        print(f"  {label:<24} n={len(sub):<5} genes={sub.gene.nunique():<4} ρ={rho:.3f} ({rho_nz:.3f} nonzero)  "
              f"mean|β| {sub.beta_base.abs().mean():.4f}→{sub.beta_immune.abs().mean():.4f}  "
              f"median Δ={d.median():+.5f}")
    print("\n  ρ≈1 ⇒ the immune term is INERT (C already sufficient) · ρ low + |β| DOWN on immune genes ⇒ "
          "OVER-CONTROL · ρ low + |β| up ⇒ genuine de-confounding")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-genes", type=int, default=200)
    run(ap.parse_args().n_genes)
