"""Nested state model Mₜ = M_H + Δ (Design §Decision H, §6 Phase 6). Is the per-site weight `M`
state-dependent, or is only abundance `x`?

Fit M_H on **healthy GTEx breast** (paired miRNA→gene, 327 donors) and M_T on **tumour** (TCGA), both on
**within-cohort z-scored abundance** so the two matrices share a gauge (Decision C: separately-gauged M are
not comparable). Then operationalise the Δ=0 test as **cross-state transfer**:

    ρ(M_H → tumour)   vs   ρ(M_T → tumour)      # do HEALTHY weights couple tumour mRNA as well as tumour weights?

If healthy-fit weights predict tumour repression as well as tumour-fit weights (ρ_H ≈ ρ_T), then Δ≈0 — the
regulatory *weight* is stable and the tumour story is **acquired abundance** (Δx), not acquired per-site
weight (the simpler world, Design §7 Q1). If ρ_T ≫ ρ_H, `M` is genuinely state-dependent (APA site loss /
AGO-capacity shift). Reports the per-arm Δ = M_T − M_H so the *direction* of any state change is inspectable.

Data: GTEx v10 breast `gene_tpm` (mRNA) + `gtex_mirna_matrix` (arms on TCGA names), donor-paired.
CLI: `python -m mirna_hallmark.learned.state`
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Lasso, LinearRegression

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM

GENE_TPM = Path("data/GTEx/gene_tpm_v10_breast.gct.gz")
_CACHE: dict = {}


def _gtex_mrna() -> pd.DataFrame:
    """GTEx breast gene TPM → log2(TPM+1), gene symbol × donor (donor = GTEX-XXXX)."""
    if "mrna" not in _CACHE:
        df = pd.read_csv(GENE_TPM, sep="\t", skiprows=2)
        df = df.drop(columns=["Name"]).groupby("Description").max()          # symbol × sample
        df.columns = ["-".join(str(c).split("-")[:2]) for c in df.columns]   # aliquot → donor
        df = df.T.groupby(level=0).mean().T                                   # dedupe donor cols
        _CACHE["mrna"] = np.log2(df + 1.0)
    return _CACHE["mrna"]


def _gtex_mirna() -> pd.DataFrame:
    """GTEx breast miRNA arm × donor on TCGA arm names → log2(TPM+1)."""
    if "mirna" not in _CACHE:
        from mirna_hallmark import gtex_mirna_matrix as G
        m = G.gtex_tcga_arm_matrix().copy()
        m.columns = ["-".join(str(c).split("-")[:2]) for c in m.columns]
        m = m.T.groupby(level=0).mean().T
        _CACHE["mirna"] = np.log2(m + 1.0)
    return _CACHE["mirna"]


def _zt(df: pd.DataFrame) -> pd.DataFrame:
    """Column-wise z-score (shared abundance gauge across cohorts)."""
    return (df - df.mean()) / df.std(ddof=0).replace(0, np.nan)


def _fit_M(Y: pd.Series, Xz: pd.DataFrame, *, alpha: float = 0.01,
           C: "pd.DataFrame | None" = None) -> pd.Series:
    """Non-negative lasso of −resid(Y | C) on z-scored abundance → per-arm repression weight M (≥0).

    ⚠ **`C` WAS MISSING ENTIRELY (fixed 2026-07-12).** This function took only (Y, Xz) and merely
    mean-centred Y — so **M_H was fit on GTEx with NO confounder block at all**, while the tumour target it
    is scored against IS C-residualised. Composition alone HALVES |β| (MH-102: no-C vs C, ρ=0.58), so roughly
    half of the "healthy wiring" this reported was composition artifact. `confounders.build_C("gtex")` now
    exists (514 donors, covering 327/327 paired) — pass it."""
    yv = Y.to_numpy(float)
    yr = yv - (LinearRegression().fit(C.to_numpy(float), yv).predict(C.to_numpy(float))
               if C is not None and C.shape[1] else yv.mean())
    Xf = Xz.fillna(0.0).to_numpy(float)
    m = Lasso(alpha=alpha, positive=True, max_iter=5000).fit(Xf, -yr)
    return pd.Series(m.coef_, index=Xz.columns)


def cross_state_transfer(gene: str, *, alpha: float = 0.01, family: bool = True) -> dict:
    """Fit M_H (GTEx) and M_T (tumour) on shared-gauge z-scored abundance; test whether healthy weights
    transfer to tumour coupling (Δ=0) and report the per-arm state change Δ = M_T − M_H."""
    Yt, Xt, Ct, _ = LD.assemble_gene(gene)
    mrna, mir = _gtex_mrna(), _gtex_mirna()
    if gene not in mrna.index:
        return {"gene": gene, "error": "gene not in GTEx mRNA"}
    donors = [d for d in mrna.columns if d in mir.columns]
    arms = [a for a in Xt.columns if a in mir.index]                         # HE arms present in both
    if len(arms) < 3 or len(donors) < 40:
        return {"gene": gene, "error": f"overlap arms={len(arms)} donors={len(donors)}"}
    Yh = mrna.loc[gene, donors]
    Xh = mir.loc[arms, donors].T                                             # donor × arm
    Xt = Xt[arms]
    if family:                                                              # collapse both to seed families
        fam = FAM.family_of(pd.Index(arms))
        Xt, _, _ = FAM.collapse_by_family(Xt, pd.Series(1.0, index=arms), fam)
        Xh, _, _ = FAM.collapse_by_family(Xh, pd.Series(1.0, index=arms), fam)
        common = [c for c in Xt.columns if c in Xh.columns]
        Xt, Xh = Xt[common], Xh[common]

    Xtz, Xhz = _zt(Xt), _zt(Xh)
    ng = len(donors)
    Cm = Ct.to_numpy(float)                                                  # tumour coupling residualises on C
    ytr = Yt.to_numpy(float) - LinearRegression().fit(Cm, Yt.to_numpy(float)).predict(Cm)
    Xtn = Xtz.fillna(0.0).to_numpy(float)

    # n-matched, held-out control: split tumour into a 327-fit set and the rest; evaluate BOTH the
    # in-state (tumour-subsample) M and the out-of-state (healthy) M on the SAME held-out tumour. Both fits
    # use n=ng, so any coupling gap between them is STATE, not estimation noise.
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(Yt))
    fit_i, test_i = perm[:ng], perm[ng:]
    from mirna_hallmark.learned import confounders as CF
    Ch = CF.build_C("gtex", list(Yh.index))                                  # ← the healthy C (was MISSING)
    Ct_df = pd.DataFrame(Cm, index=Yt.index)
    M_Tsub = _fit_M(Yt.iloc[fit_i], Xtz.iloc[fit_i], alpha=alpha, C=Ct_df.iloc[fit_i])   # in-state, n=ng
    M_H = _fit_M(Yh, Xhz, alpha=alpha, C=Ch)                                 # out-of-state, n=ng, NOW C-adjusted
    Xte = Xtn[test_i]
    rho_Tsub = spearmanr(Xte @ M_Tsub.reindex(Xtz.columns).fillna(0).to_numpy(), ytr[test_i]).correlation
    rho_H = spearmanr(Xte @ M_H.reindex(Xtz.columns).fillna(0).to_numpy(), ytr[test_i]).correlation
    M_T = _fit_M(Yt, Xtz, alpha=alpha, C=Ct_df)                              # full-n tumour M (for Δ inspection)
    rho_T = spearmanr(Xtn @ M_T.reindex(Xtz.columns).fillna(0).to_numpy(), ytr).correlation
    delta = (M_T - M_H).reindex(M_T.index.union(M_H.index)).fillna(0.0)
    top_delta = delta.reindex(delta.abs().sort_values(ascending=False).index).head(4)
    # retention = healthy coupling as a FRACTION of the n-matched in-state control (both n=ng, held-out).
    # ~1 ⇒ healthy weights transfer ⇒ Δ≈0 (only abundance is state-dependent); <0.75 ⇒ M is state-dependent.
    retention = float(rho_H / rho_Tsub) if (rho_Tsub and rho_Tsub == rho_Tsub) else np.nan
    return {"gene": gene, "n_tum": len(Yt), "n_gtex": ng, "n_pred": Xtz.shape[1],
            "rho_T_full": round(float(rho_T), 3), "rho_Tsub_ho": round(float(rho_Tsub), 3),
            "rho_H_ho": round(float(rho_H), 3), "retention": round(retention, 2),
            "state_dependent": bool(rho_Tsub < 0 and retention < 0.75),
            "top_delta": ", ".join(f"{m}={v:+.2f}" for m, v in top_delta.items())}


def run(genes=None, *, alpha: float = 0.01, family: bool = True) -> pd.DataFrame:
    genes = genes or ["PTEN", "GATA3", "ESR1", "ZEB1", "CDKN1A"]
    rows = [cross_state_transfer(g, alpha=alpha, family=family) for g in genes]
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 175, "display.max_colwidth", 58):
        print(df.to_string(index=False))
    if "retention" in df:
        ok = df.dropna(subset=["retention"])
        print(f"\nΔ=0 test (n-matched, held-out): in-state control ρ(M_Tsub)={ok['rho_Tsub_ho'].mean():.3f} "
              f"vs healthy ρ(M_H)={ok['rho_H_ho'].mean():.3f} | mean retention={ok['retention'].mean():.2f} "
              f"| state-dependent on {int(ok['state_dependent'].sum())}/{len(ok)} genes "
              f"(retention<0.75 ⇒ M state-dependent beyond n-noise; ~1 ⇒ only abundance x is)")
    return df


if __name__ == "__main__":
    run()
