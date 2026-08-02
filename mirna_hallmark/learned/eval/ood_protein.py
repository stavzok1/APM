"""Bar 5 — orthogonal protein layer (CPTAC), the anti-circularity crux (Design §5, §Decision B).

M is fit on TCGA **mRNA**. Here we score it, unrefit, on **CPTAC** at the **protein** layer it never saw
— and, most sharply, on **protein-vs-mRNA discordance** (`protein_resid` = protein beyond what its own
mRNA explains), the translational residual the mRNA fit *structurally cannot see* (Design §1b).

**TWO CPTAC cohorts, different independence — do not conflate:**
- ``prospective`` (pancan122): **independent patients**, self-contained CPTAC-2 miRNA + TMT protein →
  the true non-circular OOD (out-of-patient AND out-of-layer). **The Bar-5 claim rests on this.**
- ``tcga105``: the **same 105 patients** as the TCGA fit, patient-paired by participant id →
  **patient-circular**; only the molecular layer differs. Reported as a labeled *layer-transfer* check
  (does M carry from mRNA to protein on the same people?), NOT as independent validation.
  **Doubly weaker:** its proteome is older **iTRAQ-4** MS (shallower depth, ratio compression, lower
  quantitative accuracy) vs prospective's **TMT-11** — so the prospective cohort is superior on *both*
  counts (independent patients AND higher-resolution protein). Prospective is the primary; tcga105 is a
  caveated cross-check only.

    ρ_cptac(g,s) = Σ_{m∈R(g)} X_cptac(m,s)·M(m,g)     # TCGA-fit M applied to CPTAC arms
    rho_protein  = Spearman(ρ_cptac, protein_z[g])     # expect < 0
    rho_discord  = Spearman(ρ_cptac, protein_resid[g]) # expect < 0 (translational residual)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import regression as LR

_CACHE: dict = {}
CIRCULAR = {"prospective": False, "tcga105": True}  # tcga105 shares patients with the TCGA fit


def _cov(cohort: str, samples):
    """⚠ CELL-COMPOSITION-ADJUSTED covariate block (MH-107). The scores below were previously RAW Spearman —
    NO covariates at all — while the TCGA model being scored already conditions on the 8 Wu-major lineages.
    Bulk protein is heavily compartment-loaded (EMT proteins ↔ CAF fraction r=+0.509; ZEB1 protein +0.768) and
    epithelial miRNAs (miR-200) anti-correlate with stroma, which MANUFACTURES protein 'coupling'. Measured on
    the gold edges: **27.6% of raw protein-coupled edges SIGN-FLIP under composition and the mean |ρ| falls 39%**
    (re-derived 2026-08-01, MH-172). ⚠ this line previously read "27% … falls 60%" and cited "MH-84" — the
    **27% reproduces (71.9–72.8% keep-rate under every subset), the 60% does NOT (36–41%), and MH-84 has no
    registry row** (ids 81–85 were never written). Home for this number: **MH-172**, not a docstring.
    Returns a design matrix [1 | C] aligned to `samples`, or an intercept-only column if C is unavailable."""
    import numpy as _np
    key = f"cov_{cohort}"
    if key not in _CACHE:
        try:
            from mirna_hallmark.learned import confounders as CF
            _CACHE[key] = CF.build_C("cptac" if cohort == "prospective" else "tcga", list(samples))
        except Exception:
            _CACHE[key] = None
    C = _CACHE[key]
    if C is None:
        return _np.ones((len(samples), 1))
    Cm = C.reindex(list(samples)).apply(pd.to_numeric, errors="coerce")
    Cm = Cm.fillna(Cm.median())
    return _np.column_stack([_np.ones(len(samples)), Cm.to_numpy(float)])


def _resid(v, D):
    v = np.asarray(v, dtype=float)
    m = np.isfinite(v)
    out = np.full(v.shape, np.nan)
    if m.sum() > D.shape[1] + 2:
        Dm = D[m]
        out[m] = v[m] - Dm @ np.linalg.lstsq(Dm, v[m], rcond=None)[0]
    return out


def _cptac(cohort: str = "prospective"):
    if cohort not in _CACHE:
        from mirna_hallmark.eval import cptac_validation as V
        L = V.load_cptac_layers(cohort)
        # prospective has its OWN CPTAC-2 miRNA (independent); tcga105 reuses TCGA miRNA (same patients).
        arms = V.load_prospective_mirna_arms() if cohort == "prospective" else D.load_mirna_arms()
        _CACHE[cohort] = {"arms": arms, "protein": L["protein_z"], "resid": L["protein_resid"]}
    return _CACHE[cohort]


def _sp(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    return spearmanr(a[m], b[m]).correlation if m.sum() > 10 else np.nan


def fit_M(gene: str, *, alpha: float = 0.005) -> pd.Series:
    """Per-arm M for a gene, fit on TCGA mRNA (ledger prior + deconvolution composition in C)."""
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    return LR.fit_gene(Y, X, C, w, alpha=alpha)


def score(gene: str, M: pd.Series | None = None, *, cohort: str = "prospective", alpha: float = 0.005) -> dict:
    d = _cptac(cohort)
    arms, protein, resid = d["arms"], d["protein"], d["resid"]
    if gene not in protein.index:
        raise KeyError(f"{gene} not in CPTAC protein")
    if M is None:
        M = fit_M(gene, alpha=alpha)                                # always fit on TCGA mRNA
    regs = [a for a in M.index if a in arms.index and M[a] > 0]
    if not regs:
        raise ValueError(f"no nonzero-weight regulators of {gene} present in CPTAC arms")
    samples = [s for s in arms.columns if s in protein.columns]
    Xc = arms.loc[regs, samples].T.astype(float).fillna(0.0)
    rho = Xc.to_numpy() @ M[regs].to_numpy()                       # learned aggregate per CPTAC sample
    abund = Xc.to_numpy().mean(axis=1)                             # raw-abundance baseline
    prot = protein.loc[gene, samples].astype(float).to_numpy()
    disc = resid.loc[gene, samples].astype(float).to_numpy()
    D = _cov(cohort, samples)                                      # ⚠ MH-107: composition-adjusted (was RAW)
    rho_r, abund_r = _resid(rho, D), _resid(abund, D)
    prot_r, disc_r = _resid(prot, D), _resid(disc, D)
    r_p, r_a = _sp(rho_r, prot_r), _sp(abund_r, prot_r)
    return {
        "gene": gene, "cohort": cohort, "circular": CIRCULAR[cohort],
        "n": len(samples), "n_reg": len(regs), "n_cov": D.shape[1] - 1,
        "rho_protein": r_p, "rho_discord": _sp(rho_r, disc_r),
        "rho_protein_abund": r_a,
        "model_beats_abund": bool(np.nan_to_num(r_p, nan=0) < np.nan_to_num(r_a, nan=0)),
        # provenance: the historical, CONFOUNDED values (no covariates) — for reproducing MH-83/84
        "rho_protein_RAW": _sp(rho, prot), "rho_discord_RAW": _sp(rho, disc),
    }


HUB = ["PTEN", "ESR1", "ZEB1", "BCL2", "PDCD4", "CDKN1B", "GATA3"]


def run(genes=None, *, cohorts=("prospective", "tcga105"), alpha: float = 0.005) -> pd.DataFrame:
    genes = genes or HUB
    rows = []
    for cohort in cohorts:
        Ms = {}
        for g in genes:
            try:
                if g not in Ms:
                    Ms[g] = fit_M(g, alpha=alpha)                  # fit once per gene (TCGA), reuse per cohort
                rows.append(score(g, Ms[g], cohort=cohort))
            except Exception as e:
                rows.append({"gene": g, "cohort": cohort, "error": repr(e)[:60]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 170):
        print(df.round(3).to_string(index=False))
    for cohort in cohorts:
        ok = df[(df["cohort"] == cohort) & df.get("rho_protein").notna()] if "rho_protein" in df else df.iloc[0:0]
        if len(ok):
            tag = "CIRCULAR (same patients — layer transfer only)" if CIRCULAR[cohort] else "INDEPENDENT (non-circular OOD)"
            print(f"\n[{cohort}] {tag}: protein-coupled {int((ok['rho_protein']<0).sum())}/{len(ok)} | "
                  f"discordance-coupled {int((ok['rho_discord']<0).sum())}/{len(ok)} | "
                  f"beats abundance {int(ok['model_beats_abund'].sum())}/{len(ok)}")
    return df


if __name__ == "__main__":
    run()
