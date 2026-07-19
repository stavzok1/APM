"""State/subtype-stratified assembly + cross-state coupling (Design §Decision H; NAT/healthy/subtype plan).

Fit the learned model in **NAT** (tumour-adjacent normal, TCGA sample-type 11), a **PAM50 subtype**, or the
**NAT-paired tumour** subset — beyond the default tumour cohort. Type-filters at the BARCODE level because
`load_mirna_expression` collapses all sample types to 12-char participant (averaging tumour+NAT).

COMPOSITION CAVEAT (analysis guardrail): NAT is normal breast (adipose/stroma/epithelial mix) vs tumour
(epithelial) — cross-state coupling differences are confounded by composition, and NAT has no CIBERSORTx
deconvolution. The clean design is the **paired within-patient** difference (Phase D), which removes the
patient baseline. State-level coupling here is reported as descriptive, composition-uncontrolled.

    X_state, Y_state = state_matrices("11")           # NAT arm×participant, gene×participant
    df = cross_state_coupling(genes)                  # per-gene NAT vs tumour single-aggregate coupling
"""
from __future__ import annotations

from functools import lru_cache

from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import regression as LR


def _stype(bc: str):
    p = str(bc).split("-")
    return p[3][:2] if len(p) >= 4 and len(p[3]) >= 2 else None


@lru_cache(maxsize=1)
def _raw_mirna() -> pd.DataFrame:
    """arm × barcode log2(RPM+1) — all sample types (no participant collapse yet)."""
    from analysis.expression.mirna_target_integration import load_mimat_to_arm
    m2a = load_mimat_to_arm(C.MIRNA_MATURE_LOCI)
    df = pd.read_csv(C.MIRNA_EXPRESSION, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.rename(index=m2a).groupby(level=0).mean()


@lru_cache(maxsize=1)
def _raw_rna() -> pd.DataFrame:
    """gene × barcode log2(TPM+1) — all sample types."""
    df = pd.read_csv(C.RNA_EXPRESSION, sep="\t", index_col=0)
    return df.drop(columns=[c for c in ("Ensembl_ID",) if c in df.columns])


@lru_cache(maxsize=4)
def state_matrices(sample_type: str = "11"):
    """(arm × participant, gene × participant) for a TCGA sample type (11=NAT, 01=tumour)."""
    mi, rn = _raw_mirna(), _raw_rna()
    mic = [c for c in mi.columns if _stype(c) == sample_type]
    rnc = [c for c in rn.columns if _stype(c) == sample_type]
    X = mi[mic].copy(); X.columns = [c[:12] for c in mic]; X = X.T.groupby(level=0).mean().T
    Y = rn[rnc].copy(); Y.columns = [c[:12] for c in rnc]; Y = Y.T.groupby(level=0).mean().T
    X = X[~X.index.duplicated()]; Y = Y[~Y.index.duplicated()]   # dedupe arm/gene symbols
    return X, Y


def _he_arms(gene: str) -> list:
    he = D.high_evidence_edges()
    return sorted(set(he.loc[he["gene"] == gene, "miRNA"]))


def _resid(v, Cmat):
    return v - LinearRegression().fit(Cmat, v).predict(Cmat)


@lru_cache(maxsize=1)
def _deconv_raw():
    """The CIBERSORTx fractions file, read ONCE. Was re-read from NFS on every _cibersortx_state_cov call
    (i.e. once per gene in a batch) — the dominant per-gene I/O cost that thrashed the genome-wide card
    build under parallel workers (memory `batch-nfs-per-unit-reads`). Cache it; slice per call is cheap."""
    from mirna_hallmark.learned import data as LD
    return pd.read_csv(LD._DECONV_PATH, sep="\t").drop_duplicates("Mixture").set_index("Mixture")


def _cibersortx_state_cov(parts, sample_type: str):
    """CIBERSORTx non-malignant compartment fractions for a state — the file HAS tumour (12-char keys)
    AND NAT (`TCGA-..-....-NAT` keys); 111 paired (user, 2026-07-05). Gold-standard composition, consistent
    with the tumour C. Returns None for GTEx (not deconvolved — metagene fallback)."""
    if sample_type not in ("01", "11"):
        return None
    f = _deconv_raw()
    cols = [c for c in LD._DECONV_COLS if c in f.columns]
    key = {p: (f"{p}-NAT" if sample_type == "11" else p) for p in parts}    # NAT rows are '<participant>-NAT'
    idx = {p: k for p, k in key.items() if k in f.index}
    if len(idx) < 25:
        return None
    cov = f.loc[list(idx.values()), cols].apply(pd.to_numeric, errors="coerce")
    cov.index = list(idx.keys())                                            # back to participant keys
    return cov.reindex(parts)


def _state_metagene_cov(rna_state: pd.DataFrame) -> pd.DataFrame:
    """State-comparable composition+proliferation covariates from mRNA (precursor `_state_covariates`
    source='metagene'): proliferation + epithelial/immune/stroma marker metagenes — computable in EVERY
    state (NAT/GTEx have no CPE/CIBERSORTx). CIBERSORTx-on-NAT/GTEx is the documented future upgrade."""
    from mirna_hallmark.analyses.cross_state.cross_state_coupling import (
        EPI_MARKERS, IMMUNE_MARKERS, STROMA_MARKERS, _metagene, _prolif_metagene)
    from mirna_hallmark.hallmark_sets import HallmarkSets
    hs = HallmarkSets.load()
    return pd.DataFrame({
        "prolif": _prolif_metagene(rna_state, hs),
        "epi": _metagene(rna_state, EPI_MARKERS),
        "immune": _metagene(rna_state, IMMUNE_MARKERS),
        "stroma": _metagene(rna_state, STROMA_MARKERS),
    })


def assemble_state(gene: str, sample_type: str = "11", *, participants=None, family: bool = True,
                   metagene_cov: bool = True):
    """(Y, X, C, w) for `gene` in a state: gene mRNA + HE-arm abundance + state-comparable metagene
    composition/proliferation C (NAT/GTEx have no CPE/CIBERSORTx). Returns None if gene/regulators absent."""
    X_s, Y_s = state_matrices(sample_type)
    if gene not in Y_s.index:
        return None
    arms = [a for a in _he_arms(gene) if a in X_s.index]
    if len(arms) < 1:
        return None
    parts = list(Y_s.columns.intersection(X_s.columns))
    if participants is not None:
        parts = [p for p in parts if p in set(participants)]
    if len(parts) < 25:
        return None
    Y = Y_s.loc[gene, parts]
    X = X_s.loc[arms, parts].T                                    # participant × arm
    w = pd.Series(1.0, index=arms)
    if family:
        X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(pd.Index(arms)))
    if metagene_cov:
        C = _cibersortx_state_cov(parts, sample_type)          # gold-standard CIBERSORTx (tumour + NAT)
        if C is None:                                          # GTEx / uncovered → metagene fallback
            C = _state_metagene_cov(Y_s).reindex(parts)
        C = C.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))
    else:
        C = pd.DataFrame(1.0, index=parts, columns=["intercept"])
    return Y, X, C, w


_XSCALE_CACHE: dict = {}
_TCGA_XSCALE = 0.947          # median sd(X_fam) in tumour — the scale the historical 0.1 floor was tuned on


def x_scale(sample_type: str) -> float:
    """Median family-dose sd in a state — the scale any dose floor must be relative to.
    Measured 2026-07-12: **tumour 0.947 · NAT 0.682 · GTEx 0.235 · (CPTAC 0.861)**."""
    if sample_type in _XSCALE_CACHE:
        return _XSCALE_CACHE[sample_type]
    from mirna_hallmark.learned import structural_identity as SI
    from mirna_hallmark.learned.evidence import ledger as LG
    ed = LG.pooled_he_edges()
    sds = []
    for g in list(dict.fromkeys(ed["gene"]))[:80]:
        try:
            arms = pd.Index(sorted(SI._functional_regulators(g)))
            if arms.empty:
                continue
            d = _state_family_data(g, sample_type, FAM.family_of(arms))
            if d is None:
                continue
            sds.extend(d[1].std(ddof=0).to_numpy())
        except Exception:
            continue
    sds = [s for s in sds if np.isfinite(s) and s > 0]
    _XSCALE_CACHE[sample_type] = float(np.median(sds)) if sds else _TCGA_XSCALE
    return _XSCALE_CACHE[sample_type]


def x_floor(sample_type: Optional[str]) -> float:
    """The dose floor, TRANSLATED to the state's own scale.

    ⚠ **The historical floor was an ABSOLUTE `sd < 0.1`, and that is a SCALE BUG across states** (found via the
    CPTAC session's report of the same bug in `gauge.beta_table`, 2026-07-12). Median sd(X_fam) is **tumour
    0.947 · NAT 0.682 · GTEx 0.235** — GTEx's miRNA dose has ~4× less spread — so an absolute 0.1 zeroed
    **32% of GTEx's families vs 3% of tumour's**, silently deleting a third of the healthy regulators. Measured
    consequence: the GTEx-side M REORDERS (ρ(M_abs, M_free) = **ESR1 0.43** · PTEN 0.75 · CDK6 0.84; ZEB1 1.00),
    which reaches the published wiring verdicts.

    The floor is therefore anchored to the scale it was tuned on: `0.1 · x_scale(state) / 0.947`.
    ⇒ **tumour keeps EXACTLY 0.1 (zero ripple — `subtype.py` and every tumour result are unchanged)**; GTEx gets
    ~0.025 and NAT ~0.072, i.e. the SAME relative stringency."""
    if sample_type is None or sample_type == "01":
        return 0.1
    return float(0.1 * x_scale(sample_type) / _TCGA_XSCALE)


def _bagged_nnls(Y, Xf, Cm, *, n_boot: int = 40, seed: int = 0,
                 sample_type: Optional[str] = None) -> pd.Series:
    """THE canonical per-edge weight estimator: **bagged** non-negative least squares on C-residualised,
    **z-scored** family predictors. Fixed support (no lasso selection) kills cross-state selection
    instability; family collapse removes within-family collinearity; z-scoring conditions the NNLS
    (level-scale NNLS is ill-conditioned → recovers noise arms); BAGGING (mean over `n_boot` resamples)
    drives the estimate variance down ~1/B. Verified reproducible (corr ~0.99) and recovers curated drivers
    (miR-21/103/181 on PTEN) — replaces the single-fit lasso wherever COEFFICIENTS are read (attribution)."""
    from scipy.optimize import nnls
    yr = -_resid(Y.to_numpy(float), Cm)
    sd = Xf.std(ddof=0)
    Xz = ((Xf - Xf.mean()) / (sd + 1e-9)).fillna(0.0)
    Xz.loc[:, sd < x_floor(sample_type)] = 0.0                  # SCALE-FREE dose floor (see `x_floor`); tumour
    Xz = Xz.to_numpy(float)                                     # no information; zero it (else z-scoring by ~0 sd
    n = len(yr); rng = np.random.default_rng(seed); acc = np.zeros(Xz.shape[1])  # amplifies noise → spurious M)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        c, _ = nnls(Xz[idx], yr[idx]); acc += c
    return pd.Series(acc / n_boot, index=Xf.columns)


def _state_family_data(gene, sample_type, fam_map):
    """(Y, family-collapsed X, C) for a gene in a state — arm-level assembled, then family-collapsed onto the
    fixed `fam_map` support. sample_type: '01' tumour (deconv/cell-intrinsic) · '11' NAT · 'gtex' healthy."""
    if sample_type == "01":
        Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
        arms = list(X.columns)
    elif sample_type == "11":
        a = assemble_state(gene, "11", family=False)
        if a is None:
            return None
        Y, X, C, _ = a
        mn = Y.notna() & C.notna().all(axis=1)
        Y, X, C = Y[mn], X.loc[mn].fillna(0.0), C.loc[mn]
        arms = list(X.columns)
    else:  # GTEx healthy
        from mirna_hallmark.learned import state as STA
        Xg, Yg = STA._gtex_mirna(), STA._gtex_mrna()
        if gene not in Yg.index:
            return None
        arms = [a for a in fam_map.index if a in Xg.index]
        Cg = _state_metagene_cov(Yg)
        p = [c for c in Yg.columns if c in Xg.columns and c in Cg.index]
        if len(p) < 25 or not arms:
            return None
        X = Xg.loc[arms, p].T.fillna(0.0); Y = Yg.loc[gene, p]
        C = Cg.loc[p].apply(lambda s: s.fillna(s.median()))
        mk = Y.notna(); Y, X, C = Y[mk], X.loc[mk], C.loc[mk]
    arms = [a for a in arms if a in fam_map.index]
    Xf, _, _ = FAM.collapse_by_family(X[arms], pd.Series(1.0, index=arms), fam_map.reindex(arms))
    return Y, Xf, C


def canonical_M(gene: str, sample_type: str = "01", *, arm_level: bool = True, n_boot: int = 40) -> pd.Series:
    """Canonical per-edge weight for `gene` in a state: bagged z-scored NNLS (`_bagged_nnls`) on the fixed
    HE seed-**family** support (Design §F: family→gene is the identified estimand). `arm_level=True`
    broadcasts each family's weight to its member arms — the per-arm split is a NOMINATION (abundance
    apportions the realized contribution), not a resolved coefficient. Use everywhere M is READ for
    attribution (budget, ΔM, decompose, wiring, card) so the whole attribution layer is estimator-consistent."""
    he = D.high_evidence_edges()
    arms = pd.Index(sorted(set(he.loc[he["gene"] == gene, "miRNA"])))
    if arms.empty:
        return pd.Series(dtype=float, name="M")
    fam_map = FAM.family_of(arms)
    d = _state_family_data(gene, sample_type, fam_map)
    if d is None:
        return pd.Series(dtype=float, name="M")
    Y, Xf, C = d
    Mf = _bagged_nnls(Y, Xf, C.to_numpy(float), n_boot=n_boot, sample_type=sample_type)
    if not arm_level:
        return Mf.rename("M")
    return pd.Series({a: float(Mf.get(fam_map.get(a, a), 0.0)) for a in arms}, name="M")


def _aggregate_coupling(gene: str, sample_type: str, participants=None, *, alpha=0.01,
                        metagene_cov: bool = True) -> dict:
    a = assemble_state(gene, sample_type, participants=participants, metagene_cov=metagene_cov)
    if a is None:
        return {"rho": np.nan, "n": 0}
    Y, X, C, w = a
    Cm = C.to_numpy(float)
    Xz = ((X - X.mean()) / (X.std(ddof=0) + 1e-9)).fillna(0.0)
    M = LR.fit_gene(Y, Xz, C, w, alpha=alpha)
    agg = Xz.to_numpy(float) @ M.reindex(Xz.columns).fillna(0).to_numpy()
    return {"rho": spearmanr(_resid(agg, Cm), _resid(Y.to_numpy(float), Cm)).correlation, "n": len(Y),
            "top": ", ".join(f"{m}={v:.2f}" for m, v in M[M > 0].sort_values(ascending=False).head(3).items())}


def cross_state_coupling(genes, *, alpha: float = 0.01) -> pd.DataFrame:
    """Per-gene learned-aggregate coupling in NAT (11) vs tumour (01) — same HE regulators, z-scored,
    intercept-only C. Descriptive (composition-uncontrolled in NAT)."""
    rows = []
    for g in genes:
        nat = _aggregate_coupling(g, "11", alpha=alpha)
        tum = _aggregate_coupling(g, "01", alpha=alpha)
        rows.append({"gene": g, "n_nat": nat["n"], "n_tum": tum["n"],
                     "rho_nat": round(nat["rho"], 3) if nat["rho"] == nat["rho"] else np.nan,
                     "rho_tum": round(tum["rho"], 3) if tum["rho"] == tum["rho"] else np.nan,
                     "top_tum": tum.get("top", "")})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 160):
        print(df.to_string(index=False))
    ok = df.dropna(subset=["rho_nat", "rho_tum"])
    if len(ok):
        print(f"\nmean ρ: NAT={ok['rho_nat'].mean():.3f}  tumour={ok['rho_tum'].mean():.3f} "
              f"| coupling present in NAT (ρ<−0.15) on {int((ok['rho_nat']<-0.15).sum())}/{len(ok)} genes")
    return df


@lru_cache(maxsize=1)
def paired_delta_matrices():
    """Δ = tumour − NAT for the paired patients: (arm × patient Δmirna, gene × patient ΔmRNA, patients).
    Within-patient differencing removes the patient baseline (composition, purity, batch) — the clean
    cross-state design (Phase D)."""
    Xt, Yt = state_matrices("01")
    Xn, Yn = state_matrices("11")
    pts = sorted(set(Xt.columns) & set(Yt.columns) & set(Xn.columns) & set(Yn.columns))
    dX = Xt.reindex(columns=pts) - Xn.reindex(columns=pts)
    dY = Yt.reindex(columns=pts) - Yn.reindex(columns=pts)
    return dX, dY, pts


def realization(genes, *, alpha: float = 0.005) -> pd.DataFrame:
    """Phase D — is the NAT→tumour pressure shift REALIZED within patients? Fit M on the full tumour cohort
    (independent of the pairs), then per paired patient predict Δpressure = Σ_m M(m,g)·Δmirna(m) and correlate
    with the actual target shift ΔmRNA(g). Expect ρ<0 (miRNA gained → target lost) — the paired causal test."""
    dX, dY, pts = paired_delta_matrices()
    rows = []
    for g in genes:
        try:
            Yt, Xt, Ct, w = LD.assemble_gene(g, w_prior_source="ledger")
            M = LR.fit_gene(Yt, Xt, Ct, w, alpha=alpha)
        except Exception:
            continue
        regs = [a for a in M[M > 0].index if a in dX.index]
        if g not in dY.index or not regs:
            continue
        pred = dX.loc[regs, pts].fillna(0.0).T.to_numpy(float) @ M[regs].to_numpy()  # predicted Δpressure per patient
        dy = dY.loc[g, pts].to_numpy(float)
        msk = np.isfinite(pred) & np.isfinite(dy)
        if msk.sum() < 25:
            continue
        pred, dy_m = pred[msk], dy[msk]
        # decompose the target shift: also the mean NAT→tumour Δ of the top regulator (abundance shift)
        top = M[M > 0].sort_values(ascending=False).index[0]
        rows.append({"gene": g, "n_pairs": int(msk.sum()), "n_reg": len(regs),
                     "rho_realized": round(float(spearmanr(pred, dy_m).correlation), 3),
                     "top_reg": top, "mean_dReg": round(float(dX.loc[top, pts].mean()), 2),
                     "mean_dTarget": round(float(dy_m.mean()), 2)})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 160):
        print(df.to_string(index=False))
    ok = df.dropna(subset=["rho_realized"])
    if len(ok):
        print(f"\nrealization: mean ρ(predicted Δpressure, Δtarget) = {ok['rho_realized'].mean():.3f} "
              f"| realized (ρ<−0.1) on {int((ok['rho_realized']<-0.1).sum())}/{len(ok)} genes "
              f"(within-patient, paired — baseline-free)")
    return df


def _pam50() -> pd.Series:
    from analysis.utils.common.loaders import load_clinical
    clin = load_clinical(C.CLINICAL_UNIFIED).drop_duplicates("participant")
    return clin.set_index("participant")["PAM50_final"]


def shift_vs_weight(genes, *, alpha: float = 0.005, deconv: bool = True) -> pd.DataFrame:
    """(#1) Per edge: the WEIGHT M(m,g) vs the NAT→tumour ABUNDANCE shift Δx̄(m), and the realized
    contribution M·Δx̄ — does the model's weight, or the abundance shift, drive the target change? deconv=True
    fits M on composition-adjusted tumour (so stroma edges don't dominate)."""
    dX, dY, pts = paired_delta_matrices()
    rows = []
    for g in genes:
        try:
            Yt, Xt, Ct, w = LD.assemble_gene(g, w_prior_source="ledger", deconv=deconv)
            M = LR.fit_gene(Yt, Xt, Ct, w, alpha=alpha)
        except Exception:
            continue
        for m in M[M > 0].sort_values(ascending=False).index:
            if m not in dX.index:
                continue
            dbar = float(dX.loc[m, pts].mean())
            rows.append({"gene": g, "arm": m, "weight_M": round(float(M[m]), 3),
                         "mean_dAbund": round(dbar, 3), "contribution": round(float(M[m]) * dbar, 3)})
    df = pd.DataFrame(rows)
    if len(df):
        print("weight vs NAT→tumour abundance shift (top contributors):")
        print(df.reindex(df["contribution"].abs().sort_values(ascending=False).index).head(14).to_string(index=False))
        print(f"\nSpearman(|weight|, |Δabund|) = {spearmanr(df['weight_M'].abs(), df['mean_dAbund'].abs()).correlation:+.3f} "
              "(are high-weight edges the ones that shift? ~0 ⇒ weight and shift are independent axes)")
    return df


def realization_by_subtype(genes, *, alpha: float = 0.005, deconv: bool = True) -> pd.DataFrame:
    """(#2) Phase-D realization stratified by PAM50 subtype (of the tumour). M is fit on the WHOLE tumour
    cohort; the Δ and realization are computed per subtype → is the NAT→tumour progression subtype-specific?
    Powered mainly for LumA (n~55); others reported with n."""
    dX, dY, pts = paired_delta_matrices()
    pam = _pam50()
    subt = pd.Series({p: pam.get(p) for p in pts})
    rows = []
    for g in genes:
        try:
            Yt, Xt, Ct, w = LD.assemble_gene(g, w_prior_source="ledger", deconv=deconv)
            M = LR.fit_gene(Yt, Xt, Ct, w, alpha=alpha)
        except Exception:
            continue
        regs = [a for a in M[M > 0].index if a in dX.index]
        if g not in dY.index or not regs:
            continue
        pred_all = dX.loc[regs, pts].fillna(0.0).T.to_numpy(float) @ M[regs].to_numpy()
        dy_all = dY.loc[g, pts].to_numpy(float)
        row = {"gene": g}
        for st in ["LumA", "LumB", "Basal", "Her2"]:
            idx = [i for i, p in enumerate(pts) if subt.get(p) == st]
            if len(idx) < 10:
                row[st] = np.nan; continue
            pr, dy = pred_all[idx], dy_all[idx]
            msk = np.isfinite(pr) & np.isfinite(dy)
            row[st] = round(float(spearmanr(pr[msk], dy[msk]).correlation), 2) if msk.sum() >= 8 else np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    if len(df):
        print("realization ρ by subtype (M = whole-cohort tumour; Δ + realization per subtype):")
        print(df.to_string(index=False))
        ns = subt.value_counts()
        print(f"\npaired-n per subtype: " + ", ".join(f"{s}={int(ns.get(s,0))}" for s in ['LumA','LumB','Basal','Her2']))
    return df


def budget_shift(gene: str, *, alpha: float = 0.005, deconv: bool = True) -> pd.DataFrame:
    """(E7/G4 across states) — the WITHIN-GENE contribution ranking and how it SHIFTS NAT→tumour.
    contribution(m,state) = the arm's share of its FAMILY's realized pressure M(fam)·X_fam(state), apportioned
    by the arm's linear-RPM share of the family pool (so Σ_arms(family) = M(fam)·X_fam, NO family-size inflation:
    the family weight is broadcast to arms, and summing M·log-abundance would over-count multi-arm families by
    ~their size — see LEARNED_MODEL_ESTIMATOR_MAP §4). share = contribution / Σ; rank = within-gene rank.
    M is the CELL-INTRINSIC (deconv) weight so shares aren't stroma shadows."""
    M = canonical_M(gene, "01", arm_level=True)                          # CANONICAL bagged family weight (deconv/cell-intrinsic)
    regs = list(M[M > 0].index)
    if not regs:
        return pd.DataFrame()
    fam = FAM.family_of(pd.Index(regs)).reindex(regs)                    # arm → family (for size-correct pooling)
    Mv = M.reindex(regs).to_numpy(float)
    # state mean abundances. GTEx (healthy) uses within-GTEx rank/share only — NO QN needed: a share is a
    # within-state ratio (scale-free) and rank is ordinal, so the cross-platform TPM/RPM offset cancels.
    # MEDIAN, not mean — robust to rare-spike arms (a few very-high samples pull the mean up and overstate an arm
    # typically absent), and consistent with the card's regime/grank axes (MH-166 follow-up).
    levels = {"HLY": None, "NAT": state_matrices("11")[0].reindex(regs).median(axis=1),
              "TUM": state_matrices("01")[0].reindex(regs).median(axis=1)}
    try:
        from mirna_hallmark.learned import state as ST
        levels["HLY"] = ST._gtex_mirna().reindex(regs).median(axis=1)
    except Exception:
        levels.pop("HLY")
    df = pd.DataFrame({"arm": regs, "M": np.round(Mv, 3)})
    for st, lv in levels.items():
        a = lv.reindex(regs).fillna(0.0).to_numpy(float)                 # log2(abundance+1) per arm, this state
        df[f"abund_{st}"] = np.round(a, 2)
        a_lin = np.clip(2.0 ** a - 1.0, 0.0, None)                       # → linear RPM/TPM
        pool = pd.Series(a_lin, index=regs).groupby(fam).transform("sum").to_numpy(float)  # Σ within family
        x_fam = np.log2(1.0 + pool)                                      # pooled family abundance = the fitted predictor
        share_in_fam = np.divide(a_lin, pool, out=np.zeros_like(a_lin), where=pool > 0)    # arm's share of its family
        c = np.clip(Mv * x_fam * share_in_fam, 0.0, None)               # Σ_arms(family) = M_fam·X_fam (no inflation)
        tot = c.sum()
        df[f"share_{st}"] = np.round(c / tot, 3) if tot > 0 else 0.0
        df[f"rank_{st}"] = pd.Series(c).rank(ascending=False).astype(int).values
    df = df.sort_values("share_TUM", ascending=False)
    states_present = list(levels)
    df["d_rank_NAT_TUM"] = df["rank_NAT"] - df["rank_TUM"]                # +ve = rose NAT→tumour
    if "HLY" in states_present:
        df["d_rank_HLY_TUM"] = df["rank_HLY"] - df["rank_TUM"]            # +ve = rose healthy→tumour (acquired, dHT)
    order = ["arm", "M"] + [f"{p}_{st}" for st in states_present for p in ("share", "rank")] \
        + [c for c in ("d_rank_HLY_TUM", "d_rank_NAT_TUM") if c in df]
    print(f"\n=== {gene}: within-gene contribution budget across states {states_present} (cell-intrinsic M) ===")
    print(df[order].to_string(index=False))
    if "d_rank_HLY_TUM" in df:
        mv = df.reindex(df["d_rank_HLY_TUM"].abs().sort_values(ascending=False).index).iloc[0]
        print(f"  biggest healthy→tumour budget mover: {mv['arm']} rank {mv['rank_HLY']}→{mv['rank_TUM']} "
              f"(share {mv['share_HLY']:.2f}→{mv['share_TUM']:.2f})")
    return df


if __name__ == "__main__":
    CELL_INTRINSIC = ["PTEN", "ZEB1", "ESR1", "GATA3", "BCL2", "CDKN1A", "CDKN1B", "RB1"]
    print("=== #1: weight M vs NAT→tumour abundance shift (deconv-M, cell-intrinsic) ===")
    shift_vs_weight(CELL_INTRINSIC)
    print("\n=== #2: realization by PAM50 subtype (deconv-M) ===")
    realization_by_subtype(CELL_INTRINSIC)
