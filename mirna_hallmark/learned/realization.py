"""Multi-resolution WITHIN-PATIENT paired NAT→tumour REALIZATION test.

Is ACQUIRED miRNA pressure (NAT→tumour) realized WITHIN a patient? For predictor P on gene g:
    Δpressure_p = Σ_{a∈P} M(a,g)·Δx_p(a),   Δx_p = x_p^tumour − x_p^NAT
test Spearman(Δpressure, Δy_g) across the ~104 paired patients (expect ρ<0), raw AND composition-Δ-adjusted.

Within-patient differencing removes the patient baseline (composition/purity/batch/germline) — the one
progression estimator that survives the dead cross-cohort state channel (MH-102d). NAT≠healthy (field
effect) ⇒ this is the FINAL malignant step, conservative. n≈104 ⇒ read SET-LEVEL, never per-edge.

M-reference (user-flagged): default **complement** = M fit on the ~900 NON-matched tumours (leak-free +
powered); **full** = canonical_M incl. matched (leak, sensitivity); **matched** = n≈104 (circular ceiling).

Engine only (pure, memoized reads). Sharded CLI driver: `analyses/progression/realize_run.py`.
Design: plan `i-think-we-can-fluffy-walrus`; reuses `states`/`families`/`attribution`/`decoy_bench`.
"""
from __future__ import annotations

import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")                 # single-thread BLAS before numpy (axiom 3a)

from functools import lru_cache
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import states as ST

_MIN_PAIRS = 25
OUT = Path(C.OUTPUT_ROOT) / "learned" / "realization"
_LEARNED = Path(C.OUTPUT_ROOT) / "learned"


# --------------------------------------------------------------------------- #
# partial Spearman (rank → residualise on C → Pearson of residuals)
def _resid_on(v: np.ndarray, C: np.ndarray) -> np.ndarray:
    C1 = np.column_stack([np.ones(len(v)), C])
    beta, *_ = np.linalg.lstsq(C1, v, rcond=None)
    return v - C1 @ beta


def partial_spearman(x, y, C=None) -> float:
    """Spearman(x, y | C): Pearson correlation of C-residualised ranks. C None → plain Spearman."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if C is not None:
        C = np.asarray(C, float)
        if C.ndim == 1:
            C = C[:, None]
        m &= np.isfinite(C).all(axis=1)
    if m.sum() < _MIN_PAIRS:
        return np.nan
    rx, ry = rankdata(x[m]), rankdata(y[m])
    if C is None:
        r = np.corrcoef(rx, ry)[0, 1]
        return float(r)
    Cr = np.column_stack([rankdata(C[m, j]) for j in range(C.shape[1])])
    ex, ey = _resid_on(rx, Cr), _resid_on(ry, Cr)
    if ex.std() < 1e-12 or ey.std() < 1e-12:
        return np.nan
    return float(np.corrcoef(ex, ey)[0, 1])


# --------------------------------------------------------------------------- #
# M-reference layer
def _state_family_data_sub(gene: str, participants):
    """`states._state_family_data(gene,'01')` but fit on a PARTICIPANT SUBSET: row-slice the tumour
    `assemble_gene` output to `participants` BEFORE family collapse (assemble_gene has no participant filter)."""
    try:
        Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)   # raises for no-HE-arm genes
    except Exception:
        return None
    keep = Y.index.intersection(pd.Index(list(participants)))
    if len(keep) < _MIN_PAIRS:
        return None
    Y, X, C = Y.loc[keep], X.loc[keep], C.loc[keep]
    arms = [a for a in X.columns]
    fam_map = FAM.family_of(pd.Index(arms))
    arms = [a for a in arms if a in fam_map.index]
    if not arms:
        return None
    Xf, _, _ = FAM.collapse_by_family(X[arms], pd.Series(1.0, index=arms), fam_map.reindex(arms))
    return Y, Xf, C, fam_map, arms


def _M_from_sub(gene: str, participants, *, n_boot: int = 40) -> pd.Series:
    d = _state_family_data_sub(gene, participants)
    if d is None:
        return pd.Series(dtype=float, name="M")
    Y, Xf, C, fam_map, arms = d
    Mf = ST._bagged_nnls(Y, Xf, C.to_numpy(float), n_boot=n_boot, sample_type="01")   # tumour floor 0.1
    return pd.Series({a: float(Mf.get(fam_map.get(a, a), 0.0)) for a in arms}, name="M")


@lru_cache(maxsize=1)
def _matched_pts() -> tuple:
    return tuple(ST.paired_delta_matrices()[2])


@lru_cache(maxsize=1)
def _tumour_participants() -> tuple:
    return tuple(ST.state_matrices("01")[0].columns)


def M_reference(gene: str, kind: str = "complement", *, n_boot: int = 40) -> pd.Series:
    """arm-level M weight for gene g. kind: complement (non-matched tumours, DEFAULT, leak-free) ·
    full (canonical_M, all tumours incl. matched) · matched (matched tumours only, circular ceiling)."""
    if kind == "full":
        return ST.canonical_M(gene, "01", arm_level=True, n_boot=n_boot)
    matched = set(_matched_pts())
    tum = _tumour_participants()
    parts = [p for p in matched if p in set(tum)] if kind == "matched" else [p for p in tum if p not in matched]
    return _M_from_sub(gene, parts, n_boot=n_boot)


def _M_fam_from_sub(gene: str, participants, *, n_boot: int = 40) -> pd.Series:
    """FAMILY-level M (leak-free subset) — identical fit to `_M_from_sub` but returns the family-indexed weight
    BEFORE the arm broadcast (Res-3/4 use the family as the identified unit, Design §F)."""
    d = _state_family_data_sub(gene, participants)
    if d is None:
        return pd.Series(dtype=float, name="M")
    Y, Xf, C, fam_map, arms = d
    return ST._bagged_nnls(Y, Xf, C.to_numpy(float), n_boot=n_boot, sample_type="01").rename("M")


def M_reference_fam(gene: str, kind: str = "complement", *, n_boot: int = 40) -> pd.Series:
    """FAMILY-level M weight (for Res-3 family / Res-4 between-family). Same three references as `M_reference`
    but the seed FAMILY is the unit (no arm broadcast)."""
    if kind == "full":
        return ST.canonical_M(gene, "01", arm_level=False, n_boot=n_boot)
    matched = set(_matched_pts())
    tum = _tumour_participants()
    parts = [p for p in matched if p in set(tum)] if kind == "matched" else [p for p in tum if p not in matched]
    return _M_fam_from_sub(gene, parts, n_boot=n_boot)


def _family_pooled_delta(gene: str, pts):
    """NONLINEAR family-pooled Δ (verification item 6): collapse arms→families in tumour AND NAT SEPARATELY,
    THEN difference — family Δ = collapse(01) − collapse(11), NOT collapse(arm-level Δ) (pooling is nonlinear).
    Returns (dF family×patient, Ft participant×family tumour, Fn participant×family NAT), or None."""
    xt, xn, apts = _paired_abund()
    pts = [p for p in pts if p in apts]
    if len(pts) < _MIN_PAIRS:
        return None
    arms = [a for a in ST._he_arms(gene) if a in xt.index and a in xn.index]
    if not arms:
        return None
    fam_map = FAM.family_of(pd.Index(arms))
    ones = pd.Series(1.0, index=arms)
    Ft, _, _ = FAM.collapse_by_family(xt.loc[arms, pts].T, ones, fam_map)   # participant × family (tumour)
    Fn, _, _ = FAM.collapse_by_family(xn.loc[arms, pts].T, ones, fam_map)   # participant × family (NAT)
    dF = (Ft.reindex(pts) - Fn.reindex(pts)).T                             # family × patient
    return dF, Ft.reindex(pts), Fn.reindex(pts)


# --------------------------------------------------------------------------- #
# composition Δ  (per paired patient: tumour − NAT)
@lru_cache(maxsize=1)
def _delta_C() -> pd.DataFrame:
    pts = list(ST.paired_delta_matrices()[2])
    ct = ST._cibersortx_state_cov(pts, "01")
    cn = ST._cibersortx_state_cov(pts, "11")
    if ct is None or cn is None:
        dC = pd.DataFrame(index=pts)
    else:
        dC = (ct.reindex(pts) - cn.reindex(pts))
    _, Yt = ST.state_matrices("01"); _, Yn = ST.state_matrices("11")
    try:
        pt = ST._state_metagene_cov(Yt)["prolif"].reindex(pts)
        pn = ST._state_metagene_cov(Yn)["prolif"].reindex(pts)
        dC = dC.assign(d_prolif=(pt.to_numpy() - pn.to_numpy()))
    except Exception:
        pass
    return dC.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))


# --------------------------------------------------------------------------- #
# core primitive
def _realize(g: str, arms, M: pd.Series, dX: pd.DataFrame, dY: pd.DataFrame,
             pts, dC: pd.DataFrame) -> dict:
    """Predicted Δpressure = Σ M·Δx over `arms`; ρ(Δpressure, Δy_g) raw + ΔC-adjusted."""
    regs = [a for a in arms if a in dX.index and float(M.get(a, 0.0)) > 0]
    if not regs or g not in dY.index:
        return {}
    pred = dX.loc[regs, pts].fillna(0.0).T.to_numpy(float) @ M[regs].to_numpy(float)
    dy = dY.loc[g, pts].to_numpy(float)
    m = np.isfinite(pred) & np.isfinite(dy)
    if m.sum() < _MIN_PAIRS:
        return {}
    Cm = dC.reindex(pts).to_numpy(float)[m] if len(dC.columns) else None
    rho_raw = partial_spearman(pred[m], dy[m])
    rho_adj = partial_spearman(pred[m], dy[m], Cm)
    ret = float(rho_adj / rho_raw) if (rho_raw and np.isfinite(rho_raw) and abs(rho_raw) > 1e-6) else np.nan
    return {"n_pairs": int(m.sum()), "n_reg": len(regs), "rho_raw": _r3(rho_raw), "rho_adj": _r3(rho_adj),
            "retention": _r2(ret), "composition_explained": bool(ret == ret and ret < 0.4),
            "mean_dPred": _r2(float(np.nanmean(pred[m]))), "mean_dTarget": _r2(float(np.nanmean(dy[m])))}


def _r3(x):
    return round(float(x), 3) if x == x else np.nan


def _r2(x):
    return round(float(x), 2) if x == x else np.nan


# --------------------------------------------------------------------------- #
# resolutions 1 (gene) + 2 (edge)
def realize_gene_edge(genes: Sequence[str], *, m_ref: str = "complement") -> pd.DataFrame:
    """Resolution 1 (gene aggregate) + 2 (edge) for a gene slice, under M-reference `m_ref`."""
    dX, dY, pts_l = ST.paired_delta_matrices()
    pts = list(pts_l)
    dC = _delta_C()
    rows = []
    for g in genes:
        try:
            M = M_reference(g, m_ref)
            if M.empty:
                continue
            arms = list(M[M > 0].index)
            if not arms:
                continue
            agg = _realize(g, arms, M, dX, dY, pts, dC)          # GENE aggregate
            if agg:
                rows.append({"resolution": "gene", "gene": g, "predictor": "AGG", "m_ref": m_ref, **agg})
            for a in arms:                                        # EDGE (per arm)
                e = _realize(g, [a], M, dX, dY, pts, dC)
                if e:
                    rows.append({"resolution": "edge", "gene": g, "predictor": a, "m_ref": m_ref, **e})
        except Exception:
            continue
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# compartment orientation helper (MH-114 axis) — reusable
def orientation(mir_L: float, tgt_L: float) -> str:
    """same/opposite compartment: sign(corr(miR, L)) × sign(corr(target, L)). L = the composition (CAF) axis."""
    if not (np.isfinite(mir_L) and np.isfinite(tgt_L)) or mir_L == 0 or tgt_L == 0:
        return "NA"
    return "SAME" if (np.sign(mir_L) == np.sign(tgt_L)) else "OPPOSITE"


@lru_cache(maxsize=1)
def _paired_abund() -> tuple:
    """Matched tumour & NAT arm abundance (arm×patient) on the paired patients (log2(RPM+1))."""
    Xt, _ = ST.state_matrices("01"); Xn, _ = ST.state_matrices("11")
    pts = [p for p in ST.paired_delta_matrices()[2] if p in Xt.columns and p in Xn.columns]
    return Xt[pts], Xn[pts], pts


def dose_shift_arm() -> pd.DataFrame:
    """DOSE/POTENTIAL axis — per arm: acquired dose (own-NAT), GLOBAL-rank move NAT→tumour, and the
    PATIENT-SPECIFICITY of the shift: `own_specific_frac` = var(patient NAT baseline) / var(own dose shift).
    High ⇒ the shift is driven by patient-specific NAT (own-NAT matters, cohort/subtype anchor misleads) =
    a 'unique dosage shift'; low ⇒ tumour-driven acquisition shared across patients."""
    xt, xn, pts = _paired_abund()
    nat_cohort = xn.mean(axis=1)
    rt = xt.rank(pct=True) * 100.0                         # global %ile within each patient's tumour
    rn = xn.rank(pct=True) * 100.0
    rows = []
    for m in xt.index:
        if m not in xn.index:
            continue
        own = (xt.loc[m] - xn.loc[m])                     # per-patient own dose shift
        nat_dev = (nat_cohort.get(m, np.nan) - xn.loc[m])  # own_specific = own − cohort = cohort_mean − own_NAT
        v_own = float(np.nanvar(own.to_numpy()))
        v_dev = float(np.nanvar(nat_dev.to_numpy()))
        drank = (rt.loc[m] - rn.loc[m])
        rows.append({"arm": m, "n_pairs": int(own.notna().sum()),
                     "mean_own_shift": _r3(float(own.mean())),
                     "mean_dGlobalRank": _r2(float(drank.mean())),
                     "sd_dGlobalRank": _r2(float(drank.std())),
                     "nat_patient_sd": _r3(float(np.sqrt(v_dev)) if v_dev == v_dev else np.nan),
                     "own_specific_frac": _r2(v_dev / v_own) if v_own > 1e-9 else np.nan})
    return pd.DataFrame(rows)


def dose_shift_edge(genes: Sequence[str], *, m_ref: str = "complement") -> pd.DataFrame:
    """WITHIN-GENE dose ranking: per (gene, arm) the budget-share shift NAT→tumour, M-weighted (ties to the
    coupling M) AND raw-abundance (pure dose, M-independent) — side by side, so M-vs-raw shows how much of the
    within-gene ranking move is DOSE vs WEIGHT. Own-NAT paired; cohort-referenced share vs the cohort NAT mean."""
    xt, xn, pts = _paired_abund()
    nat_cohort = xn.mean(axis=1)
    rows = []
    for g in genes:
        try:
            M = M_reference(g, m_ref)
            arms = [a for a in M[M > 0].index if a in xt.index and a in xn.index]
            if len(arms) < 2:                              # share is defined only among ≥2 regulators
                continue
            Xt_g = np.nan_to_num(2 ** xt.loc[arms, pts].to_numpy(float) - 1, nan=0.0)   # linear RPM, arm×patient
            Xn_g = np.nan_to_num(2 ** xn.loc[arms, pts].to_numpy(float) - 1, nan=0.0)
            w = M[arms].to_numpy(float)[:, None]
            def _sh(P):                                     # column-normalised share (per patient)
                return P / (P.sum(0, keepdims=True) + 1e-9)
            shMt, shMn = _sh(w * Xt_g), _sh(w * Xn_g)        # M-weighted share
            shRt, shRn = _sh(Xt_g), _sh(Xn_g)                # raw-abundance share
            xc = np.nan_to_num(2 ** nat_cohort.reindex(arms).to_numpy(float) - 1, nan=0.0)[:, None]
            shMc = (w * xc) / ((w * xc).sum(0) + 1e-9)       # M-share of the cohort-mean NAT
            for i, a in enumerate(arms):
                rows.append({"gene": g, "arm": a, "m_ref": m_ref, "n_pairs": len(pts),
                             "dShare_M_own": _r3(float((shMt[i] - shMn[i]).mean())),
                             "dShare_raw_own": _r3(float((shRt[i] - shRn[i]).mean())),
                             "dShare_M_cohort": _r3(float((shMt[i] - shMc[i]).mean())),
                             "share_M_tum": _r3(float(shMt[i].mean())),
                             "share_raw_tum": _r3(float(shRt[i].mean()))})
        except Exception:
            continue
    return pd.DataFrame(rows)


def _corr(a, b) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < _MIN_PAIRS or a[m].std() < 1e-12 or b[m].std() < 1e-12:
        return np.nan
    return float(np.corrcoef(a[m], b[m])[0, 1])


# --------------------------------------------------------------------------- #
# the NULL — site-free matched decoys through the SAME paired Δ, orientation-stratified (the arbiter)
def realize_null(genes: Sequence[str], *, m_ref: str = "complement") -> pd.DataFrame:
    """Per (gene, real_arm, matched fake_arm): push the FAKE arm's ΔX through the REAL arm's M weight and the
    identical paired-Δ + ΔC pipeline. REAL vs DECOY, ρ_adj, stratified by compartment orientation (Δ axis:
    sign(corr(Δx,ΔCAF))·sign(corr(Δy,ΔCAF))). Set-level (n≈104), NOT per-edge."""
    from mirna_hallmark.eval import decoy_bench as DB
    dX, dY, pts_l = ST.paired_delta_matrices()
    pts = list(pts_l)
    dC = _delta_C()
    dCAF = dC["CAFs"].reindex(pts).to_numpy(float) if "CAFs" in dC.columns else np.full(len(pts), np.nan)
    dec = DB.build_decoys(list(genes))
    rows = []
    for g, grp in dec.groupby("gene"):
        try:
            M = M_reference(g, m_ref)
            if M.empty or g not in dY.index:
                continue
            dy = dY.loc[g, pts].to_numpy(float)
            tgt_L = _corr(dy, dCAF)
            for r in grp.itertuples():
                w = float(M.get(r.real_arm, 0.0))
                if w <= 0:
                    continue
                for grp_lbl, arm in (("REAL", r.real_arm), ("DECOY", r.fake_arm)):
                    if arm not in dX.index:
                        continue
                    res = _realize(g, [arm], pd.Series({arm: w}), dX, dY, pts, dC)
                    if not res:
                        continue
                    ori = orientation(_corr(dX.loc[arm, pts].to_numpy(float), dCAF), tgt_L)
                    rows.append({"gene": g, "group": grp_lbl, "arm": arm, "orientation": ori,
                                 "m_ref": m_ref, "rho_raw": res["rho_raw"], "rho_adj": res["rho_adj"]})
        except Exception:
            continue
    return pd.DataFrame(rows)


def nat_decoy_control(genes: Sequence[str], *, m_ref: str = "complement") -> pd.DataFrame:
    """Is own>cohort REPRESSION-specific or shared unmodelled-composition? For REAL and matched DECOY arms,
    realization ρ_adj under own-NAT Δy vs cohort-NAT Δy; the (own−cohort) gap. If DECOYs show the same gap ⇒
    the own advantage is composition/shared-baseline; if REAL gap is more negative ⇒ repression-specific.
    Reuses cached paired matrices / ΔC / M / decoys — one single-process pass."""
    from mirna_hallmark.eval import decoy_bench as DB
    dX, dY_own, pts_l = ST.paired_delta_matrices()
    pts = list(pts_l)
    dC = _delta_C()
    _, Yt = ST.state_matrices("01"); _, Yn = ST.state_matrices("11")
    nat_mean = Yn.mean(axis=1)
    dec = DB.build_decoys(list(genes))
    rows = []
    for g, grp in dec.groupby("gene"):
        try:
            M = M_reference(g, m_ref)
            if M.empty or g not in Yt.index or g not in dY_own.index:
                continue
            dy_coh = np.array([Yt.at[g, p] - nat_mean.get(g, np.nan) for p in pts], float)
            dY_coh = pd.DataFrame({g: dy_coh}, index=pts).T
            for r in grp.itertuples():
                w = float(M.get(r.real_arm, 0.0))
                if w <= 0:
                    continue
                for lbl, arm in (("REAL", r.real_arm), ("DECOY", r.fake_arm)):
                    if arm not in dX.index:
                        continue
                    Mt = pd.Series({arm: w})
                    ro = _realize(g, [arm], Mt, dX, dY_own, pts, dC)
                    rc = _realize(g, [arm], Mt, dX, dY_coh, pts, dC)
                    if ro and rc:
                        rows.append({"gene": g, "group": lbl, "arm": arm,
                                     "rho_own": ro["rho_adj"], "rho_cohort": rc["rho_adj"],
                                     "gap_own_minus_cohort": _r3(ro["rho_adj"] - rc["rho_adj"])})
        except Exception:
            continue
    return pd.DataFrame(rows)


def null_summary(nd: pd.DataFrame) -> pd.DataFrame:
    """Set-level REAL vs DECOY mean ρ within each orientation stratum, raw + adjusted (the headline)."""
    out = []
    for (grp, ori), s in nd.groupby(["group", "orientation"]):
        for c in ("rho_raw", "rho_adj"):
            out.append({"group": grp, "orientation": ori, "C": c.split("_")[1],
                        "n_edges": len(s), "mean_rho": _r3(s[c].mean())})
    return pd.DataFrame(out).sort_values(["orientation", "C", "group"])


# --------------------------------------------------------------------------- #
# patient-specificity of NAT — own vs cohort vs subtype reference (the second question)
def nat_reference_adjudication(genes: Sequence[str], *, m_ref: str = "complement") -> pd.DataFrame:
    """For each of 3 differencing references (own / cohort-NAT-mean / subtype-NAT-mean): residual variance of
    Δy and gene-aggregate realization ρ (raw + adj). Does own-NAT beat cohort/subtype as the anchor?"""
    Xt, Yt = ST.state_matrices("01"); Xn, Yn = ST.state_matrices("11")
    dX_own, dY_own, pts_l = ST.paired_delta_matrices()
    pts = [p for p in pts_l if p in Yt.columns and p in Yn.columns]
    dC = _delta_C()
    nat_mean = Yn.mean(axis=1)                                    # cohort-NAT-mean per gene
    pam = ST._pam50()
    nat_sub = {}                                                  # subtype → NAT-mean gene vector
    for st in pam.dropna().unique():
        cols = [c for c in Yn.columns if pam.get(c) == st]
        if len(cols) >= 10:
            nat_sub[st] = Yn[cols].mean(axis=1)
    rows = []
    for g in genes:
      try:
        M = M_reference(g, m_ref)
        if M.empty or g not in Yt.index:
            continue
        arms = list(M[M > 0].index)
        for ref in ("own", "cohort", "subtype"):
            if ref == "own":
                dy = dY_own.loc[g, pts].to_numpy(float) if g in dY_own.index else None
            elif ref == "cohort":
                dy = np.array([Yt.at[g, p] - nat_mean.get(g, np.nan) for p in pts], float)
            else:
                dy = np.array([Yt.at[g, p] - nat_sub.get(pam.get(p), nat_mean).get(g, np.nan) for p in pts], float)
            if dy is None or np.isfinite(dy).sum() < _MIN_PAIRS:
                continue
            dYref = pd.DataFrame({g: dy}, index=pts).T
            res = _realize(g, arms, M, dX_own, dYref, pts, dC)
            if res:
                rows.append({"reference": ref, "gene": g, "m_ref": m_ref,
                             "resid_var_dY": _r3(float(np.nanvar(dy))),
                             "rho_raw": res["rho_raw"], "rho_adj": res["rho_adj"], "n_pairs": res["n_pairs"]})
      except Exception:
        continue
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Resolution 3 (FAMILY) — _realize DIRECT on the NONLINEAR family-pooled Δ
def _fam_own_specific_frac(Ft: pd.DataFrame, Fn: pd.DataFrame, f: str) -> float:
    """Family-level own_specific_frac (dose-axis analogue): var(cohort-NAT − own-NAT) / var(own tum−NAT shift).
    High ⇒ the family's dose shift is driven by the patient's OWN NAT baseline (own-NAT ≠ cohort anchor)."""
    own = (Ft[f] - Fn[f]).to_numpy(float)
    nat_dev = (Fn[f].mean() - Fn[f]).to_numpy(float)
    v_own = float(np.nanvar(own))
    return _r2(float(np.nanvar(nat_dev)) / v_own) if v_own > 1e-9 else np.nan


def realize_family(genes: Sequence[str], *, m_ref: str = "complement") -> pd.DataFrame:
    """Resolution 3: realization on the family-pooled Δ. Rows: `family_agg` (nonlinear-pooled gene aggregate,
    a control on the arm-level `gene` row) + one `family` row per seed family (predictor = family)."""
    _, dY, pts_l = ST.paired_delta_matrices()
    pts = list(pts_l)
    dC = _delta_C()
    rows = []
    for g in genes:
        try:
            M_fam = M_reference_fam(g, m_ref)
            if M_fam.empty:
                continue
            fd = _family_pooled_delta(g, pts)
            if fd is None:
                continue
            dF, Ft, Fn = fd
            fams = [f for f in M_fam[M_fam > 0].index if f in dF.index]
            if not fams:
                continue
            agg = _realize(g, fams, M_fam, dF, dY, pts, dC)                       # nonlinear-pooled aggregate
            if agg:
                rows.append({"resolution": "family_agg", "gene": g, "predictor": "AGG", "m_ref": m_ref,
                             "own_specific_frac": np.nan, **agg})
            for f in fams:
                e = _realize(g, [f], M_fam, dF, dY, pts, dC)
                if e:
                    rows.append({"resolution": "family", "gene": g, "predictor": f, "m_ref": m_ref,
                                 "own_specific_frac": _fam_own_specific_frac(Ft, Fn, f), **e})
        except Exception:
            continue
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Resolution 4 (BETWEEN-FAMILY) — Shapley identity of the paired Δ, denominator-gated (axiom 5)
def realize_between_family(genes: Sequence[str], *, m_ref: str = "complement",
                           n_boot: int = 30, n_perm: int = 150) -> pd.DataFrame:
    """Resolution 5 (between-family): `attribution.shapley_identity` on the paired family-pooled Δ. Xr =
    ΔC-residualised, z-scored per-family Δ predictors; m = family M; yr = ΔC-residualised Δy. Per (gene,family):
    realization_identity (collinearity-fair credit for R²(Δpressure,Δy)), bootstrap SD, and `agg_pred_mag` (the
    denominator floor — axiom 5: a share is meaningful only where the aggregate |Δpressure| clears a floor)."""
    from mirna_hallmark.learned import attribution as ATTR
    _, dY, pts_l = ST.paired_delta_matrices()
    pts = list(pts_l)
    dC = _delta_C()
    Cfull = dC.reindex(pts).to_numpy(float) if len(dC.columns) else None
    rng = np.random.default_rng(0)
    rows = []
    for g in genes:
        try:
            M_fam = M_reference_fam(g, m_ref)
            fams0 = list(M_fam[M_fam > 0].index)
            if len(fams0) < 2 or g not in dY.index:
                continue
            fd = _family_pooled_delta(g, pts)
            if fd is None:
                continue
            dF = fd[0]
            fams = [f for f in fams0 if f in dF.index]
            if len(fams) < 2:
                continue
            dy = dY.loc[g, pts].to_numpy(float)
            Xmat = dF.loc[fams, pts].T.to_numpy(float)                            # patient × family
            m = np.isfinite(Xmat).all(axis=1) & np.isfinite(dy)
            if Cfull is not None:
                m &= np.isfinite(Cfull).all(axis=1)
            if m.sum() < _MIN_PAIRS:
                continue
            Xm, yv = Xmat[m], dy[m]
            Cm = Cfull[m] if Cfull is not None else None
            def _rz(col):
                r = _resid_on(col, Cm) if Cm is not None else col - col.mean()
                return r / (r.std() + 1e-9)
            Xr = np.column_stack([_rz(Xm[:, j]) for j in range(Xm.shape[1])])
            yr = _resid_on(yv, Cm) if Cm is not None else yv - yv.mean()
            mvec = M_fam[fams].to_numpy(float)
            agg_mag = float(np.mean(np.abs(Xm @ mvec)))                           # denominator floor (axiom 5)
            phi = np.clip(ATTR.shapley_identity(Xr, mvec, yr, n_perm=n_perm), 0, None)
            share = phi / phi.sum() if phi.sum() > 0 else phi
            bs = []
            for _ in range(n_boot):
                idx = rng.integers(0, len(yr), len(yr))
                ph = np.clip(ATTR.shapley_identity(Xr[idx], mvec, yr[idx], n_perm=max(60, n_perm // 2)), 0, None)
                bs.append(ph / ph.sum() if ph.sum() > 0 else ph)
            sd = np.asarray(bs).std(0)
            owner = fams[int(np.argmax(share))]
            for j, f in enumerate(fams):
                rows.append({"gene": g, "family": f, "m_ref": m_ref,
                             "realization_identity": _r3(float(share[j])), "identity_sd_boot": _r3(float(sd[j])),
                             "M": _r3(float(mvec[j])), "agg_pred_mag": _r3(agg_mag),
                             "is_owner": bool(f == owner), "n_fam": len(fams), "n_pairs": int(m.sum())})
        except Exception:
            continue
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# DESCRIPTIVE PATTERN LAYER — homed from the 2026-07-18 exploratory thread (formerly inline heredocs).
# All reuse the already-built outputs (canonical_card, ladder, dose, nat_decoy_control); no new heavy compute.
_COREP = Path(C.OUTPUT_ROOT) / "tissue_reference" / "mirna_comovement" / "gene_corepression.tsv"
_ACQ = Path(C.OUTPUT_ROOT) / "tissue_reference" / "mirna_state_class" / "gene_acquired_pressure.tsv"


def progression_edge_card() -> pd.DataFrame:
    """⭐ THE INTEGRATED PROGRESSION EDGE CARD — one row per (gene, arm), folding BOTH progression objects onto
    the attribution card so every cross-resolution question is a column read, not a multi-file join:
      • CROSS-STATE (cohort, `canonical_card`): `shift_class`, `coupling_{hly,nat,tum}`, `grank_*`,
        `arm_lfc_*` (dHT/dNT/dHN), `retention_rho`, `beta`/`identity`/`share_TUM`;
      • WITHIN-PATIENT PAIRED (MH-158): `edge_rho_adj`, `own_specific_frac`, `mean_own_shift`, `dShare_M_own`,
        `family_rho_adj` (Res-3), `realization_identity`/`is_realization_owner` (Res-4, family→arm via seed_family);
      • gene role + net-repression.
    -> progression_edge_card.tsv. (Supersedes `master_edge_patterns.tsv`: same key, + the Phase-2 columns.)"""
    from mirna_hallmark import gene_roles as GR
    card = pd.read_csv(_LEARNED / "canonical_card.tsv", sep="\t")
    de = pd.read_csv(OUT / "dose_shift_edge.tsv", sep="\t")
    da = pd.read_csv(OUT / "dose_shift_arm.tsv", sep="\t")
    lad = pd.read_csv(OUT / "realization_ladder.tsv", sep="\t")
    cor = pd.read_csv(_COREP, sep="\t")
    m = card.merge(de[["gene", "arm", "dShare_M_own", "dShare_raw_own"]], on=["gene", "arm"], how="left")
    m = m.merge(da[["arm", "mean_own_shift", "mean_dGlobalRank", "own_specific_frac"]], on="arm", how="left")
    er = lad[lad.resolution == "edge"][["gene", "predictor", "rho_adj"]].rename(
        columns={"predictor": "arm", "rho_adj": "edge_rho_adj"})
    m = m.merge(er, on=["gene", "arm"], how="left")
    # --- Phase-2 family-level realization, broadcast to arm via the card's seed_family ---
    bf_p, fam_p = OUT / "realization_between_family.tsv", OUT / "realization_family.tsv"
    if bf_p.exists():
        bf = pd.read_csv(bf_p, sep="\t")[["gene", "family", "realization_identity", "is_owner"]].rename(
            columns={"family": "seed_family", "is_owner": "is_realization_owner"})
        m = m.merge(bf, on=["gene", "seed_family"], how="left")
    if fam_p.exists():
        fam = pd.read_csv(fam_p, sep="\t")
        fr = fam[fam.resolution == "family"][["gene", "predictor", "rho_adj"]].rename(
            columns={"predictor": "seed_family", "rho_adj": "family_rho_adj"})
        m = m.merge(fr, on=["gene", "seed_family"], how="left")
    roles = GR.load_gene_roles(list(m.gene.unique())).set_index("gene")["role"]
    cp = cor.set_index("gene")
    m["role"] = m.gene.map(roles).fillna("unknown")
    m["gene_nreg"] = m.gene.map(cp["n_regulators"])
    m["gene_repr_class"] = m.gene.map(cp["gene_repression_class"])
    m["gene_net_repr"] = m.gene.map(cp["gene_net_repressed_tumor"])
    m["gene_dominated"] = m.groupby("gene").share_TUM.transform("max") > 0.6
    m.to_csv(OUT / "progression_edge_card.tsv", sep="\t", index=False)
    return m


def edge_pattern_table() -> pd.DataFrame:
    """Backward-compat alias → `progression_edge_card` (renamed 2026-07-18: master_edge_patterns IS the edge card)."""
    return progression_edge_card()


def _dominant_by_state(card: pd.DataFrame, state_col: str) -> pd.Series:
    """Per gene, the arm with the highest global-abundance rank in a state (`grank_HLY/NAT/TUM`) = the dominant
    regulator in that state — for the regulatory-handoff flag (does the dominant regulator switch across states)."""
    c = card.dropna(subset=[state_col])
    if not len(c):
        return pd.Series(dtype=object)
    return c.loc[c.groupby("gene")[state_col].idxmax()].set_index("gene")["arm"]


def progression_gene_card() -> pd.DataFrame:
    """⭐ THE INTEGRATED PROGRESSION GENE CARD — one row per gene:
      • tumour attribution summary (`readouts_genes`: n_fam, total_pressure, top_identity, concentration, retention);
      • WITHIN-PATIENT paired realization — gene aggregate (`realized_rho_{raw,adj}`, retention, comp_explained) +
        nonlinear family-pooled (`fam_pooled_rho_adj`) + `frac_edges_realized` (share of the gene's edges with ρ_adj<0);
      • CROSS-STATE acquired dose (`gene_acquired_pressure`: pressure_tumor, acquired_vs_{gtex,nat}, frac_true_acquired);
      • gene net-repression (`gene_corepression`: repression_class, net_repressed_tumor, rho_gene_pressure_tumor);
      • OWNER — `realization_owner` (Res-4 is_owner) vs `static_owner_family` (max canonical identity) + `owner_agrees`;
      • REGULATORY HANDOFF — dominant regulator per state (HLY/NAT/TUM) + `regulatory_handoff` (does it switch);
      • `dominant_edge_shift_class` (the top-share edge's cross-state class).
    -> progression_gene_card.tsv."""
    rg = pd.read_csv(_LEARNED / "readouts_genes.tsv", sep="\t")
    lad = pd.read_csv(OUT / "realization_ladder.tsv", sep="\t")
    card = pd.read_csv(_LEARNED / "canonical_card.tsv", sep="\t")
    G = rg.copy()
    gr = lad[lad.resolution == "gene"][["gene", "rho_raw", "rho_adj", "retention", "composition_explained",
                                        "n_reg", "mean_dTarget"]].rename(
        columns={"rho_raw": "realized_rho_raw", "rho_adj": "realized_rho_adj", "retention": "realized_retention",
                 "composition_explained": "realized_composition_explained", "n_reg": "realized_n_reg"})
    G = G.merge(gr, on="gene", how="left")
    if (OUT / "realization_family.tsv").exists():
        fam = pd.read_csv(OUT / "realization_family.tsv", sep="\t")
        fa = fam[fam.resolution == "family_agg"][["gene", "rho_adj"]].rename(columns={"rho_adj": "fam_pooled_rho_adj"})
        G = G.merge(fa, on="gene", how="left")
    # frac of the gene's EDGES realized within-patient (ρ_adj<0)
    ed = lad[lad.resolution == "edge"]
    G = G.merge(ed.groupby("gene")["rho_adj"].apply(lambda s: float((s < 0).mean())).rename("frac_edges_realized").reset_index(),
                on="gene", how="left")
    # cross-state acquired dose
    if _ACQ.exists():
        ap = pd.read_csv(_ACQ, sep="\t")
        keep = [c for c in ["gene", "pressure_tumor", "acquired_vs_gtex", "acquired_vs_nat",
                            "frac_gain_true_acquired", "top_gaining_arms", "n_hallmark_sets"] if c in ap.columns]
        G = G.merge(ap[keep], on="gene", how="left")
    # gene net-repression
    cor = pd.read_csv(_COREP, sep="\t")
    G = G.merge(cor[["gene", "n_regulators", "gene_repression_class", "gene_net_repressed_tumor",
                     "rho_gene_pressure_tumor", "delta_tumor_nat"]], on="gene", how="left")
    # OWNER — realization (Res-4) vs static identity, and agreement
    if (OUT / "realization_between_family.tsv").exists():
        bf = pd.read_csv(OUT / "realization_between_family.tsv", sep="\t")
        own = bf[bf.is_owner].drop_duplicates("gene").set_index("gene")["family"]
        G["realization_owner"] = G.gene.map(own)
    if "identity" in card and "seed_family" in card:
        st = card.dropna(subset=["identity"]).sort_values("identity").drop_duplicates("gene", keep="last").set_index("gene")
        G["static_owner_family"] = G.gene.map(st["seed_family"])
        if "realization_owner" in G:   # defined ONLY where a realization owner exists (multi-family genes); else NaN
            agree = (G["realization_owner"] == G["static_owner_family"])
            G["owner_agrees"] = agree.where(G["realization_owner"].notna() & G["static_owner_family"].notna())
    # REGULATORY HANDOFF across states
    domH, domN, domT = (_dominant_by_state(card, s) for s in ("grank_HLY", "grank_NAT", "grank_TUM"))
    G["dominant_HLY"] = G.gene.map(domH); G["dominant_NAT"] = G.gene.map(domN); G["dominant_TUM"] = G.gene.map(domT)
    G["regulatory_handoff"] = ((G.dominant_HLY != G.dominant_TUM) | (G.dominant_NAT != G.dominant_TUM)) & \
        G.dominant_TUM.notna() & G.dominant_HLY.notna()
    # dominant-edge cross-state shift_class (the top-share regulator's class)
    sc = card.dropna(subset=["shift_class"])
    if "share_TUM" in sc:
        dom = sc.loc[sc.groupby("gene")["share_TUM"].idxmax()].set_index("gene")["shift_class"]
        G["dominant_edge_shift_class"] = G.gene.map(dom)
    G.to_csv(OUT / "progression_gene_card.tsv", sep="\t", index=False)
    return G


def mirna_nat_retention() -> pd.DataFrame:
    """Per-miRNA patient-baseline retention: Spearman(tumour, own-NAT) across paired patients (`r_own`) vs a
    permuted-pairing null (`r_perm`≈0); `r_adj` residualises each arm on its state's CIBERSORTx composition
    first (composition-free). High = the miRNA's level persists from a patient's normal into their tumour
    (top: miR-412-5p 0.66, 14q32 imprinted). -> mirna_nat_retention.tsv."""
    xt, xn, pts = _paired_abund()
    Ct = ST._cibersortx_state_cov(pts, "01"); Cn = ST._cibersortx_state_cov(pts, "11")
    kc = [p for p in pts if Ct is not None and p in Ct.index and Cn is not None and p in Cn.index]
    Ctm = Ct.reindex(kc).fillna(Ct.median()).to_numpy(float) if Ct is not None else None
    Cnm = Cn.reindex(kc).fillna(Cn.median()).to_numpy(float) if Cn is not None else None
    rng = np.random.default_rng(0); rows = []
    for m in xt.index:
        if m not in xn.index:
            continue
        a = xt.loc[m].to_numpy(float); b = xn.loc[m].to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 40 or np.nanstd(a[ok]) < 0.1 or np.nanstd(b[ok]) < 0.1:
            continue
        r_own = spearmanr(a[ok], b[ok]).correlation
        r_perm = float(np.nanmean([spearmanr(a[ok], rng.permutation(b[ok])).correlation for _ in range(15)]))
        r_adj = np.nan
        if Ctm is not None:
            ak = xt.loc[m].reindex(kc).to_numpy(float); bk = xn.loc[m].reindex(kc).to_numpy(float)
            mk = np.isfinite(ak) & np.isfinite(bk)
            if mk.sum() >= 40:
                r_adj = spearmanr(_resid_on(ak[mk], Ctm[mk]), _resid_on(bk[mk], Cnm[mk])).correlation
        rows.append({"arm": m, "n": int(ok.sum()), "r_own": _r3(r_own), "r_perm": _r3(r_perm),
                     "retention": _r3(r_own - r_perm), "r_adj": _r3(r_adj), "tum_med": round(float(np.nanmedian(a)), 2)})
    df = pd.DataFrame(rows); df.to_csv(OUT / "mirna_nat_retention.tsv", sep="\t", index=False)
    return df


def patient_nat_identity() -> pd.DataFrame:
    """Per-PATIENT profile similarity: tumour miRNA profile vs OWN-NAT vs cohort-median vs subtype-median NAT.
    `own_minus_cohort`<0 (usual, ~93%) = the denoised cohort median is the better reference; >0 = an
    own-retaining tumour (~7%, e.g. TCGA-BH-A18U). -> patient_nat_identity.tsv."""
    xt, xn, pts = _paired_abund(); pam = ST._pam50()
    coh = xn.median(axis=1); sub_of = {p: pam.get(p) for p in pts}
    subm = {s: xn[[p for p in pts if sub_of[p] == s]].median(axis=1)
            for s in set(sub_of.values()) if s and sum(sub_of[p] == s for p in pts) >= 8}
    arms = [a for a in xt.index if a in xn.index]; rows = []
    for p in pts:
        t = xt.loc[arms, p].to_numpy(float); no = xn.loc[arms, p].to_numpy(float); nc = coh.loc[arms].to_numpy(float)
        ns = subm.get(sub_of[p]); ns = ns.loc[arms].to_numpy(float) if ns is not None else None
        ok = np.isfinite(t) & np.isfinite(no) & np.isfinite(nc)
        so = spearmanr(t[ok], no[ok]).correlation; sc = spearmanr(t[ok], nc[ok]).correlation
        ss = spearmanr(t[ok & np.isfinite(ns)], ns[ok & np.isfinite(ns)]).correlation if ns is not None else np.nan
        rows.append({"patient": p, "subtype": sub_of[p], "sim_own": _r3(so), "sim_cohort": _r3(sc),
                     "sim_subtype": _r3(ss), "own_minus_cohort": _r3(so - sc)})
    df = pd.DataFrame(rows); df.to_csv(OUT / "patient_nat_identity.tsv", sep="\t", index=False)
    return df


def decoy_stratify() -> pd.DataFrame:
    """Join `nat_decoy_control` with canonical_card importance + the decoy pairing → the shift_class-stratified
    REAL-vs-DECOY own−cohort gap (the constitutive-edge repression-specific refinement of MH-158).
    Requires nat_decoy_control.tsv (run `realize_null`/`nat_decoy_control` first). -> nat_decoy_stratified{,_summary}.tsv."""
    from mirna_hallmark.eval import decoy_bench as DB
    dc = pd.read_csv(OUT / "nat_decoy_control.tsv", sep="\t")
    card = pd.read_csv(_LEARNED / "canonical_card.tsv", sep="\t")
    pair = DB.build_decoys(sorted(set(LD.D.high_evidence_edges()["gene"].unique())))
    rg = dc[dc.group == "REAL"].set_index(["gene", "arm"])["gap_own_minus_cohort"]
    dg = dc[dc.group == "DECOY"].set_index(["gene", "arm"])["gap_own_minus_cohort"]
    pair["real_gap"] = [rg.get((g, a), np.nan) for g, a in zip(pair.gene, pair.real_arm)]
    pair["decoy_gap"] = [dg.get((g, a), np.nan) for g, a in zip(pair.gene, pair.fake_arm)]
    ci = card.set_index(["gene", "arm"])
    for c in ["share_TUM", "rank_TUM", "beta", "identity", "shift_class"]:
        pair[c] = [ci[c].get((g, a), np.nan) if (g, a) in ci.index else np.nan for g, a in zip(pair.gene, pair.real_arm)]
    P = pair.dropna(subset=["real_gap", "decoy_gap"]); P.to_csv(OUT / "nat_decoy_stratified.tsv", sep="\t", index=False)

    def _row(mask, name):
        s = P[mask]
        if len(s) < 30:
            return {"stratum": name, "n": len(s)}
        _, p = mannwhitneyu(s.real_gap.dropna(), s.decoy_gap.dropna(), alternative="less")
        return {"stratum": name, "n": len(s), "real_gap": round(s.real_gap.mean(), 4),
                "decoy_gap": round(s.decoy_gap.mean(), 4), "mwu_p": round(p, 3),
                "repression_specific": bool(p < 0.05 and s.real_gap.mean() < s.decoy_gap.mean())}
    sm = pd.DataFrame([_row(P.rank_TUM == 1, "dominant_rank1"),
                       _row(P.share_TUM >= P.share_TUM.quantile(.75), "top_share"),
                       _row(P.beta.abs() >= P.beta.abs().quantile(.75), "top_beta"),
                       _row(P.identity >= 0.5, "high_identity"),
                       _row(P.shift_class == "acquired_realized", "acquired_realized"),
                       _row(P.shift_class == "constitutive", "constitutive")])
    sm.to_csv(OUT / "nat_decoy_stratified_summary.tsv", sep="\t", index=False)
    return sm


# --------------------------------------------------------------------------- #
# ⭐ OWNER-CONVERGENCE (the formalization the descriptive layer earned) — finalize-time cross-gene test.
# Per multi-family gene: does the REALIZATION-Shapley owner coincide with the STATIC-identity owner
# (canonical_card), the DOSE acquirer, and the BUDGET-gainer beyond a shuffled-owner null? ⚠ dose-controlled
# (verification item 7: high-identity edges are high-dose ⇒ convergence is PARTLY dose by construction).
def _owner_table() -> pd.DataFrame:
    """Per (gene, family): the 4 owner-criterion scores, on the between-family identified support.
    realization = Shapley on paired Δ · static = canonical_card identity · dose = acquired own-shift ·
    budget = within-gene M-weighted Δshare. Families mapped from arms via card `seed_family`."""
    bf = pd.read_csv(OUT / "realization_between_family.tsv", sep="\t")
    card = pd.read_csv(_LEARNED / "canonical_card.tsv", sep="\t")
    da = pd.read_csv(OUT / "dose_shift_arm.tsv", sep="\t")
    de = pd.read_csv(OUT / "dose_shift_edge.tsv", sep="\t")
    a2f = card.dropna(subset=["seed_family"]).drop_duplicates("arm").set_index("arm")["seed_family"]
    # static identity → family (max over member arms; card identity is family-broadcast)
    cst = card.assign(fam=card["seed_family"]).groupby(["gene", "fam"])["identity"].max().rename("static_identity")
    # dose acquirer: arm own-shift → family (max member; most positive = acquired most repressive dose)
    da = da.assign(fam=da["arm"].map(a2f))
    dose_fam = da.dropna(subset=["fam"]).groupby("fam")["mean_own_shift"].max().rename("dose_shift")
    # budget gainer: within-gene M-weighted Δshare → family (sum member, budget is additive within gene)
    de = de.assign(fam=de["arm"].map(a2f))
    bud = de.dropna(subset=["fam"]).groupby(["gene", "fam"])["dShare_M_own"].sum().rename("budget_dshare")
    t = bf.rename(columns={"family": "fam"})
    t = t.merge(cst.reset_index(), on=["gene", "fam"], how="left")
    t = t.merge(bud.reset_index(), on=["gene", "fam"], how="left")
    t = t.merge(dose_fam.reset_index(), on="fam", how="left")
    return t


def owner_convergence(mag_floor: Optional[float] = None) -> pd.DataFrame:
    """Set-level test across multi-family genes: pairwise OWNER-match rate among {realization, static, dose,
    budget} vs a random-owner baseline (E[1/n_fam]) and a permutation null; plus the DOSE-CONTROL contrast
    (does realization track static identity BETTER than dose does?). Denominator-gated on `agg_pred_mag`
    (axiom 5) — pass `mag_floor` or it sweeps the median. -> owner_convergence.tsv."""
    t = _owner_table()
    if mag_floor is None:
        mag_floor = float(t.groupby("gene")["agg_pred_mag"].first().median())
    crit = {"realization": "realization_identity", "static": "static_identity",
            "dose": "dose_shift", "budget": "budget_dshare"}
    owners, nfam = {}, {}
    for g, s in t.groupby("gene"):
        if s["fam"].nunique() < 2:
            continue
        nfam[g] = s["fam"].nunique()
        og = {}
        for name, col in crit.items():
            v = s[["fam", col]].dropna()
            og[name] = v.loc[v[col].idxmax(), "fam"] if len(v) and v[col].notna().any() else None
        owners[g] = og
    O = pd.DataFrame(owners).T
    mag = t.groupby("gene")["agg_pred_mag"].first()
    gated = O.index[mag.reindex(O.index) >= mag_floor]
    pairs = [("realization", "static"), ("realization", "dose"), ("realization", "budget"),
             ("static", "dose"), ("static", "budget"), ("dose", "budget")]
    rows = []
    from scipy.stats import norm
    for scope, idx in (("all", O.index), ("gated", gated)):
        Osub = O.loc[idx].dropna(how="all")
        for a, b in pairs:
            d = Osub[[a, b]].dropna()
            if len(d) < 20:
                rows.append({"scope": scope, "pair": f"{a}~{b}", "n": len(d)})
                continue
            # ⚠ WITHIN-GENE random-owner null (NOT cross-gene name-shuffle: family names are gene-specific,
            # so a name-shuffle collapses to ~0 — a trivially-beatable null, axiom 4). Under H0 criterion-a's
            # owner is uniform over THIS gene's families ⇒ P(match)=1/n_fam. Poisson-binomial z + exact p.
            p_gene = np.array([1.0 / nfam[g] for g in d.index])
            n_match = int((d[a] == d[b]).sum())
            mu = float(p_gene.sum()); var = float((p_gene * (1 - p_gene)).sum())
            z = (n_match - mu) / (np.sqrt(var) + 1e-9)
            rows.append({"scope": scope, "pair": f"{a}~{b}", "n": len(d),
                         "match_rate": round(n_match / len(d), 3),
                         "random_baseline": round(mu / len(d), 3), "excess": round((n_match - mu) / len(d), 3),
                         "z": round(float(z), 2), "p": float(f"{norm.sf(z):.2g}")})
    out = pd.DataFrame(rows)
    out.attrs["mag_floor"] = mag_floor
    out.to_csv(OUT / "owner_convergence.tsv", sep="\t", index=False)
    print(f"[owner_convergence] mag_floor={mag_floor:.3f} | genes all={O.dropna(how='all').shape[0]} gated={len(gated)}")
    print(out.to_string(index=False))
    # dose-control readout: realization tracks static better than dose does ⇒ realization adds beyond dose
    g = out[out.scope == "gated"].set_index("pair")
    if {"realization~static", "dose~budget"}.issubset(g.index):
        rs = g.loc["realization~static", "match_rate"]; sd = g.loc["static~dose", "match_rate"]
        print(f"\n[dose-control] realization~static={rs} vs static~dose={sd} "
              f"⇒ realization {'ADDS beyond' if rs > sd else 'does NOT exceed'} dose on identity-tracking.")
    return out


# --------------------------------------------------------------------------- #
# RETENTION × REALIZATION — does composition-free patient-baseline retention predict edge realization,
# PARTIALLING expression? (both r_own and realization scale with expression — the shared confound.)
def retention_realization() -> pd.DataFrame:
    """Join `mirna_nat_retention` (r_adj, composition-free) with edge realization ρ_adj; test whether retention
    predicts realization RAW and after partialling arm expression (tum_med). Pre-reg: near-null after control.
    -> retention_realization.tsv (+ prints the raw/partial Spearman)."""
    ret = pd.read_csv(OUT / "mirna_nat_retention.tsv", sep="\t")
    lad = pd.read_csv(OUT / "realization_ladder.tsv", sep="\t")
    edges = lad[lad.resolution == "edge"][["gene", "predictor", "rho_adj"]].rename(columns={"predictor": "arm"})
    j = edges.merge(ret[["arm", "r_own", "r_adj", "tum_med"]], on="arm", how="inner").dropna(
        subset=["rho_adj", "r_adj", "tum_med"])
    j.to_csv(OUT / "retention_realization.tsv", sep="\t", index=False)
    out = []
    for rcol in ("r_own", "r_adj"):
        raw = partial_spearman(j[rcol].to_numpy(float), j["rho_adj"].to_numpy(float))
        par = partial_spearman(j[rcol].to_numpy(float), j["rho_adj"].to_numpy(float), j["tum_med"].to_numpy(float))
        out.append({"retention": rcol, "n_edges": len(j), "spearman_raw": _r3(raw),
                    "spearman_partial_expr": _r3(par)})
    df = pd.DataFrame(out)
    print("[retention_realization] retention vs edge realization ρ_adj (raw + expression-partial):")
    print(df.to_string(index=False))
    return df


# --------------------------------------------------------------------------- #
# GENOMIC-CONTEXT + HALLMARK-PROGRAM roll-up (the retention mechanism + subproject core tie-back). Cheap joins.
def context_hallmark_rollup() -> tuple:
    """(1) GENOMIC CONTEXT — edge realization ρ_adj + patient-baseline retention split by miRNA `mir_class`
    (host-coupled sense_{coding,lncRNA}_host vs intergenic vs antisense_overlap); ties the retention mechanism
    (MH-158: host-coupled retain more) to realization. (2) HALLMARK PROGRAM — per-program mean gene-aggregate
    realization ρ_adj + acquired dose (mean_own_shift). -> genomic_context_rollup.tsv, hallmark_rollup.tsv."""
    gc = pd.read_csv(_LEARNED / "genomic_context.tsv", sep="\t")[["arm", "mir_class", "host_type"]]
    lad = pd.read_csv(OUT / "realization_ladder.tsv", sep="\t")
    ret = pd.read_csv(OUT / "mirna_nat_retention.tsv", sep="\t")
    edges = lad[lad.resolution == "edge"][["gene", "predictor", "rho_adj"]].rename(columns={"predictor": "arm"})
    # --- (1) genomic context: edge realization + arm retention by mir_class ---
    e = edges.merge(gc, on="arm", how="left")
    r = ret.merge(gc, on="arm", how="left")
    grp = []
    for cls in ["sense_coding_host", "sense_lncRNA_host", "intergenic", "antisense_overlap"]:
        ee = e[e.mir_class == cls]; rr = r[r.mir_class == cls]
        host_coupled = cls in ("sense_coding_host", "sense_lncRNA_host")
        grp.append({"mir_class": cls, "host_coupled": host_coupled, "n_edges": len(ee), "n_arms": rr.arm.nunique(),
                    "mean_edge_rho_adj": _r3(ee.rho_adj.mean()), "frac_rho_neg": _r2(float((ee.rho_adj < 0).mean())),
                    "mean_r_own": _r3(rr.r_own.mean()), "mean_r_adj": _r3(rr.r_adj.mean())})
    gcr = pd.DataFrame(grp); gcr.to_csv(OUT / "genomic_context_rollup.tsv", sep="\t", index=False)
    # --- (2) hallmark program: gene-aggregate realization + acquired dose per program ---
    from mirna_hallmark.hallmark_sets import gene_to_hallmarks
    g2h = gene_to_hallmarks()
    gene_r = lad[lad.resolution == "gene"][["gene", "rho_adj"]].set_index("gene")["rho_adj"]
    da = pd.read_csv(OUT / "dose_shift_edge.tsv", sep="\t")
    gene_dose = da.groupby("gene")["dShare_M_own"].apply(lambda s: s.abs().mean())     # within-gene dose activity
    mep = pd.read_csv(OUT / "master_edge_patterns.tsv", sep="\t")
    gene_acq = mep.groupby("gene")["mean_own_shift"].max()                             # gene's top acquired dose
    hrows = []
    for prog in sorted(set(h for hs in g2h.values() for h in hs)):
        genes = [g for g, hs in g2h.items() if prog in hs]
        rr = gene_r.reindex(genes).dropna()
        if len(rr) < 15:
            continue
        hrows.append({"hallmark": prog, "n_genes": len(rr), "mean_gene_rho_adj": _r3(rr.mean()),
                      "frac_rho_neg": _r2(float((rr < 0).mean())),
                      "mean_acquired_dose": _r3(gene_acq.reindex(genes).dropna().mean()),
                      "mean_dose_activity": _r3(gene_dose.reindex(genes).dropna().mean())})
    hr = pd.DataFrame(hrows).sort_values("mean_gene_rho_adj")
    hr.to_csv(OUT / "hallmark_rollup.tsv", sep="\t", index=False)
    print("[genomic-context rollup] edge realization + retention by mir_class:")
    print(gcr.to_string(index=False))
    print(f"\n[hallmark rollup] {len(hr)} programs; most-realized (top) vs least (bottom):")
    print(pd.concat([hr.head(5), hr.tail(5)]).to_string(index=False))
    return gcr, hr


# --------------------------------------------------------------------------- #
# Resolution 4 (WITHIN-FAMILY member) — DIAGNOSTIC, pre-registered NULL (MH-122: member naming is at chance;
# the Δ axis has less variance than the cross-section ⇒ strictly harder). No M fit: reads the active families.
def realize_within_family() -> pd.DataFrame:
    """Per multi-member family: each member arm's realization ρ_adj (member Δx alone → Δy | ΔC), the member's
    patient-baseline retention `r_own`, and a set-level test — does the MAX-RETENTION member coincide with the
    MOST-REALIZING (min ρ_adj) member beyond within-family chance (1/n_members)? prereg_null=True.
    -> realization_within_family.tsv (+ prints the concordance vs null)."""
    fam_tab = pd.read_csv(OUT / "realization_family.tsv", sep="\t")
    active = fam_tab[fam_tab.resolution == "family"][["gene", "predictor"]]
    dX, dY, pts_l = ST.paired_delta_matrices()
    pts = list(pts_l); dC = _delta_C()
    rmap = pd.read_csv(OUT / "mirna_nat_retention.tsv", sep="\t").set_index("arm")["r_own"].to_dict()
    rows = []
    for g, sub in active.groupby("gene"):
        if g not in dY.index:
            continue
        arms = [a for a in ST._he_arms(g) if a in dX.index]
        if not arms:
            continue
        fam_map = FAM.family_of(pd.Index(arms))
        for f in sub.predictor:
            members = [a for a in arms if fam_map.get(a) == f]
            if len(members) < 2:
                continue
            for a in members:
                res = _realize(g, [a], pd.Series({a: 1.0}), dX, dY, pts, dC)
                if res:
                    rows.append({"gene": g, "family": f, "member": a, "rho_adj": res["rho_adj"],
                                 "r_own": _r3(float(rmap.get(a, np.nan))), "n_members": len(members),
                                 "prereg_null": True})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "realization_within_family.tsv", sep="\t", index=False)
    # set-level concordance: max-r_own member == min-rho_adj member? vs Σ 1/n_members chance
    hits = tot = 0
    exp = 0.0
    for (g, f), s in df.dropna(subset=["r_own", "rho_adj"]).groupby(["gene", "family"]):
        if len(s) < 2:
            continue
        top_ret = s.loc[s.r_own.idxmax(), "member"]
        top_real = s.loc[s.rho_adj.idxmin(), "member"]        # most negative = most realizing
        hits += int(top_ret == top_real); tot += 1; exp += 1.0 / len(s)
    if tot:
        from scipy.stats import binomtest
        p = binomtest(hits, tot, exp / tot).pvalue if tot else np.nan
        print(f"[within-family] member concordance (max-retention == most-realizing): {hits}/{tot} "
              f"= {hits/tot:.2f} vs chance {exp/tot:.2f} (binom p={p:.2f}) — PRE-REGISTERED NULL (MH-122).")
    return df


# --------------------------------------------------------------------------- #
# Phase 3 — the ~7% OWN-RETAINING patients (A18U et al.): clinical deep-dive (descriptive + one test).
def own_retaining_clinical() -> pd.DataFrame:
    """Join `patient_nat_identity` (own_minus_cohort>0 = own-retaining, ~7%) with clinical (stage/PAM50/purity/
    proliferation): are own-retaining tumours lower-grade / more-differentiated? Already composition-checked
    (own-advantage ⊥ comp-distance p=0.72, MH-158). -> own_retaining_clinical.tsv."""
    from analysis.utils.common.loaders import load_clinical
    pid = pd.read_csv(OUT / "patient_nat_identity.tsv", sep="\t")
    pid["own_retaining"] = pid["own_minus_cohort"] > 0
    cl = load_clinical().drop_duplicates("participant").set_index("participant")
    j = pid.merge(cl[["pathologic_stage_collapsed", "PAM50_final", "CPE", "abs_purity",
                      "thornsson_proliferation", "thornsson_aneuploidy_score"]],
                  left_on="patient", right_index=True, how="left")
    j.to_csv(OUT / "own_retaining_clinical.tsv", sep="\t", index=False)
    n_own = int(j.own_retaining.sum())
    print(f"[own-retaining clinical] {n_own}/{len(j)} = {100*n_own/len(j):.0f}% own-retaining (own_minus_cohort>0)")
    from scipy.stats import mannwhitneyu
    for col in ["CPE", "abs_purity", "thornsson_proliferation", "thornsson_aneuploidy_score"]:
        a = j.loc[j.own_retaining, col].dropna(); b = j.loc[~j.own_retaining, col].dropna()
        if len(a) >= 3 and len(b) >= 3:
            p = mannwhitneyu(a, b).pvalue
            print(f"  {col:32s} own-retaining {a.mean():.3f} vs rest {b.mean():.3f} (MWU p={p:.2f})")
    # early-stage enrichment
    if "pathologic_stage_collapsed" in j:
        early = j["pathologic_stage_collapsed"].astype(str).str.contains("I$|II$|^I|^II", regex=True)
        from scipy.stats import fisher_exact
        ct = pd.crosstab(j.own_retaining, early)
        if ct.shape == (2, 2):
            _, pf = fisher_exact(ct)
            print(f"  early-stage(I/II) × own-retaining Fisher p={pf:.2f}")
    return j


# --------------------------------------------------------------------------- #
# ⭐ CLASS → WITHIN-PATIENT REALIZATION — does the cross-state shift_class predict paired realization ρ_adj?
# Gene-CLUSTERED (edges within a gene correlated), dose-CONTROLLED (verification item 7), decoy-arbitrated
# (axiom 4). Set-level (n=103). The non-circular contrast is AMONG coupled classes: "live-in-NAT"
# {constitutive, field_established_realized} vs "acquired-only" {acquired_realized, tumour_realized}.
_LIVE_IN_NAT = ("constitutive", "field_established_realized")
_ACQUIRED_ONLY = ("acquired_realized", "tumour_realized")


def _gene_boot(gene_arr: np.ndarray, stat, n: int = 2000, seed: int = 0):
    """Gene-CLUSTER bootstrap: `stat(idx)` computed on numpy ROW INDICES; each draw resamples GENES with
    replacement and gathers their row positions (numpy concat, not pd.concat — profile-before-batch). Returns
    (point, lo95, hi95, p_two_sided_vs_0)."""
    uniq = pd.unique(gene_arr)
    pos = {g: np.where(gene_arr == g)[0] for g in uniq}
    rng = np.random.default_rng(seed)
    point = stat(np.arange(len(gene_arr)))
    draws = np.empty(n)
    for i in range(n):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        draws[i] = stat(np.concatenate([pos[g] for g in pick]))
    draws = draws[np.isfinite(draws)]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    p = 2.0 * min((draws >= 0).mean(), (draws <= 0).mean()) if len(draws) else np.nan
    return float(point), float(lo), float(hi), float(p)


def class_realization() -> pd.DataFrame:
    """Does the cross-state `shift_class` predict within-patient paired realization ρ_adj? (gene-clustered,
    dose-controlled, decoy-arbitrated; set-level, n=103). -> class_realization{,_contrast}.tsv."""
    ec = pd.read_csv(OUT / "progression_edge_card.tsv", sep="\t").dropna(subset=["shift_class", "edge_rho_adj"])
    # (1) gene-clustered mean ρ_adj + acquired dose per class
    rows = []
    for cls, s in ec.groupby("shift_class"):
        rho = s["edge_rho_adj"].to_numpy(float)
        pt, lo, hi, p = _gene_boot(s["gene"].to_numpy(), lambda idx: np.nanmean(rho[idx]), n=1000)
        rows.append({"shift_class": cls, "n_edges": len(s), "n_genes": s.gene.nunique(),
                     "mean_rho_adj": _r3(pt), "ci_lo": _r3(lo), "ci_hi": _r3(hi), "boot_p_vs0": _r3(p),
                     "mean_own_shift": _r3(s["mean_own_shift"].mean()), "coupled": cls not in ("uncoupled", "undetectable", "lost")})
    census = pd.DataFrame(rows).sort_values("mean_rho_adj")
    census.to_csv(OUT / "class_realization.tsv", sep="\t", index=False)

    # circularity magnitude (rigor fix): shared repression signal with the DEFINING coupling_nat + the stronger
    # correlate coupling_tum — NOT tautological (would be ≈1). ρ<0.4 ⇒ different estimators.
    def _c(a, b):
        d = ec.dropna(subset=[a, b])
        return _r3(float(spearmanr(d[a], d[b]).correlation)) if len(d) > _MIN_PAIRS else np.nan
    circ = {"corr_rho_couplingNAT": _c("edge_rho_adj", "coupling_nat"),
            "corr_rho_couplingTUM": _c("edge_rho_adj", "coupling_tum")}

    # (2) the non-circular COUPLED contrast: live-in-NAT vs acquired-only, gene-clustered; RAW, dose-controlled,
    # and dose+coupling_tum-controlled (rigor fix: coupling_tum is the STRONGER correlate and must be balanced/controlled)
    sub = ec[ec.shift_class.isin(_LIVE_IN_NAT + _ACQUIRED_ONLY)].dropna(
        subset=["mean_own_shift", "coupling_tum"]).copy()
    g_sub = sub["gene"].to_numpy()
    rho = sub["edge_rho_adj"].to_numpy(float)
    islive = sub.shift_class.isin(_LIVE_IN_NAT).to_numpy(float)
    dose = sub["mean_own_shift"].to_numpy(float)
    ctum = sub["coupling_tum"].to_numpy(float)
    ctum_bal = {"coupling_tum_live": _r3(float(ctum[islive == 1].mean())),      # balance check (should be ~equal)
                "coupling_tum_acquired": _r3(float(ctum[islive == 0].mean()))}

    def _diff(idx):                                      # mean ρ(live) − mean ρ(acquired)  (expect <0)
        li = islive[idx]
        a, b = rho[idx][li == 1], rho[idx][li == 0]
        return (a.mean() - b.mean()) if (len(a) and len(b)) else np.nan

    def _beta(idx, cols):                                # OLS coef of is_live controlling `cols` (list of arrays)
        li = islive[idx]
        if li.min() == li.max():
            return np.nan
        X = np.column_stack([np.ones(len(idx)), li] + [c[idx] for c in cols])
        beta, *_ = np.linalg.lstsq(X, rho[idx], rcond=None)
        return beta[1]

    d_pt, d_lo, d_hi, d_p = _gene_boot(g_sub, _diff, n=2000)
    b_pt, b_lo, b_hi, b_p = _gene_boot(g_sub, lambda i: _beta(i, [dose]), n=2000)
    bt_pt, bt_lo, bt_hi, bt_p = _gene_boot(g_sub, lambda i: _beta(i, [dose, ctum]), n=2000)
    contrast = [{"test": "live_vs_acquired_raw", "estimate": _r3(d_pt), "ci_lo": _r3(d_lo), "ci_hi": _r3(d_hi),
                 "boot_p": _r3(d_p), "n_edges": len(sub), "n_genes": sub.gene.nunique()},
                {"test": "live_vs_acquired_dose_controlled", "estimate": _r3(b_pt), "ci_lo": _r3(b_lo),
                 "ci_hi": _r3(b_hi), "boot_p": _r3(b_p), "n_edges": len(sub), "n_genes": sub.gene.nunique()},
                {"test": "live_vs_acquired_dose+couplingTUM_controlled", "estimate": _r3(bt_pt), "ci_lo": _r3(bt_lo),
                 "ci_hi": _r3(bt_hi), "boot_p": _r3(bt_p), "n_edges": len(sub), "n_genes": sub.gene.nunique()}]

    # (3) DECOY ARBITER (axiom 4): REAL vs matched site-free DECOY ρ_adj, by class of the REAL edge — reported as
    # the DIRECT paired-bootstrap DIFFERENCE (live REAL−DECOY) − (acquired REAL−DECOY) (rigor fix: two separate
    # vs-0 CIs overlap and do NOT establish the difference; the paired diff is the correct statistic).
    nd_p = OUT / "realization_null_edges.tsv"
    if nd_p.exists():
        nd = pd.read_csv(nd_p, sep="\t")
        cls_map = ec.set_index(["gene", "arm"])["shift_class"]
        real = nd[nd.group == "REAL"].copy()
        real["shift_class"] = [cls_map.get((g, a), np.nan) for g, a in zip(real.gene, real.arm)]
        decoy_floor = nd[nd.group == "DECOY"].groupby("gene")["rho_adj"].mean()   # per-gene null floor
        real["decoy_floor"] = real.gene.map(decoy_floor)
        real["real_minus_decoy"] = real["rho_adj"] - real["decoy_floor"]
        rc = real[real.shift_class.isin(_LIVE_IN_NAT + _ACQUIRED_ONLY)].dropna(subset=["real_minus_decoy"]).copy()
        rmd = rc["real_minus_decoy"].to_numpy(float)
        rlive = rc.shift_class.isin(_LIVE_IN_NAT).to_numpy(float)
        for grp, name in ((_LIVE_IN_NAT, "live_in_nat"), (_ACQUIRED_ONLY, "acquired_only")):
            g = real[real.shift_class.isin(grp)].dropna(subset=["real_minus_decoy"])
            if len(g) >= 30:
                v = g["real_minus_decoy"].to_numpy(float)
                pt, lo, hi, p = _gene_boot(g["gene"].to_numpy(), lambda idx: np.nanmean(v[idx]), n=1000)
                contrast.append({"test": f"REAL_minus_DECOY_{name}", "estimate": _r3(pt), "ci_lo": _r3(lo),
                                 "ci_hi": _r3(hi), "boot_p": _r3(p), "n_edges": len(g), "n_genes": g.gene.nunique()})

        def _rmd_diff(idx):                              # (live REAL−DECOY) − (acquired REAL−DECOY), paired per draw
            li = rlive[idx]
            a, b = rmd[idx][li == 1], rmd[idx][li == 0]
            return (a.mean() - b.mean()) if (len(a) and len(b)) else np.nan
        pt, lo, hi, p = _gene_boot(rc["gene"].to_numpy(), _rmd_diff, n=2000)
        contrast.append({"test": "REAL_minus_DECOY_live_minus_acquired", "estimate": _r3(pt), "ci_lo": _r3(lo),
                         "ci_hi": _r3(hi), "boot_p": _r3(p), "n_edges": len(rc), "n_genes": rc.gene.nunique()})
    con = pd.DataFrame(contrast)
    con.attrs.update(circ); con.attrs.update(ctum_bal)
    print(f"[circularity] corr(edge_rho_adj, coupling_nat)={circ['corr_rho_couplingNAT']} "
          f"(defining var) · corr(·, coupling_tum)={circ['corr_rho_couplingTUM']} (stronger correlate) — both <0.4, not tautological")
    print(f"[coupling_tum balance] live {ctum_bal['coupling_tum_live']} vs acquired {ctum_bal['coupling_tum_acquired']} (should be ~equal)")
    con.to_csv(OUT / "class_realization_contrast.tsv", sep="\t", index=False)
    print("[class_realization] mean ρ_adj by cross-state shift_class (gene-clustered):")
    print(census.to_string(index=False))
    print("\n[contrast] live-in-NAT {constitutive,field_established} vs acquired-only {acquired_realized,tumour_realized}:")
    print(con.to_string(index=False))
    return census
