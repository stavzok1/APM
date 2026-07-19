"""Per-edge (and per-gene) attribution card — the taxonomy rendered per edge on the confounder-honest,
cell-intrinsic foundation. Column/estimator spec: docs/FORMULAS.md (§7 coupling ladder, §11f `shift_class`
calibrated two-axis) + docs/LEARNED_MODEL_METHODS.md; findings: DISCOVERY_REGISTRY MH-158/160/166.
One row per regulator of a gene, joining the four layers we already compute:

  1. REGIME (range stats) — arm median RPM, %>functional-floor (RPM≥10), IQR, SPIKER flag (low median +
     high IQR = subset-driven, e.g. miR-135b 28%>floor); target IQR (a strong coupling on a low-IQR gene =
     the miRNA explains most of its variance).
  2. BUDGET share (E7/G4) — the arm's rank + share of the gene's pressure, GTEx→NAT→tumour Δ (states.budget_shift).
  3. COMPOSITION tag — deconv retention: cell-intrinsic (≥0.7) / partial / composition-explained (<0.4).
  4. SHIFT-CLASS — the calibrated two-axis progression class (MH-166): per-state CALIBRATED coupling
     (site_free_null p<0.05, replaces the −0.1 cut) × same-platform `arm_lfc_NAT_TUM` dose. New class
     `dose_acquired_uncoupled` (dose up, not calibrated-coupled) + quantified `realization_score`,
     `coupling_p_*`, `wiring_frac`. Spec: docs/FORMULAS.md §11f. Learned successor to the §6b-RETIRED
     `mirna_state_class.joint_edge_class` — do NOT read the two as the same object.

CLI: `python -m mirna_hallmark.learned.card PTEN GATA3`
"""
from __future__ import annotations

import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

warnings.filterwarnings("ignore", message="An input array is constant")   # constant-arm edges → ρ=nan (handled)

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import states as ST

_FLOOR = np.log2(11)   # RPM ≥ 10 functional floor on the log2(RPM+1) scale


def _resid(v, C):
    return v - LinearRegression().fit(C, v).predict(C)


# --------------------------------------------------------------------------------------------------------- #
# CALIBRATED coupling null (MH-123/154) — replaces the hard −0.1 repressor cut. The theoretical null is 3–4×
# too narrow; the site-free empirical null gives an honest per-state p. Fit ONCE per state (memoized + disk
# cached); the parent warms it before forking so workers inherit via COW (axiom 3a).
_NULL_CACHE: dict = {}
COUPLING_ALPHA = 0.05          # calibrated per-state "repressor in this state" call (site-free null p < 0.05)


def _null_scorers():
    if "lad" not in _NULL_CACHE:
        from mirna_hallmark.eval import site_free_null as SFN
        from mirna_hallmark.learned import state as STA
        _NULL_CACHE["lad"] = {"tum": SFN.fit_state("tumour"), "nat": SFN.fit_state("nat"),
                              "hly": SFN.fit_state("gtex")}
        _NULL_CACHE["ab"] = {"tum": LD._load()["X"].median(axis=1),
                             "nat": ST.state_matrices("11")[0].median(axis=1),
                             "hly": STA._gtex_mirna().median(axis=1)}
    return _NULL_CACHE["lad"], _NULL_CACHE["ab"]


def _healthy_leg_map():
    """arm → (healthy-leg provenance, imputed healthy abundance = repression-POTENTIAL proxy), so the `h` (GTEx)
    coupling flag is not read as evidence where the GTEx uniquely-mappable pipeline COLLAPSED the arm (raw GTEx
    coupling matrix floored ⇒ ρ≈0 mechanically, NOT biology). Cross-references the DIANA-miTED healthy-breast cohort
    (`healthy_anchor`): a floored-in-GTEx arm that miTED shows expressed is a COLLAPSE ARTIFACT, its 'uncoupled-in-
    healthy' call blind not absent (`collapse_blind`). BUT blindness only THREATENS an 'acquired' call when the arm is
    HIGHLY expressed in healthy: repression is abundance-gated (AGO/RISC loading saturates), so a LOW-abundance healthy
    arm has low repression potential whether or not coupling is measurable ⇒ its acquisition is still safe. So we carry
    a GRADED `healthy_potential` (imputed baseline, collapse-repaired) and reserve `healthy_uninformative` for
    collapse_blind AND high potential. (MH-166 follow-up.)  provenance ∈ {measured, collapse_blind, true_absent}."""
    if "hleg" not in _NULL_CACHE:
        from mirna_hallmark.learned import state as STA
        from mirna_hallmark import healthy_anchor as HA
        Xg = STA._gtex_mirna()
        gmed = Xg.median(axis=1)
        gfloor = (Xg.eq(Xg.min(axis=1), axis=0)).mean(axis=1)          # frac samples at row-min (tie/floor)
        mited = HA._mited_breast_median()
        mited = HA._reconcile_to_tcga(mited, gmed.index) if not mited.empty else pd.Series(dtype=float)
        FL = HA.ABUND_FLOOR
        prov = {}
        for a in gmed.index:
            floored = (gmed.get(a, 0) <= FL) or (gfloor.get(a, 1.0) > 0.5)
            if not floored:
                prov[a] = "measured"
            elif (mited.get(a, 0) > FL) or (gmed.get(a, 0) > FL):
                prov[a] = "collapse_blind"                             # expressed somewhere, blind in GTEx coupling
            else:
                prov[a] = "true_absent"
        pot = HA.gtex_qn_baseline()                                    # imputed healthy level (collapse-repaired)
        # "high potential" = above the project's RPM≥10 FUNCTIONAL-PRESENCE floor (_FLOOR = log2(11)). A biological
        # anchor, not a population tercile — avoids the axiom-5 pile-up a tercile cut would sit on. An arm functionally
        # present in healthy has real repression potential ⇒ its blind coupling THREATENS an acquired call.
        _NULL_CACHE["hleg"] = (prov, pot.to_dict(), float(_FLOOR))
    return _NULL_CACHE["hleg"]


SURROGATE_MIN_CORR = 0.5   # a same-seed instrument must track the collapsed arm this well (tumour, where both measured)


def _surrogate_map():
    """arm → (instrument, tumour_corr) for arms COLLAPSED in GTEx that have a SAME-SEED mate MEASURED in GTEx. The
    collapse hits multi-mapping paralogues, whose polycistronic/paralog cluster-mates ARE uniquely-mappable — and share
    the SEED, hence the SAME targets with the SAME sites. So the mate's GTEx coupling directly estimates the seed's
    healthy repression (recovers the variance the collapse destroyed; miTED-median cannot). VALIDATED where both are
    measured (tumour corr ≥ SURROGATE_MIN_CORR) before transfer. Recovers the blind healthy-coupling node for the
    ~7–30 well-instrumented flagged arms (miR-20a→miR-17, miR-181b→miR-181a). (MH-166 follow-up.)"""
    if "surr" not in _NULL_CACHE:
        from mirna_hallmark.learned import state as STA
        from mirna_hallmark.learned import families as FAM
        Xt = LD._load()["X"]
        gmed = STA._gtex_mirna().median(axis=1)
        fam = FAM.family_of(pd.Index(Xt.index))
        collapsed = [a for a in Xt.index if gmed.get(a, 0.0) <= HA_FL() and Xt.loc[a].median() > 1]
        smap = {}
        for a in collapsed:
            f = fam.get(a)
            if f is None:
                continue
            mates = [m for m in Xt.index if m != a and fam.get(m) == f and gmed.get(m, 0.0) > HA_FL()]
            if not mates:
                continue
            best = max(mates, key=lambda m: abs(spearmanr(Xt.loc[a], Xt.loc[m]).correlation))
            r = float(spearmanr(Xt.loc[a], Xt.loc[best]).correlation)
            if r >= SURROGATE_MIN_CORR:
                smap[a] = (best, round(r, 3))
        _NULL_CACHE["surr"] = smap
    return _NULL_CACHE["surr"]


def HA_FL():
    from mirna_hallmark import healthy_anchor as HA
    return HA.ABUND_FLOOR


def _coupling_calibrated(state_key: str, arm: str, rho: float):
    """(one-sided repressive p, z) of a coupling ρ against the state's site-free null, abundance-matched.
    p < COUPLING_ALPHA ⇒ significantly repressive vs pairs that CANNOT repress (calibrated, not the −0.1 cut)."""
    if rho != rho:
        return np.nan, np.nan
    lad, ab = _null_scorers()
    a = ab[state_key].get(arm, np.nan)
    if a != a:
        return np.nan, np.nan
    p = float(lad[state_key].pvalues([rho], [a])[0])
    z = float(lad[state_key].zscores([rho], [a])[0])
    return p, z


def _edge_coupling(gene: str, arm: str, *, deconv: bool):
    """Single-edge partial-Spearman of the arm vs the gene | C (core, or +composition)."""
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True, deconv=deconv)
    if arm not in X.columns:
        return np.nan
    Cm = C.to_numpy(float)
    return spearmanr(_resid(X[arm].to_numpy(float), Cm), _resid(Y.to_numpy(float), Cm)).correlation


def _gtex_coupling(gene: str, arm: str):
    """Single-edge coupling in GTEx-healthy (metagene C from GTEx mRNA) — the healthy realization axis."""
    from mirna_hallmark.learned import state as STA
    Xg, Yg = STA._gtex_mirna(), STA._gtex_mrna()
    if arm not in Xg.index or gene not in Yg.index:
        return np.nan
    Cg = ST._state_metagene_cov(Yg)
    p = [c for c in Yg.columns if c in Xg.columns and c in Cg.index]
    if len(p) < 25:
        return np.nan
    xv, yv = Xg.loc[arm, p].to_numpy(float), Yg.loc[gene, p].to_numpy(float)
    Cm = Cg.loc[p].apply(lambda s: s.fillna(s.median())).to_numpy(float)
    m = np.isfinite(xv) & np.isfinite(yv) & np.isfinite(Cm).all(axis=1)
    if m.sum() < 25:
        return np.nan
    return spearmanr(_resid(xv[m], Cm[m]), _resid(yv[m], Cm[m])).correlation


def _shift_class(row) -> str:
    """⭐ TWO-AXIS progression class (2026-07-19 rebuild): the COUPLING trajectory (CALIBRATED per-state site-free
    null p<COUPLING_ALPHA, NOT the −0.1 cut) × the SAME-PLATFORM DOSE trajectory (`arm_lfc_NAT_TUM`, NOT the soft QN
    healthy leg). Separates *when repression engages* from *when the miRNA is dosed up* — the two the old label
    conflated. NEW class `dose_acquired_uncoupled` = dose rises but no realized repression (potential w/o
    realization; the exact case the old label mislabelled `constitutive`). Healthy leg (h) is calibrated but SOFT
    (GTEx cross-platform, lowest-abundance bin unfit). Quantified companions: `realization_score`, `dose_class`,
    `wiring_frac`. Supersedes the coupling-only −0.1 label."""
    ret, spk = row["retention"], row["spiker"]
    rep = lambda p: (p == p) and p < COUPLING_ALPHA              # CALIBRATED: significantly repressive vs site-free
    h_p = row.get("coupling_p_hly_resolved", row.get("coupling_p_hly"))   # direct, else same-seed surrogate
    h, n, t = rep(h_p), rep(row.get("coupling_p_nat")), rep(row.get("coupling_p_tum"))
    dose = row.get("arm_lfc_NAT_TUM", np.nan)                    # SAME-PLATFORM NAT→tumour level shift
    dose_gain = (dose == dose) and dose > 0.3
    if not t:                                                    # not calibrated-repressor in tumour
        if spk or (row["arm_pct_floor"] == row["arm_pct_floor"] and row["arm_pct_floor"] < 30):
            return "undetectable"
        if dose_gain and not (h or n):
            return "dose_acquired_uncoupled"                    # ⭐ dose up, no realized repression (potential only)
        return "lost" if (h or n) else "uncoupled"
    if ret == ret and ret < 0.4:
        return "composition_explained"
    if h and n:
        return "constitutive"                                   # calibrated-repressor in healthy already
    if (not h) and n:
        return "field_established_realized"                     # established in NAT (field effect)
    if (not h) and (not n):
        return "acquired_realized" if dose_gain else "tumour_realized"  # tumour-specific; dose-gated by same-platform
    return "nat_decoupled" if (h and not n) else "stable"


def gene_card(gene: str, *, alpha: float = 0.005) -> pd.DataFrame:
    """One row per regulator of `gene` with regime + budget + composition + shift-class."""
    bdf = ST.budget_shift(gene, alpha=alpha)                     # rank/share GTEx→NAT→tumour (cell-intrinsic M)
    if bdf.empty:
        return bdf
    X = LD._load()["X"]; Y = LD._load()["Y"]
    yg = Y.loc[gene].dropna() if gene in Y.index else pd.Series(dtype=float)
    gene_iqr = float(yg.quantile(.75) - yg.quantile(.25)) if len(yg) else np.nan
    # --- assemble the tumour design ONCE (raw + deconv, orphans) and reuse across ALL arms; per-arm partial
    # couplings are then a cheap masked spearman, not a 0.7s re-assemble each (profile-before-batch-loops). ---
    Yr, Xr, Cr, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True)
    Yd, Xd, Cd, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True, deconv=True)
    Crm = Cr.to_numpy(float); yrr = _resid(Yr.to_numpy(float), Crm)   # assemble_gene NaN-cleans → aligned/finite
    Cdm = Cd.to_numpy(float); yrd = _resid(Yd.to_numpy(float), Cdm)

    def _edge(arm, Xa, Cm, yres):
        return spearmanr(_resid(Xa[arm].to_numpy(float), Cm), yres).correlation if arm in Xa.columns else np.nan

    # GTEx-healthy single-edge coupling — precompute the metagene C + per-gene residual target ONCE
    from mirna_hallmark.learned import state as STA
    Xg, Yg = STA._gtex_mirna(), STA._gtex_mrna()
    if gene in Yg.index:
        Cg = ST._state_metagene_cov(Yg)
        pg = [c for c in Yg.columns if c in Xg.columns and c in Cg.index]
        Cgm = Cg.loc[pg].apply(lambda s: s.fillna(s.median())).to_numpy(float) if len(pg) >= 25 else None
        ygv = Yg.loc[gene, pg].to_numpy(float) if len(pg) >= 25 else None
    else:
        pg, Cgm, ygv = [], None, None

    def _gtex_edge(arm):
        if Cgm is None or arm not in Xg.index:
            return np.nan
        xv = Xg.loc[arm, pg].to_numpy(float)
        m = np.isfinite(xv) & np.isfinite(ygv) & np.isfinite(Cgm).all(axis=1)
        return spearmanr(_resid(xv[m], Cgm[m]), _resid(ygv[m], Cgm[m])).correlation if m.sum() >= 25 else np.nan

    # NAT single-edge coupling (state-comparable CIBERSORTx/metagene C) — the realization axis in NAT
    Xn, Yn = ST.state_matrices("11")
    Cn = ST._cibersortx_state_cov(list(Yn.columns), "11")
    if Cn is None:
        Cn = ST._state_metagene_cov(Yn).reindex(Yn.columns)
    Cn = Cn.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))

    def _nat_coupling(arm):
        if arm not in Xn.index or gene not in Yn.index:
            return np.nan
        p = [c for c in Yn.columns if c in Xn.columns and c in Cn.index]
        if len(p) < 25:
            return np.nan
        xv = Xn.loc[arm, p].to_numpy(float); yv = Yn.loc[gene, p].to_numpy(float)
        Cm = Cn.loc[p].to_numpy(float)
        m = np.isfinite(xv) & np.isfinite(yv) & np.isfinite(Cm).all(axis=1)
        if m.sum() < 25:
            return np.nan
        return spearmanr(_resid(xv[m], Cm[m]), _resid(yv[m], Cm[m])).correlation

    # state mean levels for logFC + QN'd healthy baseline (TCGA scale) for the healthy→tumour logFC
    Xt2, Yt2 = ST.state_matrices("01")
    from mirna_hallmark import healthy_anchor as HA
    try:
        hbase = HA.load_baseline()
    except Exception:
        hbase = pd.Series(dtype=float)
    tg = float(Yt2.loc[gene].mean()) if gene in Yt2.index else np.nan
    ng = float(Yn.loc[gene].mean()) if gene in Yn.index else np.nan
    gene_lfc = round(tg - ng, 2) if (tg == tg and ng == ng) else np.nan   # target NAT→tumour logFC (direct)
    # GLOBAL abundance rank (percentile among ALL miRNAs, per state) — the mirna_state_class level axis,
    # complementary to the gene-centric budget rank (share among the GENE's regulators). QN-safe (percentile).
    # ⭐ rank on the MEDIAN, not the mean — the mean is pulled up by rare-spike arms (miR-122-3p rank 15→81 on the
    # mean), overstating arms typically absent but occasionally high; median is robust and consistent with the card's
    # regime/spiker axes (`arm_med_rpm`). (MH-166 follow-up; ~15% of arms shift >10 pctile, corr 0.96.)
    gr_tum = Xt2.median(axis=1).rank(pct=True) * 100
    gr_nat = Xn.median(axis=1).rank(pct=True) * 100
    try:
        from mirna_hallmark.learned import state as STA
        gmean = STA._gtex_mirna().mean(axis=1)                 # raw GTEx per-arm mean — kept (mean) for the raw logFC
        gr_hly = STA._gtex_mirna().median(axis=1).rank(pct=True) * 100   # rank on median (robust to rare spikes)
    except Exception:
        gmean, gr_hly = pd.Series(dtype=float), pd.Series(dtype=float)

    # WIRING axis (SOFT — ΔM is n-limited at NAT): canonical bagged-family M per state → Δ(M·a) decomposition.
    try:
        Mt_w = ST.canonical_M(gene, "01", arm_level=True); Mn_w = ST.canonical_M(gene, "11", arm_level=True)
    except Exception:
        Mt_w = Mn_w = pd.Series(dtype=float)

    rows = []
    for r in bdf.itertuples():
        arm = r.arm
        xa = X.loc[arm].dropna() if arm in X.index else pd.Series(dtype=float)
        med = float(xa.median()) if len(xa) else np.nan
        pct = 100 * float((xa > _FLOOR).mean()) if len(xa) else np.nan
        iqr = float(xa.quantile(.75) - xa.quantile(.25)) if len(xa) else np.nan
        spiker = bool(pct < 40 and iqr > 1.5) if pct == pct else False
        ct = _edge(arm, Xr, Crm, yrr)
        cd = _edge(arm, Xd, Cdm, yrd)
        ret = float(cd / ct) if (ct and ct == ct) else np.nan
        cn = _nat_coupling(arm)
        ch = _gtex_edge(arm)
        tm = float(Xt2.loc[arm].mean()) if arm in Xt2.index else np.nan   # arm logFCs
        nm = float(Xn.loc[arm].mean()) if arm in Xn.index else np.nan
        hm = float(hbase.get(arm)) if arm in getattr(hbase, "index", []) else np.nan
        gm = float(gmean.get(arm, np.nan)) if arm in getattr(gmean, "index", []) else np.nan
        lfc_nt = round(tm - nm, 2) if (tm == tm and nm == nm) else np.nan          # NAT→tumour direct
        lfc_ht = round(tm - hm, 2) if (tm == tm and hm == hm) else np.nan          # healthy→tumour QN
        # RAW GTEx→tumour/NAT logFC (no QN) — for miRNAs (near-equal length, library dominated by ~200-300
        # arms) the cross-platform RPM/TPM difference is roughly interpretable; report alongside the QN.
        lfc_ht_raw = round(tm - gm, 2) if (tm == tm and gm == gm) else np.nan
        lfc_hn_raw = round(nm - gm, 2) if (nm == nm and gm == gm) else np.nan
        grt = float(gr_tum.get(arm, np.nan)); grn = float(gr_nat.get(arm, np.nan))
        grh = float(gr_hly.get(arm, np.nan)) if arm in getattr(gr_hly, "index", []) else np.nan
        # CALIBRATED coupling (site-free null p per state — replaces the −0.1 cut)
        cp_t, cz_t = _coupling_calibrated("tum", arm, ct)
        cp_n, _ = _coupling_calibrated("nat", arm, cn)
        cp_h, _ = _coupling_calibrated("hly", arm, ch)
        # SURROGATE healthy coupling for GTEx-collapsed arms via a same-seed MEASURED instrument (same seed ⇒ same
        # targets/sites ⇒ the mate's GTEx coupling estimates the seed's healthy repression) — RESOLVES the blind leg.
        ch_sur = cp_h_sur = sur_corr = np.nan; sur_inst = ""
        if _healthy_leg_map()[0].get(arm, "true_absent") == "collapse_blind" and arm in _surrogate_map():
            sur_inst, sur_corr = _surrogate_map()[arm]
            ch_sur = _gtex_edge(sur_inst)
            cp_h_sur, _ = _coupling_calibrated("hly", sur_inst, ch_sur)
        # WIRING decomposition Δ(M·a) = M_NAT·Δa (DOSE) + a_NAT·ΔM (WIRING) + Δa·ΔM (INTERACT)  [SOFT]
        mt_w, mn_w = float(Mt_w.get(arm, np.nan)), float(Mn_w.get(arm, np.nan))
        if mt_w == mt_w and mn_w == mn_w and tm == tm and nm == nm:
            da_w, dM_w = tm - nm, mt_w - mn_w
            t_ab, t_wi, t_in = mn_w * da_w, nm * dM_w, da_w * dM_w
            wf = abs(t_wi) / (abs(t_ab) + abs(t_wi) + abs(t_in) + 1e-9)
        else:
            t_ab = t_wi = t_in = wf = np.nan
        rows.append({"arm": arm, "share_TUM": r.share_TUM, "rank_TUM": r.rank_TUM,       # gene-centric budget
                     "d_rank_HLY_TUM": getattr(r, "d_rank_HLY_TUM", np.nan),
                     "grank_TUM": round(grt, 0) if grt == grt else np.nan,               # GLOBAL abundance %ile
                     "grank_NAT": round(grn, 0) if grn == grn else np.nan,
                     "grank_HLY": round(grh, 0) if grh == grh else np.nan,
                     "dGlobal_HLY_TUM": round(grt - grh, 0) if (grt == grt and grh == grh) else np.nan,
                     "arm_lfc_NAT_TUM": lfc_nt, "arm_lfc_HLY_TUM_QN": lfc_ht,
                     "arm_lfc_HLY_TUM_raw": lfc_ht_raw, "arm_lfc_HLY_NAT_raw": lfc_hn_raw,
                     "gene_lfc_NAT_TUM": gene_lfc,
                     "coupling_hly": round(ch, 3) if ch == ch else np.nan,
                     "arm_med_rpm": round(med, 1), "arm_pct_floor": round(pct, 0) if pct == pct else np.nan,
                     "arm_iqr": round(iqr, 1), "spiker": spiker,
                     "coupling_tum": round(ct, 3) if ct == ct else np.nan,
                     "coupling_nat": round(cn, 3) if cn == cn else np.nan,
                     "coupling_p_tum": round(cp_t, 4) if cp_t == cp_t else np.nan,   # CALIBRATED site-free-null p
                     "coupling_p_nat": round(cp_n, 4) if cp_n == cp_n else np.nan,
                     "coupling_p_hly": round(cp_h, 4) if cp_h == cp_h else np.nan,
                     "coupling_hly_surrogate": round(ch_sur, 3) if ch_sur == ch_sur else np.nan,  # same-seed instrument
                     "coupling_p_hly_surrogate": round(cp_h_sur, 4) if cp_h_sur == cp_h_sur else np.nan,
                     "surrogate_instrument": sur_inst, "surrogate_corr": sur_corr,
                     "coupling_z_tum": round(cz_t, 2) if cz_t == cz_t else np.nan,
                     "term_ABUND": round(t_ab, 3) if t_ab == t_ab else np.nan,       # dose-vs-wiring (SOFT)
                     "term_WIRING": round(t_wi, 3) if t_wi == t_wi else np.nan,
                     "term_INTERACT": round(t_in, 3) if t_in == t_in else np.nan,
                     "wiring_frac": round(wf, 2) if wf == wf else np.nan,
                     "healthy_leg": _healthy_leg_map()[0].get(arm, "true_absent"),  # GTEx-collapse provenance (miTED-aware)
                     "healthy_potential": round(_healthy_leg_map()[1].get(arm, np.nan), 2),  # imputed healthy level (potential)
                     "retention": round(ret, 2) if ret == ret else np.nan,
                     "gene_iqr": round(gene_iqr, 2)})
    card = pd.DataFrame(rows)
    # RESOLVED healthy coupling p: direct where measured; the same-seed surrogate where GTEx collapsed the arm. This is
    # what the class's healthy (h) leg reads — so a surrogate-recovered edge is classed on real evidence, not a blind NaN.
    card["coupling_p_hly_resolved"] = card["coupling_p_hly"].fillna(card["coupling_p_hly_surrogate"])
    card["shift_class"] = card.apply(_shift_class, axis=1)
    # ⚠ the healthy (h) leg is uninformative ONLY where GTEx collapsed the arm, the arm is HIGHLY expressed in healthy
    # (miTED, real repression potential), AND no same-seed surrogate recovered the coupling. A LOW-abundance collapse_blind
    # arm has abundance-bounded potential ⇒ safe; a surrogate-RESOLVED arm is no longer blind (MH-166 follow-up; graded).
    _hi = _healthy_leg_map()[2]
    card["healthy_uninformative"] = (card["healthy_leg"].eq("collapse_blind") & (card["healthy_potential"] > _hi)
                                     & card["coupling_p_hly_surrogate"].isna())
    # QUANTIFIED companions (the class is categorical; these are the continuous axes)
    card["realization_score"] = -card["coupling_z_tum"]         # calibrated tumour repression strength (>0 = repressive)
    card["dose_class"] = card["arm_lfc_NAT_TUM"].apply(          # SAME-PLATFORM NAT→tumour dose direction
        lambda d: "gain" if (d == d and d > 0.3) else ("loss" if (d == d and d < -0.3) else "flat"))
    card = card.sort_values("share_TUM", ascending=False)
    with pd.option_context("display.width", 200, "display.max_colwidth", 30):
        print(f"\n=== {gene} — per-edge attribution card (cell-intrinsic, confounder-honest) ===")
        print("gene-centric budget rank (rank_TUM/share) + GLOBAL abundance %ile (grank) + logFC + 3-state coupling:")
        print(card[["arm", "rank_TUM", "share_TUM", "grank_HLY", "grank_TUM", "dGlobal_HLY_TUM",
                    "arm_lfc_HLY_TUM_QN", "arm_lfc_HLY_TUM_raw", "spiker", "coupling_hly", "coupling_nat",
                    "coupling_tum", "retention", "shift_class"]].to_string(index=False))
    return card


def stable_wiring(gene: str) -> pd.DataFrame:
    """STABLE ΔM (wiring) — fixed FAMILY support across states + canonical bagged NNLS (no selection). ΔM is
    comparable because the estimand (the same families) is fixed in every state (fixes the cross-state lasso
    instability). Family-level view of states.canonical_M (arm_level=False)."""
    Ms = {}
    for st, lab in [("01", "T"), ("11", "NAT"), ("gtex", "HLY")]:
        m = ST.canonical_M(gene, st, arm_level=False)
        if len(m):
            Ms[lab] = m
    shared = sorted(set.intersection(*[set(m.index) for m in Ms.values()])) if len(Ms) > 1 else []
    df = pd.DataFrame({f"M_{s}": Ms[s].reindex(shared) for s in Ms}).round(3)
    if "M_T" in df and "M_HLY" in df:
        df["dM_HLY_TUM"] = (df["M_T"] - df["M_HLY"]).round(3)
    if "M_T" in df and "M_NAT" in df:
        df["dM_NAT_TUM"] = (df["M_T"] - df["M_NAT"]).round(3)
    df = df.sort_values("M_T", ascending=False) if "M_T" in df else df
    print(f"\n=== {gene} — STABLE ΔM (fixed family support + NNLS; {len(shared)} shared families) ===")
    print(df.to_string())
    return df


def state_wiring(gene: str, *, alpha: float = 0.01) -> pd.DataFrame:
    """The WEIGHT (M) shift itself, per edge, across states — the CANONICAL bagged family M (arm-broadcast)
    per state, compared. ΔM = the WIRING change (a·ΔM component: APA site-loss / AGO capacity), distinct from
    the abundance shift (M·Δa in the card). NOISY for small-n states (NAT 104, GTEx 327) — read ΔM only on
    arms with real M in ≥2 states. (Arm-level view of stable_wiring, with undetectable-state flags.)"""
    from mirna_hallmark.learned import state as STA
    Ms = {}
    for st, lab in [("01", "T"), ("11", "NAT"), ("gtex", "HLY")]:
        m = ST.canonical_M(gene, st, arm_level=True)
        if len(m) and m.abs().sum() > 0:
            Ms[lab] = m
    Xg, Yg = STA._gtex_mirna(), STA._gtex_mrna()
    idx = sorted(set().union(*[set(m.index) for m in Ms.values()]))
    df = pd.DataFrame({f"M_{s}": Ms[s].reindex(idx) for s in Ms}).fillna(0.0)
    if "M_T" in df and "M_HLY" in df:
        df["dM_HLY_TUM"] = df["M_T"] - df["M_HLY"]
    if "M_T" in df and "M_NAT" in df:
        df["dM_NAT_TUM"] = df["M_T"] - df["M_NAT"]
    # split M_state=0 into UNDETECTABLE (low per-state IQR = range-restricted) vs true-zero (user 2026-07-05)
    def _iqr(s):
        s = s.dropna(); return float(s.quantile(.75) - s.quantile(.25)) if len(s) else np.nan
    Xn2, Yn2 = ST.state_matrices("11")
    gih = _iqr(Yg.loc[gene]) if gene in Yg.index else np.nan
    gin = _iqr(Yn2.loc[gene]) if gene in Yn2.index else np.nan
    LOW = 1.0
    if "M_HLY" in df:
        df["undet_HLY"] = [bool(df.at[a, "M_HLY"] == 0 and ((a in Xg.index and _iqr(Xg.loc[a]) < LOW)
                           or (gih == gih and gih < LOW))) for a in df.index]
    if "M_NAT" in df:
        df["undet_NAT"] = [bool(df.at[a, "M_NAT"] == 0 and ((a in Xn2.index and _iqr(Xn2.loc[a]) < LOW)
                           or (gin == gin and gin < LOW))) for a in df.index]
    df = df.round(3).sort_values("M_T", ascending=False)
    print(f"\n=== {gene} — per-edge WEIGHT (M) shift across states (WIRING axis, z-scored abundance) ===")
    print(df[df.abs().sum(axis=1) > 0].to_string())
    return df


def decompose(gene: str, *, alpha: float = 0.005) -> pd.DataFrame:
    """NAT→tumour abundance/wiring decomposition per edge — both TCGA (clean gauge, no cross-platform issue):
        Δ(M·a) = M_NAT·Δa (ABUNDANCE) + a_NAT·ΔM (WIRING) + Δa·ΔM (INTERACTION / co-acquisition).
    M = the CANONICAL bagged family weight per state (states.canonical_M) so ΔM is stable; M·a = contribution.
    Surfaces per-state IQR of arm + gene so a low-variance (undetectable / range-restricted) state is visible."""
    Yt, Xt, Ct, wt = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    an = ST.assemble_state(gene, "11", family=False)
    if an is None:
        print(f"{gene}: no NAT assembly"); return pd.DataFrame()
    Yn, Xn, Cn, wn = an
    mn = Yn.notna() & Cn.notna().all(axis=1)                      # assemble_state doesn't NaN-clean; do it here
    Yn, Xn, Cn = Yn[mn], Xn.loc[mn].fillna(0.0), Cn.loc[mn]
    arms = [a for a in Xt.columns if a in Xn.columns]
    Xt_a, Xn_a = Xt[arms].fillna(0.0), Xn[arms].fillna(0.0)
    Mt = ST.canonical_M(gene, "01", arm_level=True).reindex(arms).fillna(0)   # CANONICAL bagged family M per state
    Mn = ST.canonical_M(gene, "11", arm_level=True).reindex(arms).fillna(0)   # → stable ΔM (was single-fit lasso)
    at, anm = Xt_a.mean(), Xn_a.mean()
    da, dM = at - anm, Mt - Mn

    def iqr(s):
        s = s.dropna(); return round(float(s.quantile(.75) - s.quantile(.25)), 2) if len(s) else np.nan
    df = pd.DataFrame({
        "M_NAT": Mn.round(3), "M_TUM": Mt.round(3), "dabund": da.round(2), "dM": dM.round(3),
        "term_ABUND": (Mn * da).round(3), "term_WIRING": (anm * dM).round(3),
        "term_INTERACT": (da * dM).round(3),
        "arm_iqr_NAT": [iqr(Xn[a]) for a in arms], "arm_iqr_TUM": [iqr(Xt[a]) for a in arms],
    })
    df["dContrib"] = (df["term_ABUND"] + df["term_WIRING"] + df["term_INTERACT"]).round(3)
    gN, gT = iqr(Yn), iqr(Yt)
    df = df[(df[["M_NAT", "M_TUM"]].abs().sum(axis=1) > 0)].sort_values("dContrib")
    print(f"\n=== {gene} — NAT→tumour ABUNDANCE / WIRING / INTERACTION decomposition (gene IQR NAT={gN} TUM={gT}) ===")
    print(df[["M_NAT", "M_TUM", "dabund", "dM", "term_ABUND", "term_WIRING", "term_INTERACT", "dContrib",
              "arm_iqr_NAT", "arm_iqr_TUM"]].to_string())
    tot = df[["term_ABUND", "term_WIRING", "term_INTERACT"]].sum()
    print(f"  Σ terms: ABUNDANCE {tot['term_ABUND']:+.2f} | WIRING {tot['term_WIRING']:+.2f} | "
          f"INTERACTION {tot['term_INTERACT']:+.2f}  (which axis drives the gene's NAT→tumour pressure shift)")
    return df


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    fn = (stable_wiring if "--stable" in sys.argv else state_wiring if "--wiring" in sys.argv
          else decompose if "--decompose" in sys.argv else gene_card)
    for g in (args or ["GATA3", "PTEN"]):
        fn(g)
