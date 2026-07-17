"""Per-edge (and per-gene) attribution card — the taxonomy rendered per edge on the confounder-honest,
cell-intrinsic foundation (docs/ATTRIBUTION_CONTEXT_AXIS.md). One row per regulator of a gene, joining the
four layers we already compute:

  1. REGIME (range stats) — arm median RPM, %>functional-floor (RPM≥10), IQR, SPIKER flag (low median +
     high IQR = subset-driven, e.g. miR-135b 28%>floor); target IQR (a strong coupling on a low-IQR gene =
     the miRNA explains most of its variance).
  2. BUDGET share (E7/G4) — the arm's rank + share of the gene's pressure, GTEx→NAT→tumour Δ (states.budget_shift).
  3. COMPOSITION tag — deconv retention: cell-intrinsic (≥0.7) / partial / composition-explained (<0.4).
  4. SHIFT-CLASS — a learned-model class (level × realization × composition × detectability), compatible in
     spirit with the precursor `mirna_state_class.joint_edge_class` {acquired_realized/lost/stable/…}.

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
    """Full H→NAT→T annotation (level × 3-state realization × composition × detectability), compatible with
    the precursor `mirna_state_class.joint_edge_class`."""
    ret, spk = row["retention"], row["spiker"]
    rep = lambda c: (c == c) and c <= -0.1                       # "repressor in this state"
    h, n, t = rep(row["coupling_hly"]), rep(row["coupling_nat"]), rep(row["coupling_tum"])
    if not t:
        return "undetectable" if (spk or (row["arm_pct_floor"] == row["arm_pct_floor"] and row["arm_pct_floor"] < 30)) \
            else ("lost" if (h or n) else "uncoupled")
    if ret == ret and ret < 0.4:
        return "composition_explained"
    acq = (row.get("arm_lfc_HLY_TUM_QN", 0) or 0) > 0.3            # acquired elevation over healthy (dHT logFC)
    if h and n:
        return "constitutive"                                   # present+realized in healthy already
    if (not h) and n:
        return "field_established_realized"                     # established in NAT (field effect, dHN)
    if (not h) and (not n):
        return "acquired_realized" if acq else "tumour_realized"  # tumour-specific (dNT)
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
    gr_tum = Xt2.mean(axis=1).rank(pct=True) * 100
    gr_nat = Xn.mean(axis=1).rank(pct=True) * 100
    try:
        from mirna_hallmark.learned import state as STA
        gmean = STA._gtex_mirna().mean(axis=1)                 # raw GTEx per-arm mean (log2(TPM+1))
        gr_hly = gmean.rank(pct=True) * 100
    except Exception:
        gmean, gr_hly = pd.Series(dtype=float), pd.Series(dtype=float)

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
                     "retention": round(ret, 2) if ret == ret else np.nan,
                     "gene_iqr": round(gene_iqr, 2)})
    card = pd.DataFrame(rows)
    card["shift_class"] = card.apply(_shift_class, axis=1)
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
