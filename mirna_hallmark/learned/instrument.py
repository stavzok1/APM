"""CN-locus instrument — the causal test (Design §5 Bar 4, §Decision J). Bar-4 evaluation.

miRNA-locus copy number is a *natural genetic perturbation* of arm dose that is (largely) not caused by
the target's expression — so it breaks the confounding loop observational coupling cannot. Per edge:

    first stage :  x_m  = γ·CN_locus(m) + C·δ + u      (does locus CN move the arm? F > 10 ⇒ usable)
    reduced form:  y_g  = π·CN_locus(m) + C·β + e      (does the genetic dose move the target? π < 0 ⇒ repression)
    conditioned :  partial-Spearman(x_m, y_g | C)      (the observational coupling, for comparison)

A well-instrumented edge whose reduced-form coupling is **negative and sign-concordant** with the
observational coupling is causal evidence, not covariation. Weak instruments (low F) are reported as
INCONCLUSIVE, not failed (Design §J). Pleiotropy caveat: locus CN can hit neighbours — restrict to focal
loci / check the arm locus does not overlap the target gene. Data: mirna_locus_cnv (ASCAT3 → locus CN).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.linear_model import LinearRegression

from mirna_hallmark import data_loaders as D
from mirna_hallmark.coupling_inference import benjamini_yekutieli
from mirna_hallmark.learned import data as LD

_CACHE: dict = {}

_ADMIT_SOFT_C = 10.0   # CN-admission soft null-corrector: A_cn·(1−c/F)₊ (C4 default; cleaner than the hard F>10 cliff — §C, `soft_fweight_sweep`)

_ENTITY_CNV = "mirna_hallmark/output/mirna_locus_cnv/tables/sample_entity_cnv.tsv.gz"
_PARALOG_MAP = "mirna_hallmark/output/mirna_locus_cnv/maps/mirna_cnv_mimat_paralog_map.tsv"


def _entity_long() -> pd.DataFrame:
    """Cached read of the full participant × entity CNV long table (locus + arm + family)."""
    if "long" not in _CACHE:
        from pathlib import Path
        _CACHE["long"] = pd.read_csv(Path(_ENTITY_CNV), sep="\t", low_memory=False)
    return _CACHE["long"]


def arm_cn() -> pd.DataFrame:
    """participant × arm copy number (ASCAT3). **LEGACY** paralog-AGGREGATED CN (w_norm composite of an
    arm's loci) — kept for comparison only; the instrument now sources per-locus CN (`locus_cn` +
    `_arm_focal_locus`) so paralog arms stop contaminating (CN_INSTRUMENT.md §9, gap 1)."""
    if "cn" not in _CACHE:
        from mirna_hallmark import mirna_locus_cnv as MC
        cn = MC._participant_arm_cn(_entity_long())
        _CACHE["cn"] = cn.pivot_table(index="participant", columns="entity_label",
                                      values="copy_number", aggfunc="first")
    return _CACHE["cn"]


def locus_cn() -> pd.DataFrame:
    """participant × pre-miRNA hairpin locus (`MI*`) copy number — the per-locus ASCAT3 CN, un-aggregated.
    This is the instrument's exogenous source (CN_INSTRUMENT.md §6/§9 gap 1): CN is observed per genomic
    segment, so each locus is a separable genetic perturbation; a mature arm's aggregated CN is not."""
    if "locus_cn" not in _CACHE:
        lg = _entity_long()
        loc = lg[lg["entity_level"] == "locus"].copy()
        loc["copy_number"] = pd.to_numeric(loc["copy_number"], errors="coerce")
        _CACHE["locus_cn"] = loc.pivot_table(index="participant", columns="entity_id",
                                             values="copy_number", aggfunc="first")
    return _CACHE["locus_cn"]


def arm_loci_map() -> dict:
    """arm label → its contributing hairpin loci as a DataFrame [locus_id, chrom, midpoint, w_norm],
    sorted by w_norm (mature-dose contribution) descending. Source: the MIMAT paralog map (authoritative —
    keeps the full multi-locus structure + per-locus dose weights); arms absent from it (the single-locus
    tail the map doesn't enumerate) fall back to the arm entity row's own `pre_gene_id`."""
    if "arm_loci" not in _CACHE:
        from pathlib import Path
        pm = pd.read_csv(Path(_PARALOG_MAP), sep="\t")
        arm_rows = (_entity_long().loc[lambda d: d["entity_level"] == "arm",
                    ["mirbase_mature_id", "mature_accession", "pre_gene_id"]]
                    .drop_duplicates("mirbase_mature_id"))
        mimat2arm = dict(zip(arm_rows["mature_accession"], arm_rows["mirbase_mature_id"]))
        m: dict = {}
        for mimat, g in pm.groupby("mimat"):
            arm = mimat2arm.get(mimat)
            if arm is None:
                continue
            m[arm] = (g[["locus_id", "chrom", "midpoint", "w_norm"]]
                      .sort_values("w_norm", ascending=False).reset_index(drop=True))
        for _, r in arm_rows.iterrows():                                   # single-locus fallback
            arm = r["mirbase_mature_id"]
            if arm in m or pd.isna(r["pre_gene_id"]):
                continue
            m[arm] = pd.DataFrame({"locus_id": [r["pre_gene_id"]], "chrom": [np.nan],
                                   "midpoint": [np.nan], "w_norm": [1.0]})
        _CACHE["arm_loci"] = m
    return _CACHE["arm_loci"]


def _arm_focal_locus(arm: str) -> dict:
    """(memoised per arm — pure function of the paralog map / locus_origin.) The arm's focal instrument locus for
    the *single-instrument* rank test (`edge_instrument`). Prefers an
    **active_source** hairpin (`locus_origin`: the arm is actually transcribed there) over the mere dominant-
    w_norm one — so a high-w_norm but SILENT locus is never picked as the instrument; falls back to max-w_norm
    when origin is unavailable. `n_active_source` is exposed because with ≥2 active sources a single focal is a
    simplification: the proper estimate is `multi_iv` over ALL of them (single-locus arm ⇒ CN ≡ aggregated exactly)."""
    memo = _CACHE.setdefault("focal_locus", {})
    if arm in memo:
        return memo[arm]
    g = arm_loci_map().get(arm)
    if g is None or not len(g):
        return memo.setdefault(arm, {"locus_id": None, "n_loci": 0, "multilocus": False, "n_active_source": 0})
    n_loci = len(g)
    focal, n_active = str(g.iloc[0]["locus_id"]), 0                                 # default: dominant w_norm
    if n_loci > 1:                                                                  # single-locus is trivially focal
        origin = locus_origin(arm)
        if not origin.empty and "status" in origin.columns:
            active = origin.loc[origin["status"] == "active_source", "locus_id"].tolist()
            n_active = len(active)
            if active:
                wmap = dict(zip(g["locus_id"], g["w_norm"]))
                focal = str(max(active, key=lambda L: wmap.get(L, 0.0)))            # dominant AMONG active sources
    return memo.setdefault(arm, {"locus_id": focal, "n_loci": n_loci, "multilocus": bool(n_loci > 1),
                                 "n_active_source": n_active})


def _arm_active_loci(arm: str) -> list:
    """The arm's ACTIVE-SOURCE instrument loci (w_norm-desc, focal first) — the COMPLETE valid instrument set for the
    arm. Silent-source loci are NOT valid instruments (relevance fails: their CN doesn't move the mature dose), so they
    are excluded — hence for a SINGLE-active-source arm this returns just the focal locus and the focal simplification is
    **EXACT, not an approximation** (783/785 HE arms). With ≥2 active sources it returns all of them (the `multi_iv`-
    complete set), so `exclusion` can instrument X_fam with each separately + run a within-arm over-ID."""
    memo = _CACHE.setdefault("active_loci", {})
    if arm in memo:
        return memo[arm]
    foc = _arm_focal_locus(arm)
    out = [foc["locus_id"]] if foc.get("locus_id") else []
    if foc.get("n_active_source", 0) >= 2:
        try:
            o = locus_origin(arm); g = arm_loci_map().get(arm)
            wmap = dict(zip(g["locus_id"].astype(str), g["w_norm"])) if g is not None else {}
            act = sorted({str(L) for L in o.loc[o["status"] == "active_source", "locus_id"]},
                         key=lambda L: -wmap.get(L, 0.0))               # w_norm desc → focal first
            if act:
                out = act
        except Exception:
            pass
    return memo.setdefault(arm, out)


def _resid(v, Cmat):
    return v - LinearRegression().fit(Cmat, v).predict(Cmat)


def _first_stage_F(x, instr, Cmat) -> float:
    """F for adding the instrument to x ~ C (first-stage strength)."""
    n = len(x)
    Xf = np.column_stack([np.ones(n), Cmat, instr]); Xr = np.column_stack([np.ones(n), Cmat])
    br, *_ = np.linalg.lstsq(Xr, x, rcond=None); bf, *_ = np.linalg.lstsq(Xf, x, rcond=None)
    rss_r = float(((x - Xr @ br) ** 2).sum()); rss_f = float(((x - Xf @ bf) ** 2).sum())
    df2 = n - Xf.shape[1]
    return float(((rss_r - rss_f) / 1) / (rss_f / df2)) if rss_f > 0 and df2 > 0 else np.nan


def edge_instrument(arm: str, gene: str, *, locus: str | None = None) -> dict:
    """CN instrument for one edge, sourced from **per-locus** CN (gap 1). The instrument is the arm's focal
    hairpin-locus CN (or an explicit `locus`, the hook multi-IV/gap 2 calls per-locus), not the paralog-
    aggregated arm CN — so a multi-locus arm's chr-mixed dose no longer contaminates the exclusion."""
    Y, X, C, _ = LD.assemble_gene(gene, he_only=False)
    if arm not in X.columns:
        raise KeyError(f"{arm} not a candidate regulator of {gene}")
    foc = _arm_focal_locus(arm)
    lid = locus or foc["locus_id"]
    if lid is None:
        raise KeyError(f"{arm} has no locus CN")
    lcn = locus_cn()
    if lid not in lcn.columns:
        raise KeyError(f"locus {lid} ({arm}) has no CN")
    parts = X.index.intersection(lcn.index)
    x = X.loc[parts, arm].to_numpy(float)
    y = Y.loc[parts].to_numpy(float)
    Cmat = C.loc[parts].to_numpy(float)
    z = pd.to_numeric(lcn.loc[parts, lid], errors="coerce").to_numpy(float)
    ok = np.isfinite(z)
    x, y, Cmat, z = x[ok], y[ok], Cmat[ok], z[ok]
    F = _first_stage_F(x, z, Cmat)
    rho_fs = spearmanr(_resid(z, Cmat), _resid(x, Cmat)).correlation      # CN → arm expr (expect > 0)
    _rf = spearmanr(_resid(z, Cmat), _resid(y, Cmat))                     # CN → target  (expect < 0)
    rho_rf = _rf.correlation
    p_rf = (_rf.pvalue / 2) if (rho_rf == rho_rf and rho_rf < 0) else (1 - _rf.pvalue / 2)   # one-sided rho_reduced<0
    rho_obs = spearmanr(_resid(x, Cmat), _resid(y, Cmat)).correlation     # observational coupling
    return {"arm": arm, "gene": gene, "n": int(len(x)), "locus_id": lid,
            "n_loci": foc["n_loci"], "n_active_source": foc.get("n_active_source", 0),
            "multilocus": foc["multilocus"], "first_stage_F": F,
            "rho_firststage": rho_fs, "rho_reduced": rho_rf, "p_reduced": float(p_rf), "rho_observational": rho_obs,
            "usable": bool(F is not None and F > 10),
            "causal_concordant": bool(rho_rf < 0 and rho_obs < 0 and F is not None and F > 10)}


def _mirna_coords():
    """arm → (chrom, midpoint, precursor). For the pleiotropy check."""
    if "coords" not in _CACHE:
        from mirna_hallmark import config as C
        m = pd.read_csv(C.MIRNA_MATURE_LOCI)
        m = m.dropna(subset=["mirbase_mature_id"]).drop_duplicates("mirbase_mature_id")
        m["mid"] = (m["start"] + m["end"]) / 2
        _CACHE["coords"] = m.set_index("mirbase_mature_id")[["chrom", "mid", "pre_gene_name"]]
    return _CACHE["coords"]


def pleiotropy(arm: str, gene: str, *, window: float = 1e6, edges=None) -> dict:
    """Horizontal-pleiotropy check for the CN instrument (Design §5 caveat). The locus CN of `arm` also
    changes any OTHER miRNA at the same locus — so if a co-located miRNA ALSO targets `gene`, the
    reduced-form (CN→target) is the *locus/cluster* effect, not `arm`'s alone. Returns the co-located
    miRNAs that target `gene` (the pleiotropy sources) + same-precursor (polycistron) members."""
    co = _mirna_coords()
    edges = edges if edges is not None else D.high_evidence_edges()
    reg = set(edges.loc[edges["gene"] == gene, "miRNA"])
    if arm not in co.index:
        return {"arm": arm, "gene": gene, "coords": None}
    c = co.loc[arm]
    near = co[(co["chrom"] == c["chrom"]) & ((co["mid"] - c["mid"]).abs() <= window)]
    near_arms = set(near.index) - {arm}
    polycistron = sorted(set(co[(co["pre_gene_name"] == c["pre_gene_name"])].index) - {arm})
    coloc_targeting = sorted(near_arms & reg)                   # co-located AND regulate the same gene
    return {"arm": arm, "gene": gene, "chrom": c["chrom"], "n_colocated": len(near_arms),
            "coloc_targeting_gene": coloc_targeting, "polycistron_members": polycistron,
            "pleiotropic": len(coloc_targeting) > 0}


def cluster_attribution(arm: str, gene: str, *, window: float = 1e6) -> dict:
    """Rescue the CN instrument's exclusion restriction at the cluster level (user 2026-07-05). The locus CN
    moves the whole cluster, so if a co-located miRNA ALSO targets `gene`, the reduced-form is the cluster's
    effect, not `arm`'s. Two outcomes:
      • CLUSTER-CLEAN — no co-located co-targeter ⇒ exclusion holds, the CN edge is `arm`-specific.
      • otherwise ATTRIBUTE within the cluster by EXPRESSION: is `arm` the UNIQUE member of the co-locating
        co-targeter set whose anti-corr with the target SURVIVES conditioning on the others (+C)? If yes, the
        cluster's causal effect is carried by `arm`. Bounded-circular: CN supplies the exogeneity; expression
        only resolves WHICH cluster member carries it (same logic as the seed-family driver nomination, §4b)."""
    p = pleiotropy(arm, gene, window=window)
    Y, X, C, _ = LD.assemble_gene(gene, he_only=False)
    co_t = [m for m in (p.get("coloc_targeting_gene") or []) if m in X.columns]
    Cm = C.to_numpy(float); yv = Y.to_numpy(float)
    S = [arm] + co_t

    def _cond(m, others):
        Z = np.column_stack([Cm] + ([X[others].to_numpy(float)] if others else []))
        xr = _resid(X[m].to_numpy(float), Z)
        if float(np.std(xr)) < 1e-6:
            return np.nan
        return float(spearmanr(xr, _resid(yv, Z)).correlation)

    part = {m: _cond(m, [o for o in S if o != m]) for m in S if m in X.columns}
    our = part.get(arm, np.nan)
    others_surv = [m for m in co_t if part.get(m) == part.get(m) and part[m] < -0.1]
    unique = bool(our == our and our < -0.1 and not others_surv)
    return {"arm": arm, "gene": gene, "n_cotargeters": len(co_t),
            "cluster_clean": len(co_t) == 0,
            "arm_cluster_partial": round(our, 3) if our == our else np.nan,
            "cotargeters": ",".join(co_t) if co_t else "",
            "cotargeter_partials": ";".join(f"{m.replace('hsa-','')}={part[m]:+.2f}" if part.get(m) == part.get(m)
                                            else f"{m.replace('hsa-','')}=collinear" for m in co_t),
            "arm_unique_anticorr": unique,
            "cluster_attributable": bool(len(co_t) == 0 or unique)}


def _genomic_clusters(arms, window: float = 5e4) -> dict:
    """arm → genomic-cluster id (same chrom, single-linkage within `window` bp). The CLUSTER is what a
    copy-number event moves together — the correct aggregation unit for the CN instrument (a seed family
    can span several clusters on different loci, each with independent CN)."""
    co = _mirna_coords()
    a = [x for x in arms if x in co.index]
    d = co.loc[a].sort_values(["chrom", "mid"])
    cid, k, prev = {}, 0, None
    for arm, row in d.iterrows():
        if prev is None or row["chrom"] != prev[0] or abs(row["mid"] - prev[1]) > window:
            k += 1
        cid[arm] = k
        prev = (row["chrom"], row["mid"])
    return cid


def _shapley_r2(Xr: np.ndarray, yr: np.ndarray) -> np.ndarray:
    """Exact Shapley split of R²(yr ~ Xr) among predictors — variation-based marginal contributions."""
    import itertools
    from math import factorial
    k = Xr.shape[1]
    yss = float((yr ** 2).sum()) or 1.0

    def r2(S):
        if not S:
            return 0.0
        Xs = Xr[:, list(S)]
        b, *_ = np.linalg.lstsq(Xs, yr, rcond=None)
        return 1 - float(((yr - Xs @ b) ** 2).sum()) / yss

    phi = np.zeros(k)
    for i in range(k):
        others = [c for c in range(k) if c != i]
        for r in range(len(others) + 1):
            w = factorial(r) * factorial(k - r - 1) / factorial(k)
            for S in itertools.combinations(others, r):
                phi[i] += w * (r2(set(S) | {i}) - r2(set(S)))
    return np.clip(phi, 0.0, None)


def family_multi_iv(gene: str, mem, foc=None, *, min_n: int = 60, assembled=None) -> dict:
    """Family-level multi-IV — the §7/§8 synthesis (replaces the mean-CN collapse). Instrument the family
    aggregate abundance `X_fam = log2(Σ_m 2^x_m)` by EACH member's focal-locus CN (a separate instrument),
    residualised on C. Yields, per member m:
        gamma_m  — first stage (CN_m → X_fam): member m's DOSE delivery to the family pool;
        pi_m     — reduced form (CN_m → target): member m's genetic-dose effect on the TARGET (the causal leg);
        wald_m   — pi_m/gamma_m: target effect per unit dose (same-seed ⇒ should be a shared W).
    Group: 2SLS `beta` (family effect per unit X_fam), joint first-stage F, Hansen-J over-ID (tests the wald_m
    equal = same-seed interchangeability). CAUSAL attribution to member m = pi_m (share = gamma_m share, since W
    is shared). **CAVEAT:** pi_m is confounded if OTHER co-located seed families also target `gene` (the cluster CN
    moves them too) — that exclusion is the between-family Kind-B step (parked); over-ID/pleiotropy partly flags it."""
    from scipy.stats import chi2
    Y, X, C, _ = assembled if assembled is not None else LD.assemble_gene(gene, he_only=False)
    if foc is None:
        foc = {m: _arm_focal_locus(m)["locus_id"] for m in mem}
    lcn = locus_cn()
    mem = [m for m in mem if m in X.columns and foc.get(m) in lcn.columns]
    if len(mem) < 2:
        return {}
    parts = X.index.intersection(lcn.index)
    Xm = X.loc[parts, mem].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    CNm = np.column_stack([pd.to_numeric(lcn.loc[parts, foc[m]], errors="coerce").to_numpy(float) for m in mem])
    y, Cm = Y.loc[parts].to_numpy(float), C.loc[parts].to_numpy(float)
    ok = np.isfinite(Xm).all(1) & np.isfinite(CNm).all(1) & np.isfinite(y) & np.isfinite(Cm).all(1)
    Xm, CNm, y, Cm = Xm[ok], CNm[ok], y[ok], Cm[ok]
    n, k = len(y), len(mem)
    if n < min_n:
        return {}
    xr = _resid(np.log2((2.0 ** Xm).sum(1)), Cm)                                    # X_fam residualised on C
    yr = _resid(y, Cm)
    Zr = np.column_stack([_resid(CNm[:, j], Cm) for j in range(k)])
    tss = float((xr ** 2).sum())
    solo = {}
    for j, m in enumerate(mem):
        zj = Zr[:, j]
        gj, pj = float((zj @ xr) / (zj @ zj)), float((zj @ yr) / (zj @ zj))         # solo first stage / reduced form
        sf = float(((xr - zj * gj) ** 2).sum())
        Fj = float(((tss - sf) / 1) / (sf / (n - 1))) if sf > 0 else np.nan
        solo[m] = {"F": Fj, "gamma": gj, "pi": pj, "wald": (pj / gj) if abs(gj) > 1e-12 else np.nan,
                   "sd_cn": float(np.std(zj))}
    contrib = {m: abs(solo[m]["gamma"]) * solo[m]["sd_cn"] for m in mem}            # dose share (cn_copy)
    ctot = sum(contrib.values()) or 1.0
    cn_copy = {m: contrib[m] / ctot for m in mem}                                   # THE attribution share (§8)
    strong = [m for m in mem if solo[m]["F"] == solo[m]["F"] and solo[m]["F"] > 10]
    ks = len(strong)
    beta = F = overid_J = overid_p = np.nan
    n_indep = 0
    if ks:
        Zs = Zr[:, [mem.index(m) for m in strong]]
        g = np.linalg.lstsq(Zs, xr, rcond=None)[0]
        xhat = Zs @ g
        ssr = float(((xr - xhat) ** 2).sum())
        F = float(((tss - ssr) / ks) / (ssr / (n - ks))) if ssr > 0 and n > ks else np.nan
        beta = float((xhat @ yr) / (xhat @ xhat)) if (xhat @ xhat) > 0 else np.nan  # family effect per unit X_fam
        # over-ID needs INDEPENDENT strong instruments — same-CN members are one direction; keep a full-rank subset
        indep, kept = [], []
        for m in sorted(strong, key=lambda a: -solo[a]["F"]):
            col = Zr[:, mem.index(m)]
            if all(abs(np.corrcoef(col, Zr[:, mem.index(o)])[0, 1]) < 0.999 for o in kept):
                indep.append(mem.index(m)); kept.append(m)
        n_indep = len(kept)
        if n_indep > 1:                                                             # Hansen-J over-ID (same-seed check)
            Zi = Zr[:, indep]
            e = yr - beta * xr
            gm = Zi.T @ e / n
            S = (Zi * (e ** 2)[:, None]).T @ Zi / n
            try:
                overid_J = float(n * gm @ np.linalg.solve(S, gm))
                overid_p = float(chi2.sf(overid_J, n_indep - 1))
            except np.linalg.LinAlgError:
                pass
    return {"gene": gene, "n": n, "members": mem, "k": k, "n_strong": ks, "n_indep_cn": n_indep,
            "first_stage_F": F, "beta_family": beta, "overid_J": overid_J, "overid_p": overid_p,
            "over_id_powered": bool(n_indep > 1),
            "recurse_flag": bool(n_indep > 1 and overid_p == overid_p and overid_p < 0.05),
            "solo": solo, "cn_copy": cn_copy}


def family_causal_attribution(gene: str, mem, foc=None, *, min_n: int = 60, min_F: float = 10.0) -> pd.DataFrame:
    """**Ring-1 (CN_INSTRUMENT §7): the hierarchical within-family causal attribution.** Composes the two things
    only CN+expression together can do:
        portion(arm) = p_S(its segment)  ×  functional-dose share WITHIN that segment
    where **p_S** = the segment's EXOGENOUS CN-causal reduced-form share (CN separates independent-CN segments),
    and **within a segment** (CN-blind, same-seed members) the split is **dose-delivery** (§8, L2): K_D is CONSTANT
    within a family (same seed → same affinity) so it can't split members — use the MEASURED arm-resolved signal
    (chimeric binding, `l2_src=chimeric`) where present, else abundance (loading resolved-INERT 2026-07-10).
    Reports the isolation regime (n_segments, per-segment first-stage F). CAVEAT (§ exclusion): p_S is confounded by
    co-located OTHER-seed families that also hit `gene` — the parked between-family (Kind-B) step."""
    from mirna_hallmark.learned import chimeric_evidence as CE
    Y, X, C, _ = LD.assemble_gene(gene, he_only=False)
    if foc is None:
        foc = {m: _arm_focal_locus(m)["locus_id"] for m in mem}
    lcn = locus_cn()
    mem = [m for m in mem if m in X.columns and foc.get(m) in lcn.columns]
    if len(mem) < 2:
        return pd.DataFrame()
    parts = X.index.intersection(lcn.index)
    Cm = C.loc[parts].to_numpy(float)
    yv = Y.loc[parts].to_numpy(float)
    CN = {m: pd.to_numeric(lcn.loc[parts, foc[m]], errors="coerce").to_numpy(float) for m in mem}
    segs = []                                                                       # group members by CN collinearity (|ρ|>0.999)
    for m in mem:
        for s in segs:
            c = np.corrcoef(np.nan_to_num(CN[m]), np.nan_to_num(CN[s[0]]))[0, 1]
            if abs(c) > 0.999:
                s.append(m); break
        else:
            segs.append([m])
    import warnings
    seg_info = []
    for s in segs:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)                        # all-NaN CN row → NaN (masked below)
            zc = np.nanmean(np.column_stack([CN[m] for m in s]), axis=1)
        xf = X.loc[parts, s].mean(axis=1).to_numpy(float)
        ok = np.isfinite(zc) & np.isfinite(yv) & np.isfinite(xf)
        if ok.sum() < min_n:
            seg_info.append((s, np.nan, np.nan)); continue
        F = _first_stage_F(xf[ok], zc[ok], Cm[ok])
        pi = float(spearmanr(_resid(zc[ok], Cm[ok]), _resid(yv[ok], Cm[ok])).correlation)   # reduced form (CN_seg→target)
        seg_info.append((s, pi, F))
    tot = sum(abs(pi) for _, pi, F in seg_info if pi == pi and F == F and F > min_F) or 1.0
    rows = []
    for si, (s, pi, F) in enumerate(seg_info):
        pS = (abs(pi) / tot) if (pi == pi and F == F and F > min_F) else 0.0        # exogenous segment share
        # L2 within-segment dose-delivery: K_D is CONSTANT within a family (same seed → same affinity), so it can't
        # split members — use the MEASURED arm-resolved signal (chimeric binding) where present, else abundance.
        chim = CE.evidence(gene, s)
        if float(chim.sum()) > 0:
            w = {m: float(chim.get(m, 0.0)) for m in s}; src = "chimeric"
        else:
            w = {m: max(float(np.median(2.0 ** X.loc[parts, m].to_numpy(float))), 0.0) for m in s}; src = "abund"
        wt = sum(w.values()) or 1.0
        for m in s:
            rows.append({"gene": gene, "arm": m, "segment": si, "seg_pi": round(pi, 3) if pi == pi else np.nan,
                         "seg_F": round(F, 1) if F == F else np.nan, "seg_pS": round(pS, 3), "l2_src": src,
                         "within_dose_share": round(w[m] / wt, 3), "portion": round(pS * w[m] / wt, 4)})
    # FAMILY-LEVEL: absolute magnitude (2SLS β) + CN-vs-expression fraction (how much of the family's antagonism is
    # CN-driven/exogenous vs total observational — the family-level `cn_attribution`, MH-43 idea at the family unit).
    fr = family_multi_iv(gene, mem, foc)
    beta = float(fr.get("beta_family", np.nan)) if fr else np.nan
    Xm = X.loc[parts, mem].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    CNm = np.column_stack([CN[m] for m in mem])
    okf = np.isfinite(Xm).all(1) & np.isfinite(CNm).all(1) & np.isfinite(yv) & np.isfinite(Cm).all(1)
    obs = cnfrac = np.nan
    if okf.sum() >= min_n:
        xfam = np.log2((2.0 ** Xm[okf]).sum(1)); Cf = Cm[okf]; yr = _resid(yv[okf], Cf)
        obs = float(spearmanr(_resid(xfam, Cf), yr).correlation)                    # observational family antagonism
        Zc = np.column_stack([_resid(CNm[okf][:, j], Cf) for j in range(len(mem))])
        b = np.linalg.lstsq(Zc, _resid(xfam, Cf), rcond=None)[0]                    # CN-driven part of X_fam
        rcn = float(spearmanr(Zc @ b, yr).correlation)
        cnfrac = round(rcn / obs, 3) if obs and abs(obs) > 1e-6 else np.nan
    df = pd.DataFrame(rows)
    if len(df):
        df["portion_abs"] = (df["portion"] * beta).round(4) if beta == beta else np.nan   # = share × family 2SLS β
        df["family_beta"] = round(beta, 3) if beta == beta else np.nan
        df["family_obs_coupling"] = round(obs, 3) if obs == obs else np.nan               # total observational antagonism
        df["cn_driven_fraction"] = cnfrac                                                  # CN-driven / observational
    return df.sort_values("portion", ascending=False) if len(df) else df


def run_family_causal_attribution(genes=None, *, min_members: int = 2, progress: int = 100,
                                  out: str = "mirna_hallmark/output/learned/instrument/family_causal_attribution.tsv") -> pd.DataFrame:
    """Batch **Ring-1**: `family_causal_attribution` over every seed family (≥`min_members` regulators) of each gene
    in `genes` (default = the HE gene universe). Persists per-(gene, family, arm) portions + family-level magnitude/
    CN-driven fraction. CLI: `python -m mirna_hallmark.learned.instrument --ring1`."""
    from pathlib import Path
    from mirna_hallmark.learned import families as FAM
    if genes is None:
        genes = sorted(D.high_evidence_edges()["gene"].dropna().astype(str).unique())
    rows = []
    for gi, g in enumerate(genes):
        if progress and gi % progress == 0:
            print(f"[ring1] {gi}/{len(genes)} genes, {len(rows)} family-tables", flush=True)
        try:
            _, X, _, _ = LD.assemble_gene(g, he_only=False)
        except Exception:
            continue
        fo = FAM.family_of(list(X.columns))
        fams: dict = {}
        for a in X.columns:
            f = fo.get(a)
            if f:
                fams.setdefault(f, []).append(a)
        for f, mem in fams.items():
            if len(mem) < min_members:
                continue
            try:
                df = family_causal_attribution(g, mem)
            except Exception:
                continue
            if len(df):
                df.insert(1, "family", f)
                rows.append(df)
    if not rows:
        return pd.DataFrame()
    out_df = pd.concat(rows, ignore_index=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out, sep="\t", index=False)
    caus = out_df[(out_df["seg_pS"] > 0) & (out_df["portion"] > 0)]
    print(f"[ring1] {len(out_df)} arm-rows | {out_df['gene'].nunique()} genes × {out_df['family'].nunique()} families | "
          f"{caus['gene'].nunique()} genes with ≥1 causal segment → {out}")
    return out_df


_GENCODE_CSV = f"{__import__('os').environ.get('APM_ROOT', '/sci/labs/michall/stavzok/APM')}/data/gencode.v49.annotation.gtf.csv"
_GENE_CN_CACHE = "mirna_hallmark/output/matrices/cnv_target_genes.tsv.gz"
_PROLIF_GENES = ["MKI67", "PCNA", "TOP2A", "CCNB1", "CCNB2", "AURKA", "AURKB", "BUB1", "CDK1", "CDC20",
                 "CENPF", "MCM2", "MCM6", "PLK1", "TYMS", "RRM2", "E2F1"]


def _coding_coords() -> pd.DataFrame:
    """protein-coding gene → (chrom, mid) from GENCODE (cached). For the coding-pleiotropy co-location scan."""
    if "coding_coords" not in _CACHE:
        g = pd.read_csv(_GENCODE_CSV, usecols=["seqname", "start", "end", "feature", "gene_type", "gene_name"])
        g = g[(g["feature"] == "gene") & (g["gene_type"] == "protein_coding")]
        g["mid"] = (g["start"] + g["end"]) / 2
        _CACHE["coding_coords"] = g.groupby("gene_name").agg(chrom=("seqname", "first"), mid=("mid", "mean")).reset_index()
    return _CACHE["coding_coords"]


_GENE_CN_GENOMEWIDE = "mirna_hallmark/output/matrices/gene_cn_genomewide.tsv.gz"


def _gene_cn() -> pd.DataFrame:
    """gene × participant copy number for the co-amplification filter. Prefers the **genome-wide** coding-gene CN
    (20k genes, ASCAT-segment × GENCODE intersection, `build_genomewide_gene_cn.py`, MH-99 TODO #1); falls back to the
    4363-gene Hallmark cache if the genome-wide matrix isn't built."""
    if "gene_cn" not in _CACHE:
        import os
        path = _GENE_CN_GENOMEWIDE if os.path.exists(_GENE_CN_GENOMEWIDE) else _GENE_CN_CACHE
        _CACHE["gene_cn"] = pd.read_csv(path, sep="\t", index_col=0)
    return _CACHE["gene_cn"]


def _as_np(df: pd.DataFrame, key: str):
    """(values, row_name→i, col_name→j) for a gene×participant frame — CACHED. Lets `coding_pleiotropy` pull gene/participant
    blocks by NUMPY positional indexing instead of pandas `.loc[list, list]` fancy-reindex (which was 77% of the scan)."""
    if key not in _CACHE:
        _CACHE[key] = (df.to_numpy(float),
                       {n: i for i, n in enumerate(df.index)},
                       {c: j for j, c in enumerate(df.columns)})
    return _CACHE[key]


def coding_pleiotropy(gene: str, mem, *, window: float = 3e6, coamp: float = 0.4, min_corr: float = 0.3,
                      surv: float = 0.15, dmin: float = 0.05, topk: int = 3, assembled=None,
                      host_genes: dict | None = None, scan: bool = True) -> pd.DataFrame:
    """SAFE coding-gene + HOST pleiotropy gate (CHANNEL_FUSION §H escalation step 2, MH-99; host fold-in MH-101). A CN
    event moves co-amplified CODING genes too; if one regulates the target, its path Z_s→CN_c→X_c→Y confounds the
    exclusion δ_s. Per segment:
      1. candidate co-located coding genes = **CN co-amplified** with the miRNA-locus CN (`coamp`, the correct
         'co-located' criterion — NOT genomic distance; recovers RB1@13q14 at 1.7 Mb) AND Y-expression-correlated
         (`min_corr`: partial **Spearman**(coding, target | C), z-scored — aligned with the host gate's relatedness
         criterion, MH-101, though the conditioning sets differ ON PURPOSE: coding removes proliferation via C + the
         surv guard, the host gate keeps it because a known co-transcription host makes prolif the genuine mediator);
      2. SPECIFICITY guard — keep only those whose Y-confound SURVIVES proliferation-residualisation (a proliferation
         PC1), so proliferation *markers* (POLR2A/AURKB) on proliferation-*effector* targets are down-selected;
      2b. **`host_genes` (MH-101, per-locus dict): the intronic HOST genes are added as MANDATORY, GUARD-EXEMPT
         candidates** — bypassing the coamp/min_corr/surv scan filters because they carry their own validated
         `host_target_relatedness` gate (partial ρ|CPE,miR≥0.3) + annotation. This UNIFIES the coding & host gates into
         ONE joint conditioning instead of two separate down-weights combined via MAX (which assumed they explain the
         same slice of δ_s; the joint captures the case where a host and a coding neighbour explain DIFFERENT slices);
      3. condition δ_s on {kept coding} ∪ {host} JOINTLY; **ACT ONLY ON REDUCTIONS** — a valid backdoor-blocking
         conditioning REDUCES the confounded δ_s (→ a `pleio_down_weight`); an INCREASE is a collider artifact → IGNORE.
    `scan=False` skips the window scan (cheap host-only path — condition on `host_genes` alone, no coords/CN/prolif load).
    Emits per segment `pleio_down_weight` ∈ [0,1] = the CN-channel down-weight (inflate `s²_π`), NOT a subtraction, +
    `n_host`/`pleio_sources`. Anchor: RB1@13q14 → miR-16→CCND1 42%. Returns per-(gene,family,segment) rows."""
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import families as FAM
    host_genes = host_genes or {}
    asm = assembled if assembled is not None else LD.assemble_gene(gene, he_only=False)
    Y, X, C, _ = asm
    fam_name = str(FAM.family_of([mem[0]]).get(mem[0], mem[0])) if mem else ""
    G = LD._load()["Y"]; lcn = locus_cn()
    coords = _coding_coords() if scan else None
    mc = _mirna_coords() if scan else None
    if scan:                                                               # numpy views + name→index maps (was pandas .loc reindex)
        G_np, g_rows, g_cols = _as_np(G, "rna_np")
        CN_np, cn_rows, cn_cols = _as_np(_gene_cn(), "gene_cn_np")
    foc = {m: _arm_focal_locus(m)["locus_id"] for m in mem}
    segs = {}
    for m in mem:
        if foc.get(m) in lcn.columns and (not scan or m in mc.index):
            segs.setdefault(foc[m], (mc.loc[m, "chrom"], mc.loc[m, "mid"]) if scan else (None, None))
    if scan:
        parts = [p for p in X.index if p in lcn.index and p in g_cols and p in cn_cols]
    else:
        parts = [p for p in X.index if p in lcn.index and p in G.columns]
    if len(parts) < 60 or not segs:
        return pd.DataFrame()
    Xm = X.loc[parts, mem].apply(pd.to_numeric, errors="coerce")
    xfam_raw = np.log2((2.0 ** Xm.to_numpy(float)).sum(1)); y_raw = Y.loc[parts].to_numpy(float)
    Cm0 = C.loc[parts].to_numpy(float)
    if scan:
        pcolG = np.fromiter((g_cols[p] for p in parts), int, len(parts))   # participant→column indices (once per gene)
        pcolCN = np.fromiter((cn_cols[p] for p in parts), int, len(parts))
        _pc = _CACHE.setdefault("prolif_pc", {})                           # PERF: the prolif PC is GENE-level (function of `parts`
        _hit = _pc.get(gene)                                               # only), but was re-SVD'd per FAMILY — cache it per gene
        if _hit is not None and _hit[0] == parts:
            prolif = _hit[1]
        else:
            prow = np.fromiter((g_rows[g2] for g2 in _PROLIF_GENES if g2 in g_rows), int, -1)
            P = G_np[prow][:, pcolG].T; P = (P - P.mean(0)) / (P.std(0) + 1e-9)
            prolif = np.linalg.svd(P - P.mean(0), full_matrices=False)[0][:, 0]
            _pc[gene] = (parts, prolif)
    # PERF (max-speed rewrite): residualisation via PRECOMPUTED lstsq projectors (was N× sklearn LinearRegression fits);
    # co-amp filter VECTORISED over all near genes; CN/G pulled as single blocks; yr/resy hoisted (were per-gene). Numerically
    # identical to the sklearn version (OLS residual on [1, C] / [1, C, prolif] / [1, xf_p]).
    rows = []
    for loc, (chrom, mid) in segs.items():
        Zc_raw = pd.to_numeric(lcn.loc[parts, loc], errors="coerce").to_numpy(float)
        base = np.isfinite(xfam_raw) & np.isfinite(y_raw) & np.isfinite(Zc_raw) & np.isfinite(Cm0).all(1)
        nb = int(base.sum())
        if nb < 60:
            continue
        A0 = np.column_stack([np.ones(nb), Cm0[base]]); P0 = np.linalg.pinv(A0)      # residualiser on [1, C]

        def _r0(vb):
            r = vb - A0 @ (P0 @ vb); s = r.std(); return r / (s or 1.0)

        yb = y_raw[base]
        xf, Zr, yr = _r0(xfam_raw[base]), _r0(Zc_raw[base]), _r0(yb)
        host_set = set(host_genes.get(loc, []))                            # this locus's related hosts (MH-101)
        keep = []                                                          # (gene, resid, is_host)
        if scan:                                                           # window scan for co-amplified coding neighbours
            Ap = np.column_stack([np.ones(nb), Cm0[base], prolif[base]]); Pp = np.linalg.pinv(Ap)  # residualiser on [1, C, prolif]

            def _rp(vb):
                r = vb - Ap @ (Pp @ vb); s = r.std(); return r / (s or 1.0)

            xf_p, yr_p = _rp(xfam_raw[base]), _rp(yb)
            B = np.column_stack([np.ones(nb), xf_p]); PB = np.linalg.pinv(B)          # surv-guard residualiser on [1, xf_p]
            resy = yr_p - B @ (PB @ yr_p)                                             # CONSTANT — hoisted out of the gene loop
            near = list(dict.fromkeys(g2 for g2 in coords[(coords.chrom == chrom) & ((coords.mid - mid).abs() < window)].gene_name
                                      if g2 in g_rows and g2 in cn_rows and g2 != gene))
            if near:
                CNb = CN_np[np.fromiter((cn_rows[g2] for g2 in near), int, len(near))][:, pcolCN][:, base]  # m × nb (numpy, was .loc)
                zc = Zc_raw[base]; zc0 = zc - zc.mean(); zcn = float(np.sqrt(zc0 @ zc0)) or 1.0
                CNc = CNb - CNb.mean(1, keepdims=True); CNn = np.sqrt((CNc ** 2).sum(1)); CNn[CNn == 0] = 1.0
                coamp_r = np.abs(CNc @ zc0) / (CNn * zcn)                             # VECTORISED |corr(CN_i, Z_s)|
                idx = np.where(np.isfinite(CNb).all(1) & (CNb.std(1) > 1e-9) & (coamp_r >= coamp))[0]
                if len(idx):
                    gk = [near[i] for i in idx]
                    Gk = G_np[np.fromiter((g_rows[g2] for g2 in gk), int, len(gk))][:, pcolG][:, base]      # m' × nb (numpy, was .loc)
                    Xc = Gk - (A0 @ (P0 @ Gk.T)).T                                    # residual | C for ALL candidates (matches _r0)
                    Xc = Xc / (Xc.std(1, keepdims=True) + 1e-12)                      # z per gene
                    # VECTORISED partial |Spearman(Xc_j, yr | C)| = |Pearson(rank Xc_j, rank yr)| — one rank + matmul,
                    # replaces the per-candidate scipy `spearmanr` (which re-did cov/pvalue 30k× = the scan's hot loop).
                    ry = rankdata(yr); ry = (ry - ry.mean()) / (ry.std() + 1e-12)
                    Rx = np.apply_along_axis(rankdata, 1, Xc)
                    Rx = (Rx - Rx.mean(1, keepdims=True)) / (Rx.std(1, keepdims=True) + 1e-12)
                    sp = (Rx @ ry) / len(ry)
                    for j in np.where(np.abs(sp) >= min_corr)[0]:                     # candidates passing the Y-corr filter
                        resxc = (xcp := _rp(Gk[j])) - B @ (PB @ xcp)
                        if abs(np.corrcoef(resxc, resy)[0, 1]) >= surv:              # survives proliferation-residualisation
                            keep.append((gk[j], Xc[j], gk[j] in host_set))           # tag scan-found hosts too
        for h in sorted(host_set):                                         # MANDATORY guard-exempt HOST candidates the scan MISSED
            if h in G.index and h != gene and h not in [k[0] for k in keep]:
                keep.append((h, _r0(G.loc[h, parts].to_numpy(float)[base]), True))
        d0 = float(np.linalg.lstsq(np.column_stack([xf, Zr]), yr, rcond=None)[0][1])
        d1 = float(np.linalg.lstsq(np.column_stack([xf, Zr] + [k[1] for k in keep]), yr, rcond=None)[0][1]) if keep else d0
        # ACT ONLY ON REDUCTIONS: valid backdoor block reduces |δ_s|; inflation = collider → ignore
        dw = max(0.0, 1.0 - abs(d1) / abs(d0)) if abs(d0) >= dmin else 0.0
        n_host = sum(1 for k in keep if k[2])
        verdict = ("pleio_confounded" if dw > 0 else ("d~0" if abs(d0) < dmin else "collider_ignore"))
        rows.append({"gene": gene, "family": fam_name, "locus": loc, "delta_s": round(d0, 3),
                     "delta_cond": round(d1, 3), "pleio_down_weight": round(dw, 2), "n_host": n_host, "verdict": verdict,
                     "pleio_sources": ", ".join(f"{k[0]}*" if k[2] else k[0]                          # host tagged *, shown first
                                                for k in sorted(keep, key=lambda k: not k[2])[:topk])})
    return pd.DataFrame(rows)


_GENOMIC_CONTEXT = "mirna_hallmark/output/learned/genomic_context.tsv"
_HOST_REL_TSV = "mirna_hallmark/output/learned/host_target_relatedness.tsv"


def _sense_coding_hosts() -> dict:
    """arm (full `hsa-` name) → [host gene(s)] for every sense_coding_host arm (genomic_context, MH-101)."""
    if "sense_coding_hosts" not in _CACHE:
        import os
        m: dict = {}
        if os.path.exists(_GENOMIC_CONTEXT):
            gc = pd.read_csv(_GENOMIC_CONTEXT, sep="\t")
            for _, r in gc[gc["mir_class"] == "sense_coding_host"].iterrows():
                if isinstance(r["host"], str):
                    m[r["arm"]] = [h.strip() for h in r["host"].split(",") if h.strip()]
        _CACHE["sense_coding_hosts"] = m
    return _CACHE["sense_coding_hosts"]


def _related_host_set() -> set:
    """{(arm, host, target)} where the host RELATES to the target (host_target_relatedness partial ρ≥0.3, MH-101).
    Now PER-LOCUS complete (2026-07-11): host_target_relatedness enumerates every distinct coding host across an arm's
    CN loci (not just the coding-first representative), so a per-locus coding host (CTDSPL alongside CTDSP2, SLIT3
    alongside SLIT2) has its own (arm, host, target) tuple here."""
    if "related_host_set" not in _CACHE:
        import os
        s: set = set()
        if os.path.exists(_HOST_REL_TSV):
            rel = pd.read_csv(_HOST_REL_TSV, sep="\t")
            rel = rel[rel["related"] == True]                               # noqa: E712 (pandas mask)
            s = set(zip(rel["arm"], rel["host"], rel["target"]))
        _CACHE["related_host_set"] = s
    return _CACHE["related_host_set"]


_LOCUS_CTX = "mirna_hallmark/output/learned/locus_context.tsv"


def _locus_coding_host() -> dict:
    """CN locus_id (`MI*`) → coding host gene, for loci classified `sense_coding_host` (PER-LOCUS map, MH-101 fix). Lets
    `exclusion` condition each ACTIVE instrument locus on the host it ACTUALLY sits in — not the arm's coding-first
    representative, which for a multi-locus arm can belong to a DIFFERENT (or SILENT, non-instrumented) locus. Source:
    `genomic_context.locus_host_map()` persisted to `locus_context.tsv`."""
    if "locus_coding_host" not in _CACHE:
        import os
        m: dict = {}
        if os.path.exists(_LOCUS_CTX):
            lc = pd.read_csv(_LOCUS_CTX, sep="\t")
            for _, r in lc[lc["mir_class"] == "sense_coding_host"].iterrows():
                if isinstance(r["host"], str):
                    m[str(r["locus_id"])] = r["host"]
        _CACHE["locus_coding_host"] = m
    return _CACHE["locus_coding_host"]


def host_pleiotropy(gene: str, mem, *, assembled=None, dmin: float = 0.05, topk: int = 3) -> pd.DataFrame:
    """**Host-gene exclusion gate for intronic miRNAs (MH-101 → CN exclusion, design axis BP).** The a-priori extreme
    of `coding_pleiotropy`: a SENSE-intronic miRNA's locus CN co-amplifies its HOST gene by construction (same locus,
    perfect co-amplification), so if the host RELATES to the target the path Z_s→CN_host→host protein→Y bypasses the
    miRNA and violates exclusion. **COMPLEMENTARY to `coding_pleiotropy` — NOT a replacement (measured):** with the
    genome-wide CN matrix the host IS a candidate in the window scan too, but the scan frequently misses it — of 29
    scoped host-confounds, the scan down-weights only 14 and lists the HOST itself in just 5; **15/29 the scan misses
    entirely.** The clearest gap is proliferation-gene hosts (MCM7 for miR-17~92/106b~25): the scan's prolif-survival
    guard deliberately drops them as 'prolif markers', but here MCM7 is the a-priori co-transcription source, so the host
    gate retains it. Conversely the scan catches strong NON-host neighbours the host gate can't (PIK3CG: host_dw 0.09 vs
    coding_dw 0.90) — hence the two are combined via MAX per segment in `exclusion()`. This gate: host from ANNOTATION
    (`genomic_context`, sense_coding_host), gated by the validated `host_target_relatedness` (partial ρ(host,target|CPE)
    ≥0.3). Per focal-locus segment carrying a related-host member: condition δ_s on the host expression, **ACT ONLY ON
    REDUCTIONS** (a valid backdoor block shrinks |δ_s|; inflation = collider → ignore). X_fam is held in the regression,
    so the LEGITIMATE co-transcription path (host→miR→Y) is already accounted — conditioning blocks only the residual
    host-protein path. Emits `host_down_weight` ∈ [0,1] = the CN-channel down-weight (inflate s²_π), NOT a subtraction.
    ⚠ For PROLIFERATION-gene hosts (MCM7) the flag overlaps the MH-100 confound-vs-mechanism question — cross-reference
    `prolif_verdict`, don't blind-apply; the clean additive value is the ~11 NON-prolif hosts (WDR82/CTDSP2/SMC4/AOPEP)."""
    from mirna_hallmark.learned import families as FAM
    asm = assembled if assembled is not None else LD.assemble_gene(gene, he_only=False)
    Y, X, C, _ = asm
    fam_name = str(FAM.family_of([mem[0]]).get(mem[0], mem[0])) if mem else ""
    G = LD._load()["Y"]; lcn = locus_cn()
    lch = _locus_coding_host(); rel = _related_host_set()                   # PER-LOCUS host (MH-101 fix; consistent with exclusion)
    seg_hosts: dict = {}                                                    # ACTIVE locus → {related host genes present in Y}
    for m in mem:
        for loc in _arm_active_loci(m):                                     # each active locus → its OWN coding host
            if loc not in lcn.columns:
                continue
            h = lch.get(loc)
            if h and (m, h, gene) in rel and h in G.index:
                seg_hosts.setdefault(loc, set()).add(h)
    if not seg_hosts:
        return pd.DataFrame()
    parts = [p for p in X.index if p in lcn.index and p in G.columns]
    Xm = X.loc[parts, mem].apply(pd.to_numeric, errors="coerce")
    xfam_raw = np.log2((2.0 ** Xm.to_numpy(float)).sum(1)); y_raw = Y.loc[parts].to_numpy(float)
    Cm0 = C.loc[parts].to_numpy(float)
    rows = []
    for loc, hosts in seg_hosts.items():
        Zc_raw = pd.to_numeric(lcn.loc[parts, loc], errors="coerce").to_numpy(float)
        base = np.isfinite(xfam_raw) & np.isfinite(y_raw) & np.isfinite(Zc_raw) & np.isfinite(Cm0).all(1)
        if base.sum() < 60:
            continue
        Cm = Cm0[base]

        def rz(v):
            r = v[base] - LinearRegression().fit(Cm, v[base]).predict(Cm); return r / (r.std() or 1)

        xf, Zr, yr = rz(xfam_raw), rz(Zc_raw), rz(y_raw)
        hc = [(h, rz(G.loc[h, parts].to_numpy(float))) for h in sorted(hosts)]
        zz = float(Zr @ Zr); gamma = float((Zr @ xf) / zz) if zz > 0 else 0.0                 # first stage γ (CN→X_fam)
        d0 = float(np.linalg.lstsq(np.column_stack([xf, Zr]), yr, rcond=None)[0][1])          # δ_s (as in exclusion T1)
        d1 = float(np.linalg.lstsq(np.column_stack([xf, Zr] + [h[1] for h in hc]), yr, rcond=None)[0][1])
        dw = max(0.0, 1.0 - abs(d1) / abs(d0)) if abs(d0) >= dmin else 0.0    # ACT ONLY ON REDUCTIONS
        verdict = ("host_confounded" if dw > 0 else ("d~0" if abs(d0) < dmin else "collider_ignore"))
        # CONFIDENCE (δ_s sign): a residual-miR leak into δ_s has sign −sign(γ) (the miR REPRESSES ⇒ β_mir<0). So δ_s of the
        # SAME sign as γ CANNOT be lost-miR signal ⇒ UNAMBIGUOUS host pleiotropy (clean_antirepression). Opposite sign
        # (repression_direction) COULD be residual-miR (X_fam mismeasures the dose) → grade via the miR-partialled relatedness + prolif_verdict.
        confidence = "clean_antirepression" if d0 * gamma > 0 else "repression_direction"
        rows.append({"gene": gene, "family": fam_name, "locus": loc, "gamma_s": round(gamma, 3),
                     "delta_s": round(d0, 3), "delta_cond": round(d1, 3), "host_down_weight": round(dw, 2),
                     "verdict": verdict, "confidence": confidence, "hosts": ", ".join(h[0] for h in hc[:topk])})
    return pd.DataFrame(rows)


def between_family_exclusion(gene: str, mem, foc=None, *, window: float = 1e6, min_n: int = 60,
                             min_F: float = 10.0, edges=None, assembled=None, fam_of=None) -> pd.DataFrame:
    """**Gap 2B (CN_INSTRUMENT §9): between-family exclusion — condition `p_S` on co-located OTHER-seed co-targeters.**
    Ring-1's `seg_pS` (per-segment CN reduced-form share, `family_causal_attribution`) is confounded when a
    co-located OTHER-seed family also targets `gene`: the cluster CN moves *them* too, so `seg_pi` (CN_seg→target)
    mixes our family's causal leg with theirs. This is the frequentist interim of the exclusion step:

      1. **enumerate** co-located other-seed co-targeter families via `pleiotropy()` — union of each member's
         co-located gene-targeting arms, mapped to seed families, our own family removed;
      2. **condition** each segment's reduced form on those families' pooled abundance → `seg_pi_cond`
         (blocks the CN_seg → other-family-abundance → target mediation path), reporting
         `pi_shift = |seg_pi_cond| − |seg_pi_raw|` and the re-normalised `seg_pS_cond`;
      3. **attribute** the coupling among {our family} ∪ {co-located other-seed families} with the between-family
         exoneration gate + Shapley (`attribution._exonerate_between` + `shapley_identity`) → `between_family_share`
         (our family's collinearity-fair ownership) and `family_exonerated` (our family is itself a collinear rider).

    A segment is EXCLUSION-ROBUST when conditioning barely moves `pi` and our family survives exoneration; CONFOUNDED
    when `pi` collapses/flips under conditioning or our family is exonerated. `between_clean` = no co-located other-seed
    co-targeter (exclusion holds trivially — the common case). Returns per-segment rows + family-level columns.
    CAVEAT: conditions on the other families' *abundance* (mediation adjustment); the exogenous CN-anchor form of this
    (cross-family posterior covariance) is the parked Bayes-unification item (§9)."""
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned.attribution import _exonerate_between, shapley_identity
    Y, X, C, _ = assembled if assembled is not None else LD.assemble_gene(gene, he_only=False)
    if foc is None:
        foc = {m: _arm_focal_locus(m)["locus_id"] for m in mem}
    lcn = locus_cn()
    mem = [m for m in mem if m in X.columns and foc.get(m) in lcn.columns]
    if len(mem) < 2:
        return pd.DataFrame()
    edges = edges if edges is not None else D.high_evidence_edges()
    fo = fam_of if fam_of is not None else FAM.family_of(list(X.columns))
    our_fam = str(fo.get(mem[0], mem[0]))
    # (1) co-located OTHER-seed co-targeters of `gene` (pleiotropy union over our members → other seed families)
    coloc: set = set()
    for a in mem:
        coloc |= set(pleiotropy(a, gene, window=window, edges=edges).get("coloc_targeting_gene") or [])
    other_fams: dict = {}
    for c in sorted(coloc):                                                 # sorted → deterministic family order (Shapley)
        if c not in X.columns:
            continue
        f = str(fo.get(c, c))
        if f != our_fam:                                                    # only OTHER-seed families confound
            other_fams.setdefault(f, []).append(c)
    other_fams = {f: other_fams[f] for f in sorted(other_fams)}
    parts = X.index.intersection(lcn.index)
    Cm = C.loc[parts].to_numpy(float)
    yv = Y.loc[parts].to_numpy(float)

    def _pool(arms):                                                        # log2 Σ 2^x pooled abundance
        arms = [a for a in arms if a in X.columns]
        return np.log2((2.0 ** X.loc[parts, arms].to_numpy(float)).sum(1))

    Xother = np.column_stack([_pool(a) for a in other_fams.values()]) if other_fams else None
    # (2) per-segment reduced form, raw vs conditioned on the other-seed pool
    CN = {m: pd.to_numeric(lcn.loc[parts, foc[m]], errors="coerce").to_numpy(float) for m in mem}
    segs = []
    for m in mem:                                                          # group members by CN collinearity (|ρ|>0.999)
        for s in segs:
            if abs(np.corrcoef(np.nan_to_num(CN[m]), np.nan_to_num(CN[s[0]]))[0, 1]) > 0.999:
                s.append(m); break
        else:
            segs.append([m])
    import warnings
    seg_info = []                                                          # (members, pi_raw, pi_cond, F)
    for s in segs:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            zc = np.nanmean(np.column_stack([CN[m] for m in s]), axis=1)
        xf = X.loc[parts, s].mean(axis=1).to_numpy(float)
        ok = np.isfinite(zc) & np.isfinite(yv) & np.isfinite(xf) & np.isfinite(Cm).all(1)
        if ok.sum() < min_n:
            seg_info.append((s, np.nan, np.nan, np.nan)); continue
        F = _first_stage_F(xf[ok], zc[ok], Cm[ok])
        pi_raw = float(spearmanr(_resid(zc[ok], Cm[ok]), _resid(yv[ok], Cm[ok])).correlation)
        pi_cond = pi_raw
        if Xother is not None:
            okc = ok & np.isfinite(Xother).all(1)
            if okc.sum() >= min_n:
                Ca = np.column_stack([Cm[okc], Xother[okc]])              # C + other-seed pool → mediation-adjusted
                pi_cond = float(spearmanr(_resid(zc[okc], Ca), _resid(yv[okc], Ca)).correlation)
        seg_info.append((s, pi_raw, pi_cond, F))

    def _tot(idx):                                                        # denom for pS shares over F-strong segments
        return sum(abs(t[idx]) for t in seg_info if t[idx] == t[idx] and t[3] == t[3] and t[3] > min_F) or 1.0

    tot_raw, tot_cond = _tot(1), _tot(2)
    # (3) between-family attribution: exonerate collinear riders, then Shapley split among {our_fam} ∪ other_fams
    supp = [our_fam] + list(other_fams)
    bshare, family_exon = 1.0, False
    if other_fams:
        Xf = pd.DataFrame(index=parts)
        Xf[our_fam] = _pool(mem)
        for f, arms in other_fams.items():
            Xf[f] = _pool(arms)
        okf = np.isfinite(Xf.to_numpy(float)).all(1) & np.isfinite(yv) & np.isfinite(Cm).all(1)
        if okf.sum() >= min_n:
            Xf_ok, yv_ok, Cm_ok = Xf.loc[okf], yv[okf], Cm[okf]
            survivors, exon = _exonerate_between(Xf_ok, supp, yv_ok, Cm_ok)
            family_exon = our_fam in exon
            if not survivors:
                survivors = supp
            yr = _resid(yv_ok, Cm_ok)

            def _rz(col):                                                 # C-residualised, z-scored (identity_vs_magnitude style)
                v = _resid(Xf_ok[col].to_numpy(float), Cm_ok)
                return (v - v.mean()) / (v.std(ddof=0) + 1e-9)

            phi = np.clip(shapley_identity(np.column_stack([_rz(f) for f in survivors]),
                                           np.ones(len(survivors)), yr), 0, None)
            if our_fam in survivors and phi.sum() > 0:
                bshare = float(phi[survivors.index(our_fam)] / phi.sum())
            elif family_exon:
                bshare = 0.0
    # family-level verdict
    good = [(pr, pc) for (_, pr, pc, F) in seg_info if pr == pr and F == F and F > min_F]
    if not other_fams:
        verdict = "between_clean"
    elif not good:
        verdict = "weak_instrument"
    elif family_exon:
        verdict = "confounded_rider"
    else:
        pr, pc = max(good, key=lambda t: abs(t[0]))                       # driving (max |pi_raw|) segment
        atten = max(0.0, abs(pr) - (abs(pc) if pc == pc else 0.0)) / (abs(pr) or 1.0)
        flip = pc == pc and pr < 0 <= pc
        verdict = "confounded" if (flip or atten > 0.5) else ("partly_attenuated" if atten > 0.2 else "robust")
    rows = []
    for si, (s, pi_raw, pi_cond, F) in enumerate(seg_info):
        strong = pi_raw == pi_raw and F == F and F > min_F
        rows.append({
            "gene": gene, "family": our_fam, "segment": si, "members": ",".join(m.replace("hsa-", "") for m in s),
            "seg_F": round(F, 1) if F == F else np.nan,
            "seg_pi_raw": round(pi_raw, 3) if pi_raw == pi_raw else np.nan,
            "seg_pi_cond": round(pi_cond, 3) if pi_cond == pi_cond else np.nan,
            "pi_shift": round((abs(pi_cond) - abs(pi_raw)), 3) if (pi_raw == pi_raw and pi_cond == pi_cond) else np.nan,
            "seg_pS_raw": round(abs(pi_raw) / tot_raw, 3) if strong else 0.0,
            "seg_pS_cond": round(abs(pi_cond) / tot_cond, 3) if (strong and pi_cond == pi_cond) else 0.0,
        })
    df = pd.DataFrame(rows)
    if len(df):
        df["n_other_seed_families"] = len(other_fams)
        df["other_seed_families"] = ";".join(sorted(other_fams)) if other_fams else ""
        df["between_family_share"] = round(bshare, 3)
        df["family_exonerated"] = family_exon
        df["exclusion_verdict"] = verdict
    return df


def run_between_family_exclusion(genes=None, *, min_members: int = 2, window: float = 1e6, progress: int = 100,
                                 out: str = "mirna_hallmark/output/learned/instrument/between_family_exclusion.tsv"
                                 ) -> pd.DataFrame:
    """Batch **Gap 2B**: `between_family_exclusion` over every seed family (≥`min_members` regulators) of each gene in
    `genes` (default = HE gene universe). Persists per-(gene, family, segment) raw-vs-conditioned `pi`/`pS` + the
    family-level between-family Shapley share and exclusion verdict. CLI: `python -m mirna_hallmark.learned.instrument --gap2b`."""
    from pathlib import Path
    from mirna_hallmark.learned import families as FAM
    if genes is None:
        genes = sorted(D.high_evidence_edges()["gene"].dropna().astype(str).unique())
    edges = D.high_evidence_edges()
    rows = []
    for gi, g in enumerate(genes):
        if progress and gi % progress == 0:
            print(f"[gap2b] {gi}/{len(genes)} genes, {len(rows)} family-tables", flush=True)
        try:
            asm = LD.assemble_gene(g, he_only=False)
        except Exception:
            continue
        X = asm[1]
        fo = FAM.family_of(list(X.columns))
        fams: dict = {}
        for a in X.columns:
            f = fo.get(a)
            if f:
                fams.setdefault(f, []).append(a)
        for f, mem in fams.items():
            if len(mem) < min_members:
                continue
            try:
                df = between_family_exclusion(g, mem, window=window, edges=edges, assembled=asm, fam_of=fo)
            except Exception:
                continue
            if len(df):
                rows.append(df)
    if not rows:
        return pd.DataFrame()
    out_df = pd.concat(rows, ignore_index=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out, sep="\t", index=False)
    fam_lvl = out_df.drop_duplicates(["gene", "family"])
    vc = fam_lvl["exclusion_verdict"].value_counts().to_dict()
    conf = fam_lvl[fam_lvl["n_other_seed_families"] > 0]
    print(f"[gap2b] {len(out_df)} segment-rows | {fam_lvl.shape[0]} (gene,family) | "
          f"{int((fam_lvl['n_other_seed_families'] > 0).sum())} with co-located other-seed co-targeter | "
          f"verdicts={vc} → {out}")
    if len(conf):
        print(f"[gap2b] among confoundable families: median |pi_shift| (driving seg) reported per-row; "
              f"robust={int((conf['exclusion_verdict'] == 'robust').sum())} "
              f"partly={int((conf['exclusion_verdict'] == 'partly_attenuated').sum())} "
              f"confounded={int(conf['exclusion_verdict'].isin(['confounded', 'confounded_rider']).sum())}")
    return out_df


def exclusion(gene: str, mem, foc=None, *, min_n: int = 60, min_F: float = 10.0,
              delta_thr: float = 0.05, assembled=None, n_boot: int = 0, boot_seed: int = 0,
              coding: bool = False, host: bool = False, shuffle_cn: bool = False, shuffle_seed: int = 0) -> pd.DataFrame:
    """**Scope-3 CN-instrument exclusion gate (CHANNEL_FUSION_DESIGN §2.2, MH-89).** Two COMPLEMENTARY tests on the
    **family exposure** `X_fam = log2 Σ_m 2^x` (the coupling unit, §8; = `family_multi_iv`'s exposure — off-segment
    same-family arms are constituents of X_fam, NOT controls):
      • **T1 — runs on every edge.** Per independent-CN segment s, `δ_s` = the direct effect of the segment CN on Y that
        does NOT flow through X_fam: regress `(Y|C)` on `[X_fam|C, Z_s|C]`, read the **`Z_s` coefficient δ_s + SE** (NOT
        the ratio `π_cond/π_raw` — that inherits the Wald small-`π_raw` instability). `δ_s ≠ 0` ⇒ pleiotropy (catches
        HOMOGENEOUS pleiotropy, over-ID's blind spot). CONSERVATIVE (collider if a shared program ∉ C ⇒ over-flags).
      • **over-ID (Hansen-J) — where ≥2 independent-CN segments.** Reused from `family_multi_iv` (`overid_p`, collider-
        free, catches HETEROGENEOUS pleiotropy).
    Admission = `F>10 ∧ T1 δ_s≈0 ∧ over-ID-clean(where powered)`. Arm-level 'which arm' is ORTHOGONAL (§8, Ring-1 L2),
    not part of exclusion. Returns per-segment rows + family verdict/weight."""
    from scipy.stats import t as tdist, chi2
    from mirna_hallmark.learned import families as FAM
    asm = assembled if assembled is not None else LD.assemble_gene(gene, he_only=False)
    Y, X, C, _ = asm
    if foc is None:
        foc = {m: _arm_focal_locus(m)["locus_id"] for m in mem}
    lcn = locus_cn()
    mem = [m for m in mem if m in X.columns and foc.get(m) in lcn.columns]
    if len(mem) < 1:
        return pd.DataFrame()
    fam_name = str(FAM.family_of([mem[0]]).get(mem[0], mem[0]))
    parts = X.index.intersection(lcn.index)
    # INSTRUMENT columns = per (arm, ACTIVE-SOURCE locus), DECOUPLED from the exposure: X_fam stays per-ARM (below), but a
    # multi-active-source arm contributes one CN instrument PER active locus (multi_iv-complete; the extra columns feed the
    # segment grouping + over-ID). For single-active-source arms (783/785) this is exactly the old focal-per-arm behaviour.
    inst = [(m, L) for m in mem for L in _arm_active_loci(m) if L in lcn.columns]
    if not inst:
        return pd.DataFrame()
    Xm = X.loc[parts, mem].apply(pd.to_numeric, errors="coerce").to_numpy(float)      # per ARM → X_fam (the exposure)
    CNm = np.column_stack([pd.to_numeric(lcn.loc[parts, L], errors="coerce").to_numpy(float) for (_, L) in inst])
    y, Cm = Y.loc[parts].to_numpy(float), C.loc[parts].to_numpy(float)
    ok = np.isfinite(Xm).all(1) & np.isfinite(CNm).all(1) & np.isfinite(y) & np.isfinite(Cm).all(1)
    Xm, CNm, y, Cm = Xm[ok], CNm[ok], y[ok], Cm[ok]
    if shuffle_cn:                                                                    # SHUFFLED-CN NULL (§C): permute the
        CNm = CNm[np.random.default_rng(shuffle_seed).permutation(len(CNm))]          # instrument rows → dead instrument
    n, k = len(y), len(inst)                                                          # k = # instrument columns (≥ #arms)
    if n < min_n:
        return pd.DataFrame()
    xfam_raw = np.log2((2.0 ** Xm).sum(1))                          # X_fam (raw) — sum over ARMS, not instrument columns
    # --- vectorised residualisation: ONE lstsq residualises X_fam, Y and every member CN on [1, C] (was k+2 sklearn fits) ---
    Cd = np.column_stack([np.ones(n), Cm])
    V = np.column_stack([xfam_raw, y, CNm])
    R = V - Cd @ np.linalg.lstsq(Cd, V, rcond=None)[0]             # residuals | C
    xfam_r, y_r, CN_r = R[:, 0], R[:, 1], R[:, 2:]                 # RAW residuals (over-ID + first-stage F)

    def _zc(v):
        sd = float(np.std(v)); return v / sd if sd > 1e-12 else v

    xr, yr = _zc(xfam_r), _zc(y_r)                                 # z-scored (T1, standardised scale)
    Zr = np.column_stack([_zc(CN_r[:, j]) for j in range(k)])
    tss = float(xfam_r @ xfam_r) or 1.0
    df2 = n - Cm.shape[1] - 2                                      # first-stage F denominator dof (= _first_stage_F)
    ncond = Cm.shape[1] + 1
    # ONE k×k CN correlation matrix (was ~2k² pairwise corrcoef calls — segment grouping + over-ID dedup reuse it)
    Rcc = np.abs(np.nan_to_num(np.corrcoef(CN_r, rowvar=False))) if k > 1 else np.array([[1.0]])
    segs = []                                                     # group members by CN collinearity |ρ|>0.999
    for j in range(k):
        for s in segs:
            if Rcc[j, s[0]] > 0.999:
                s.append(j); break
        else:
            segs.append([j])
    # --- over-ID (Hansen-J) INLINE on raw residuals — EXACT match to family_multi_iv (verified) — collider-free
    #     HETEROGENEOUS test; over-powered at large n ⇒ family-level CONTEXT flag, not a per-edge verdict ---
    zznorm = (CN_r ** 2).sum(0)                                   # vectorised solo first-stage F (was a k-loop)
    numj = CN_r.T @ xfam_r
    explained = np.where(zznorm > 0, numj ** 2 / np.where(zznorm > 0, zznorm, 1.0), 0.0)
    sfj = tss - explained
    with np.errstate(divide="ignore", invalid="ignore"):
        soloF = np.where(sfj > 0, explained / (sfj / (n - 1)), np.nan)
    strong_idx = [j for j in range(k) if soloF[j] == soloF[j] and soloF[j] > 10]
    over_id_powered = False; overid_p = fam_F = np.nan; n_indep = 0
    if strong_idx:
        Zs = CN_r[:, strong_idx]; g = np.linalg.lstsq(Zs, xfam_r, rcond=None)[0]
        xhat = Zs @ g; ssr = float(((xfam_r - xhat) ** 2).sum()); ks = len(strong_idx)
        fam_F = float(((tss - ssr) / ks) / (ssr / (n - ks))) if ssr > 0 and n > ks else np.nan
        beta = float((xhat @ y_r) / (xhat @ xhat)) if (xhat @ xhat) > 0 else np.nan
        indep, kept = [], []
        for j in sorted(strong_idx, key=lambda a: -soloF[a]):
            if all(Rcc[j, o] < 0.999 for o in kept):
                indep.append(j); kept.append(j)
        n_indep = len(kept)
        if n_indep > 1:
            Zi = CN_r[:, indep]; e = y_r - beta * xfam_r
            gm = Zi.T @ e / n; S = (Zi * (e ** 2)[:, None]).T @ Zi / n
            try:
                overid_p = float(chi2.sf(float(n * gm @ np.linalg.solve(S, gm)), n_indep - 1))
            except np.linalg.LinAlgError:
                pass
        over_id_powered = bool(n_indep > 1)
    dof_fs, dof_t1 = n - ncond - 1, n - ncond - 2                  # residual dof (C already partialled)
    pdw = {}; psrc = {}                                            # JOINT coding+host pleiotropy down-weight per locus (MH-99/101)
    if coding or host:                                             # ONE joint conditioning (host = guard-exempt candidate), no MAX
        hg = {}                                                    # per-locus related host genes (from genomic_context + relatedness)
        if host:
            lch = _locus_coding_host(); rel = _related_host_set()      # PER-LOCUS host (MH-101 fix): each ACTIVE locus's OWN host
            for m in mem:
                for lc in _arm_active_loci(m):                          # silent loci aren't instruments → not candidates
                    h = lch.get(lc) if lc in lcn.columns else None     # THIS locus's coding host (not the arm's representative)
                    if h and (m, h, gene) in rel:
                        hg.setdefault(lc, set()).add(h)
            hg = {k: sorted(v) for k, v in hg.items()}
        # `coding=True` ⇒ FULL-universe window scan on EVERY family (the generic MH-99 coding-pleiotropy gate); host-only
        # (coding=False) ⇒ no scan, just the guard-exempt host candidates. The old host-family restriction was a speed
        # workaround — now that the scan is ~7 ms/family (numpy pulls + cached prolif PC) it's affordable everywhere.
        do_scan = coding
        try:
            cp = coding_pleiotropy(gene, mem, assembled=asm, host_genes=hg, scan=do_scan)
            if len(cp):
                pdw = dict(zip(cp["locus"], cp["pleio_down_weight"]))
                psrc = dict(zip(cp["locus"], cp["pleio_sources"]))
        except Exception:
            pdw = {}; psrc = {}
    rows = []
    for si, s in enumerate(segs):
        zs = Zr[:, s].mean(1)                                      # segment CN | C, z
        zz = float(zs @ zs) or 1.0
        pi_raw = float((zs @ yr) / np.sqrt(zz * float(yr @ yr)))   # standardised reduced form (context/decomposition check)
        # first stage  a = Z_s → X_fam  (both z, C-residualised) + SE_a  (the channel's γ_s)
        a = float((zs @ xr) / zz)
        rss_a = float(((xr - a * zs) ** 2).sum())
        se_a = float(np.sqrt((rss_a / dof_fs) / zz)) if dof_fs > 0 else np.nan
        F_s = float((float(xr @ xr) - rss_a) / (rss_a / dof_fs)) if rss_a > 0 and dof_fs > 0 else np.nan
        # T1  yr ~ [X_fam, Z_s] :  b = X_fam→Y|Z_s ,  δ_s = direct effect (pleiotropy)
        Mt = np.column_stack([xr, zs])
        coef, *_ = np.linalg.lstsq(Mt, yr, rcond=None)
        b, delta = float(coef[0]), float(coef[1])
        resid = yr - Mt @ coef
        sig2 = float(resid @ resid) / dof_t1 if dof_t1 > 0 else np.nan
        try:
            covm = sig2 * np.linalg.inv(Mt.T @ Mt); se_b = float(np.sqrt(covm[0, 0])); se_d = float(np.sqrt(covm[1, 1]))
        except np.linalg.LinAlgError:
            se_b = se_d = np.nan
        # ⛔⛔ `pi_causal` IS **NOT** AN INSTRUMENTAL-VARIABLE ESTIMATE — RETRACTED (MH-124r, 2026-07-13).
        # It is the PRODUCT-OF-COEFFICIENTS **MEDIATION** estimator: a = the exogenous first stage (CN_s → X_fam),
        # but **b is the OBSERVATIONAL OLS SLOPE of Y on X_fam given [C, Z_s]** — the instrument-ORTHOGONAL
        # (endogenous) variation, carrying every confound. It reproduces the plain no-copy-number partial
        # correlation at ρ≈0.93, and the whole "3× HE-vs-decoy separation" was REPRODUCED FROM A SIMULATION
        # WITH ZERO CAUSAL EFFECT. The only genuinely exogenous factor, `a`, is SITE-BLIND (HE +0.1993 vs
        # DECOY +0.1984, ratio 1.004×, p=0.20).
        # ⇒ **DO NOT USE `pi_causal` AS A CAUSAL QUANTITY OR AS AN ATTRIBUTION ARBITER.** It is a channel-fusion
        #   (mediation) statistic, which is what it was originally built for. The only quantity here that
        #   reaches Y solely through copy number is the REDUCED FORM `pi_raw` — and that is SITE-BLIND too
        #   (arm-clustered p=0.115; the within-arm FE replacement also failed, MH-133).
        # See MH124_ANTICOUPLING_VALIDITY.md (both retraction banners), registry MH-124r / MH-133.
        pi_causal = a * b                                         # ⚠ MEDIATION (a × an OBSERVATIONAL slope), NOT an IV
        # s²_π = SAMPLING variance of the product a·b + δ_s² BIAS bound. The δ_s² term is essential (parameter-free): the
        # pleiotropy correction δ_s is precisely estimated but COLLIDER-biased, so a big correction must down-weight the edge —
        # without it a pleiotropic edge gets HIGHER precision than a clean one (verified backwards on miR-16→CCND1). If the whole
        # δ_s were collider artifact the error in pi_causal is δ_s, so δ_s² is the honest bias-variance floor.
        if n_boot > 0:
            # NONPARAMETRIC BOOTSTRAP of the product a·b (WHATS_NEXT ①): the analytic Sobel `b²se_a²+a²se_b²` is a delta-method
            # approx that is ANTI-CONSERVATIVE when the first stage γ=a is weak/near-0 (the small-γ regime that makes this
            # channel weak). Resampling patients captures the product's true (skewed) sampling variance + the (a,b) covariance,
            # honestly INFLATING s²_π on weak-γ segments → lower admission weight (the F>10 gate's job, made continuous).
            brng = np.random.default_rng(boot_seed + si); pib = np.empty(n_boot)
            for bi in range(n_boot):
                ix = brng.integers(0, n, n); zsi = zs[ix]; zzi = float(zsi @ zsi) or 1.0
                ai = float((zsi @ xr[ix]) / zzi)
                cbi, *_ = np.linalg.lstsq(np.column_stack([xr[ix], zsi]), yr[ix], rcond=None)
                pib[bi] = ai * float(cbi[0])
            samp_var = float(np.var(pib, ddof=1))
        else:
            samp_var = b ** 2 * se_a ** 2 + a ** 2 * se_b ** 2           # analytic Sobel (delta-method) sampling variance
        s2_pi = samp_var + delta ** 2                                    # → the channel s²_π (sampling + pleiotropy-bias)
        se_pi = float(np.sqrt(s2_pi)) if s2_pi == s2_pi and s2_pi >= 0 else np.nan
        prec = float(1.0 / s2_pi) if s2_pi == s2_pi and s2_pi > 1e-12 else 0.0   # CONTINUOUS admission weight (no threshold)
        flag = ("weak" if not (F_s == F_s and F_s > min_F) else    # DESCRIPTIVE label only — NOT used to filter
                "pleiotropic" if abs(delta) > delta_thr else "clean")
        seg_loc = inst[s[0]][1] if s else None                     # this segment's locus (JOINT coding+host down-weight lookup, MH-99/101)
        pw = float(pdw.get(seg_loc, 0.0))
        s2_pi_pleio = s2_pi / max(1e-6, 1.0 - pw)                  # pleiotropy → inflate s²_π (weaker CN pull), NOT a subtraction
        # δ_s-sign CONFIDENCE: a residual-miR leak has sign −sign(γ), so δ_s SAME sign as the first-stage γ (=a) can't be
        # lost-miR ⇒ clean_antirepression (unambiguous); opposite ⇒ repression_direction (could be miR-leak or co-targeting family).
        conf = ("clean_antirepression" if delta * a > 0 else "repression_direction") if pw > 0 else ""
        rows.append({"gene": gene, "family": fam_name, "segment": si,
                     "seg_members": ",".join(dict.fromkeys(inst[j][0].replace("hsa-", "") for j in s)),  # arms (deduped; multi-locus)
                     "seg_F": round(F_s, 1) if F_s == F_s else np.nan,
                     "gamma_s": round(a, 3), "b_fam": round(b, 3),
                     "pi_raw": round(pi_raw, 3), "delta_s": round(delta, 3),
                     "delta_se": round(se_d, 3) if se_d == se_d else np.nan,
                     "pi_causal": round(pi_causal, 4),             # channel π̂_s (pleiotropy-corrected mediated effect)
                     "se_pi_causal": round(se_pi, 4) if se_pi == se_pi else np.nan,
                     "s2_pi": s2_pi, "pleio_down_weight": round(pw, 2), "s2_pi_pleio": s2_pi_pleio,
                     "pleio_confidence": conf, "pleio_sources": psrc.get(seg_loc, "") if pw > 0 else "",
                     "precision": round(prec, 1), "seg_flag": flag})
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    df["n_segments"] = len(segs)
    df["n_indep_cn"] = n_indep
    df["over_id_powered"] = over_id_powered
    df["overid_p"] = round(overid_p, 4) if overid_p == overid_p else np.nan   # heterogeneity CONTEXT (over-powered)
    df["family_F"] = round(fam_F, 1) if fam_F == fam_F else np.nan
    # Channel weight = `A_cn = γ_s²/s²_π` (§1, the contribution to β's precision) with the SOFT null-correction (C4, default
    # 2026-07-12): `A_cn·(1−c/F)₊` = `(γ²−c·se_a²)₊/s²_π` (since se_a²=γ²/F). This ramps from 0 at F=c to 1 as F→∞ — 0 for
    # F≤c (same admitted SET as the old hard F>10 cliff) but it down-weights the borderline F∈[c,2c] band the cliff
    # over-credited. VERIFIED a strictly cleaner gate than C5: same ranking (ρ=0.995) + 4.4× less shuffled-CN null leakage
    # (SNR 43→124, `soft_fweight_sweep`, §C). `_ADMIT_SOFT_C` is the F-scale null-corrector.
    _soft = np.clip(1.0 - _ADMIT_SOFT_C / pd.to_numeric(df["seg_F"], errors="coerce"), 0.0, None).fillna(0.0)  # (1−c/F)₊
    df["beta_weight"] = np.where(df["seg_flag"] != "weak", (df["gamma_s"] ** 2 / df["s2_pi"]) * _soft, 0.0).round(1)
    # pleiotropy-adjusted admission weight: A_cn on the JOINT coding+host-inflated s²_π (one gate, no MAX; MH-99/101)
    df["beta_weight_pleio"] = np.where(df["seg_flag"] != "weak", (df["gamma_s"] ** 2 / df["s2_pi_pleio"]) * _soft, 0.0).round(1)
    df["exclusion_weight"] = df["beta_weight"]
    # seg_dose_share = this segment's ESTIMATOR-grade share of the family's CN identification (β-precision) — the locus-level
    # dose-leverage the abundance `fam_dominant_arm` proxies (MH-101 (a)↔(c) connection). A host down-weight on a low-share
    # segment is correct-but-low-leverage (the family's dose is instrumented by another locus); makes the abundance dominant-arm
    # annotation honest (they disagree ~50% on multi-arm families) WITHOUT touching the arithmetic (dose is already in gamma_s).
    _tot_bw = float(df.loc[df["seg_flag"] != "weak", "beta_weight"].sum())
    df["seg_dose_share"] = np.where(df["seg_flag"] != "weak", df["beta_weight"] / (_tot_bw or 1.0), 0.0).round(3)
    strong = df[df["seg_flag"] != "weak"]
    df["family_verdict"] = ("weak" if not len(strong) else
                            "has_clean" if (strong["seg_flag"] == "clean").any() else "pleiotropic")
    return df


_SOFT_SWEEP_OUT = "mirna_hallmark/output/learned/instrument/soft_fweight_sweep.tsv"


def soft_fweight_sweep(cs=(5, 10, 20),
                       real_path: str = "mirna_hallmark/output/learned/instrument/exclusion.tsv",
                       null_path: str = "mirna_hallmark/output/learned/instrument/exclusion_nullCN.tsv",
                       out: str = _SOFT_SWEEP_OUT) -> pd.DataFrame:
    """**C-axis admission-weight diagnostic (CHANNEL_FUSION §C, board ④).** Compares the CURRENT hard F>10 gate (C5,
    `beta_weight = γ²/s²_π · 1[F>10]`) with the SOFT null-corrected weight (C4): `w_soft(c) = (γ²−c·se_a²)₊ / s²_π`.
    Since the first-stage `se_a² = γ²/F`, this equals `(γ²/s²_π)·(1−c/F)₊` — the raw weight times a soft ramp that is 0
    below F=c and →1 as F→∞ (STRICTER than the hard cliff just above threshold). Evaluated on the REAL table and a
    SHUFFLED-CN NULL (`exclusion(shuffle_cn=True)` = dead instrument): a good gate suppresses the null. Reports, per gate,
    the real vs null total admitted weight + the real/null signal-to-noise ratio, the rank-agreement with the hard gate,
    and the borderline mass (fraction of real weight in F∈[c,2c] where soft and hard diverge). Pure arithmetic on the two
    persisted tables. Deliverable: `soft_fweight_sweep.tsv`."""
    import os
    if not (os.path.exists(real_path) and os.path.exists(null_path)):
        return pd.DataFrame()
    from scipy.stats import spearmanr
    real = pd.read_csv(real_path, sep="\t"); null = pd.read_csv(null_path, sep="\t")

    def raw(df):
        return (pd.to_numeric(df["gamma_s"], errors="coerce") ** 2 / pd.to_numeric(df["s2_pi"], errors="coerce")).fillna(0.0)

    def soft(df, c):
        F = pd.to_numeric(df["seg_F"], errors="coerce")
        return (raw(df) * np.clip(1.0 - c / F, 0.0, None)).fillna(0.0)

    hard_real = pd.to_numeric(real["beta_weight"], errors="coerce").fillna(0.0)   # C5: γ²/s²_π gated at F>10
    hard_null = pd.to_numeric(null["beta_weight"], errors="coerce").fillna(0.0)
    Fr = pd.to_numeric(real["seg_F"], errors="coerce")
    rows = [dict(gate="hard_F>10", real_total=round(float(hard_real.sum()), 1), null_total=round(float(hard_null.sum()), 2),
                 snr=round(float(hard_real.sum()) / max(float(hard_null.sum()), 1e-9), 1),
                 real_admit=int((hard_real > 0).sum()), null_admit=int((hard_null > 0).sum()),
                 rho_vs_hard=1.0, borderline_massfrac=np.nan)]
    for c in cs:
        sr, sn = soft(real, c), soft(null, c)
        bl = float(sr[(Fr >= c) & (Fr < 2 * c)].sum()) / max(float(sr.sum()), 1e-9)   # F∈[c,2c] mass (where soft≠hard most)
        rows.append(dict(gate=f"soft_c={c}", real_total=round(float(sr.sum()), 1), null_total=round(float(sn.sum()), 2),
                         snr=round(float(sr.sum()) / max(float(sn.sum()), 1e-9), 1),
                         real_admit=int((sr > 0).sum()), null_admit=int((sn > 0).sum()),
                         rho_vs_hard=round(float(spearmanr(hard_real, sr).correlation), 3),
                         borderline_massfrac=round(bl, 3)))
    res = pd.DataFrame(rows)
    from pathlib import Path
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    return res


_EXR: dict = {"min_members": 1, "coding": False, "host": False, "shuffle_cn": False, "shuffle_seed": 0}


def _exclusion_gene(g):
    """One gene's family tables (module-level so ProcessPoolExecutor can fork it). Config in `_EXR`."""
    from mirna_hallmark.learned import families as FAM
    try:
        asm = LD.assemble_gene(g, he_only=False)
    except Exception:
        return []
    X = asm[1]; fo = FAM.family_of(list(X.columns)); fams: dict = {}
    for a in X.columns:
        f = fo.get(a)
        if f:
            fams.setdefault(f, []).append(a)
    out = []
    for f, m in fams.items():
        if len(m) < _EXR["min_members"]:
            continue
        try:
            import zlib
            df = exclusion(g, m, assembled=asm, coding=_EXR["coding"], host=_EXR["host"],
                           shuffle_cn=_EXR["shuffle_cn"],
                           shuffle_seed=_EXR["shuffle_seed"] + zlib.crc32(g.encode()))   # stable per-gene null permutation
        except Exception:
            continue
        if len(df):
            out.append(df)
    return out


_FAM_CTX = "mirna_hallmark/output/learned/family_context.tsv"


def _attach_family_context(df: pd.DataFrame) -> pd.DataFrame:
    """MH-101 (c) → (a): attach the per-family DOSE-WEIGHTED host-exposure map (`family_context`) to each exclusion row,
    so the CN-admission table carries WHICH arm carries any host confound (dose-dominant) and how host-exposed the family
    is — the arm-source attribution that ties the CN gate (a) to (c)/MH-100 (e.g. a miR-17~92 down-weight is tagged
    `fam_dose_cotx_frac`=0.88 / `fam_dominant_arm`=miR-93-5p@MCM7). INTERPRETATION layer (annotation), NOT an estimator
    change — the measured `pleio_down_weight` already captures the effect; this just labels its family context."""
    import os
    if not os.path.exists(_FAM_CTX) or "family" not in df.columns:
        return df
    fc = pd.read_csv(_FAM_CTX, sep="\t")[
        ["family", "dose_cotranscribed_frac", "dominant_arm", "dominant_host", "context_homogeneous"]].rename(
        columns={"dose_cotranscribed_frac": "fam_dose_cotx_frac", "dominant_arm": "fam_dominant_arm",
                 "dominant_host": "fam_dominant_host", "context_homogeneous": "fam_context_homogeneous"})
    return df.merge(fc, on="family", how="left")


_DOSE_ATTR_OUT = "mirna_hallmark/output/learned/instrument/host_dose_attribution.tsv"


def dose_attribute_host_downweights(exclusion_df=None, *, run_delta: bool = True,
                                    excl_path: str = "mirna_hallmark/output/learned/instrument/exclusion.tsv",
                                    out: str = _DOSE_ATTR_OUT) -> pd.DataFrame:
    """**MH-101 (a)↔(c)↔δ-pooling CONNECTION — estimator-grade DOSE-ATTRIBUTION of the CN host down-weights.**
    Closes parked thread (4) ("the confound follows the dose-dominant arm") by wiring `family_context` (c, abundance)
    and `within_family.delta_pooling` (MH-98, estimator) onto every host down-weight from `exclusion(host=True)`.

    **Deliberately does NOT change the down-weight ARITHMETIC — that is already dose-gated, latently.** The first stage
    `gamma_s` (`Z_s → X_fam`) regresses the host-locus CN on the POOLED family dose, so a dose-minor host weakly moves
    X_fam ⇒ low `F_s` ⇒ `weak`/low `beta_weight`; and `pleio_down_weight` is correctly the dose-INDEPENDENT fraction of
    δ_s the host protein explains. Scaling the down-weight by dose-share would UNDER-correct a real confound on an
    admitted segment. So this layer makes the leverage CHECKABLE, fusing three dose readouts:
      • `host_beta_share` = Σ`beta_weight`(host segments) / Σ`beta_weight`(all strong segments) — the EXACT
        estimator-grade BETWEEN-segment dose-leverage of the host locus, straight from the exclusion table (no re-fit);
        `host_is_dose_dominant` = share ≥ 0.5.
      • `abund_dominant_arm` = family_context's ABUNDANCE dose-dominant arm; `abund_host_dominant` = is a host arm it.
      • `est_dominant_arm`/`est_host_share` = δ-pooling's fused dose-share — resolves WITHIN-segment collinear arms via
        chimeric (the one thing `beta_weight` can't, since collinear arms share a segment). `run_delta=False` skips it.
    A down-weight on a dose-DOMINANT host = real leverage (correct to down-weight); on a dose-MINOR host = conservative /
    low-leverage (the family's identification lives on another locus, already reflected by `beta_weight`). The
    abundance↔estimator agreement (`abund_est_agree`) validates family_context's abundance dominant-arm. Deliverable:
    `host_dose_attribution.tsv`. NOT an estimator change (annotation/validation)."""
    import os
    from pathlib import Path
    if exclusion_df is None:
        if not os.path.exists(excl_path):
            return pd.DataFrame()
        exclusion_df = pd.read_csv(excl_path, sep="\t")
    df = exclusion_df.copy()
    if "pleio_sources" not in df.columns or "beta_weight" not in df.columns:
        return pd.DataFrame()
    df["pleio_sources"] = df["pleio_sources"].fillna("")
    df["_is_host_seg"] = (pd.to_numeric(df["pleio_down_weight"], errors="coerce").fillna(0) > 0) & \
                         df["pleio_sources"].str.contains(r"\*", regex=True)                       # host tagged `*`
    host_fams = df.loc[df["_is_host_seg"], ["gene", "family"]].drop_duplicates()
    if not len(host_fams):
        return pd.DataFrame()
    lch = _locus_coding_host()                                                                     # MI* locus → coding host (per-locus)
    fam_dom = {}
    if os.path.exists(_FAM_CTX):
        fc = pd.read_csv(_FAM_CTX, sep="\t")
        fam_dom = dict(zip(fc["family"], fc["dominant_arm"]))
    rows = []
    for _, hf in host_fams.iterrows():
        gene, fam = hf["gene"], hf["family"]
        sub = df[(df["gene"] == gene) & (df["family"] == fam)]
        strong = sub[sub["seg_flag"] != "weak"]
        tot_bw = float(pd.to_numeric(strong["beta_weight"], errors="coerce").fillna(0).sum())
        hseg = sub[sub["_is_host_seg"] & (sub["seg_flag"] != "weak")]
        host_bw = float(pd.to_numeric(hseg["beta_weight"], errors="coerce").fillna(0).sum())
        host_share = host_bw / tot_bw if tot_bw > 0 else 0.0
        host_genes = set()                                                                         # starred host genes
        for s in sub.loc[sub["_is_host_seg"], "pleio_sources"]:
            host_genes |= {p.strip().rstrip("*") for p in str(s).split(",") if p.strip().endswith("*")}
        members = []
        for sm in sub["seg_members"]:
            members += [m.strip() for m in str(sm).split(",") if m.strip()]
        members = list(dict.fromkeys(members))                                                     # short (hsa- stripped)
        host_arms = [m for m in members                                                            # per-locus: arm whose ACTIVE locus's host is starred
                     if any(lch.get(l) in host_genes for l in _arm_active_loci("hsa-" + m))]
        abund_dom = str(fam_dom.get(fam, "")).replace("hsa-", "")
        max_dw = float(pd.to_numeric(sub.loc[sub["_is_host_seg"], "pleio_down_weight"], errors="coerce").max())
        est_dom = ""; est_host_share = np.nan; est_host_dom = None
        if run_delta and len(members) >= 2:
            try:
                from mirna_hallmark.learned import within_family as WF
                dp, _ = WF.delta_pooling(gene, fam)
                dp = dp.dropna(subset=["fused_share"]) if "fused_share" in dp.columns else pd.DataFrame()
                if len(dp):
                    est_dom = str(dp.iloc[0]["arm"]).replace("hsa-", "")
                    est_host_share = float(dp[dp["arm"].str.replace("hsa-", "", regex=False).isin(host_arms)]["fused_share"].sum())
                    est_host_dom = est_dom in host_arms
            except Exception:
                pass
        rows.append({"gene": gene, "family": fam, "n_strong_seg": int(len(strong)),
                     "host_arm": ",".join(host_arms), "host_gene": ",".join(sorted(host_genes)),
                     "host_down_weight": round(max_dw, 2),
                     "host_beta_share": round(host_share, 3), "host_is_dose_dominant": host_share >= 0.5,
                     "abund_dominant_arm": abund_dom, "abund_host_dominant": abund_dom in host_arms,
                     "est_dominant_arm": est_dom, "est_host_share": round(est_host_share, 3) if est_host_share == est_host_share else np.nan,
                     "est_host_dominant": est_host_dom,
                     "abund_est_agree": (est_dom == abund_dom) if est_dom else None})
    res = pd.DataFrame(rows).sort_values(["host_beta_share", "host_down_weight"], ascending=False).reset_index(drop=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    return res


_HOST_LENS_TSV = "mirna_hallmark/output/learned/host_confound_lens.tsv"
_GRADE_OUT = "mirna_hallmark/output/learned/instrument/host_downweight_grades.tsv"


def grade_host_downweights(exclusion_df=None, *, excl_path: str = "mirna_hallmark/output/learned/instrument/exclusion.tsv",
                           lens_path: str = _HOST_LENS_TSV, out: str = _GRADE_OUT) -> pd.DataFrame:
    """**Grade the ambiguous `repression_direction` host down-weights (MH-101 tidy-up, parked thread 3).** A host
    down-weight tagged `clean_antirepression` is unambiguous (δ_s same-sign as γ ⇒ can't be a lost-miR leak ⇒ genuine
    host pleiotropy). A `repression_direction` one is AMBIGUOUS: δ_s opposite-sign to γ ⇒ the down-weight could be a
    real host confound OR an artifact of residual-miR leak / a co-targeting family. This resolves each ambiguous one by
    joining the **`host_confound_lens` OOF verdict** (does conditioning on the SPECIFIC host HELP or HURT held-out
    coupling — the non-circular 2×2 from `prolif_verdict`, already miR-partial-gated via `rel`), keyed on
    (target, host, host-arm). Graded on the OOF `dC` MAGNITUDE (|dC|≤0.02 ⇒ immaterial, the weak-channel default):
      • `dC > +0.02` (host materially de-confounds) ⇒ **`confound_confirmed`** — the down-weight is JUSTIFIED;
      • `dC < −0.02` (conditioning on the host materially HURTS held-out coupling — host is the mechanism / the edge is
        real) ⇒ **`over_control`** — the down-weight is discounting a real edge (candidate to RELAX);
      • `|dC| ≤ 0.02` ⇒ **`host_immaterial`** (host control neither helps nor hurts — common on this weak channel);
      • lens `ambiguous` / no row ⇒ **`unresolved`**; `clean_antirepression` stays **`clean_pleiotropy`**.
    Cheap (no re-fit — reads the two persisted tables). Deliverable: `host_downweight_grades.tsv`. RESULT (2026-07-11, 26
    repression_direction): confound_confirmed 6 · over_control 6 · host_immaterial 11 · unresolved 3 (+10 clean)."""
    import os
    if exclusion_df is None:
        if not os.path.exists(excl_path):
            return pd.DataFrame()
        exclusion_df = pd.read_csv(excl_path, sep="\t")
    df = exclusion_df.copy()
    df["pleio_sources"] = df.get("pleio_sources", "").fillna("")
    hd = df[(pd.to_numeric(df["pleio_down_weight"], errors="coerce").fillna(0) > 0)
            & (df["seg_flag"] != "weak") & df["pleio_sources"].str.contains(r"\*", regex=True)].copy()
    if not len(hd):
        return pd.DataFrame()
    lens = pd.read_csv(lens_path, sep="\t") if os.path.exists(lens_path) else pd.DataFrame()
    lv = {}                                                                # (target, host, arm-short) → (edge_verdict, dC)
    if len(lens):
        for _, r in lens.iterrows():
            lv[(str(r["target"]), str(r["host"]), str(r["arm"]).replace("hsa-", ""))] = (r.get("edge_verdict"), r.get("dC"))
    lch = _locus_coding_host()
    eps = 0.02                                                             # |dC| below this = host-immaterial to held-out coupling
    rows = []
    for _, r in hd.iterrows():
        gene, conf = r["gene"], str(r.get("pleio_confidence", ""))
        hosts = {p.strip().rstrip("*") for p in str(r["pleio_sources"]).split(",") if p.strip().endswith("*")}
        members = [m.strip() for m in str(r["seg_members"]).split(",") if m.strip()]
        for h in sorted(hosts):
            harm = next((m for m in members if any(lch.get(l) == h for l in _arm_active_loci("hsa-" + m))), members[0] if members else "")
            ev, dc = lv.get((str(gene), h, harm), (None, np.nan))
            has_dc = (dc == dc and dc is not None)
            if conf == "clean_antirepression":                             # δ_s same-sign as γ ⇒ already unambiguous host pleiotropy
                grade = "clean_pleiotropy"
            elif not has_dc or ev == "ambiguous":                          # no OOF signal / lens couldn't resolve
                grade = "unresolved"
            elif float(dc) > eps:                                          # host MATERIALLY de-confounds ⇒ down-weight JUSTIFIED
                grade = "confound_confirmed"
            elif float(dc) < -eps:                                         # host control MATERIALLY hurts ⇒ down-weight OVER-CONTROLS a real edge
                grade = "over_control"
            else:                                                          # |dC|≤eps ⇒ host immaterial to coupling (weak-channel default)
                grade = "host_immaterial"
            rows.append({"gene": gene, "family": r["family"], "host_arm": harm, "host_gene": h,
                         "host_down_weight": round(float(r["pleio_down_weight"]), 2),
                         "seg_dose_share": round(float(r.get("seg_dose_share", np.nan)), 3) if r.get("seg_dose_share") == r.get("seg_dose_share") else np.nan,
                         "pleio_confidence": conf, "lens_edge_verdict": ev,
                         "lens_dC": round(float(dc), 3) if dc == dc and dc is not None else np.nan, "grade": grade})
    res = pd.DataFrame(rows).sort_values(["grade", "host_down_weight"], ascending=[True, False]).reset_index(drop=True)
    from pathlib import Path
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    return res


def run_exclusion(genes=None, *, min_members: int = 1, progress: int = 100, coding: bool = False, host: bool = False,
                  workers: int = 1, shuffle_cn: bool = False, shuffle_seed: int = 0,
                  out: str = "mirna_hallmark/output/learned/instrument/exclusion.tsv") -> pd.DataFrame:
    """Batch **scope-3 exclusion gate**: `exclusion` over every seed family (≥`min_members`) of each gene in `genes`
    (default = HE gene universe). Persists per-(gene,family,segment) T1 `δ_s` + over-ID + verdict/weight — the CN
    channel's admission list. `host=True` adds the intronic-host BP down-weight (MH-101); `coding=True` adds the generic
    coding-pleiotropy scan (MH-99, slower — window scan per family). `workers>1` parallelises over genes (fork-COW: the
    parent warms caches, workers inherit). CLI: `python -m mirna_hallmark.learned.instrument --exclusion`."""
    import time
    from pathlib import Path
    from mirna_hallmark.learned import families as FAM
    if genes is None:
        genes = sorted(D.high_evidence_edges()["gene"].dropna().astype(str).unique())
    _EXR.update(min_members=min_members, coding=coding, host=host, shuffle_cn=shuffle_cn, shuffle_seed=shuffle_seed)
    rows, t0 = [], time.time()
    if workers > 1:
        # warm every read-only cache in the PARENT before forking so workers share it copy-on-write
        LD._load(); locus_cn(); _sense_coding_hosts(); _related_host_set(); _locus_coding_host()
        FAM.family_of(["hsa-miR-21-5p"])
        if coding:
            _coding_coords(); _gene_cn(); _mirna_coords()
            _as_np(LD._load()["Y"], "rna_np"); _as_np(_gene_cn(), "gene_cn_np")   # warm numpy views for COW-shared worker access
        from concurrent.futures import ProcessPoolExecutor
        print(f"[exclusion] PARALLEL {len(genes)} genes · {workers} workers · coding={coding} host={host}", flush=True)
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for gi, dfs in enumerate(ex.map(_exclusion_gene, genes, chunksize=4)):
                rows.extend(dfs)
                if progress and gi % progress == 0:
                    print(f"[exclusion] {gi}/{len(genes)} genes, {len(rows)} family-tables, {(time.time()-t0)/60:.1f}min", flush=True)
    else:
        for gi, g in enumerate(genes):
            if progress and gi % progress == 0:
                print(f"[exclusion] {gi}/{len(genes)} genes, {len(rows)} family-tables", flush=True)
            rows.extend(_exclusion_gene(g))
    if not rows:
        return pd.DataFrame()
    out_df = pd.concat(rows, ignore_index=True)
    out_df = _attach_family_context(out_df)                             # MH-101 (c): family dose-weighted host exposure + arm-source
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out, sep="\t", index=False)
    fam_lvl = out_df.drop_duplicates(["gene", "family"])
    strong_seg = out_df[out_df["seg_flag"] != "weak"]
    vc = fam_lvl["family_verdict"].value_counts().to_dict()
    tot_prec = float(out_df["precision"].sum()) or 1.0
    print(f"[exclusion] {len(out_df)} segment-rows | {fam_lvl.shape[0]} (gene,family) | "
          f"{len(strong_seg)} strong segments (F>10) | descriptive: clean {int((strong_seg['seg_flag']=='clean').sum())} "
          f"pleiotropic {int((strong_seg['seg_flag']=='pleiotropic').sum())} | "
          f"CONTINUOUS weight: top-decile precision holds {out_df.nlargest(max(1,len(out_df)//10),'precision')['precision'].sum()/tot_prec*100:.0f}% of total → {out}")
    print(f"[exclusion] NOTE: admission is CONTINUOUS (exclusion_weight = 1/s²_π); seg_flag/family_verdict are DESCRIPTIVE only.")
    return out_df


def _family_copy_shares(mem, X, Cm, yraw, parts, foc) -> dict:
    """Within-family (Kind-A, §8) dose-delivery split by THREE lenses, for comparison (never a repression-identity
    claim — same seed ⇒ identical per-molecule repression). Returns per-member normalized shares + `n`:
      • `abundance` — level share `2^a_m / Σ 2^a` (budget_shift, §5a);
      • `cn_copy`   — the EXOGENOUS leg: first-stage contribution `|γ_m|·sd(CN_m)` of each member's focal-locus CN
        to the family aggregate abundance `X_fam = log2 Σ 2^x` — which genomic COPY delivered the CN-driven dose;
      • `shapley`   — variation share: Shapley of R²(y ~ members | C)."""
    lcn = locus_cn()
    memc = [m for m in mem if foc.get(m) in lcn.columns]
    Xm = X.loc[parts, mem].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    CNm = (np.column_stack([pd.to_numeric(lcn.loc[parts, foc[m]], errors="coerce").to_numpy(float) for m in memc])
           if memc else np.empty((len(parts), 0)))
    mask = np.isfinite(Xm).all(axis=1) & np.isfinite(yraw) & np.isfinite(Cm).all(axis=1)
    if CNm.shape[1]:
        mask &= np.isfinite(CNm).all(axis=1)
    Xm, CNm, Cmk, yk = Xm[mask], CNm[mask], Cm[mask], yraw[mask]
    if mask.sum() < 40:
        u = {m: 1.0 / len(mem) for m in mem}
        return {"abundance": u, "cn_copy": u, "shapley": u, "n": int(mask.sum())}
    lvl = np.median(2.0 ** Xm, axis=0)
    ab = {m: float(lvl[i] / (lvl.sum() or 1.0)) for i, m in enumerate(mem)}
    contrib = {m: 0.0 for m in mem}                                                  # CN-copy: instrument X_fam by member CN
    if CNm.shape[1]:
        xfam = np.log2((2.0 ** Xm).sum(axis=1))
        Zc = np.column_stack([_resid(CNm[:, i], Cmk) for i in range(CNm.shape[1])])
        g = np.linalg.lstsq(Zc, _resid(xfam, Cmk), rcond=None)[0]
        for i, m in enumerate(memc):
            contrib[m] = abs(float(g[i])) * float(np.std(Zc[:, i]))
    ctot = sum(contrib.values()) or 1.0
    cn_copy = {m: contrib[m] / ctot for m in mem}
    yr = _resid(yk, Cmk)                                                             # Shapley-by-variation → target
    Xr = np.column_stack([_resid(Xm[:, i], Cmk) for i in range(Xm.shape[1])])
    phi = _shapley_r2(Xr, yr)
    ptot = float(phi.sum()) or 1.0
    sh = {m: float(phi[i] / ptot) for i, m in enumerate(mem)}
    return {"abundance": ab, "cn_copy": cn_copy, "shapley": sh, "n": int(mask.sum())}


def cluster_instrument(gene: str, *, unit: str = "cluster", window: float = 5e4,
                       min_members: int = 2, attribute: bool = False) -> pd.DataFrame:
    """CN instrument at the **genomic-cluster** level (what CN moves), then break down by identifiability
    (user synthesis, 2026-07-04). `unit='cluster'` groups co-located arms; `unit='family'` groups by seed
    (for contrast — wrong for the instrument, since family clusters have independent CN).

    Per group: family-CN = mean member locus CN; first-stage F (CN→group abundance); reduced-form
    (CN→target | C) < 0 ⇒ the CLUSTER causally represses the gene. Then within-group CN collinearity:
    low ⇒ a member/sub-locus has independent CN and CAN be resolved; high ⇒ irreducibly cluster-level (§F).

    `attribute=True` (family unit, well-instrumented groups) adds the **within-family dose-delivery split** (§8,
    gap 2 within-family): per-member `cn_copy` (exogenous — which copy delivered the CN-driven dose) vs
    `abundance` (level) vs `shapley` (variation), with their rank concordances — to test whether the exogenous
    CN-copy leg agrees with or adds to the observational splits. NEVER a repression-identity claim (same seed)."""
    from collections import defaultdict
    from mirna_hallmark.learned import families as FAM
    lcn = locus_cn()
    Y, X, C, _ = LD.assemble_gene(gene, he_only=False)
    # per-arm focal-locus CN (gap 1): each arm sourced from its dominant hairpin, not the paralog composite
    foc = {a: _arm_focal_locus(a)["locus_id"] for a in X.columns}
    reg_loci = {a: foc[a] for a in X.columns if foc.get(a) in lcn.columns}
    cn = pd.DataFrame({a: lcn[lid] for a, lid in reg_loci.items()})
    parts = X.index.intersection(cn.index)
    Cm = C.loc[parts].to_numpy(float)
    y = _resid(Y.loc[parts].to_numpy(float), Cm)
    regs = [a for a in X.columns if a in cn.columns]
    grp = _genomic_clusters(regs, window) if unit == "cluster" else FAM.family_of(regs).to_dict()
    fam_of = FAM.family_of(regs)
    groups = defaultdict(list)
    for a, gid in grp.items():
        groups[gid].append(a)
    rows = []
    co = _mirna_coords()
    for gid, mem in groups.items():
        mem = [a for a in mem if pd.to_numeric(cn[a], errors="coerce").loc[parts].notna().any()]
        if len(mem) < min_members:
            continue
        cc = cn.loc[parts, mem].apply(pd.to_numeric, errors="coerce").corr(method="spearman").abs()
        maxc = float(cc.where(~np.eye(len(mem), dtype=bool)).max().max())
        base = {"gene": gene, "unit": unit, "chrom": co.loc[mem[0], "chrom"] if mem[0] in co.index else "?",
                "n_members": len(mem), "seed_families": ", ".join(sorted(set(fam_of.reindex(mem).dropna()))),
                "cn_collinearity": round(maxc, 2), "members_resolvable": bool(maxc < 0.7), "members": ", ".join(mem)}
        if unit == "family":                                                         # proper family multi-IV (§7/§8)
            fr = family_multi_iv(gene, mem, foc)
            if not fr:
                continue
            F, beta = fr["first_stage_F"], fr["beta_family"]
            row = {**base, "first_stage_F": round(F, 1) if F == F else F,
                   "beta_family": round(beta, 3) if beta == beta else beta,
                   "cluster_causal": bool(F == F and F > 10 and beta == beta and beta < 0),
                   "n_strong": fr["n_strong"], "n_indep_cn": fr["n_indep_cn"],
                   "overid_p": round(fr["overid_p"], 3) if fr["overid_p"] == fr["overid_p"] else np.nan,
                   "recurse_flag": fr["recurse_flag"]}
            if attribute and F == F and F > 10:                                       # §8 within-family dose-delivery
                sh = _family_copy_shares(mem, X, Cm, Y.loc[parts].to_numpy(float), parts, foc)
                cnc = fr["cn_copy"]                                                   # cn_copy from the SAME estimator
                order = sorted(mem)
                cnv = np.array([cnc[m] for m in order]); abv = np.array([sh["abundance"][m] for m in order])
                shv = np.array([sh["shapley"][m] for m in order])
                _rho = lambda u, v: (float(spearmanr(u, v).correlation) if len(u) > 2 and np.std(u) > 0 and np.std(v) > 0
                                     else np.nan)
                row.update({
                    "member_cn_copy": ";".join(f"{m.replace('hsa-','')}={cnc[m]:.2f}" for m in order),
                    "member_pi": ";".join(f"{m.replace('hsa-','')}={fr['solo'][m]['pi']:+.3f}" for m in order),
                    "member_abundance": ";".join(f"{m.replace('hsa-','')}={sh['abundance'][m]:.2f}" for m in order),
                    "member_shapley": ";".join(f"{m.replace('hsa-','')}={sh['shapley'][m]:.2f}" for m in order),
                    "cn_copy_top": order[int(np.argmax(cnv))].replace("hsa-", ""),
                    "abundance_top": order[int(np.argmax(abv))].replace("hsa-", ""),
                    "shapley_top": order[int(np.argmax(shv))].replace("hsa-", ""),
                    "rho_cncopy_abund": round(_rho(cnv, abv), 2) if _rho(cnv, abv) == _rho(cnv, abv) else np.nan,
                    "rho_cncopy_shapley": round(_rho(cnv, shv), 2) if _rho(cnv, shv) == _rho(cnv, shv) else np.nan,
                    "l1_cncopy_abund": round(float(np.abs(cnv - abv).sum()), 2),
                    "l1_cncopy_shapley": round(float(np.abs(cnv - shv).sum()), 2),
                    "l1_abund_shapley": round(float(np.abs(abv - shv).sum()), 2)})
        else:                                                                        # cluster unit: mean-collapse screen
            xF = X.loc[parts, mem].mean(axis=1).to_numpy(float)
            zF = pd.to_numeric(cn.loc[parts, mem].mean(axis=1), errors="coerce")
            ok = np.isfinite(zF.to_numpy())
            if ok.sum() < 40:
                continue
            F = _first_stage_F(xF[ok], zF.to_numpy(float)[ok], Cm[ok])
            rho_rf = spearmanr(_resid(zF.to_numpy(float)[ok], Cm[ok]), y[ok]).correlation
            row = {**base, "first_stage_F": round(F, 1) if F == F else F,
                   "rho_clusterCN_target": round(rho_rf, 3), "cluster_causal": bool(F > 10 and rho_rf < 0)}
        rows.append(row)
    df = pd.DataFrame(rows)
    return df.sort_values("first_stage_F", ascending=False) if len(df) else df


def focal_locus_concordance(*, force: bool = False,
                            out: str = "mirna_hallmark/output/learned/instrument/focal_locus_concordance.tsv") -> pd.DataFrame:
    """Per-arm CN→expr concordance on the arm's **focal locus** CN (partial-Spearman | CONFOUNDER_NUMERIC),
    the instrument-admission analogue of `mirna_cnv_expr_concordance` but paralog-decontaminated (gap-1
    side-finding). Single-locus arms are bit-identical to the aggregated table (arm CN ≡ locus CN); multi-
    locus arms are scored on their dominant (w_norm) hairpin, un-diluted — which recovers paralog arms the
    aggregated CN buries below threshold (e.g. miR-92a-3p 0.15→0.31, miR-19b-3p 0.13→0.24). Cached to `out`."""
    from pathlib import Path
    from mirna_hallmark import config as C, stats as S
    p = Path(out)
    if p.exists() and not force:
        return pd.read_csv(p, sep="\t")
    expr = D.load_mirna_arms()
    clin = D.load_clinical_strata().set_index("participant")
    conf = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    lcn, amap = locus_cn(), arm_loci_map()
    rows = []
    for arm in expr.index:
        foc = _arm_focal_locus(arm)
        lid = foc["locus_id"]
        if lid is None or lid not in lcn.columns:
            continue
        m = pd.concat([expr.loc[arm].rename("expr"), lcn[lid].rename("cn")], axis=1).dropna()
        if len(m) < 40:
            continue
        st = S.correlation_pair(m["expr"], m["cn"], clin.reindex(m.index)[conf] if conf else None)
        rows.append({"arm": arm, "focal_locus": lid, "n_loci": foc["n_loci"], "multilocus": foc["multilocus"],
                     "n": int(st["n"]), "partial_rho": float(st["partial_rho"]),
                     "partial_p": float(st["partial_p"]), "spearman_rho": float(st["spearman_rho"])})
    df = pd.DataFrame(rows).sort_values("partial_rho", ascending=False)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep="\t", index=False)
    return df


def locus_origin(arm: str) -> pd.DataFrame:
    """Per-locus FIRST-STAGE / **locus-of-origin** assay for a mature arm — **gene-independent** (the mature
    product is locus-agnostic, §6; measured once per (arm,locus), never per edge). For each contributing
    hairpin: does its CN move the mature-arm abundance (partial-Spearman + first-stage F | CONFOUNDER_NUMERIC)?
    Separates AVAILABILITY (CN variation `frac_nondiploid`) from ORIGIN (first-stage F/ρ). Classifies:
      • `active_source` — CN varies AND drives the mature arm (F>10, ρ>0): valid instrument + transcriptional source;
      • `silent_source` — CN varies but does NOT drive the mature arm ⇒ the arm is not (meaningfully) transcribed
        from this copy (e.g. miR-92a/19b **chrX 106a~363** vs **chr13 17~92**) — a positive locus-of-origin finding,
        and NOT a valid instrument (using it violates relevance, and its target-associations are pleiotropy);
      • `cn_invariant` — too little CN variation to instrument at all (an uninformative weak instrument)."""
    key = ("origin", arm)
    if key in _CACHE:
        return _CACHE[key]
    from mirna_hallmark import config as C, stats as S  # noqa: F401
    if "lo_data" not in _CACHE:                               # CACHE the exact loaders ONCE (was a fresh load_mirna_arms +
        _CACHE["lo_data"] = (D.load_mirna_arms(),             # clinical read PER new-arm call — byte-identical, just cached)
                             D.load_clinical_strata().set_index("participant"))
    expr, clin = _CACHE["lo_data"]
    conf = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    lcn, g = locus_cn(), arm_loci_map().get(arm)
    if g is None or arm not in expr.index:
        _CACHE[key] = pd.DataFrame(); return _CACHE[key]
    ex, rows = expr.loc[arm], []
    for _, r in g.iterrows():
        L = r["locus_id"]
        if L not in lcn.columns:
            continue
        m = pd.concat([ex.rename("e"), pd.to_numeric(lcn[L], errors="coerce").rename("cn")], axis=1).dropna()
        cov = clin.reindex(m.index)[conf]
        keep = cov.notna().all(axis=1).to_numpy()
        m, Cm = m[keep], cov[keep].to_numpy(float)
        if len(m) < 40:
            continue
        F = _first_stage_F(m["e"].to_numpy(float), m["cn"].to_numpy(float), Cm)
        rho = float(spearmanr(_resid(m["cn"].to_numpy(float), Cm), _resid(m["e"].to_numpy(float), Cm)).correlation)
        frac = float((np.round(m["cn"].to_numpy()) != 2).mean())
        status = ("cn_invariant" if frac < 0.10
                  else "active_source" if (F == F and F > 10 and rho > 0)
                  else "silent_source")
        rows.append({"arm": arm, "locus_id": L, "chrom": r["chrom"], "w_norm": round(float(r["w_norm"]), 3),
                     "n": len(m), "frac_nondiploid": round(frac, 3),
                     "first_stage_F": round(F, 1) if F == F else F, "first_stage_rho": round(rho, 3),
                     "status": status})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("first_stage_F", ascending=False)
    _CACHE[key] = df
    return df


def multi_iv(arm: str, gene: str, *, min_n: int = 60) -> dict:
    """Multi-instrument 2SLS + **Sargan over-ID** for a multi-locus arm (CN_INSTRUMENT.md §5, gap 2A). Instruments
    = the arm's **active-source** loci only (`locus_origin`): a silent-source copy doesn't perturb the mature dose,
    so it is not a valid instrument (relevance fails). On the active sources: first-stage `x ~ γ₁CN_L1+…+C`,
    reduced-form `y ~ π₁CN_L1+…+C`, 2SLS β (dose→target). With ≥2 active sources the Wald ratios must agree (same
    molecule) — a **Hansen J** over-ID test (heteroskedasticity-robust) tests it; rejection ⇒ a source-locus is
    pleiotropic (its co-amplified cluster also hits `gene`) → `recurse_flag` (§5/§7). Excluded **silent-source**
    loci whose CN still moves the target are reported
    as `silent_pleiotropy` (horizontal-pleiotropy QC — exactly why they're not instruments). OLS/2SLS (linear)."""
    from scipy.stats import chi2
    Y, X, C, _ = LD.assemble_gene(gene, he_only=False)
    if arm not in X.columns:
        raise KeyError(f"{arm} not a candidate regulator of {gene}")
    origin = locus_origin(arm)
    if origin.empty:
        raise KeyError(f"{arm} has no locus-origin assay")
    lcn = locus_cn()
    active = [L for L in origin.loc[origin["status"] == "active_source", "locus_id"] if L in lcn.columns]
    silent = [L for L in origin.loc[origin["status"] == "silent_source", "locus_id"] if L in lcn.columns]
    if not active:
        raise ValueError(f"{arm}: no active-source locus (origins: {origin['status'].tolist()})")
    parts = X.index.intersection(lcn.index)
    Z = lcn.loc[parts, active].apply(pd.to_numeric, errors="coerce")
    x, y, Cd = X.loc[parts, arm], Y.loc[parts], C.loc[parts]
    ok = Z.notna().all(axis=1) & x.notna() & y.notna() & Cd.notna().all(axis=1)
    Z, x, y, Cm = Z[ok].to_numpy(float), x[ok].to_numpy(float), y[ok].to_numpy(float), Cd[ok].to_numpy(float)
    n, ka = len(x), Z.shape[1]
    if n < min_n:
        raise ValueError(f"insufficient data (n={n})")
    xr, yr = _resid(x, Cm), _resid(y, Cm)                                          # FWL: partial out C
    Zr = np.column_stack([_resid(Z[:, j], Cm) for j in range(ka)])
    tss = float((xr ** 2).sum())
    solo = {}                                                                       # per active-source Wald ratio
    for j, L in enumerate(active):
        zj = Zr[:, j]
        gj, pj = float((zj @ xr) / (zj @ zj)), float((zj @ yr) / (zj @ zj))
        s_ssr = float(((xr - zj * gj) ** 2).sum())
        solo[L] = {"F": float(((tss - s_ssr) / 1) / (s_ssr / (n - 1))) if s_ssr > 0 else np.nan,
                   "wald": (pj / gj) if abs(gj) > 1e-12 else np.nan}
    g = np.linalg.lstsq(Zr, xr, rcond=None)[0]                                      # joint 2SLS on active sources
    xhat = Zr @ g
    ssr = float(((xr - xhat) ** 2).sum())
    F = float(((tss - ssr) / ka) / (ssr / (n - ka))) if ssr > 0 and n > ka else np.nan
    beta = float((xhat @ yr) / (xhat @ xhat)) if (xhat @ xhat) > 0 else np.nan
    overid_J = overid_p = np.nan
    dev = None
    if ka > 1:                                                                      # Hansen J (heteroskedasticity-
        e = yr - beta * xr                                                          # robust over-ID) ~ χ²(ka-1)
        gmom = Zr.T @ e / n                                                         # sample moment vector (H0: →0)
        S = (Zr * (e ** 2)[:, None]).T @ Zr / n                                     # robust GMM weight (not σ²·Z'Z)
        try:
            J = float(n * gmom @ np.linalg.solve(S, gmom))
            overid_J, overid_p = J, float(chi2.sf(J, ka - 1))
        except np.linalg.LinAlgError:
            pass
        dev = max(active, key=lambda L: abs((solo[L]["wald"] if solo[L]["wald"] == solo[L]["wald"] else 0) - beta))
    pleio = []                                                                       # silent-source horizontal pleiotropy
    for L in silent:
        zz = pd.to_numeric(lcn.loc[parts, L], errors="coerce")[ok].to_numpy(float)
        okz = np.isfinite(zz)
        if okz.sum() < min_n:
            continue
        rr = spearmanr(_resid(zz[okz], Cm[okz]), _resid(y[okz], Cm[okz]))
        if rr.pvalue < 0.05:
            pleio.append(f"{L}={rr.correlation:+.2f}")
    powered = bool(ka > 1)
    return {"arm": arm, "gene": gene, "n": n, "n_active_source": ka, "n_silent_source": len(silent),
            "active_loci": ",".join(active), "silent_loci": ",".join(silent),
            "first_stage_F": F, "beta_2sls": beta,
            "wald_by_locus": ";".join(f"{L}={solo[L]['wald']:+.3f}(F{solo[L]['F']:.0f})" for L in active),
            "overid_J": overid_J, "overid_p": overid_p, "over_id_powered": powered,
            "over_id_ok": (bool(overid_p > 0.05) if overid_p == overid_p else np.nan),
            "deviating_locus": dev, "recurse_flag": bool(powered and overid_p == overid_p and overid_p < 0.05),
            "silent_pleiotropy": ";".join(pleio), "has_silent_pleiotropy": bool(pleio),
            "multilocus": bool(len(origin) > 1)}


def run_family_copy_attribution(genes=None, *, min_members: int = 2,
                                out: str = "mirna_hallmark/output/learned/instrument/family_copy_attribution.tsv") -> pd.DataFrame:
    """Gap 2 within-family (§8): run `cluster_instrument(unit='family', attribute=True)` across `genes` and
    quantify whether the exogenous **CN-copy** dose-delivery split agrees with / adds to the observational
    **abundance** and **shapley** splits — reported separately for CN-RESOLVABLE families (members with
    independent CN, where CN-copy CAN split) vs non-resolvable (where CN correctly abstains ~uniform)."""
    from pathlib import Path
    if genes is None:
        genes = sorted(D.high_evidence_edges()["gene"].unique())
    rows = []
    for g in genes:
        try:
            df = cluster_instrument(g, unit="family", attribute=True, min_members=min_members)
        except Exception:
            continue
        if "member_cn_copy" in df.columns:
            rows.append(df[df["member_cn_copy"].notna()])
    if not rows:
        return pd.DataFrame()
    fam = pd.concat(rows, ignore_index=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fam.to_csv(out, sep="\t", index=False)
    res = fam[fam["members_resolvable"]]                                            # independent-CN members
    nonres = fam[~fam["members_resolvable"]]
    print(f"attributed families: {len(fam)} over {len(genes)} genes | CN-resolvable (independent CN): {len(res)} | "
          f"non-resolvable (shared CN): {len(nonres)}")
    print("  [top-match + median L1 distance between share vectors; L1∈[0,2], 0=identical splits]")
    for lab, sub in [("CN-resolvable", res), ("non-resolvable", nonres)]:
        if not len(sub):
            continue
        print(f"  {lab} (n={len(sub)}): cn_top==ab_top {int((sub['cn_copy_top']==sub['abundance_top']).sum())}/{len(sub)} | "
              f"cn_top==shap_top {int((sub['cn_copy_top']==sub['shapley_top']).sum())}/{len(sub)} | "
              f"L1(cn,ab)={sub['l1_cncopy_abund'].median():.2f} L1(cn,shap)={sub['l1_cncopy_shapley'].median():.2f} "
              f"L1(ab,shap)={sub['l1_abund_shapley'].median():.2f}")
    print(f"  → {out}")
    return fam


def well_instrumented_arms(min_rho: float = 0.2, min_n: int = 500, *, source: str = "focal_locus") -> list[str]:
    """Arms whose CN is a strong-enough first-stage instrument. `source='focal_locus'` (default, gap-1) scores
    each arm on its focal-hairpin CN so paralog arms aren't diluted out; `source='aggregated'` = the legacy
    mature-aggregated concordance table."""
    if source == "focal_locus":
        c = focal_locus_concordance()
    else:
        c = pd.read_csv("mirna_hallmark/output/mirna_locus_cnv/mirna_cnv_expr_concordance.tsv", sep="\t")
    sub = c[(c["partial_rho"] > min_rho) & (c["n"] >= min_n)]
    return sub.sort_values("partial_rho", ascending=False)["arm"].tolist()


def run(n_edges: int = 12) -> pd.DataFrame:
    """Report the instrument on HE edges whose regulator arm is well-instrumented."""
    arms = set(well_instrumented_arms())
    edges = D.high_evidence_edges()
    cand = edges[edges["miRNA"].isin(arms)][["miRNA", "gene"]].drop_duplicates()
    rows = []
    for _, r in cand.head(200).iterrows():
        try:
            rows.append(edge_instrument(r["miRNA"], r["gene"]))
        except Exception:
            continue
        if len(rows) >= n_edges:
            break
    df = pd.DataFrame(rows)
    if len(df):
        with pd.option_context("display.width", 170):
            print(df.round(3).to_string(index=False))
        u = df[df["usable"]]
        print(f"\nusable instruments (F>10): {len(u)}/{len(df)} | "
              f"causal-concordant (reduced<0 AND obs<0): {int(df['causal_concordant'].sum())}")
    return df


def run_multi_iv(*, min_rho: float = 0.2, min_n: int = 500,
                 out: str = "mirna_hallmark/output/learned/instrument/cn_multi_iv.tsv") -> pd.DataFrame:
    """Gap 2A: multi-IV + Sargan over-ID over every HE edge of a **multi-locus** well-instrumented arm. Each
    arm's contributing hairpin loci are separate instruments; the estimate is 2SLS on the strong ones and, where
    ≥2 are strong, the over-ID test flags loci whose Wald ratio disagrees (a co-amplified-cluster effect, not the
    single arm) for recursion (§5/§7). Empirically most edges reduce to the one strong (focal) locus — the value
    is the over-ID DIAGNOSIS where both paralog loci are live."""
    from pathlib import Path
    foc = focal_locus_concordance()
    ml = foc[foc["multilocus"] & (foc["partial_rho"] > min_rho) & (foc["n"] >= min_n)]["arm"].tolist()
    edges = D.high_evidence_edges()
    rows = []
    for a in ml:
        for g in edges.loc[edges["miRNA"] == a, "gene"].unique():
            try:
                r = multi_iv(a, g)
            except Exception:
                continue
            rows.append(r)
    df = pd.DataFrame(rows)
    if len(df):
        pw = df[df["over_id_powered"]]
        if len(pw):
            df["overid_q_by"] = np.nan
            df.loc[df["over_id_powered"], "overid_q_by"] = benjamini_yekutieli(pw["overid_p"].values)
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False)
        print(f"multi-locus arms: {ml} | edges scored: {len(df)} | "
              f"single active-source (reduces to it): {int((df['n_active_source'] == 1).sum())} | "
              f"≥2 active (over-ID powered): {int(df['over_id_powered'].sum())} | "
              f"over-ID REJECT → recurse: {int(df['recurse_flag'].sum())} | "
              f"edges w/ silent-source pleiotropy: {int(df['has_silent_pleiotropy'].sum())}  → {out}")
    return df


def run_clean(n_edges: int = 60, *, window: float = 1e6,
              out: str = "mirna_hallmark/output/learned/instrument/cn_instrument_clean.tsv") -> pd.DataFrame:
    """CN instrument + cluster-cleanness + within-cluster attribution (user 2026-07-05). The **strong causal**
    set = causal-concordant AND (cluster-clean OR `arm` uniquely carries the within-cluster anti-corr) — i.e.
    the reduced-form effect can be attributed to `arm`, not just its genomic cluster."""
    from pathlib import Path
    arms = set(well_instrumented_arms())
    edges = D.high_evidence_edges()
    cand = edges[edges["miRNA"].isin(arms)][["miRNA", "gene"]].drop_duplicates()
    rows = []
    for _, r in cand.head(400).iterrows():
        try:
            e = edge_instrument(r["miRNA"], r["gene"])
            a = cluster_attribution(r["miRNA"], r["gene"], window=window)
        except Exception:
            continue
        e.update({k: a[k] for k in ("n_cotargeters", "cluster_clean", "arm_unique_anticorr",
                                    "cluster_attributable", "cotargeter_partials")})
        e["strong_causal"] = bool(e["causal_concordant"] and e["cluster_attributable"])
        rows.append(e)
        if len(rows) >= n_edges:
            break
    df = pd.DataFrame(rows)
    if len(df):
        df["q_reduced_by"] = benjamini_yekutieli(df["p_reduced"].values)   # dependence-robust FDR on the causal (reduced-form) test
        df["strong_causal_fdr"] = df["strong_causal"] & (df["q_reduced_by"] < 0.05)
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False)
        cc = df["causal_concordant"]
        print(f"usable F>10: {int(df['usable'].sum())}/{len(df)} | causal-concordant: {int(cc.sum())} | "
              f"cluster-clean: {int(df['cluster_clean'].sum())} | concordant&clean: {int((cc & df['cluster_clean']).sum())} | "
              f"concordant&pleiotropic-but-arm-unique: {int((cc & ~df['cluster_clean'] & df['arm_unique_anticorr']).sum())} | "
              f"STRONG causal: {int(df['strong_causal'].sum())} | STRONG & FDR-sig (q_BY<0.05): {int(df['strong_causal_fdr'].sum())}  → {out}")
        with pd.option_context("display.width", 200):
            print(df[cc][["arm", "gene", "rho_reduced", "n_cotargeters", "cluster_clean",
                          "arm_unique_anticorr", "strong_causal", "cotargeter_partials"]].round(3).to_string(index=False))
    return df


if __name__ == "__main__":
    _argv = __import__("sys").argv
    if "--ring1" in _argv:                                    # batch Ring-1 within-family causal attribution
        run_family_causal_attribution()
    elif "--gap2b" in _argv:                                  # batch Gap-2B between-family exclusion (condition p_S)
        run_between_family_exclusion()
    elif "--exclusion" in _argv:                              # scope-3 CN-instrument exclusion gate (over-ID + T1 δ_s)
        run_exclusion()
    elif "--multi-iv" in _argv:                               # gap-2A multi-locus over-ID table
        run_multi_iv()
    elif "--clean" in _argv:
        run_clean()
    else:
        run()
