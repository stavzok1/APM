"""CN channel builder — the multi-family segment terms for `spike_slab._gibbs_posterior(channels=)`.

CHANNEL_FUSION_DESIGN §1/§2: a CN segment (co-amplified genomic locus) `s` is one exogenous shock whose reduced
form on the target Y is `π̂_s = Σ_{f on s} γ_{f,s}·β_f` (sum over the co-located FAMILIES that regulate g). The
Gibbs term `N(π̂_s; Σ_f γ_{f,s} β_f, s²_π)` constrains that LINEAR COMBINATION, so when ≥2 co-located families share
a segment the posterior can only split π̂_s between their β's by trading them off → an induced **cross-family posterior
covariance** = the between-family exclusion, graded with width (replaces MH-88's point Shapley).

This unifies the two cases:
  • 1-member term  = the MH-91 single-family MVP (one family owns the segment).
  • ≥2-member term = the cross-family covariance case (co-located other-seed families on the same amplicon).

`multi_family_terms(gene, cols, members, ...)` returns the `channels=` list keyed to the dense fit's column order
(`cols` = family names in fit order; `members` = family→arms from `families.collapse_by_family`). Coding-gene
pleiotropy (H2/H3) and δ-pooling (F2) are separate follow-ups — this builder is the family-β covariance core.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as IN


def _resid(V: np.ndarray, Cd: np.ndarray) -> np.ndarray:
    """Residualise every column of V on [1,C]."""
    return V - Cd @ np.linalg.lstsq(Cd, V, rcond=None)[0]


def _z(v: np.ndarray) -> np.ndarray:
    sd = float(v.std())
    return v / sd if sd > 1e-12 else v


def multi_family_terms(gene: str, cols, members: dict, *, assembled=None,
                       window: float = 5e4, min_F: float = 10.0, min_n: int = 60,
                       shuffle_cn: bool = False, seed: int = 0) -> list[dict]:
    """Build the CN channel `channels=` list for `gene`.

    cols     : family names in the dense fit's COLUMN order (channel `members` index into this).
    members  : family → [arms]  (from `families.collapse_by_family`).
    Groups all families' member arms into genomic clusters (`_genomic_clusters`, single-linkage `window` bp) — each
    cluster = one shared CN segment Z_c. Per cluster carrying ≥1 family: γ_{f,c} = Z_c→X_fam,f, π̂_c = Z_c→Y, and
    s²_π = se(π̂_c)² + δ_c² (δ_c = the joint direct effect of Z_c on Y not through any family on c = the pleiotropy
    floor, §2.2). Emits one term per cluster with F(Z_c→any family)>min_F.
    """
    asm = assembled if assembled is not None else LD.assemble_gene(gene, he_only=False)
    Y, X, C, _ = asm
    lcn = IN.locus_cn()
    col_ix = {f: i for i, f in enumerate(cols)}

    # arm → focal locus, keep families whose members have a usable locus
    fam_arms: dict = {}
    for f in cols:
        arms = [m for m in members.get(f, []) if m in X.columns]
        loci = {m: IN._arm_focal_locus(m)["locus_id"] for m in arms}
        arms = [m for m in arms if loci.get(m) in lcn.columns]
        if arms:
            fam_arms[f] = arms
    if not fam_arms:
        return []

    all_arms = sorted({m for arms in fam_arms.values() for m in arms})
    clusters = IN._genomic_clusters(all_arms, window=window)         # arm → genomic cluster id
    # cluster → the participating (family, its arms in this cluster)
    clu_fams: dict = {}
    for f, arms in fam_arms.items():
        by_c: dict = {}
        for m in arms:
            by_c.setdefault(clusters.get(m), []).append(m)
        for c, ms in by_c.items():
            if c is not None:
                clu_fams.setdefault(c, {})[f] = ms

    parts = X.index.intersection(lcn.index)
    Cd = np.column_stack([np.ones(len(parts)), C.loc[parts].apply(pd.to_numeric, errors="coerce").to_numpy(float)])
    y = pd.to_numeric(Y.loc[parts], errors="coerce").to_numpy(float)
    rng = np.random.default_rng(seed)
    terms: list[dict] = []

    for c, fam_ms in clu_fams.items():
        # shared segment CN Z_c = mean locus CN over all arms in the cluster (co-amplified → move together)
        loci = sorted({IN._arm_focal_locus(m)["locus_id"] for ms in fam_ms.values() for m in ms})
        Zc = lcn.loc[parts, [l for l in loci if l in lcn.columns]].apply(pd.to_numeric, errors="coerce").mean(1).to_numpy(float)
        # family pooled dose X_fam,f (true RPM pool, re-logged) for each family on this cluster (ALL its arms, not just in-cluster)
        Xf = {}
        for f in fam_ms:
            am = fam_arms[f]
            xr = X.loc[parts, am].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            Xf[f] = np.log2((2.0 ** xr).sum(1) + 1.0)
        ok = np.isfinite(y) & np.isfinite(Zc) & np.all([np.isfinite(v) for v in Xf.values()], 0) & np.isfinite(Cd).all(1)
        if ok.sum() < min_n:
            continue
        Cd_o, y_o, Zc_o = Cd[ok], y[ok], Zc[ok]
        if shuffle_cn:
            Zc_o = rng.permutation(Zc_o)                             # shuffled-CN null (§3.4): break the instrument
        Xmat = np.column_stack([Xf[f][ok] for f in fam_ms])
        # GAUGE (§E1): match the dense fit EXACTLY — r = -resid(Y|C) RAW (NOT z-scored, as `_prep`); Z and each family
        # dose z-scored (as the dense `Xz`). π̂ on the raw-r scale ⇒ Σγβ = π̂ holds with the dense β (verified: IV-β=π̂/γ
        # lands on the dense scale). z-scoring r (an earlier bug) put π̂ ~6× off and silently killed the channel.
        R = _resid(np.column_stack([y_o, Zc_o, Xmat]), Cd_o)
        rr = -R[:, 0]                                             # raw dense response r = -resid(Y|C)
        zr = _z(R[:, 1])                                          # Z_c z-scored
        Xr = np.column_stack([_z(R[:, 2 + j]) for j in range(R.shape[1] - 2)])   # each family dose z-scored
        zz = float(zr @ zr) or 1.0
        # first stage γ_{f,c} = Z_c → X_fam,f  and its F (strongest family on the segment)
        gammas, Fs = {}, []
        for j, f in enumerate(fam_ms):
            a = float((zr @ Xr[:, j]) / zz)
            rss = float(((Xr[:, j] - a * zr) ** 2).sum()); dof = ok.sum() - Cd.shape[1] - 1
            Fj = float((float(Xr[:, j] @ Xr[:, j]) - rss) / (rss / dof)) if rss > 0 and dof > 0 else 0.0
            gammas[f] = a; Fs.append(Fj)
        if not Fs or max(Fs) < min_F:                                # weak-instrument cut (F>10 on the strongest family)
            continue
        # segment reduced form π̂_c = Z_c → r  + its SE  (r scale = dense fit)
        pihat = float((zr @ rr) / zz)
        rss_y = float(((rr - pihat * zr) ** 2).sum()); dofy = ok.sum() - Cd.shape[1] - 1
        se_pi = float(np.sqrt((rss_y / dofy) / zz)) if dofy > 0 else np.nan
        # δ_c = joint direct effect of Z_c on r NOT through ANY family on the segment (pleiotropy floor, §2.2)
        M = np.column_stack([Xr, zr])
        coef, *_ = np.linalg.lstsq(M, rr, rcond=None)
        delta = float(coef[-1])
        s2 = (se_pi ** 2 if se_pi == se_pi else 0.0) + delta ** 2
        mem = [(col_ix[f], gammas[f]) for f in fam_ms if f in col_ix]
        if not mem or s2 <= 0:
            continue
        terms.append({"members": mem, "pihat": pihat, "s2": float(s2),
                      "families": list(fam_ms), "F": round(max(Fs), 1), "delta": round(delta, 3),
                      "n_fam": len(mem)})
    return terms


def between_family_bayes(gene: str, *, min_fam: int = 2, seed: int = 0,
                         n_iter: int = 1600, burn: int = 500, n_boot: int = 20) -> pd.DataFrame:
    """CONTINUOUS-INCLUSION between-family exclusion (G2, MH-94) — the Bayesian replacement for MH-88's binary
    'does the family enter the coalition' gate. Runs the dense sampler in **INCLUSION mode** (evidence-π) with vs
    without the CN channel and reports, per co-located family on a multi-family segment:
      • `pip_mrna` / `pip_cn` = posterior P(family enters) — the CONTINUOUS entry; ΔPIP = the CN existence-evidence.
      • `shap_share` ± `shap_sd` = **Bayesian Shapley** of the segment's reduced form `π=Σγβ`: since π is LINEAR in β,
        Shapley of family f is exact = `γ_f β_f`, so the per-draw share `γ_fβ_f/Σγβ` → posterior mean ± sd (fair
        collinear credit WITH width, vs MH-88's ±0.01 MC point). CN folded in (channel-informed posterior).
    General beyond CN: the same (inclusion z → PIP-entry, per-draw Shapley → width) reframe applies to any Shapley.

    ⚠ **`shap_sd` (the POSTERIOR width) IS NOT THE HONEST WIDTH — use `shap_sd_boot`** (MH-102e, 2026-07-12).
    Re-derived by bootstrapping participants and refitting (`n_boot`, default 20), the posterior width is:
      • **CONSERVATIVE where the share is contested** — PTEN miR-141/200a is 0.735 ± **0.441** posterior vs
        **0.243** bootstrap (1.8× too WIDE). The MH-94 non-identifiability claim therefore STANDS and is if
        anything under-sold. (Two earlier assertions that these widths were *understated* are RETRACTED.)
      • **CATASTROPHICALLY OVERCONFIDENT AT THE BOUNDARY** — when one family's PIP≈1 and its co-located
        families' PIP≈0, every posterior draw gives the same split, so `shap_sd` collapses to **0.000** and the
        segment reports a share of **1.000 ± 0.000**. The bootstrap shows the share genuinely moves by
        **0.04–0.13** (PTEN miR-23-3p / miR-26-5p / miR-29-3p). **The posterior is claiming certainty it does
        not have** — it is overconfident about the *inclusion* (z collapses to 0 in every draw), which makes
        the share degenerate. `width_flag` names this; `shap_sd_boot` measures it.

    NOTE the calibration does NOT transfer from β: a Shapley share is a **RATIO** (`γ_fβ_f/Σγβ`), so a common
    inflation of all β's CANCELS. β's 1.73× overconfidence (`calibration.posterior_calibration`) therefore says
    NOTHING about the share's width — which is exactly why it had to be re-derived, not rescaled.
    `width_flag`: `ok` · `posterior_degenerate` (0-width claim the data does not support — DO NOT publish
    `shap_sd`) · `posterior_conservative` (posterior > 2× the sampling width).
    """
    from mirna_hallmark.learned import families as FAM, attribution_eb as AE, spike_slab as SS
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = FAM.family_of(pd.Index(X.columns)); Xf, wf, members = FAM.collapse_by_family(X, w, fam)
    yr, Xz, cols = AE._prep(Y, Xf, C)
    pi = np.clip(AE._evidence_pi(gene, cols), 0.0, 1.0)
    terms = multi_family_terms(gene, cols, members, assembled=(Y, X, C, w))
    multi = [t for t in terms if t["n_fam"] >= min_fam]
    if not multi:
        return pd.DataFrame()
    _, _, pip0 = SS._gibbs_posterior(Xz, yr, pi, n_iter=n_iter, burn=burn, seed=seed)
    _, _, pipC, samp = SS._gibbs_posterior(Xz, yr, pi, n_iter=n_iter, burn=burn, seed=seed,
                                           channels=terms, return_samples=True)
    rows = []
    for si, t in enumerate(multi):
        mem = t["members"]                                        # [(col, γ)]
        contrib = np.column_stack([g * samp[:, c] for (c, g) in mem])   # per-draw γ_f·β_f (Shapley of the linear π)
        tot = contrib.sum(1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(np.abs(tot) > 1e-9, contrib / tot, np.nan)
        for j, (c, g) in enumerate(mem):
            sj = share[:, j][np.isfinite(share[:, j])]
            rows.append({"gene": gene, "segment": si, "F": t["F"], "delta": t["delta"],
                         "family": cols[c], "gamma": round(g, 3), "prior_pi": round(float(pi[c]), 3),
                         "pip_mrna": round(float(pip0[c]), 3), "pip_cn": round(float(pipC[c]), 3),
                         "d_pip": round(float(pipC[c] - pip0[c]), 3),
                         "shap_share": round(float(sj.mean()), 3) if sj.size else np.nan,
                         "shap_sd": round(float(sj.std()), 3) if sj.size else np.nan})
    df = pd.DataFrame(rows)
    if not n_boot or df.empty:
        df["shap_sd_boot"] = np.nan
        df["width_flag"] = "not_computed"
        return df

    # ---- HONEST width: bootstrap the participants, refit, and take the SD of the share across resamples.
    # This is the ONLY valid correction (the β calibration factor cannot be transferred to a ratio).
    n = len(yr)
    rng = np.random.default_rng(seed)
    acc: dict = {}
    for b in range(n_boot):
        bi = rng.integers(0, n, n)
        _, _, _, sb = SS._gibbs_posterior(Xz[bi], yr[bi], pi, n_iter=max(n_iter // 2, 600),
                                          burn=max(burn // 2, 200), seed=b, channels=terms,
                                          return_samples=True)
        bm = sb.mean(0)                                        # posterior-mean β for this resample
        for si, t in enumerate(multi):
            ctr = np.array([g * bm[c] for (c, g) in t["members"]])
            tot = ctr.sum()
            if abs(tot) < 1e-9:
                continue
            for j, (c, _g) in enumerate(t["members"]):
                acc.setdefault((si, cols[c]), []).append(ctr[j] / tot)

    def _boot_sd(r):
        v = acc.get((r["segment"], r["family"]))
        return float(np.std(v)) if v and len(v) > 2 else np.nan

    df["shap_sd_boot"] = df.apply(_boot_sd, axis=1).round(3)

    def _flag(r):
        ps, bs = r["shap_sd"], r["shap_sd_boot"]
        if not np.isfinite(bs):
            return "not_computed"
        if ps < 0.01 and bs > 0.02:                            # 1.000 ± 0.000 — certainty it does not have
            return "posterior_degenerate"
        if bs > 1e-6 and ps > 2.0 * bs:
            return "posterior_conservative"
        return "ok"

    df["width_flag"] = df.apply(_flag, axis=1)
    return df
