"""Guide/passenger CN reduced-form contrast — POWERED UP with honest inference (SPEED-OPTIMIZED null).

Adds to cn_falsify.py F3:
  (1) HAIRPIN-level clustering (guide & passenger share the identical CN column ⇒ perfectly dependent
      within a (hairpin,gene) cell; original clustered only arm+gene).
  (2) ROLE-SWAP PERMUTATION NULL: flip guide<->passenger per hairpin p=0.5, recompute gap tau_G-tau_P.
      One-sided p = frac(null gap <= observed). Tests "the functional (high-abundance guide) arm represses
      more than its passenger", preserving CN/abundance/sites.
  (3) Reduced-form is the GENEROUS scale (structural Wald = reduced/first-stage divides by rho~0.2 and a
      near-zero passenger first stage ⇒ strictly larger SE). n.s. here ⇒ n.s. structural.

SPEED (axiom 3): the FE within-transform is LINEAR, so demean each arm's site-indicator vector ONCE
(DV, n x n_arms), then every draw's demeaned regressors are matvecs: Xr_G = DV @ s_guide,
Xr_P = DVsum - Xr_G. Turns 3000 full FE demeans into 3 demeans + cheap matvecs.
"""
import time, numpy as np, pandas as pd
from scipy import stats
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as INS

import os; SP = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20260718)
t0 = time.time()

R = pd.read_csv(f"{SP}/within_arm_iv_nogate.tsv", sep="\t").dropna(subset=["reduced"])
TS = pd.read_csv("data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt", sep="\t",
                 usecols=["Gene Symbol", "miRNA", "Site Type"], low_memory=False)
TS = TS.rename(columns={"Gene Symbol": "gene", "miRNA": "arm", "Site Type": "stype"})
TS = TS[TS.stype.isin([1, 2, 3])].groupby(["arm", "gene"], as_index=False).first()
KD = pd.concat([pd.read_csv(f"data/external_cache/scanmir/{f}", sep="\t", usecols=["arm", "gene"])
                for f in ["genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz"]],
               ignore_index=True).drop_duplicates()
KD["kd"] = 1
M = R.merge(TS[["arm", "gene"]].assign(ts=1), on=["arm", "gene"], how="left").merge(KD, on=["arm", "gene"], how="left")
M["ts"] = M.ts.fillna(0).astype(int)
M["kd"] = M.kd.fillna(0).astype(int)
M["site_free"] = ((M.ts == 0) & (M.kd == 0)).astype(int)
M["w"] = np.clip(M.F, 0, None)

XALL = LD._load()["X"]
AB = XALL.median(axis=1)
loc = {a: tuple(INS._arm_active_loci(a)) for a in M.arm.unique()}
byl = {}
for a, l in loc.items():
    if l:
        byl.setdefault(l, []).append(a)
prs = {l: v for l, v in byl.items() if len(v) == 2}
role, asym, hp = {}, {}, {}
for hid, (l, (a1, a2)) in enumerate(prs.items()):
    e1, e2 = float(AB.get(a1, np.nan)), float(AB.get(a2, np.nan))
    if not (np.isfinite(e1) and np.isfinite(e2)):
        continue
    g_, p_ = (a1, a2) if e1 >= e2 else (a2, a1)
    role[g_], role[p_] = "guide", "passenger"
    asym[a1] = asym[a2] = abs(e1 - e2)
    hp[a1] = hp[a2] = hid
M["role"] = M.arm.map(role); M["asym"] = M.arm.map(asym); M["hp"] = M.arm.map(hp)
P = M[M.role.notna() & ((M.ts == 1) | (M.site_free == 1))].copy()
print(f"[{time.time()-t0:.0f}s] hairpins both arms: {len(prs)} | {P.arm.nunique()} arms, {len(P):,} rows "
      f"| sites G {int((P.ts*(P.role=='guide')).sum()):,} P {int((P.ts*(P.role=='passenger')).sum()):,}")


def wtransform(Mraw, W, ai, gi, iters=20):
    """weighted 2-way (arm,gene) FE within-transform applied to every column of Mraw (n x k)."""
    X = np.asarray(Mraw, float).copy()
    na, ng = int(ai.max()) + 1, int(gi.max()) + 1
    for _ in range(iters):
        for idxs, m in ((ai, na), (gi, ng)):
            sw = np.bincount(idxs, weights=W, minlength=m); sw[sw == 0] = 1.0
            GS = np.zeros((m, X.shape[1]))
            np.add.at(GS, idxs, W[:, None] * X)
            X -= (GS / sw[:, None])[idxs]
    return X


def analytic(dd, W, cl_a, cl_g):
    """observed gap + multiway-clustered SE (cl_a,cl_g). Design = [ts_G, ts_P]."""
    ts = dd.ts.to_numpy(float)
    Xr = wtransform(np.column_stack([ts * (dd.role.to_numpy() == "guide"),
                                     ts * (dd.role.to_numpy() == "passenger")]), W,
                    dd.arm.astype("category").cat.codes.to_numpy(),
                    dd.gene.astype("category").cat.codes.to_numpy())
    yr = wtransform(dd.reduced.to_numpy(float)[:, None], W,
                    dd.arm.astype("category").cat.codes.to_numpy(),
                    dd.gene.astype("category").cat.codes.to_numpy())[:, 0]
    Xw = Xr * W[:, None]
    Xi = np.linalg.inv(Xr.T @ Xw)
    b = Xi @ (Xw.T @ yr)
    e = yr - Xr @ b
    S = Xw * e[:, None]
    def meat(ids):
        i = pd.factorize(ids)[0]; k = 2; Mt = np.zeros((k, k))
        for p in range(k):
            bp = np.bincount(i, weights=S[:, p])
            for q in range(k):
                Mt[p, q] = float(bp @ np.bincount(i, weights=S[:, q]))
        return Mt
    both = pd.Series(cl_a).astype(str) + "|" + pd.Series(cl_g).astype(str)
    V = Xi @ (meat(cl_a) + meat(cl_g) - meat(both.to_numpy())) @ Xi
    se = np.sqrt(max(V[0, 0] + V[1, 1] - 2 * V[0, 1], 1e-18))
    return b[0], b[1], (b[0] - b[1]) / se


def perm_null(dd, W, B=2000):
    """role-swap null via precomputed per-arm demeaned site vectors (linear FE ⇒ matvec per draw)."""
    arms = dd.arm.to_numpy()
    uarms = pd.unique(arms)
    ai_arm = {a: k for k, a in enumerate(uarms)}
    col = np.array([ai_arm[a] for a in arms])
    ts = dd.ts.to_numpy(float)
    # DV_raw: column per arm = that arm's site indicator (ts on its rows)
    DVraw = np.zeros((len(dd), len(uarms)))
    DVraw[np.arange(len(dd)), col] = ts
    aic = dd.arm.astype("category").cat.codes.to_numpy()
    gic = dd.gene.astype("category").cat.codes.to_numpy()
    DV = wtransform(DVraw, W, aic, gic)
    yr = wtransform(dd.reduced.to_numpy(float)[:, None], W, aic, gic)[:, 0]
    DVsum = DV.sum(1)
    Wy = W * yr
    # pairs (arm indices) + observed guide selector
    hpv = dd.hp.to_numpy()
    pair_arms = {}
    guide_obs = np.zeros(len(uarms), bool)
    for a in uarms:
        rows = arms == a
        pair_arms.setdefault(int(hpv[rows][0]), []).append(ai_arm[a])
        if dd.role.to_numpy()[rows][0] == "guide":
            guide_obs[ai_arm[a]] = True
    pairs = [tuple(v) for v in pair_arms.values() if len(v) == 2]

    def gap(sg):
        XG = DV @ sg
        XP = DVsum - XG
        WG, WP = W * XG, W * XP
        # 2x2 normal equations for [XG, XP]
        A = np.array([[XG @ WG, XG @ WP], [XP @ WG, XP @ WP]])
        r = np.array([XG @ Wy, XP @ Wy])
        b = np.linalg.solve(A, r)
        return b[0] - b[1]

    obs = gap(guide_obs.astype(float))
    null = np.empty(B)
    for bb in range(B):
        sg = guide_obs.copy()
        flip = rng.integers(0, 2, len(pairs)).astype(bool)
        for (i1, i2), f in zip(pairs, flip):
            if f:
                sg[i1] = not sg[i1]; sg[i2] = not sg[i2]
        null[bb] = gap(sg.astype(float))
    p = (np.sum(null <= obs) + 1) / (B + 1)
    return obs, p, null.mean(), null.std()


print(f"[{time.time()-t0:.0f}s] running specs...")
for lab, S, wt in [("unweighted", P, None), ("F-weighted", P, "w"),
                   ("confident guide |dAB|>=2, F-wtd", P[P.asym >= 2.0], "w")]:
    dd = S.dropna(subset=["reduced"]).reset_index(drop=True)
    if int((dd.ts * (dd.role == "guide")).sum()) < 50 or int((dd.ts * (dd.role == "passenger")).sum()) < 50:
        print(f"{lab}: skipped"); continue
    W = np.ones(len(dd)) if wt is None else dd[wt].to_numpy(float)
    tG, tP, z_arm = analytic(dd, W, dd.arm.to_numpy(), dd.gene.to_numpy())
    _, _, z_hp = analytic(dd, W, dd.hp.to_numpy(), dd.gene.to_numpy())
    obs, p_perm, nm, ns = perm_null(dd, W)
    print(f"\n### {lab}  (n={len(dd):,})   [{time.time()-t0:.0f}s]")
    print(f"  tau_G {tG:+.5f}  tau_P {tP:+.5f}  gap {tG-tP:+.5f}  (perm obs {obs:+.5f})")
    print(f"  arm+gene cluster : one-sided p={stats.norm.cdf(z_arm):.4f}   (original)")
    print(f"  HAIRPIN+gene clu : one-sided p={stats.norm.cdf(z_hp):.4f}   <- honest analytic")
    print(f"  role-swap NULL   : one-sided p={p_perm:.4f}   (null mean {nm:+.5f} sd {ns:.5f})")
print(f"\n[{time.time()-t0:.0f}s] DONE")
