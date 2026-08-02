"""⭐ MH-129e — THE FALSIFICATION PROGRAM FOR THE CN INSTRUMENT. (User: "the genuinely site-free control was
never built — then it needs to. go on on falsification.")

THE ESTIMATOR (established):  reduced(m,g) = alpha_m + gamma_g + tau*site(m,g) + eps
  alpha_m absorbs ALL locus pleiotropy (a FE, so no observational slope is ever touched — unlike exclusion()).
  gamma_g absorbs the gene's CN sensitivity.  C_g = [CPE, target_cn(g), mal_prolif, PLOIDY], rebuilt per gene.
  tau is identified by SITE PRESENCE = SEQUENCE = EXOGENOUS TO COPY NUMBER.
  CURRENT: tau = -0.0039 (p=0.029, 2-way clustered), placebo at 0, Wald -0.019 [-0.034, -0.003].

THREE FIXES / FALSIFICATIONS:

  F1 ⭐ THE GENUINELY SITE-FREE CONTROL. My control was "not curated" (and then "no TargetScan"). But scanMiR
     duplexes cover most of the panel, so the control class was full of REAL-but-unannotated binders => tau was
     DILUTED. The genome-wide scanMiR scan EXISTS (data/external_cache/scanmir/genomewide_kd*.tsv.gz — the same
     source MH-127's build_fakes.py uses). Build: SITE-FREE = no TargetScan context++ site AND no scanMiR
     duplex. PREDICTION: |tau| INCREASES.

  F2 ⭐⭐ THE SITE-TYPE EFFICACY LADDER — THE REAL FALSIFICATION.
     TargetScan Site Type is the canonical efficacy ordering (Bartel): 8mer(3) > 7mer-m8(2) > 7mer-A1(1).
     Median context++ in this very file: -0.263 / -0.177 / -0.134 — the ladder is right there.
     **COPY NUMBER CANNOT SEE A SITE TYPE.** If the CN-mediated effect is genuine site-mediated repression it
     MUST follow the ladder. If it does not, the effect is NOT repression and the causal reading dies.
     (The CONTINUOUS context++ gradient already FAILED among site-bearing pairs, p=0.55. The site-type ladder
     is coarser and far more robust — it is the fair version of that test.)
     PREDICTION: tau_8mer < tau_7mer-m8 < tau_7mer-A1 < 0, monotone.

  F3 GUIDE vs PASSENGER, POWERED. The F>10 gate keeps only 231/515 arms and starves the passenger side
     (1,011 predicted sites vs the guide's 7,025) => the contrast was underpowered (p=0.25), NOT refuted.
     FIX: drop the hard gate; keep ALL arms with a usable instrument and WEIGHT by first-stage precision.
     PREDICTION: tau_GUIDE < 0, tau_PASSENGER ~ 0, and now the DIFFERENCE resolves.

⚠ This script REBUILDS the reduced-form panel WITHOUT the F>10 gate (F kept as a column for weighting).
"""
import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as INS

import os; OUT = os.path.dirname(os.path.abspath(__file__))
LD.assemble_gene("PTEN", w_prior_source="ledger")
d = LD._load()
XALL, YALL = d["X"], d["Y"]
LCN = INS.locus_cn().apply(pd.to_numeric, errors="coerce")
PLOIDY = LCN.mean(axis=1)

E = D.load_hallmark_edges()
HE = {(r.miRNA, r.gene) for r in E[E.high_evidence].itertuples()}
ARMS = sorted({a for a, _ in HE if a in XALL.index})
GENES = sorted({g for _, g in HE if g in YALL.index})
_tcn = LD._target_cn(GENES)    # dict {gene: Series(participant)}
TCN = pd.DataFrame({g: s for g, s in _tcn.items() if s is not None}).T  # gene rows, participant cols

_, _, C0, _ = LD.assemble_gene(GENES[0], w_prior_source="ledger")
IDX = C0.index
n = len(IDX)
Zc, Za = [], []
for a in ARMS:
    li = [c for c in INS._arm_active_loci(a) if c in LCN.columns]
    if not li:
        continue
    z = LCN[li].mean(axis=1).reindex(IDX)
    if z.isna().mean() > 0.2 or float(np.nanstd(z)) < 1e-6:
        continue
    Zc.append(z.fillna(z.median()).to_numpy(float)); Za.append(a)
Z = np.column_stack(Zc)
Xm = XALL.reindex(Za).T.reindex(IDX).apply(pd.to_numeric, errors="coerce")
Xm = Xm.fillna(Xm.median()).to_numpy(float)
pl = PLOIDY.reindex(IDX).fillna(PLOIDY.median()).to_numpy(float)
cpe = pd.to_numeric(C0["CPE"], errors="coerce").fillna(0).to_numpy(float)
prl = pd.to_numeric(C0["mal_prolif"], errors="coerce").fillna(0).to_numpy(float)

_res = lambda V, A: V - A @ np.linalg.lstsq(A, V, rcond=None)[0]
A0 = np.column_stack([np.ones(n), cpe, prl, pl])
Zr0, Xr0 = _res(Z, A0), _res(Xm, A0)
zz = (Zr0 ** 2).sum(0)
a_fs = (Zr0 * Xr0).sum(0) / np.where(zz > 0, zz, 1)
rss = ((Xr0 - Zr0 * a_fs) ** 2).sum(0)
Fs = ((Xr0 ** 2).sum(0) - rss) / (rss / (n - 5))
print(f"arms with locus CN: {len(Za)} | F>10: {(Fs>10).sum()} | ⭐ NO GATE — all {len(Za)} kept, F is a WEIGHT")

rows = []
for j, g in enumerate(GENES):
    y = pd.to_numeric(YALL.loc[g].reindex(IDX), errors="coerce")
    if y.isna().mean() > 0.2 or float(np.nanstd(y)) < 1e-6:
        continue
    y = y.fillna(y.median()).to_numpy(float)
    tc = pd.to_numeric(TCN.loc[g], errors="coerce").reindex(IDX) if g in TCN.index else None
    tc = tc.fillna(tc.median()).to_numpy(float) if tc is not None and tc.notna().any() else np.zeros(n)
    Ag = np.column_stack([np.ones(n), cpe, prl, pl, tc])
    yr = _res(y[:, None], Ag)[:, 0]
    Zr = _res(Z, Ag)
    ry = stats.rankdata(yr); ry = (ry - ry.mean()) / (ry.std() + 1e-12)
    rZ = np.apply_along_axis(stats.rankdata, 0, Zr)
    rZ = (rZ - rZ.mean(0)) / (rZ.std(0) + 1e-12)
    red = (rZ * ry[:, None]).mean(0)
    for i, a in enumerate(Za):
        rows.append((a, g, float(red[i]), int((a, g) in HE), float(Fs[i])))
    if (j + 1) % 400 == 0:
        print(f"   {j+1}/{len(GENES)}")
R = pd.DataFrame(rows, columns=["arm", "gene", "reduced", "curated", "F"]).dropna(subset=["reduced"])
R.to_csv(f"{OUT}/within_arm_iv_nogate.tsv", sep="\t", index=False)
print(f"{len(R):,} arm-gene reduced forms | {R.arm.nunique()} arms x {R.gene.nunique()} genes")

# ---------------- SITE ANNOTATIONS ----------------
TS = pd.read_csv("data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt", sep="\t",
                 usecols=["Gene Symbol", "miRNA", "Site Type", "context++ score"], low_memory=False)
TS = TS.rename(columns={"Gene Symbol": "gene", "miRNA": "arm", "Site Type": "stype",
                        "context++ score": "ctx"})
TS = TS[TS.stype.isin([1, 2, 3])]
TS = TS.sort_values("stype", ascending=False).groupby(["arm", "gene"], as_index=False).first()  # BEST site type
print(f"TargetScan: {len(TS):,} pairs | 8mer {int((TS.stype==3).sum()):,} "
      f"7mer-m8 {int((TS.stype==2).sum()):,} 7mer-A1 {int((TS.stype==1).sum()):,}")

KD = pd.concat([pd.read_csv(f"data/external_cache/scanmir/{f}", sep="\t", usecols=["arm", "gene"])
                for f in ["genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz"]],
               ignore_index=True).drop_duplicates()
KD["kd"] = 1
print(f"⭐ genome-wide scanMiR duplexes: {len(KD):,} arm-gene pairs")

M = R.merge(TS[["arm", "gene", "stype", "ctx"]], on=["arm", "gene"], how="left") \
     .merge(KD, on=["arm", "gene"], how="left")
M["ts"] = M.stype.notna().astype(int)
M["kd"] = M.kd.fillna(0).astype(int)
M["site_free"] = ((M.ts == 0) & (M.kd == 0)).astype(int)
print(f"\nOF {len(M):,} PANEL PAIRS:")
print(f"  TargetScan site        : {int(M.ts.sum()):,} ({M.ts.mean():.1%})")
print(f"  scanMiR duplex         : {int(M.kd.sum()):,} ({M.kd.mean():.1%})   <- the DILUTION")
print(f"  ⭐ GENUINELY SITE-FREE  : {int(M.site_free.sum()):,} ({M.site_free.mean():.1%})")
print(f"  curated (HE)           : {int(M.curated.sum()):,} ({M.curated.mean():.2%})")


def fe(df, xcols, w=None):
    dd = df.dropna(subset=["reduced"] + xcols).reset_index(drop=True)
    Y = dd["reduced"].to_numpy(float).copy()
    X = dd[xcols].to_numpy(float).copy()
    W = np.ones(len(dd)) if w is None else dd[w].to_numpy(float)
    ai = dd.arm.astype("category").cat.codes.to_numpy()
    gi = dd.gene.astype("category").cat.codes.to_numpy()
    for _ in range(15):
        for idxs in (ai, gi):
            sw = np.bincount(idxs, weights=W); sw[sw == 0] = 1
            Y -= (np.bincount(idxs, weights=W * Y) / sw)[idxs]
            for j in range(X.shape[1]):
                X[:, j] -= (np.bincount(idxs, weights=W * X[:, j]) / sw)[idxs]
    Xw = X * W[:, None]
    Xi = np.linalg.inv(X.T @ Xw)
    b = Xi @ (Xw.T @ Y)
    e = Y - X @ b
    S = Xw * e[:, None]
    def meat(i):
        k = X.shape[1]; Mt = np.zeros((k, k))
        for p in range(k):
            for q in range(k):
                Mt[p, q] = float(np.bincount(i, weights=S[:, p]) @ np.bincount(i, weights=S[:, q]))
        return Mt
    both = pd.factorize(pd.Series(ai).astype(str) + "|" + pd.Series(gi).astype(str))[0]
    V = Xi @ (meat(ai) + meat(gi) - meat(both)) @ Xi
    return b, np.sqrt(np.maximum(np.diag(V), 0)), V, len(dd)


print("\n" + "=" * 96)
print("=== F1 — THE GENUINELY SITE-FREE CONTROL (no TargetScan site AND no scanMiR duplex) ===")
M["w"] = np.clip(M.F, 0, None)                              # first-stage precision weight (no hard gate)
for lab, sub, wt in [
    ("curated vs NOT-CURATED (old, diluted)", M, None),
    ("curated vs no-TargetScan", M[(M.curated == 1) | (M.ts == 0)], None),
    ("⭐ curated vs GENUINELY SITE-FREE", M[(M.curated == 1) | (M.site_free == 1)], None),
    ("⭐ same, F-WEIGHTED (no gate)", M[(M.curated == 1) | (M.site_free == 1)], "w"),
    ("⭐ same, F>10 only (the old gate)", M[((M.curated == 1) | (M.site_free == 1)) & (M.F > 10)], None),
]:
    b, se, _, nn = fe(sub, ["curated"], wt)
    z = b[0] / se[0]
    print(f"  {lab:40s} tau={b[0]:+.5f}  SE={se[0]:.5f}  z={z:+.2f}  "
          f"p={2*stats.norm.sf(abs(z)):.4f}  n={nn:,}")

print("\n" + "=" * 96)
print("=== F2 ⭐⭐ THE SITE-TYPE EFFICACY LADDER — COPY NUMBER CANNOT SEE A SITE TYPE ===")
print("    canonical (Bartel): 8mer > 7mer-m8 > 7mer-A1.  tau must follow it, or the effect is NOT repression.")
L = M[(M.site_free == 1) | (M.ts == 1)].copy()             # site-free controls + TargetScan-site pairs
for t, nm in [(3, "8mer"), (2, "7mer-m8"), (1, "7mer-A1")]:
    L[f"s{t}"] = ((L.ts == 1) & (L.stype == t)).astype(int)
b, se, V, nn = fe(L, ["s3", "s2", "s1"])
print(f"    n={nn:,}  (baseline = GENUINELY SITE-FREE)")
for i, (t, nm) in enumerate([(3, "8mer  "), (2, "7mer-m8"), (1, "7mer-A1")]):
    z = b[i] / se[i]
    print(f"      tau[{nm}] = {b[i]:+.5f}   SE {se[i]:.5f}   z {z:+.2f}   p={2*stats.norm.sf(abs(z)):.4f}"
          f"   (n_sites {int(L[f's{t}'].sum()):,})")
dv = V[0, 0] + V[2, 2] - 2 * V[0, 2]
zl = (b[0] - b[2]) / np.sqrt(max(dv, 1e-18))
mono = (b[0] < b[1] < b[2])
print(f"    LADDER: 8mer vs 7mer-A1  diff {b[0]-b[2]:+.5f}  z={zl:+.2f}  one-sided p={stats.norm.cdf(zl):.4f}")
print(f"    MONOTONE (8mer < 7mer-m8 < 7mer-A1)? {mono}   "
      f"{'** LADDER CONFIRMED — this is site-mediated repression' if mono and stats.norm.cdf(zl)<0.05 else '   NOT CONFIRMED'}")

print("\n" + "=" * 96)
print("=== F3 — GUIDE vs PASSENGER, POWERED (no F gate; TargetScan-predicted sites) ===")
AB = XALL.median(axis=1)
loc = {a: tuple(INS._arm_active_loci(a)) for a in M.arm.unique()}
byl = {}
for a, l in loc.items():
    if l:
        byl.setdefault(l, []).append(a)
prs = {l: v for l, v in byl.items() if len(v) == 2}
role, asym = {}, {}
for l, (a1, a2) in prs.items():
    e1, e2 = float(AB.get(a1, np.nan)), float(AB.get(a2, np.nan))
    if not (np.isfinite(e1) and np.isfinite(e2)):
        continue
    g_, p_ = (a1, a2) if e1 >= e2 else (a2, a1)
    role[g_], role[p_] = "guide", "passenger"
    asym[a1] = asym[a2] = abs(e1 - e2)
M["role"] = M.arm.map(role); M["asym"] = M.arm.map(asym)
P = M[M.role.notna() & ((M.ts == 1) | (M.site_free == 1))].copy()
P["ts_G"] = P.ts * (P.role == "guide"); P["ts_P"] = P.ts * (P.role == "passenger")
print(f"  hairpins with BOTH arms: {len(prs)}  |  {P.arm.nunique()} arms, {len(P):,} rows")
print(f"  predicted sites — GUIDE {int(P.ts_G.sum()):,}  PASSENGER {int(P.ts_P.sum()):,}  "
      f"(was 7,025 / 1,011 under the F>10 gate)")
for lab, S, wt in [("unweighted", P, None), ("F-WEIGHTED", P, "w"),
                   ("confident guide |Δ|>=2, F-wtd", P[P.asym >= 2.0], "w")]:
    if S.ts_G.sum() < 50 or S.ts_P.sum() < 50:
        continue
    b, se, V, nn = fe(S, ["ts_G", "ts_P"], wt)
    dv = V[0, 0] + V[1, 1] - 2 * V[0, 1]
    zd = (b[0] - b[1]) / np.sqrt(max(dv, 1e-18))
    print(f"\n  {lab:30s} n={nn:,}")
    print(f"    tau_GUIDE     = {b[0]:+.5f}  SE {se[0]:.5f}  z {b[0]/se[0]:+.2f}  p={2*stats.norm.sf(abs(b[0]/se[0])):.4f}")
    print(f"    tau_PASSENGER = {b[1]:+.5f}  SE {se[1]:.5f}  z {b[1]/se[1]:+.2f}  p={2*stats.norm.sf(abs(b[1]/se[1])):.4f}")
    print(f"    H1 (G < P): diff {b[0]-b[1]:+.5f}  z={zd:+.2f}  one-sided p={stats.norm.cdf(zd):.4f}"
          f"  {'** CONFIRMED' if stats.norm.cdf(zd)<0.05 else '   not resolved'}")
