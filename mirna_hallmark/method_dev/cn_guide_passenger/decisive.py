"""Disentangle CO-CLUSTER LEAK from ABUNDANCE (user caution). solo-vs-clustered is confounded with abundance
(clustered = high-abundance oncomiR polycistrons). CLEAN test: hold the arm set fixed (clustered hairpins),
split the GUIDE's own target genes into co-targeted-by-a-co-cluster-member vs NOT, and estimate the guide
site effect for each. Arm FE absorbs abundance level, so this isolates the co-targeting channel from abundance.
Also: profile abundance+variance of solo vs clustered guides; and the guide's OWN effect (tau_G, non-cotgt,
site vs site-free) under hairpin+gene clustering.
"""
import time, numpy as np, pandas as pd
from scipy import stats
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as INS

import os; SP = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20260718); t0 = time.time(); FLOOR = np.log2(11)
R = pd.read_csv(f"{SP}/within_arm_iv_nogate.tsv", sep="\t").dropna(subset=["reduced"])
TSraw = pd.read_csv("data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt", sep="\t",
                    usecols=["Gene Symbol", "miRNA", "Site Type"], low_memory=False).rename(
                    columns={"Gene Symbol": "gene", "miRNA": "arm", "Site Type": "stype"})
TSraw = TSraw[TSraw.stype.isin([1, 2, 3])].groupby(["arm", "gene"], as_index=False).first()
KD = pd.concat([pd.read_csv(f"data/external_cache/scanmir/{f}", sep="\t", usecols=["arm", "gene"])
                for f in ["genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz"]],
               ignore_index=True).drop_duplicates(); KD["kd"] = 1
M = R.merge(TSraw[["arm", "gene"]].assign(ts=1), on=["arm", "gene"], how="left").merge(KD, on=["arm", "gene"], how="left")
M["ts"] = M.ts.fillna(0).astype(int); M["kd"] = M.kd.fillna(0).astype(int)
M["site_free"] = ((M.ts == 0) & (M.kd == 0)).astype(int)
XALL = LD._load()["X"]; Xn = XALL.apply(pd.to_numeric, errors="coerce"); AB = Xn.median(axis=1); VAR = Xn.var(axis=1)
loc = {a: tuple(INS._arm_active_loci(a)) for a in M.arm.unique()}
byl = {}
for a, l in loc.items():
    if l: byl.setdefault(l, []).append(a)
prs = {l: v for l, v in byl.items() if len(v) == 2}
hp_of, hp_loci, role, hp_arms, guide_ab, guide_dab = {}, {}, {}, {}, {}, {}
for hid, (l, (a1, a2)) in enumerate(prs.items()):
    if not (np.isfinite(AB.get(a1, np.nan)) and np.isfinite(AB.get(a2, np.nan))): continue
    hp_of[a1] = hp_of[a2] = hid; hp_loci[hid] = list(l); hp_arms[hid] = (a1, a2)
    g_, p_ = (a1, a2) if AB[a1] >= AB[a2] else (a2, a1); role[g_], role[p_] = "guide", "passenger"
    guide_ab[hid] = AB[g_]; guide_dab[hid] = abs(AB[a1] - AB[a2])
M["hp"] = M.arm.map(hp_of); M["role"] = M.arm.map(role)
detarms = set(AB[AB > FLOOR].index)
CL = INS._genomic_clusters(sorted(detarms | set(M.arm.dropna().unique())))
clm = {}
for a, c in CL.items(): clm.setdefault(c, []).append(a)
hp_co = {h: [a for a in clm.get(CL.get(a1), []) if a not in (a1, a2) and a in detarms] for h, (a1, a2) in hp_arms.items()}
ts_by_arm = TSraw.groupby("arm").gene.apply(set).to_dict()
cotgt = {h: set().union(*[ts_by_arm.get(a, set()) for a in co]) if co else set() for h, co in hp_co.items()}
M["solo"] = M.hp.map({h: (len(hp_co[h]) == 0) for h in hp_arms})
LC = pd.read_csv("mirna_hallmark/output/learned/locus_context.tsv", sep="\t").set_index("locus_id")["chrom"]
hp_chrom = {h: next((LC.get(c) for c in li if c in LC.index), None) for h, li in hp_loci.items()}
gc = pd.read_csv("data/gencode.v49.annotation.gtf.csv", usecols=["seqname", "feature", "gene_name"])
gc = gc[gc.feature == "gene"].drop_duplicates("gene_name").set_index("gene_name")["seqname"]
M["hp_chr"] = M.hp.map(hp_chrom); M["g_chr"] = M.gene.map(gc)
TR = M[(M.hp.notna()) & ((M.ts == 1) | (M.site_free == 1)) & M.hp_chr.notna() & M.g_chr.notna()
       & (M.hp_chr != M.g_chr)].copy()
TR["co_tgt"] = [int(g in cotgt.get(int(h), set())) for h, g in zip(TR.hp, TR.gene)]

# --- abundance/variance profile: solo vs clustered (is solo just lower-abundance?) ---
soloh = [h for h in hp_arms if not hp_co[h]]; clh = [h for h in hp_arms if hp_co[h]]
print(f"[{time.time()-t0:.0f}s] guide ABUNDANCE (median log2RPM): solo {np.median([guide_ab[h] for h in soloh]):.2f} "
      f"vs clustered {np.median([guide_ab[h] for h in clh]):.2f}  | guide |dAB|: solo "
      f"{np.median([guide_dab[h] for h in soloh]):.2f} vs clustered {np.median([guide_dab[h] for h in clh]):.2f}  "
      f"| guide VARIANCE: solo {np.median([VAR.get(hp_arms[h][0] if AB[hp_arms[h][0]]>=AB[hp_arms[h][1]] else hp_arms[h][1]) for h in soloh]):.2f} "
      f"vs clustered {np.median([VAR.get(hp_arms[h][0] if AB[hp_arms[h][0]]>=AB[hp_arms[h][1]] else hp_arms[h][1]) for h in clh]):.2f}")


def demean(V, W, ai, gi, it=20):
    X = np.asarray(V, float).copy(); na, ng = int(ai.max()) + 1, int(gi.max()) + 1
    for _ in range(it):
        for idxs, m in ((ai, na), (gi, ng)):
            sw = np.bincount(idxs, weights=W, minlength=m); sw[sw == 0] = 1.0
            if X.ndim == 1: X -= (np.bincount(idxs, weights=W * X, minlength=m) / sw)[idxs]
            else:
                GS = np.zeros((m, X.shape[1])); np.add.at(GS, idxs, W[:, None] * X); X -= (GS / sw[:, None])[idxs]
    return X


def clustered_se(dd, Xcols):
    """FE(arm,gene) WLS of reduced ~ Xcols; hairpin+gene clustered SE. Returns b, se per col."""
    dd = dd.reset_index(drop=True); W = np.ones(len(dd))
    ai = dd.arm.astype("category").cat.codes.to_numpy(); gi = dd.gene.astype("category").cat.codes.to_numpy()
    yr = demean(dd.reduced.to_numpy(float), W, ai, gi)
    Xr = demean(dd[Xcols].to_numpy(float), W, ai, gi)
    Xi = np.linalg.inv(Xr.T @ Xr); b = Xi @ (Xr.T @ yr); e = yr - Xr @ b; S = Xr * e[:, None]; k = len(Xcols)
    def meat(ids):
        i = pd.factorize(ids)[0]; Mt = np.zeros((k, k))
        for p in range(k):
            bp = np.bincount(i, weights=S[:, p])
            for q in range(k): Mt[p, q] = float(bp @ np.bincount(i, weights=S[:, q]))
        return Mt
    both = (pd.Series(dd.hp.to_numpy()).astype(str) + "|" + pd.Series(gi).astype(str)).to_numpy()
    V = Xi @ (meat(dd.hp.to_numpy()) + meat(gi) - meat(both)) @ Xi
    return b, np.sqrt(np.maximum(np.diag(V), 0))


# --- DECISIVE: within CLUSTERED hairpins, guide site split by co-targeting; passenger as one term ---
CLU = TR[TR.solo == False].copy()
CLU["gG0"] = ((CLU.role == "guide") & (CLU.ts == 1) & (CLU.co_tgt == 0)).astype(float)  # guide, NOT co-targeted
CLU["gG1"] = ((CLU.role == "guide") & (CLU.ts == 1) & (CLU.co_tgt == 1)).astype(float)  # guide, co-targeted
CLU["gP"] = ((CLU.role == "passenger") & (CLU.ts == 1)).astype(float)
b, se = clustered_se(CLU, ["gG0", "gG1", "gP"])
print(f"\n[{time.time()-t0:.0f}s] === CLUSTERED hairpins (abundance-matched): guide effect split by co-targeting ===")
for nm, i in [("guide, NON-co-targeted (its OWN targets)", 0), ("guide, CO-targeted (leak-prone)", 1),
              ("passenger", 2)]:
    print(f"    tau[{nm:38s}] = {b[i]:+.5f}  SE {se[i]:.5f}  p={2*stats.norm.sf(abs(b[i]/se[i])):.3f}  "
          f"(n_site={int(CLU[['gG0','gG1','gP'][i]].sum()):,})")
print(f"    => if guide-OWN (non-cotgt) ~ 0 and guide-CO << 0, the signal is co-cluster leak, NOT guide repression.")

# same split on SOLO for reference (guide only has non-cotgt by definition)
SOLO = TR[TR.solo == True].copy()
SOLO["gG"] = ((SOLO.role == "guide") & (SOLO.ts == 1)).astype(float); SOLO["gP"] = ((SOLO.role == "passenger") & (SOLO.ts == 1)).astype(float)
bs, ses = clustered_se(SOLO, ["gG", "gP"])
print(f"    [ref] SOLO hairpins: guide tau={bs[0]:+.5f} (SE {ses[0]:.5f}, p={2*stats.norm.sf(abs(bs[0]/ses[0])):.3f}), "
      f"passenger {bs[1]:+.5f}")
print(f"[{time.time()-t0:.0f}s] DONE")
