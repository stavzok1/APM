"""Sharpen in TCGA: is the guide's CLEAN own effect (non-co-targeted genes, leak-free) genuine site-mediated
repression? Decisive test = does it SCALE WITH SITE STRENGTH (context++ magnitude; 8mer>7mer-m8>7mer-A1)?
Real repression must; a confound won't. Plus: profile WHY the passenger is positive.
Everything TRANS + arm/gene FE (absorbs abundance) + hairpin/gene-clustered SE.
"""
import time, numpy as np, pandas as pd
from scipy import stats
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as INS

import os; SP = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20260718); t0 = time.time(); FLOOR = np.log2(11)
R = pd.read_csv(f"{SP}/within_arm_iv_nogate.tsv", sep="\t").dropna(subset=["reduced"])
TS = pd.read_csv("data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt", sep="\t",
                 usecols=["Gene Symbol", "miRNA", "Site Type", "context++ score"], low_memory=False).rename(
                 columns={"Gene Symbol": "gene", "miRNA": "arm", "Site Type": "stype", "context++ score": "ctx"})
TS = TS[TS.stype.isin([1, 2, 3])].sort_values("stype", ascending=False).groupby(["arm", "gene"], as_index=False).first()
KD = pd.concat([pd.read_csv(f"data/external_cache/scanmir/{f}", sep="\t", usecols=["arm", "gene"])
                for f in ["genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz"]],
               ignore_index=True).drop_duplicates(); KD["kd"] = 1
M = R.merge(TS, on=["arm", "gene"], how="left").merge(KD, on=["arm", "gene"], how="left")
M["ts"] = M.stype.notna().astype(int); M["kd"] = M.kd.fillna(0).astype(int)
M["site_free"] = ((M.ts == 0) & (M.kd == 0)).astype(int)
XALL = LD._load()["X"]; AB = XALL.apply(pd.to_numeric, errors="coerce").median(axis=1)
loc = {a: tuple(INS._arm_active_loci(a)) for a in M.arm.unique()}
byl = {}
for a, l in loc.items():
    if l: byl.setdefault(l, []).append(a)
prs = {l: v for l, v in byl.items() if len(v) == 2}
hp_of, hp_loci, role, hp_arms = {}, {}, {}, {}
for hid, (l, (a1, a2)) in enumerate(prs.items()):
    if not (np.isfinite(AB.get(a1, np.nan)) and np.isfinite(AB.get(a2, np.nan))): continue
    hp_of[a1] = hp_of[a2] = hid; hp_loci[hid] = list(l); hp_arms[hid] = (a1, a2)
    g_, p_ = (a1, a2) if AB[a1] >= AB[a2] else (a2, a1); role[g_], role[p_] = "guide", "passenger"
M["hp"] = M.arm.map(hp_of); M["role"] = M.arm.map(role)
detarms = set(AB[AB > FLOOR].index)
CL = INS._genomic_clusters(sorted(detarms | set(M.arm.dropna().unique())))
clm = {}
for a, c in CL.items(): clm.setdefault(c, []).append(a)
hp_co = {h: [a for a in clm.get(CL.get(a1), []) if a not in (a1, a2) and a in detarms] for h, (a1, a2) in hp_arms.items()}
ts_by_arm = TS.groupby("arm").gene.apply(set).to_dict()
cotgt = {h: set().union(*[ts_by_arm.get(a, set()) for a in co]) if co else set() for h, co in hp_co.items()}
LC = pd.read_csv("mirna_hallmark/output/learned/locus_context.tsv", sep="\t").set_index("locus_id")["chrom"]
hp_chrom = {h: next((LC.get(c) for c in li if c in LC.index), None) for h, li in hp_loci.items()}
gc = pd.read_csv("data/gencode.v49.annotation.gtf.csv", usecols=["seqname", "feature", "gene_name"])
gc = gc[gc.feature == "gene"].drop_duplicates("gene_name").set_index("gene_name")["seqname"]
M["hp_chr"] = M.hp.map(hp_chrom); M["g_chr"] = M.gene.map(gc)
TR = M[(M.hp.notna()) & ((M.ts == 1) | (M.site_free == 1)) & M.hp_chr.notna() & M.g_chr.notna()
       & (M.hp_chr != M.g_chr)].copy()
TR["co_tgt"] = [int(g in cotgt.get(int(h), set())) for h, g in zip(TR.hp, TR.gene)]


def demean(V, W, ai, gi, it=20):
    X = np.asarray(V, float).copy(); na, ng = int(ai.max()) + 1, int(gi.max()) + 1
    for _ in range(it):
        for idxs, m in ((ai, na), (gi, ng)):
            sw = np.bincount(idxs, weights=W, minlength=m); sw[sw == 0] = 1.0
            if X.ndim == 1: X -= (np.bincount(idxs, weights=W * X, minlength=m) / sw)[idxs]
            else:
                GS = np.zeros((m, X.shape[1])); np.add.at(GS, idxs, W[:, None] * X); X -= (GS / sw[:, None])[idxs]
    return X


def fe(dd, Xcols):
    dd = dd.reset_index(drop=True); W = np.ones(len(dd))
    ai = dd.arm.astype("category").cat.codes.to_numpy(); gi = dd.gene.astype("category").cat.codes.to_numpy()
    yr = demean(dd.reduced.to_numpy(float), W, ai, gi); Xr = demean(dd[Xcols].to_numpy(float), W, ai, gi)
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


# ===== clean guide-own subset: non-co-targeted guide sites + site-free =====
GO = TR[((TR.role == "guide") & (TR.ts == 1) & (TR.co_tgt == 0)) | (TR.site_free == 1)].copy()
print(f"[{time.time()-t0:.0f}s] clean guide-own (non-co-targeted) sites: {int(((GO.role=='guide')&(GO.ts==1)).sum()):,}")

# (1) SITE-TYPE ladder on the clean guide-own effect
for t in (3, 2, 1): GO[f"s{t}"] = ((GO.ts == 1) & (GO.stype == t)).astype(float)
b, se = fe(GO, ["s3", "s2", "s1"])
print("  SITE-TYPE ladder (baseline site-free):")
for i, nm in [(0, "8mer  "), (1, "7mer-m8"), (2, "7mer-A1")]:
    print(f"    tau[{nm}] = {b[i]:+.5f}  SE {se[i]:.5f}  p={2*stats.norm.sf(abs(b[i]/se[i])):.3f}")
print(f"    monotone 8>7m8>7A1? {b[0]<b[1]<b[2]}   8mer vs 7mer-A1 diff {b[0]-b[2]:+.5f}")

# (2) CONTINUOUS context++ strength: among clean guide-own SITE rows, does stronger site => more neg reduced?
GS = GO[(GO.role == "guide") & (GO.ts == 1) & GO.ctx.notna()].copy()
GS["ctxz"] = (GS.ctx - GS.ctx.mean()) / (GS.ctx.std() + 1e-9)   # ctx more NEGATIVE = stronger repression
bc, sec = fe(GS, ["ctxz"])
print(f"  CONTINUOUS context++ slope (among clean guide sites, n={len(GS):,}): "
      f"beta={bc[0]:+.5f} SE {sec[0]:.5f} p={2*stats.norm.sf(abs(bc[0]/sec[0])):.3f}"
      f"  (>0 => stronger/more-neg site -> more-neg reduced = repression signature)")

# ===== why is the passenger positive? profile passenger sites vs guide sites =====
print(f"\n[{time.time()-t0:.0f}s] === passenger positive: profiling ===")
gsite = TR[(TR.role == "guide") & (TR.ts == 1)]; psite = TR[(TR.role == "passenger") & (TR.ts == 1)]
for lab, s in [("GUIDE sites", gsite), ("PASSENGER sites", psite)]:
    st = s.stype.value_counts(normalize=True)
    print(f"  {lab:16s} n={len(s):,} | 8mer {st.get(3,0):.0%} 7m8 {st.get(2,0):.0%} 7A1 {st.get(1,0):.0%} | "
          f"mean ctx {s.ctx.mean():+.3f} | mean arm abund {AB.reindex(s.arm).mean():.2f} | mean reduced {s.reduced.mean():+.4f}")
# passenger effect vs site strength (does the positive concentrate in WEAK passenger sites = noise?)
PO = TR[((TR.role == "passenger") & (TR.ts == 1)) | (TR.site_free == 1)].copy()
for t in (3, 2, 1): PO[f"s{t}"] = ((PO.ts == 1) & (PO.stype == t)).astype(float)
bp, sep = fe(PO, ["s3", "s2", "s1"])
print("  passenger tau by site type (baseline site-free):")
for i, nm in [(0, "8mer  "), (1, "7mer-m8"), (2, "7mer-A1")]:
    print(f"    tau[{nm}] = {bp[i]:+.5f}  SE {sep[i]:.5f}  p={2*stats.norm.sf(abs(bp[i]/sep[i])):.3f}")
print(f"[{time.time()-t0:.0f}s] DONE")
