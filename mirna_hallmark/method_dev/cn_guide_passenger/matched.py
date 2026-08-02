"""User's refinement: make the guide/passenger comparison FAIR by matching gene-sets on position + leak, and
comparing arms at a comparable HIGH position within the gene's regulatory profile.
 - arm_rank = the arm's abundance rank among the gene's detectable site-bearing regulators (1=top).
 - HIGH position = rank<=2; leak-free = not co-targeted by a co-cluster member.
 - Compare guide-high-leakfree vs passenger-high-leakfree (matched position + leak), + site-strength test.
 - Also COMPARABLE-ABUNDANCE hairpins (|dAB| small) where both arms are functional.
TRANS, arm/gene FE, hairpin/gene-clustered SE + role-swap null.
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
XALL = LD._load()["X"]; AB = XALL.apply(pd.to_numeric, errors="coerce").median(axis=1); RPM = np.power(2.0, AB) - 1
loc = {a: tuple(INS._arm_active_loci(a)) for a in M.arm.unique()}
byl = {}
for a, l in loc.items():
    if l: byl.setdefault(l, []).append(a)
prs = {l: v for l, v in byl.items() if len(v) == 2}
hp_of, hp_loci, role, hp_arms, dab = {}, {}, {}, {}, {}
for hid, (l, (a1, a2)) in enumerate(prs.items()):
    if not (np.isfinite(AB.get(a1, np.nan)) and np.isfinite(AB.get(a2, np.nan))): continue
    hp_of[a1] = hp_of[a2] = hid; hp_loci[hid] = list(l); hp_arms[hid] = (a1, a2)
    g_, p_ = (a1, a2) if AB[a1] >= AB[a2] else (a2, a1); role[g_], role[p_] = "guide", "passenger"; dab[hid] = abs(AB[a1]-AB[a2])
M["hp"] = M.arm.map(hp_of); M["role"] = M.arm.map(role); M["dab"] = M.hp.map(dab)
detarms = set(AB[AB > FLOOR].index)
CL = INS._genomic_clusters(sorted(detarms | set(M.arm.dropna().unique())))
clm = {}
for a, c in CL.items(): clm.setdefault(c, []).append(a)
hp_co = {h: [a for a in clm.get(CL.get(a1), []) if a not in (a1, a2) and a in detarms] for h, (a1, a2) in hp_arms.items()}
ts_by_arm = TS.groupby("arm").gene.apply(set).to_dict()
cotgt = {h: set().union(*[ts_by_arm.get(a, set()) for a in co]) if co else set() for h, co in hp_co.items()}
TSd = TS[TS.arm.isin(detarms)].copy(); TSd["rpm"] = TSd.arm.map(RPM).fillna(0.0)
TSd["rk"] = TSd.groupby("gene").rpm.rank(ascending=False, method="min")
rank_map = {(a, g): rk for a, g, rk in zip(TSd.arm, TSd.gene, TSd.rk)}
M["arm_rank"] = [rank_map.get((a, g), np.nan) for a, g in zip(M.arm, M.gene)]
LC = pd.read_csv("mirna_hallmark/output/learned/locus_context.tsv", sep="\t").set_index("locus_id")["chrom"]
hp_chrom = {h: next((LC.get(c) for c in li if c in LC.index), None) for h, li in hp_loci.items()}
gc = pd.read_csv("data/gencode.v49.annotation.gtf.csv", usecols=["seqname", "feature", "gene_name"])
gc = gc[gc.feature == "gene"].drop_duplicates("gene_name").set_index("gene_name")["seqname"]
M["hp_chr"] = M.hp.map(hp_chrom); M["g_chr"] = M.gene.map(gc)
M["co_tgt"] = [int(g in cotgt.get(int(h), set())) if pd.notna(h) else 0 for h, g in zip(M.hp, M.gene)]
TR = M[(M.hp.notna()) & ((M.ts == 1) | (M.site_free == 1)) & M.hp_chr.notna() & M.g_chr.notna() & (M.hp_chr != M.g_chr)].copy()

# profile: how often is each arm HIGH-position on its target genes?
for lab, s in [("guide sites", TR[(TR.role=="guide")&(TR.ts==1)]), ("passenger sites", TR[(TR.role=="passenger")&(TR.ts==1)])]:
    print(f"[{time.time()-t0:.0f}s] {lab}: rank1 {int((s.arm_rank==1).sum()):,} rank2 {int((s.arm_rank==2).sum()):,} "
          f"rank>=3 {int((s.arm_rank>=3).sum()):,} | leak-free&rank<=2 {int(((s.arm_rank<=2)&(s.co_tgt==0)).sum()):,}")


def demean(V, W, ai, gi, it=20):
    X = np.asarray(V, float).copy(); na, ng = int(ai.max()) + 1, int(gi.max()) + 1
    for _ in range(it):
        for idxs, m in ((ai, na), (gi, ng)):
            sw = np.bincount(idxs, weights=W, minlength=m); sw[sw == 0] = 1.0
            if X.ndim == 1: X -= (np.bincount(idxs, weights=W * X, minlength=m) / sw)[idxs]
            else:
                GS = np.zeros((m, X.shape[1])); np.add.at(GS, idxs, W[:, None] * X); X -= (GS / sw[:, None])[idxs]
    return X


def contrast(dd, B=2000):
    dd = dd.reset_index(drop=True); W = np.ones(len(dd))
    ng = int((dd.ts*(dd.role=="guide")).sum()); npg = int((dd.ts*(dd.role=="passenger")).sum())
    if ng < 20 or npg < 20: return None
    arms = dd.arm.to_numpy(); uarms = pd.unique(arms); idx = {a: k for k, a in enumerate(uarms)}
    col = np.array([idx[a] for a in arms]); ts = dd.ts.to_numpy(float)
    DVraw = np.zeros((len(dd), len(uarms))); DVraw[np.arange(len(dd)), col] = ts
    aic = dd.arm.astype("category").cat.codes.to_numpy(); gic = dd.gene.astype("category").cat.codes.to_numpy()
    DV = demean(DVraw, W, aic, gic); yr = demean(dd.reduced.to_numpy(float), W, aic, gic); DVsum = DV.sum(1)
    hpv = dd.hp.to_numpy(); pair = {}; g0 = np.zeros(len(uarms), bool)
    for a in uarms:
        rr = arms == a; pair.setdefault(int(hpv[rr][0]), []).append(idx[a])
        if dd.role.to_numpy()[rr][0] == "guide": g0[idx[a]] = True
    pairs = [tuple(v) for v in pair.values() if len(v) == 2]
    if not pairs: return None
    def coef(sg):
        XG = DV @ sg; XP = DVsum - XG
        return np.linalg.solve(np.array([[XG@XG, XG@XP], [XP@XG, XP@XP]]), np.array([XG@yr, XP@yr]))
    bG, bP = coef(g0.astype(float)); obs = bG - bP
    null = np.empty(B)
    for bb in range(B):
        sg = g0.copy(); fl = rng.integers(0, 2, len(pairs)).astype(bool)
        for (i1, i2), f in zip(pairs, fl):
            if f: sg[i1] = not sg[i1]; sg[i2] = not sg[i2]
        c = coef(sg.astype(float)); null[bb] = c[0] - c[1]
    return dict(gap=obs, tG=bG, tP=bP, p=(np.sum(null <= obs)+1)/(B+1), nhp=len(pairs), nGs=ng, nPs=npg)


def show(lab, o):
    print(f"  {lab:46s} " + ("skipped(sparse)" if o is None else
          f"hp={o['nhp']:>3} nG={o['nGs']:>5,} nP={o['nPs']:>5,} tG={o['tG']:+.5f} tP={o['tP']:+.5f} gap={o['gap']:+.5f} p={o['p']:.4f}"))


print(f"\n[{time.time()-t0:.0f}s] === matched-position, leak-free contrasts (TRANS) ===")
sf = TR[TR.site_free == 1]
hi_lf = TR[((TR.ts == 1) & (TR.arm_rank <= 2) & (TR.co_tgt == 0))]           # high-position, leak-free SITES
show("HIGH-position (rank<=2) + leak-free", contrast(pd.concat([hi_lf, sf])))
show("  rank==1 + leak-free", contrast(pd.concat([TR[(TR.ts==1)&(TR.arm_rank==1)&(TR.co_tgt==0)], sf])))
# comparable-abundance hairpins (both arms functional): |dAB| small
for thr in (1.0, 1.5, 2.0):
    show(f"comparable-abundance hairpins |dAB|<{thr}", contrast(TR[TR.dab < thr]))
show("comparable-abund |dAB|<2 + HIGH-pos + leak-free",
     contrast(pd.concat([TR[(TR.dab<2)&(TR.ts==1)&(TR.arm_rank<=2)&(TR.co_tgt==0)], TR[(TR.dab<2)&(TR.site_free==1)]])))
print(f"[{time.time()-t0:.0f}s] DONE")
