"""Co-cluster leak (Q2) + strongest cases + canonical positive controls + seedless negative control.
Optimized: co-targeting precomputed per-hairpin (159 sets), mapped over TR only."""
import time, numpy as np, pandas as pd
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
clmembers = {}
for a, c in CL.items(): clmembers.setdefault(c, []).append(a)
hp_co = {h: [a for a in clmembers.get(CL.get(a1), []) if a not in (a1, a2) and a in detarms]
         for h, (a1, a2) in hp_arms.items()}
# genes co-targeted by a hairpin's co-cluster members (precompute per hairpin)
ts_by_arm = TSraw.groupby("arm").gene.apply(set).to_dict()
cotgt_genes = {h: set().union(*[ts_by_arm.get(a, set()) for a in co]) if co else set() for h, co in hp_co.items()}
nsolo = sum(1 for h in hp_arms if not hp_co[h])
LC = pd.read_csv("mirna_hallmark/output/learned/locus_context.tsv", sep="\t").set_index("locus_id")["chrom"]
hp_chrom = {h: next((LC.get(c) for c in li if c in LC.index), None) for h, li in hp_loci.items()}
gc = pd.read_csv("data/gencode.v49.annotation.gtf.csv", usecols=["seqname", "feature", "gene_name"])
gc = gc[gc.feature == "gene"].drop_duplicates("gene_name").set_index("gene_name")["seqname"]
M["hp_chr"] = M.hp.map(hp_chrom); M["g_chr"] = M.gene.map(gc)
M["solo"] = M.hp.map({h: (len(hp_co[h]) == 0) for h in hp_arms})
TR = M[(M.hp.notna()) & ((M.ts == 1) | (M.site_free == 1)) & M.hp_chr.notna() & M.g_chr.notna()
       & (M.hp_chr != M.g_chr)].copy()
TR["co_tgt"] = [int(g in cotgt_genes.get(int(h), set())) for h, g in zip(TR.hp, TR.gene)]
gco = TR[(TR.role == "guide") & (TR.ts == 1)].co_tgt.mean(); pco = TR[(TR.role == "passenger") & (TR.ts == 1)].co_tgt.mean()
print(f"[{time.time()-t0:.0f}s] hairpins {len(hp_arms)} | SOLO {nsolo} clustered {len(hp_arms)-nsolo} | "
      f"co-cluster co-targeting of site genes: guide {gco:.1%} vs passenger {pco:.1%}")


def demean(V, W, ai, gi, it=20):
    X = np.asarray(V, float).copy(); na, ng = int(ai.max()) + 1, int(gi.max()) + 1
    for _ in range(it):
        for idxs, m in ((ai, na), (gi, ng)):
            sw = np.bincount(idxs, weights=W, minlength=m); sw[sw == 0] = 1.0
            if X.ndim == 1: X -= (np.bincount(idxs, weights=W * X, minlength=m) / sw)[idxs]
            else:
                GS = np.zeros((m, X.shape[1])); np.add.at(GS, idxs, W[:, None] * X); X -= (GS / sw[:, None])[idxs]
    return X


def contrast(dd, B=1000):
    dd = dd.reset_index(drop=True); W = np.ones(len(dd))
    if int((dd.ts*(dd.role=="guide")).sum()) < 25 or int((dd.ts*(dd.role=="passenger")).sum()) < 25: return None
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
    return dict(gap=obs, tG=bG, tP=bP, p=(np.sum(null <= obs)+1)/(B+1), nhp=len(pairs), nGs=int((dd.ts*(dd.role=="guide")).sum()))


def show(lab, o):
    print(f"  {lab:44s} " + ("skipped(sparse)" if o is None else
          f"hp={o['nhp']:>3} nGs={o['nGs']:>6,} tG={o['tG']:+.5f} tP={o['tP']:+.5f} gap={o['gap']:+.5f} p={o['p']:.4f}"))


print(f"\n[{time.time()-t0:.0f}s] === Q2 co-cluster leak (TRANS, role-swap null) ===")
show("baseline (all)", contrast(TR))
show("SOLO hairpins only", contrast(TR[TR.solo]))
show("EXCLUDE co-cluster-co-targeted genes", contrast(TR[TR.co_tgt == 0]))
show("clustered hairpins only", contrast(TR[TR.solo == False]))

print(f"\n[{time.time()-t0:.0f}s] === strongest cases: guide arms by mean reduced on their TRANS site genes (>=20 sites) ===")
gs = TR[(TR.role == "guide") & (TR.ts == 1)]
prof = gs.groupby("arm").agg(n=("reduced", "size"), mean_red=("reduced", "mean")).query("n>=20").sort_values("mean_red")
print(prof.head(12).to_string(float_format=lambda x: f"{x:+.4f}"))

print(f"\n[{time.time()-t0:.0f}s] === canonical POSITIVE controls (is the miRNA a guide? reduced on target?) ===")
CANON = [("hsa-miR-16-5p", "BCL2"), ("hsa-miR-16-5p", "CCND1"), ("hsa-miR-200c-3p", "ZEB1"),
         ("hsa-miR-93-5p", "PTEN"), ("hsa-miR-106b-5p", "PTEN"), ("hsa-miR-21-5p", "PTEN"),
         ("hsa-miR-155-5p", "SOCS1"), ("hsa-miR-17-5p", "CDKN1A")]
for a, g in CANON:
    r = TR[(TR.arm == a) & (TR.gene == g)]
    rr = r.reduced.iloc[0] if len(r) else np.nan
    rl = role.get(a, "n/a"); inpanel = g in set(TR.gene)
    print(f"  {a:18s}->{g:8s} role={rl:9s} reduced={rr if len(r) else float('nan'):+.4f} "
          f"(in panel gene set: {inpanel}, edge rows: {len(r)})")
print(f"[{time.time()-t0:.0f}s] DONE")
