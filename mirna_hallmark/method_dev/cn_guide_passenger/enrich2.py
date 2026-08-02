"""Corrected context readouts.
 A) STRATIFY the guide/passenger contrast by the guide's DOMINANCE in the gene's regulatory profile
    (guide = the gene's TOP regulator vs a MINOR one). Full contrast + role-swap null per stratum.
 B) FAMILY-level aggregation with the family defined by the GUIDE arm (the two arms of a hairpin are
    different seed families, so the hairpin is assigned to its guide's family), grouping hairpins.
Loads the enriched table columns rebuilt here. TRANS only.
"""
import time, numpy as np, pandas as pd
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as INS
from mirna_hallmark.learned import families as FAM

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
XALL = LD._load()["X"]; Xn = XALL.apply(pd.to_numeric, errors="coerce")
AB = Xn.median(axis=1); RPM = np.power(2.0, AB) - 1
loc = {a: tuple(INS._arm_active_loci(a)) for a in M.arm.unique()}
byl = {}
for a, l in loc.items():
    if l: byl.setdefault(l, []).append(a)
prs = {l: v for l, v in byl.items() if len(v) == 2}
hp_of, hp_loci, role, hp_guide = {}, {}, {}, {}
for hid, (l, (a1, a2)) in enumerate(prs.items()):
    if not (np.isfinite(AB.get(a1, np.nan)) and np.isfinite(AB.get(a2, np.nan))): continue
    hp_of[a1] = hp_of[a2] = hid; hp_loci[hid] = list(l)
    g_, p_ = (a1, a2) if AB[a1] >= AB[a2] else (a2, a1)
    role[g_], role[p_] = "guide", "passenger"; hp_guide[hid] = g_
M["hp"] = M.arm.map(hp_of); M["role"] = M.arm.map(role)
FAMILY = FAM.family_of(list(M.arm.dropna().unique())).to_dict()
hp_gfam = {h: FAMILY.get(g) for h, g in hp_guide.items()}
M["hp_gfam"] = M.hp.map(hp_gfam)
LC = pd.read_csv("mirna_hallmark/output/learned/locus_context.tsv", sep="\t").set_index("locus_id")["chrom"]
hp_chrom = {h: next((LC.get(c) for c in li if c in LC.index), None) for h, li in hp_loci.items()}
gc = pd.read_csv("data/gencode.v49.annotation.gtf.csv", usecols=["seqname", "feature", "gene_name"])
gc = gc[gc.feature == "gene"].drop_duplicates("gene_name").set_index("gene_name")["seqname"]
M["hp_chr"] = M.hp.map(hp_chrom); M["g_chr"] = M.gene.map(gc)
# guide dominance rank within a gene's detectable site-bearing regulators
detset = set(AB[AB > FLOOR].index)
TSd = TSraw[TSraw.arm.isin(detset)].copy(); TSd["rpm"] = TSd.arm.map(RPM).fillna(0.0)
TSd["rk"] = TSd.groupby("gene").rpm.rank(ascending=False, method="min")
rank_map = {(a, g): rk for a, g, rk in zip(TSd.arm, TSd.gene, TSd.rk)}
M["guide_rank"] = [rank_map.get((a, g), np.nan) for a, g in zip(M.arm, M.gene)]
TR = M[(M.hp.notna()) & ((M.ts == 1) | (M.site_free == 1)) & M.hp_chr.notna() & M.g_chr.notna()
       & (M.hp_chr != M.g_chr)].copy()


def demean(V, W, ai, gi, it=20):
    X = np.asarray(V, float).copy(); na, ng = int(ai.max()) + 1, int(gi.max()) + 1
    for _ in range(it):
        for idxs, m in ((ai, na), (gi, ng)):
            sw = np.bincount(idxs, weights=W, minlength=m); sw[sw == 0] = 1.0
            if X.ndim == 1: X -= (np.bincount(idxs, weights=W * X, minlength=m) / sw)[idxs]
            else:
                GS = np.zeros((m, X.shape[1])); np.add.at(GS, idxs, W[:, None] * X); X -= (GS / sw[:, None])[idxs]
    return X


def contrast(dd, B=2000, minsite=25):
    dd = dd.reset_index(drop=True); W = np.ones(len(dd))
    if int((dd.ts*(dd.role=="guide")).sum()) < minsite or int((dd.ts*(dd.role=="passenger")).sum()) < minsite:
        return None
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
    return dict(gap=obs, tG=bG, tP=bP, p=(np.sum(null <= obs)+1)/(B+1), nhp=len(pairs),
                nGs=int((dd.ts*(dd.role=="guide")).sum()))


def show(lab, o):
    print(f"  {lab:40s} " + ("skipped" if o is None else
          f"hp={o['nhp']:>3} nGsite={o['nGs']:>6,} tG={o['tG']:+.5f} tP={o['tP']:+.5f} gap={o['gap']:+.5f} p={o['p']:.4f}"))


# ===== Readout A: contrast stratified by guide dominance in the gene profile =====
print(f"[{time.time()-t0:.0f}s] === A: guide/passenger contrast vs guide DOMINANCE (gene's top vs minor regulator) ===")
base_pass = TR[(TR.role == "passenger")] ; base_sf = TR[TR.site_free == 1]
for lab, gmask in [("guide = gene's TOP regulator (rank1)", TR.guide_rank == 1),
                   ("guide = rank 2", TR.guide_rank == 2),
                   ("guide = MINOR (rank>=3)", TR.guide_rank >= 3)]:
    gs = TR[(TR.role == "guide") & (TR.ts == 1) & gmask]
    sub = pd.concat([gs, base_pass, base_sf]).drop_duplicates()
    show(lab, contrast(sub))

# ===== Readout B: family aggregation (family = GUIDE arm's family), grouping hairpins =====
print(f"\n[{time.time()-t0:.0f}s] === B: FAMILY-level aggregation (family = guide's seed family) ===")
famhp = TR.groupby("hp_gfam").hp.nunique().sort_values(ascending=False)
big = famhp[famhp >= 3]
print(f"  guide-families with >=3 both-arm hairpins: {len(big)}")
rows = []
for fam, nhp in big.items():
    o = contrast(TR[TR.hp_gfam == fam], minsite=15)
    if o: rows.append((fam, o["nhp"], o["gap"], o["tG"], o["tP"], o["p"]))
FB = pd.DataFrame(rows, columns=["guide_family", "nhp", "gap", "tG", "tP", "p"]).sort_values("gap")
print(FB.to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
if len(FB): print(f"  families gap<0: {int((FB.gap<0).sum())}/{len(FB)} | p<0.05: {int((FB.p<0.05).sum())}")
print(f"[{time.time()-t0:.0f}s] DONE")
