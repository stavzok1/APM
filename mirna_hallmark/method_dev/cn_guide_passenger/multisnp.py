"""Stronger instruments WITHOUT imputation: cross-fitted multi-SNP cis allelic score per miRNA arm.
Per arm: y = miRNA residualized on ancestry PCs; 5-fold cross-fitted Ridge(y ~ cis-SNP dosages) -> OOF score
(unbiased instrument, no winner's curse). OOF first-stage F. Then reduced form + target-specificity control
on arms reaching F>10. Tests: do STRONGER instruments reveal target-specificity (curated vs non-target)?
"""
import os, time, numpy as np, pandas as pd
os.environ["OMP_NUM_THREADS"] = "1"
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from scipy import stats
from mirna_hallmark.learned import data as LD
from mirna_hallmark import data_loaders as D

HERE = os.path.dirname(os.path.abspath(__file__)); t0 = time.time(); rng = np.random.default_rng(3); WIN = 1_000_000
G = np.load(f"{HERE}/geno_cis.npy"); parts = pd.read_csv(f"{HERE}/geno_participants.txt")["case"].tolist()
ann = pd.read_csv(f"{HERE}/geno_cis_annot.tsv", sep="\t", dtype={"chr": str, "pos": "Int64"})
dos = np.where(G >= 0, G, np.nan).astype(np.float32); cm = np.nanmean(dos, 0)
cr = np.mean(G >= 0, 0); maf = np.nanmean(dos, 0)/2; maf = np.minimum(maf, 1-maf)
ii = np.where(np.isnan(dos)); dos[ii] = np.take(cm, ii[1])
okS = (cr >= 0.95) & (maf >= 0.05) & ann.pos.notna().values
dos = dos[:, okS]; ann = ann[okS].reset_index(drop=True)
Z = (dos-dos.mean(0))/(dos.std(0)+1e-6); PC = TruncatedSVD(10, random_state=0).fit_transform(Z); PC=(PC-PC.mean(0))/(PC.std(0)+1e-6)
print(f"[{time.time()-t0:.0f}s] geno {dos.shape}, PCs done")

d = LD._load(); X, Y, C = d["X"], d["Y"], d["C"]; norm = lambda s: "-".join(str(s).split("-")[:3])
pmap = {p: i for i, p in enumerate(parts)}
cpe_by = {norm(i): v for i, v in zip(C.index, pd.to_numeric(C["CPE"], errors="coerce").values)}
E = D.load_hallmark_edges(); he = E[E.high_evidence]; tgt_of = he.groupby("miRNA").gene.apply(lambda s: set(s)).to_dict()
LCX = pd.read_csv(f"{HERE.rsplit('/method_dev',1)[0]}/output/learned/locus_context.tsv", sep="\t")
LCX["chr"] = LCX.chrom.str.replace("chr","",regex=False)
from mirna_hallmark.learned import instrument as INS
loc_mid = {r.locus_id:(str(r.chr),(int(r.start)+int(r.end))//2) for r in LCX.itertuples()}
Yc = {norm(c): j for j, c in enumerate(Y.columns)}; Xc = {norm(c): j for j, c in enumerate(X.columns)}
allP = [p for p in parts if p in Yc and p in Xc]
gi = np.array([pmap[p] for p in allP]); yj=np.array([Yc[p] for p in allP]); xj=np.array([Xc[p] for p in allP])
cpe = np.nan_to_num(np.array([cpe_by.get(p,np.nan) for p in allP]), nan=0.0)
PCg = PC[gi]; n=len(allP); Yv,Xv=Y.values,X.values; allg=list(Y.index)
A_pc = np.column_stack([np.ones(n), PCg]); A_c = np.column_stack([np.ones(n), PCg, cpe])
bychr = {c:(np.array(g.index), g.pos.astype(int).values) for c,g in ann.groupby("chr")}
def resid(y,A):
    m=np.isfinite(y); Ai=A[m]; b=np.linalg.lstsq(Ai,y[m],rcond=None)[0]; r=np.full_like(y,np.nan); r[m]=y[m]-Ai@b; return r
kf = KFold(5, shuffle=True, random_state=0)
HE_arms = sorted({r.miRNA for r in he.itertuples()} & set(X.index))

def oof_score(arm):
    loci=[l for l in INS._arm_active_loci(arm) if l in loc_mid]
    if not loci: return None
    ch,mid=loc_mid[loci[0]]
    if ch not in bychr: return None
    idx,pos=bychr[ch]; cis=idx[(pos>=mid-WIN)&(pos<=mid+WIN)]
    if len(cis)<5: return None
    xe=pd.to_numeric(pd.Series(Xv[X.index.get_loc(arm)][xj]),errors="coerce").values.astype(float)
    y=resid(xe,A_pc); m=np.isfinite(y)
    if m.sum()<400 or np.nanstd(y[m])<1e-6: return None
    Xg=dos[gi][:,cis]; y2=np.nan_to_num(y,nan=np.nanmedian(y[m]))
    score=np.zeros(n)
    for tr,te in kf.split(Xg):
        r=RidgeCV(alphas=[1,10,100,1000]).fit(Xg[tr]-Xg[tr].mean(0), y2[tr])
        score[te]=(Xg[te]-Xg[tr].mean(0))@r.coef_
    # OOF first-stage F: y ~ score | (already PC-resid)
    sc=resid(score,A_pc); yy=resid(xe,A_pc); mm=np.isfinite(sc)&np.isfinite(yy)
    r=np.corrcoef(sc[mm],yy[mm])[0,1]; F=(mm.sum()-2)*r**2/max(1-r**2,1e-9)
    return score, float(F), int(len(cis))

# scan all HE arms -> OOF F
rows=[]; scores={}
for i,arm in enumerate(HE_arms):
    o=oof_score(arm)
    if o is None: continue
    scores[arm]=o[0]; rows.append((arm,o[2],o[1]))
    if (i+1)%80==0: print(f"[{time.time()-t0:.0f}s] {i+1}/{len(HE_arms)}")
FS=pd.DataFrame(rows,columns=["arm","n_cis","oof_F"]).sort_values("oof_F",ascending=False)
FS.to_csv(f"{HERE}/multisnp_firststage.tsv",sep="\t",index=False)
print(f"\n[{time.time()-t0:.0f}s] === MULTI-SNP OOF first stage ({len(FS)} arms) ===")
print(f"  arms OOF-F>10: {int((FS.oof_F>10).sum())} | >20: {int((FS.oof_F>20).sum())} | >30: {int((FS.oof_F>30).sum())}")
print("  (single-SNP had 38 perm-sig F>10, only 1-2 F>20)")
print(FS.head(12).to_string(index=False,float_format=lambda x:f"{x:.1f}"))

# specificity on strong arms
strong=FS[FS.oof_F>10].arm.tolist()
def aligned(sc, gene):
    yv=pd.to_numeric(pd.Series(Yv[Y.index.get_loc(gene)][yj]),errors="coerce").values.astype(float)
    yr=resid(yv,A_c); sr=resid(sc,A_c); m=np.isfinite(yr)&np.isfinite(sr)
    if m.sum()<300 or np.nanstd(yr[m])<1e-9: return np.nan
    rho,_=stats.spearmanr(sr[m],yr[m]); return -rho
spec=[]
for arm in strong:
    if arm not in tgt_of: continue
    sc=scores[arm]; tg=[x for x in tgt_of[arm] if x in Y.index]
    if len(tg)<4: continue
    cur=np.array([aligned(sc,x) for x in tg]); cur=cur[np.isfinite(cur)]
    nt=[x for x in rng.choice(allg,size=min(80,len(allg)),replace=False) if x not in tgt_of[arm]][:len(tg)*3]
    ct=np.array([aligned(sc,x) for x in nt]); ct=ct[np.isfinite(ct)]
    if len(cur)<4 or len(ct)<4: continue
    spec.append((arm,len(cur),cur.mean(),ct.mean(),cur.mean()-ct.mean()))
SP=pd.DataFrame(spec,columns=["arm","n_tgt","cur","ctl","spec"])
SP.to_csv(f"{HERE}/multisnp_specificity.tsv",sep="\t",index=False)
print(f"\n[{time.time()-t0:.0f}s] === TARGET SPECIFICITY with STRONG multi-SNP instruments ({len(SP)} arms) ===")
if len(SP)>=4:
    t,p=stats.ttest_rel(SP.cur,SP.ctl); s=stats.binomtest((SP.spec>0).sum(),len(SP),0.5,alternative='greater').pvalue
    print(f"  curated aligned-pi {SP.cur.mean():+.4f} vs control {SP.ctl.mean():+.4f} | Δ={SP.spec.mean():+.4f} "
          f"paired t={t:+.2f} p={p:.3f} sign p={s:.3f} | curated>0 in {(SP.cur>0).mean():.0%} arms")
    print(SP.sort_values("spec",ascending=False).to_string(index=False,float_format=lambda x:f"{x:+.4f}"))
print(f"[{time.time()-t0:.0f}s] DONE")
