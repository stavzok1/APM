"""DECISIVE control: are the germline instrument reduced forms TARGET-SPECIFIC (repression-directed for the
miRNA's CURATED targets) or indiscriminate (pleiotropy)? Per instrumented arm, compare aligned-pi on curated
targets vs matched RANDOM non-target genes (same SNP, same sign-alignment). Paired arm-level test.
If curated >> control -> specific (supports miRNA-mediated). If curated ~ control -> pleiotropy/confound.
"""
import os, time, numpy as np, pandas as pd
os.environ["OMP_NUM_THREADS"] = "1"
from sklearn.decomposition import TruncatedSVD
from scipy import stats
from mirna_hallmark.learned import data as LD
from mirna_hallmark import data_loaders as D

HERE = os.path.dirname(os.path.abspath(__file__)); t0 = time.time(); rng = np.random.default_rng(11)
G = np.load(f"{HERE}/geno_cis.npy"); parts = pd.read_csv(f"{HERE}/geno_participants.txt")["case"].tolist()
ann = pd.read_csv(f"{HERE}/geno_cis_annot.tsv", sep="\t"); pid2col = {p: i for i, p in enumerate(ann.probeid)}
FS = pd.read_csv(f"{HERE}/firststage_cis_eqtl.tsv", sep="\t")
sig = FS[(FS.F_best > 10) & (FS.perm_p < 0.05)].copy()
sig["imprinted"] = (sig.chr.astype(str) == "14") & (sig.locus_mid.between(100_700_000, 101_700_000))
dos = np.where(G >= 0, G, np.nan).astype(np.float32); cm = np.nanmean(dos, 0)
ii = np.where(np.isnan(dos)); dos[ii] = np.take(cm, ii[1])
Z = (dos-dos.mean(0))/(dos.std(0)+1e-6); PC = TruncatedSVD(10, random_state=0).fit_transform(Z); PC=(PC-PC.mean(0))/(PC.std(0)+1e-6)
d = LD._load(); X, Y, C = d["X"], d["Y"], d["C"]; norm = lambda s: "-".join(str(s).split("-")[:3])
pmap = {p: i for i, p in enumerate(parts)}
cpe_by = {norm(i): v for i, v in zip(C.index, pd.to_numeric(C["CPE"], errors="coerce").values)}
E = D.load_hallmark_edges(); he = E[E.high_evidence]; tgt_of = he.groupby("miRNA").gene.apply(lambda s: set(s)).to_dict()
Yc = {norm(c): j for j, c in enumerate(Y.columns)}; Xc = {norm(c): j for j, c in enumerate(X.columns)}
allP = [p for p in parts if p in Yc and p in Xc]
gi = np.array([pmap[p] for p in allP]); yj = np.array([Yc[p] for p in allP]); xj = np.array([Xc[p] for p in allP])
cpe = np.array([cpe_by.get(p, np.nan) for p in allP]); cpe = np.nan_to_num(cpe, nan=np.nanmedian(cpe))
PCg = PC[gi]; n = len(allP); Yv, Xv = Y.values, X.values; allg = list(Y.index)
A0 = np.column_stack([np.ones(n), PCg, cpe])
def resid(y, A):
    m = np.isfinite(y); Ai = A[m]; b = np.linalg.lstsq(Ai, y[m], rcond=None)[0]; r = np.full_like(y, np.nan); r[m]=y[m]-Ai@b; return r
print(f"[{time.time()-t0:.0f}s] n={n}, {len(sig)} instrumented arms")

def aligned_pi(g_snp_r, gene, sgn):
    yv = pd.to_numeric(pd.Series(Yv[Y.index.get_loc(gene)][yj]), errors="coerce").values.astype(float)
    yr = resid(yv, A0); m = np.isfinite(g_snp_r) & np.isfinite(yr)
    if m.sum() < 300 or np.nanstd(yr[m]) < 1e-9: return np.nan
    rho, _ = stats.spearmanr(g_snp_r[m], yr[m]); return -sgn*rho

rows = []
for r in sig.itertuples():
    arm = r.arm
    if arm not in X.index or arm not in tgt_of: continue
    g = dos[:, pid2col[r.best_snp]][gi]; gr = resid(g, A0)
    xe = pd.to_numeric(pd.Series(Xv[X.index.get_loc(arm)][xj]), errors="coerce").values.astype(float)
    xr = resid(xe, A0); mm = np.isfinite(gr) & np.isfinite(xr)
    gamma = np.polyfit(gr[mm], xr[mm], 1)[0]; sgn = np.sign(gamma) if gamma != 0 else 1.0
    tg = [x for x in tgt_of[arm] if x in Y.index]
    if len(tg) < 4: continue
    cur = np.array([aligned_pi(gr, x, sgn) for x in tg]); cur = cur[np.isfinite(cur)]
    nont = [x for x in rng.choice(allg, size=min(60, len(allg)), replace=False) if x not in tgt_of[arm]][:len(tg)*3]
    ctl = np.array([aligned_pi(gr, x, sgn) for x in nont]); ctl = ctl[np.isfinite(ctl)]
    if len(cur) < 4 or len(ctl) < 4: continue
    rows.append((arm, bool(r.imprinted), len(cur), cur.mean(), len(ctl), ctl.mean(), cur.mean()-ctl.mean()))
S = pd.DataFrame(rows, columns=["arm","imprinted","n_tgt","cur_mean","n_ctl","ctl_mean","specificity"])
S.to_csv(f"{HERE}/mr_specificity.tsv", sep="\t", index=False)
cl = S[~S.imprinted]
print(f"\n[{time.time()-t0:.0f}s] === TARGET SPECIFICITY (curated aligned-pi vs matched non-target control) ===")
print(f"  clean arms: {len(cl)}")
print(f"  curated mean aligned-pi:  {cl.cur_mean.mean():+.4f}  (frac arms>0 {(cl.cur_mean>0).mean():.0%})")
print(f"  control mean aligned-pi:  {cl.ctl_mean.mean():+.4f}  (frac arms>0 {(cl.ctl_mean>0).mean():.0%})")
t,p = stats.ttest_rel(cl.cur_mean, cl.ctl_mean); sgn = stats.binomtest((cl.specificity>0).sum(), len(cl), 0.5, alternative='greater').pvalue
print(f"  PAIRED curated vs control: Δ={cl.specificity.mean():+.4f}  t={t:+.2f} p={p:.3f}  sign-test p={sgn:.3f}")
print(f"  ⇒ if curated>>control: target-specific (miRNA-mediated). if ≈: pleiotropy.")
print(S.sort_values("specificity",ascending=False).to_string(index=False, float_format=lambda x:f"{x:+.4f}"))
