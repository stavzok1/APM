"""Reduced-form MR: for arms with a genuine germline cis-eQTL, does the instrument SNP move the arm's
CURATED targets in the repression direction? Sign-align: SNP up -> (if gamma>0) miRNA up -> target DOWN,
so a real repressive edge has sign(pi_reduced) = -sign(gamma). Test curated edges vs matched control genes.
C = [PCs, CPE, target_cn]. Flag DLK1-DIO3 imprinted-cluster arms (exclusion compromised by co-transcription).
"""
import os, time, numpy as np, pandas as pd
os.environ["OMP_NUM_THREADS"] = "1"
from sklearn.decomposition import TruncatedSVD
from scipy import stats
from mirna_hallmark.learned import data as LD
from mirna_hallmark import data_loaders as D

HERE = os.path.dirname(os.path.abspath(__file__)); t0 = time.time(); rng = np.random.default_rng(7)
G = np.load(f"{HERE}/geno_cis.npy"); parts = pd.read_csv(f"{HERE}/geno_participants.txt")["case"].tolist()
ann = pd.read_csv(f"{HERE}/geno_cis_annot.tsv", sep="\t"); pid2col = {p: i for i, p in enumerate(ann.probeid)}
FS = pd.read_csv(f"{HERE}/firststage_cis_eqtl.tsv", sep="\t")
sig = FS[(FS.F_best > 10) & (FS.perm_p < 0.05)].copy()
# DLK1-DIO3 imprinted cluster (chr14q32 ~100.8-101.6Mb)
sig["imprinted"] = (sig.chr.astype(str) == "14") & (sig.locus_mid.between(100_700_000, 101_700_000))
print(f"[{time.time()-t0:.0f}s] instrumented arms: {len(sig)} | DLK1-DIO3 imprinted: {int(sig.imprinted.sum())} | clean solo: {int((~sig.imprinted).sum())}")

dos = np.where(G >= 0, G, np.nan).astype(np.float32)
cm = np.nanmean(dos, 0); ii = np.where(np.isnan(dos)); dos[ii] = np.take(cm, ii[1])
Z = (dos - dos.mean(0)) / (dos.std(0) + 1e-6)
PC = TruncatedSVD(n_components=10, random_state=0).fit_transform(Z); PC = (PC-PC.mean(0))/(PC.std(0)+1e-6)

d = LD._load(); X, Y = d["X"], d["Y"]; C = d["C"]
norm = lambda s: "-".join(str(s).split("-")[:3])
xcol = [norm(c) for c in X.columns]; ycol = [norm(c) for c in Y.columns]
pmap = {p: i for i, p in enumerate(parts)}
# CPE aligned to parts
cpe_by = {norm(i): v for i, v in zip(C.index, pd.to_numeric(C["CPE"], errors="coerce").values)}
E = D.load_hallmark_edges(); he = E[E.high_evidence]
tgt_of = he.groupby("miRNA").gene.apply(lambda s: sorted(set(s))).to_dict()

def geno_col(snp): return dos[:, pid2col[snp]]
def resid(y, A):
    m = np.isfinite(y); Ai = A[m]; b = np.linalg.lstsq(Ai, y[m], rcond=None)[0]; r = y.copy(); r[m] = y[m]-Ai@b; r[~m]=np.nan; return r

# participant alignment: geno rows
Yc = {norm(c): j for j, c in enumerate(Y.columns)}
Xc = {norm(c): j for j, c in enumerate(X.columns)}
allP = [p for p in parts if p in Yc and p in Xc]
gi = np.array([pmap[p] for p in allP]); yj = np.array([Yc[p] for p in allP]); xj = np.array([Xc[p] for p in allP])
cpe = np.array([cpe_by.get(p, np.nan) for p in allP]); cpe = np.nan_to_num(cpe, nan=np.nanmedian(cpe))
PCg = PC[gi]; n = len(allP); print(f"[{time.time()-t0:.0f}s] MR n={n}")
Yv = Y.values; Xv = X.values

def target_cn_vec(gene):
    from mirna_hallmark.learned import data as LDx
    tc = LDx._target_cn([gene]).get(gene)
    if tc is None: return np.zeros(n)
    m = {norm(i): v for i, v in zip(tc.index, pd.to_numeric(tc, errors="coerce").values)}
    return np.nan_to_num(np.array([m.get(p, np.nan) for p in allP]), nan=0.0)

rows = []
for r in sig.itertuples():
    arm = r.arm; snp = r.best_snp
    if arm not in X.index or arm not in tgt_of: continue
    g_snp = geno_col(snp)[gi]
    xe = pd.to_numeric(pd.Series(Xv[X.index.get_loc(arm)][xj]), errors="coerce").values.astype(float)
    A_fs = np.column_stack([np.ones(n), PCg]);
    gamma = np.polyfit(resid(g_snp, A_fs), resid(xe, A_fs), 1)[0]  # first-stage sign
    sgn = np.sign(gamma) if gamma != 0 else 1.0
    for gene in tgt_of[arm]:
        if gene not in Y.index: continue
        yv = pd.to_numeric(pd.Series(Yv[Y.index.get_loc(gene)][yj]), errors="coerce").values.astype(float)
        A = np.column_stack([np.ones(n), PCg, cpe, target_cn_vec(gene)])
        gr = resid(g_snp, A); yr = resid(yv, A); m = np.isfinite(gr) & np.isfinite(yr)
        if m.sum() < 300: continue
        rho, p = stats.spearmanr(gr[m], yr[m])
        rows.append((arm, gene, snp, float(gamma), rho, -sgn*rho, p, bool(r.imprinted)))
M = pd.DataFrame(rows, columns=["arm","gene","snp","gamma","pi_rho","pi_aligned","p","imprinted"])
M.to_csv(f"{HERE}/mr_reduced_form.tsv", sep="\t", index=False)
print(f"[{time.time()-t0:.0f}s] curated instrumented edges tested: {len(M)}")

def summ(lab, sub):
    if len(sub) < 5: print(f"  {lab:34s} n={len(sub)} (too few)"); return
    t, pt = stats.ttest_1samp(sub.pi_aligned, 0)
    frac = (sub.pi_aligned > 0).mean()
    print(f"  {lab:34s} n={len(sub):>4}  mean aligned-pi {sub.pi_aligned.mean():+.4f}  "
          f"frac repression-dir {frac:.2%}  t={t:+.2f} p={pt:.3f}")
print(f"\n[{time.time()-t0:.0f}s] === reduced-form MR: is SNP->target repression-directed for curated edges? ===")
print("  (aligned-pi>0 = repression-consistent: instrument that raises the miRNA lowers the target)")
summ("ALL instrumented curated edges", M)
summ("CLEAN (non-imprinted) arms", M[~M.imprinted])
summ("DLK1-DIO3 imprinted arms", M[M.imprinted])
# control: same SNPs vs RANDOM (non-target) genes
allg = list(Y.index)
ctrl = []
for r in sig.itertuples():
    if r.arm not in tgt_of: continue
    g_snp = geno_col(r.best_snp)[gi]; sgn = 1.0
    for gene in rng.choice(allg, size=min(8, len(tgt_of.get(r.arm,[]))), replace=False):
        yv = pd.to_numeric(pd.Series(Yv[Y.index.get_loc(gene)][yj]), errors="coerce").values.astype(float)
        A = np.column_stack([np.ones(n), PCg, cpe]); gr=resid(g_snp,A); yr=resid(yv,A); m=np.isfinite(gr)&np.isfinite(yr)
        if m.sum()<300: continue
        rho,_ = stats.spearmanr(gr[m], yr[m]); ctrl.append(abs(rho))
print(f"  CONTROL |pi| (random genes): mean {np.mean(ctrl):.4f} vs curated |pi| mean {M.pi_rho.abs().mean():.4f}")
print(f"[{time.time()-t0:.0f}s] DONE")
