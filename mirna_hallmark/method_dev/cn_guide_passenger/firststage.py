"""FEASIBILITY GATE: do HE miRNA arms have a strong germline cis-eQTL (first-stage F>10)?
For each arm: cis-SNPs within ±1Mb of its locus; residualize arm tumor expression + each SNP on ancestry
PCs; best-cis-SNP F. Winner's-curse guard: per-arm permutation null (shuffle expression, recompute max-F)
=> empirical p. Report arms with F>10 AND perm-p<0.05 = genuine strong instruments.
Ancestry PCs computed from the (genome-wide-distributed) cis-SNP matrix itself.
"""
import os, time, numpy as np, pandas as pd
os.environ["OMP_NUM_THREADS"] = "1"
from sklearn.decomposition import TruncatedSVD
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import instrument as INS
from mirna_hallmark import data_loaders as D

HERE = os.path.dirname(os.path.abspath(__file__)); t0 = time.time()
rng = np.random.default_rng(20260718); WIN = 1_000_000

G = np.load(f"{HERE}/geno_cis.npy")                                    # participants × cis-SNP int8 (-1 miss)
parts = pd.read_csv(f"{HERE}/geno_participants.txt")["case"].tolist()
ann = pd.read_csv(f"{HERE}/geno_cis_annot.tsv", sep="\t", dtype={"chr": str, "pos": "Int64"})
print(f"[{time.time()-t0:.0f}s] geno {G.shape} | {len(parts)} participants")

# QC + mean-impute dosage
dos = np.where(G >= 0, G, np.nan).astype(np.float32)
callrate = np.mean(G >= 0, 0); maf = np.nanmean(dos, 0) / 2; maf = np.minimum(maf, 1 - maf)
ok = (callrate >= 0.95) & (maf >= 0.05) & ann.pos.notna().values
dos = dos[:, ok]; ann = ann[ok].reset_index(drop=True)
colmean = np.nanmean(dos, 0); inds = np.where(np.isnan(dos)); dos[inds] = np.take(colmean, inds[1])
print(f"[{time.time()-t0:.0f}s] QC'd SNPs: {dos.shape[1]:,} (MAF>=0.05, call>=0.95)")

# ancestry PCs from standardized dosages (genome-wide-distributed)
Z = (dos - dos.mean(0)) / (dos.std(0) + 1e-6)
PC = TruncatedSVD(n_components=10, random_state=0).fit_transform(Z)
PC = (PC - PC.mean(0)) / (PC.std(0) + 1e-6)
A = np.column_stack([np.ones(len(parts)), PC])                        # PC covariate block
Ainv = np.linalg.pinv(A)
def resid(M): return M - A @ (Ainv @ M)
dosr = resid(dos)                                                     # SNPs residualized on PCs (once)
dosr /= (np.linalg.norm(dosr, axis=0) + 1e-9)
print(f"[{time.time()-t0:.0f}s] PCs computed; SNPs residualized")

# arm -> locus coords (chr,pos midpoint) via active loci + locus_context
LCX = pd.read_csv(f"{HERE.rsplit('/method_dev',1)[0]}/output/learned/locus_context.tsv", sep="\t")
LCX["chr"] = LCX.chrom.str.replace("chr", "", regex=False)
loc_mid = {r.locus_id: (str(r.chr), (int(r.start)+int(r.end))//2) for r in LCX.itertuples()}
X = LD._load()["X"]; pmap = {p: i for i, p in enumerate(parts)}
# align arm-expression columns to geno participants
Xcols = ["-".join(str(c).split("-")[:3]) for c in X.columns]
gi = np.array([pmap.get(c, -1) for c in Xcols]); keepP = gi >= 0
gi = gi[keepP]                                                        # geno-row index per kept X column
E = D.load_hallmark_edges(); HE_arms = sorted({r.miRNA for r in E[E.high_evidence].itertuples()} & set(X.index))
n = len(gi)

# per-chrom SNP index for fast cis lookup
bychr = {c: (np.array(g.index), g.pos.astype(int).values) for c, g in ann.groupby("chr")}
NPERM = 200
rows = []
for arm in HE_arms:
    loci = [l for l in INS._arm_active_loci(arm) if l in loc_mid]
    if not loci: continue
    ch, mid = loc_mid[loci[0]]
    if ch not in bychr: continue
    idx, pos = bychr[ch]
    m = (pos >= mid-WIN) & (pos <= mid+WIN)
    cis = idx[m]
    if len(cis) < 3: continue
    y = pd.to_numeric(X.loc[arm].values[keepP], errors="coerce").astype(np.float32)
    if np.isnan(y).mean() > 0.2 or np.nanstd(y) < 1e-6: continue
    y = np.nan_to_num(y, nan=np.nanmedian(y))
    yr = y[None, :] - 0                                                # residualize y on PCs (aligned to geno rows)
    # build y residualized on PCs (using geno-row-aligned A)
    Ai = A[gi]; yr = y - Ai @ (np.linalg.pinv(Ai) @ y)
    Xc = dosr[gi][:, cis]                                             # SNPs (already PC-residualized, per-full-cohort) at cis, geno rows
    Xc = Xc - Xc.mean(0); Xc /= (np.linalg.norm(Xc, axis=0) + 1e-9)
    yn = (yr - yr.mean()); yn /= (np.linalg.norm(yn) + 1e-9)
    r = Xc.T @ yn                                                     # partial corr per SNP
    F = (n - 12) * r**2 / np.maximum(1 - r**2, 1e-9)
    fmax = float(F.max()); best = cis[int(F.argmax())]
    # permutation null: shuffle y, recompute max-F
    permmax = np.empty(NPERM)
    for b in range(NPERM):
        yp = yn[rng.permutation(n)]
        rp = Xc.T @ yp; permmax[b] = ((n-12)*rp**2/np.maximum(1-rp**2,1e-9)).max()
    pval = (np.sum(permmax >= fmax) + 1) / (NPERM + 1)
    rows.append((arm, ch, mid, len(cis), fmax, ann.loc[best,"probeid"], int(ann.loc[best,"pos"]), pval))

R = pd.DataFrame(rows, columns=["arm","chr","locus_mid","n_cis","F_best","best_snp","best_pos","perm_p"])
R.to_csv(f"{HERE}/firststage_cis_eqtl.tsv", sep="\t", index=False)
sig = R[(R.F_best>10) & (R.perm_p<0.05)]
print(f"\n[{time.time()-t0:.0f}s] === FIRST-STAGE FEASIBILITY ({len(R)} HE arms scanned, n={n}) ===")
print(f"  arms with best-cis F>10:              {int((R.F_best>10).sum())} ({(R.F_best>10).mean():.0%})")
print(f"  arms with perm-p<0.05 (genuine eQTL): {int((R.perm_p<0.05).sum())} ({(R.perm_p<0.05).mean():.0%})")
print(f"  ⭐ arms with F>10 AND perm-p<0.05:     {len(sig)} ({len(sig)/len(R):.0%})  <- usable strong instruments")
print(f"  F>20 & perm-sig: {int(((R.F_best>20)&(R.perm_p<0.05)).sum())} | F>30: {int(((R.F_best>30)&(R.perm_p<0.05)).sum())}")
print("  top instruments:"); print(R.sort_values("F_best",ascending=False).head(12).to_string(index=False))
