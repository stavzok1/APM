"""Extract cis-SNP germline genotypes (Affy SNP6 birdseed) for the miRNA cohort.
Probe order is identical across all birdseed files ⇒ subset by precomputed line indices (_cis_idx.npy).
QC: birdseed confidence > 0.1 -> missing (-1). Output: participant × cis-SNP int8 matrix + SNP annot.
"""
import os, glob, time, numpy as np, pandas as pd
os.environ["OMP_NUM_THREADS"] = "1"
from mirna_hallmark.learned import data as LD

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = "data/SNV/germline_array/raw"
t0 = time.time()

# file -> Case ID (12-char TCGA barcode)
ss = pd.read_csv("data/SNV/germline_array/gdc_sample_sheet.2026-06-07.tsv", sep="\t")
file2case = dict(zip(ss["File Name"], ss["Case ID"]))
# miRNA cohort participants
xcols = list(LD._load()["X"].columns)
norm = lambda s: "-".join(str(s).split("-")[:3])
cohort = set(norm(c) for c in xcols)

files = sorted(glob.glob(f"{RAW}/*.birdseed.data.txt"))
# keep cohort files, dedupe by case (first)
keep, seen = [], set()
for f in files:
    case = file2case.get(os.path.basename(f))
    if case in cohort and case not in seen:
        keep.append((f, case)); seen.add(case)
print(f"[{time.time()-t0:.0f}s] cohort birdseed files: {len(keep)} (of {len(files)})")

cis_idx = np.load(f"{HERE}/_cis_idx.npy")
probes = pd.read_csv(keep[0][0], sep="\t", skiprows=2, header=None, usecols=[0])[0].values[cis_idx]
G = np.full((len(keep), len(cis_idx)), -1, dtype=np.int8)
cases = []
for i, (f, case) in enumerate(keep):
    d = pd.read_csv(f, sep="\t", skiprows=2, header=None, usecols=[1, 2],
                    names=["call", "conf"], dtype={"call": "int8", "conf": "float32"})
    call = d["call"].values[cis_idx].copy(); conf = d["conf"].values[cis_idx]
    call[conf > 0.1] = -1            # QC: low-confidence -> missing
    G[i] = call; cases.append(case)
    if (i + 1) % 200 == 0: print(f"[{time.time()-t0:.0f}s]  {i+1}/{len(keep)}")

# annotation for the cis SNPs (aligned to columns)
ann = pd.read_csv("data/SNV/germline_array/snp6.na35.remap.hg38.subset.txt.gz", sep="\t",
                  dtype={"chr": str, "pos": "Int64"}).set_index("probeid")
cis_ann = ann.reindex(probes)[["chr", "pos"]].reset_index()
cis_ann.columns = ["probeid", "chr", "pos"]
cis_ann["chr"] = cis_ann.chr.str.replace("chr", "", regex=False)

np.save(f"{HERE}/geno_cis.npy", G)
pd.Series(cases, name="case").to_csv(f"{HERE}/geno_participants.txt", index=False)
cis_ann.to_csv(f"{HERE}/geno_cis_annot.tsv", sep="\t", index=False)
# QC summary
callrate = (G >= 0).mean(0)
dos = np.where(G >= 0, G, np.nan)
maf = np.nanmean(dos, 0) / 2.0; maf = np.minimum(maf, 1 - maf)
print(f"[{time.time()-t0:.0f}s] DONE. G={G.shape} | median call-rate {np.median(callrate):.3f} | "
      f"SNPs MAF>=0.05 & callrate>=0.95: {int(((maf>=0.05)&(callrate>=0.95)).sum()):,}")
PY = None
