import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
e = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
t = pd.to_numeric(g.top_identity, errors="coerce").dropna()
print(f"top_identity: n={len(t)}")
for q in [0,.01,.05,.25,.5,.75,.95,.99,1]:
    print(f"   p{int(q*100):3d} = {t.quantile(q):+12.4f}")
print(f"\n   OUTSIDE [0,1]:  >1 -> {(t>1).sum()} ({(t>1).mean():.2%})   <0 -> {(t<0).sum()} ({(t<0).mean():.2%})")
print(f"   >2  -> {(t>2).sum()}     >10 -> {(t>10).sum()}    max {t.max():.1f}")
i = pd.to_numeric(e.identity, errors="coerce").dropna()
print(f"\nedge-card `identity`: n={len(i)}  min {i.min():+.4f}  max {i.max():+.4f}  >1 -> {(i>1).sum()}  <0 -> {(i<0).sum()}")
print("   -> per-gene sums of edge identity:")
s = e.assign(idn=pd.to_numeric(e.identity, errors='coerce')).groupby("gene").idn.sum().dropna()
print(f"      median {s.median():.4f}  p95 {s.quantile(.95):.4f}  max {s.max():.4f}  frac>1.01 {(s>1.01).mean():.1%}")
bad = g[pd.to_numeric(g.top_identity, errors="coerce")>1]
print(f"\ngenes with top_identity>1: {len(bad)}")
print(f"   their ceiling: median {pd.to_numeric(bad.ctx_ceiling,errors='coerce').median():+.4f}  "
      f"vs others {pd.to_numeric(g[pd.to_numeric(g.top_identity,errors='coerce')<=1].ctx_ceiling,errors='coerce').median():+.4f}")
print(f"   frac with ceiling<=0: {(pd.to_numeric(bad.ctx_ceiling,errors='coerce')<=0).mean():.1%} "
      f"vs {(pd.to_numeric(g[pd.to_numeric(g.top_identity,errors='coerce')<=1].ctx_ceiling,errors='coerce')<=0).mean():.1%}")
print("\n=== producer ===")
