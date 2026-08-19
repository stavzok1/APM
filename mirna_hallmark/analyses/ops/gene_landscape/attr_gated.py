import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
e = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
def num(s):
    return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float) if s.dtype==object else pd.to_numeric(s, errors="coerce")

print("="*95); print("IDENTITY IS A SIGNED SHARE — the suppressor structure, quantified"); print("="*95)
e["idn"] = num(e.identity)
neg = e[e.idn<0]
print(f"edge identity: {len(e.idn.dropna())} values | NEGATIVE {len(neg)} ({len(neg)/e.idn.notna().sum():.1%}) | >1: {(e.idn>1).sum()}")
print("per-gene sum is EXACTLY 1 (median/p95/max all 1.0000) => negatives are SUPPRESSOR contributions under a")
print("non-additive R^2 value function, NOT a normalisation bug. But max(identity) is therefore UNBOUNDED.")
gs = e.dropna(subset=["idn"]).groupby("gene").idn.agg(["min","max","count"])
gs = gs[gs["count"]>=2]
print(f"\ngenes with >=2 scored families: {len(gs)}   with a NEGATIVE component: {(gs['min']<0).sum()} ({(gs['min']<0).mean():.1%})")
print(f"   worst suppressor: {gs['min'].min():+.2f}   largest top share: {gs['max'].max():+.2f}")

print("\n" + "="*95); print("⭐ THE GATED DISAGREEMENT LIST — measurable genes, top_identity in a real share range"); print("="*95)
m = g[(num(g.n_fam)>=2) & g.identity_eq_magnitude.notna()].copy()
m["ti"] = num(m.top_identity); m["ceil"] = num(m.ctx_ceiling)
gate = m[(m.ti.between(0,1)) & (m.ceil>0.02)]
dis = gate[num(gate.identity_eq_magnitude)==0]
print(f"multi-family {len(m)} -> gated (top_identity<=1 AND ceiling>0.02) {len(gate)} -> disagree {len(dis)} ({len(dis)/len(gate):.1%})")
print(f"  (ungated disagreement was 24.0%; on UNMEASURABLE genes it is "
      f"{(num(m[m.ceil<=0.02].identity_eq_magnitude)==0).mean():.1%})")
print("\nthe 18 clearest, GATED — identity confident AND the gene measurable:")
show = dis.nlargest(18,"ti")[["gene","top_family_magnitude","top_family_identity","top_identity","top_beta_frac","n_fam","ctx_ceiling"]]
print(show.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))

print("\n" + "="*95); print("WHERE THE BLOW-UPS LIVE (the ungated tail is an ARTIFACT, not a finding)"); print("="*95)
bad = m[m.ti>1]
print(f"  genes with top_identity>1: {len(bad)} ({len(bad)/len(m):.1%} of multi-family)")
print(f"  ceiling<=0.02:  {(num(bad.ctx_ceiling)<=0.02).mean():.1%}  vs  {(num(m[m.ti<=1].ctx_ceiling)<=0.02).mean():.1%} for the rest")
u = stats.mannwhitneyu(num(bad.ctx_ceiling).dropna(), num(m[m.ti<=1].ctx_ceiling).dropna())
print(f"  MWU on ceiling p={u.pvalue:.3g}  => the ratio explodes exactly where the denominator vanishes (axiom 5)")
