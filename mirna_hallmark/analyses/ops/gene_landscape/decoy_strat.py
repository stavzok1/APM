"""DECOY at the gene rung: the stratification and the TAILS, with named genes."""
import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
def num(s):
    if s.dtype==object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")
for c in ["ctx_gap_core","ctx_gap_deconv","ctx_ceiling","n_fam","n_arms","ctx_gap_retention",
          "median_retention","total_pressure","ctx_dose_max","greal_med_coupling","top_beta"]:
    g[c] = num(g[c])
g["measurable"] = num(g["ctx_measurable"])
d = g.dropna(subset=["ctx_gap_core"]).copy()
print(f"genes with a gap: {len(d)}   median {d.ctx_gap_core.median():+.4f}   mean {d.ctx_gap_core.mean():+.4f}")
print(f"  frac gap<0 (real beats fake): {(d.ctx_gap_core<0).mean():.1%}")

print("\n" + "="*90); print("STRATUM 1 — a-priori competence class"); print("="*90)
for k, s in d.groupby("ctx_apriori_class"):
    t = stats.wilcoxon(s.ctx_gap_core) if len(s)>10 else None
    print(f"  {k:14s} n={len(s):4d}  median={s.ctx_gap_core.median():+.4f}  mean={s.ctx_gap_core.mean():+.4f}"
          f"  frac<0={(s.ctx_gap_core<0).mean():5.1%}  p={t.pvalue:.2e}" if t else f"  {k} n={len(s)}")

print("\n" + "="*90); print("STRATUM 2 — design width (the degeneracy split)"); print("="*90)
d["wbin"] = pd.cut(d.n_fam, [0,1,2,3,5,99], labels=["1 (β INERT)","2","3","4-5","6+"])
for k, s in d.groupby("wbin"):
    if len(s)<10: continue
    t = stats.wilcoxon(s.ctx_gap_core)
    print(f"  n_fam={str(k):12s} n={len(s):4d}  median={s.ctx_gap_core.median():+.4f}  frac<0={(s.ctx_gap_core<0).mean():5.1%}  p={t.pvalue:.2e}")

print("\n" + "="*90); print("STRATUM 3 — measurability (the sharpest separator; partly definitional)"); print("="*90)
d["cbin"] = pd.cut(d.ctx_ceiling, [-9,0,0.02,0.05,0.15,9], labels=["<=0","0-0.02","0.02-0.05","0.05-0.15",">0.15"])
for k, s in d.groupby("cbin"):
    if len(s)<10: continue
    t = stats.wilcoxon(s.ctx_gap_core)
    print(f"  ceiling {str(k):11s} n={len(s):4d}  median={s.ctx_gap_core.median():+.4f}  frac<0={(s.ctx_gap_core<0).mean():5.1%}  p={t.pvalue:.2e}")

print("\n" + "="*90); print("⭐ THE TAILS — the 20 genes where curation beats the fake MOST"); print("="*90)
cols = ["gene","ctx_gap_core","ctx_gap_deconv","ctx_ceiling","n_fam","n_arms","top_family_magnitude",
        "median_retention","heur_repression_class"]
lo = d.nsmallest(20, "ctx_gap_core")[cols]
print(lo.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
print("\n" + "="*90); print("⛔ THE OTHER TAIL — the 20 genes where the FAKE beats curation"); print("="*90)
hi = d.nlargest(20, "ctx_gap_core")[cols]
print(hi.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

print("\n" + "="*90); print("WHAT SEPARATES THE TAILS (top-100 vs bottom-100, MWU)"); print("="*90)
A = d.nsmallest(100,"ctx_gap_core"); Z = d.nlargest(100,"ctx_gap_core")
for c in ["ctx_ceiling","n_fam","n_arms","median_retention","total_pressure","ctx_dose_max",
          "greal_med_coupling","top_beta","measurable"]:
    a, z = A[c].dropna(), Z[c].dropna()
    if len(a)<20 or len(z)<20: continue
    u = stats.mannwhitneyu(a, z)
    print(f"  {c:22s} real-favoured med={a.median():+.4f}   fake-favoured med={z.median():+.4f}   p={u.pvalue:.2e}")
