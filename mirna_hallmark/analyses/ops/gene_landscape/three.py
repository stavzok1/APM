import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
e = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
def num(s):
    if s.dtype==object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")

print("="*78); print("(a) ARM vs FAMILY — where does resolving arms WITHIN a family buy anything?"); print("="*78)
sub = e.copy()
sub["n_arm_in_cell"] = num(sub["n_arm_in_cell"]); sub["oof_drho"]=num(sub["oof_drho"])
sub["arm_resolvable"]=num(sub["arm_resolvable"]); sub["arm_sep_z"]=num(sub["arm_sep_z"])
sub["oof_rho_arm"]=num(sub["oof_rho_arm"]); sub["oof_rho_fam"]=num(sub["oof_rho_fam"])
multi = sub[sub["n_arm_in_cell"]>1]
print(f"edges total {len(sub)} | in MULTI-ARM family cells {len(multi)} ({len(multi)/len(sub):.1%}) "
      f"-- outside these, arm==family BY CONSTRUCTION")
d = multi["oof_drho"].dropna()
print(f"\noof_drho (arm-level OOF rho MINUS family-level), multi-arm cells: n={len(d)}")
print(f"  median {d.median():+.4f} | mean {d.mean():+.4f} | IQR [{d.quantile(.25):+.4f},{d.quantile(.75):+.4f}]")
print(f"  frac ARM BETTER (<0 = more negative rho = better repression fit): {(d<0).mean():.1%}")
w = stats.wilcoxon(d) if len(d)>10 else None
print(f"  Wilcoxon vs 0: stat={w.statistic:.0f} p={w.pvalue:.3g}" if w else "")
print(f"  ** the TAIL (never quote the median): p05={d.quantile(.05):+.4f}  p95={d.quantile(.95):+.4f} "
      f" max|gain|={d.abs().max():.4f}")
r = multi["arm_resolvable"].dropna()
print(f"\narm_resolvable in multi-arm cells: {r.mean():.1%} of {len(r)} edges")
for lo,hi,lab in [(1,1,"1 arm"),(2,2,"2 arms"),(3,3,"3 arms"),(4,99,"4+ arms")]:
    m = sub[(sub["n_arm_in_cell"]>=lo)&(sub["n_arm_in_cell"]<=hi)]
    dd = m["oof_drho"].dropna(); rr = m["arm_resolvable"].dropna()
    print(f"  {lab:8s} n={len(m):5d}  oof_drho med={dd.median() if len(dd) else float('nan'):+.4f}"
          f"  resolvable={rr.mean() if len(rr) else float('nan'):.1%}")

print()
print("="*78); print("(b) PROTEIN — what CPTAC cohorts exist and how much is actually populated"); print("="*78)
for pre,label in [("cptac_prosp","PROSPECTIVE (independent patients, TMT-11)"),
                  ("cptac_t105","TCGA-105 (same patients as TCGA => layer-only)")]:
    ne = e[f"{pre}_n"].notna().sum() if f"{pre}_n" in e else 0
    ng = g[f"{pre}_n_arms"].notna().sum() if f"{pre}_n_arms" in g else 0
    nmed = num(e[f"{pre}_n"]).median() if f"{pre}_n" in e else np.nan
    print(f"\n{label}")
    print(f"  edge card: {ne} edges carry it | gene card: {ng} genes | median patients per edge n={nmed:.0f}")
    for c in [f"{pre}_rho_rna", f"{pre}_rho_prot", f"{pre}_rho_disc", f"{pre}_ret_prot"]:
        if c in e.columns:
            s = num(e[c]).dropna()
            print(f"    {c:26s} n={len(s):5d}  med={s.median():+.4f}  frac<0={(s<0).mean():.1%}")
print("\nRPPA (MH-203, n=866 same patients):", "output/learned/rppa/", 
      [p.name for p in (B/"rppa").iterdir()] if (B/"rppa").exists() else "MISSING")
print("Is RPPA on any card?", [c for c in list(e.columns)+list(g.columns) if "rppa" in c.lower()] or "NO -- not carded")

print()
print("="*78); print("(c) STATE HANDOFF — who hands off to whom"); print("="*78)
h = g[num(g["regulatory_handoff"])==1]
print(f"genes with a regulatory handoff: {len(h)} / {len(g)} ({len(h)/len(g):.1%})")
pair = h.groupby(["dominant_HLY","dominant_TUM"]).size().sort_values(ascending=False)
print("\ntop HEALTHY -> TUMOUR dominant-regulator transitions:")
print(pair.head(12).to_string())
print("\narms GAINING dominance (tumour) vs LOSING it (healthy), among handoff genes:")
gain = h["dominant_TUM"].value_counts().head(8)
lose = h["dominant_HLY"].value_counts().head(8)
net = (h["dominant_TUM"].value_counts() - h["dominant_HLY"].value_counts()).dropna().sort_values()
print("  top TUMOUR-dominant:", dict(gain))
print("  top HEALTHY-dominant:", dict(lose))
print("  biggest NET GAINERS:", dict(net.tail(6)))
print("  biggest NET LOSERS :", dict(net.head(6)))
print("\ndominant_edge_shift_class among handoff genes:")
print(h["dominant_edge_shift_class"].value_counts().head(8).to_string())
