import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
e = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
def num(s):
    if s.dtype==object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")

print("="*95); print("1. ATTRIBUTION — where LMG and magnitude NAME DIFFERENT REGULATORS"); print("="*95)
m = g[(num(g.n_fam)>=2) & g.identity_eq_magnitude.notna()].copy()
dis = m[num(m.identity_eq_magnitude)==0]
print(f"multi-family genes {len(m)} | disagree {len(dis)} ({len(dis)/len(m):.1%})")
dis = dis.assign(ti=num(dis.top_identity), tb=num(dis.top_beta_frac))
show = dis.nlargest(18,"ti")[["gene","top_family_magnitude","top_family_identity","top_identity","top_beta_frac","n_fam","ctx_ceiling"]]
print("\n⭐ The clearest disagreements (identity most confident that magnitude names the WRONG family):")
print(show.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
print("\n⭐ Genes where identity is UNDEFINED (NNLS zeroed EVERY family — beta cannot say this):")
und = g[g.top_family_identity.isna() & num(g.n_fam).ge(2)]
print(f"   {len(und)} multi-family genes. Examples: {', '.join(und.gene.head(14))}")
print(f"   their median ceiling {num(und.ctx_ceiling).median():+.4f} vs defined {num(g[g.top_family_identity.notna()].ctx_ceiling).median():+.4f}")
print("\n⭐ Genes where identity AGREES WITH THE LITERATURE but magnitude does NOT (identity's win cases):")
w = g[(num(g.lit_agrees_identity)==1) & (num(g.lit_agrees_magnitude)==0)]
print(f"   n={len(w)}: {', '.join(w.gene.head(20))}")
l = g[(num(g.lit_agrees_identity)==0) & (num(g.lit_agrees_magnitude)==1)]
print(f"   reverse (magnitude right, identity wrong) n={len(l)}: {', '.join(l.gene.head(20))}")
print(f"   => identity net {len(w)-len(l):+d} genes over magnitude against the literature")

print("\n" + "="*95); print("2. ARM-IN-FAMILY — the genes where resolving the arm matters MOST"); print("="*95)
mm = e[num(e.n_arm_in_cell)>1].copy(); mm["dr"] = num(mm.oof_drho)
top = mm.nsmallest(15,"dr")[["gene","arm","seed_family","oof_drho","cell_arm_sep_z","cell_arms_resolvable","n_arm_in_cell","coupling_tum"]]
print(top.to_string(index=False))

print("\n" + "="*95); print("3. PROTEIN — the strongest protein-coupled genes, and prospective-vs-t105 agreement"); print("="*95)
p = g.dropna(subset=["cptac_prosp_agg_rho_prot"]).copy()
p["pr"] = num(p.cptac_prosp_agg_rho_prot); p["t1"] = num(p.cptac_t105_agg_rho_prot); p["rna"] = num(p.tcga_agg_rho_rna)
print("⭐ top 15 by prospective protein coupling (the INDEPENDENT cohort):")
print(p.nsmallest(15,"pr")[["gene","cptac_prosp_agg_rho_prot","cptac_t105_agg_rho_prot","tcga_agg_rho_rna",
                            "cptac_prosp_agg_ret_prot","top_family_magnitude"]].to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
both = p.dropna(subset=["t1"])
r,pv = stats.spearmanr(both.pr, both.t1)
print(f"\ncross-cohort agreement prospective vs t105: n={len(both)} rho={r:+.4f} p={pv:.3g}")
r2,p2 = stats.spearmanr(both.pr, both.rna.reindex(both.index))
print(f"protein(prosp) vs mRNA(TCGA):               n={both.rna.notna().sum()} rho={r2:+.4f} p={p2:.3g}")
rep = both[(both.pr<-0.10)&(both.t1<-0.10)]
print(f"⭐ genes protein-coupled BELOW -0.10 in BOTH cohorts: n={len(rep)} -> {', '.join(rep.gene.head(18))}")

print("\n" + "="*95); print("4. STATE — the genes miR-21 takes over, and the biggest handoffs"); print("="*95)
h = g[(num(g.regulatory_handoff)==1) & (g.dominant_HLY!=g.dominant_TUM)]
t21 = h[h.dominant_TUM=="hsa-miR-21-5p"]
print(f"genes where miR-21-5p becomes dominant in tumour: {len(t21)}")
print(f"  they were dominated in HEALTHY by: {t21.dominant_HLY.value_counts().head(8).to_dict()}")
print(f"  named: {', '.join(t21.gene.head(24))}")
print(f"\n  their repression class: {t21.heur_repression_class.value_counts().to_dict()}")
print(f"  vs all handoff genes:   {h.heur_repression_class.value_counts().to_dict()}")
lost = h[h.dominant_HLY=="hsa-miR-26a-5p"]
print(f"\ngenes miR-26a-5p LOSES: {len(lost)} -> {', '.join(lost.gene.head(18))}")
print(f"  handed to: {lost.dominant_TUM.value_counts().head(6).to_dict()}")
