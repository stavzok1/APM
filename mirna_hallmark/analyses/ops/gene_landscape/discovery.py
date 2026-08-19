"""DISCOVERY + the CONVERGENT EVIDENCE LADDER — is the ladder actually monotone?"""
import pandas as pd, numpy as np, pathlib, warnings, json
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
gs = pd.read_csv(B/"discovery_gold_set.tsv", sep="\t", low_memory=False)
e  = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
def num(s):
    if s.dtype==object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")

print("="*92); print("A. THE GOLD SET — what it actually is"); print("="*92)
print(f"rows {len(gs)} | genes {gs.gene.nunique()} | arms {gs.arm.nunique()} | seed families {gs.seed_family.nunique()}")
fam = gs.seed_family.value_counts()
print(f"\n⭐ CONCENTRATION — always quote the families beside the edges:")
for k,v in fam.head(6).items(): print(f"   {v:4d} ({v/len(gs):5.1%})  {k}")
print(f"   top family alone = {fam.iloc[0]/len(gs):.1%} of all 'discoveries'")
print(f"\narms carrying them:"); 
for k,v in gs.arm.value_counts().head(8).items(): print(f"   {v:4d}  {k}")
print(f"\nnovelty guards:  no_he_gene={num(gs.no_he_gene).sum():.0f}/{len(gs)} genes have NO curated regulator at all")
print(f"                 same_seed_he={num(gs.same_seed_he).sum():.0f}/{len(gs)} are same-seed paralogues of an already-curated regulator")
print(f"\nsite/evidence composition:")
for c in ["site_manakov","site_clip_any","has_pred_site"]:
    if c in gs: print(f"   {c:16s} {num(gs[c]).sum():.0f}/{len(gs)} ({num(gs[c]).mean():.1%})")
print(f"   best_type: {gs.best_type.value_counts().head(5).to_dict()}")
print(f"   partial_coupling median {num(gs.partial_coupling).median():+.4f}   null_z median {num(gs.null_z).median():+.3f}")
print(f"   q_bh: min {num(gs.q_bh).min():.3g}  frac<0.05 {(num(gs.q_bh)<0.05).mean():.1%}  <- per-edge FDR is EMPTY by construction")

print("\n" + "="*92); print("B. ⭐ THE CONVERGENT EVIDENCE LADDER — rebuilt on the EDGE CARD (n=5,649)"); print("="*92)
e["cpl"] = num(e["coupling_tum"]); e["site"] = num(e["adm_has_site"]); e["chim"] = num(e["echim_any"])
e["kd"] = num(e["kd_affinity_pct"]); e["nsrc"] = num(e["echim_n_sources"]); e["abund"] = num(e["arm_med_rpm"])
print("\n(1) SITE PRESENT vs ABSENT")
for lab, s in [("has site", e[e.site==1]), ("seedless", e[e.site==0])]:
    c = s.cpl.dropna()
    print(f"   {lab:10s} n={len(c):5d}  median coupling={c.median():+.4f}  frac<-0.10={(c<-0.10).mean():5.1%}")
u = stats.mannwhitneyu(e[e.site==1].cpl.dropna(), e[e.site==0].cpl.dropna()); print(f"   MWU p={u.pvalue:.3g}")
print("\n(2) CHIMERIC DUPLEX — the physical rung")
for lab, s in [("chimeric", e[e.chim==1]), ("no chimeric", e[e.chim==0])]:
    c = s.cpl.dropna()
    print(f"   {lab:12s} n={len(c):5d}  median coupling={c.median():+.4f}  frac<-0.10={(c<-0.10).mean():5.1%}")
u = stats.mannwhitneyu(e[e.chim==1].cpl.dropna(), e[e.chim==0].cpl.dropna()); print(f"   MWU p={u.pvalue:.3g}")
print("\n(3) ⭐ THE LADDER — is coupling MONOTONE in accumulated evidence?")
def rung(r):
    if not r.site: return "0 seedless"
    if r.chim != 1: return "1 site only"
    if r.nsrc == 1: return "2 site+1 chimeric src"
    return "3 site+2 chimeric srcs"
e["rung"] = e.apply(rung, axis=1)
order = ["0 seedless","1 site only","2 site+1 chimeric src","3 site+2 chimeric srcs"]
prev = None
for r in order:
    s = e[e.rung==r].cpl.dropna()
    if len(s)<10: continue
    arrow = "" if prev is None else ("  ↓ stronger" if s.median()<prev else "  ↑ WEAKER — ladder breaks here")
    print(f"   {r:24s} n={len(s):5d}  median={s.median():+.4f}  frac<-0.10={(s<-0.10).mean():5.1%}{arrow}")
    prev = s.median()
ok = e.dropna(subset=["cpl"]).copy(); ok["rn"] = ok.rung.map({v:i for i,v in enumerate(order)})
rho,p = stats.spearmanr(ok.rn, ok.cpl); print(f"\n   Spearman(rung, coupling) = {rho:+.4f}  p={p:.3g}   (negative = ladder climbs)")

print("\n(4) ⚠ IS THE LADDER JUST ABUNDANCE? — abundance-quintile matched")
ok["aq"] = pd.qcut(ok.abund, 5, labels=False, duplicates="drop")
rows=[]
for q in sorted(ok.aq.dropna().unique()):
    s = ok[ok.aq==q]
    a = s[s.rung=="0 seedless"].cpl; b = s[s.rung.isin(["2 site+1 chimeric src","3 site+2 chimeric srcs"])].cpl
    if len(a)>=15 and len(b)>=15:
        rows.append((int(q), len(a), a.median(), len(b), b.median(), stats.mannwhitneyu(a,b).pvalue))
print("   abundance   seedless            site+chimeric        MWU p")
for q,na,ma,nb,mb,pv in rows:
    print(f"   Q{q+1}         n={na:4d} med={ma:+.4f}   n={nb:4d} med={mb:+.4f}   {pv:.3g}")
print(f"   -> gap holds in {sum(1 for r in rows if r[4] < r[2])}/{len(rows)} quintiles")

print("\n(5) K_D AFFINITY — a continuous rung, coupling vs predicted affinity percentile")
s = ok.dropna(subset=["kd"])
rho,p = stats.spearmanr(s.kd, s.cpl)
print(f"   n={len(s)}  Spearman(kd_affinity_pct, coupling) = {rho:+.4f}  p={p:.3g}")
s2 = s[s.site==1]; rho2,p2 = stats.spearmanr(s2.kd, s2.cpl)
print(f"   within SEEDED edges only: n={len(s2)}  rho={rho2:+.4f}  p={p2:.3g}")
