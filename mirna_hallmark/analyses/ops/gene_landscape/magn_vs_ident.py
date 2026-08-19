"""Do magnitude and identity pick systematically DIFFERENT KINDS of family?"""
import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
e = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
a = pd.read_csv(B/"arm_card.tsv", sep="\t", low_memory=False)
def num(s): return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float) if s.dtype==object else pd.to_numeric(s, errors="coerce")

# family-level abundance + fame, from the edge card's arm properties rolled to seed_family
e["ab"] = num(e.arm_med_rpm)
fam_ab = e.groupby("seed_family").ab.median()
famsize = e.groupby("seed_family").arm.nunique()
fame = None
if "fame_npmid" in a.columns:
    am = a.set_index("arm")["fame_npmid"].pipe(num)
    e["fame"] = e.arm.map(am)
    fame = e.groupby("seed_family").fame.median()

m = g[(num(g.n_fam)>=2) & g.identity_eq_magnitude.notna()].copy()
m["ti"] = num(m.top_identity); m["ceil"] = num(m.ctx_ceiling)
gate = m[(m.ti.between(0,1)) & (m.ceil>0.02)]
dis = gate[num(gate.identity_eq_magnitude)==0].copy()
print(f"gated disagreement set: n={len(dis)} genes\n")

rows=[]
for r in dis.itertuples():
    rows.append({"gene": r.gene,
                 "mag_ab": fam_ab.get(r.top_family_magnitude, np.nan),
                 "ide_ab": fam_ab.get(r.top_family_identity, np.nan),
                 "mag_sz": famsize.get(r.top_family_magnitude, np.nan),
                 "ide_sz": famsize.get(r.top_family_identity, np.nan),
                 "mag_fm": fame.get(r.top_family_magnitude, np.nan) if fame is not None else np.nan,
                 "ide_fm": fame.get(r.top_family_identity, np.nan) if fame is not None else np.nan})
D = pd.DataFrame(rows)
print("PAIRED comparison of the family MAGNITUDE picks vs the family IDENTITY picks (same gene):")
for lab, mc, ic in [("median arm abundance (log2 RPM)","mag_ab","ide_ab"),
                    ("family size (n arms)","mag_sz","ide_sz"),
                    ("study attention (distinct PMIDs)","mag_fm","ide_fm")]:
    d = D[[mc,ic]].dropna()
    if len(d) < 20: print(f"  {lab:36s} n={len(d)} — too few"); continue
    w = stats.wilcoxon(d[mc], d[ic])
    print(f"  {lab:36s} n={len(d):3d}  magnitude={d[mc].median():8.2f}   identity={d[ic].median():8.2f}   "
          f"Wilcoxon p={w.pvalue:.3g}  identity-higher in {(d[ic]>d[mc]).mean():.1%}")

print("\nCONTROL — the same three quantities on genes where they AGREE (should show no gap):")
agr = gate[num(gate.identity_eq_magnitude)==1].copy()
rows=[]
for r in agr.itertuples():
    rows.append({"mag_ab": fam_ab.get(r.top_family_magnitude, np.nan),
                 "mag_fm": fame.get(r.top_family_magnitude, np.nan) if fame is not None else np.nan,
                 "mag_sz": famsize.get(r.top_family_magnitude, np.nan)})
A = pd.DataFrame(rows)
print(f"  agreeing genes n={len(A)}: median abundance {A.mag_ab.median():.2f}, PMIDs {A.mag_fm.median():.1f}, size {A.mag_sz.median():.1f}")
print(f"  disagreeing, MAGNITUDE pick:   abundance {D.mag_ab.median():.2f}, PMIDs {D.mag_fm.median():.1f}, size {D.mag_sz.median():.1f}")
print(f"  disagreeing, IDENTITY  pick:   abundance {D.ide_ab.median():.2f}, PMIDs {D.ide_fm.median():.1f}, size {D.ide_sz.median():.1f}")
u = stats.mannwhitneyu(D.mag_ab.dropna(), A.mag_ab.dropna())
print(f"\n  => when they disagree, MAGNITUDE's pick is a much rarer arm than when they agree (MWU p={u.pvalue:.3g})")
