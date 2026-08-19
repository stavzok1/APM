"""Condition the four outcomes on DESIGN WIDTH, which several of them are monotone in by construction."""
import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
from scipy import stats
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False).set_index("gene")
A = pd.read_csv(B/"gene_axes_matrix.tsv", sep="\t").set_index("gene")

def num(s):
    if s.dtype == object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")

nfam = num(g["n_fam"]); narm = num(g["n_arms"])

print("="*78); print("A. HOW MECHANICAL IS EACH OUTCOME? (is it even DEFINED at n_fam==1?)"); print("="*78)
for nm, o in [("identity!=magnitude", 1.0-num(g["identity_eq_magnitude"])),
              ("regulatory_handoff",  num(g["regulatory_handoff"])),
              ("owner_agrees",        num(g["owner_agrees"]))]:
    one = o[nfam <= 1]; multi = o[nfam >= 2]
    print(f"{nm:22s} n_fam==1: n={one.notna().sum():4d} mean={one.mean():.4f}   "
          f"n_fam>=2: n={multi.notna().sum():4d} mean={multi.mean():.4f}")

print()
print("="*78); print("B. PARTIAL ON DESIGN WIDTH — does anything survive?"); print("="*78)
def partial_spearman(y, x, ctrl):
    d = pd.concat([y.rename("y"), x.rename("x")] + [c.rename(f"c{i}") for i,c in enumerate(ctrl)], axis=1).dropna()
    if len(d) < 60: return np.nan, np.nan, len(d)
    R = d.rank()
    C = np.column_stack([np.ones(len(R))] + [R[f"c{i}"].to_numpy() for i in range(len(ctrl))])
    ry = R["y"].to_numpy() - C @ np.linalg.lstsq(C, R["y"].to_numpy(), rcond=None)[0]
    rx = R["x"].to_numpy() - C @ np.linalg.lstsq(C, R["x"].to_numpy(), rcond=None)[0]
    rho, p = stats.spearmanr(ry, rx)
    return rho, p, len(d)

OUT = {
 "identity!=magnitude": 1.0-num(g["identity_eq_magnitude"]),
 "lit_agrees_identity": num(g["lit_agrees_identity"]),
 "DECOY gap_core":      num(g["ctx_gap_core"]),
 "PROT rho_prot":       num(g["cptac_prosp_agg_rho_prot"]),
 "STATE handoff":       num(g["regulatory_handoff"]),
 "STATE owner_agrees":  num(g["owner_agrees"]),
}
CTRL = [nfam, narm]
CIRC = {"card_n_identified","card_frac_identified","ident_n_def","card_top_beta","card_concentration",
        "card_top_beta_frac","card_n_fam","card_n_arms","card_n_regulators","card_n_dense_included",
        "reg_n","card_realized_n_reg","card_greal_n_scored","card_greal_med_coupling",
        "card_greal_frac_c10","card_total_pressure"}
for nm, o in OUT.items():
    rows = []
    for ax in A.columns:
        if ax in CIRC: continue
        rho, p, n = partial_spearman(o, A[ax], CTRL)
        if np.isfinite(rho): rows.append((ax, n, rho, p))
    if not rows: continue
    d = pd.DataFrame(rows, columns=["axis","n","rho_partial","p"]).sort_values("p")
    m = len(d); raw = d["p"].to_numpy() * m / np.arange(1, m+1)
    d["q"] = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1)   # BH step-up: reverse cumulative min
    sig = d[d["q"] < 0.05]
    print(f"\n--- {nm}  (controls: n_fam + n_arms) --- {len(sig)}/{m} survive q<0.05")
    print(d.head(6).to_string(index=False))
