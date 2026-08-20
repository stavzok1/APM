"""Assemble per-gene axes from the cards and scan the four outcomes the user asked about."""
import pandas as pd, numpy as np, pathlib, json, warnings
warnings.filterwarnings("ignore")
from mirna_hallmark.learned import gene_axes as GA

B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False).set_index("gene")
e = pd.read_csv(B/"realization/edge_card.tsv", sep="\t", low_memory=False)
print(f"genes {len(g)} | edges {len(e)}")

# ---- reg_* : the REGULATOR ENSEMBLE, built per gene from the edge card ----
# level  = arm_med_rpm (per-arm median over samples)  == regulator_axes' `med`
# range  = arm_iqr     (per-arm IQR over samples)     -- DEVIATION: gene_axes uses sd; iqr is what the card carries
P = GA._promiscuity_map()
rows = []
for gene, sub in e.groupby("gene"):
    lev = sub["arm_med_rpm"].to_numpy(float)
    rng = sub["arm_iqr"].to_numpy(float)
    arms = sub["arm"].tolist()
    p = np.array([P.get(a, np.nan) for a in arms], dtype=float)
    pf = p[np.isfinite(p)]
    d = {"gene": gene, "reg_n": float(len(sub))}
    if np.isfinite(lev).sum() >= 2:
        d.update({"reg_dose_med": np.nanmedian(lev), "reg_dose_max": np.nanmax(lev),
                  "reg_dose_hhi": GA.hhi(lev[np.isfinite(lev)])})
    if np.isfinite(rng).sum() >= 2:
        fr = rng[np.isfinite(rng)]
        cut = np.nanpercentile(fr, 10) + 1e-9
        d.update({"reg_var_med": np.nanmedian(fr), "reg_var_max": np.nanmax(fr),
                  "reg_var_min": np.nanmin(fr), "reg_var_hhi": GA.hhi(fr),
                  "reg_frac_flat": float(np.mean(fr < cut))})
    d["reg_promisc_cov"] = float(np.mean(np.isfinite(p))) if len(p) else np.nan
    if len(pf) >= 2:
        d.update({"reg_promisc_med": np.median(pf), "reg_promisc_max": np.max(pf),
                  "reg_promisc_min": np.min(pf), "reg_promisc_sd": np.std(pf),
                  "reg_promisc_hhi": GA.hhi(np.expm1(pf))})
    # arm-within-family identifiability (the user's explicit sub-question)
    # ⛔ RIPPLE FIX 2026-08-19 (unit 25): `echim_any` was PRUNED (unit 21 — a literal `True` carrying
    # 0 bits) and `n_regulators`/`n_dense_included` were renamed/dropped. Both call sites here are
    # GUARDED (`if c in sub.columns` / a `[c for c in CARD if ...]` filter), so nothing raised — the
    # axes simply stopped being built. **A guarded reference converts a prune into a silently
    # missing AXIS**, which in a scan reads as 'tested and null'.
    for c, nm in [("cell_arms_resolvable","armres_frac"),("adm_has_site","site_frac"),
                  ("adm_admissible","adm_frac"),("echim_n_sources","chim_frac")]:
        if c in sub.columns:
            v = sub[c]
            v = v.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float) if v.dtype==object else v.astype(float)
            d[nm] = v.mean(skipna=True)
    for c, nm in [("oof_drho","armgain_med"),("arm_sep_z","armsep_med"),
                  ("n_arm_in_cell","narm_cell_max"),("kd_affinity_pct","kd_pct_med"),
                  ("arm_share_of_family_dose","famdose_share_max")]:
        if c in sub.columns:
            s = pd.to_numeric(sub[c], errors="coerce")
            d[nm] = (s.max() if nm.endswith("_max") else s.median())
    rows.append(d)
R = pd.DataFrame(rows).set_index("gene")

# ---- card_* / ident_* / self_* from the gene card ----
CARD = ["n_arms","n_fam","n_identified","frac_identified","concentration",
        "top_beta_frac","total_pressure","top_beta","w_max","ctx_ceiling","ctx_dose_max",
        "ctx_dose_med","ctx_n_abund","ctx_frac_abund","ctx_n_fam_multi","ctx_d_collin",
        "median_retention","realized_retention","realized_n_reg","greal_n_scored",
        # ⛔ `n_dense_included` DROPPED (bit-identical to `n_fam`, already in this list) and
        # `n_regulators` -> `heur_n_regulators` (the §6b-RETIRED lane; `n_arms` is the fit's own width and
        # is already here). Both were silently absent from the axis set until 2026-08-19.
        "greal_med_coupling","greal_frac_c10","n_hallmark_sets", "heur_n_regulators"]
IDENT = ["top_identity","max_beta_frac_sd"]
A = g[[c for c in CARD if c in g.columns]].add_prefix("card_")
A = A.join(g[[c for c in IDENT if c in g.columns]].add_prefix("ident_"))
A["ident_n_def"] = g["top_family_identity"].notna().astype(float)
A = A.join(R)
A = A.apply(pd.to_numeric, errors="coerce")
A = A.loc[:, A.notna().sum() >= 100]
print(f"axes built: {A.shape[1]} axes on {A.shape[0]} genes")
A.to_csv(B/"gene_axes_matrix.tsv", sep="\t")

# ---- outcomes ----
def as_num(s):
    if s.dtype == object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")
OUT = {
 "ATTR identity!=magnitude": 1.0 - as_num(g["identity_eq_magnitude"]),
 "ATTR agrees with literature (identity)": as_num(g["lit_agrees_identity"]),
 "DECOY gap_core (gene)": pd.to_numeric(g["ctx_gap_core"], errors="coerce"),
 "PROTEIN cptac_prosp rho_prot": pd.to_numeric(g["cptac_prosp_agg_rho_prot"], errors="coerce"),
 "PROTEIN prot-beyond-mRNA (disc)": pd.to_numeric(g["cptac_prosp_agg_rho_disc"], errors="coerce"),
 "STATE regulatory_handoff": as_num(g["regulatory_handoff"]),
 "STATE owner_agrees": as_num(g["owner_agrees"]),
 "REALIZED rho_adj": pd.to_numeric(g["realized_rho_adj"], errors="coerce"),
}
CIRC = {"card_n_identified","card_frac_identified","ident_n_def","card_top_beta",
        "card_concentration","card_top_beta_frac"}
res = {}
for name, o in OUT.items():
    o = o.reindex(A.index)
    s = GA.scan(o, A, circular=CIRC)
    res[name] = s
    n = int(o.notna().sum())
    print(f"\n===== {name}   (n={n}) =====")
    sig = s[s["q"] < 0.05] if "q" in s.columns else s.head(0)
    print(f"  axes scanned: {len(s)} | FDR q<0.05: {len(sig)}")
    print(s.head(12).to_string())
json.dump({k: v.to_dict("records") for k, v in res.items()}, open(f"{__import__('os').environ['SP']}/scans.json","w"), default=str)
