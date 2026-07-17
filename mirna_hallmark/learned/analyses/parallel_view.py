"""Parallel view — the design model (learned M) vs the abundance/rank infrastructure.

The abundance stack (`arm_expression` global rank, `mirna_state_class` trajectory, `budget_shift` within-gene
rank, `mirna_comovement` grip) owns **X** (abundance and its cross-state change). The design model adds **M**
(the learned per-edge repression WEIGHT). Realized force = `M·X`. This module puts the two views side by side:

(a) RANK CONCORDANCE (`within_gene_ranks`) — within a gene, rank its regulator families by ABUNDANCE (X̄) vs the
    design ranks: budget (`M·X̄`, force), identity (Shapley coupling ownership), structural (designed owner).
    Spearman concordance + the divergent families (loud passenger / quiet owner). This is a RANK lens.

(b) CHANGE TRAJECTORY (`change_trajectory`) — per regulator across the TRIO GTEx→NAT→tumour: expression
    **logFC** (ΔX magnitude, not just rank) vs wiring change **ΔM**, and the EXACT decomposition of the realized-
    force change `Δ(M·X) = M̄·ΔX (abundance-driven) + X̄·ΔM (wiring-driven)`. This is a MAGNITUDE lens.

Gauge: M fit on z-scored abundance per state (shared gauge, Design §C). Abundance logFC in log2-level, cohort-
median-aligned (GTEx miRNA is TPM vs TCGA RPM) → the trio levels are comparable *relative* abundances; the
GTEx→NAT leg is cross-platform (flagged), NAT→tumour is within-TCGA.

CLI: `python -m mirna_hallmark.learned.analyses.parallel_view PTEN ESR1 ZEB1`
"""
from __future__ import annotations

import contextlib
import io
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import state as ST
from mirna_hallmark.learned import states as STS


# ============================ (a) RANK CONCORDANCE ============================

def within_gene_ranks(gene: str, *, n_perm: int = 400) -> pd.DataFrame:
    """Rank a gene's regulator families by abundance vs the design views (budget / identity / structural),
    with Spearman rank-concordance and the divergent families."""
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import structural_identity as SI
    from mirna_hallmark.learned.attribution import identity_vs_magnitude
    with contextlib.redirect_stdout(io.StringIO()):
        im = identity_vs_magnitude(gene, n_perm=n_perm)      # M, indiv_rho, identity_share, magnitude_share, gap
        st = SI.structural_identity(gene)                    # structural_share (designed owner)
        _, X, _, _ = LD.assemble_gene(gene, w_prior_source="ledger")
    if im.empty:
        print(f"{gene}: identity undefined (<2 nonzero families)"); return pd.DataFrame()
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, _, _ = FAM.collapse_by_family(X, pd.Series(1.0, index=X.columns), fam)
    df = im.rename(columns={"magnitude_share": "budget", "identity_share": "identity"}).copy()
    df["abundance"] = Xf.mean().reindex(df.index).round(2)   # mean log2(RPM+1) per family
    df["structural"] = (st["structural_share"].reindex(df.index).fillna(0.0) if not st.empty else 0.0)
    for src, r in [("abundance", "r_abund"), ("budget", "r_budget"), ("identity", "r_ident"), ("structural", "r_struct")]:
        df[r] = df[src].rank(ascending=False, method="min").astype(int)
    df = df.sort_values("abundance", ascending=False)

    with pd.option_context("display.width", 150):
        print(f"\n=== {gene} — within-gene RANK: abundance vs design views ===")
        print(df[["abundance", "budget", "identity", "structural", "r_abund", "r_budget", "r_ident", "r_struct"]].round(3).to_string())
    # FULL pairwise Spearman-rank concordance matrix (all 6 pairs) — each disagreement names a regulation mode
    keys = ["abundance", "budget", "identity", "structural"]
    cm = pd.DataFrame({a: {b: round(float(spearmanr(df[a], df[b]).correlation), 2) for b in keys} for a in keys})
    print("  pairwise rank concordance:")
    print("   " + cm.to_string().replace("\n", "\n   "))
    meaning = {("budget", "identity"): "force vs coupling-ownership (loud-passenger / quiet-owner gap)",
               ("identity", "structural"): "REALIZED owner vs DESIGNED owner (opportunistic / silenced-specialist)",
               ("abundance", "structural"): "current abundance vs design (anti ⇒ designed owner is silent)",
               ("abundance", "identity"): "does coupling-ownership track loudness"}
    for (a, b), m in meaning.items():
        print(f"   ρ({a},{b})={cm.loc[a,b]:+.2f}  — {m}")
    quiet = df.loc[df["gap"].idxmax()]; loud = df.loc[df["gap"].idxmin()]
    print(f"  ⚑ quiet owner (identity≫budget): {quiet.name} (abund {quiet['abundance']}, budget {quiet['budget']:.2f}, identity {quiet['identity']:.2f})")
    print(f"  ⚑ loud passenger (budget≫identity): {loud.name} (abund {loud['abundance']}, budget {loud['budget']:.2f}, identity {loud['identity']:.2f})")
    # KIND-B collinearity, RESOLVED (not just named): Shapley SPLITS collinear identity credit; the cluster
    # resolver says who'd OWN it under conditioning (COLLINEARITY_AND_IDENTIFIABILITY.md §2). Join the role so a
    # "shared" identity is annotated with its conditioned-coupling owner / co-transcribed rider.
    with contextlib.redirect_stdout(io.StringIO()):
        cr = FAM.cluster_resolution(gene)
    if not cr.empty:
        df["cluster_role"] = cr.set_index("family")["role"].reindex(df.index)
        for ci, sub in cr.groupby("cluster"):
            own = list(sub.loc[sub["role"] == "owner", "family"])
            rid = list(sub.loc[sub["role"] == "rider", "family"])
            print(f"  ⚑ collinear cluster {list(sub['family'])} → owner(conditioned coupling): "
                  f"{own or '— (none survives — genuinely shared)'}" + (f"; rider(no-site): {rid}" if rid else ""))
    return df


# ============================ (b) CHANGE TRAJECTORY (the trio) ============================

def _trio(gene: str, *, n_boot: int = 40):
    """Per regulator-family, across GTEx→NAT→tumour: level abundance + the CANONICAL bagged-NNLS M per state
    (fixed family support → no cross-state selection flips), converted to a level-scale weight M_lvl = M_z/sd
    so `force = M_lvl·X̄` is on one gauge. GTEx abundance quantile-mapped onto the NAT gauge (TPM↔RPM)."""
    from mirna_hallmark.learned import structural_identity as SI
    arms = pd.Index(sorted(SI._functional_regulators(gene)))
    if arms.empty:
        return None
    fam_map = FAM.family_of(arms)
    cols = {}
    for name, st in [("gtex", "gtex"), ("nat", "11"), ("tum", "01")]:
        d = STS._state_family_data(gene, st, fam_map)
        if d is None:
            return None
        Y, Xf, _C = d
        M = STS._bagged_nnls(Y, Xf, _C.to_numpy(float), n_boot=n_boot, sample_type=st)   # z-scale weight, bagged
        sd = Xf.std(ddof=0).clip(lower=0.1)
        cols[f"x_{name}"] = Xf.mean()                                    # level abundance (log2)
        cols[f"M_{name}"] = (M / sd).reindex(Xf.columns)                # level-scale weight M_lvl = M_z/sd
    t = pd.DataFrame(cols).dropna(how="any")
    if len(t) < 3:
        return None
    t["x_gtex_raw"] = t["x_gtex"]                        # naive log2(TPM+1) level (kept for the raw RPM−TPM logFC)
    # GTEx (log2 TPM) → quantile-map onto NAT (log2 RPM): `x_gtex` = QN-into-TCGA gauge (relative change, platform
    # scale + any global shift removed). The naive `x_gtex_raw` is ALSO kept: for miRNAs RPM≈TPM (near-uniform
    # mature length ⇒ length-normalisation ≈ constant; similar count-dominating miRNA sets across libraries), so
    # the raw RPM−TPM difference is a defensible absolute logFC too — reported alongside so both are inspectable.
    ref = np.sort(t["x_nat"].to_numpy(float))
    q = (t["x_gtex"].rank(method="first") - 1) / max(len(t) - 1, 1)
    t["x_gtex"] = np.interp(q, np.linspace(0, 1, len(ref)), ref)
    return t


def change_trajectory(gene: str, *, n_boot: int = 40, top: int = 12) -> pd.DataFrame:
    """ΔlogFC(abundance) vs ΔM(wiring) across GTEx→NAT→tumour + the exact force-change decomposition
    Δ(M·X) = M̄·ΔX (abundance-driven) + X̄·ΔM (wiring-driven), per family."""
    t = _trio(gene, n_boot=n_boot)
    if t is None:
        print(f"{gene}: insufficient trio overlap (GTEx mRNA / NAT / tumour arms)"); return pd.DataFrame()
    # ABUNDANCE change (the robust per-family signal). NAT→tumour = clean within-TCGA logFC; GTEx→NAT reported
    # two ways: relFC_GN = QN-into-TCGA (relative), rawFC_GN = naive RPM−TPM (defensible absolute for miRNAs).
    t["relFC_GN"] = (t["x_nat"] - t["x_gtex"]).round(2)
    t["rawFC_GN"] = (t["x_nat"] - t["x_gtex_raw"]).round(2)
    t["logFC_NT"] = (t["x_tum"] - t["x_nat"]).round(2)
    # WIRING (level-weight M per state). Per-family ΔM tabulated for inspection, but NOT the trusted split —
    # NAT n≈104 vs tumour n≈1000, so ΔM here is partly estimation-power, not biology (§ caveat below). The
    # trusted wiring verdict is the n-matched GTEx→tumour transfer from `state.cross_state_transfer`.
    t["dM_GN"] = (t["M_nat"] - t["M_gtex"]).round(3)
    t["dM_NT"] = (t["M_tum"] - t["M_nat"]).round(3)
    t = t.reindex(t["logFC_NT"].abs().sort_values(ascending=False).index)
    with pd.option_context("display.width", 165):
        print(f"\n=== {gene} — CHANGE trajectory GTEx→NAT→tumour ===")
        print(t[["x_gtex_raw", "x_gtex", "x_nat", "x_tum", "rawFC_GN", "relFC_GN", "logFC_NT", "M_gtex", "M_nat", "M_tum", "dM_GN", "dM_NT"]].head(top).round(3).to_string())
    up = t.nlargest(3, "logFC_NT"); dn = t.nsmallest(3, "logFC_NT")
    print(f"  ABUNDANCE (clean, NAT→tumour logFC): ↑ {', '.join(f'{i}(+{v:.1f})' for i,v in up['logFC_NT'].items())}"
          f"  |  ↓ {', '.join(f'{i}({v:.1f})' for i,v in dn['logFC_NT'].items())}")
    try:                                                                # WIRING verdict — n-matched, aggregate
        cst = ST.cross_state_transfer(gene)
        if "retention" in cst:
            verdict = "STATE-DEPENDENT (wiring rewires healthy→tumour)" if cst.get("state_dependent") else "CONSERVED (only abundance shifts)"
            print(f"  WIRING (n-matched GTEx→tumour, `cross_state_transfer`): retention={cst['retention']} → {verdict}; top ΔM: {cst.get('top_delta','')}")
    except Exception:
        pass
    print("  [relFC_GN cross-platform TPM↔RPM (relative); logFC_NT clean; per-family ΔM is n-limited at NAT — "
          "use the n-matched retention for the wiring verdict]")
    return t


# ============================ (c) STATE PRESSURE ATTRIBUTION (progression axis) ============================

def state_pressure_attribution(gene: str, *, n_boot: int = 40) -> pd.DataFrame:
    """Pinpoint the tumour PRESSURE REALIZER via the progression axis — supply-shift matched to target-response
    shift — for the case where within-tumour variance is uninformative (flat/collinear/power-starved).

    Per regulator family across GTEx→NAT→tumour: realized SUPPLY two ways — `budget_M = X̄·M` (data-force) and
    `budget_pot = X̄·potential` (design-force, potential state-invariant) — the target's own level, and the
    within-state spread. Then the difference-in-differences claim: a family is the tumour pressure realizer iff
    its Δbudget is the CONCENTRATED, DIRECTION-MATCHED supply-shift and the target truly moves (target-CN-adjusted).
    Primary leg NAT→tumour (both TCGA, clean); GTEx (QN'd) is the healthy anchor. REGIME tag (not a hard filter):
      confirmed        — shifts across states AND varies within a state (within-state coupling can corroborate);
      cross-state-only — flat within every state but shifts across (the state contrast is the ONLY signal);
      constitutive     — flat within AND ~no shift across ⇒ unidentifiable by variation (design-only).
    """
    from mirna_hallmark.learned import structural_identity as SI
    arms = pd.Index(sorted(SI._functional_regulators(gene)))
    if arms.empty:
        print(f"{gene}: no functional regulators"); return pd.DataFrame()
    fam_map = FAM.family_of(arms)
    strength, _, _ = SI._potential_matrix("fused")
    per, ylvl, ycn = {}, {}, {}
    for name, st in [("gtex", "gtex"), ("nat", "11"), ("tum", "01")]:
        d = STS._state_family_data(gene, st, fam_map)
        if d is None:
            print(f"{gene}: state {name} unavailable"); return pd.DataFrame()
        Y, Xf, C = d
        M = STS._bagged_nnls(Y, Xf, C.to_numpy(float), n_boot=n_boot, sample_type=st)
        sd = Xf.std(ddof=0)
        per[name] = pd.DataFrame({"xbar": Xf.mean(), "Mlvl": M / sd.clip(lower=0.1), "xstd": sd})
        ylvl[name] = float(Y.mean())
        ycn[name] = float(C["target_cn"].mean()) if "target_cn" in C.columns else np.nan
    fams = sorted(set(per["gtex"].index) & set(per["nat"].index) & set(per["tum"].index))
    if len(fams) < 2:
        print(f"{gene}: <2 families across all three states"); return pd.DataFrame()
    pot = {f: (float(strength.loc[[m for m in fam_map.index[fam_map == f] if m in strength.index], gene].max())
               if gene in strength.columns else 0.0) for f in fams}
    df = pd.DataFrame(index=fams)
    for s in ("gtex", "nat", "tum"):
        df[f"x_{s}"] = per[s].loc[fams, "xbar"]; df[f"xsd_{s}"] = per[s].loc[fams, "xstd"]; df[f"M_{s}"] = per[s].loc[fams, "Mlvl"]
    ref = np.sort(df["x_nat"].to_numpy(float))                            # QN GTEx abundance onto the NAT gauge
    q = (df["x_gtex"].rank(method="first") - 1) / max(len(df) - 1, 1)
    df["x_gtex"] = np.interp(q, np.linspace(0, 1, len(ref)), ref)
    df["potential"] = [round(pot[f], 2) for f in fams]
    # supply (budget) per state — two gauges; NAT→tumour is the clean primary leg
    for s in ("nat", "tum"):
        df[f"budM_{s}"] = (df[f"x_{s}"] * df[f"M_{s}"]).round(3)
        df[f"budP_{s}"] = (df[f"x_{s}"] * df["potential"]).round(2)
    df["dBudM"] = (df["budM_tum"] - df["budM_nat"]).round(3)              # supply shift (data-force)
    df["dBudP"] = (df["budP_tum"] - df["budP_nat"]).round(2)             # supply shift (design-force)
    dY = round(ylvl["tum"] - ylvl["nat"], 2)                             # target move (NAT→tumour)
    dYcn = round(ycn["tum"] - 2.0, 2) if ycn["tum"] == ycn["tum"] else np.nan  # target CN vs diploid (NAT≈2): its own-CN confound
    # concentration = share of the total |ΔbudP| supply-shift; direction-match = supply↑ ⇒ target↓
    tot = df["dBudP"].abs().sum() or 1.0
    df["concentration"] = (df["dBudP"].abs() / tot).round(2)
    df["dir_match"] = np.sign(df["dBudP"]) != np.sign(dY)                 # opposite signs = supply explains the move
    # regime: does the family shift across (|ΔbudP| notable) and/or vary within a state
    within_var = df[["xsd_gtex", "xsd_nat", "xsd_tum"]].max(axis=1)
    shifts = df["dBudP"].abs() >= 0.1 * tot
    df["regime"] = np.where(~shifts, "constitutive",
                            np.where(within_var >= 0.3, "confirmed-able", "cross-state-only"))
    df = df.sort_values("concentration", ascending=False)
    show = ["potential", "x_nat", "x_tum", "budP_nat", "budP_tum", "dBudP", "dBudM", "concentration", "dir_match", "regime"]
    with pd.option_context("display.width", 165):
        print(f"\n=== {gene} — STATE pressure attribution (NAT→tumour; target ΔY={dY:+.2f}, target ΔCN={dYcn:+.2f}) ===")
        print(df[show].head(10).to_string())
    matched = df[df["dir_match"]]
    if len(matched) and dY == dY:
        top = matched["dBudP"].abs().idxmax(); r = df.loc[top]; uniform = 1.0 / len(df)
        conc = "CONCENTRATED" if r["concentration"] >= 1.5 * uniform else "diffuse (shift spread across regulators)"
        cnflag = (" ⚠ target move aligns with its OWN CN change — confounded" if
                  (dYcn == dYcn and np.sign(dYcn) == np.sign(dY) and abs(dYcn) > 0.2) else "")
        print(f"  ⚑ tumour pressure-realizer nominee: {top} (ΔbudP {r['dBudP']:+.2f} = {int(r['concentration']*100)}% of "
              f"the supply-shift vs {int(uniform*100)}% uniform → {conc}; regime={r['regime']}){cnflag}")
    else:
        print(f"  no direction-matched supply-shift (target ΔY={dY:+.2f})")
    return df


if __name__ == "__main__":
    genes = [a for a in sys.argv[1:] if not a.startswith("-")] or ["PTEN", "ESR1", "ZEB1"]
    for g in genes:
        within_gene_ranks(g)
        change_trajectory(g)
        state_pressure_attribution(g)
