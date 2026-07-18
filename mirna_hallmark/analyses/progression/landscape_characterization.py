"""BIOLOGICAL CHARACTERIZATION of the progression/state landscape — the descriptive map the axis never got.

NOT an estimand test: this DESCRIBES the biology already encoded in the annotations + cards (shift_class, acquired
dose, realization, roles, host/genomic context, seed family, PAM50 subtype, arm trajectory). Four threads:
  1. CONVERGENCE HUBS   — the TSGs that acquire the most repressive edges (PTEN pile-up), coordination.
  2. CLUSTER-AS-UNIT    — do genomic clusters (shared host) / seed families co-acquire / co-lose as units?
  3. SUBTYPE trajectory — which arms/edges acquire subtype-specifically (Basal vs Luminal), from PAM50 + tau.
  4. DRIVER/LOSER roster — named, dose-ranked acquired oncomiRs (→ TSG) and silenced TSG-miRs (→ oncogene).

⚠ DESCRIPTIVE, on annotations with known caveats: shift_class per-edge FDR null is 3–4× too narrow
(LANDSCAPE_REPORT banner debt); role annotations incomplete (many 'unknown'). Read the PATTERNS (cluster
co-movement, TSG-miR loss, convergence hubs), NOT per-edge significance. Reuses built outputs; cheap.
⚠ CROSS-PLATFORM (GTEx TPM → TCGA RPM): TRUST the RANK axis (dHT/grank percentile deltas) + Shapley over the QN
magnitude bridge (log2fc/arm_lfc_HLY_TUM_QN) — QN is correct+settled but a softer assumption. The healthy→tumour
'acquired' calls lead with rank; QN-magnitude-only gainers are the lower-confidence layer. The NAT→tumour paired
dose/realization (①–④) is TCGA same-platform ⇒ NOT a QN concern.

  python -m mirna_hallmark.analyses.progression.landscape_characterization
"""
from __future__ import annotations

import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
from mirna_hallmark import config as C

OUT = Path(C.OUTPUT_ROOT) / "learned" / "realization"
_SC = Path(C.OUTPUT_ROOT) / "tissue_reference" / "mirna_state_class"
_GX = Path(C.OUTPUT_ROOT) / "learned" / "genomic_context.tsv"
_GAINED = ("acquired_realized", "tumour_realized", "field_established_realized")   # gained-repression classes


def _load():
    ec = pd.read_csv(OUT / "progression_edge_card.tsv", sep="\t")
    gx = pd.read_csv(_GX, sep="\t")[["arm", "mir_class", "host", "clustered", "promoter_gene"]]
    ec = ec.merge(gx, on="arm", how="left")
    sc = pd.read_csv(_SC / "mirna_state_class.tsv", sep="\t").rename(columns={"miRNA": "arm"})
    return ec, sc


# --------------------------------------------------------------------------- #
def convergence_hubs(ec: pd.DataFrame) -> pd.DataFrame:
    """TSG/oncogene targets ranked by how many ACQUIRED repressive edges converge on them — the pile-up map."""
    g = ec[ec.shift_class.isin(_GAINED)].copy()
    rows = []
    for gene, s in g.groupby("gene"):
        arms = s.arm.nunique(); fams = s.seed_family.nunique()
        rows.append({"gene": gene, "role": s.role.iloc[0], "n_acquired_edges": arms, "n_families": fams,
                     "coordinated": fams >= 3,                                  # many independent families converge
                     "mean_dose": round(float(s.mean_own_shift.mean()), 2),
                     "mean_realization_rho": round(float(s.edge_rho_adj.mean()), 3),
                     "gene_net_repressed": bool(s.gene_net_repr.iloc[0]) if s.gene_net_repr.notna().any() else None,
                     "top_arms": ", ".join(s.sort_values("mean_own_shift", ascending=False).arm.str.replace("hsa-", "", regex=False).head(4))})
    df = pd.DataFrame(rows).sort_values("n_acquired_edges", ascending=False)
    df.to_csv(OUT / "landscape_convergence_hubs.tsv", sep="\t", index=False)
    return df


def cluster_units(ec: pd.DataFrame) -> pd.DataFrame:
    """Do genomic clusters move as UNITS? Two co-transcription proxies: shared HOST gene and shared SEED FAMILY.
    Per unit (≥2 landscape arms): mean acquired dose, dose concordance (all same sign), dominant shift_class."""
    arm = ec.dropna(subset=["mean_own_shift"]).drop_duplicates("arm")
    rows = []
    for key, col in (("host", "host"), ("seed_family", "seed_family")):
        for unit, s in arm.dropna(subset=[col]).groupby(col):
            if s.arm.nunique() < 2 or unit in ("", "nan"):
                continue
            dose = s.mean_own_shift
            same_dir = float(np.sign(dose).replace(0, np.nan).dropna().pipe(lambda x: (x == x.mode().iloc[0]).mean())) if dose.notna().any() else np.nan
            sc = ec[ec.arm.isin(s.arm)].shift_class.dropna()
            rows.append({"unit_type": key, "unit": unit, "n_arms": s.arm.nunique(),
                         "mean_dose": round(float(dose.mean()), 2), "dose_sd": round(float(dose.std()), 2),
                         "dose_concordance": round(same_dir, 2) if same_dir == same_dir else np.nan,
                         "direction": "acquire" if dose.mean() > 0.2 else ("lose" if dose.mean() < -0.2 else "mixed"),
                         "dominant_shift_class": sc.mode().iloc[0] if len(sc) else None,
                         "members": ", ".join(s.sort_values("mean_own_shift", ascending=False).arm.str.replace("hsa-", "", regex=False).head(6))})
    df = pd.DataFrame(rows).sort_values(["unit_type", "n_arms"], ascending=[True, False])
    df.to_csv(OUT / "landscape_cluster_units.tsv", sep="\t", index=False)
    return df


def subtype_trajectories(sc: pd.DataFrame) -> pd.DataFrame:
    """Subtype-specific ACQUISITION: arms that gain (acquired_gainer) and are subtype-specific (high subtype_tau),
    by dominant subtype — who acquires in Basal vs Luminal vs Her2."""
    g = sc[sc.acquired_gainer == True].copy() if "acquired_gainer" in sc else sc.copy()
    keep = [c for c in ["arm", "dHT", "log2fc_tumor_vs_healthy", "subtype_tau", "dominant_subtype",
                        "secondary_class", "acquired_axis"] if c in g.columns]
    g = g[keep].copy()
    g["arm"] = g["arm"].str.replace("hsa-", "", regex=False)
    g = g.sort_values("subtype_tau", ascending=False)
    g.to_csv(OUT / "landscape_subtype_trajectories.tsv", sep="\t", index=False)
    return g


def driver_loser_rosters(ec: pd.DataFrame):
    """Named catalogs: (drivers) acquired oncomiRs — dose>0, target=TSG, realized (ρ<0); (losers) silenced
    TSG-miRs — dose<0 / lost, target=oncogene."""
    drv = ec[(ec.mean_own_shift > 0.5) & (ec.role == "tsg") & (ec.edge_rho_adj < 0)].copy()
    drv = drv.sort_values("mean_own_shift", ascending=False)[
        ["arm", "gene", "mean_own_shift", "edge_rho_adj", "coupling_tum", "seed_family", "host", "own_specific_frac"]]
    drv["arm"] = drv["arm"].str.replace("hsa-", "", regex=False)
    drv.to_csv(OUT / "landscape_driver_roster.tsv", sep="\t", index=False)
    los = ec[(ec.shift_class == "lost") & (ec.role == "oncogene")].copy()
    los = los.sort_values("mean_own_shift")[["arm", "gene", "mean_own_shift", "coupling_nat", "coupling_tum", "seed_family", "host"]]
    los["arm"] = los["arm"].str.replace("hsa-", "", regex=False)
    los.to_csv(OUT / "landscape_loser_roster.tsv", sep="\t", index=False)
    return drv, los


# --------------------------------------------------------------------------- #
# WIDENING threads (2026-07-18) — more descriptive dimensions across the landscape.
def trajectory_archetypes(sc: pd.DataFrame) -> pd.DataFrame:
    """Census of per-arm TRAJECTORY archetypes (`primary_class` = progressive/tumor_acquired/field_effect × gain/loss,
    + `acquired_axis` = rank vs magnitude) with named exemplars — the miRNA-level trajectory taxonomy."""
    rows = []
    for cls, s in sc.groupby("primary_class"):
        ex = s.reindex(s.get("dHT", pd.Series(0, index=s.index)).abs().sort_values(ascending=False).index)
        rows.append({"primary_class": cls, "n_arms": len(s),
                     "frac_acquired_gainer": round(float((s.get("acquired_gainer") == True).mean()), 2) if "acquired_gainer" in s else np.nan,
                     "exemplars": ", ".join(ex["arm"].str.replace("hsa-", "", regex=False).head(5)) if "arm" in ex else ""})
    df = pd.DataFrame(rows).sort_values("n_arms", ascending=False)
    if "acquired_axis" in sc:
        df.attrs["axis_census"] = sc[sc.get("acquired_gainer") == True]["acquired_axis"].value_counts().to_dict()
    df.to_csv(OUT / "landscape_trajectory_archetypes.tsv", sep="\t", index=False)
    return df


def regulatory_handoffs(gc: pd.DataFrame) -> pd.DataFrame:
    """Genes whose DOMINANT regulator SWITCHES across states (healthy→tumour) — named handoffs (who takes over)."""
    h = gc[gc.regulatory_handoff == True].dropna(subset=["dominant_HLY", "dominant_TUM"]).copy()
    for c in ["dominant_HLY", "dominant_NAT", "dominant_TUM"]:
        h[c] = h[c].str.replace("hsa-", "", regex=False)
    h["handoff"] = h["dominant_HLY"] + " → " + h["dominant_TUM"]
    keep = [c for c in ["gene", "dominant_HLY", "dominant_NAT", "dominant_TUM", "realized_rho_adj",
                        "acquired_vs_nat", "gene_repression_class"] if c in h.columns]
    out = h[keep].sort_values("acquired_vs_nat", ascending=False)
    out.to_csv(OUT / "landscape_regulatory_handoffs.tsv", sep="\t", index=False)
    return out


def dose_realization_quadrants(gc: pd.DataFrame) -> pd.DataFrame:
    """The FUNCTIONAL 2×2 — gene acquired dose (acquired_vs_nat) × realized (realized_rho_adj<−0.1): DRIVER
    (acquired+realized) · BUFFERED (acquired+unrealized) · PRE-SET (not-acquired+realized) · INERT."""
    g = gc.dropna(subset=["acquired_vs_nat", "realized_rho_adj"]).copy()
    acq = g.acquired_vs_nat > 0; real = g.realized_rho_adj < -0.1
    g["quadrant"] = np.select([acq & real, acq & ~real, ~acq & real], ["DRIVER", "BUFFERED", "PRE-SET"], "INERT")
    cen = g.groupby("quadrant").agg(n_genes=("gene", "size"), mean_dose=("acquired_vs_nat", "mean"),
                                    mean_realization=("realized_rho_adj", "mean"),
                                    frac_net_repressed=("gene_net_repressed_tumor", "mean")).round(3)
    g[["gene", "quadrant", "acquired_vs_nat", "realized_rho_adj", "n_regulators", "gene_net_repressed_tumor"]].to_csv(
        OUT / "landscape_dose_realization_quadrants.tsv", sep="\t", index=False)
    cen.attrs["drivers_top"] = ", ".join(g[g.quadrant == "DRIVER"].sort_values("realized_rho_adj").gene.head(10))
    cen.attrs["buffered_top"] = ", ".join(g[g.quadrant == "BUFFERED"].sort_values("acquired_vs_nat", ascending=False).gene.head(10))
    return cen


def subtype_driver_rosters(sc: pd.DataFrame) -> pd.DataFrame:
    """Per-PAM50-subtype ACQUIRED drivers — the top subtype-specific gainers per dominant_subtype."""
    if "dominant_subtype" not in sc:
        return pd.DataFrame()
    g = sc[(sc.get("acquired_gainer") == True) & (sc.subtype_tau > 0.5)].copy()
    g["arm"] = g["arm"].str.replace("hsa-", "", regex=False)
    rows = []
    for st, s in g.groupby("dominant_subtype"):
        top = s.sort_values("subtype_tau", ascending=False).head(6)
        rows.append({"subtype": st, "n_specific_gainers": len(s),
                     "top_arms": ", ".join(top["arm"] + "(" + top["subtype_tau"].round(2).astype(str) + ")")})
    df = pd.DataFrame(rows).sort_values("n_specific_gainers", ascending=False)
    df.to_csv(OUT / "landscape_subtype_driver_rosters.tsv", sep="\t", index=False)
    return df


def genomic_context_census(ec: pd.DataFrame) -> pd.DataFrame:
    """Trajectory by GENOMIC CONTEXT — host-coupled (sense coding/lncRNA host) vs intergenic vs antisense, and the
    imprinted 14q32 locus. Uses the paired NAT→tumour dose (`mean_own_shift`, same-platform) — trusted, not QN."""
    arm = ec.dropna(subset=["mean_own_shift"]).drop_duplicates("arm")
    rows = []
    for cls, s in arm.dropna(subset=["mir_class"]).groupby("mir_class"):
        d = s.mean_own_shift
        sc_ = ec[ec.arm.isin(s.arm)].shift_class.dropna()
        rows.append({"mir_class": cls, "n_arms": s.arm.nunique(), "mean_dose": round(float(d.mean()), 2),
                     "frac_acquire": round(float((d > 0.2).mean()), 2), "frac_lose": round(float((d < -0.2).mean()), 2),
                     "dominant_shift_class": sc_.mode().iloc[0] if len(sc_) else None})
    # imprinted 14q32 / DLK1-DIO3 (MEG3/8/9, MIR493HG hosts)
    imp = arm[arm.host.astype(str).str.contains("MEG|MIR493HG|DLK1|DIO3", na=False)]
    if len(imp):
        rows.append({"mir_class": "IMPRINTED_14q32", "n_arms": imp.arm.nunique(),
                     "mean_dose": round(float(imp.mean_own_shift.mean()), 2),
                     "frac_acquire": round(float((imp.mean_own_shift > 0.2).mean()), 2),
                     "frac_lose": round(float((imp.mean_own_shift < -0.2).mean()), 2),
                     "dominant_shift_class": ec[ec.arm.isin(imp.arm)].shift_class.dropna().mode().iloc[0]
                     if ec[ec.arm.isin(imp.arm)].shift_class.notna().any() else None})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "landscape_genomic_context_census.tsv", sep="\t", index=False)
    return df


def program_acquired_dose(gc: pd.DataFrame) -> pd.DataFrame:
    """Which Hallmark PROGRAMS' genes acquire the most repressive pressure — mean gene `acquired_vs_nat` (paired
    NAT→tumour, SAME-PLATFORM, trusted) + realization. (Realization is NOT program-specific per MH-163; this is the
    DOSE description, distinct — which programs the acquired pressure LANDS on.)"""
    from mirna_hallmark.hallmark_sets import gene_to_hallmarks
    g2h = gene_to_hallmarks()
    ap = gc.set_index("gene")["acquired_vs_nat"] if "acquired_vs_nat" in gc else pd.Series(dtype=float)
    rr = gc.set_index("gene")["realized_rho_adj"] if "realized_rho_adj" in gc else pd.Series(dtype=float)
    rows = []
    for prog in sorted(set(h for hs in g2h.values() for h in hs)):
        genes = [g for g, hs in g2h.items() if prog in hs]
        a = ap.reindex(genes).dropna()
        if len(a) < 15:
            continue
        rows.append({"hallmark": prog, "n_genes": len(a), "mean_acquired_dose_NAT_TUM": round(float(a.mean()), 3),
                     "frac_gained": round(float((a > 0).mean()), 2),
                     "mean_realization": round(float(rr.reindex(genes).dropna().mean()), 3)})
    df = pd.DataFrame(rows).sort_values("mean_acquired_dose_NAT_TUM", ascending=False)
    df.to_csv(OUT / "landscape_program_acquired_dose.tsv", sep="\t", index=False)
    return df


def trajectory_two_step(sc: pd.DataFrame) -> pd.DataFrame:
    """The FULL 3-state trajectory via the HEALTHY ANCHOR, RANK-based (trusted): decompose each arm's healthy→tumour
    move into the FIELD step dHN (healthy→NAT, GTEx-anchored) and the MALIGNANT step dNT (NAT→tumour). Where does
    the change happen — in the pre-malignant field or the malignant transition? Uses rank deltas, NOT QN magnitude."""
    d = sc.dropna(subset=["dHN", "dNT"]).copy()
    THR = 0.1
    field, malig = d.dHN.abs() > THR, d.dNT.abs() > THR
    d["change_locus"] = np.select([field & ~malig, ~field & malig, field & malig],
                                  ["field_only", "malignant_only", "both"], "stable")
    from scipy.stats import spearmanr
    cont = spearmanr(d.dHN, d.dNT, nan_policy="omit").correlation      # does the field step CONTINUE into malignancy?
    rows = []
    for loc, s in d.groupby("change_locus"):
        ex = s.reindex(s.dHT.abs().sort_values(ascending=False).index)
        rows.append({"change_locus": loc, "n_arms": len(s), "mean_dHN_field": round(float(s.dHN.mean()), 3),
                     "mean_dNT_malignant": round(float(s.dNT.mean()), 3),
                     "exemplars": ", ".join(ex["arm"].str.replace("hsa-", "", regex=False).head(5))})
    df = pd.DataFrame(rows).sort_values("n_arms", ascending=False)
    df.attrs["field_continues_into_malignant_rho"] = round(float(cont), 3) if cont == cont else np.nan
    df.attrs["field_dominates"] = float((d.dHN.abs() > d.dNT.abs()).mean())
    df.to_csv(OUT / "landscape_trajectory_two_step.tsv", sep="\t", index=False)
    return df


def healthy_anchored_drivers(sc: pd.DataFrame) -> pd.DataFrame:
    """RANK-supported (dHT>0.15, NOT QN-magnitude-only) healthy→tumour acquired gainers, split by whether they are
    acquired DE NOVO (silent in healthy, rank_gtex low) vs AMPLIFIED from an already-present healthy baseline."""
    g = sc[(sc.get("dHT", pd.Series(0, index=sc.index)) > 0.15)].copy()
    if "rank_gtex" in g:
        g["healthy_status"] = np.where(g.rank_gtex < 0.5, "acquired_de_novo", "amplified_from_healthy")
    g["arm"] = g["arm"].str.replace("hsa-", "", regex=False)
    keep = [c for c in ["arm", "dHN", "dNT", "dHT", "rank_gtex", "rank_tumor", "healthy_status",
                        "primary_class", "dominant_subtype"] if c in g.columns]
    out = g[keep].sort_values("dHT", ascending=False)
    out.to_csv(OUT / "landscape_healthy_anchored_drivers.tsv", sep="\t", index=False)
    return out


def cohort_view(ec: pd.DataFrame) -> pd.DataFrame:
    """COHORT contrast (all ~1000 tumours vs ~104 NAT, UNPAIRED) via `arm_lfc_NAT_TUM` — better-powered than the 103
    pairs, no within-patient control. (1) validates the paired view generalises (paired `mean_own_shift` vs cohort
    LFC agreement); (2) cohort driver/loser roster; (3) `coupling_tum` (~1000-tumour) is the well-powered cohort
    coupling vs the n=103 paired `edge_rho_adj`."""
    arm = ec.dropna(subset=["arm_lfc_NAT_TUM", "mean_own_shift"]).drop_duplicates("arm")
    rho = spearmanr(arm.mean_own_shift, arm.arm_lfc_NAT_TUM).correlation
    sign_agree = float((np.sign(arm.mean_own_shift) == np.sign(arm.arm_lfc_NAT_TUM)).mean())
    div = arm[(np.sign(arm.mean_own_shift) != np.sign(arm.arm_lfc_NAT_TUM)) & (arm.mean_own_shift.abs() > 1)]
    ac = arm.copy(); ac["arm"] = ac["arm"].str.replace("hsa-", "", regex=False)
    out = ac[["arm", "mean_own_shift", "arm_lfc_NAT_TUM", "coupling_nat", "coupling_tum"]].sort_values(
        "arm_lfc_NAT_TUM", ascending=False)
    out.to_csv(OUT / "landscape_cohort_view.tsv", sep="\t", index=False)
    out.attrs.update({"n_arms": len(arm), "paired_vs_cohort_rho": round(float(rho), 3),
                      "sign_agreement": round(sign_agree, 2), "n_big_divergences": int(len(div))})
    return out


_CN = Path(C.OUTPUT_ROOT) / "mirna_locus_cnv"


def _norm(a):
    return str(a).replace("hsa-", "").lower()


def direct_gtex_tumor(sc: pd.DataFrame) -> pd.DataFrame:
    """DIRECT healthy(GTEx)→tumour comparison via dHT rank (bypasses the field-contaminated NAT intermediate — the
    'what's different in cancer vs TRULY healthy' question). Is the direct change field-mediated (tracks dHN) or
    malignant (tracks dNT)? Named direct gainers/losers."""
    d = sc.dropna(subset=["dHT"]).copy()
    d["arm"] = d["arm"].str.replace("hsa-", "", regex=False)
    has = d.dropna(subset=["dHN", "dNT"])
    cf = spearmanr(has.dHT, has.dHN, nan_policy="omit").correlation
    cm = spearmanr(has.dHT, has.dNT, nan_policy="omit").correlation
    top = d.sort_values("dHT", ascending=False)
    out = top[[c for c in ["arm", "dHT", "dHN", "dNT", "primary_class", "dominant_subtype"] if c in top.columns]]
    out.to_csv(OUT / "landscape_direct_gtex_tumor.tsv", sep="\t", index=False)
    out.attrs.update({"dHT_tracks_field_dHN_rho": round(float(cf), 3), "dHT_tracks_malignant_dNT_rho": round(float(cm), 3),
                      "n_direct_gain": int((d.dHT > 0.15).sum()), "n_direct_lose": int((d.dHT < -0.15).sum())})
    return out


def cn_layer(ec: pd.DataFrame) -> pd.DataFrame:
    """Is ACQUISITION copy-number-driven? Join per-arm cohort tumour CN + CN↔expr concordance with acquired dose
    (`arm_lfc_NAT_TUM`, cohort). Correlate; flag CN-driven acquirers (e.g. miR-17~92 at 13q31)."""
    cn = pd.read_csv(_CN / "mirna_cnv_by_stratum.tsv", sep="\t")
    cn = cn[(cn.entity_level == "arm") & (cn.stratum == "cohort_all")].copy()
    cn["k"] = cn.entity_label.map(_norm); cn_med = cn.groupby("k")["median_copy_number"].mean()
    conc = pd.read_csv(_CN / "mirna_cnv_expr_concordance.tsv", sep="\t")
    conc["k"] = conc.arm.map(_norm); cn_drv = conc.set_index("k")["partial_rho"]
    arm = ec.dropna(subset=["arm_lfc_NAT_TUM"]).drop_duplicates("arm").copy()
    arm["k"] = arm.arm.map(_norm)
    arm["tumour_cn"] = arm.k.map(cn_med); arm["cn_expr_partial_rho"] = arm.k.map(cn_drv)
    arm["arm2"] = arm.arm.str.replace("hsa-", "", regex=False)
    m = arm.dropna(subset=["tumour_cn"])
    rho_dose_cn = spearmanr(m.arm_lfc_NAT_TUM, m.tumour_cn, nan_policy="omit").correlation
    out = m[["arm2", "arm_lfc_NAT_TUM", "tumour_cn", "cn_expr_partial_rho"]].sort_values("tumour_cn", ascending=False)
    out.to_csv(OUT / "landscape_cn_layer.tsv", sep="\t", index=False)
    out.attrs.update({"n_arms": len(m), "rho_acquired_dose_vs_tumourCN": round(float(rho_dose_cn), 3),
                      "median_cn_expr_partial_rho": round(float(conc.partial_rho.median()), 3)})
    return out


def ago_layer(ec: pd.DataFrame) -> pd.DataFrame:
    """Is acquired pressure REALIZED where RISC/AGO loading is present? Stratify edge realization (`edge_rho_adj`)
    by `ago_loading` tertile, among ACQUIRED-dose edges (mean_own_shift>0)."""
    e = ec.dropna(subset=["edge_rho_adj", "ago_loading", "mean_own_shift"]).copy()
    acq = e[e.mean_own_shift > 0].copy()
    acq["ago_tertile"] = pd.cut(acq.ago_loading, [-0.01, 0.33, 0.67, 1.01], labels=["low", "mid", "high"])
    out = acq.groupby("ago_tertile").agg(n_edges=("edge_rho_adj", "size"),
                                         mean_realization=("edge_rho_adj", "mean"),
                                         frac_realized=("edge_rho_adj", lambda s: float((s < 0).mean()))).round(3)
    out.to_csv(OUT / "landscape_ago_layer.tsv", sep="\t")
    rho = spearmanr(acq.ago_loading, acq.edge_rho_adj, nan_policy="omit").correlation
    out.attrs["rho_agoloading_vs_realization"] = round(float(rho), 3)
    return out


def characterize():
    ec, sc = _load()
    gc = pd.read_csv(OUT / "progression_gene_card.tsv", sep="\t")
    ctx = genomic_context_census(ec); pad = program_acquired_dose(gc)
    twostep = trajectory_two_step(sc); hdrv = healthy_anchored_drivers(sc)
    coh = cohort_view(ec); dgt = direct_gtex_tumor(sc); cnl = cn_layer(ec); agol = ago_layer(ec)
    hubs = convergence_hubs(ec); clusters = cluster_units(ec)
    subt = subtype_trajectories(sc); drv, los = driver_loser_rosters(ec)
    arch = trajectory_archetypes(sc); hand = regulatory_handoffs(gc)
    quad = dose_realization_quadrants(gc); subr = subtype_driver_rosters(sc)
    print("=" * 100)
    print("PROGRESSION LANDSCAPE — BIOLOGICAL CHARACTERIZATION (descriptive; annotation caveats apply)")
    print("=" * 100)
    print("\n[1] CONVERGENCE HUBS — TSGs/oncogenes by # acquired repressive edges (coordinated = ≥3 seed families):")
    print(hubs[hubs.role.isin(["tsg", "oncogene", "oncogene/tsg"])].head(14).to_string(index=False))
    print("\n[2] CLUSTER-AS-UNIT — co-transcribed (host) / seed-family units that move together (top by n_arms):")
    cu = clusters[(clusters.n_arms >= 3) & (clusters.direction != "mixed")].sort_values("n_arms", ascending=False)
    print(cu.head(14).to_string(index=False))
    print("\n[3] SUBTYPE-SPECIFIC ACQUIRERS — subtype_tau-ranked acquired gainers:")
    print(subt.head(12).to_string(index=False))
    print("\n[4a] DRIVER ROSTER — acquired oncomiRs repressing TSGs (dose-ranked):")
    print(drv.head(14).to_string(index=False))
    print("\n[4b] LOSER ROSTER — silenced TSG-miRs de-repressing oncogenes:")
    print(los.head(12).to_string(index=False))
    print("\n[5] TRAJECTORY ARCHETYPES — per-arm trajectory taxonomy (primary_class):")
    print(arch.to_string(index=False))
    print(f"    acquired-gainer axis census (rank vs magnitude): {arch.attrs.get('axis_census')}")
    print("\n[6] REGULATORY HANDOFFS — genes whose dominant regulator SWITCHES healthy→tumour (top by acquired pressure):")
    print(hand.head(14).to_string(index=False))
    print(f"    ({len(hand)} genes switch dominant regulator across states)")
    print("\n[7] FUNCTIONAL 2×2 — acquired dose × realized:")
    print(quad.to_string())
    print(f"    DRIVER (acquired+realized): {quad.attrs.get('drivers_top')}")
    print(f"    BUFFERED (acquired+unrealized): {quad.attrs.get('buffered_top')}")
    print("\n[8] PER-SUBTYPE ACQUIRED DRIVERS:")
    print(subr.to_string(index=False))
    print("\n[9] GENOMIC-CONTEXT trajectory census (paired NAT→tumour dose, same-platform):")
    print(ctx.to_string(index=False))
    print("\n[10] PROGRAM ACQUIRED DOSE — which Hallmark programs' genes acquire the most repressive pressure (top/bottom):")
    print(pd.concat([pad.head(6), pad.tail(4)]).to_string(index=False))
    print("\n[11] TWO-STEP TRAJECTORY (healthy anchor, RANK-based) — field step dHN vs malignant step dNT:")
    print(twostep.to_string(index=False))
    print(f"    field CONTINUES into malignant (ρ dHN,dNT)={twostep.attrs.get('field_continues_into_malignant_rho')} · "
          f"field-step-dominates in {100*twostep.attrs.get('field_dominates',0):.0f}% of arms")
    print("\n[12] HEALTHY-ANCHORED acquired gainers (rank dHT>0.15, trusted) — de novo vs amplified:")
    if "healthy_status" in hdrv:
        print("    " + str(hdrv.healthy_status.value_counts().to_dict()))
    print(hdrv.head(12).to_string(index=False))
    print(f"\n[13] COHORT VIEW (~1000 tumour vs ~104 NAT, unpaired) vs PAIRED — {coh.attrs}")
    print("     ⇒ paired generalises to the cohort (validates ①–④); coupling_tum is the well-powered cohort coupling.")
    print(f"\n[14] DIRECT GTEx→TUMOUR (rank dHT, bypasses NAT) — {dgt.attrs}")
    print(dgt.head(6).to_string(index=False))
    print(f"\n[15] CN LAYER — is acquisition copy-number-driven? {cnl.attrs}")
    print(cnl.head(6).to_string(index=False))
    print(f"\n[16] AGO LAYER — is acquired pressure realized where RISC loading is present? {agol.attrs}")
    print(agol.to_string())
    print(f"\n-> {OUT}/ : landscape_{{convergence_hubs, cluster_units, subtype_trajectories, driver_roster, loser_roster, "
          f"trajectory_archetypes, regulatory_handoffs, dose_realization_quadrants, subtype_driver_rosters}}.tsv")
    return hubs, clusters, subt, drv, los, arch, hand, quad, subr


if __name__ == "__main__":
    characterize()
