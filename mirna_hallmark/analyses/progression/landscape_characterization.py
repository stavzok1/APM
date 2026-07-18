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


def characterize():
    ec, sc = _load()
    hubs = convergence_hubs(ec); clusters = cluster_units(ec)
    subt = subtype_trajectories(sc); drv, los = driver_loser_rosters(ec)
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
    print(f"\n-> {OUT}/ : landscape_{{convergence_hubs, cluster_units, subtype_trajectories, driver_roster, loser_roster}}.tsv")
    return hubs, clusters, subt, drv, los


if __name__ == "__main__":
    characterize()
