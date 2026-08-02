"""RUNG PROVENANCE — for every card column, WHICH UNIT was it computed on?

    .venv/bin/python3 -m mirna_hallmark.learned.card_rungs [--check]

WHY (user-asked, 2026-08-01). `edge_card.tsv` is keyed `(gene, arm)`, so **every row is an ARM** — but the
VALUES on that row come from four different units, and nothing on the card says which:

    beta, identity        estimated on the SEED FAMILY, then BROADCAST to each member arm
    coupling_tum          a single-ARM partial Spearman
    ctx_ceiling, comp_*   a property of the GENE, broadcast
    ctx_arm_dose          a property of the ARM alone

That silent mixing is not hypothetical — it caused three defects in one day: `cptac_card` applied a
FAMILY-estimated β to RAW ARM abundance (MH-179); the card asserted a FAMILY `beta` beside an ARM
`coupling_tum` with nothing marking it (MH-187); and a within-cell Shapley was compared against a
GENE-level OOF statistic (MH-188). **Every one was a unit mismatch, i.e. axiom 6.**

⚠ **AND THE NAMING IS NOW WORSE THAN UNIFORM SILENCE.** The columns added on 2026-08-01 carry their rung
(`coupling_fam`, `beta_arm`, `arm_credit_share`) while the canonical ones do not (`beta` is family,
`coupling_tum` is arm) — so a reader who trusts the naming convention will infer the OPPOSITE of the truth
for the oldest and most-used columns. This module fixes that by RECORD rather than by renaming 160
columns: renaming would break every consumer and every historical ledger/registry row.

THE FOUR RUNGS
    key      the join keys
    gene     a property of the gene, broadcast to all its rows
    family   computed with the SEED FAMILY as the unit, broadcast to member arms
    arm      a property of the arm alone
    edge     computed on the (gene, arm) PAIR with the arm as the unit

⚠ `family` vs `edge` is the distinction that keeps biting: both look like "per-edge" on the card, but a
`family` value is SHARED by every same-seed arm of that gene while an `edge` value is not. If you group by
arm and average a `family` column, you are averaging duplicates.
"""
from __future__ import annotations

import argparse

import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
CARD = OUT / "realization/edge_card.tsv"
DEST = OUT / "edge_card_rungs.tsv"

# explicit column -> rung. Prefix rules follow; explicit wins.
EXPLICIT = {
    "gene": "key", "arm": "key", "seed_family": "key",
    # ── FAMILY rung: estimated on the seed family, BROADCAST to member arms
    **{c: "family" for c in (
        "beta", "beta_sd", "z", "identified", "identity", "identity_deconv", "identity_eq_magnitude",
        "beta_frac", "beta_frac_sd", "m_nnls", "pip_dense", "pip_discovery", "prior_pi",
        "beta_deconv", "retention_beta", "composition_class", "net_pressure", "family_size",
        "coupling_fam", "share_TUM", "rank_TUM", "share_HLY", "share_NAT", "rank_HLY", "rank_NAT",
        "identity_allocated", "family_rho_adj", "realization_identity", "is_realization_owner")},
    # ── ARM rung: a property of the arm alone
    **{c: "arm" for c in (
        "arm_med_rpm", "arm_pct_floor", "arm_iqr", "arm_id_status", "ago_loading", "detection",
        "arm_lfc_HLY_NAT_raw", "arm_lfc_HLY_TUM_QN", "arm_lfc_HLY_TUM_raw", "arm_lfc_NAT_TUM",
        "ctx_arm_dose", "ctx_arm_abundant", "spiker", "grank_HLY", "grank_NAT", "grank_TUM",
        "mean_own_shift", "mean_dGlobalRank", "own_specific_frac",
        "beta_arm", "sd_arm", "z_arm")},
    # ── EDGE rung: the (gene, arm) pair, arm as the unit
    **{c: "edge" for c in (
        "coupling_hly", "coupling_nat", "coupling_tum", "retention_rho", "shift_class",
        "edge_rho_adj", "dShare_M_own", "dShare_raw_own", "realization_score", "dose_class",
        "coupling_p_tum", "coupling_p_nat", "coupling_z_tum", "wiring_frac", "acquisition_score")},
    # ── residual sweep (2026-08-01): every column that was UNASSIGNED on first run, classified by
    # what UNIT actually produced it. An unassigned column is where a unit mismatch hides, so this
    # list is kept exhaustive and `--check` fails the build if a new one appears.
    **{c: "gene" for c in ("n", "p_fam", "role", "family_role")},
    **{c: "family" for c in (
        "beta_sd_deconv", "beta_frac_deconv", "beta_frac_abs", "beta_frac_reliable",
        "retention_reliable", "family_dose_share",
        "d_rank_HLY_NAT", "d_rank_NAT_TUM", "d_rank_HLY_TUM")},
    **{c: "arm" for c in (
        "dose_rank_HLY", "dose_rank_NAT", "dose_rank_TUM",
        "dGlobal_HLY_NAT", "dGlobal_NAT_TUM", "dGlobal_HLY_TUM",
        "healthy_leg", "healthy_potential", "healthy_uninformative",
        "surrogate_instrument", "surrogate_corr", "coupling_hly_surrogate",
        "dose_comp_retention", "dose_prolif_retention", "dose_confounded")},
    **{c: "edge" for c in (
        "coupling_p_hly", "coupling_p_hly_resolved", "coupling_p_hly_surrogate",
        "coupling_z_hly", "coupling_z_nat",
        "term_ABUND", "term_WIRING", "term_INTERACT")},
    # ── WITHIN-FAMILY: an arm's share of its family — neither purely arm nor family
    **{c: "arm-in-family" for c in (
        "arm_credit_share", "arm_dbeta", "arm_sep_z", "arm_resolvable", "n_arm_in_cell",
        "oof_rho_arm", "oof_rho_fam", "oof_drho")},
}
PREFIX = (("ctx_", "gene"), ("comp_", "gene"), ("cptac_", "edge"), ("tcga_", "gene"),
          ("cal_", "family"), ("gene_", "gene"))   # cal_* derive from beta => family rung


def rung_of(col: str) -> str:
    if col in EXPLICIT:
        return EXPLICIT[col]
    for p, r in PREFIX:
        if col.startswith(p):
            return r
    return "UNASSIGNED"


def build() -> pd.DataFrame:
    card = pd.read_csv(CARD, sep="\t", low_memory=False, nrows=5)
    prov = {}
    pf = OUT / "edge_card_base_provenance.tsv"
    if pf.exists():
        p = pd.read_csv(pf, sep="\t")
        prov = dict(zip(p["column"], p.get("estimator", pd.Series(dtype=str))))
    rows = [{"column": c, "rung": rung_of(c), "estimator": prov.get(c, "")} for c in card.columns]
    R = pd.DataFrame(rows)
    R.to_csv(DEST, sep="\t", index=False)
    return R


def report(R: pd.DataFrame) -> None:
    print(f"=== RUNG PROVENANCE — {len(R)} columns on edge_card.tsv ===\n")
    for r, s in R.groupby("rung"):
        print(f"  {r:14s} {len(s):3d}")
    print()
    fam = R[R.rung == "family"]["column"].tolist()
    print("⚠ FAMILY-rung columns — SHARED by every same-seed arm of the gene. Averaging over arms")
    print("  double-counts them; a groupby(arm) on these is a groupby over duplicates:")
    print("   ", ", ".join(fam[:14]), "..." if len(fam) > 14 else "")
    bad = R[R.rung == "UNASSIGNED"]
    print(f"\n{'⛔' if len(bad) else '✅'} UNASSIGNED: {len(bad)}")
    if len(bad):
        print("   ", ", ".join(bad["column"].tolist()[:20]))
        print("    add them to EXPLICIT/PREFIX — an unassigned column is exactly how a unit mismatch hides.")
    print("\n⚠ NAMING IS NOT A GUIDE: `beta` is FAMILY, `coupling_tum` is EDGE, and neither says so.")
    print("  The columns added 2026-08-01 DO carry their rung, which makes the older ones look")
    print("  arm-level by contrast. Read THIS table, not the column name.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="exit non-zero if any column is UNASSIGNED")
    a = ap.parse_args()
    R = build()
    report(R)
    if a.check and (R.rung == "UNASSIGNED").any():
        raise SystemExit(1)
