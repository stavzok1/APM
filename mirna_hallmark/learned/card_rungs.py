"""CANONICAL CARD REGISTRY — for every (card, column): which UNIT is it defined on, and does the data agree?

    .venv/bin/python3 -m mirna_hallmark.learned.card_rungs            # full report
    .venv/bin/python3 -m mirna_hallmark.learned.card_rungs --check    # non-zero exit on any gap

WHY (user-asked, 2026-08-01). Three cards exist at three grains and a column's meaning depends on WHICH
card it is on. Nothing recorded that, and the resulting unit mismatches produced four defects in one day
— MH-179 (a FAMILY-estimated beta applied to RAW ARM abundance), MH-187 (a family weight beside an arm
correlation, unmarked), MH-188 (a within-cell Shapley compared against a GENE-level OOF statistic), and
MH-191 (beta labelled `family` when `readouts.run(level="arm")` fits it PER ARM). All four are the same
error: **not knowing what unit a number lives on.**

THE TRAP THIS REGISTRY EXISTS FOR — `beta` MEANS A DIFFERENT THING ON EACH CARD:
      edge_card.tsv    beta is EDGE rung   — readouts.run(level="arm"),    the SS8 collapse REMOVED
      family_card.tsv  beta is FAMILY rung — readouts.run(level="family"), the SS8 collapse APPLIED
    Same name, same estimator, DIFFERENT unit. => THE RUNG IS A PROPERTY OF (CARD, COLUMN), NEVER OF THE
    COLUMN NAME. A single shared prefix map cannot express that, which is why each card gets its own.

TWO ORTHOGONAL FACTS PER COLUMN — do not conflate them:
    rung    the unit the value is DEFINED on: key / gene / family / arm / edge / arm-in-family
    agg_of  if the value SUMMARISES a lower unit, which one (else "")
            e.g. on the gene card `cptac_prosp_agg_rho_prot` is rung=gene, agg_of=arm — one row per gene,
            but computed by summing beta*X over the gene's ARMS, so it inherits the arm rung's caveats.

THE LABELS ARE TESTED, NOT ASSERTED. Each rung implies an exact invariance GIVEN THE CARD'S KEY, and
`verify()` checks it against the data:
    on a (gene, arm)-keyed card:   gene => constant within gene · family => constant within
                                   (gene, seed_family) · arm => constant within arm ACROSS genes ·
                                   edge => free
    on a (gene, family)-keyed card: gene => constant within gene · family => free (it IS the grain)
    on a (gene)-keyed card:        one row per gene => nothing is checkable; `agg_of` carries the meaning
A declared rung the data contradicts is a mislabel — that is how MH-191's wrong labels and the
`dose_rank_*` / `family_role` errors were caught.
"""
from __future__ import annotations

import argparse

import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"

# the readouts attribution block. Its RUNG depends on the card: the same estimator runs at two levels.
_READOUTS = ("beta", "beta_sd", "z", "identified", "identity", "identity_deconv",
             "identity_eq_magnitude", "beta_frac", "beta_frac_sd", "beta_frac_abs",
             "beta_frac_reliable", "beta_sd_deconv", "beta_frac_deconv", "retention_reliable",
             "m_nnls", "pip_dense", "pip_discovery", "prior_pi", "beta_deconv", "retention_beta",
             "composition_class", "net_pressure", "identity_allocated")
_GENE_DESC = ("n", "p_fam", "role", "n_arms", "n_fam")
_ARM_DESC = ("arm_med_rpm", "arm_pct_floor", "arm_iqr", "arm_id_status", "ago_loading", "detection",
             "spiker", "healthy_leg", "healthy_potential", "healthy_uninformative",
             "arm_id_status_UNUSED", "surrogate_instrument", "surrogate_corr",
             "dose_comp_retention", "dose_prolif_retention", "dose_confounded",
             "arm_lfc_HLY_NAT_raw", "arm_lfc_HLY_TUM_QN", "arm_lfc_HLY_TUM_raw", "arm_lfc_NAT_TUM",
             "grank_HLY", "grank_NAT", "grank_TUM", "dGlobal_HLY_NAT", "dGlobal_NAT_TUM",
             "dGlobal_HLY_TUM", "mean_own_shift", "mean_dGlobalRank", "own_specific_frac",
             "ctx_arm_dose", "ctx_arm_abundant")
_EDGE_ONLY = ("coupling_hly", "coupling_nat", "coupling_tum", "retention_rho", "shift_class",
              "edge_rho_adj", "dShare_M_own", "dShare_raw_own", "realization_score", "dose_class",
              "coupling_p_tum", "coupling_p_nat", "coupling_z_tum", "coupling_z_hly", "coupling_z_nat",
              "coupling_p_hly", "coupling_p_hly_resolved", "coupling_p_hly_surrogate",
              "wiring_frac", "acquisition_score", "term_ABUND", "term_WIRING", "term_INTERACT",
              "share_TUM", "share_HLY", "share_NAT", "rank_TUM", "rank_HLY", "rank_NAT",
              "d_rank_HLY_NAT", "d_rank_NAT_TUM", "d_rank_HLY_TUM",
              "dose_rank_HLY", "dose_rank_NAT", "dose_rank_TUM",
              "realization_identity", "is_realization_owner",
              # caught empirically: these are fit/computed per (gene, arm), not per arm alone
              "beta_arm", "sd_arm", "z_arm", "coupling_hly_surrogate")
_FAMILY_ATTR = ("family_size", "coupling_fam", "family_rho_adj", "seed_family", "arms")
_ARM_IN_FAM = ("arm_credit_share", "arm_dbeta", "arm_sep_z", "arm_resolvable", "n_arm_in_cell",
               "oof_rho_arm", "oof_rho_fam", "oof_drho",
               # caught empirically: these vary INSIDE a (gene, family) cell, so they describe
               # an ARM WITHIN ITS FAMILY, not the family as a whole (MH-166 allocation work)
               "arm_id_status", "family_dose_share", "family_role")
_GENE_ROLLUP = ("gene_dominated", "gene_net_repr", "gene_nreg", "gene_repr_class",
                "gene_lfc_NAT_TUM", "gene_iqr", "gene_med_rpm", "gene_pct_floor")

CARDS = {
    "edge": dict(
        path=OUT / "realization/edge_card.tsv", key=["gene", "arm"],
        explicit={**{c: "edge" for c in _READOUTS + _EDGE_ONLY},
                  **{c: "arm" for c in _ARM_DESC},
                  **{c: "family" for c in _FAMILY_ATTR},
                  **{c: "arm-in-family" for c in _ARM_IN_FAM},
                  **{c: "gene" for c in _GENE_DESC + _GENE_ROLLUP},
                  "gene": "key", "arm": "key"},
        prefix=(("ctx_", "gene"), ("comp_", "gene"), ("cptac_", "edge"), ("tcga_", "gene"),
                ("cal_", "edge"))),
    "family": dict(
        path=OUT / "family_card.tsv", key=["gene", "family"],
        # the SAME readouts columns are FAMILY rung here — level="family", the collapse applied
        explicit={**{c: "family" for c in _READOUTS + _FAMILY_ATTR},
                  **{c: "gene" for c in _GENE_DESC + _GENE_ROLLUP},
                  "retention": "family",     # beta_deconv/beta_core at the family grain
                  "n_arms": "family",        # arms IN THIS FAMILY (caught empirically: varies within gene)
                  "w_max": "family",         # max curated evidence weight over the family's arms
                  "gene": "key", "family": "key"},
        prefix=(("ctx_", "gene"), ("comp_", "gene"), ("cal_", "family"))),
    "gene": dict(
        path=OUT / "realization/gene_card.tsv", key=["gene"],
        explicit={"gene": "key"},
        prefix=(("", "gene"),)),          # one row per gene => everything is gene rung by construction
    # ⭐ THE FOURTH RUNG (MH-209). Until 2026-08-03 there was NO arm card, so every per-arm property
    # (promiscuity, K_D, host locus, AGO loading, seed rarity) lived scattered in its own module — which
    # is how MH-208's promiscuity annotation stayed mis-specified for five weeks with one consumer.
    "arm": dict(
        path=OUT / "arm_card.tsv", key=["arm"],
        explicit={"arm": "key"},
        prefix=(("", "arm"),)),           # one row per arm => everything is arm rung by construction
}
AGG_OF = {"gene": {**{c: "arm" for c in ("cptac_prosp_agg_rho_rna", "cptac_prosp_agg_rho_prot",
                                         "cptac_prosp_agg_rho_disc", "cptac_t105_agg_rho_rna",
                                         "cptac_t105_agg_rho_prot", "cptac_t105_agg_rho_disc",
                                         "tcga_agg_rho_rna", "cptac_prosp_abund_rho_rna",
                                         "cptac_prosp_abund_rho_prot", "cptac_prosp_abund_rho_disc",
                                         "tcga_abund_rho_rna", "n_arms")},
                   **{c: "family" for c in ("top_family_magnitude", "top_family_identity", "top_beta",
                                            "top_identity", "top_beta_frac", "n_fam",
                                            "identity_eq_magnitude")}}}


# ── DOMAIN: the condition under which a column is DEFINED. A column NaN on 80% of rows is not
# "missing data" if it is only defined on 20% of rows BY DESIGN — and the difference decides whether a
# NaN is a warning or a fact. VERIFIED (2026-08-01): the stated domain matches the observed non-NaN
# pattern at 98.6-100% for every entry below, so these NaNs are STRUCTURAL.
DOMAIN = {
    "multi-arm cells only (n_arm_in_cell > 1) — 20.3% of edges": (
        "arm_dbeta", "arm_sep_z", "arm_resolvable", "n_arm_in_cell", "oof_rho_arm", "oof_rho_fam",
        "oof_drho", "beta_arm", "sd_arm", "z_arm", "arm_credit_share"),
    "arms with a same-seed SURROGATE (MH-166) — 3.6% of edges": (
        "surrogate_instrument", "surrogate_corr", "coupling_hly_surrogate",
        "coupling_p_hly_surrogate", "coupling_p_hly_resolved"),
    "genes/edges present in the CPTAC cohort": ("cptac_",),          # prefix
    "MULTI-FAMILY genes only (an owner needs >1 family to choose between)": (
        "realization_owner", "static_owner_family", "owner_agrees"),
    "genes with a family-level realization fit": ("family_rho_adj", "realization_identity",
                                                  "is_realization_owner"),
    # ── STATE AVAILABILITY. The progression block is defined only where a state was measured, so its
    # NaN rate tracks cohort coverage, not data loss: GTEx-healthy is the thinnest leg (~38-46% absent),
    # NAT next (~34-40%), the tumour/budget block ~30-34%.
    "arms with a GTEx-HEALTHY leg (the thinnest state — many arms are multi-mapping-collapsed, MH-166)": (
        "coupling_hly", "coupling_p_hly", "coupling_z_hly", "grank_HLY", "dGlobal_HLY_NAT",
        "dGlobal_HLY_TUM", "arm_lfc_HLY_TUM_raw", "arm_lfc_HLY_NAT_raw", "arm_lfc_HLY_TUM_QN",
        "healthy_potential", "healthy_leg", "healthy_uninformative", "share_HLY", "rank_HLY",
        "dose_rank_HLY", "d_rank_HLY_NAT", "d_rank_HLY_TUM", "acquisition_score"),
    "edges with a matched NAT leg (n~104 paired participants)": (
        "coupling_nat", "coupling_p_nat", "coupling_z_nat", "grank_NAT", "dGlobal_NAT_TUM",
        "arm_lfc_NAT_TUM", "share_NAT", "rank_NAT", "dose_rank_NAT", "d_rank_NAT_TUM",
        "dose_comp_retention", "dose_prolif_retention", "dose_confounded"),
    "edges present in the progression panel (arm expressed + state measured)": (
        "grank_TUM", "share_TUM", "rank_TUM", "dose_rank_TUM", "arm_med_rpm", "arm_pct_floor",
        "arm_iqr", "spiker", "arm_id_status", "family_dose_share", "family_role", "wiring_frac",
        "term_ABUND", "term_WIRING", "term_INTERACT", "realization_score", "dose_class",
        "shift_class", "coupling_tum", "coupling_p_tum", "coupling_z_tum", "retention_rho"),
    "genes in the CPTAC-PROTEIN compartment analysis (protein layer, n~101)": ("comp_cptac_prot_",),
    "genes with a gated compartment DRIVER (|rho_raw| >= 0.05 — axiom 5's denominator gate)": (
        "comp_tcga_mrna_driver", "comp_tcga_mrna_driver_ret", "comp_tcga_mrna_driver_share"),
    "genes with a NAT-TUM gene-level leg": ("gene_lfc_NAT_TUM", "gene_iqr", "gene_med_rpm",
                                            "gene_pct_floor", "frac_gain_true_acquired"),
    "genes with >=2 seed families (a collinearity measure needs 2 columns)": ("ctx_d_collin",),
    "the WITHIN-PATIENT paired subset (n=103 matched tumour/NAT pairs, MH-158)": (
        "edge_rho_adj", "dShare_M_own", "dShare_raw_own", "mean_own_shift", "mean_dGlobalRank",
        "own_specific_frac"),
    # ── ARM CARD (MH-209). Each block covers a DIFFERENT scan, so its NaNs mean UNSCANNED, never zero.
    # Measured on the denominator that matters (803 model arms, ≥1 HE edge): seq 65.8% · site 96.0% ·
    # ts 33.7% · cur 99.0%. On the union of all 3,241 arms the same numbers read 23.0/23.8/9.9/29.7%.
    "arms in the scanMiR GENOME-WIDE K_D scan (746 arms; breast-expressed genes)": ("seq_",),
    "arms in the scanMiR site-type table (771 arms) — ⚠ HALLMARK-SCOPED, 1,432 genes only": ("site_",),
    "arms in TargetScan's human default predictions (⚠ only 321 arms)": ("ts_",),
    "arms with miRTarBase curation / ledger PMIDs — ⛔ FAME axes, not targetome (MH-208)": (
        "cur_", "fame_"),
    "arms detected in the TCGA miRNA matrix (2,236 of 3,241)": ("abund_",),
    "arms carrying >=1 curated HE edge in the model design (803 arms)": ("model_",),
    "arms with a resolvable precursor locus / RISC-loading measurement": (
        "host_class", "clustered", "mirgenedb", "ago_reads"),
}


def domain_of(col: str) -> str:
    for dom, cols in DOMAIN.items():
        for c in cols:
            if col == c or (c.endswith("_") and col.startswith(c)):
                return dom
    return ""


def rung_of(card: str, col: str) -> str:
    spec = CARDS[card]
    if col in spec["explicit"]:
        return spec["explicit"][col]
    for p, r in spec["prefix"]:
        if col.startswith(p):
            return r
    return "UNASSIGNED"


def build(card: str) -> pd.DataFrame:
    spec = CARDS[card]
    cols = pd.read_csv(spec["path"], sep="\t", low_memory=False, nrows=3).columns
    R = pd.DataFrame({"card": card, "column": cols,
                      "rung": [rung_of(card, c) for c in cols],
                      "agg_of": [AGG_OF.get(card, {}).get(c, "") for c in cols],
                      "domain": [domain_of(c) for c in cols]})
    R.to_csv(OUT / f"{card}_card_rungs.tsv", sep="\t", index=False)
    return R


def verify(card: str, sample_genes: int = 900) -> pd.DataFrame:
    """Test each declared rung against the data's invariance, GIVEN THIS CARD'S KEY."""
    spec = CARDS[card]
    d = pd.read_csv(spec["path"], sep="\t", low_memory=False)
    if "gene" in d and d.gene.nunique() > sample_genes:
        d = d[d.gene.isin(d.gene.drop_duplicates().sample(sample_genes, random_state=0))]
    R = build(card)
    grp = {"gene": ["gene"],
           "family": ["gene", "seed_family"] if card == "edge" else None,
           "arm": ["arm"] if card == "edge" else None}
    out = []
    for _, r in R.iterrows():
        g = grp.get(r["rung"])
        if not g or r["column"] not in d.columns or any(k not in d.columns for k in g):
            continue
        # NaN MUST NOT COUNT AS A VALUE. `nunique(dropna=False)` made every partially-covered column
        # look like it violated its rung: `healthy_leg`, `spiker`, `gene_iqr`, `gene_lfc_NAT_TUM` all
        # showed ~38% "violation" purely because they are NaN on 30% of rows (a left-join gap), and all
        # four drop to EXACTLY 0% once NaN is excluded. The labels were right; the CHECK was wrong.
        # Coverage is reported separately — a partially-covered column is worth knowing about, but it
        # is not a unit mismatch.
        nu = d.groupby(g, dropna=False)[r["column"]].nunique(dropna=True)
        v = int((nu > 1).sum())
        out.append({"card": card, "column": r["column"], "rung": r["rung"],
                    "groups": int(len(nu)), "violations": v, "rate": v / max(len(nu), 1),
                    "na_rate": float(d[r["column"]].isna().mean())})
    return pd.DataFrame(out)


def main() -> int:
    allR, gaps = [], 0
    for card in CARDS:
        R, V = build(card), verify(card)
        allR.append(R)
        n_un = int((R.rung == "UNASSIGNED").sum())
        nv = int((V.rate > 0).sum()) if len(V) else 0
        rows = sum(1 for _ in open(CARDS[card]["path"])) - 1
        print(f"\n=== {card.upper()} CARD — {rows:,} rows x {len(R)} cols · key={CARDS[card]['key']} ===")
        print("   rung:", dict(R.rung.value_counts()))
        if (R.agg_of != "").any():
            print("   agg_of (row lives here, but SUMMARISES a lower unit):",
                  dict(R[R.agg_of != ""].agg_of.value_counts()))
        print(f"   UNASSIGNED {n_un} · invariance violations {nv}"
              + (f" of {len(V)} checkable" if len(V) else " (none checkable on this grain)"))
        # ⚠ `V` is EMPTY when nothing is invariance-checkable on this grain — true for any single-key card
        # whose columns are all its own rung (the ARM card, MH-209, is the first such card). An empty frame
        # has no `.rate`, so the guard is required, not defensive: without it `--check` died with
        # `AttributeError: 'DataFrame' object has no attribute 'rate'` AFTER printing a clean-looking
        # per-card summary — a checker that exits 1 mid-report reads as "checked and fine" to a skimmer.
        if len(V):
            for _, b in V[V.rate > 0].sort_values("rate", ascending=False).head(8).iterrows():
                print(f"      X {b.column:26s} {b.rung:7s} varies in {b.rate:.0%} of groups")
        # ⚠ coverage is computed over EVERY column, not just the rung-checkable ones. It used to read
        # off `V`, which only holds gene/family/arm rungs — so it reported 25 thin columns on the edge
        # card where the true number is 89. A coverage report that silently skips two thirds of the
        # table is worse than none.
        _d = pd.read_csv(CARDS[card]["path"], sep="\t", low_memory=False)
        _na = _d.isna().mean()
        thin = [c for c in _d.columns if _na[c] > 0.25]
        if thin:
            R2 = R.set_index("column")
            named = [c for c in thin if R2.loc[c, "domain"]]
            print(f"   coverage: {len(thin)} of {len(R)} column(s) NaN on >25% of rows — "
                  f"{len(named)} carry a DECLARED, VERIFIED DOMAIN (structural, not missing data), "
                  f"{len(thin) - len(named)} do not")
            for c in thin:
                if not R2.loc[c, "domain"]:
                    print(f"      ? {c:28s} {_na[c]:.0%} NaN — no domain declared")
        gaps += n_un + nv
    A = pd.concat(allR, ignore_index=True)
    A.to_csv(OUT / "card_registry.tsv", sep="\t", index=False)
    dup = A.groupby("column").rung.nunique()
    shared = sorted(dup[dup > 1].index)
    print(f"\n{'=' * 78}\nSAME NAME, DIFFERENT RUNG ACROSS CARDS — {len(shared)} column(s)")
    print("   legitimate, but NEVER join these across cards without re-checking the unit:")
    for c in shared[:14]:
        print("     ", c, dict(A[A.column == c].set_index("card").rung))
    print(f"\n-> {OUT / 'card_registry.tsv'}   ({len(A)} card-column rows)")
    print(f"{'GAPS' if gaps else 'CLEAN'}: total unassigned + invariance violations = {gaps}")
    return gaps


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    n = main()
    if a.check and n:
        raise SystemExit(1)
