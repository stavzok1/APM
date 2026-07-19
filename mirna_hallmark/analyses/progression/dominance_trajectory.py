"""Per-TARGET dominance trajectory (MH-166 follow-up) — for every edge (arm→gene), WHEN across GTEx→NAT→TUM does the
arm become a DOMINANT realized repressor of THAT gene, and in which step (FIELD = GTEx→NAT vs MALIGNANT = NAT→TUM)?

⭐ 'how it enters' (the discipline this arc was built under):
  • REALIZED dominance = calibrated per-state coupling (site-free-null p, sign-aware) — NOT dose-budget. The GTEx (hly)
    node reads `coupling_p_hly_resolved` = direct where measured, SAME-SEED SURROGATE where the arm was GTEx-collapsed
    (card.py::_surrogate_map). Where no surrogate exists the hly node is BLIND (unknown), not False — reported as such.
  • NAT and TUM nodes are TCGA/uncollapsed ⇒ available for ALL arms. So the MALIGNANT (NAT→TUM) step is always clean;
    the FIELD (GTEx→NAT) step is clean for measured+surrogate arms, blind for the still-blind minority.
  • dose-budget rank is deliberately NOT the dominance axis here (fixed_M×level is a dose placement, not fitted
    dominance; imputing a constant kills the variance a re-fit needs) — see FORMULAS §11g.

Reads `output/learned/realization/progression_edge_card.tsv`; writes `dominance_trajectory.tsv` + prints the summary.
Run:  .venv/bin/python3 -m mirna_hallmark.analyses.progression.dominance_trajectory
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from mirna_hallmark.learned.realization import OUT

ALPHA = 0.05          # calibrated-repressor call (site-free null p)
TOPK = 2              # 'dominant' = a top-K realized repressor of the gene in that state


def _state_rho_p(ec: pd.DataFrame):
    """resolved per-state (rho, p): hly reads direct-or-surrogate; nat/tum direct."""
    hly_rho = ec["coupling_hly"].fillna(ec.get("coupling_hly_surrogate"))
    hly_p = ec["coupling_p_hly_resolved"] if "coupling_p_hly_resolved" in ec else ec["coupling_p_hly"]
    return {"hly": (hly_rho, hly_p, ec["healthy_leg"]),
            "nat": (ec["coupling_nat"], ec["coupling_p_nat"], None),
            "tum": (ec["coupling_tum"], ec["coupling_p_tum"], None)}


def dominance_trajectory() -> pd.DataFrame:
    ec = pd.read_csv(OUT / "progression_edge_card.tsv", sep="\t")
    sr = _state_rho_p(ec)
    for st, (rho, p, _leg) in sr.items():
        # realized-repressor score: -log10(p) where the arm is a calibrated repressor (rho<0), else NaN
        rep = (rho < 0) & (p < ALPHA)
        score = pd.Series(np.where(rep, -np.log10(p.clip(lower=1e-12)), np.nan), index=ec.index)
        ec[f"coup_score_{st}"] = score.round(2)
        ec[f"coup_rank_{st}"] = score.groupby(ec["gene"]).rank(ascending=False, method="min")   # 1 = top repressor
        ec[f"dom_{st}"] = rep & (ec[f"coup_rank_{st}"] <= TOPK)
    # hly node availability (blind = collapse_blind AND no resolved p)
    hly_p = sr["hly"][1]
    ec["hly_blind"] = ec["healthy_leg"].eq("collapse_blind") & hly_p.isna()

    def _timing(r):
        if not r["dom_tum"]:
            return "not_dominant"
        if r["hly_blind"]:                                   # can't see the field node
            return "malignant_dom" if not r["dom_nat"] else "field_or_earlier_dom(hly-blind)"
        if r["dom_hly"]:
            return "constitutive_dom"                        # dominant already in true healthy
        if r["dom_nat"]:
            return "field_acquired_dom"                      # became dominant in GTEx→NAT (field)
        return "malignant_acquired_dom"                      # became dominant NAT→TUM
    ec["dominance_timing"] = ec.apply(_timing, axis=1)

    out = OUT / "dominance_trajectory.tsv"
    keep = ["gene", "arm", "shift_class", "healthy_leg", "surrogate_instrument",
            "coup_score_hly", "coup_score_nat", "coup_score_tum",
            "coup_rank_hly", "coup_rank_nat", "coup_rank_tum",
            "dom_hly", "dom_nat", "dom_tum", "hly_blind", "dominance_timing"]
    ec[keep].to_csv(out, sep="\t", index=False)

    dom = ec[ec.dom_tum]
    print(f"[dominance_trajectory] {len(ec)} edges | tumour-dominant repressors: {len(dom)} across {dom.gene.nunique()} genes")
    print("\ntiming of tumour-dominance (WHEN the arm became a top realized repressor of the gene):")
    print(dom.dominance_timing.value_counts().to_string())
    # of the CLEAN (non-blind) tumour-dominant edges, field vs malignant vs constitutive
    clean = dom[~dom.hly_blind]
    print(f"\nclean (hly node readable) tumour-dominant: {len(clean)}")
    for t in ["constitutive_dom", "field_acquired_dom", "malignant_acquired_dom"]:
        print(f"  {t:24s} {int((clean.dominance_timing==t).sum())}")
    print("\ntop FIELD-acquired dominant edges (became a top repressor already by NAT):")
    fa = clean[clean.dominance_timing == "field_acquired_dom"].sort_values("coup_score_tum", ascending=False)
    print(fa[["gene", "arm", "coup_score_hly", "coup_score_nat", "coup_score_tum", "shift_class"]].head(12).to_string(index=False))
    print(f"\n-> {out}")
    return ec[keep]


if __name__ == "__main__":
    dominance_trajectory()
