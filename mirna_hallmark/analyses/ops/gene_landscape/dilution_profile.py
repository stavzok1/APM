"""⭐ WHICH ADDED MATES DILUTE THE FAMILY AGGREGATE, AND WHY? (user-asked 2026-08-19)

    PYTHONPATH=. .venv/bin/python3 -m mirna_hallmark.analyses.ops.gene_landscape.dilution_profile

`excluded_mates_evidence.py` PART B found that adding a gene's excluded same-seed mates back into the
family aggregate is a WASH pooled — 56.1% improve, 40.9% dilute. That number hides the mechanism. This
asks what the DILUTING arms are, on the two axes the user named: ABUNDANCE and EVIDENCE.

⚠ THE UNIT IS THE ADDED ARM, NOT THE CELL. A cell can dilute because of one bad mate among five, so a
per-cell contrast cannot say what a diluting arm IS. Every row here is one (gene, family, added arm), and
the cell's Δ is attributed to it — with the per-arm coupling measured directly so the mechanism is
testable rather than assumed.

THE HYPOTHESIS TESTED, AND ⛔ HALF OF IT WAS WRONG. I predicted "ABUNDANT but UNCOUPLED dilutes", on the
reasoning that the aggregate is a dose SUM so an abundant non-repressor enters the predictor with weight.
**The coupling half holds and the ABUNDANCE half does not.**

RESULT (660 cells / 1,420 added-arm rows, 2026-08-19)
-----------------------------------------------------
**HOW WIDESPREAD — much less than the raw 40.9% suggests.** Most of that is noise around zero:
**material dilution (Δ > +0.05) is 67 of 660 cells = 10.2%**, against 98 (14.8%) material gain; at
|Δ| > 0.20 it is 8 (1.2%) vs 13 (2.0%). ⇒ quote **~10%**, not 41%.

**WHAT A DILUTING ARM IS — one axis, and it is not abundance:**

    its own coupling rho     diluting −0.0020  vs improving −0.0870   p=1.4e-22   <- THE axis
    abundance (log2 RPM)              7.10             7.70           p=0.61      <- NO effect
    study attention (PMIDs)             54               54           p=0.37      <- NO effect
    evidence state           HT_ONLY 43.0% · NEVER_LOOKED 43.2% · WEAK_ONLY 48.3%  chi2 p=0.72

The 2×2 settles it — dilution tracks coupling and is FLAT in abundance:

                       coupled   UNCOUPLED
        sparse           31.8%      55.2%
        abundant         32.0%      52.5%

⇒ **a mate that does not itself couple dilutes; whether it is abundant is irrelevant.** My dose-sum story
predicted an abundance gradient and there is none.

⛔⛔ **THE CONSEQUENCE FOR THE REFIT (board 6b leg c): EVIDENCE CANNOT PRE-SCREEN WHICH MATES TO ADD.**
Evidence state does not predict dilution at all (p=0.72), and neither does fame or abundance. The only
axis that separates them is **the arm's own coupling — which is the outcome**, so screening on it is
selection-on-outcome and would manufacture the very gain the refit is meant to test. ⇒ either add all
mates and accept the wash, or screen on an **out-of-fold** coupling estimate. Do not screen on evidence
class and expect it to help.

**BY FAMILY** — dilution is family-structured: miR-196-5p (4/4 cells dilute), miR-302 (58%, median 8 mates
added), miR-29-3p (58%), miR-34/449 (65%) vs the gainers miR-99/100 (median Δ −0.099), miR-146-5p (12.5%
dilute), miR-141/200a (25%), miR-200bc/429 (7.7%).
"""
from __future__ import annotations

import pathlib
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

OUT = pathlib.Path(__file__).resolve().parents[3] / "output" / "learned"


def main() -> None:
    A = pd.read_csv(OUT / "family_aggregate_with_mates.tsv", sep="\t")
    EV = pd.read_csv(OUT / "excluded_mates_evidence.tsv", sep="\t")
    a = pd.read_csv(OUT / "arm_card.tsv", sep="\t", low_memory=False)
    num = lambda s: pd.to_numeric(s, errors="coerce")

    print("=" * 96)
    print("1. HOW WIDESPREAD IS DILUTION, AND HOW BIG?")
    print("=" * 96)
    d = A.delta
    print(f"  cells scored {len(A)}   improved {(d<0).mean():.1%}   diluted {(d>0).mean():.1%}   exactly 0 {(d==0).mean():.1%}")
    for cut in (0.01, 0.05, 0.10, 0.20):
        print(f"    |Δ| > {cut:<5}:  improved {(d < -cut).sum():4d} ({(d < -cut).mean():5.1%})   "
              f"diluted {(d > cut).sum():4d} ({(d > cut).mean():5.1%})")
    print(f"  ⇒ MATERIAL dilution (Δ > +0.05) affects {(d>0.05).sum()} of {len(A)} cells ({(d>0.05).mean():.1%});")
    print(f"    material gain    (Δ < −0.05) affects {(d<-0.05).sum()} ({(d<-0.05).mean():.1%}).")

    print("\n  BY FAMILY — which families dilute (cells with >=4 scored):")
    fam = A.groupby("seed_family").agg(n=("delta", "size"), med=("delta", "median"),
                                       dil=("delta", lambda s: (s > 0).mean()),
                                       med_added=("n_added", "median"))
    fam = fam[fam.n >= 4].sort_values("med", ascending=False)
    print(fam.head(10).to_string())
    print("  ... and the families that GAIN:")
    print(fam.tail(8).to_string())

    # ---- per added arm ----
    EV = EV.merge(a[["arm", "arm_med_rpm", "fame_npmid"]].drop_duplicates("arm"), on="arm", how="left")
    EV["rpm"] = num(EV.arm_med_rpm)
    M = EV.merge(A[["gene", "seed_family", "delta", "n_in", "n_added", "rho_in"]],
                 on=["gene", "seed_family"], how="inner")
    M["dilutes"] = M.delta > 0
    print("\n" + "=" * 96)
    print(f"2. WHAT CHARACTERISES A DILUTING ADDED ARM?  (n={len(M)} added-arm rows in scored cells)")
    print("=" * 96)
    print(f"{'axis':<26}{'DILUTING':>14}{'IMPROVING':>14}{'MWU p':>12}")
    for lab, col in [("its own coupling rho", "rho"), ("abundance (log2 RPM)", "rpm"),
                     ("study attention (PMIDs)", "fame_npmid"), ("mates added to the cell", "n_added"),
                     ("in-design arms in cell", "n_in")]:
        x, y = num(M.loc[M.dilutes, col]).dropna(), num(M.loc[~M.dilutes, col]).dropna()
        if len(x) < 20 or len(y) < 20:
            continue
        p = stats.mannwhitneyu(x, y).pvalue
        print(f"{lab:<26}{x.median():>14.4f}{y.median():>14.4f}{p:>12.3g}")

    print("\n  EVIDENCE STATE of the added arm:")
    ct = pd.crosstab(M.ev_state, M.dilutes, normalize="index")
    n_by = M.ev_state.value_counts()
    for st in ct.index:
        print(f"    {st:14s} n={n_by[st]:5d}   dilutes {ct.loc[st, True]:.1%}")
    if {True, False} <= set(ct.columns):
        chi = stats.chi2_contingency(pd.crosstab(M.ev_state, M.dilutes))[1]
        print(f"    chi2 p={chi:.3g}")

    print("\n" + "=" * 96)
    print("3. ⭐ THE HYPOTHESIS TEST — abundant-but-uncoupled is the diluting profile?")
    print("=" * 96)
    M["ab_hi"] = M.rpm > M.rpm.median()
    M["cpl_flat"] = num(M.rho) > -0.05          # essentially uncoupled
    tab = M.groupby(["ab_hi", "cpl_flat"]).agg(n=("dilutes", "size"), dil=("dilutes", "mean"))
    tab.index = [f"{'abundant' if a else 'sparse':9s} × {'UNCOUPLED' if c else 'coupled':9s}" for a, c in tab.index]
    print(tab.to_string())
    hi_flat = M[(M.ab_hi) & (M.cpl_flat)].dilutes.mean()
    hi_cpl = M[(M.ab_hi) & (~M.cpl_flat)].dilutes.mean()
    sp_flat = M[(~M.ab_hi) & (M.cpl_flat)].dilutes.mean()
    print(f"\n  uncoupled mates dilute at {hi_flat:.1%} (abundant) and {sp_flat:.1%} (sparse);")
    print(f"  coupled mates dilute at {hi_cpl:.1%} (abundant). ⇒ THE AXIS IS COUPLING, NOT ABUNDANCE.")
    print("  ⛔ My dose-sum prediction ('abundant non-repressors dilute') expected an abundance")
    print("     gradient. There is none (p=0.61) — a mate that does not couple dilutes at any abundance.")
    print("\n  ⛔⛔ CONSEQUENCE FOR THE REFIT: evidence state does NOT predict dilution (chi2 p=0.72),")
    print("     nor does fame or abundance. The only separating axis is the arm's OWN COUPLING — the")
    print("     outcome — so screening on it is selection-on-outcome. Use an OUT-OF-FOLD estimate,")
    print("     or add all mates and accept the wash.")


if __name__ == "__main__":
    main()
