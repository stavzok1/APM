"""MH-142 — is E6's self-partialling bug LIVE in production? NO. Verified, not argued.

MH-124's E6 row says: "SELF-PARTIALLING BUG (live): `discovery`/`dossier` partial every edge on
`C + he_agg` ... an HE arm is partialled on a covariate containing itself". Its evidence column, however,
reads "class-matched estimator comparison" -- i.e. the effect was measured in an experiment that deliberately
pushed HE arms through the estimator to class-match them against candidates and fakes.

Production never does that: every call site scores ORPHANS, which are disjoint from `he_agg` by construction.
This script is the check. Run: .venv/bin/python3 -m mirna_hallmark.eval._e6_production_check
"""
import pandas as pd
from mirna_hallmark.learned.eval import dossier as D
from mirna_hallmark.learned.evidence import ledger as LG

c = pd.read_csv(D.COMBINED, sep="\t")
he = LG.pooled_he_edges()
hepairs = set(zip(he["gene"], he["miRNA"]))
cp = set(zip(c["gene"], c["arm"]))
inter = cp & hepairs
print(f"HE edges (pooled_he_edges)         : {len(hepairs)}")
print(f"dossier candidate (gene,arm) pairs : {len(cp)}")
print(f"HE arms scored through he_agg      : {len(inter)} ({len(inter)/max(len(cp),1):.2%})")
assert not inter, f"E6 IS LIVE: {len(inter)} HE arms scored against an he_agg containing them"
print("\nOK — production scores only orphans; he_agg cannot contain them; E6 is NOT live.")
print("E6 remains a REAL trap for any CLASS-MATCHED comparison that scores HE arms through this estimator")
print("(as MH-124's own E6 experiment did). Do not remove `he_agg`: it is what makes an orphan's coupling")
print("PARTIAL on what curation already captures, which is the whole point of the discovery lane.")
