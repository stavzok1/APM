"""WHAT IS THE UNIVERSE OF ARMS? Measure it end to end, from source to model."""
import pathlib

import numpy as np
import pandas as pd

B = pathlib.Path("/sci/labs/michall/stavzok/APM/mirna_hallmark/output/learned")
a = pd.read_csv(B / "arm_card.tsv", sep="\t", low_memory=False)
e = pd.read_csv(B / "realization/edge_card.tsv", sep="\t", low_memory=False)
sf = pd.read_csv(B / "seed_family_card.tsv", sep="\t", low_memory=False)

b = lambda c: a[c].map({True: 1, False: 0, "True": 1, "False": 0}) if c in a.columns else None
n = lambda c: pd.to_numeric(a[c], errors="coerce") if c in a.columns else None

print("THE FUNNEL — every arm the system knows, down to the ones that carry a claim\n")
rows = [
    ("arms on the arm card", len(a), "the universe: one row per mature miRNA arm"),
    ("… with any curated literature", int(n("fame_npmid").notna().sum()), "fame_npmid present"),
    ("… scanned for expression", int(b("cov_expr").sum()), "cov_expr True"),
    ("… with an abundance value", int(n("abund_med").notna().sum()), "abund_med present"),
    ("… EXPRESSED above floor", int(b("tier_expressed").sum()), "tier_expressed True"),
    ("… in a seed family with >1 member", int((n("famrole_n_members") > 1).sum()), "famrole_n_members > 1"),
    ("… with a genomic-context call", int(n("gctx_n_loci").notna().sum()), "gctx_ block present"),
    ("… in MirGeneDB", int(b("gctx_mirgenedb").sum()) if b("gctx_mirgenedb") is not None else 0,
     "high-confidence miRNA set"),
    ("… with AGO-CLIP reads", int(n("ago_reads").notna().sum()), "ago_reads present"),
    ("… with ANY model edge", int(n("arb_n_edges").notna().sum()), "arb_n_edges present"),
    ("… on the delivered EDGE card", e.arm.nunique(), "distinct arms with a delivered edge"),
    ("… with an admissible edge", int((n("adm_n_admissible") > 0).sum()), "adm_n_admissible > 0"),
]
for label, v, how in rows:
    print(f"   {label:<38}{v:>7,}   {v/len(a):>6.1%}   {how}")

print(f"\nSEED FAMILIES: {len(sf):,} rows")
fn = pd.to_numeric(sf.get("fam_n_members"), errors="coerce")
print(f"   singletons                          {int((fn == 1).sum()):>7,}   {(fn==1).mean():>6.1%}")
print(f"   multi-member                        {int((fn > 1).sum()):>7,}   {(fn>1).mean():>6.1%}")
print(f"   largest family                      {int(fn.max()):>7,} members")
print(f"   arms accounted for by families      {int(fn.sum()):>7,}")

print("\nTHE TWO NUMBERS THAT MATTER MOST")
mod = int(n("arb_n_edges").notna().sum())
exp = int(b("tier_expressed").sum())
print(f"   only {mod:,} of {len(a):,} arms ({mod/len(a):.0%}) carry ANY model edge")
print(f"   only {exp:,} of {len(a):,} arms ({exp/len(a):.0%}) are expressed above the floor")
both = int(((n("arb_n_edges").notna()) & (b("tier_expressed") == 1)).sum())
print(f"   and {both:,} are BOTH — the arms a claim can actually rest on")
print(f"   ⇒ every arm-level statistic over the whole card is dominated by the "
      f"{len(a)-exp:,} arms that are not expressed")
