"""Site-primed cooperative product terms (Design §Decision G). STUB.

Add + ψ_{g,(m,m')}·occ_m·occ_m' ONLY for arm-pairs with DISTINCT predicted sites 8–40 nt apart on the
same 3'UTR (Saetrom/Grimson/Broderick window), support primed by TargetScan/TarBase site coords.
Strong spike-and-slab, low base rate; promote a product only if it beats its mains out-of-fold.
Seed redundancy (same/overlapping site, sub-additive) is handled by the family grouping (§F), not here.
"""
raise_note = "Phase-G: cooperative products primed by 8–40 nt distinct-site pairs (site coords)."
