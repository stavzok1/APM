"""MH-138 — is `readouts.share` (=beta_f/Sum beta) contaminated by half-normal prior drift?

WHY: METHODS §5 ("canonical attribution = bagged NNLS") and SYNTHESIS §6b ("dense posterior ->
attribution") appeared to conflict. RATIONALE §2e predicted the mechanism that resolves it: the
half-normal slab N+(0,nu^2) has a STRICTLY POSITIVE mean, so an un-informed family cannot be zeroed
-- it relaxes to the prior -- whereas NNLS returns exact zeros. This script measures whether that
drift is real in the SHIPPED table, and how much of `share`'s denominator it occupies.

Run: .venv/bin/python3 -m mirna_hallmark.eval._share_drift_check   (from repo root)
"""
import pandas as pd, numpy as np

d = pd.read_csv("mirna_hallmark/output/learned/readouts_edges.tsv", sep="\t").dropna(subset=["beta", "z"])
gcol = "gene" if "gene" in d.columns else [c for c in d.columns if "gene" in c.lower()][0]
d["unid"] = d["z"].abs() <= 2

print(f"n edges = {len(d)}")
print(f"beta == 0 exactly        : {(d.beta == 0).sum()}   <- NNLS would zero the un-informed; Gibbs cannot")
print(f"beta  > 0                : {(d.beta > 0).sum()} ({(d.beta > 0).mean():.1%})")
print(f"POOLED unidentified mass : {d.loc[d.unid,'beta'].sum() / d.beta.sum():.1%}")

g = d.groupby(gcol).apply(
    lambda x: x.loc[x.unid, "beta"].sum() / x.beta.sum() if x.beta.sum() > 0 else np.nan,
    include_groups=False).dropna()
nfam = d.groupby(gcol).size()
print(f"per-gene median          : {g.median():.1%}")
print(f"genes at 100% drift      : {(g > 0.999).sum()} ({(g > 0.999).mean():.1%}); "
      f"single-family among them: {(nfam[g[g > 0.999].index] == 1).sum()}")
print(f"n_fam>=3 genes median    : {g[nfam.reindex(g.index) >= 3].median():.1%}")
