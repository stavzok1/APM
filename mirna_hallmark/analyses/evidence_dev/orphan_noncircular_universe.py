"""TASK 2 — build the NON-CIRCULAR orphan universe and its site-free abundance-matched null.

Nomination trace (measured, see module docstring of the runner):
  every persisted "orphan" (discovery_dossier.tsv, 6744 edges) passed an EXPRESSION-derived selection
  (partial-Spearman on TCGA mRNA, or a spike-slab PIP / dense-ridge |z| from a Gibbs posterior on the
  TCGA-residualised design). in_targetscan / in_kd are only the CANDIDACY SOURCE, not an independent route.

Here we rebuild the SEQUENCE-ONLY candidacy universes exactly as `data.assemble_gene(orphans=True, ...)`
constructs them -- but WITHOUT the expression filter:
  TS route : TargetScan context++ site (ts_mag > 0)                         [pure sequence]
  KD route : scanMiR predicted repression < 0, top-80 per gene              [pure sequence/biophysics]
restricted to arms in the miRNA matrix, genes in the mRNA matrix and the pooled-HE gene universe
(the genes the discovery scan actually visited), EXCLUDING every curated pooled-HE edge.

NULL: degree-preserving rewiring of the orphan universe onto SITE-FREE pairs (no TS site, no scanMiR K_D
entry at all, not curated). Rewiring by arm keeps each arm's edge count EXACTLY -> the arm ABUNDANCE and
DETECTION marginals are matched by construction (they are arm properties).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import kd as KD
from mirna_hallmark.learned.evidence import ledger as LG

OUT = Path("mirna_hallmark/output/learned/orphan_noncircular")


def build_universes(kd_top: int = 80):
    d = LD._load()
    X, Y = d["X"], d["Y"]
    arms_x = set(X.index)
    genes_y = set(Y.index)

    he = LG.pooled_he_edges()
    he_genes = sorted(set(he["gene"].dropna().astype(str)) & genes_y)  # the scanned gene universe
    hg = set(he_genes)
    curated = set(zip(he["gene"].astype(str), he["miRNA"].astype(str)))

    # ---- SEQUENCE route 1: TargetScan context++ sites
    ts = LD._targetscan_context()
    ts = ts[(ts["ts_mag"] > 0) & ts["arm"].isin(arms_x) & ts["gene"].isin(hg)]
    ts_pairs = set(zip(ts["gene"].astype(str), ts["arm"].astype(str)))

    # ---- SEQUENCE route 2: scanMiR K_D (same rule assemble_gene uses: repression<0, top-80 per gene)
    aff = KD.affinity()
    aff = aff[aff["arm"].isin(arms_x) & aff["gene"].isin(hg) & (aff["repression"] < 0)]
    aff = aff.drop_duplicates(["gene", "arm"])
    kd = aff.sort_values("repression").groupby("gene", sort=False).head(kd_top)
    kd_pairs = set(zip(kd["gene"].astype(str), kd["arm"].astype(str)))

    seq = (ts_pairs | kd_pairs) - curated                       # NON-CIRCULAR orphan universe
    # ---- SITE-FREE pool: pairs with NO sequence support of ANY kind (genome-wide K_D, not the top-80)
    ga = KD.genome_affinity()
    ga = ga[ga["arm"].isin(arms_x) & ga["gene"].isin(hg)]
    kd_any = set(zip(ga["gene"].astype(str), ga["arm"].astype(str)))   # ANY scanMiR entry (even weak)
    ts_any_full = LD._targetscan_context()
    ts_any = set(zip(ts_any_full["gene"].astype(str), ts_any_full["arm"].astype(str)))
    has_site = ts_any | kd_any | curated                        # anything with sequence support or curation

    return dict(X=X, Y=Y, he_genes=he_genes, curated=curated,
                ts_pairs=ts_pairs, kd_pairs=kd_pairs, seq=seq, has_site=has_site)


def matched_null(seq: set, has_site: set, he_genes: list, arms_x: set, *, seed: int = 0) -> set:
    """Degree-preserving null: keep each arm's edge COUNT, resample its gene partners from site-free genes.
    Arm-level abundance/detection marginals are therefore matched EXACTLY."""
    rng = np.random.default_rng(seed)
    by_arm: dict = {}
    for g, a in seq:
        by_arm.setdefault(a, []).append(g)
    genes = np.array(he_genes)
    out = set()
    for a, gs in by_arm.items():
        need = len(gs)
        pool = [g for g in genes if (g, a) not in has_site]
        if not pool:
            continue
        k = min(need, len(pool))
        picks = rng.choice(len(pool), size=k, replace=False)
        for i in picks:
            out.add((pool[i], a))
    return out


if __name__ == "__main__":
    u = build_universes()
    OUT.mkdir(parents=True, exist_ok=True)
    seq = u["seq"]
    null = matched_null(seq, u["has_site"], u["he_genes"], set(u["X"].index))
    print(f"HE genes scanned      : {len(u['he_genes'])}")
    print(f"curated pooled-HE     : {len(u['curated'])}")
    print(f"TS-site pairs         : {len(u['ts_pairs'])}")
    print(f"scanMiR-KD pairs      : {len(u['kd_pairs'])}")
    print(f"NON-CIRCULAR universe : {len(seq)}  (TS-only {len(u['ts_pairs']-u['kd_pairs']-u['curated'])}, "
          f"KD-only {len(u['kd_pairs']-u['ts_pairs']-u['curated'])}, both {len(u['ts_pairs']&u['kd_pairs']-u['curated'])})")
    print(f"SITE-FREE matched null: {len(null)}")
    pd.DataFrame(sorted(seq), columns=["gene", "arm"]).to_csv(OUT / "universe_seq.tsv", sep="\t", index=False)
    pd.DataFrame(sorted(null), columns=["gene", "arm"]).to_csv(OUT / "universe_null.tsv", sep="\t", index=False)
    pd.DataFrame(sorted(u["curated"]), columns=["gene", "arm"]).to_csv(OUT / "universe_he.tsv", sep="\t", index=False)
