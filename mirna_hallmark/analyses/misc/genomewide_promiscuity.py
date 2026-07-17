"""Genome-wide, breast-expression-filtered miRNA promiscuity (the corrected attribution denominator).

The pressure attribution's specificity weight `w(m,g)=log_ev(m,g)/D(m)` discounts an arm by its
promiscuity D(m). Until now D(m) was computed from the **Hallmark-scoped** edge set (`mirna_target_degree`
on `load_mirtar_edges(hs.universe)`) — i.e. only the arm's targets *within the 50 hallmarks*. That
under-counts the true division-of-attention, which is a property of the arm's whole targetome.

This module computes D(m) over the arm's **genome-wide targetome**, with two corrections agreed in
`docs/ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md` §11 (items 2b, 2c):
  - **HE edges only** (`Support Type == 'Functional MTI'`, the strong functional set — NOT the 3.98M
    `Functional MTI (Weak)` records, which swamp everything and are not the targetome we model).
  - **breast-tumor expression filter**: a target that isn't expressed can't be competed for. Keep targets
    with TCGA-tumor median log2>floor (NOT GTEx-only — that would drop cancer-specific genes).

Empirical: HE genome-wide degree is ~1.5-2× the Hallmark-scoped degree and **rank-faithful to it
(Spearman ≈0.92)** — so this is a modest, rank-preserving rescale, not a reordering. (The earlier
all-evidence count was ~84× larger and rank-uncorrelated — an artifact of the weak edges.)

Run:  .venv/bin/python3 -m mirna_hallmark.analyses.misc.genomewide_promiscuity            # build + cache + assess
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D

CACHE = C.OUTPUT_ROOT / "genomewide_promiscuity" / "he_targetome_promiscuity.tsv"
HE_SUPPORT = "Functional MTI"   # strong functional MTI only (HE); excludes the weak/non-functional sets
EXPR_FLOOR = 1.0                # TCGA-tumor median log2(x+1) floor for "expressed in breast"


def _breast_expressed_genes(floor: float = EXPR_FLOOR) -> set:
    rna = D.load_rna()
    med = rna.median(axis=1)
    return set(med[med > floor].index)


def build_targetome_promiscuity(expr_filter: bool = True, floor: float = EXPR_FLOOR) -> pd.DataFrame:
    """Per-arm genome-wide HE promiscuity (cached). Columns: degree, degree_expr, log1p variants."""
    raw = pd.read_csv(C.MIRTARBASE_CSV, usecols=["miRNA", "Target Gene", "Support Type"])
    he = raw[(raw["miRNA"].str.startswith("hsa-", na=False)) & (raw["Support Type"] == HE_SUPPORT)]
    he = he.rename(columns={"Target Gene": "gene"})[["miRNA", "gene"]].drop_duplicates()
    deg = he.groupby("miRNA")["gene"].nunique().rename("he_degree")
    out = deg.to_frame()
    if expr_filter:
        expr = _breast_expressed_genes(floor)
        deg_x = he[he["gene"].isin(expr)].groupby("miRNA")["gene"].nunique().rename("he_degree_expr")
        out = out.join(deg_x, how="left")
        out["he_degree_expr"] = out["he_degree_expr"].fillna(0).astype(int)
    else:
        out["he_degree_expr"] = out["he_degree"]
    out["promisc_he"] = np.log1p(out["he_degree"])
    out["promisc_he_expr"] = np.log1p(out["he_degree_expr"])
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    out.reset_index().rename(columns={"miRNA": "arm"}).to_csv(CACHE, sep="\t", index=False)
    return out


def load_promiscuity(arms: Optional[Iterable[str]] = None, *, col: str = "promisc_he_expr",
                     rebuild: bool = False) -> pd.Series:
    """Per-arm promiscuity discount (log1p of HE [+expr-filtered] genome-wide degree). Cached.

    Use as the attribution specificity denominator in place of the Hallmark-scoped `_promisc_discount`.
    Arms with no HE targetome entry get the median (a neutral, non-zero discount).
    """
    if rebuild or not CACHE.exists():
        build_targetome_promiscuity()
    df = pd.read_csv(CACHE, sep="\t").set_index("arm")
    s = df[col]
    if arms is not None:
        arms = list(arms)
        s = s.reindex(arms)
        s = s.fillna(s.median() if s.notna().any() else 0.0)
    return s


def main() -> None:
    out = build_targetome_promiscuity()
    print(f"[gw-promisc] HE targetome promiscuity for {len(out)} arms -> {CACHE}")
    from mirna_hallmark.coupling_predictor_comparison import _load_grid_inputs
    from mirna_hallmark.pressure_engine import mirna_target_degree
    _hs, edges, *_ = _load_grid_inputs(None)
    deg_hm = mirna_target_degree(edges)
    shared = deg_hm.index.intersection(out.index)
    rho = pd.Series(deg_hm[shared]).corr(pd.Series(out.loc[shared, "he_degree"]), method="spearman")
    print(f"[gw-promisc] Spearman(hallmark degree, HE genome-wide degree) over {len(shared)} arms = {rho:.2f} "
          f"(rank-faithful) | median HE/hallmark ratio = "
          f"{(out.loc[shared,'he_degree']/deg_hm[shared].replace(0,np.nan)).median():.1f}x")
    print("\n  arm                HE_deg  HE_expr  promisc_he_expr")
    for a in ["hsa-miR-33b-5p", "hsa-miR-137", "hsa-miR-17-5p", "hsa-miR-21-5p", "hsa-miR-370-3p", "hsa-miR-375"]:
        if a in out.index:
            r = out.loc[a]
            print(f"  {a:18s} {int(r.he_degree):5d}  {int(r.he_degree_expr):6d}  {r.promisc_he_expr:.3f}")


if __name__ == "__main__":
    main()
