"""Per-arm miRNA promiscuity D(m) — TWO definitions, and ⛔ **the original one is a FAME axis.**

⛔⛔ **READ THIS FIRST (MH-208, 2026-08-03). `he_degree*` / `promisc_he*` — everything built below by
`build_targetome_promiscuity` — is a CURATION COUNT, and it measures how well-studied an arm is, not how
broad its targetome is.** Measured over the 572 arms carrying both definitions:

| | vs FAME (#distinct PMIDs) | vs arm abundance | vs the SEQUENCE targetome |
|---|---|---|---|
| **curated `he_degree_expr`** | **ρ=+0.736** (p=1e-163) | +0.556 | **+0.124** (any-site: **−0.011**) |
| **sequence `seq_degree_strong`** | +0.197 | **−0.004** (p=0.91) | — |

The two definitions' **top-10 lists do not overlap at all**: curated = miR-21 · miR-125b · miR-155 ·
miR-29a/b · miR-34a · miR-145 · miR-27a · miR-17 · miR-26a (the cancer-miRNA hall of fame); sequence =
miR-3662 · miR-149-3p · miR-4428 · miR-153-5p · miR-590-3p · miR-205-3p · miR-4784 · miR-30c-2-3p …
And the scales are incommensurate: the curated **median arm has 2 targets**, the sequence median is **3,634**
— curated captures ~0.05% of the actual targetome.

⚠ **"Genome-wide" in the original docstring meant "not Hallmark-scoped", i.e. curated targets across the
genome — NOT the genome.** The name describes the SCOPE of the curation, not the quantity (axiom 7:
classify the quantity, not its container). ⇒ the retired pressure heuristic's specificity weight
`w(m,g)=log_ev(m,g)/D(m)` was **dividing by fame**, which is part of why it measured ≈ abundance ≈ null
(MH-115).

**USE `promisc_seq_strong` (→ `build_sequence_promiscuity`) as the promiscuity axis.** It is curation-free
(scanMiR RBNS K_D on sequence) and — the property that makes it worth having — **orthogonal to abundance
(ρ=−0.004)**, the confound that dominates every other regulator axis. ⚠ Its coverage is the **746
K_D-scanned detectable arms**, not all ~2,600 (see PROGRAM_FORWARD_BOARD §E, "K_D scan scope").

--- original docstring, for the curated definition only ---
Genome-wide, breast-expression-filtered miRNA promiscuity (the corrected attribution denominator).

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


SEQ_CACHE = C.OUTPUT_ROOT / "genomewide_promiscuity" / "seq_targetome_promiscuity.tsv"
SCANMIR_FILES = ("genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz")
_SEQ_COLS = {"seq_degree_strong", "seq_degree_any", "promisc_seq_strong", "promisc_seq_any"}


def build_sequence_promiscuity(floor: float = EXPR_FLOOR, chunksize: int = 2_000_000) -> pd.DataFrame:
    """⭐ The CURATION-FREE targetome size per arm — the promiscuity definition to actually use (MH-208).

    Counts the breast-expressed genes carrying a scanMiR RBNS-K_D site for each arm: `seq_degree_strong`
    (a STRONG site, `repression_strong` non-null) and `seq_degree_any`. Sequence + biochemistry only, so
    unlike `he_degree*` it does not inherit the literature's attention (fame ρ=+0.197 vs +0.736) and is
    **abundance-orthogonal (ρ=−0.004)**.

    ⚠ Coverage is the K_D-scanned set (**746 arms**), not the full ~2,600 — an arm absent here is
    UNSCANNED, not un-promiscuous. Callers must not read a missing value as zero breadth.
    """
    expr = _breast_expressed_genes(floor)
    strong: dict = {}
    anys: dict = {}
    for f in SCANMIR_FILES:
        path = Path("data/external_cache/scanmir") / f
        if not path.exists():
            print(f"  ⚠ missing {path} — sequence promiscuity will UNDER-count")
            continue
        for chunk in pd.read_csv(path, sep="\t", usecols=["arm", "gene", "repression_strong"],
                                 chunksize=chunksize):
            chunk = chunk[chunk["gene"].isin(expr)]
            for a, g in chunk[chunk["repression_strong"].notna()].groupby("arm")["gene"].apply(set).items():
                strong.setdefault(a, set()).update(g)
            for a, g in chunk.groupby("arm")["gene"].apply(set).items():
                anys.setdefault(a, set()).update(g)
    out = pd.DataFrame({"arm": sorted(anys),
                        "seq_degree_strong": [len(strong.get(a, ())) for a in sorted(anys)],
                        "seq_degree_any": [len(anys[a]) for a in sorted(anys)]})
    out["promisc_seq_strong"] = np.log1p(out["seq_degree_strong"])
    out["promisc_seq_any"] = np.log1p(out["seq_degree_any"])
    SEQ_CACHE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(SEQ_CACHE, sep="\t", index=False)
    return out.set_index("arm")


def load_promiscuity(arms: Optional[Iterable[str]] = None, *, col: str = "promisc_seq_strong",
                     rebuild: bool = False, fill: str = "median") -> pd.Series:
    """Per-arm promiscuity. ⭐ DEFAULT IS NOW THE SEQUENCE DEFINITION (`promisc_seq_strong`, MH-208).

    ⛔ `col="promisc_he*"` selects the CURATED count, which is a FAME axis (ρ=+0.736 with #PMIDs) — see the
    module docstring. It is kept only for reproducing the retired pressure-arc results.

    `fill`: "median" (a neutral non-zero discount, the legacy behaviour) or "nan" — ⭐ prefer "nan" for the
    sequence column, where a missing arm means UNSCANNED, and imputing the median invents breadth for it.
    """
    seq = col in _SEQ_COLS
    cache = SEQ_CACHE if seq else CACHE
    if rebuild or not cache.exists():
        build_sequence_promiscuity() if seq else build_targetome_promiscuity()
    df = pd.read_csv(cache, sep="\t").set_index("arm")
    s = df[col]
    if arms is not None:
        arms = list(arms)
        s = s.reindex(arms)
        if fill == "median":
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
