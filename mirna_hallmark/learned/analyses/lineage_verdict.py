"""Per-edge LINEAGE (hormone-receptor / subtype) confound-vs-mechanism verdict — the MH-100 pattern,
applied to the axis the model does NOT control for.

The question
------------
Is a coupling just a BETWEEN-SUBTYPE artifact — "basal tumours happen to have high miR-X and low gene-Y"
— rather than a within-tumour regulatory relationship? That is a legitimate worry and it deserves an
ANSWER PER EDGE. It does **not** deserve a blanket adjustment, for two reasons, both measured:

⚠ **1. PAM50 MUST NEVER GO IN C — it is computed FROM the mRNA matrix.** 27 of the 50 PAM50 classifier
genes (ESR1, ERBB2, MYC, BCL2, EGFR, CCNE1, MKI67, PGR, …) are **targets in our own edge universe**, so
conditioning on the PAM50 label conditions on a *function of Y*. The damage is exactly where the theory
says it should be (measured, 2026-07-12, mean |β| change):

      covariate              all      PAM50-classifier genes    non-PAM50 genes
      + PAM50 dummies       −21%              −36%                  −18%
      + IHC ER/HER2         −11%              −27%                   −8%

  ESR1/miR-18 collapses **−79%** under PAM50 (ESR1 *is* the gene that defines luminal-vs-basal);
  BCL2/miR-224 **−68%**. This reproduces the documented `data._latent` failure — 10 global expression PCs
  killed miR-206→ESR1 (−0.32 → −0.03) — because **PAM50 is a coarse discretization of those same PCs.**

⚠ **2. Even the NON-CIRCULAR covariate is not safe as a global control.** `IHC ER/HER2` is measured by
protein staining (independent of the RNA matrix; 95%/66% coverage) and costs only −8% on non-PAM50
targets — so it is **the right lineage instrument if one is used at all**. But an −8% shift carries exactly
the ambiguity MH-100 resolved for proliferation, where a global fuller control was **REJECTED** because it
over-controlled the majority. Lineage is *more* entangled than proliferation: the miRNAs that define the
luminal/basal axis (miR-200 family, miR-21, miR-18) **are the regulators being studied**.

⇒ **C is unchanged. This module emits a FLAG.** (The sibling of `prolif_verdict`, `coding_pleiotropy`'s
"trust only reductions", and `host_confound_lens` — the program's standing "flag, don't subtract" rule.)

The test (identical, non-circular OOF 2×2 — reused from `prolif_verdict`)
-------------------------------------------------------------------------
    train {C, C+L} × evaluate {C-space, C+L-space}   →   does controlling lineage help HELD-OUT coupling?
      helps in both spaces   → lineage is a genuine CONFOUND      (the coupling was a subtype artifact)
      hurts in both spaces   → lineage is the MECHANISM           (controlling = over-control)
      winner FLIPS           → FRAGILE — OOF cannot settle it

`L` = IHC ER + HER2 (**never PAM50**). Targets that are themselves lineage markers (ESR1, PGR, ERBB2 …)
are flagged `circular=True`: for them **no** lineage covariate is valid — ER status *is* ESR1's own state —
so their verdict is definitionally over-control and must not be read as biology.

CLI: `python -m mirna_hallmark.learned.analyses.lineage_verdict [GENE ...]`
"""
from __future__ import annotations

import sys
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned.prolif_verdict import (TOL, TOL_E, _classify, _edge_verdict, _oof_pred, _resid)

OUT = "mirna_hallmark/output/learned/lineage_verdict.tsv"
EDGE_OUT = "mirna_hallmark/output/learned/lineage_verdict_edges.tsv"
_ANN = "annotations/BRCA_clinical_immune_unified.tsv"

#: The PAM50 classifier genes. A target in this set CANNOT be lineage-adjusted (the covariate and the
#: outcome are the same object); its verdict is definitionally over-control ⇒ `circular=True`.
PAM50_GENES = {"ESR1", "ERBB2", "MYC", "CCNE1", "FOXA1", "MKI67", "BCL2", "EGFR", "CDH3", "KRT5",
               "KRT14", "KRT17", "MMP11", "GRB7", "CEP55", "BIRC5", "CCNB1", "MYBL2", "RRM2", "UBE2C",
               "AURKA", "TYMS", "ANLN", "MELK", "EXO1", "CDC20", "KIF2C", "NDC80", "PTTG1", "UBE2T",
               "ACTR3B", "BAG1", "BLVRA", "CDC6", "CXXC5", "FGFR4", "GPR160", "MAPT", "MDM2", "MIA",
               "NAT1", "PGR", "PHGDH", "SFRP1", "SLC39A6", "TMEM45B", "MLPH", "NUF2"}
_LINEAGE_MARKERS = {"ESR1", "PGR", "ERBB2"}          # the strictest core: the covariate IS the target

_CACHE: dict = {}


def lineage_covariates(parts) -> Optional[pd.DataFrame]:
    """`L` — the NON-CIRCULAR lineage block: IHC ER + HER2, measured by protein staining, independent of
    the RNA matrix the model is fit on. **Deliberately NOT PAM50** (see the module docstring)."""
    if "ihc" not in _CACHE:
        a = pd.read_csv(_ANN, sep="\t").drop_duplicates("participant").set_index("participant")
        mp = {"pos": 1.0, "neg": 0.0}
        _CACHE["ihc"] = pd.DataFrame({"ihc_er": a["ihc_er_status"].map(mp),
                                      "ihc_her2": a["ihc_her2_status"].map(mp)})
    L = _CACHE["ihc"].reindex(parts)
    if L["ihc_er"].notna().mean() < 0.5:
        return None
    return L.apply(lambda s: s.fillna(s.median())).astype(float)


def gene_verdict(gene: str, *, folds: int = 5, n_iter: int = 900) -> Optional[dict]:
    """OOF 2×2: does controlling LINEAGE help this gene's held-out coupling?"""
    Yg, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    if X.shape[1] < 1:
        return None
    L = lineage_covariates(Yg.index)
    if L is None:
        return None
    Caug = C.copy()
    for c in L.columns:
        Caug[c] = L[c].to_numpy()

    yv = Yg.to_numpy(float)
    Cc, Ca = C.to_numpy(float), Caug.to_numpy(float)
    # how much lineage signal is still in the target AFTER C? (the residual crack, MH-100's rho_resid)
    er_resid = _resid(L["ihc_er"].to_numpy(float), Cc)
    rho_resid = spearmanr(er_resid, _resid(yv, Cc)).correlation

    pC = _oof_pred(Yg, X, C, w, folds, n_iter)
    pCR = _oof_pred(Yg, X, Caug, w, folds, n_iter)
    # REPRESSION-DIRECTED coupling = −ρ (the MH-100 abs()→−ρ fix: abs conflated repression-strengthening
    # with sign-FLIPS to anti-repression). Same convention here, deliberately.
    def rho(pred, S):
        return -spearmanr(_resid(pred, S), _resid(yv, S)).correlation

    dC = rho(pCR, Cc) - rho(pC, Cc)
    dCR = rho(pCR, Ca) - rho(pC, Ca)
    cls = _classify(dC, dCR)

    from mirna_hallmark.learned import spike_slab as SS
    M0, _ = SS.fit_gene_ss(Yg, X, C, w, n_iter=n_iter, burn=n_iter // 3)
    M1, _ = SS.fit_gene_ss(Yg, X, Caug, w, n_iter=n_iter, burn=n_iter // 3)
    rel = {f: (float(M1[f] - M0[f]) / abs(float(M0[f]))) if abs(float(M0[f])) > 1e-6 else 0.0
           for f in X.columns}
    circ = gene in PAM50_GENES
    return {"gene": gene, "n": len(Yg), "rho_resid_lineage": round(float(rho_resid), 3),
            "delta_Cspace": round(float(dC), 4), "delta_CLspace": round(float(dCR), 4),
            "class": cls, "circular": circ,
            "circular_reason": ("PAM50 classifier gene — the lineage covariate IS (a function of) this "
                                "target; verdict is definitionally over-control, NOT biology" if circ else ""),
            "mean_rel_dbeta": round(float(np.mean(list(rel.values()))), 3),
            "_rel": rel, "_class": cls}


def run(genes: Optional[Sequence[str]] = None, *, folds: int = 5) -> pd.DataFrame:
    genes = list(genes) if genes else ["PTEN", "ZEB1", "CDK6", "BCL2", "ESR1", "EGFR", "TP63",
                                       "CDKN1A", "MYB", "THBS1", "SERPINB5", "KIT"]
    rows, erows = [], []
    for g in genes:
        try:
            v = gene_verdict(g, folds=folds)
        except Exception as e:
            print(f"  {g}: skipped ({type(e).__name__})"); continue
        if v is None:
            continue
        rel, cls = v.pop("_rel"), v.pop("_class")
        rows.append(v)
        for f, r in rel.items():
            # `_edge_verdict` is reused verbatim from prolif_verdict; only its "robust" label is
            # proliferation-specific, so rename it for this axis (the semantics are identical).
            ev = _edge_verdict(cls, r).replace("prolif_robust", "lineage_robust")
            erows.append({"gene": v["gene"], "family": f, "rel_dbeta": round(r, 3),
                          "gene_class": cls, "circular": v["circular"], "edge_verdict": ev})
    df, ed = pd.DataFrame(rows), pd.DataFrame(erows)
    if df.empty:
        print("no verdicts"); return df
    from pathlib import Path
    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    ed.to_csv(EDGE_OUT, sep="\t", index=False)

    with pd.option_context("display.width", 170):
        print("\n=== LINEAGE verdict (OOF 2×2; L = IHC ER/HER2 — NON-circular, never PAM50) ===")
        print(df[["gene", "n", "rho_resid_lineage", "delta_Cspace", "delta_CLspace", "class",
                  "circular", "mean_rel_dbeta"]].to_string(index=False))
    clean = df[~df["circular"]]
    print(f"\n  GENE classes (excluding the {int(df['circular'].sum())} circular/PAM50 targets): "
          f"{dict(clean['class'].value_counts())}")
    if len(ed):
        ce = ed[~ed["circular"]]
        print(f"  EDGE verdicts (non-circular): {dict(ce['edge_verdict'].value_counts())}")
    print(f"\n  ⚠ `circular=True` rows (PAM50 classifier genes) are NOT interpretable as biology — for them the "
          f"lineage covariate IS the target.\n  → {OUT} · {EDGE_OUT}")
    return df


if __name__ == "__main__":
    run([a for a in sys.argv[1:] if not a.startswith("-")] or None)
