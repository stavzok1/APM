"""CPTAC prospective data assembly for the PROTEIN axis (Model 3) — the CPTAC twin of `learned/data.py`.

The mediation system (CPTAC_PROTEIN_CHANNEL_PLAN §2), per gene g:

    M = −X·β  + e_M                 β  = mRNA destabilization   (the §6b TCGA latent)
    P = a·M − X·βᵗ + e_P            βᵗ = TRANSLATIONAL repression (mRNA-invisible)

  Model 1  `M ~ X + C`      → β         (the current model, on TCGA)
  Model 2  `P ~ X + C`      → τ = a·β + βᵗ
  Model 3  `P ~ X + M + C`  → βᵗ (coef on X) AND a (coef on M)   ← THE BUILD

**RNA SOURCE — `linkedomics` is the DEFAULT (an earlier `star` default was WRONG; measured, MH-104).**
The two CPTAC mRNA matrices are only **r=0.81** correlated per-gene across samples (1.8× SD difference), so the
choice is a real researcher degree of freedom, not a detail. **Two independent, non-circular measurements both
favour LinkedOmics:**
  1. agreement with the *independent protein assay* (which RNA tracks protein better): LinkedOmics **0.3713** vs
     STAR 0.3678;
  2. **Bar-5 rung (ii)** — how much of the TCGA-learned M's out-of-sample strength survives the cohort jump:
     LinkedOmics retains **0.193** of the TCGA-held-out ceiling vs STAR's **0.148** (binomial p 1.8e-3 vs 1.5e-2,
     581 genes).
The earlier `star` default rested on a *gauge-commensurability* argument (same transform as TCGA's Y ⇒ β_CPTAC
natively poolable with β_TCGA). That argument is now **moot**: the gauge standardizes Y anyway (`zscore_y=True`),
and the cross-cohort β pooling it was for is cancelled (β does not transport — MH-104).
  * ``linkedomics`` — HS_CPTAC_BRCA_2018_RNA_GENE.cct. Gene-median-centered log2 (NEGATIVE values; absolute TPM is
                      NOT recoverable — which is why CIBERSORTx cannot use it). **DEFAULT: the better measurement.**
  * ``star``        — GDC CPTAC-2 STAR, log2(TPM+1). Same transform as TCGA's Y; the **only** usable source for
                      CIBERSORTx (needs absolute linear TPM). Keep for the mandatory both-ways robustness check.

Confounders: `confounders.build_C("cptac", parts)` (8 Wu-major lineages + purity + prolif; REF-A, one builder,
four cohorts) + the gene's own `target_cn`. **PAM50 is NEVER added** (`lineage_verdict`: it is computed FROM the
mRNA — 27/50 of its classifier genes are our targets — and costs −36% of |β|). TMT `plex` is a PROTEIN-ASSAY
artifact and belongs to the PROTEIN equation only (axis N2).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM

ROOT = Path(__file__).resolve().parents[2]
STAR_TPM = ROOT / "output/cibersortx_transfer/mixture_cptac_brca.txt"     # genes × X<case>, LINEAR TPM
LINKEDOMICS_RNA = ROOT / "data/CPTAC/HS_CPTAC_BRCA_2018_RNA_GENE.cct"
CPTAC_GENE_CN = ROOT / "data/CPTAC/HS_CPTAC_BRCA_2018_CNA.cct"
LOCUS_CN = ROOT / "mirna_hallmark/output/matrices/mirna_locus_cn_cptac.tsv.gz"

RNA_SOURCES = ("star", "linkedomics")
_CACHE: dict = {}


# --------------------------------------------------------------------------------------- layers
def arms() -> pd.DataFrame:
    """arm × participant  log2(RPM+1)  (CPTAC-2 prospective miRNA-seq, GDC)."""
    if "arms" not in _CACHE:
        from mirna_hallmark.eval import cptac_validation as CV
        X = CV.load_prospective_mirna_arms()
        _CACHE["arms"] = X[~X.index.duplicated(keep="first")]
    return _CACHE["arms"]


def mrna(source: str = "linkedomics") -> pd.DataFrame:
    """gene × participant mRNA. See the module docstring — the source is a REAL fork, not a detail."""
    if source not in RNA_SOURCES:
        raise ValueError(f"rna_source must be one of {RNA_SOURCES}, got {source!r}")
    key = f"mrna_{source}"
    if key not in _CACHE:
        if source == "star":
            if not STAR_TPM.exists():
                raise FileNotFoundError(f"{STAR_TPM} missing — run scripts/cptac/build_cptac_cibersortx_mixture.py")
            M = np.log2(pd.read_csv(STAR_TPM, sep="\t", index_col=0) + 1.0)   # SAME transform as TCGA's Y
        else:
            M = pd.read_csv(LINKEDOMICS_RNA, sep="\t", index_col=0)
            M = M[~M.index.isin(["IDX"])].apply(pd.to_numeric, errors="coerce")
        _CACHE[key] = M[~M.index.duplicated(keep="first")]
    return _CACHE[key]


def protein() -> pd.DataFrame:
    """gene × participant protein (cohort z-scored TMT-11, LinkedOmics prospective)."""
    if "protein" not in _CACHE:
        from mirna_hallmark.eval import cptac_validation as CV
        P = CV.load_cptac_layers("prospective")["protein_z"]
        _CACHE["protein"] = P[~P.index.duplicated(keep="first")]
    return _CACHE["protein"]


def gene_cn() -> pd.DataFrame:
    """gene × participant ABSOLUTE copy number (CPTAC CNA is a log2 ratio → CN = 2·2^log2r)."""
    if "gene_cn" not in _CACHE:
        c = pd.read_csv(CPTAC_GENE_CN, sep="\t", index_col=0)
        c = c[~c.index.isin(["IDX"])].apply(pd.to_numeric, errors="coerce")
        c = c[~c.index.duplicated(keep="first")]
        if float(c.min().min()) < 0:
            c = 2.0 * np.power(2.0, c)
        _CACHE["gene_cn"] = c
    return _CACHE["gene_cn"]


def locus_cn() -> pd.DataFrame:
    """participant × MI* locus ABSOLUTE CN — the CPTAC exogenous instrument (`build_locus_cn_from_gene_cn.py`;
    nearest-5-gene proxy, validated on TCGA vs ASCAT at median r=0.997)."""
    if "locus_cn" not in _CACHE:
        if not LOCUS_CN.exists():
            raise FileNotFoundError(f"{LOCUS_CN} missing — run scripts/analysis/build_locus_cn_from_gene_cn.py")
        _CACHE["locus_cn"] = pd.read_csv(LOCUS_CN, sep="\t", index_col=0)
    return _CACHE["locus_cn"]


def plex() -> Optional[pd.Series]:
    """TMT plex per participant — a PROTEIN-ASSAY artifact (protein equation only; NOT in the mRNA equation)."""
    if "plex" not in _CACHE:
        try:
            from mirna_hallmark import cptac_batch as CB
            _CACHE["plex"] = CB.plex_map("prospective") if hasattr(CB, "plex_map") else None
        except Exception:
            _CACHE["plex"] = None
    return _CACHE["plex"]


# --------------------------------------------------------------------------------------- assembly
def participants(rna_source: str = "linkedomics") -> list[str]:
    """The analysis cohort: participants with miRNA AND mRNA AND protein AND a confounder block."""
    X, M, P = arms(), mrna(rna_source), protein()
    dec = CF.deconv("cptac")
    return [p for p in P.columns if p in X.columns and p in M.columns and (dec is None or p in dec.index)]


def assemble_gene(gene: str, *, rna_source: str = "linkedomics", he_only: bool = False,
                  target_cn: bool = True, parts: Optional[Sequence[str]] = None):
    """Per-gene CPTAC design for Model 3.

    Returns ``(M, P, X, C)``:
      M : Series[participant]           the gene's own mRNA        (the MEDIATOR)
      P : Series[participant]           the gene's own protein     (the OUTCOME)
      X : DataFrame[participant × arm]  regulator arm log2(RPM+1)  (columns = canonical EDGE arm names)
      C : DataFrame[participant × conf] build_C("cptac") + target_cn   (NO PAM50, NO plex)

    Arm names are resolved through `data._arm_name_map` (the `.N`/case/suffixless-guide fix): a bare
    `m in X.index` test drops 58 CPTAC edge arms.
    """
    Xall, Mall, Pall = arms(), mrna(rna_source), protein()
    if gene not in Mall.index:
        raise KeyError(f"{gene} not in CPTAC mRNA ({rna_source})")
    if gene not in Pall.index:
        raise KeyError(f"{gene} not in CPTAC protein")

    pk = list(parts) if parts is not None else participants(rna_source)

    from mirna_hallmark.learned.evidence import ledger as LG
    ed = LG.pooled_he_edges() if he_only else LG.pooled_he_edges()
    nm = LD._arm_name_map(Xall)
    regs = [m for m in ed.loc[ed["gene"] == gene, "miRNA"].unique()
            if (nm.get(m) or nm.get(str(m).lower()))]
    if len(regs) < 2:
        raise ValueError(f"{gene}: <2 regulator arms present in the CPTAC arm matrix")
    regs_x = [nm.get(m) or nm.get(str(m).lower()) or m for m in regs]

    M = pd.to_numeric(Mall.loc[gene, pk], errors="coerce")
    P = pd.to_numeric(Pall.loc[gene, pk], errors="coerce")
    keep = [p for p, ok in zip(pk, (M.notna() & P.notna())) if ok]        # per-gene missingness (CPTAC is sparse)
    M, P = M.loc[keep], P.loc[keep]

    X = Xall.loc[regs_x, keep].T.astype(float).fillna(0.0)
    X.columns = list(regs)                                                # canonical edge names (matches TCGA)

    C = CF.build_C("cptac", keep)
    if target_cn:
        gc = gene_cn()
        if gene in gc.index:
            tcn = pd.to_numeric(gc.loc[gene].reindex(keep), errors="coerce")
            C = C.join(tcn.rename("target_cn"))
    C = C.apply(pd.to_numeric, errors="coerce")
    C = C.fillna(C.median())                                              # preserve n (as build_C does)
    return M, P, X, C


def family_design(gene: str, **kw):
    """`assemble_gene` + the §8 seed-family collapse → (M, P, X_fam, C, members). X_fam is THE predictor."""
    M, P, X, C = assemble_gene(gene, **kw)
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, wf, members = FAM.collapse_by_family(X, pd.Series(1.0, index=X.columns), fam)
    return M, P, Xf, C, members
