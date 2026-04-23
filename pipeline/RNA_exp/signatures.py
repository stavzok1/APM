from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


BUFFA_HYPOXIA_METAGENE: Tuple[str, ...] = (
    # MSigDB: BUFFA_HYPOXIA_METAGENE (Buffa et al., PLoS One 2010; PubMed 20087356).
    # Stored as upper-case gene symbols for matching the TPM loader normalization.
    "NDRG1",
    "LDHA",
    "PGK1",
    "TPI1",
    "PGAM1",
    "PFKP",
    "ENO1",
    "SLC2A1",
    "PNP",
    "BNIP3",
    "PSRC1",
    "SLC16A1",
    "DDIT4",
    "ADM",
    "HK2",
    "SEC61G",
    "AK4",
    "CA9",
    "P4HA1",
    "ACOT7",
    "GPI",
    "TUBB6",
    "TUBA1C",
    "CDKN3",
    "CTSV",
    "VEGFA",
    "TUBA1A",
    "LRRC42",
    "PSMA7",
    "GAPDH",
    "CHCHD2",
    "YKT6",
    "MIF",
    "MAP7D1",
    "MRPL15",
    "MRPL13",
    "MCTS1",
    "UTP11",
    "KIF4A",
    "HILPDA",
    "MRGBP",
    "KIF20A",
    "MRPS17",
    "ESRP1",
    "SHCBP1",
    "SLC25A32",
    "CORO1C",
    "ANLN",
    "MAD2L2",
    "ANKRD37",
)


AUTOPHAGY_CORE_GENE_SET: Tuple[str, ...] = (
    # Pragmatic “autophagy state” proxy gene set (bulk RNA cannot measure flux directly).
    # Chosen to cover initiation/nucleation/elongation + selective autophagy adaptors.
    "ULK1",
    "ULK2",
    "BECN1",
    "ATG14",
    "PIK3C3",
    "ATG5",
    "ATG7",
    "ATG12",
    "ATG3",
    "ATG4B",
    "WIPI1",
    "WIPI2",
    "MAP1LC3B",
    "GABARAP",
    "GABARAPL1",
    "SQSTM1",
    "NBR1",
    "TFEB",
    "LAMP1",
    "LAMP2",
)


MHCII_READOUT_GENE_SET: Tuple[str, ...] = (
    "HLA-DRA",
    "CD74",
    "CIITA",
)


def compute_mean_signature(expr: pd.DataFrame) -> pd.Series:
    """
    Compute per-sample mean across rows (genes).

    ``expr`` is genes x samples.
    """
    if expr.empty:
        return pd.Series(dtype=float)
    return expr.mean(axis=0, skipna=True)


def compute_cytolytic_score(expr: pd.DataFrame) -> pd.Series:
    """
    CYT = geometric_mean(GZMA, PRF1).

    Accepts a 2-row DataFrame (genes x samples) or a wider DataFrame that contains
    rows named "GZMA" and "PRF1".
    """
    if expr.empty:
        return pd.Series(dtype=float)
    if "GZMA" not in expr.index or "PRF1" not in expr.index:
        missing = [g for g in ("GZMA", "PRF1") if g not in expr.index]
        raise KeyError(f"CYT requires genes {missing} in expression index")

    a = pd.to_numeric(expr.loc["GZMA"], errors="coerce")
    b = pd.to_numeric(expr.loc["PRF1"], errors="coerce")
    x = np.sqrt(a.astype(float) * b.astype(float))
    return pd.Series(x, index=expr.columns, name="CYT")


@dataclass(frozen=True)
class GeneSets:
    apm_class_i: Tuple[str, ...] = (
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "TAP1",
        "TAP2",
        "PSMB8",
        "PSMB9",
        "NLRC5",
    )
    cd8: Tuple[str, ...] = ("CD8A", "CD8B")
    nk: Tuple[str, ...] = ("NKG7", "GNLY", "PRF1", "KLRD1")
    ifng: Tuple[str, ...] = ("IFNG", "CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1")


def default_gene_sets() -> Dict[str, Sequence[str]]:
    """
    Default RNA signature sets used by cohort post-processing.

    Notes:
    - Buffa hypoxia is typically scored as a mean in tumors (see PubMed 39892389 evaluation).
    - Autophagy score is a proxy for transcriptional state, not flux.
    """
    g = GeneSets()
    return {
        "APM_classI_mean": list(g.apm_class_i),
        "CD8_mean": list(g.cd8),
        "NK_mean": list(g.nk),
        "IFNG_mean": list(g.ifng),
        "HYPOXIA_BUFFA_mean": list(BUFFA_HYPOXIA_METAGENE),
        "AUTOPHAGY_core_mean": list(AUTOPHAGY_CORE_GENE_SET),
        "MHCII_readout_mean": list(MHCII_READOUT_GENE_SET),
    }


def _normalize_gene_symbol(x: str) -> str:
    return str(x).strip().upper()


def read_tpm_wide_subset(
    tpm_path: Path,
    *,
    genes: Sequence[str],
    sample_cols: Sequence[str],
    gene_col: str = "gene_symbol",
    sep: str = "\t",
    chunksize: int = 10_000,
) -> pd.DataFrame:
    """
    Read a *wide* TPM table (rows=genes, cols=samples) but only for a small set of genes + samples.

    Returns genes x samples with index = gene symbol (uppercased).
    """
    tpm_path = Path(tpm_path)
    want_genes = {_normalize_gene_symbol(g) for g in genes}
    usecols = [gene_col] + list(sample_cols)
    out_rows: List[pd.DataFrame] = []

    for chunk in pd.read_csv(
        tpm_path, sep=sep, usecols=usecols, chunksize=chunksize, low_memory=False
    ):
        if gene_col not in chunk.columns:
            raise KeyError(f"TPM table missing gene_col={gene_col!r}: {tpm_path}")
        chunk[gene_col] = chunk[gene_col].astype(str).map(_normalize_gene_symbol)
        hit = chunk[chunk[gene_col].isin(want_genes)]
        if not hit.empty:
            out_rows.append(hit)

    if not out_rows:
        return pd.DataFrame(index=[], columns=list(sample_cols), dtype=float)

    df = pd.concat(out_rows, axis=0, ignore_index=True)
    df = df.drop_duplicates(subset=[gene_col], keep="first")
    df = df.set_index(gene_col)

    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.reindex(index=sorted(df.index), columns=list(sample_cols))


def compute_gene_set_scores_from_tpm(
    tpm_path: Path,
    *,
    sample_ids: Sequence[str],
    gene_col: str = "gene_symbol",
    gene_sets: Optional[Mapping[str, Sequence[str]]] = None,
) -> pd.DataFrame:
    """
    Compute simple expression-defined signatures from a wide TPM matrix.

    Output: samples x scores.
    """
    sets = dict(gene_sets or {})
    if not sets:
        sets = default_gene_sets()
    genes_all: List[str] = sorted(
        {_normalize_gene_symbol(x) for xs in sets.values() for x in xs} | {"GZMA", "PRF1"}
    )

    expr = read_tpm_wide_subset(
        tpm_path,
        genes=genes_all,
        sample_cols=list(sample_ids),
        gene_col=gene_col,
    )

    scores: Dict[str, pd.Series] = {}
    for name, gs in sets.items():
        idx = [_normalize_gene_symbol(g) for g in gs]
        present = [g for g in idx if g in expr.index]
        if not present:
            scores[name] = pd.Series(
                [np.nan] * len(sample_ids), index=sample_ids, name=name
            )
        else:
            scores[name] = compute_mean_signature(expr.loc[present]).reindex(sample_ids)

    try:
        scores["CYT"] = compute_cytolytic_score(expr).reindex(sample_ids)
    except KeyError:
        scores["CYT"] = pd.Series([np.nan] * len(sample_ids), index=sample_ids, name="CYT")

    out = pd.DataFrame(scores)
    out.index.name = "sample_id"
    return out

