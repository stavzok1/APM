from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from ..config import UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE


def normalize_expression_mat(
    rna: pd.DataFrame,
    gene_col: str = "gene_symbol",
    mapping: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Normalize expression matrix gene symbols (UCSC fixes, HGNC aliases, etc.).

    Replaces gene symbols in ``gene_col`` according to ``mapping``,
    leaving all other symbols unchanged.

    Args:
        rna: RNA expression DataFrame
        gene_col: Column containing gene symbols
        mapping: old_symbol -> new_symbol. If None, uses only
            ``UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE`` from config.

    Returns:
        DataFrame with normalized gene symbols
    """
    if mapping is None:
        mapping = dict(UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE)
    if isinstance(gene_col, Path):
        gene_col = gene_col.as_posix()
    else:
        gene_col = str(gene_col)
    print(f'gene column: {gene_col}')
    print(f'matrix columns: {rna.columns[:2].tolist()}')
    if gene_col not in rna.columns:
        if "sample" in rna.columns:
            rna.rename(columns={"sample": gene_col}, inplace=True)
        else:
            raise ValueError(f"Column '{gene_col}' not found in DataFrame")
    
    rna = rna.copy()
    
    rna[gene_col] = rna[gene_col].replace(mapping)
    
    return rna
