"""GTEx breast miRNA abundance on **TCGA arm names** (MIMAT-primary join).

GTEx small-RNA features are keyed by RNAcentral URS with MIMAT + mirbase_name
annotations. TCGA GDC miRNA-seq is keyed by MIMAT → miRBase mature ID via our
loci file. Coupling and landscape code must join on the **same TCGA arm index**
used for tumor/NAT, not on GTEx display names alone (e.g. ``hsa-miR-375`` on
TCGA vs ``hsa-miR-375-3p`` in GTEx mirbase_name).

Join order per TCGA arm:
  1. **MIMAT** — aggregate all GTEx rows whose MIMAT maps to the arm (primary).
  2. **mirbase_name fallback** — exact / alias / locus-suffix strip on the legacy
     mirbase_name matrix (for arms without a resolved MIMAT in loci).

``healthy_anchor`` already uses MIMAT medians; this module exposes the full
donor×arm matrix for coupling.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.family_normal_reference import (
    GTEX_MIRNA_TPM,
    GTEX_SMALLRNA_ANNOT,
    _gtex_available,
    _gtex_breast_samples,
    _gtex_donor,
    _gtex_urs_to_arm,
)
from mirna_hallmark.mirna_arm_resolve import mirna_arm_alias_candidates

_LOCUS_SUFFIX = re.compile(r"\.\d+$")

_ARM_MAP: Optional[pd.DataFrame] = None


def _strip_locus(arm: str) -> str:
    return _LOCUS_SUFFIX.sub("", str(arm))


@lru_cache(maxsize=1)
def _mimat_to_arm() -> Dict[str, str]:
    from analysis.expression.mirna_target_integration import load_mimat_to_arm

    return load_mimat_to_arm(C.MIRNA_MATURE_LOCI)


@lru_cache(maxsize=1)
def _load_gtex_mimat_donor_matrix() -> pd.DataFrame:
    """MIMAT × GTEx breast donor, log2(TPM+1)."""
    if not _gtex_available():
        return pd.DataFrame()
    annot = pd.read_csv(
        GTEX_SMALLRNA_ANNOT, sep="\t", usecols=["id", "MIMAT"], low_memory=False,
    )
    annot = annot.dropna(subset=["MIMAT"])
    annot = annot.loc[annot["MIMAT"].astype(str).str.startswith("MIMAT")]
    id_to_mimat = dict(zip(annot["id"].astype(str), annot["MIMAT"].astype(str)))
    breast = set(_gtex_breast_samples("SMLRNA"))
    rows: List[pd.DataFrame] = []
    for chunk in pd.read_csv(GTEX_MIRNA_TPM, sep="\t", index_col=0, chunksize=400, low_memory=False):
        keep = chunk.loc[chunk.index.astype(str).isin(id_to_mimat)]
        if not keep.empty:
            rows.append(keep)
    if not rows:
        return pd.DataFrame()
    mat = pd.concat(rows)
    mat.index = mat.index.astype(str).map(id_to_mimat)
    cols = [c for c in mat.columns if c in breast]
    mat = mat[cols].groupby(level=0).mean()
    mat.columns = [_gtex_donor(c) for c in mat.columns]
    mat = mat.loc[:, ~mat.columns.duplicated()]
    return np.log2(mat.clip(lower=0) + 1)


@lru_cache(maxsize=1)
def _load_gtex_mirbase_donor_matrix() -> pd.DataFrame:
    """mirbase_name × GTEx breast donor, log2(TPM+1) — homolog fallback."""
    if not _gtex_available():
        return pd.DataFrame()
    urs_to_arm = _gtex_urs_to_arm()
    arm_urs = set(urs_to_arm)
    breast = set(_gtex_breast_samples("SMLRNA"))
    rows: List[pd.DataFrame] = []
    for chunk in pd.read_csv(GTEX_MIRNA_TPM, sep="\t", index_col=0, chunksize=400, low_memory=False):
        keep = chunk.loc[chunk.index.astype(str).isin(arm_urs)]
        if not keep.empty:
            rows.append(keep)
    if not rows:
        return pd.DataFrame()
    mat = pd.concat(rows)
    mat.index = mat.index.astype(str).map(urs_to_arm)
    cols = [c for c in mat.columns if c in breast]
    mat = mat[cols].groupby(level=0).mean()
    mat.columns = [_gtex_donor(c) for c in mat.columns]
    mat = mat.loc[:, ~mat.columns.duplicated()]
    return np.log2(mat.clip(lower=0) + 1)


def _resolve_mirbase_feature(tcga_arm: str, mirbase_index: set[str]) -> Tuple[Optional[str], str]:
    for i, cand in enumerate(mirna_arm_alias_candidates(tcga_arm)):
        if cand in mirbase_index:
            if i == 0:
                return cand, "mirbase_exact"
            if cand.endswith(".1") or cand.endswith(".2") or cand.endswith(".3"):
                return cand, "mirbase_isoform"
            return cand, "mirbase_alias"
    base = _strip_locus(tcga_arm)
    if base != tcga_arm and base in mirbase_index:
        return base, "locus_strip"
    return None, "missing"


def gtex_arm_map_table() -> pd.DataFrame:
    """Per TCGA arm: ``gtex_feature`` + ``gtex_arm_map_kind`` (mimat / mirbase_* / missing)."""
    global _ARM_MAP
    if _ARM_MAP is not None:
        return _ARM_MAP.copy()
    gtex_tcga_arm_matrix()  # builds cache
    return _ARM_MAP.copy() if _ARM_MAP is not None else pd.DataFrame()


def gtex_tcga_arm_matrix() -> pd.DataFrame:
    """TCGA miRBase arm × GTEx breast donor, log2(TPM+1)."""
    global _ARM_MAP
    if not _gtex_available():
        _ARM_MAP = pd.DataFrame(columns=["gtex_feature", "gtex_arm_map_kind"])
        return pd.DataFrame()

    mimat_mat = _load_gtex_mimat_donor_matrix()
    mirbase_mat = _load_gtex_mirbase_donor_matrix()
    m2a = _mimat_to_arm()

    from mirna_hallmark.data_loaders import load_mirna_arms

    tcga_arms = list(load_mirna_arms().index)
    out: Dict[str, pd.Series] = {}
    audit_rows: List[dict] = []

    arm_to_mimats: Dict[str, List[str]] = {}
    for mimat in mimat_mat.index:
        arm = m2a.get(str(mimat))
        if not arm or not str(arm).startswith("hsa-"):
            continue
        arm_to_mimats.setdefault(str(arm), []).append(str(mimat))

    for arm, mimats in arm_to_mimats.items():
        sub = mimat_mat.loc[mimats]
        out[arm] = sub.mean(axis=0) if len(mimats) > 1 else sub.iloc[0]
        feat = mimats[0] if len(mimats) == 1 else "|".join(mimats)
        audit_rows.append({"tcga_arm": arm, "gtex_feature": feat, "gtex_arm_map_kind": "mimat"})

    mirbase_idx = set(mirbase_mat.index.astype(str))
    for arm in tcga_arms:
        if arm in out:
            continue
        hit, kind = _resolve_mirbase_feature(str(arm), mirbase_idx)
        if hit is None:
            audit_rows.append({"tcga_arm": arm, "gtex_feature": "", "gtex_arm_map_kind": "missing"})
            continue
        out[arm] = mirbase_mat.loc[hit]
        audit_rows.append({"tcga_arm": arm, "gtex_feature": hit, "gtex_arm_map_kind": kind})

    mat = pd.DataFrame(out).T if out else pd.DataFrame()
    _ARM_MAP = (
        pd.DataFrame(audit_rows).set_index("tcga_arm")
        if audit_rows
        else pd.DataFrame(columns=["gtex_feature", "gtex_arm_map_kind"])
    )
    return mat


def gtex_arm_map_for(tcga_arm: str) -> Tuple[str, str]:
    """Return (gtex_feature, kind) for one TCGA arm."""
    tbl = gtex_arm_map_table()
    if tcga_arm not in tbl.index:
        return "", "missing"
    row = tbl.loc[tcga_arm]
    return str(row.get("gtex_feature", "") or ""), str(row.get("gtex_arm_map_kind", "missing") or "missing")
