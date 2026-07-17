"""Shared data loaders for the miRNA x Hallmark subproject.

Thin wrappers that reuse the canonical APM loaders (so we do not duplicate
parsing logic) and add the subproject-specific joins:

- Hallmark miRNA->target edges (built by ``build_edges.py``)
- clinical strata (PAM50 / TNBC / immune / stage) + Tier-A confounders
- miRNA arm x participant log2(RPM+1)         (reuses mirna_target_integration)
- RNA gene x participant log2(TPM+1)          (reuses utils.common.loaders)
- AGO / RISC expression subset
- CNV (ASCAT3 integer copy number) per target gene x participant, cached

All matrices are keyed on the 12-char ``participant`` and restricted to primary
tumors upstream; multi-vial collapse uses mean (expression) / median (CNV).
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C


# --------------------------------------------------------------------------- #
# Hallmark edges
# --------------------------------------------------------------------------- #
_EDGES_RAW: dict = {}    # cache the parsed TSV per path — the read is ~0.7s and is on every assemble_gene


def load_hallmark_edges(path: Optional[Path] = None) -> pd.DataFrame:
    """Load the (miRNA, gene, hallmark_set) edge table from build_edges.py. The parsed TSV is cached per
    path (deterministic); a fresh copy is returned so callers may mutate freely."""
    p = Path(path or C.HALLMARK_EDGES_TSV)
    key = str(p)
    if key not in _EDGES_RAW:
        if not p.exists():
            raise FileNotFoundError(
                f"Hallmark edges missing: {p}. Build with "
                f".venv/bin/python3 -m mirna_hallmark.build_edges"
            )
        _EDGES_RAW[key] = pd.read_csv(p, sep="\t", low_memory=False)
    return _EDGES_RAW[key].copy()


def high_evidence_edges(edges: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Subset to high-evidence rows (Functional MTI + low-throughput validation)."""
    edges = edges if edges is not None else load_hallmark_edges()
    return edges.loc[edges["high_evidence"].astype(bool)].copy()


# --------------------------------------------------------------------------- #
# Clinical strata + confounders
# --------------------------------------------------------------------------- #
def load_clinical_strata() -> pd.DataFrame:
    """Participant-keyed clinical frame with strata + Tier-A confounders.

    Columns (when available): participant, PAM50_final, tnbc_subtype_4,
    thornsson_immune_subtype, pathologic_stage_collapsed, CPE, thornsson_hrd_score.
    """
    from analysis.utils.common.loaders import load_clinical

    clin = load_clinical(C.CLINICAL_UNIFIED)
    clin = clin.drop_duplicates("participant").copy()

    # Join TNBC (Lehmann TNBCtype-4) subtype on participant
    if C.TNBC_SUBTYPES.exists():
        tnbc = pd.read_csv(C.TNBC_SUBTYPES, sep="\t", low_memory=False)
        tnbc = tnbc[["participant", "tnbc_subtype_4"]].drop_duplicates("participant")
        clin = clin.merge(tnbc, on="participant", how="left")

    wanted = [
        "participant",
        "PAM50_final",
        "tnbc_subtype_4",
        "thornsson_immune_subtype",
        "pathologic_stage_collapsed",
        "CPE",
        "thornsson_hrd_score",
        # Broad genome-wide aneuploidy / CIN measures (Taylor et al. arm-level
        # aneuploidy score + fraction-genome-altered) — used as confounders to
        # separate locus-specific CN effects from the genome-wide aneuploidy tide.
        "thornsson_aneuploidy_score",
        "thornsson_fraction_altered",
    ]
    cols = [c for c in wanted if c in clin.columns]
    out = clin[cols].copy()
    # Global cohort switch: discard PAM50 Normal-like everywhere (config.EXCLUDE_NORMAL_LIKE).
    if getattr(C, "EXCLUDE_NORMAL_LIKE", False) and "PAM50_final" in out.columns:
        n0 = len(out)
        out = out.loc[out["PAM50_final"] != C.NORMAL_LIKE_LABEL].copy()
        print(f"[data_loaders] EXCLUDE_NORMAL_LIKE: dropped {n0 - len(out)} Normal-like participants "
              f"({len(out)} remain)")
    return out


def augment_tcga_batch(cov: pd.DataFrame, kind: Optional[str] = None) -> pd.DataFrame:
    """Append TCGA sequencing-batch (analyte-plate) dummies to a **participant-indexed**
    covariate frame, per ``config.TCGA_BATCH_KIND`` (no-op when kind=='none').

    Centralizes the batch-correction wiring so every confounder-adjusted coupling test
    adds batch identically (covariate adjustment, not ComBat). The cov frame must be
    indexed by the 12-char participant id.
    """
    k = C.TCGA_BATCH_KIND if kind is None else kind
    if not k or k == "none":
        return cov
    from mirna_hallmark.tcga_batch import augment_cov as _augment_cov

    out, _added = _augment_cov(cov, k)
    return out


def normal_like_participants() -> set:
    """Set of 12-char participants whose PAM50 call is Normal-like (read from the RAW
    clinical, i.e. *before* the EXCLUDE_NORMAL_LIKE filter). Used to drop Normal-like
    primary-tumour columns from the barcode-keyed cross-state matrices (which split by
    sample type, not PAM50, so they would otherwise retain them)."""
    from analysis.utils.common.loaders import load_clinical

    clin = load_clinical(C.CLINICAL_UNIFIED).drop_duplicates("participant")
    if "PAM50_final" not in clin.columns:
        return set()
    return set(clin.loc[clin["PAM50_final"] == C.NORMAL_LIKE_LABEL, "participant"])


# --------------------------------------------------------------------------- #
# Expression matrices
# --------------------------------------------------------------------------- #
def load_mirna_arms() -> pd.DataFrame:
    """arm_name x participant log2(RPM+1) (reuses mirna_target_integration loaders)."""
    from analysis.expression.mirna_target_integration import (
        load_mimat_to_arm,
        load_mirna_expression,
    )

    mimat_to_arm = load_mimat_to_arm(C.MIRNA_MATURE_LOCI)
    return load_mirna_expression(C.MIRNA_EXPRESSION, mimat_to_arm=mimat_to_arm)


def load_rna() -> pd.DataFrame:
    """gene_symbol x participant log2(TPM+1) (reuses utils.common.loaders)."""
    from analysis.utils.common.loaders import load_rna_log2tpm

    return load_rna_log2tpm(C.RNA_EXPRESSION)


def load_ago_expression(
    rna: Optional[pd.DataFrame] = None,
    *,
    include_tnrc6: bool = False,
) -> pd.DataFrame:
    """Subset RNA log2(TPM+1) to AGO (+ optional TNRC6) genes x participant."""
    rna = rna if rna is not None else load_rna()
    genes: List[str] = list(C.AGO_GENES) + (list(C.RISC_EXTRA_GENES) if include_tnrc6 else [])
    present = [g for g in genes if g in rna.index]
    missing = [g for g in genes if g not in rna.index]
    if missing:
        print(f"[data_loaders] AGO/RISC genes absent from RNA matrix: {missing}")
    return rna.loc[present].copy()


# --------------------------------------------------------------------------- #
# CNV (target genes)
# --------------------------------------------------------------------------- #
def _combined_ascat3_path() -> Path:
    return C.CNV_GENE_TABLES_DIR / "cnv_gene_calls_all_samples_ascat3.csv"


def load_cnv_target_genes(
    genes: Sequence[str],
    *,
    cache: bool = True,
    cache_path: Optional[Path] = None,
    chunksize: int = 1_000_000,
) -> pd.DataFrame:
    """gene_name x participant integer copy number (ASCAT3), median across vials.

    The combined ``cnv_gene_calls_all_samples_ascat3.csv`` is large, so we cache
    the gene-subset matrix under ``output/matrices/`` and reuse it.
    """
    cache_path = Path(cache_path or (C.OUTPUT_ROOT / "matrices" / "cnv_target_genes.tsv.gz"))
    gene_set = set(genes)

    if cache and cache_path.exists():
        cached = pd.read_csv(cache_path, sep="\t", index_col=0)
        if gene_set.issubset(set(cached.index)):
            return cached.loc[sorted(gene_set & set(cached.index))]

    src = _combined_ascat3_path()
    if not src.exists():
        raise FileNotFoundError(f"ASCAT3 combined CNV table missing: {src}")

    parts: List[pd.DataFrame] = []
    usecols = ["gene_name", "copy_number", "participant"]
    for chunk in pd.read_csv(src, usecols=lambda c: c in usecols, chunksize=chunksize):
        sub = chunk.loc[chunk["gene_name"].isin(gene_set)].dropna(subset=["copy_number"])
        if not sub.empty:
            parts.append(sub)
    if not parts:
        raise ValueError("No CNV rows matched the requested genes")

    long = pd.concat(parts, ignore_index=True)
    long["copy_number"] = long["copy_number"].round().astype(int)
    mat = long.groupby(["gene_name", "participant"])["copy_number"].median().unstack()

    if cache:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        mat.to_csv(cache_path, sep="\t")
    return mat


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def common_participants(*frames: Iterable[str]) -> List[str]:
    """Intersection of participant id collections, sorted."""
    sets = [set(f) for f in frames if f is not None]
    if not sets:
        return []
    inter = set.intersection(*sets)
    return sorted(inter)


def classify_cn(series: pd.Series) -> pd.Series:
    """Integer CN -> {deep_del, loss, neutral, gain, amp} (matches CNV module)."""
    def _c(v: float) -> str:
        if pd.isna(v):
            return "NA"
        v = int(round(v))
        if v <= 0:
            return "deep_del"
        if v == 1:
            return "loss"
        if v == 2:
            return "neutral"
        if v == 3:
            return "gain"
        return "amp"

    return series.map(_c)
