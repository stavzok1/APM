"""
miRTarBase experimentally-validated miRNA–gene interaction processing.

Loads miRTarBase data, normalizes support types and experiment classes,
collapses to study-level interactions, and produces four summary tables:
    - interaction_study:   (miRNA, gene, study) level with binary flags
    - interaction_summary: (miRNA, gene) level with study counts and evidence score
    - gene_summary:        (gene) level with per-miRNA detail dicts
    - mirna_summary:       (miRNA) level with per-gene detail dicts
    - family_summary:      (miRNA family) level counts

Optionally integrates TargetScan miR family annotations for seed-based
family assignment, falling back to a simple name-based proxy.
"""

import os
import re
import json
from pathlib import Path
from typing import List, Optional, Dict, Any

import numpy as np
import pandas as pd

from ..config import PATHS, PIPELINE_GENE_PANEL


# =============================================================================
# CONSTANTS
# =============================================================================

# Normalized support types used everywhere
SUPPORT_TYPES: List[str] = [
    "Functional MTI",
    "Functional MTI (Weak)",
    "Non-Functional MTI",
    "Non-Functional MTI (Weak)",
]

# Experiment classes for stratification
EXPERIMENT_CLASSES: List[str] = [
    "reporter",
    "binding",
    "protein",
    "rna",
    "perturbation",
    "other",
]


# =============================================================================
# COLUMN DETECTION
# =============================================================================

def _find_first_existing_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first column name from *candidates* that exists in *df*."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _infer_mirtarbase_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    """
    Infer miRTarBase column names from common header variants.

    Returns a dict mapping logical names to actual column names.
    Raises ValueError if required columns (mirna, gene, experiments,
    support_type) cannot be found.
    """
    colmap = {
        "mirtarbase_id": _find_first_existing_column(df, [
            "miRTarBase ID", "miRTarBase_ID", "MIRTarBase ID", "MIRTarBase_ID",
        ]),
        "mirna": _find_first_existing_column(df, [
            "miRNA", "miRNA Name", "Mature miRNA", "miRNA_name",
        ]),
        "mirna_species": _find_first_existing_column(df, [
            "Species (miRNA)", "miRNA Species", "Species_miRNA",
        ]),
        "gene": _find_first_existing_column(df, [
            "Target Gene", "Target gene", "Gene", "Target Symbol", "Target",
        ]),
        "gene_entrez": _find_first_existing_column(df, [
            "Target Gene (Entrez ID)", "Entrez ID", "Target Entrez ID",
        ]),
        "gene_species": _find_first_existing_column(df, [
            "Species (Target Gene)", "Target Gene Species", "Species_Target",
        ]),
        "experiments": _find_first_existing_column(df, [
            "Experiments", "Experiment", "Experimental Methods",
        ]),
        "support_type": _find_first_existing_column(df, [
            "Support Type", "SupportType",
        ]),
        "references": _find_first_existing_column(df, [
            "References", "Reference", "PMID", "PMIDs", "PubMed ID", "Source",
        ]),
    }

    required = ["mirna", "gene", "experiments", "support_type"]
    missing = [k for k in required if colmap[k] is None]
    if missing:
        raise ValueError(
            f"Could not infer required miRTarBase columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    return colmap


# =============================================================================
# NORMALIZATION HELPERS
# =============================================================================

def _normalize_support_type(x) -> str:
    """Normalize miRTarBase support-type strings to canonical form."""
    if pd.isna(x):
        return np.nan
    x_low = re.sub(r"\s+", " ", str(x).strip().lower().replace("_", " "))

    mapping = {
        "functional mti": "Functional MTI",
        "functional mti (weak)": "Functional MTI (Weak)",
        "non-functional mti": "Non-Functional MTI",
        "non-functional mti (weak)": "Non-Functional MTI (Weak)",
    }
    return mapping.get(x_low, str(x).strip())


def _normalize_gene_symbol(x) -> str:
    if pd.isna(x):
        return np.nan
    return str(x).strip()


def _normalize_mirna_name(x) -> str:
    if pd.isna(x):
        return np.nan
    return re.sub(r"\s+", "", str(x).strip())


def _simple_mirna_family(mirna: str) -> str:
    """
    Simple family proxy from miRNA name.

    hsa-miR-122-5p  -> miR-122
    hsa-let-7a-5p   -> let-7a

    NOT seed-family exact like TargetScan, but a usable fallback.
    """
    if pd.isna(mirna):
        return np.nan
    m = str(mirna).strip()
    m = re.sub(r"^[A-Za-z]{3}-", "", m)       # remove species prefix
    m = re.sub(r"-(5p|3p)$", "", m, flags=re.IGNORECASE)  # remove arm
    return m


def _split_experiments(exp_str) -> List[str]:
    """Split miRTarBase experiment strings (// ; , | delimited)."""
    if pd.isna(exp_str):
        return []
    parts = re.split(r"//|;|,|\|", str(exp_str).strip())
    return [p.strip() for p in parts if p.strip()]


def _classify_experiment(exp: str) -> str:
    """Collapse detailed experiment names into broad classes."""
    e = str(exp).strip().lower()

    reporter_kw = ["luciferase", "reporter assay", "reporter", "gfp reporter"]
    binding_kw = [
        "clip", "hits-clip", "par-clip", "clash", "rip", "rip-chip",
        "pull-down", "pulldown", "immunoprecipitation", "binding assay",
    ]
    protein_kw = [
        "western blot", "elisa", "flow cytometry", "immunoblot",
        "protein assay", "immunohistochemistry", "ihc",
    ]
    rna_kw = [
        "qrt-pcr", "qpcr", "rt-pcr", "northern blot", "microarray",
        "rna-seq", "sequencing", "expression profiling", "pcr",
    ]
    perturbation_kw = [
        "overexpression", "knockdown", "transfection", "sirna", "shrna",
        "mimic", "inhibitor", "antisense", "gene silencing",
    ]

    if any(k in e for k in reporter_kw):
        return "reporter"
    if any(k in e for k in binding_kw):
        return "binding"
    if any(k in e for k in protein_kw):
        return "protein"
    if any(k in e for k in rna_kw):
        return "rna"
    if any(k in e for k in perturbation_kw):
        return "perturbation"
    return "other"


def _extract_pmids(ref_value) -> List[str]:
    """Extract PMIDs (5–9 digit numbers) from a reference string."""
    if pd.isna(ref_value):
        return []
    s = str(ref_value)
    pmids = list(dict.fromkeys(re.findall(r"\b\d{5,9}\b", s)))
    if pmids:
        return pmids
    s = s.strip()
    return [f"REF::{s}"] if s else []


def _support_slug(support_type: str) -> str:
    """Convert support type to a column-name-safe slug."""
    return (
        support_type.lower()
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("-", "")
    )


def _dict_to_json(d: dict) -> str:
    return json.dumps(d, ensure_ascii=False, sort_keys=True)


# =============================================================================
# LOAD & NORMALIZE
# =============================================================================

def load_mirtarbase(
    mirtarbase_csv: Optional[Path] = None,
    family_info_tsv: Optional[Path] = None,
    gene_panel: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load miRTarBase CSV, apply TargetScan family annotations, normalize,
    and optionally filter to a gene panel.

    Parameters
    ----------
    mirtarbase_csv : Path, optional
        Path to miRTarBase CSV. Defaults to PATHS.mirtarbase_csv.
    family_info_tsv : Path, optional
        Path to TargetScan miR_Family_Info.txt. Defaults to PATHS.mir_family_info.
    gene_panel : list[str], optional
        Gene symbols to retain. Defaults to ``PIPELINE_GENE_PANEL``.

    Returns
    -------
    pd.DataFrame
        Normalized, filtered miRTarBase table with columns:
        miRNA_norm, gene_norm, support_type_norm, miRNA_family,
        experiment_list, experiment_classes, study_ids.
    """
    mirtarbase_csv = Path(mirtarbase_csv or PATHS.mirtarbase_csv)
    family_info_tsv = Path(family_info_tsv or PATHS.mir_family_info)
    gene_panel = gene_panel or PIPELINE_GENE_PANEL

    print(f"Loading miRTarBase from: {mirtarbase_csv}")
    mirs = pd.read_csv(mirtarbase_csv)
    cols = _infer_mirtarbase_columns(mirs)
    print("  Detected columns:")
    for k, v in cols.items():
        print(f"    {k}: {v}")

    # --- TargetScan family lookup ---
    mirna_to_family: Dict[str, str] = {}
    if family_info_tsv.exists():
        print(f"Loading miR family info from: {family_info_tsv}")
        families = pd.read_csv(family_info_tsv, sep="\t")
        families = families[families["Species ID"] == 9606].copy()
        families = families.rename(columns={
            "miR family": "mir_family",
            "MiRBase ID": "mirna",
        })
        families["mirna_norm"] = families["mirna"].map(_normalize_mirna_name)
        mirna_to_family = (
            families
            .drop_duplicates("mirna_norm")
            .set_index("mirna_norm")["mir_family"]
            .to_dict()
        )
        print(f"  Loaded {len(mirna_to_family):,} miRNA → family mappings")
    else:
        print(f"  miR family file not found, using name-based proxy only")

    # --- Normalize ---
    df = mirs.copy()
    df["miRNA_norm"] = df[cols["mirna"]].map(_normalize_mirna_name)
    df["gene_norm"] = df[cols["gene"]].map(_normalize_gene_symbol)

    import os

    if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    ):
        from .symbol_normalization import apply_symbol_mapping_series, default_symbol_mapping

        df["gene_norm"] = apply_symbol_mapping_series(
            df["gene_norm"], default_symbol_mapping()
        )

    # Filter to gene panel
    if gene_panel is not None:
        gene_set = set(gene_panel)
        before = len(df)
        df = df[df["gene_norm"].isin(gene_set)].copy()
        print(f"  Filtered to gene panel: {before:,} → {len(df):,} rows "
              f"({df['gene_norm'].nunique():,} genes)")

    df["support_type_norm"] = df[cols["support_type"]].map(_normalize_support_type)

    # Family assignment: TargetScan first, simple proxy fallback
    df["miRNA_family_targetscan"] = df["miRNA_norm"].map(mirna_to_family)
    df["miRNA_family"] = df["miRNA_family_targetscan"].fillna(
        df["miRNA_norm"].map(_simple_mirna_family)
    )

    df["experiment_list"] = df[cols["experiments"]].map(_split_experiments)
    df["experiment_classes"] = df["experiment_list"].map(
        lambda xs: sorted(set(_classify_experiment(x) for x in xs)) if xs else []
    )

    # Study IDs from references/PMID column
    if cols["references"] is not None:
        df["study_ids"] = df[cols["references"]].map(_extract_pmids)
    else:
        fallback_col = cols["mirtarbase_id"]
        if fallback_col is None:
            raise ValueError(
                "No references/PMID column and no miRTarBase ID fallback found."
            )
        df["study_ids"] = df[fallback_col].astype(str).map(lambda x: [f"MIRTAR::{x}"])

    # Drop incomplete rows
    df = df.dropna(subset=["miRNA_norm", "gene_norm", "support_type_norm"]).copy()
    df = df[df["support_type_norm"].isin(SUPPORT_TYPES)].copy()

    print(f"  Normalized rows retained: {len(df):,}")
    return df


# =============================================================================
# EXPLODE TO STUDY LEVEL
# =============================================================================

def _explode_to_study_level(df: pd.DataFrame) -> pd.DataFrame:
    """Explode study_ids so each row is one (miRNA, gene, study) record."""
    df_study = df.explode("study_ids").rename(columns={"study_ids": "study_id"}).copy()
    df_study["study_id"] = df_study["study_id"].astype(str).str.strip()
    df_study = df_study[df_study["study_id"] != ""].copy()
    print(f"  Study-exploded rows: {len(df_study):,}")
    return df_study


# =============================================================================
# COLLAPSE TO (miRNA, gene, study) WITH BINARY FLAGS
# =============================================================================

def _collapse_interaction_study(group: pd.DataFrame) -> pd.Series:
    """One row per (miRNA, gene, study_id) with support/experiment flags."""
    out = {
        "miRNA": group["miRNA_norm"].iloc[0],
        "miRNA_family": group["miRNA_family"].iloc[0],
        "gene": group["gene_norm"].iloc[0],
        "study_id": group["study_id"].iloc[0],
    }

    support_set = set(group["support_type_norm"].dropna())
    for st in SUPPORT_TYPES:
        out[f"has_{_support_slug(st)}"] = int(st in support_set)

    exp_classes = set()
    for xs in group["experiment_classes"]:
        exp_classes.update(xs)
    for ex in EXPERIMENT_CLASSES:
        out[f"has_{ex}"] = int(ex in exp_classes)

    # experiment × support cross flags
    for ex in EXPERIMENT_CLASSES:
        for st in SUPPORT_TYPES:
            mask = (
                group["support_type_norm"].eq(st)
                & group["experiment_classes"].map(lambda xs, _ex=ex: _ex in xs)
            )
            out[f"has_{ex}__{_support_slug(st)}"] = int(mask.any())

    return pd.Series(out)


def build_interaction_study_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build the (miRNA, gene, study) collapsed table with binary flags.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_mirtarbase().

    Returns
    -------
    pd.DataFrame
        One row per (miRNA, gene, study) with has_* flag columns.
    """
    df_study = _explode_to_study_level(df)

    interaction_study = (
        df_study
        .groupby(["miRNA_norm", "gene_norm", "study_id"], dropna=False)
        .apply(_collapse_interaction_study)
        .reset_index(drop=True)
    )
    print(f"  Interaction-study rows: {len(interaction_study):,}")
    return interaction_study


# =============================================================================
# TABLE 1: INTERACTION SUMMARY  (miRNA × gene)
# =============================================================================

def _summarize_interaction(group: pd.DataFrame) -> pd.Series:
    out = {
        "miRNA": group["miRNA"].iloc[0],
        "miRNA_family": group["miRNA_family"].iloc[0],
        "gene": group["gene"].iloc[0],
        "n_studies": group["study_id"].nunique(),
    }

    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        out[f"n_{slug}_studies"] = int(group[f"has_{slug}"].sum())

    for ex in EXPERIMENT_CLASSES:
        out[f"n_{ex}_studies"] = int(group[f"has_{ex}"].sum())

    exp_support_nested = {}
    for ex in EXPERIMENT_CLASSES:
        exp_support_nested[ex] = {}
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            val = int(group[f"has_{ex}__{slug}"].sum())
            out[f"n_{ex}__{slug}_studies"] = val
            exp_support_nested[ex][st] = val
    out["experiment_support_counts_json"] = _dict_to_json(exp_support_nested)

    # evidence score (tweak weights as needed)
    out["evidence_score"] = (
        3 * int(group["has_reporter"].sum())
        + 3 * int(group["has_binding"].sum())
        + 2 * int(group["has_protein"].sum())
        + 1 * int(group["has_rna"].sum())
        + 1 * int(group["has_perturbation"].sum())
    )

    return pd.Series(out)


def build_interaction_summary(interaction_study: pd.DataFrame) -> pd.DataFrame:
    """Row = (miRNA, gene). Study counts and evidence score."""
    result = (
        interaction_study
        .groupby(["miRNA", "gene"], dropna=False)
        .apply(_summarize_interaction)
        .reset_index(drop=True)
    )
    print(f"  Interaction summary rows: {len(result):,}")
    return result


# =============================================================================
# TABLE 2: GENE SUMMARY  (row = gene)
# =============================================================================

def _summarize_gene(group: pd.DataFrame) -> pd.Series:
    out = {
        "gene": group["gene"].iloc[0],
        "n_unique_miRNAs": group["miRNA"].nunique(),
        "n_unique_families": group["miRNA_family"].nunique(),
        "n_total_studies": int(group["study_id"].nunique()),
    }

    # --- study counts per support type and experiment class ---
    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        out[f"n_{slug}_studies"] = int(group[f"has_{slug}"].sum())

    for ex in EXPERIMENT_CLASSES:
        out[f"n_{ex}_studies"] = int(group[f"has_{ex}"].sum())

    exp_support_nested = {}
    for ex in EXPERIMENT_CLASSES:
        exp_support_nested[ex] = {}
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            val = int(group[f"has_{ex}__{slug}"].sum())
            out[f"n_{ex}__{slug}_studies"] = val
            exp_support_nested[ex][st] = val
    out["experiment_support_counts_json"] = _dict_to_json(exp_support_nested)

    # --- per-miRNA / per-family study count dicts ---
    out["mirna_study_counts_json"] = _dict_to_json(
        group.groupby("miRNA")["study_id"].nunique()
        .sort_values(ascending=False).to_dict()
    )
    out["family_study_counts_json"] = _dict_to_json(
        group.groupby("miRNA_family")["study_id"].nunique()
        .sort_values(ascending=False).to_dict()
    )

    # --- membership dicts: support_type → list of miRNAs ---
    support_to_mirnas = {}
    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        support_to_mirnas[st] = sorted(
            group.loc[group[f"has_{slug}"] == 1, "miRNA"].dropna().unique().tolist()
        )
    out["support_to_mirnas_json"] = _dict_to_json(support_to_mirnas)

    experiment_to_mirnas = {}
    for ex in EXPERIMENT_CLASSES:
        experiment_to_mirnas[ex] = sorted(
            group.loc[group[f"has_{ex}"] == 1, "miRNA"].dropna().unique().tolist()
        )
    out["experiment_to_mirnas_json"] = _dict_to_json(experiment_to_mirnas)

    exp_support_to_mirnas = {}
    for ex in EXPERIMENT_CLASSES:
        exp_support_to_mirnas[ex] = {}
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            exp_support_to_mirnas[ex][st] = sorted(
                group.loc[group[f"has_{ex}__{slug}"] == 1, "miRNA"]
                .dropna().unique().tolist()
            )
    out["experiment_support_to_mirnas_json"] = _dict_to_json(exp_support_to_mirnas)

    # --- same for families ---
    support_to_families = {}
    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        support_to_families[st] = sorted(
            group.loc[group[f"has_{slug}"] == 1, "miRNA_family"]
            .dropna().unique().tolist()
        )
    out["support_to_families_json"] = _dict_to_json(support_to_families)

    experiment_to_families = {}
    for ex in EXPERIMENT_CLASSES:
        experiment_to_families[ex] = sorted(
            group.loc[group[f"has_{ex}"] == 1, "miRNA_family"]
            .dropna().unique().tolist()
        )
    out["experiment_to_families_json"] = _dict_to_json(experiment_to_families)

    exp_support_to_families = {}
    for ex in EXPERIMENT_CLASSES:
        exp_support_to_families[ex] = {}
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            exp_support_to_families[ex][st] = sorted(
                group.loc[group[f"has_{ex}__{slug}"] == 1, "miRNA_family"]
                .dropna().unique().tolist()
            )
    out["experiment_support_to_families_json"] = _dict_to_json(exp_support_to_families)

    # --- detailed per-miRNA dict ---
    mirna_detail = {}
    for mir, sub in group.groupby("miRNA"):
        entry = {
            "n_studies": int(sub["study_id"].nunique()),
            "support_counts": {},
            "experiment_support_counts": {},
        }
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            entry["support_counts"][st] = int(sub[f"has_{slug}"].sum())
        for ex in EXPERIMENT_CLASSES:
            entry["experiment_support_counts"][ex] = {}
            for st in SUPPORT_TYPES:
                slug = _support_slug(st)
                entry["experiment_support_counts"][ex][st] = int(
                    sub[f"has_{ex}__{slug}"].sum()
                )
        mirna_detail[mir] = entry
    out["mirna_detail_json"] = _dict_to_json(mirna_detail)

    return pd.Series(out)


def build_gene_summary(interaction_study: pd.DataFrame) -> pd.DataFrame:
    """Row = gene. Aggregated miRNA evidence with nested detail dicts."""
    result = (
        interaction_study
        .groupby("gene", dropna=False)
        .apply(_summarize_gene)
        .reset_index(drop=True)
    )
    print(f"  Gene summary rows: {len(result):,}")
    return result


# =============================================================================
# TABLE 3: miRNA SUMMARY  (row = miRNA)
# =============================================================================

def _summarize_mirna(group: pd.DataFrame) -> pd.Series:
    out = {
        "miRNA": group["miRNA"].iloc[0],
        "miRNA_family": group["miRNA_family"].iloc[0],
        "n_unique_targets": group["gene"].nunique(),
        "n_total_studies": int(group["study_id"].nunique()),
    }

    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        out[f"n_{slug}_studies"] = int(group[f"has_{slug}"].sum())

    for ex in EXPERIMENT_CLASSES:
        out[f"n_{ex}_studies"] = int(group[f"has_{ex}"].sum())

    exp_support_nested = {}
    for ex in EXPERIMENT_CLASSES:
        exp_support_nested[ex] = {}
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            val = int(group[f"has_{ex}__{slug}"].sum())
            out[f"n_{ex}__{slug}_studies"] = val
            exp_support_nested[ex][st] = val
    out["experiment_support_counts_json"] = _dict_to_json(exp_support_nested)

    out["target_gene_study_counts_json"] = _dict_to_json(
        group.groupby("gene")["study_id"].nunique()
        .sort_values(ascending=False).to_dict()
    )

    # membership dicts: support/experiment → list of genes
    support_to_genes = {}
    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        support_to_genes[st] = sorted(
            group.loc[group[f"has_{slug}"] == 1, "gene"].dropna().unique().tolist()
        )
    out["support_to_genes_json"] = _dict_to_json(support_to_genes)

    experiment_to_genes = {}
    for ex in EXPERIMENT_CLASSES:
        experiment_to_genes[ex] = sorted(
            group.loc[group[f"has_{ex}"] == 1, "gene"].dropna().unique().tolist()
        )
    out["experiment_to_genes_json"] = _dict_to_json(experiment_to_genes)

    exp_support_to_genes = {}
    for ex in EXPERIMENT_CLASSES:
        exp_support_to_genes[ex] = {}
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            exp_support_to_genes[ex][st] = sorted(
                group.loc[group[f"has_{ex}__{slug}"] == 1, "gene"]
                .dropna().unique().tolist()
            )
    out["experiment_support_to_genes_json"] = _dict_to_json(exp_support_to_genes)

    # detailed per-gene dict
    target_detail = {}
    for gene, sub in group.groupby("gene"):
        entry = {
            "n_studies": int(sub["study_id"].nunique()),
            "support_counts": {},
            "experiment_support_counts": {},
        }
        for st in SUPPORT_TYPES:
            slug = _support_slug(st)
            entry["support_counts"][st] = int(sub[f"has_{slug}"].sum())
        for ex in EXPERIMENT_CLASSES:
            entry["experiment_support_counts"][ex] = {}
            for st in SUPPORT_TYPES:
                slug = _support_slug(st)
                entry["experiment_support_counts"][ex][st] = int(
                    sub[f"has_{ex}__{slug}"].sum()
                )
        target_detail[gene] = entry
    out["target_gene_detail_json"] = _dict_to_json(target_detail)

    return pd.Series(out)


def build_mirna_summary(interaction_study: pd.DataFrame) -> pd.DataFrame:
    """Row = miRNA. Aggregated gene-target evidence with nested detail dicts."""
    result = (
        interaction_study
        .groupby("miRNA", dropna=False)
        .apply(_summarize_mirna)
        .reset_index(drop=True)
    )
    print(f"  miRNA summary rows: {len(result):,}")
    return result


# =============================================================================
# TABLE 4: FAMILY SUMMARY  (row = miRNA family)
# =============================================================================

def _summarize_family(group: pd.DataFrame) -> pd.Series:
    out = {
        "miRNA_family": group["miRNA_family"].iloc[0],
        "n_unique_miRNAs": group["miRNA"].nunique(),
        "n_unique_targets": group["gene"].nunique(),
        "n_total_studies": int(group["study_id"].nunique()),
    }
    for st in SUPPORT_TYPES:
        slug = _support_slug(st)
        out[f"n_{slug}_studies"] = int(group[f"has_{slug}"].sum())
    for ex in EXPERIMENT_CLASSES:
        out[f"n_{ex}_studies"] = int(group[f"has_{ex}"].sum())
    return pd.Series(out)


def build_family_summary(interaction_study: pd.DataFrame) -> pd.DataFrame:
    """Row = miRNA family. Simple counts."""
    result = (
        interaction_study
        .groupby("miRNA_family", dropna=False)
        .apply(_summarize_family)
        .reset_index(drop=True)
    )
    print(f"  Family summary rows: {len(result):,}")
    return result


# =============================================================================
# SAVE
# =============================================================================

def _save_tables(
    interaction_study: pd.DataFrame,
    interaction_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    mirna_summary: pd.DataFrame,
    family_summary: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Write all five tables to *output_dir* as CSV."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "mirtar_interaction_study_collapsed.csv": interaction_study,
        "mirtar_interaction_summary.csv": interaction_summary,
        "mirtar_gene_summary.csv": gene_summary,
        "mirtar_mirna_summary.csv": mirna_summary,
        "mirtar_family_summary.csv": family_summary,
    }
    print(f"\n  Saving miRTarBase tables to: {output_dir}")
    for fname, table in files.items():
        path = output_dir / fname
        table.to_csv(path, index=False)
        print(f"    {fname}  ({len(table):,} rows)")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def get_mirtarbase_targets(
    mirtarbase_csv: Optional[Path] = None,
    family_info_tsv: Optional[Path] = None,
    gene_panel: Optional[List[str]] = None,
    output_dir: Optional[Path] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Run the full miRTarBase processing pipeline.

    Loads data, normalizes, collapses to study level, builds all summary
    tables, and saves to *output_dir*.

    Parameters
    ----------
    mirtarbase_csv : Path, optional
        Defaults to PATHS.mirtarbase_csv.
    family_info_tsv : Path, optional
        Defaults to PATHS.mir_family_info.
    gene_panel : list[str], optional
        Defaults to ``PIPELINE_GENE_PANEL``.
    output_dir : Path, optional
        Defaults to PATHS.working_dir / OUTPUT_SUBDIRS["mirna"] / "mirtarbase".

    Returns
    -------
    dict
        Keys: "interaction_study", "interaction_summary", "gene_summary",
              "mirna_summary", "family_summary".
    """
    from ..config import OUTPUT_SUBDIRS

    output_dir = Path(
        output_dir
        or PATHS.working_dir / OUTPUT_SUBDIRS["mirna"] / "mirtarbase"
    )

    # Load & normalize
    df = load_mirtarbase(
        mirtarbase_csv=mirtarbase_csv,
        family_info_tsv=family_info_tsv,
        gene_panel=gene_panel,
    )

    # Build tables
    interaction_study = build_interaction_study_table(df)
    interaction_summary = build_interaction_summary(interaction_study)
    gene_summary = build_gene_summary(interaction_study)
    mirna_summary = build_mirna_summary(interaction_study)
    family_summary = build_family_summary(interaction_study)

    # Save
    _save_tables(
        interaction_study,
        interaction_summary,
        gene_summary,
        mirna_summary,
        family_summary,
        output_dir,
    )

    return {
        "interaction_study": interaction_study,
        "interaction_summary": interaction_summary,
        "gene_summary": gene_summary,
        "mirna_summary": mirna_summary,
        "family_summary": family_summary,
    }
