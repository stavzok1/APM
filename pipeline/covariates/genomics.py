from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Sequence, Tuple

import pandas as pd

from pipeline.sample_ids import normalize_tcga_id


MAF_REQUIRED_COLS = ("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")


def _is_missense(vc: str) -> bool:
    return str(vc) in {"Missense_Mutation"}


def _is_truncating(vc: str) -> bool:
    return str(vc) in {
        "Nonsense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Nonstop_Mutation",
        "Splice_Site",
        "Translation_Start_Site",
    }


def _is_other_nonsyn(vc: str) -> bool:
    return str(vc) in {
        "In_Frame_Del",
        "In_Frame_Ins",
        "Splice_Region",
    }


def classify_tp53_from_maf(maf: pd.DataFrame) -> pd.DataFrame:
    """
    Rule-based TP53 status from a MAF-like table.

    Output is keyed by TCGA sample_vial (preferred) with columns:
      - TP53_has_mutation (bool)
      - TP53_mutation_class (WT/missense/truncating/other_nonsyn/unknown)
    """
    missing = [c for c in MAF_REQUIRED_COLS if c not in maf.columns]
    if missing:
        raise KeyError(f"MAF missing required columns: {missing}")

    df = maf.copy()
    df = df[df["Hugo_Symbol"].astype(str).str.upper() == "TP53"].copy()
    if df.empty:
        return pd.DataFrame(
            {
                "TP53_has_mutation": pd.Series(dtype="boolean"),
                "TP53_mutation_class": pd.Series(dtype="string"),
            }
        )

    # Normalize sample ids.
    raw = df["Tumor_Sample_Barcode"].astype(str).map(lambda x: normalize_tcga_id(x).raw)
    ids = df["Tumor_Sample_Barcode"].astype(str).map(normalize_tcga_id)
    df["sample_vial"] = ids.map(lambda x: x.sample_vial or x.sample or x.participant or x.raw)
    df["variant_class"] = df["Variant_Classification"].astype(str)

    def classify_group(g: pd.DataFrame) -> Tuple[bool, str]:
        vcs = list(g["variant_class"].dropna().astype(str))
        if any(_is_truncating(vc) for vc in vcs):
            return True, "truncating"
        if any(_is_missense(vc) for vc in vcs):
            return True, "missense"
        if any(_is_other_nonsyn(vc) for vc in vcs):
            return True, "other_nonsyn"
        return True, "unknown"

    rows = []
    for sid, grp in df.groupby("sample_vial", dropna=False):
        has_mut, cls = classify_group(grp)
        rows.append({"sample_vial": str(sid), "TP53_has_mutation": bool(has_mut), "TP53_mutation_class": cls})

    out = pd.DataFrame(rows).set_index("sample_vial")
    out.index.name = "sample_vial"
    return out


def flag_hotspot_mutations_from_maf(
    maf: pd.DataFrame,
    *,
    gene: str,
    aa_hotspots: Sequence[str],
    protein_change_col_candidates: Sequence[str] = ("HGVSp_Short", "Protein_Change", "HGVSp"),
) -> pd.DataFrame:
    """
    Generic hotspot flagger for MAF-like tables.

    Returns sample_vial-keyed booleans:
      - {gene}_any_mutation
      - {gene}_hotspot_mutation
    """
    missing = [c for c in MAF_REQUIRED_COLS if c not in maf.columns]
    if missing:
        raise KeyError(f"MAF missing required columns: {missing}")

    df = maf.copy()
    df = df[df["Hugo_Symbol"].astype(str).str.upper() == gene.upper()].copy()
    if df.empty:
        return pd.DataFrame(
            {
                f"{gene}_any_mutation": pd.Series(dtype="boolean"),
                f"{gene}_hotspot_mutation": pd.Series(dtype="boolean"),
            }
        )

    prot_col = next((c for c in protein_change_col_candidates if c in df.columns), None)
    if prot_col is None:
        df[prot_col] = ""
    ids = df["Tumor_Sample_Barcode"].astype(str).map(normalize_tcga_id)
    df["sample_vial"] = ids.map(lambda x: x.sample_vial or x.sample or x.participant or x.raw)
    df["prot_change"] = df[prot_col].astype(str)

    hs = set(str(x).strip() for x in aa_hotspots)

    rows = []
    for sid, grp in df.groupby("sample_vial", dropna=False):
        any_mut = True
        hotspot = any(str(x).strip() in hs for x in grp["prot_change"].dropna().astype(str))
        rows.append(
            {
                "sample_vial": str(sid),
                f"{gene}_any_mutation": bool(any_mut),
                f"{gene}_hotspot_mutation": bool(hotspot),
            }
        )

    out = pd.DataFrame(rows).set_index("sample_vial")
    out.index.name = "sample_vial"
    return out

