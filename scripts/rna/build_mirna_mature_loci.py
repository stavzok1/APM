#!/usr/bin/env python3
"""
Build a mature-arm miRNA loci table (MIMAT...) with genomic coordinates from PATHS.mirna_gff.

Output is used when we need *arm-resolution by coordinates*:
- an arm is reported as a hit only if a variant/SV/CNV overlaps (or is within window of) the arm interval itself.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402


def _parse_attrs(s: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in str(s).strip().split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def _load_mimat_to_mirbase_id_map() -> Dict[str, str]:
    """
    Map MIMAT accession -> miRBase mature ID (e.g. MIMAT0000062 -> hsa-let-7a-5p).

    Source: PATHS.mir_family_info (TargetScan family table), filtered to human (Species ID 9606).
    """
    p = getattr(PATHS, "mir_family_info", None)
    if p is None:
        return {}
    p = Path(p)
    if not p.exists():
        return {}
    try:
        df = pd.read_csv(p, sep="\t", low_memory=False)
    except Exception:
        df = pd.read_csv(p, low_memory=False)
    want = {"Species ID", "MiRBase ID", "MiRBase Accession"}
    if not want.issubset(set(df.columns)):
        return {}
    sub = df.loc[df["Species ID"].astype(str) == "9606", ["MiRBase Accession", "MiRBase ID"]].dropna()
    out: Dict[str, str] = {}
    for _, r in sub.iterrows():
        acc = str(r["MiRBase Accession"]).strip()
        mid = str(r["MiRBase ID"]).strip()
        if acc and mid and acc.upper().startswith("MIMAT"):
            out[acc.upper()] = mid
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Build mature-arm miRNA loci (MIMAT) table from hsa.gff.")
    ap.add_argument("--gff", type=Path, default=None, help="Input GFF (default: PATHS.mirna_gff)")
    ap.add_argument("--out", type=Path, default=None, help="Output CSV (default: PATHS.mirna_mature_loci_csv)")
    args = ap.parse_args()

    gff_path = args.gff or getattr(PATHS, "mirna_gff", (PATHS.working_dir / "miRNA" / "hsa.gff"))
    out_path = args.out or getattr(PATHS, "mirna_mature_loci_csv", (PATHS.working_dir / "miRNA" / "mirna_mature_loci.csv"))

    mimat_to_mirbase = _load_mimat_to_mirbase_id_map()

    cols = ["chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=cols, low_memory=False)
    gff["attrs"] = gff["attributes"].apply(_parse_attrs)
    attrs = gff["attrs"].apply(pd.Series)
    gff = pd.concat([gff.drop(columns=["attributes", "attrs"]), attrs], axis=1)

    pre = gff.loc[gff["type"] == "pre_miRNA", ["chrom", "start", "end", "strand", "ID", "Alias"]].copy()
    pre = pre.rename(columns={"ID": "pre_gene_name", "Alias": "pre_gene_id"})
    pre["pre_base"] = (
        pre["pre_gene_name"]
        .astype(str)
        .str.replace(r"\r", "", regex=True)
        .str.strip()
        .str.replace(r"_pre$", "", regex=True)
    )

    mat = gff.loc[gff["type"] == "miRNA", ["chrom", "start", "end", "strand", "ID", "Alias"]].copy()
    mat = mat.rename(columns={"ID": "mature_name", "Alias": "mature_accession"})
    mat["mature_base"] = (
        mat["mature_name"]
        .astype(str)
        .str.replace(r"\r", "", regex=True)
        .str.strip()
        .str.replace("_5p*", "", regex=False)
        .str.replace("_3p*", "", regex=False)
        .str.replace("_5p", "", regex=False)
        .str.replace("_3p", "", regex=False)
    )

    # Join mature → pre by base name (MirGeneDB-style IDs in this GFF)
    merged = mat.merge(
        pre[["pre_base", "pre_gene_name", "pre_gene_id"]],
        left_on="mature_base",
        right_on="pre_base",
        how="left",
    ).drop(columns=["pre_base"])

    out = merged.rename(
        columns={
            "chrom": "chrom",
            "start": "start",
            "end": "end",
            "strand": "strand",
        }
    )

    # Canonical columns used by overlap mappers
    out["gene_name"] = out["mature_name"]
    out["gene_id"] = out["mature_accession"]

    # Canonical miRBase mature names (arm-specific) for easy joins/reporting.
    # This is only populated when the GFF provides a MIMAT accession and it is present in miR_Family_Info.
    out["mirbase_mature_id"] = out["mature_accession"].apply(
        lambda x: mimat_to_mirbase.get(str(x).strip().upper(), "") if pd.notna(x) else ""
    )

    # Keep only rows with a mature accession when possible; but allow missing aliases.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} mature rows)")


if __name__ == "__main__":
    main()

