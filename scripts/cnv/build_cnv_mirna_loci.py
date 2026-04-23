#!/usr/bin/env python3
"""
Build CNV miRNA loci table with mature-arm mapping from MirGeneDB-style GFF.

Input:  PATHS.working_dir/miRNA/hsa.gff
Output: PATHS.mirna_path (default: data/miRNA/cnv_miRNA.csv)

The output is used by pipeline/CNV to annotate CNV segments with nearby/overlapping
pre-miRNA loci. We augment each hairpin locus with:
  - mature_names:      comma-separated mature arm IDs (e.g. Hsa-Mir-8-P1a_3p)
  - mature_accessions: comma-separated MIMAT accessions when present
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Tuple

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


def main() -> None:
    ap = argparse.ArgumentParser(description="Build CNV miRNA loci table (pre-miRNA + mature-arm mapping).")
    ap.add_argument("--gff", type=Path, default=None, help="Input GFF (default: data/miRNA/hsa.gff)")
    ap.add_argument("--out", type=Path, default=None, help="Output CSV (default: PATHS.mirna_path)")
    args = ap.parse_args()

    gff_path = args.gff or getattr(PATHS, "mirna_gff", (PATHS.working_dir / "miRNA" / "hsa.gff"))
    out_path = args.out or PATHS.mirna_path

    cols = ["chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=cols, low_memory=False)
    gff["attrs"] = gff["attributes"].apply(_parse_attrs)
    attrs = gff["attrs"].apply(pd.Series)
    gff = pd.concat([gff.drop(columns=["attributes", "attrs"]), attrs], axis=1)

    # pre-miRNA loci: ID (name) + Alias (MI...)
    pre = gff.loc[gff["type"] == "pre_miRNA"].copy()
    if pre.empty:
        raise SystemExit(f"No pre_miRNA rows found in {gff_path}")

    pre = pre.rename(
        columns={
            "chrom": "chrom",
            "start": "start",
            "end": "end",
            "strand": "strand",
            "ID": "gene_name",
            "Alias": "gene_id",
        }
    )
    pre = pre[["chrom", "start", "end", "strand", "gene_name", "gene_id"]].dropna()

    # Join key: mature rows use <base>_5p/_3p; pre rows use <base>_pre
    pre["pre_base"] = (
        pre["gene_name"]
        .astype(str)
        .str.replace(r"\r", "", regex=True)
        .str.strip()
        .str.replace(r"_pre$", "", regex=True)
    )

    # mature arms: ID (name) + Alias (MIMAT...) — map back to hairpin by trimming suffix.
    # In this GFF, mature IDs look like: <pre_id>_5p*, <pre_id>_3p
    mat = gff.loc[gff["type"] == "miRNA"].copy()
    mat["mature_name"] = mat.get("ID")
    mat["mature_accession"] = mat.get("Alias")
    mat = mat[["mature_name", "mature_accession"]].dropna(subset=["mature_name"])

    def pre_key(mature_id: str) -> str:
        s = str(mature_id)
        s = s.replace("_5p*", "").replace("_3p*", "").replace("_5p", "").replace("_3p", "")
        return s.replace("\r", "").strip()

    mat["pre_id"] = mat["mature_name"].map(pre_key)

    names_by_pre: DefaultDict[str, List[str]] = defaultdict(list)
    acc_by_pre: DefaultDict[str, List[str]] = defaultdict(list)
    for _, r in mat.iterrows():
        pid = str(r["pre_id"])
        mn = str(r["mature_name"])
        names_by_pre[pid].append(mn)
        ma = r.get("mature_accession")
        if isinstance(ma, str) and ma.strip():
            acc_by_pre[pid].append(str(ma).strip())

    pre["mature_names"] = pre["pre_base"].map(lambda x: ",".join(sorted(set(names_by_pre.get(str(x), [])))))
    pre["mature_accessions"] = pre["pre_base"].map(lambda x: ",".join(sorted(set(acc_by_pre.get(str(x), [])))))

    pre = pre.drop(columns=["pre_base"])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pre.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(pre)} pre-miRNA loci)")


if __name__ == "__main__":
    main()

