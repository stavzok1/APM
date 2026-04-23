#!/usr/bin/env python3
"""
Compare miRNA IDs in TCGA miRNA expression matrix vs miRTarBase.

Outputs a small summary + overlap tables to an output directory.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, Set, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402


def _norm_mir(s: str) -> str:
    s = str(s)
    # miRTarBase sometimes has trailing tabs/spaces or other whitespace characters.
    # Remove *all* whitespace and normalize separators/case.
    s = re.sub(r"\s+", "", s, flags=re.UNICODE).strip()
    s = s.replace("_", "-")
    return s.lower()


def _strip_species(s: str) -> str:
    s = _norm_mir(s)
    return re.sub(r"^(hsa|mmu|rno|gga|dme)-", "", s, flags=re.IGNORECASE)


def _strip_arm(s: str) -> str:
    s = _strip_species(s)
    # remove -5p/-3p
    return re.sub(r"-(5p|3p)$", "", s, flags=re.IGNORECASE)


def _read_expression_ids(path: Path) -> pd.Series:
    # Many TCGA miRNA matrices have first col = miRNA ID
    df = pd.read_csv(path, sep="\t", low_memory=False, usecols=[0])
    col = df.columns[0]
    return df[col].dropna().astype(str)


def _read_mirtar_ids(path: Path) -> pd.Series:
    df = pd.read_csv(path, low_memory=False, usecols=["miRNA"])
    return df["miRNA"].dropna().astype(str)


def _maybe_load_mimat_to_name_map() -> dict:
    """
    Map miRBase Accession (MIMAT...) -> miRBase ID (e.g. hsa-let-7a-5p) for human.
    Uses PATHS.mir_family_info when available (seed family table).
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
    # Species ID 9606 = human
    sub = df.loc[df["Species ID"].astype(str) == "9606", ["MiRBase ID", "MiRBase Accession"]].dropna()
    out = {}
    for _, r in sub.iterrows():
        acc = _norm_mir(r["MiRBase Accession"])
        mid = _norm_mir(r["MiRBase ID"])
        if acc and mid and acc.startswith("mimat"):
            out[acc] = mid
    return out


def _looks_like_mimat(ids: Set[str]) -> bool:
    if not ids:
        return False
    n = len(ids)
    k = sum(1 for x in ids if re.fullmatch(r"mimat\d+", x))
    return (k / n) >= 0.8


def _map_mimat_to_names(ids: Set[str], mimat_map: dict) -> Set[str]:
    out: Set[str] = set()
    for x in ids:
        if re.fullmatch(r"mimat\d+", x):
            y = mimat_map.get(x)
            if y:
                out.add(y)
        else:
            out.add(x)
    return out


def _set(series: Iterable[str], fn) -> Set[str]:
    out: Set[str] = set()
    for x in series:
        y = fn(x)
        if y and y.lower() not in ("nan", "none"):
            out.add(y)
    return out


def _overlap(a: Set[str], b: Set[str]) -> Tuple[int, int, int]:
    return (len(a), len(b), len(a & b))


def main() -> None:
    ap = argparse.ArgumentParser(description="miRNA expression vs miRTarBase overlap.")
    ap.add_argument("--expr", type=Path, default=None, help="Expression TSV (default: PATHS.mirna_expression_tsv)")
    ap.add_argument("--mirtar", type=Path, default=None, help="miRTarBase CSV (default: PATHS.mirtarbase_csv)")
    ap.add_argument("--out", type=Path, required=True, help="Output directory")
    args = ap.parse_args()

    expr_path = args.expr or PATHS.mirna_expression_tsv
    mirtar_path = args.mirtar or PATHS.mirtarbase_csv
    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    expr_ids = _read_expression_ids(expr_path)
    mirtar_ids = _read_mirtar_ids(mirtar_path)

    expr_raw = _set(expr_ids, _norm_mir)
    mirtar_raw = _set(mirtar_ids, _norm_mir)

    mimat_map = _maybe_load_mimat_to_name_map()
    expr_is_mimat = _looks_like_mimat(expr_raw)
    if expr_is_mimat and mimat_map:
        expr_raw = _map_mimat_to_names(expr_raw, mimat_map)

    # Ensure derived views are consistent with any MIMAT->name mapping we applied.
    expr_no_species = {_strip_species(x) for x in expr_raw}
    mirtar_no_species = _set(mirtar_ids, _strip_species)

    expr_no_arm = {_strip_arm(x) for x in expr_raw}
    mirtar_no_arm = _set(mirtar_ids, _strip_arm)

    # quick “granularity” signal: do expression IDs carry arm suffix?
    n_arm = sum(1 for x in expr_raw if re.search(r"-(5p|3p)$", x, flags=re.IGNORECASE))
    arm_pct = 100.0 * n_arm / max(len(expr_raw), 1)

    rows = []
    for label, a, b in (
        ("raw (exact string)", expr_raw, mirtar_raw),
        ("strip species prefix", expr_no_species, mirtar_no_species),
        ("strip species + arm (-5p/-3p)", expr_no_arm, mirtar_no_arm),
    ):
        na, nb, ni = _overlap(a, b)
        rows.append(
            {
                "mode": label,
                "n_expr_unique": na,
                "n_mirtar_unique": nb,
                "n_overlap": ni,
                "pct_expr_covered_by_mirtar": round(100.0 * ni / max(na, 1), 3),
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "mirna_mirtarbase_overlap_summary.tsv", sep="\t", index=False)

    # dump explicit overlaps (best-effort: arm-stripped tends to match if expr lacks -5p/-3p)
    inter = sorted(expr_no_arm & mirtar_no_arm)
    pd.DataFrame({"mirna_id_norm": inter}).to_csv(out_dir / "mirna_overlap_ids.tsv", sep="\t", index=False)

    missing = sorted(expr_no_arm - mirtar_no_arm)
    pd.DataFrame({"mirna_id_norm": missing}).to_csv(out_dir / "mirna_missing_from_mirtarbase_ids.tsv", sep="\t", index=False)

    if expr_is_mimat:
        print("Detected expression IDs as MIMAT accessions; mapped to miRBase mature names using miR_Family_Info.txt.")
        if not mimat_map:
            print("Warning: could not load MIMAT->name map (PATHS.mir_family_info missing/unreadable).")
    print(f"Expression IDs with explicit -5p/-3p arm suffix: {n_arm}/{len(expr_raw)} ({arm_pct:.2f}%)")
    print(summary.to_string(index=False))
    print(f"Wrote: {out_dir}")


if __name__ == "__main__":
    main()

