#!/usr/bin/env python3
"""
Map Ensembl gene IDs in a wide RNA expression matrix to gene symbols.

Expected input format (TSV):
- First column: Ensembl gene id (often versioned like ENSG... .5)
- Remaining columns: samples

Uses GENCODE probemap (id -> gene symbol) at:
  annotations/RNA/gencode.v36.annotation.gtf.gene.probemap

Writes:
- Mapped matrix (symbols as first column)
- A small mapping report (counts) to stdout
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def _strip_ens_version(s: str) -> str:
    s = (s or "").strip()
    if not s:
        return s
    # ENSG0000...(.version)
    return s.split(".", 1)[0]


def load_probemap(probemap_path: Path) -> Dict[str, str]:
    """
    Returns mapping from *unversioned* ENSG -> gene symbol.
    """
    m: Dict[str, str] = {}
    with probemap_path.open("r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if "id" not in r.fieldnames or "gene" not in r.fieldnames:
            raise ValueError(f"Unexpected probemap columns: {r.fieldnames}")
        for row in r:
            eid = _strip_ens_version(row["id"])
            sym = (row["gene"] or "").strip()
            if not eid or not sym:
                continue
            # Keep first; duplicates are usually identical symbols.
            m.setdefault(eid, sym)
    return m


def _agg_init(n: int) -> List[float]:
    return [0.0] * n


def _agg_apply(acc: List[float], row: List[float], *, mode: str, counts: List[int]) -> None:
    if mode in ("sum", "mean"):
        for i, v in enumerate(row):
            acc[i] += v
        counts[0] += 1
        return
    if mode == "max":
        for i, v in enumerate(row):
            if v > acc[i]:
                acc[i] = v
        counts[0] += 1
        return
    if mode == "first":
        if counts[0] == 0:
            for i, v in enumerate(row):
                acc[i] = v
        counts[0] += 1
        return
    raise ValueError(f"Unknown collapse mode: {mode}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-tsv", required=True, help="Wide RNA matrix with Ensembl IDs in first column")
    ap.add_argument("--out-tsv", required=True, help="Output TSV with gene symbols in first column")
    ap.add_argument(
        "--probemap",
        default=None,
        help="GENCODE probemap path. Default: pipeline.config.PATHS.rna_gencode_gene_probemap",
    )
    ap.add_argument(
        "--collapse",
        choices=["mean", "sum", "max", "first"],
        default="mean",
        help="How to collapse multiple Ensembl IDs mapping to same symbol (default: mean)",
    )
    ap.add_argument(
        "--keep-unmapped",
        action="store_true",
        help="If set, keep unmapped Ensembl IDs as their own row key (ENSG...)",
    )
    args = ap.parse_args()

    in_path = Path(args.in_tsv)
    out_path = Path(args.out_tsv)

    if args.probemap is None:
        try:
            # Allow running without installation
            root = Path(__file__).resolve().parents[2]
            if str(root) not in sys.path:
                sys.path.insert(0, str(root))
            from pipeline.config import PATHS  # noqa: E402

            probemap_path = Path(getattr(PATHS, "rna_gencode_gene_probemap"))
        except Exception:
            probemap_path = Path("annotations/RNA/gencode.v36.annotation.gtf.gene.probemap")
    else:
        probemap_path = Path(args.probemap)

    ens_to_sym = load_probemap(probemap_path)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_mapped = 0
    n_unmapped = 0
    n_collapsed_multi = 0

    # symbol -> (accumulator vector, [count])
    agg: Dict[str, Tuple[List[float], List[int]]] = {}

    with in_path.open("r", newline="") as fin:
        r = csv.reader(fin, delimiter="\t")
        header = next(r)
        if not header:
            raise SystemExit("Empty input")
        id_col = header[0]
        sample_cols = header[1:]
        n = len(sample_cols)

        for row in r:
            if not row:
                continue
            n_in += 1
            raw_id = row[0]
            eid = _strip_ens_version(raw_id)
            sym = ens_to_sym.get(eid)
            if sym:
                key = sym
                n_mapped += 1
            else:
                n_unmapped += 1
                if not args.keep_unmapped:
                    continue
                key = eid if eid else raw_id

            # parse numeric values, tolerate blanks
            vals = []
            for x in row[1 : 1 + n]:
                try:
                    vals.append(float(x))
                except Exception:
                    vals.append(0.0)

            if key not in agg:
                agg[key] = (_agg_init(n), [0])
            acc, cnt = agg[key]
            if cnt[0] >= 1:
                n_collapsed_multi += 1
            _agg_apply(acc, vals, mode=args.collapse, counts=cnt)

    # write
    with out_path.open("w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["gene"] + sample_cols)
        for key in sorted(agg.keys()):
            acc, cnt = agg[key]
            if args.collapse == "mean" and cnt[0] > 0:
                out_vals = [v / float(cnt[0]) for v in acc]
            else:
                out_vals = acc
            w.writerow([key] + [f"{v:.6g}" for v in out_vals])

    print(
        "[RNA map] rows_in=%d mapped=%d unmapped=%d kept_unmapped=%s collapsed_duplicates=%d out_rows=%d"
        % (n_in, n_mapped, n_unmapped, str(args.keep_unmapped), n_collapsed_multi, len(agg))
    )
    print(f"[RNA map] wrote -> {out_path}")


if __name__ == "__main__":
    main()

