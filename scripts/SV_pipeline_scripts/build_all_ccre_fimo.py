#!/usr/bin/env python3
"""
Build ``all_ccre_fimo.tsv`` by running MEME FIMO on every cCRE interval from ``GRCh38-cCREs.csv``.

This is **separate** from the SV per-sample FIMO steps (``05_fimo_tsv``). Output is formatted for
``pipeline.SV.motif_scanning.load_all_ccre_fimo``:

- ``sequence_name`` = chromosome **without** the ``chr`` prefix (e.g. ``1``, ``X``).
- ``start`` / ``stop`` = **genomic** coordinates in the same **1-based inclusive** convention as
  standard FIMO output (equivalently: positions on a 1-based chromosome model; see MEME FIMO docs).

Defaults match the SV motif leg (``run_sv_motif_scanning`` / ``run_sv_pipeline``):

- MEME file: ``PATHS.sv_meme_file``
- Reference FASTA: ``PATHS.sv_reference_fasta``
- FIMO ``--thresh``: ``THRESHOLDS.fimo_pvalue_threshold``
- cCRE table: ``PATHS.ccre_csv``

Requirements: ``bedtools`` and ``fimo`` on ``PATH``.

Example::

    cd /path/to/APM
    python scripts/SV_pipeline_scripts/build_all_ccre_fimo.py \\
        --out data/SV/build_ccre_fimo/all_ccre_fimo.tsv \\
        --chunk-size 3000

Override inputs::

    python scripts/SV_pipeline_scripts/build_all_ccre_fimo.py \\
        --ccre-csv data/GRCh38-cCREs.csv \\
        --meme data/SV/motifs/custom.meme \\
        --ref-fasta data/genome_assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
        --thresh 1e-4 \\
        --out /scratch/$USER/all_ccre_fimo.tsv
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd

from pipeline.config import PATHS, THRESHOLDS
from pipeline.utils import harmonize_chrom_column


def _require_bin(name: str) -> str:
    p = shutil.which(name)
    if not p:
        raise SystemExit(f"ERROR: `{name}` not found on PATH.")
    return p


def _strip_chr(chrom: str) -> str:
    s = str(chrom).strip()
    if s.lower().startswith("chr"):
        return s[3:]
    if s.lower().startswith("chrom"):
        return s[5:]
    return s


def _load_ccre(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df, _ = harmonize_chrom_column(df, "chrom")
    need = {"chrom", "start", "end"}
    if not need.issubset(df.columns):
        raise SystemExit(f"cCRE CSV must contain columns {sorted(need)}; got {list(df.columns)}")
    id_col = "cCRE_id" if "cCRE_id" in df.columns else "elem_id" if "elem_id" in df.columns else None
    if id_col is None:
        raise SystemExit("cCRE CSV must contain `cCRE_id` or `elem_id`.")
    if id_col != "cCRE_id":
        df = df.rename(columns={id_col: "cCRE_id"})
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["cCRE_id"] = df["cCRE_id"].astype(str)
    bad = (df["end"] <= df["start"]) | (df["start"] < 0)
    if bad.any():
        n = int(bad.sum())
        df = df.loc[~bad].copy()
        print(f"[WARN] Dropped {n} cCRE rows with invalid coordinates.", flush=True)
    return df


def _run_getfasta(bedtools: str, ref: Path, bed: Path, fasta: Path) -> None:
    cmd = [
        bedtools,
        "getfasta",
        "-fi",
        str(ref),
        "-bed",
        str(bed),
        "-name",
        "-fo",
        str(fasta),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def _run_fimo(fimo: str, thresh: float, meme: Path, fasta: Path, oc_dir: Path) -> Path:
    oc_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        fimo,
        "--oc",
        str(oc_dir),
        "--thresh",
        str(thresh),
        str(meme),
        str(fasta),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    out = oc_dir / "fimo.tsv"
    if not out.exists():
        raise RuntimeError(f"FIMO did not write {out}")
    return out


def _read_fimo(path: Path) -> pd.DataFrame:
    # Reuse robust parser from motif_scanning (handles MEME 5.5+ headers).
    from pipeline.SV.motif_scanning import _load_fimo_tsv_as_dataframe

    return _load_fimo_tsv_as_dataframe(path)


def _remap_chunk(
    fimo_df: pd.DataFrame,
    bounds: Dict[str, Tuple[str, int, int]],
) -> pd.DataFrame:
    """Map per-cCRE FASTA hits to genomic coordinates + chrom as sequence_name."""
    if fimo_df.empty:
        return fimo_df
    if "sequence_name" not in fimo_df.columns:
        raise ValueError("FIMO output missing sequence_name")
    rows: List[dict] = []
    for _, r in fimo_df.iterrows():
        sid = str(r["sequence_name"]).strip()
        if sid not in bounds:
            continue
        chrom, rs, re = bounds[sid]
        try:
            fs = int(pd.to_numeric(r["start"], errors="raise"))
            fe = int(pd.to_numeric(r["stop"], errors="raise"))
        except Exception:
            continue
        # FIMO coordinates are 1-based inclusive on the extracted sequence.
        # BED interval is [rs, re) 0-based half-open in ``bedtools getfasta`` input.
        g_start = rs + fs
        g_stop = rs + fe
        if g_stop < g_start:
            continue
        chrom_out = _strip_chr(chrom)
        out = {
            "motif_id": r.get("motif_id", ""),
            "motif_alt_id": r.get("motif_alt_id", ""),
            "sequence_name": chrom_out,
            "start": g_start,
            "stop": g_stop,
            "strand": r.get("strand", "+"),
            "score": r.get("score", ""),
            "p-value": r.get("p-value", r.get("p_value", "")),
            "q-value": r.get("q-value", r.get("q_value", "")),
            "matched_sequence": r.get("matched_sequence", ""),
        }
        rows.append(out)
    return pd.DataFrame(rows)


def main(argv: List[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help=f"Output TSV (default: {{PATHS.working_dir}}/SV/build_ccre_fimo/all_ccre_fimo.tsv)",
    )
    p.add_argument("--ccre-csv", type=Path, default=None, help="cCRE table (default: PATHS.ccre_csv)")
    p.add_argument("--meme", type=Path, default=None, help="MEME motif file (default: PATHS.sv_meme_file)")
    p.add_argument("--ref-fasta", type=Path, default=None, help="Reference FASTA (default: PATHS.sv_reference_fasta)")
    p.add_argument(
        "--thresh",
        type=float,
        default=None,
        help=f"FIMO --thresh p-value cap (default: THRESHOLDS.fimo_pvalue_threshold = {THRESHOLDS.fimo_pvalue_threshold})",
    )
    p.add_argument("--chunk-size", type=int, default=3000, help="cCRE rows per bedtools/FIMO batch")
    p.add_argument("--max-chunks", type=int, default=None, help="Stop after N chunks (debug)")
    p.add_argument("--work-dir", type=Path, default=None, help="Temp directory (default: system temp)")
    args = p.parse_args(argv)

    ccre_path = Path(args.ccre_csv or PATHS.ccre_csv).expanduser()
    meme_path = Path(args.meme or PATHS.sv_meme_file).expanduser()
    ref_path = Path(args.ref_fasta or PATHS.sv_reference_fasta).expanduser()
    thresh = float(args.thresh if args.thresh is not None else THRESHOLDS.fimo_pvalue_threshold)
    out_path = args.out
    if out_path is None:
        out_path = PATHS.working_dir / "SV" / "build_ccre_fimo" / "all_ccre_fimo.tsv"
    out_path = Path(out_path).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for label, path in ("cCRE CSV", ccre_path), ("MEME", meme_path), ("reference FASTA", ref_path):
        if not path.is_file():
            raise SystemExit(f"ERROR: {label} not found: {path}")

    bedtools = _require_bin("bedtools")
    fimo = _require_bin("fimo")

    print(f"[INFO] cCRE table: {ccre_path}", flush=True)
    print(f"[INFO] MEME:       {meme_path}", flush=True)
    print(f"[INFO] FASTA:      {ref_path}", flush=True)
    print(f"[INFO] FIMO --thresh: {thresh}", flush=True)
    print(f"[INFO] Output:     {out_path}", flush=True)

    ccre = _load_ccre(ccre_path)
    n = len(ccre)
    print(f"[INFO] Loaded {n} cCRE intervals.", flush=True)
    if n == 0:
        raise SystemExit("No cCRE rows after filtering.")

    chunk = max(1, int(args.chunk_size))
    bedtools_bin = bedtools
    fimo_bin = fimo

    all_parts: List[pd.DataFrame] = []
    tmp_root = Path(args.work_dir).expanduser() if args.work_dir else None

    n_chunks = (n + chunk - 1) // chunk
    done_chunks = 0
    with tempfile.TemporaryDirectory(dir=str(tmp_root) if tmp_root else None) as td:
        tdir = Path(td)
        for i in range(0, n, chunk):
            if args.max_chunks is not None and done_chunks >= int(args.max_chunks):
                break
            sub = ccre.iloc[i : i + chunk].copy()
            bounds: Dict[str, Tuple[str, int, int]] = {}
            bed_lines: List[str] = []
            for _, row in sub.iterrows():
                cid = str(row["cCRE_id"])
                chrom = str(row["chrom"])
                rs = int(row["start"])
                re = int(row["end"])
                bounds[cid] = (chrom, rs, re)
                # BED 0-based; name column for getfasta -name → FASTA header = cCRE_id
                bed_lines.append(f"{chrom}\t{rs}\t{re}\t{cid}\n")
            bed_path = tdir / f"chunk_{i:09d}.bed"
            fa_path = tdir / f"chunk_{i:09d}.fa"
            fimo_oc = tdir / f"chunk_{i:09d}_fimo"
            bed_path.write_text("".join(bed_lines), encoding="utf-8")

            done_chunks += 1
            print(
                f"[INFO] Chunk {done_chunks}/{n_chunks} rows {i}-{i + len(sub) - 1} "
                f"({len(sub)} cCREs) …",
                flush=True,
            )
            _run_getfasta(bedtools_bin, ref_path, bed_path, fa_path)
            if not fa_path.is_file() or fa_path.stat().st_size == 0:
                print(f"[WARN] Empty FASTA for chunk starting {i}; skipping FIMO.", flush=True)
                continue
            fimo_tsv = _run_fimo(fimo_bin, thresh, meme_path, fa_path, fimo_oc)
            raw = _read_fimo(fimo_tsv)
            mapped = _remap_chunk(raw, bounds)
            if not mapped.empty:
                all_parts.append(mapped)
            shutil.rmtree(fimo_oc, ignore_errors=True)

    if not all_parts:
        raise SystemExit("ERROR: No FIMO hits after processing all chunks (check MEME/FASTA/chrom names).")

    out_df = pd.concat(all_parts, ignore_index=True)
    # Stable order for diffs / merges
    out_df.sort_values(["sequence_name", "start", "stop", "motif_id"], inplace=True, kind="mergesort")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Wrote {len(out_df)} rows → {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
