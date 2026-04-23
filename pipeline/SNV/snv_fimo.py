"""
MEME FIMO on fixed-width hg38 windows around each SNV (reference allele).

This is **orthogonal** to VEP ``motif_hits`` (consequence on alternate allele):
here we scan the **reference** sequence with the same motif database used for
the SV pipeline (``PATHS.sv_meme_file`` by default).

Requires ``bedtools``. ``fimo`` is resolved via ``resolve_fimo_argv()`` (``APM_FIMO_BIN``,
PATH, or ``VEP_ENV`` / ``vep_env`` conda roots).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from ..SV.motif_scanning import _load_fimo_tsv_as_dataframe
from ..utils import harmonize_chrom_column


def resolve_fimo_argv() -> List[str]:
    """
    Return argv prefix to run FIMO (``[fimo_bin]``).

    Resolution order:
    1. ``APM_FIMO_BIN`` env (full path to ``fimo``)
    2. ``fimo`` on ``PATH``
    3. ``$VEP_ENV/bin/fimo`` or ``$vep_env/bin/fimo`` when set to a conda env root
    4. Common conda locations for env ``vep_env`` (including Miniforge)
    """
    apm = os.environ.get("APM_FIMO_BIN", "").strip()
    if apm and Path(apm).is_file():
        return [apm]
    which = shutil.which("fimo")
    if which:
        return [which]
    for key in ("VEP_ENV", "vep_env", "VENV_ENV"):
        root = os.environ.get(key, "").strip()
        if not root:
            continue
        p = Path(root).expanduser() / "bin" / "fimo"
        if p.is_file():
            return [str(p)]
    for home_sub in (
        Path.home() / "miniforge3/envs/vep_env",
        Path.home() / "miniconda3/envs/vep_env",
        Path.home() / "anaconda3/envs/vep_env",
        Path.home() / "mambaforge/envs/vep_env",
        Path.home() / "micromamba/envs/vep_env",
    ):
        p = home_sub / "bin" / "fimo"
        if p.is_file():
            return [str(p)]
    return ["fimo"]


def _tf_from_motif_id(motif_id: str) -> str:
    s = str(motif_id)
    return s.split(".", 1)[0] if "." in s else s


def variant_reference_window_0based(pos: int, flank_bp: int) -> tuple[int, int]:
    """1-based VCF ``POS`` -> 0-based half-open BED ``[start, end)`` centered on that base."""
    p = int(pos)
    start0 = max(0, p - flank_bp - 1)
    end_excl = p + flank_bp
    return start0, end_excl


def build_snv_fimo_bed(df: pd.DataFrame, bed_out: Path, flank_bp: int) -> tuple[List[int], List[int]]:
    """
    Write a 4-column BED: chrom, start0, end, name with ``name`` = 0..n-1 row order.

    Returns parallel lists ``window_starts0``, ``variant_positions`` (1-based POS per row).
    """
    work = df.copy()
    work, _ = harmonize_chrom_column(work, "chrom")

    lines: List[str] = []
    window_starts0: List[int] = []
    variant_positions: List[int] = []

    for j in range(len(work)):
        row = work.iloc[j]
        chrom = str(row["chrom"]).strip()
        pos = int(row["pos"])
        s0, e_excl = variant_reference_window_0based(pos, flank_bp)
        window_starts0.append(s0)
        variant_positions.append(pos)
        lines.append(f"{chrom}\t{s0}\t{e_excl}\t{j}")

    bed_out.parent.mkdir(parents=True, exist_ok=True)
    bed_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return window_starts0, variant_positions


def _run_bedtools_getfasta(bed: Path, ref_fasta: Path, fasta_out: Path) -> None:
    cmd = [
        "bedtools",
        "getfasta",
        "-fi",
        str(ref_fasta),
        "-bed",
        str(bed),
        "-name",
        "-fo",
        str(fasta_out),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def _run_fimo(meme_file: Path, fasta: Path, out_dir: Path, thresh: float) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    fimo_bin = resolve_fimo_argv()
    cmd = [
        *fimo_bin,
        "--oc",
        str(out_dir),
        "--thresh",
        str(thresh),
        str(meme_file),
        str(fasta),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    fimo_tsv = out_dir / "fimo.tsv"
    if not fimo_tsv.is_file():
        raise FileNotFoundError(f"FIMO did not write {fimo_tsv}")
    return fimo_tsv


def aggregate_fimo_hits_per_variant(
    fimo_df: pd.DataFrame,
    window_starts0: List[int],
    variant_positions: List[int],
    max_hits_per_variant: int,
) -> List[List[Dict[str, Any]]]:
    """
    Map FIMO rows (``sequence_name`` = integer row index) to per-variant hit lists.

    FIMO coordinates are **1-based inclusive** within each extracted window; ``window_starts0[i]``
    is the hg38 0-based start of that window's first base.
    """
    n = len(window_starts0)
    raw: Dict[int, List[Dict[str, Any]]] = defaultdict(list)

    if (
        fimo_df.empty
        or "start" not in fimo_df.columns
        or "stop" not in fimo_df.columns
        or "sequence_name" not in fimo_df.columns
    ):
        return [[] for _ in range(n)]

    for _, row in fimo_df.iterrows():
        try:
            idx = int(str(row["sequence_name"]).strip().split()[0])
        except (ValueError, TypeError, AttributeError):
            continue
        if idx < 0 or idx >= n:
            continue

        motif_id = row.get("motif_id")
        if motif_id is None or (isinstance(motif_id, float) and pd.isna(motif_id)):
            continue

        try:
            rel_start = int(row["start"])
            rel_stop = int(row["stop"])
        except (ValueError, TypeError):
            continue

        g0 = window_starts0[idx]
        g_start0 = g0 + rel_start - 1
        g_end_excl = g0 + rel_stop

        pval = row.get("p-value")
        qval = row.get("q-value")
        score = row.get("score")
        strand = row.get("strand")
        matched = row.get("matched_sequence")

        def _maybe_float(x: object) -> Optional[float]:
            if x is None:
                return None
            try:
                if isinstance(x, float) and pd.isna(x):
                    return None
                return float(x)
            except (TypeError, ValueError):
                return None

        hit: Dict[str, Any] = {
            "TF": _tf_from_motif_id(str(motif_id)),
            "motif_id": str(motif_id),
            "rel_start": rel_start,
            "rel_stop": rel_stop,
            "start": int(g_start0),
            "end": int(g_end_excl),
            "strand": None
            if strand is None or (isinstance(strand, float) and pd.isna(strand))
            else str(strand),
            "p_value": _maybe_float(pval),
            "q_value": _maybe_float(qval),
            "score": _maybe_float(score),
            "variant_pos": int(variant_positions[idx]),
            "matched_sequence": None
            if matched is None or (isinstance(matched, float) and pd.isna(matched))
            else str(matched),
        }
        raw[idx].append(hit)

    out: List[List[Dict[str, Any]]] = []
    for i in range(n):
        hs = raw.get(i, [])
        hs.sort(
            key=lambda h: (h["p_value"] is None, h["p_value"] if h["p_value"] is not None else 1.0)
        )
        out.append(hs[:max_hits_per_variant])
    return out


def annotate_snvs_with_fimo(
    df: pd.DataFrame,
    *,
    ref_fasta: Path,
    meme_file: Path,
    work_dir: Path,
    flank_bp: int = 30,
    fimo_threshold: float = 1e-4,
    max_hits_per_variant: int = 50,
    fimo_tsv_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Add a ``fimo_hits`` column: list of motif dicts per row (possibly empty).

    Unless ``fimo_tsv_path`` is set, runs ``bedtools getfasta`` then ``fimo`` under ``work_dir``.
    ``fimo_tsv_path`` is intended for tests or for re-attaching after a manual FIMO run.

    Args:
        df: Filtered SNV table; must include ``chrom``, ``pos``, ``ref``, ``alt``.
        ref_fasta: hg38 reference (same convention as SV FIMO).
        meme_file: MEME motif file (e.g. ``PATHS.sv_meme_file``).
        work_dir: Writable directory for BED/FASTA/FIMO outputs.
        flank_bp: Half-width around ``POS`` (inclusive span ``2 * flank_bp + 1`` bp).
        fimo_threshold: Passed to ``fimo --thresh``.
        max_hits_per_variant: Cap hits per variant (best by ``p_value``).
        fimo_tsv_path: If provided, skip scanning and parse this TSV (must use BED names ``0``..``n-1``
            matching the row order of ``df`` after ``reset_index(drop=True)``).

    Returns:
        Copy of ``df`` with ``fimo_hits`` added (row order preserved; index reset to RangeIndex).
    """
    out = df.copy().reset_index(drop=True)
    if out.empty:
        out["fimo_hits"] = []
        return out

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    ref_fasta = Path(ref_fasta).expanduser()
    meme_file = Path(meme_file)
    if not ref_fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")
    if not meme_file.is_file():
        raise FileNotFoundError(f"MEME motif file not found: {meme_file}")

    bed_path = work_dir / "snv_fimo_windows.bed"
    fa_path = work_dir / "snv_fimo_windows.fa"
    fimo_run_dir = work_dir / "fimo_run"

    window_starts0, variant_positions = build_snv_fimo_bed(out, bed_path, flank_bp)

    if fimo_tsv_path is not None:
        tsv_path = Path(fimo_tsv_path)
        if not tsv_path.is_file():
            raise FileNotFoundError(f"fimo_tsv_path not found: {tsv_path}")
    else:
        shutil.rmtree(fimo_run_dir, ignore_errors=True)
        _run_bedtools_getfasta(bed_path, ref_fasta, fa_path)
        tsv_path = _run_fimo(meme_file, fa_path, fimo_run_dir, fimo_threshold)

    fimo_df = _load_fimo_tsv_as_dataframe(tsv_path)
    hits_cols = aggregate_fimo_hits_per_variant(
        fimo_df, window_starts0, variant_positions, max_hits_per_variant
    )
    out["fimo_hits"] = hits_cols
    return out
