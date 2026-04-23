"""
VCF loading and parsing for Manta SV calls.

Functions for:
- Loading VCF files with vcfpy
- Extracting sample information
- Parsing BND ALT alleles for remote coordinates
"""

import gzip
import re
from pathlib import Path
from typing import Tuple, Optional, List

import numpy as np
import pandas as pd
import vcfpy


# =============================================================================
# VCF LOADING
# =============================================================================


def _open_manta_sv_reader(vcf_path: str) -> vcfpy.Reader:
    """
    Open a Manta VCF for vcfpy. Uses gzip only when content has gzip magic bytes.
    GDC-style names often end in ``.vcf.gz`` even when the body is uncompressed VCF;
    ``vcfpy.Reader.from_path`` would treat the suffix as gzip and fail.
    """
    path = Path(vcf_path)
    with path.open("rb") as bf:
        magic = bf.read(2)
    is_gzip = magic == b"\x1f\x8b"
    lower = path.name.lower()
    looks_compressed_suffix = lower.endswith(".gz") or lower.endswith(".bgz")

    if is_gzip:
        stream = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    elif looks_compressed_suffix and not is_gzip:
        stream = path.open("rt", encoding="utf-8", errors="replace")
    else:
        stream = path.open("rt", encoding="utf-8", errors="replace")

    return vcfpy.Reader.from_stream(stream, path=str(path))


def load_manta_sv_vcf(vcf_path: str) -> Tuple[pd.DataFrame, str, str]:
    """
    Load a Manta SV VCF file into a DataFrame.
    
    Extracts all INFO fields and sample-level FORMAT data (PR, SR).
    Automatically detects normal and tumor samples.
    
    Args:
        vcf_path: Path to VCF file
    
    Returns:
        Tuple of (DataFrame, normal_sample_name, tumor_sample_name)
    """
    reader = _open_manta_sv_reader(vcf_path)
    try:
        samples = reader.header.samples.names
        print(f"Samples in VCF: {samples}")

        # Identify normal/tumor samples
        normal_sample = None
        tumor_sample = None
        for s in samples:
            su = s.upper()
            if "NORMAL" in su:
                normal_sample = s
            if "TUMOR" in su:
                tumor_sample = s

        # Fallback to positional
        if normal_sample is None and len(samples) >= 1:
            normal_sample = samples[0]
        if tumor_sample is None and len(samples) >= 2:
            tumor_sample = samples[-1]

        print(f"Assuming normal sample: {normal_sample}")
        print(f"Assuming tumor sample: {tumor_sample}")

        rows = []

        for rec in reader:
            row = {
                "chrom": rec.CHROM,
                "pos": rec.POS,
                "id": rec.ID[0] if rec.ID else None,
                "ref": rec.REF,
                "alt": ",".join(
                    getattr(a, "value", str(a)) for a in rec.ALT
                ),
                "qual": rec.QUAL,
                "filter": ";".join(rec.FILTER),
            }

            # INFO fields: flatten single-element lists
            for key, value in rec.INFO.items():
                if isinstance(value, list) and len(value) == 1:
                    row[key] = value[0]
                else:
                    row[key] = value

            # FORMAT / calls: PR and SR for each sample
            for sample in samples:
                call = rec.call_for_sample[sample]
                pr = call.data.get("PR")
                sr = call.data.get("SR")

                row[f"{sample}_PR_ref"] = pr[0] if pr else np.nan
                row[f"{sample}_PR_alt"] = pr[1] if pr else np.nan
                row[f"{sample}_SR_ref"] = sr[0] if sr else np.nan
                row[f"{sample}_SR_alt"] = sr[1] if sr else np.nan

            rows.append(row)

        df = pd.DataFrame(rows)

        # Ensure numeric fields
        numeric_cols = [
            "SVLEN", "SOMATICSCORE",
            f"{normal_sample}_PR_ref", f"{normal_sample}_PR_alt",
            f"{normal_sample}_SR_ref", f"{normal_sample}_SR_alt",
            f"{tumor_sample}_PR_ref", f"{tumor_sample}_PR_alt",
            f"{tumor_sample}_SR_ref", f"{tumor_sample}_SR_alt",
        ]
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # Convenience columns: total alt counts
        df["normal_alt"] = (
            df.get(f"{normal_sample}_PR_alt", 0).fillna(0) +
            df.get(f"{normal_sample}_SR_alt", 0).fillna(0)
        )
        df["tumor_alt"] = (
            df.get(f"{tumor_sample}_PR_alt", 0).fillna(0) +
            df.get(f"{tumor_sample}_SR_alt", 0).fillna(0)
        )

        df["normal_sr_alt"] = df.get(f"{normal_sample}_SR_alt", 0).fillna(0)
        df["tumor_sr_alt"] = df.get(f"{tumor_sample}_SR_alt", 0).fillna(0)
        df["tumor_pr_alt"] = df.get(f"{tumor_sample}_PR_alt", 0).fillna(0)

        return df, normal_sample, tumor_sample
    finally:
        try:
            reader.stream.close()
        except Exception:
            pass


# =============================================================================
# BND PARSING
# =============================================================================

def parse_breakend_string(s: str) -> Optional[Tuple[str, int, str, str]]:
    """
    Parse a BreakEnd string representation.
    
    Parses: BreakEnd('chrX', 101138003, '-', '-', 'T', True)
    
    Returns:
        Tuple of (chrom2, pos2, orientation_self, orientation_remote)
        or None if parsing fails.
    """
    if not isinstance(s, str):
        return None

    m = re.search(
        r"BreakEnd\('([^']+)',\s*([0-9]+),\s*'([+-])',\s*'([+-])'",
        s
    )
    if not m:
        return None

    chrom2 = m.group(1)
    pos2 = int(m.group(2))
    orientation_self = m.group(3)
    orientation_remote = m.group(4)

    return chrom2, pos2, orientation_self, orientation_remote


def get_bnd_remote_coords(alt_value) -> Tuple[Optional[str], Optional[int]]:
    """
    Extract remote chromosome + position from a Manta BND ALT.

    Handles both BreakEnd object repr and VCF bracket notation.
    
    Returns:
        Tuple of (chrom2, pos2) or (None, None).
    """
    if alt_value is None or (isinstance(alt_value, float) and np.isnan(alt_value)):
        return None, None

    # If list/tuple, take first
    if isinstance(alt_value, (list, tuple)):
        if len(alt_value) == 0:
            return None, None
        alt_value = alt_value[0]

    s = str(alt_value)

    # Pattern 1: BreakEnd('chr5', 16982822, ...
    m = re.search(r"BreakEnd\('([^']+)',\s*([0-9]+)", s)
    if m:
        chrom2 = m.group(1)
        pos2 = int(m.group(2))
        return chrom2, pos2

    # Pattern 2: VCF bracket notation like ]chr5:12345] or [chr5:12345[
    m = re.search(r"[\[\]](chr[0-9XYM]+):([0-9]+)[\[\]]", s)
    if m:
        chrom2 = m.group(1)
        pos2 = int(m.group(2))
        return chrom2, pos2

    # Fallback: chrXXX:12345 or chrXXX, 12345
    m = re.search(r"(chr[0-9XYM]+)[^\d]+([0-9]+)", s)
    if m:
        chrom2 = m.group(1)
        pos2 = int(m.group(2))
        return chrom2, pos2

    return None, None


def breakend_to_vcf_alt(be_str: str) -> Optional[str]:
    """
    Convert BreakEnd string to VCF ALT format.
    
    Args:
        be_str: String like BreakEnd('chrX', 101138003, '-', '-', 'T', True)
    
    Returns:
        VCF-style ALT string like N]chrX:101138003] or None if parsing fails.
    """
    parsed = parse_breakend_string(be_str)
    if parsed is None:
        return None
    
    chrom2, pos2, orientation_self, orientation_remote = parsed

    # VCF BND construction rules
    if orientation_self == "+" and orientation_remote == "+":
        return f"N]{chrom2}:{pos2}]"
    elif orientation_self == "+" and orientation_remote == "-":
        return f"N[{chrom2}:{pos2}["
    elif orientation_self == "-" and orientation_remote == "+":
        return f"]{chrom2}:{pos2}]N"
    else:  # "-" and "-"
        return f"[{chrom2}:{pos2}[N"


def add_bnd_remote_coords(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add bnd_remote_chrom and bnd_remote_pos columns for BND variants.
    
    Args:
        df: DataFrame with SVTYPE and alt columns
    
    Returns:
        DataFrame with new columns added
    """
    df = df.copy()

    bnd_mask = df["SVTYPE"] == "BND"

    # Nullable string / int columns (float64 NaN placeholders caused LossySetitemError
    # when assigning chr* strings from BND remote coords).
    n = len(df)
    df["bnd_remote_chrom"] = pd.Series([pd.NA] * n, index=df.index, dtype="string")
    df["bnd_remote_pos"] = pd.Series([pd.NA] * n, index=df.index, dtype="Int64")

    if not bnd_mask.any():
        return df

    coords = df.loc[bnd_mask, "alt"].apply(get_bnd_remote_coords)

    coords_df = pd.DataFrame(
        coords.tolist(),
        index=coords.index,
        columns=["bnd_remote_chrom", "bnd_remote_pos"],
    )

    df.loc[coords_df.index, "bnd_remote_chrom"] = coords_df["bnd_remote_chrom"].astype("string")
    df.loc[coords_df.index, "bnd_remote_pos"] = pd.to_numeric(
        coords_df["bnd_remote_pos"], errors="coerce"
    ).astype("Int64")

    return df
