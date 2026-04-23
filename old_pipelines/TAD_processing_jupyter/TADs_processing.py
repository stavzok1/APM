from __future__ import annotations

import gzip
import subprocess
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd


# -------------------------
# Core loaders / converters
# -------------------------

def read_arrowhead_domainlist(path_gz: str | Path) -> pd.DataFrame:
    """
    Read Arrowhead domain list (gzipped TSV) like:
    chr1 x1 x2 chr2 y1 y2 color f1 f2 f3 f4 f5

    Returns a DataFrame with raw columns.
    """
    path_gz = Path(path_gz)
    with gzip.open(path_gz, "rt") as f:
        df = pd.read_csv(f, sep="\t")

    # Normalize column names just in case
    df.columns = [c.strip() for c in df.columns]

    required = {"chr1", "x1", "x2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in Arrowhead file: {sorted(missing)}")

    return df


def _ensure_chr_prefix(chrom: pd.Series) -> pd.Series:
    # Input sometimes "1" instead of "chr1"
    chrom = chrom.astype(str)
    return chrom.where(chrom.str.startswith("chr"), "chr" + chrom)


def arrowhead_to_domains_bed(
    arrowhead_df: pd.DataFrame,
    add_score: bool = False,
    score_mode: str = "mean",  # "mean" or "max" over f1..f5
    name_prefix: str = "HMEC_Rao_Arrowhead_domain",
) -> pd.DataFrame:
    """
    Convert Arrowhead table to BED-like domain table:
      chrom, start, end, name[, score]

    If add_score=True, computes score from f1..f5.
    """
    df = arrowhead_df.copy()

    df["chrom"] = _ensure_chr_prefix(df["chr1"])
    df["start"] = df["x1"].astype(np.int64)
    df["end"] = df["x2"].astype(np.int64)

    # sanity: enforce start<end
    bad = (df["end"] <= df["start"])
    if bad.any():
        raise ValueError(f"Found {bad.sum()} rows with end <= start")

    df = df[["chrom", "start", "end"]].drop_duplicates().sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    df["name"] = [f"{name_prefix}_{i}" for i in range(len(df))]

    if add_score:
        fcols = [c for c in ["f1", "f2", "f3", "f4", "f5"] if c in arrowhead_df.columns]
        if len(fcols) < 1:
            raise ValueError("add_score=True but no f1..f5 columns were found.")
        vals = arrowhead_df.loc[df.index, fcols].astype(float)
        if score_mode == "mean":
            df["score"] = vals.mean(axis=1).values
        elif score_mode == "max":
            df["score"] = vals.max(axis=1).values
        else:
            raise ValueError("score_mode must be 'mean' or 'max'")

    return df


def domains_to_boundaries(
    domains_bed: pd.DataFrame,
    flank_bp: int = 20000,
    merge_within_bp: int = 20000,
    name_prefix: str = "HMEC_Rao_boundary",
) -> pd.DataFrame:
    """
    From domains (chrom,start,end), make boundaries:
      take start and end positions for each domain
      convert each pos -> window [pos-flank_bp, pos+flank_bp]
      then merge boundary windows within merge_within_bp per chromosome.

    Returns BED-like boundaries:
      chrom, start, end, name, n_support (how many raw boundaries collapsed here)
    """
    d = domains_bed.copy()

    # raw boundary positions (2 per domain)
    left = d[["chrom", "start"]].rename(columns={"start": "pos"})
    right = d[["chrom", "end"]].rename(columns={"end": "pos"})
    raw = pd.concat([left, right], ignore_index=True)

    # window them
    raw["start"] = (raw["pos"] - flank_bp).clip(lower=0).astype(np.int64)
    raw["end"] = (raw["pos"] + flank_bp).astype(np.int64)
    raw = raw[["chrom", "start", "end", "pos"]].sort_values(["chrom", "pos"]).reset_index(drop=True)

    # merge per chrom using a tolerance on positions (not windows) OR on starts/ends
    merged_rows = []
    for chrom, sub in raw.groupby("chrom", sort=False):
        sub = sub.sort_values("pos").reset_index(drop=True)
        # cluster boundary positions within merge_within_bp
        cluster_id = np.zeros(len(sub), dtype=np.int64)
        cid = 0
        last_pos = None
        for i, p in enumerate(sub["pos"].values):
            if last_pos is None or (p - last_pos) > merge_within_bp:
                cid += 1
            cluster_id[i] = cid
            last_pos = p

        sub["cluster"] = cluster_id
        for _, g in sub.groupby("cluster", sort=False):
            # merged window = min start, max end
            mstart = int(g["start"].min())
            mend = int(g["end"].max())
            n_support = int(len(g))
            # representative position (median)
            rep_pos = int(np.median(g["pos"].values))
            merged_rows.append((chrom, mstart, mend, rep_pos, n_support))

    out = pd.DataFrame(merged_rows, columns=["chrom", "start", "end", "pos", "n_support"])
    out = out.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    out["name"] = [f"{name_prefix}_{i}" for i in range(len(out))]
    return out[["chrom", "start", "end", "name", "pos", "n_support"]]





# -------------------------
# Liftover options
# -------------------------

def liftover_bed_ucsc(
    bed_df: pd.DataFrame,
    liftover_bin: str | Path,
    chain_file: str | Path,
    tmp_dir: str | Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Liftover using UCSC liftOver binary.

    bed_df must have: chrom,start,end,(name optional)
    Returns: (lifted_df, unmapped_df)
    """
    liftover_bin = str(Path(liftover_bin))
    chain_file = str(Path(chain_file))
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    in_bed = tmp_dir / "in.bed"
    out_bed = tmp_dir / "out.bed"
    un_bed = tmp_dir / "unmapped.bed"

    cols = bed_df.columns.tolist()
    name_col = "name" if "name" in cols else None

    # write BED4 (or BED3)
    to_write = bed_df.copy()
    if name_col is None:
        to_write = to_write[["chrom", "start", "end"]]
    else:
        to_write = to_write[["chrom", "start", "end", "name"]]
    to_write.to_csv(in_bed, sep="\t", header=False, index=False)

    cmd = [liftover_bin, str(in_bed), chain_file, str(out_bed), str(un_bed)]
    subprocess.run(cmd, check=True)

    # read results
    if out_bed.exists() and out_bed.stat().st_size > 0:
        if name_col is None:
            lifted = pd.read_csv(out_bed, sep="\t", header=None, names=["chrom", "start", "end"])
        else:
            lifted = pd.read_csv(out_bed, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    else:
        lifted = pd.DataFrame(columns=(["chrom", "start", "end", "name"] if name_col else ["chrom", "start", "end"]))

    if un_bed.exists() and un_bed.stat().st_size > 0:
        # unmapped lines have original coords + reason lines; keep raw text
        unmapped_text = un_bed.read_text()
        unmapped = pd.DataFrame({"raw": unmapped_text.splitlines()})
    else:
        unmapped = pd.DataFrame(columns=["raw"])

    return lifted, unmapped


def liftover_bed_pyliftover(
    bed_df: pd.DataFrame,
    chain_file: str | Path,
    drop_ambiguous: bool = True,
) -> pd.DataFrame:
    """
    Liftover using pyliftover (pip install pyliftover).
    Maps each interval by lifting its start and end independently.
    """
    from pyliftover import LiftOver

    lo = LiftOver(str(chain_file))

    rows = []
    for _, r in bed_df.iterrows():
        chrom, start, end = r["chrom"], int(r["start"]), int(r["end"])
        name = r["name"] if "name" in bed_df.columns else None

        s = lo.convert_coordinate(chrom, start)
        e = lo.convert_coordinate(chrom, end)

        if not s or not e:
            continue

        # s/e can have multiple mappings
        if drop_ambiguous and (len(s) != 1 or len(e) != 1):
            continue

        s0 = s[0]  # (newchrom, newpos, strand, score)
        e0 = e[0]

        if s0[0] != e0[0]:
            continue  # interval split across chroms after lift

        newchrom = s0[0]
        newstart = int(s0[1])
        newend = int(e0[1])

        if newend <= newstart:
            continue

        out = {"chrom": newchrom, "start": newstart, "end": newend}
        if name is not None:
            out["name"] = name
        rows.append(out)

    return pd.DataFrame(rows).sort_values(["chrom", "start", "end"]).reset_index(drop=True)


# -------------------------
# End-to-end runner
# -------------------------

def build_hmec_tads_from_arrowhead(
    arrowhead_gz: str | Path,
    out_dir: str | Path,
    flank_bp: int = 20000,
    merge_within_bp: int = 20000,
    liftover: bool = True,
    liftover_method: str = "ucsc",  # "ucsc" or "pyliftover"
    liftover_bin: Optional[str | Path] = None,  # required if method="ucsc"
    chain_file: Optional[str | Path] = None,     # required if liftover=True
) -> None:
    """
    Writes:
      - hmec_domains_hg19.bed
      - hmec_boundaries_hg19.bed
      - hmec_domains_hg38.bed (if liftover)
      - hmec_boundaries_hg38.bed (if liftover)

    Plus unmapped logs for UCSC method.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = read_arrowhead_domainlist(arrowhead_gz)
    domains_hg19 = arrowhead_to_domains_bed(raw, add_score=False)
    boundaries_hg19 = domains_to_boundaries(domains_hg19, flank_bp=flank_bp, merge_within_bp=merge_within_bp)

    # write hg19
    domains_hg19[["chrom", "start", "end", "name"]].to_csv(out_dir / "hmec_domains_hg19.bed", sep="\t", header=False, index=False)
    boundaries_hg19[["chrom", "start", "end", "name", "pos", "n_support"]].to_csv(
        out_dir / "hmec_boundaries_hg19.bed", sep="\t", header=False, index=False
    )

    if not liftover:
        return

    if chain_file is None:
        raise ValueError("liftover=True requires chain_file (hg19ToHg38.over.chain.gz)")

    # liftover domains and boundaries
    if liftover_method == "ucsc":
        if liftover_bin is None:
            raise ValueError("liftover_method='ucsc' requires liftover_bin path to UCSC liftOver executable")

        lifted_domains, unm_domains = liftover_bed_ucsc(
            domains_hg19[["chrom", "start", "end", "name"]],
            liftover_bin=liftover_bin,
            chain_file=chain_file,
            tmp_dir=out_dir / "tmp_liftover_domains",
        )
        lifted_boundaries, unm_bounds = liftover_bed_ucsc(
            boundaries_hg19[["chrom", "start", "end", "name"]],
            liftover_bin=liftover_bin,
            chain_file=chain_file,
            tmp_dir=out_dir / "tmp_liftover_boundaries",
        )

        

        lifted_domains.to_csv(out_dir / "hmec_domains_hg38.bed", sep="\t", header=False, index=False)
        lifted_boundaries.to_csv(out_dir / "hmec_boundaries_hg38.bed", sep="\t", header=False, index=False)

        unm_domains.to_csv(out_dir / "hmec_domains_unmapped.txt", sep="\t", index=False)
        unm_bounds.to_csv(out_dir / "hmec_boundaries_unmapped.txt", sep="\t", index=False)

    elif liftover_method == "pyliftover":
        lifted_domains = liftover_bed_pyliftover(domains_hg19[["chrom", "start", "end", "name"]], chain_file=chain_file)
        lifted_boundaries = liftover_bed_pyliftover(boundaries_hg19[["chrom", "start", "end", "name"]], chain_file=chain_file)

        lifted_domains.to_csv(out_dir / "hmec_domains_hg38.bed", sep="\t", header=False, index=False)
        lifted_boundaries.to_csv(out_dir / "hmec_boundaries_hg38.bed", sep="\t", header=False, index=False)

    else:
        raise ValueError("liftover_method must be 'ucsc' or 'pyliftover'")


# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    build_hmec_tads_from_arrowhead(
        arrowhead_gz=r"C:\Users\stavz\Desktop\masters\APM\data\TADs\HMEC\GSE63525_HMEC_Arrowhead_domainlist.txt.gz",
        out_dir=r"C:\Users\stavz\Desktop\masters\APM\data\TADs\HMEC",
        flank_bp=20000,         # boundary window half-width (20kb => 40kb region)
        merge_within_bp=20000,  # merge boundaries within 20kb
        liftover=True,
        liftover_method="ucsc",  # or "pyliftover"
        liftover_bin=r"liftOver",  # required for ucsc method
        chain_file=r"\\wsl.localhost\Ubuntu\home\stavz\refs\ucsc_chains\hg19ToHg38.over.chain.gz",
    )
