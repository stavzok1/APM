import os
import gzip
import subprocess
from pathlib import Path

import pandas as pd

########################
# CONFIG
########################

LIFTOVER_BIN = "liftOver"  # name or full path to UCSC liftOver
CHAIN_FILE   = "~/refs/ucsc_chains/hg19ToHg38.over.chain.gz"

# Filters
MIN_PETS   = 2
MIN_DIST   = 5_000        # 5 kb
MAX_DIST   = 2_000_000    # 2 Mb


########################
# 1. Load + basic filtering in hg19
########################

def load_loop_counts_bedpe(path,
                           min_pets=MIN_PETS,
                           min_dist=MIN_DIST,
                           max_dist=MAX_DIST):
    """
    Load a *.intra.loop_counts.bedpe.gz file in hg19 and filter:

      - keep only intra-chrom loops (chr1 == chr2)
      - PETs >= min_pets
      - min_dist <= distance <= max_dist

    Returns a DataFrame with:
      ["chr1","start1","end1","chr2","start2","end2","PETs"]
    """
    path = Path(path)
    if path.suffix == ".gz":
        opener = gzip.open
    else:
        opener = open

    cols = [
        "chr1","start1","end1","chr2","start2","end2",
        "FDR","PETs"
    ]

    df = pd.read_csv(path, sep="\t", header=None,
                     names=cols, comment="#", compression="infer")

    # Keep intra only
    df = df[df["chr1"] == df["chr2"]].copy()

    # Ensure numeric
    df["start1"] = df["start1"].astype(int)
    df["end1"]   = df["end1"].astype(int)
    df["start2"] = df["start2"].astype(int)
    df["end2"]   = df["end2"].astype(int)
    df["PETs"]   = df["PETs"].astype(int)

    # Distance between loop centers
    center1 = (df["start1"] + df["end1"]) // 2
    center2 = (df["start2"] + df["end2"]) // 2
    dist = (center2 - center1).abs()

    df = df[(df["PETs"] >= min_pets) &
            (dist >= min_dist) &
            (dist <= max_dist)].copy()

    df = df[["chr1","start1","end1","chr2","start2","end2","PETs"]]
    df.reset_index(drop=True, inplace=True)
    df["loop_id"] = df.index  # unique id within this replicate

    return df


########################
# 2. Write anchors to BED, run liftOver, read back
########################

def write_anchor_bed(df, side, out_bed):
    """
    side = 1 or 2 for anchor1 / anchor2
    Creates a BED:
      chr, start, end, loop_id
    """
    assert side in (1, 2)
    cols = [f"chr{side}", f"start{side}", f"end{side}", "loop_id"]
    bed = df[cols].copy()
    bed.columns = ["chr", "start", "end", "name"]
    bed.to_csv(out_bed, sep="\t", header=False, index=False)


def run_liftover(in_bed, out_bed, unmapped_bed, chain_file=CHAIN_FILE):
    """
    Call UCSC liftOver via subprocess.
    """
    cmd = [
        LIFTOVER_BIN,
        str(in_bed),
        str(chain_file),
        str(out_bed),
        str(unmapped_bed),
    ]
    subprocess.run(cmd, check=True)


def read_lifted_bed(in_bed):
    """
    Read lifted BED:
      chr, start, end, name
    Returns DataFrame with columns:
      ["chr", "start", "end", "loop_id"]
    """
    df = pd.read_csv(in_bed, sep="\t", header=None,
                     names=["chr","start","end","loop_id"])
    return df


########################
# 3. Liftover one replicate hg19 -> hg38
########################

def liftover_bedpe_df(df_hg19, tmp_dir, cell_line, rep_label):
    """
    Given a filtered hg19 BEDPE DataFrame with columns:
      ["chr1","start1","end1","chr2","start2","end2","PETs","loop_id"]
    Perform liftOver on each anchor separately and return a hg38 loop DF.

    tmp_dir is used for intermediate BED files.
    """
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # ---- anchor1 ----
    a1_in  = tmp_dir / f"{cell_line}_{rep_label}_a1_hg19.bed"
    a1_out = tmp_dir / f"{cell_line}_{rep_label}_a1_hg38.bed"
    a1_unm = tmp_dir / f"{cell_line}_{rep_label}_a1_unmapped.bed"

    write_anchor_bed(df_hg19, side=1, out_bed=a1_in)
    run_liftover(a1_in, a1_out, a1_unm)
    a1_hg38 = read_lifted_bed(a1_out)
    a1_hg38.rename(columns={
        "chr": "chr1_hg38",
        "start": "start1_hg38",
        "end": "end1_hg38",
        "loop_id": "loop_id"
    }, inplace=True)

    # ---- anchor2 ----
    a2_in  = tmp_dir / f"{cell_line}_{rep_label}_a2_hg19.bed"
    a2_out = tmp_dir / f"{cell_line}_{rep_label}_a2_hg38.bed"
    a2_unm = tmp_dir / f"{cell_line}_{rep_label}_a2_unmapped.bed"

    write_anchor_bed(df_hg19, side=2, out_bed=a2_in)
    run_liftover(a2_in, a2_out, a2_unm)
    a2_hg38 = read_lifted_bed(a2_out)
    a2_hg38.rename(columns={
        "chr": "chr2_hg38",
        "start": "start2_hg38",
        "end": "end2_hg38",
        "loop_id": "loop_id"
    }, inplace=True)

    # Merge anchors with original PETs
    merged = df_hg19[["loop_id","PETs"]].merge(
        a1_hg38, on="loop_id", how="inner"
    ).merge(
        a2_hg38, on="loop_id", how="inner"
    )

    # keep only intra-chrom after liftOver
    merged = merged[merged["chr1_hg38"] == merged["chr2_hg38"]].copy()

    # Final schema
    merged = merged[[
        "chr1_hg38","start1_hg38","end1_hg38",
        "chr2_hg38","start2_hg38","end2_hg38",
        "PETs"
    ]].copy()

    merged.rename(columns={
        "chr1_hg38": "chr1",
        "start1_hg38": "start1",
        "end1_hg38": "end1",
        "chr2_hg38": "chr2",
        "start2_hg38": "start2",
        "end2_hg38": "end2",
    }, inplace=True)

    return merged


########################
# 4. Merge two replicates into one loop set
########################

def merge_replicates(df1, df2):
    """
    df1, df2 are hg38 DataFrames from liftover with columns:
      ["chr1","start1","end1","chr2","start2","end2","PETs"]

    Returns merged DF with:
      ["chr1","start1","end1","chr2","start2","end2","weight","n_reps"]
    where:
      - weight = sum PETs across replicates
      - n_reps = in how many reps the loop appeared (1 or 2)
    """
    df1 = df1.copy()
    df2 = df2.copy()
    df1["rep"] = "r1"
    df2["rep"] = "r2"

    df_all = pd.concat([df1, df2], ignore_index=True)

    group_cols = ["chr1","start1","end1","chr2","start2","end2"]

    agg = df_all.groupby(group_cols).agg(
        weight=("PETs","sum"),
        n_reps=("rep","nunique")
    ).reset_index()

    return agg


########################
# 5. Convenience driver
########################

def process_cell_line(
    cell_line,
    rep1_path,
    rep2_path,
    tmp_dir,
    out_path,
    min_pets=MIN_PETS,
    min_dist=MIN_DIST,
    max_dist=MAX_DIST
):
    """
    Full pipeline for one cell line (e.g., MCF7 or ZR751):

      1. Load + filter rep1 & rep2 hg19
      2. Liftover each to hg38
      3. Merge replicates
      4. Save hg38 BEDPE-like file with columns:
         chr1, start1, end1, chr2, start2, end2, weight, n_reps
    """
    print(f"[{cell_line}] Loading rep1:", rep1_path)
    df1_hg19 = load_loop_counts_bedpe(rep1_path,
                                      min_pets=min_pets,
                                      min_dist=min_dist,
                                      max_dist=max_dist)

    print(f"[{cell_line}] Loading rep2:", rep2_path)
    df2_hg19 = load_loop_counts_bedpe(rep2_path,
                                      min_pets=min_pets,
                                      min_dist=min_dist,
                                      max_dist=max_dist)

    print(f"[{cell_line}] Liftover rep1 -> hg38")
    df1_hg38 = liftover_bedpe_df(df1_hg19, tmp_dir, cell_line, "rep1")

    print(f"[{cell_line}] Liftover rep2 -> hg38")
    df2_hg38 = liftover_bedpe_df(df2_hg19, tmp_dir, cell_line, "rep2")

    print(f"[{cell_line}] Merging replicates")
    merged = merge_replicates(df1_hg38, df2_hg38)

    # Sort for sanity
    merged.sort_values(["chr1","start1","chr2","start2"], inplace=True)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[{cell_line}] Saving merged loops to:", out_path)
    merged.to_csv(out_path, sep="\t", header=True, index=False)

    print(f"[{cell_line}] Done. N loops:", len(merged))


########################
# Example usage
########################

if __name__ == "__main__":
    # Example paths – adjust to your actual file locations

    # -------- MCF7 --------
    process_cell_line(
        cell_line="MCF7",
        rep1_path="~/masters/gdc/HiCHIP/MCF7/MCF7_H3K27ac_HiChIP_rep1.intra.loop_counts.bedpe.gz",
        rep2_path="~/masters/gdc/HiCHIP/MCF7/MCF7_H3K27ac_HiChIP_rep2.intra.loop_counts.bedpe.gz",
        tmp_dir="/tmp/hichip_MCF7",
        out_path="~/masters/gdc/HiCHIP/MCF7/MCF7_H3K27ac_loops_hg38.bedpe"
    )

    # -------- ZR751 --------
    process_cell_line(
        cell_line="ZR751",
        rep1_path="~/masters/gdc/HiCHIP/ZR751/ZR751_H3K27ac_HiChIP_rep1.intra.loop_counts.bedpe.gz",
        rep2_path="~/masters/gdc/HiCHIP/ZR751/ZR751_H3K27ac_HiChIP_rep2.intra.loop_counts.bedpe.gz",
        tmp_dir="/tmp/hichip_ZR751",
        out_path="~/masters/gdc/HiCHIP/ZR751/ZR751_H3K27ac_loops_hg38.bedpe"
    )
