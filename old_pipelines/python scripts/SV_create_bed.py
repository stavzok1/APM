import ast
import pandas as pd
import numpy as np
import os
import sys

def build_sv_flanks_and_overlaps_bed(
    sv_sample_path: str,
    output_dir: str,
    flank: int = 150
) -> pd.DataFrame:

    s = pd.read_csv(sv_sample_path)
    file_name = os.path.basename(sv_sample_path).split(".")[0]

    # output file
    os.makedirs(output_dir, exist_ok=True)
    out_bed = os.path.join(output_dir, f"{file_name}_flanks_and_overlaps.bed")

    mask = s["SVTYPE"].isin(["DEL", "DUP", "INS", "BND"]) & (s["elem_hits"] != "[]")
    df = s[mask].copy()

    if "chr" in str(df["chrom"].iloc[0]):
        df["chrom"] = df["chrom"].str.replace("chr", "", regex=False)

    records = []

    for idx, row in df.iterrows():
        chrom = str(row["chrom"])
        pos   = int(row["pos"])
        svt   = str(row["SVTYPE"])
        sv_id = str(row.get("id", idx))

        end_val = row.get("END", np.nan)
        has_end = not pd.isna(end_val)
        end = int(end_val) if has_end else None

        # ---- Flanks ----
        if svt == "DEL":
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank_left"))

            if end is not None:
                s2 = max(1, end - flank)
                e2 = end + flank
                records.append((chrom, s2 - 1, e2, f"{sv_id}|{svt}_flank_right"))

        elif svt == "INS":
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank_left"))

            center2 = pos + 1
            s2 = max(1, center2 - flank)
            e2 = center2 + flank
            records.append((chrom, s2 - 1, e2, f"{sv_id}|{svt}_flank_right"))

        elif svt == "BND":
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank"))

        # ---- Element overlaps ----
        try:
            hits = ast.literal_eval(row["elem_hits"])
        except Exception:
            continue

        if not isinstance(hits, list):
            continue

        for h in hits:
            if not isinstance(h, dict):
                continue
            if not h.get("overlaps_flag", 0):
                continue

            e_chrom = str(h["chrom"]).replace("chr", "")
            o_start = int(h["overlap_start"])
            o_end   = int(h["overlap_end"])

            elem_id   = str(h["elem_id"])
            elem_type = str(h.get("elem_type", ""))
            hit_side  = str(h.get("hit_side", ""))

            name = f"{sv_id}|{svt}|elem:{elem_id}|type:{elem_type}|hit_side:{hit_side}"

            records.append((e_chrom, o_start - 1, o_end, name))

    bed = pd.DataFrame.from_records(records, columns=["chrom", "start", "end", "name"])
    bed.to_csv(out_bed, sep="\t", header=False, index=False)

    print(f"Wrote: {out_bed}")
    return bed


# ==================================================
# MAIN — process a whole directory
# ==================================================
if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    flank = int(sys.argv[3])

    os.makedirs(output_dir, exist_ok=True)

    # Process all .csv files in directory
    for fname in os.listdir(input_dir):
        if fname.endswith(".csv") or fname.endswith(".tsv"):
            fpath = os.path.join(input_dir, fname)
            print(f"Processing: {fpath}")
            build_sv_flanks_and_overlaps_bed(fpath, output_dir, flank)
