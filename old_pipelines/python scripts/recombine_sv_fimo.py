#!/usr/bin/env python3
import os
import sys
import ast
from collections import defaultdict

import pandas as pd
import numpy as np


def parse_elem_hits_cell(x):
    """Parse elem_hits cell into a list of dicts."""
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    if x == "" or x == "[]":
        return []
    try:
        return ast.literal_eval(x)
    except Exception:
        return []


def parse_flank_side(name_tail: str) -> str:
    """
    Example name_tail: 'DEL_flank_left' or 'INS_flank_right'
    We want 'left' or 'right'.
    """
    parts = name_tail.split("_")
    if len(parts) >= 3 and parts[-2] == "flank":
        return parts[-1]
    # fallback: last part
    return parts[-1]


def extract_elem_id_from_name(name: str) -> str | None:
    """
    Extract elem id from something like:
    'MantaDUP:...|DUP|elem:EH38D2925526|type:dELS|hit_side:span'
    """
    for part in name.split("|"):
        if part.startswith("elem:"):
            return part[len("elem:"):]
    return None


def recombine_one(merged_bed_path: str, sv_csv_path: str, out_csv_path: str):
    print(f"[INFO] Recombining:\n  BED: {merged_bed_path}\n  CSV: {sv_csv_path}\n  OUT: {out_csv_path}")

    if not os.path.isfile(merged_bed_path):
        print(f"[WARN] merged bed file not found, skipping: {merged_bed_path}")
        return
    if not os.path.isfile(sv_csv_path):
        print(f"[WARN] SV csv file not found, skipping: {sv_csv_path}")
        return

    # Load SV CSV
    s = pd.read_csv(sv_csv_path)

    if "id" not in s.columns:
        print(f"[ERROR] CSV {sv_csv_path} has no 'id' column. Skipping.")
        return
    if "elem_hits" not in s.columns:
        print(f"[ERROR] CSV {sv_csv_path} has no 'elem_hits' column. Skipping.")
        return
    if "pos" not in s.columns:
        print(f"[ERROR] CSV {sv_csv_path} has no 'pos' column (needed for distance). Skipping.")
        return

    # Prepare flank_motif_hits column
    if "flank_motif_hits" not in s.columns:
        s["flank_motif_hits"] = [ [] for _ in range(len(s)) ]
    else:
        s["flank_motif_hits"] = s["flank_motif_hits"].apply(
            lambda x: [] if (pd.isna(x) or x == "" or x == "[]") else (
                x if isinstance(x, list) else ast.literal_eval(x)
            )
        )

    # Parse elem_hits to list-of-dicts
    s["elem_hits_parsed"] = s["elem_hits"].apply(parse_elem_hits_cell)

    # Build an index from SV id -> row indices
    id_to_indices = defaultdict(list)
    for idx, sv_id in s["id"].items():
        id_to_indices[str(sv_id)].append(idx)

    # Load merged bed
    # Columns:
    # 0 sv_chrom, 1 sv_start, 2 sv_end, 3 sv_name,
    # 4 motif_chrom, 5 motif_start, 6 motif_end,
    # 7 motif_id, 8 p_value, 9 q_value, 10 score, 11 strand
    if os.path.getsize(merged_bed_path) == 0:
        print(f"[INFO] BED {merged_bed_path} is empty, nothing to add.")
        s.drop(columns=["elem_hits_parsed"], inplace=True)
        s.to_csv(out_csv_path, index=False)
        return

    mcols = [
        "sv_chrom", "sv_start", "sv_end", "sv_name",
        "motif_chrom", "motif_start", "motif_end",
        "motif_id", "p_value", "q_value", "score", "strand"
    ]
    merged = pd.read_csv(
        merged_bed_path,
        sep="\t",
        header=None,
        names=mcols,
        dtype={"sv_chrom": str, "motif_chrom": str}
    )

    if merged.empty:
        print(f"[INFO] BED {merged_bed_path} has no rows, nothing to add.")
        s.drop(columns=["elem_hits_parsed"], inplace=True)
        s.to_csv(out_csv_path, index=False)
        return

    # Separate flank vs regulatory hits
    is_flank = merged["sv_name"].str.contains("flank", na=False)
    flank_rows = merged[is_flank].copy()
    reg_rows   = merged[~is_flank].copy()

    # ------------------------------------------------------------------
    # Collect flank motif hits per SV id (without distance yet)
    # ------------------------------------------------------------------
    flank_hits_by_sv = defaultdict(list)

    for _, row in flank_rows.iterrows():
        name = row["sv_name"]
        sv_id = name.split("|")[0]               # before first '|'
        tail = name.split("|", 1)[1] if "|" in name else ""
        flank_side = parse_flank_side(tail)

        motif_start = int(row["motif_start"])
        motif_end   = int(row["motif_end"])
        motif_id    = str(row["motif_id"])

        # parse stats
        try:
            p_val = float(row["p_value"])
        except Exception:
            p_val = np.nan
        try:
            q_val = float(row["q_value"])
        except Exception:
            q_val = np.nan
        try:
            score = float(row["score"])
        except Exception:
            score = np.nan

        strand     = str(row["strand"])
        tf_short   = motif_id.split(".", 1)[0] if "." in motif_id else motif_id

        hit_dict = {
            "flank_side": flank_side,
            "start": motif_start,
            "end": motif_end,
            "TF": tf_short,
            "motif_id": motif_id,
            "score": score,
            "p_value": p_val,
            "q_value": q_val,
            "strand": strand,
            # distance_to_pos will be added per-row when attaching
        }
        flank_hits_by_sv[sv_id].append(hit_dict)

    # ------------------------------------------------------------------
    # Collect regulatory motif hits per (SV id, elem_id) (without distance yet)
    # ------------------------------------------------------------------
    elem_motif_hits = defaultdict(lambda: defaultdict(list))

    for _, row in reg_rows.iterrows():
        name = row["sv_name"]
        sv_id = name.split("|")[0]
        elem_id = extract_elem_id_from_name(name)
        if elem_id is None:
            continue

        motif_start = int(row["motif_start"])
        motif_end   = int(row["motif_end"])
        motif_id    = str(row["motif_id"])

        try:
            p_val = float(row["p_value"])
        except Exception:
            p_val = np.nan
        try:
            q_val = float(row["q_value"])
        except Exception:
            q_val = np.nan
        try:
            score = float(row["score"])
        except Exception:
            score = np.nan

        strand     = str(row["strand"])
        tf_short   = motif_id.split(".", 1)[0] if "." in motif_id else motif_id

        hit_dict = {
            "start": motif_start,
            "end": motif_end,
            "TF": tf_short,
            "motif_id": motif_id,
            "score": score,
            "p_value": p_val,
            "q_value": q_val,
            "strand": strand,
            # distance_to_pos will be added per-row when attaching
        }
        elem_motif_hits[sv_id][elem_id].append(hit_dict)

    # ------------------------------------------------------------------
    # Attach flank hits into s["flank_motif_hits"], add distance_to_pos
    # ------------------------------------------------------------------
    for sv_id, hits in flank_hits_by_sv.items():
        idxs = id_to_indices.get(sv_id, [])
        if not idxs:
            continue
        for idx in idxs:
            try:
                bp_pos = int(s.at[idx, "pos"])
            except Exception:
                bp_pos = None

            existing = s.at[idx, "flank_motif_hits"]
            if not isinstance(existing, list):
                if pd.isna(existing) or existing == "" or existing == "[]":
                    existing = []
                else:
                    try:
                        existing = ast.literal_eval(existing)
                    except Exception:
                        existing = []

            hits_with_dist = []
            for h in hits:
                h_copy = dict(h)
                if bp_pos is not None:
                    center = (h_copy["start"] + h_copy["end"]) // 2
                    h_copy["distance_to_pos"] = center - bp_pos
                else:
                    h_copy["distance_to_pos"] = None
                hits_with_dist.append(h_copy)

            existing.extend(hits_with_dist)
            s.at[idx, "flank_motif_hits"] = existing

    # ------------------------------------------------------------------
    # Attach regulatory motif hits into s["elem_hits_parsed"], add distance_to_pos
    # ------------------------------------------------------------------
    for sv_id, elem_dict in elem_motif_hits.items():
        idxs = id_to_indices.get(sv_id, [])
        if not idxs:
            continue

        for idx in idxs:
            try:
                bp_pos = int(s.at[idx, "pos"])
            except Exception:
                bp_pos = None

            elem_list = s.at[idx, "elem_hits_parsed"]
            if not isinstance(elem_list, list):
                continue

            for elem in elem_list:
                if not isinstance(elem, dict):
                    continue
                eid = elem.get("elem_id")
                if not eid:
                    continue
                hits_for_elem = elem_dict.get(eid)
                if not hits_for_elem:
                    continue

                existing = elem.get("motif_hits", [])
                if not isinstance(existing, list):
                    existing = []

                hits_with_dist = []
                for h in hits_for_elem:
                    h_copy = dict(h)
                    if bp_pos is not None:
                        center = (h_copy["start"] + h_copy["end"]) // 2
                        h_copy["distance_to_pos"] = center - bp_pos
                    else:
                        h_copy["distance_to_pos"] = None
                    hits_with_dist.append(h_copy)

                existing.extend(hits_with_dist)
                elem["motif_hits"] = existing

    # ------------------------------------------------------------------
    # Finalize and write
    # ------------------------------------------------------------------
    # convert parsed elem_hits back to stringified form
    s["elem_hits"] = s["elem_hits_parsed"].apply(lambda x: repr(x))
    s.drop(columns=["elem_hits_parsed"], inplace=True)

    # ensure flank_motif_hits is stringified as well
    s["flank_motif_hits"] = s["flank_motif_hits"].apply(
        lambda x: "[]" if (x is None or (isinstance(x, float) and np.isnan(x)) or x == []) else repr(x)
    )

    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)
    s.to_csv(out_csv_path, index=False)
    print(f"[INFO] Wrote updated CSV: {out_csv_path}")


def main():
    if len(sys.argv) != 4:
        print("Usage: python recombine_sv_fimo.py <merged_bed_dir> <sv_csv_dir> <output_csv_dir>", file=sys.stderr)
        sys.exit(1)

    merged_bed_dir = sys.argv[1]
    sv_csv_dir     = sys.argv[2]
    out_csv_dir    = sys.argv[3]

    os.makedirs(out_csv_dir, exist_ok=True)

    for fname in os.listdir(merged_bed_dir):
        if not fname.endswith("_fimo_merged.bed"):
            continue

        merged_bed_path = os.path.join(merged_bed_dir, fname)

        # Example:
        #   TCGA-A8-A09K-01A_strict_sv_set_flanks_and_overlaps_fimo_merged.bed
        # → base_csv = TCGA-A8-A09K-01A_strict_sv_set
        base = fname[:-len("_fimo_merged.bed")]  # remove suffix
        if base.endswith("_flanks_and_overlaps"):
            base_csv = base[:-len("_flanks_and_overlaps")]
        else:
            base_csv = base

        csv_name = base_csv + ".csv"
        sv_csv_path = os.path.join(sv_csv_dir, csv_name)
        out_csv_path = os.path.join(out_csv_dir, csv_name)

        recombine_one(merged_bed_path, sv_csv_path, out_csv_path)


if __name__ == "__main__":
    main()
