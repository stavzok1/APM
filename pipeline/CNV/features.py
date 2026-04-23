"""Basic CNV feature derivation."""

import numpy as np
import pandas as pd


def _classify_cn_state(cn):
    if pd.isna(cn):
        return "NA"
    cn = float(cn)
    if cn == 0: return "deep_del"
    elif cn == 1: return "loss"
    elif cn == 2: return "neutral"
    elif cn == 3: return "gain"
    elif cn >= 4: return "amp"
    return "other"


def add_basic_cnv_features(cnv_df: pd.DataFrame) -> pd.DataFrame:
    """Derive CNV features and attach atlas-aligned lowercase aliases.

    The raw ASCAT3 GDC tables use ``Chromosome/Start/End/Copy_Number`` etc.
    The Data Structure Atlas uses lowercase ``chrom/start/end/num_probes/
    segment_mean/segment_length``. We keep the original columns intact (other
    code relies on them) and add the atlas aliases plus an approximate
    ``segment_mean = log2((cn_total + eps) / 2)`` so the table is usable by
    tools that expect log-copy-ratio segmenter output.
    """
    df = cnv_df.copy()
    df["cn_total"] = df["Copy_Number"].astype(float)
    df["cn_major"] = df["Major_Copy_Number"].astype(float)
    df["cn_minor"] = df["Minor_Copy_Number"].astype(float)
    df["loh_flag"] = (df["cn_minor"] == 0).astype(int)
    df["cn_state"] = df["cn_total"].map(_classify_cn_state)
    df["segment_len"] = df["End"] - df["Start"] + 1

    if "chrom" not in df.columns:
        df["chrom"] = df["Chromosome"].astype(str)
    if "start" not in df.columns:
        df["start"] = df["Start"].astype("Int64")
    if "end" not in df.columns:
        df["end"] = df["End"].astype("Int64")
    if "segment_length" not in df.columns:
        df["segment_length"] = df["segment_len"]
    if "num_probes" not in df.columns:
        # ASCAT/GDC exports vary: sometimes probe count exists, sometimes not.
        # Populate when we can; otherwise keep as NA (atlas-aligned column still exists).
        src = None
        for cand in (
            "Num_Probes",
            "NumProbes",
            "num_probes",
            "N_PROBES",
            "n_probes",
            "probes",
            "Nprobes",
            "nProbes",
            "num.mark",
            "num_mark",
            "MARKERS",
            "markers",
        ):
            if cand in df.columns:
                src = cand
                break
        if src is not None:
            df["num_probes"] = pd.to_numeric(df[src], errors="coerce").astype("Int64")
        else:
            df["num_probes"] = pd.array([pd.NA] * len(df), dtype="Int64")
    if "segment_mean" not in df.columns:
        eps = 1e-3
        df["segment_mean"] = np.log2((df["cn_total"].astype(float) + eps) / 2.0)

    return df
