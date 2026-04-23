"""CNV segment annotation with regulatory element (cCRE) hits."""
import pandas as pd
from .geometry import compute_signed_distance_interval, classify_elem_hit

def annotate_cnv_with_elem_hits(cnv_df, elements_df, pw=5000):
    df = cnv_df.copy(); df["Chromosome"] = df["Chromosome"].astype(str)
    el = elements_df.copy(); el["chrom"] = el["chrom"].astype(str)
    hc = [[] for _ in range(len(df))]
    for ch, idxs in df.groupby("Chromosome").indices.items():
        segs = df.loc[list(idxs)]
        ec = el[el["chrom"] == ch]
        if ec.empty: continue
        ec = ec.sort_values(["start", "end"])
        for si, sr in segs.iterrows():
            ss, se = int(sr["Start"]), int(sr["End"])
            ca = ec[(ec["end"]>=ss-pw)&(ec["start"]<=se+pw)]
            if ca.empty: continue
            sh = []
            for _, e in ca.iterrows():
                es_, ee_ = int(e["start"]), int(e["end"])
                sd = compute_signed_distance_interval(ss, se, es_, ee_)
                ann = classify_elem_hit(ss, se, es_, ee_, sd, pw)
                if ann["region_hit"] == "distal": continue
                sh.append({"element_id": e["elem_id"],
                    "element_type": e.get("element_type"),
                    "signed_dist": sd, **ann})
            hc[si] = sh
    df["elem_hits"] = hc
    return df
