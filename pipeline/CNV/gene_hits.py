"""CNV segment annotation with gene-level hits."""
import pandas as pd
from .geometry import compute_signed_distance_with_overlap, classify_gene_hit

def _annotate_segs(cnv_df, feat_gene, feat_all, col, sub=True, pu=2000, pd_=500, pw=5000):
    df = cnv_df.copy()
    df["Chromosome"] = df["Chromosome"].astype(str)
    g = feat_gene.copy(); g["chrom"] = g["chrom"].astype(str)
    af = feat_all.copy() if feat_all is not None else None
    if af is not None: af["chrom"] = af["chrom"].astype(str)
    hc = [[] for _ in range(len(df))]
    for ch, idxs in df.groupby("Chromosome").indices.items():
        segs = df.loc[list(idxs)]
        gc = g[g["chrom"] == ch]
        if gc.empty: continue
        fc = af[af["chrom"] == ch] if af is not None else None
        gc = gc.sort_values(["start", "end"])
        for si, sr in segs.iterrows():
            ss, se = int(sr["Start"]), int(sr["End"])
            ca = gc[(gc["end"] >= ss - pw) & (gc["start"] <= se + pw)]
            if ca.empty: continue
            sh = []
            for _, r in ca.iterrows():
                gs_, ge_ = int(r["start"]), int(r["end"])
                sd, _, _ = compute_signed_distance_with_overlap(ss, se, gs_, ge_, r["strand"])
                gf = fc[fc["gene_name"]==r["gene_name"]] if sub and fc is not None else None
                ann = classify_gene_hit(ss, se, gs_, ge_, r["strand"], sd, gf, pu, pd_, pw)
                if ann["region_hit"] == "intergenic": continue
                hit = {"gene_name": r["gene_name"], "gene_id": r.get("gene_id"), "strand": r["strand"], "signed_dist": sd, **ann}
                # CNV miRNA loci table can optionally carry mature-arm mappings / canonical names.
                for extra in (
                    "mature_names",
                    "mature_accessions",
                    "mirbase_mature_id",
                    "pre_gene_name",
                    "pre_gene_id",
                ):
                    if extra in r and pd.notna(r[extra]):
                        hit[extra] = r[extra]
                sh.append(hit)
            hc[si] = sh
    df[col] = hc
    return df

def annotate_cnv_with_gene_hits(cnv_df, genes_all_features_df, pu=2000, pd_=500, pw=5000):
    gg = genes_all_features_df[genes_all_features_df["feature"]=="gene"].copy()
    return _annotate_segs(cnv_df, gg, genes_all_features_df, "gene_hits", True, pu, pd_, pw)

def annotate_cnv_with_lncrna_hits(cnv_df, lnc_all, pu=2000, pd_=500, pw=5000):
    """
    ``lnc_all`` may be GENCODE-style (many ``feature`` rows per gene) or gene-centric
    (one row per lncRNA, e.g. ``lncRNAs_with_genes_1000000bp.csv``).
    """
    if "feature" in lnc_all.columns:
        lg = lnc_all[lnc_all["feature"] == "gene"].copy()
        return _annotate_segs(cnv_df, lg, lnc_all, "lncRNAs_hits", True, pu, pd_, pw)
    lg = lnc_all.copy()
    return _annotate_segs(cnv_df, lg, None, "lncRNAs_hits", False, pu, pd_, pw)

def annotate_cnv_with_mirna_hits(cnv_df, mirna_df, pu=2000, pd_=500, pw=5000):
    return _annotate_segs(cnv_df, mirna_df, None, "mir_hits", False, pu, pd_, pw)
