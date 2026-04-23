import vcfpy
import pandas as pd
import numpy as np
import os, sys
import re
def load_manta_sv_vcf(vcf_path: str) -> pd.DataFrame:
    reader = vcfpy.Reader.from_path(vcf_path)
    samples = reader.header.samples.names  # e.g. ["NORMAL", "TUMOR"]
    print("Samples in VCF:", samples)

    # Try to identify normal/tumor explicitly, otherwise fallback to order
    normal_sample = None
    tumor_sample = None
    for s in samples:
        su = s.upper()
        if "NORMAL" in su:
            normal_sample = s
        if "TUMOR" in su:
            tumor_sample = s
    if normal_sample is None and len(samples) >= 1:
        normal_sample = samples[0]
    if tumor_sample is None and len(samples) >= 2:
        tumor_sample = samples[-1]

    print("Assuming normal sample:", normal_sample)
    print("Assuming tumor sample :", tumor_sample)

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

        # INFO fields: take everything, flatten single-element lists
        for key, value in rec.INFO.items():
            if isinstance(value, list) and len(value) == 1:
                row[key] = value[0]
            else:
                row[key] = value

        # FORMAT / calls: PR and SR for each sample
        for sample in samples:
            call = rec.call_for_sample[sample]   # dict: sample -> Call
            pr = call.data.get("PR")
            sr = call.data.get("SR")

            row[f"{sample}_PR_ref"] = pr[0] if pr else np.nan
            row[f"{sample}_PR_alt"] = pr[1] if pr else np.nan
            row[f"{sample}_SR_ref"] = sr[0] if sr else np.nan
            row[f"{sample}_SR_alt"] = sr[1] if sr else np.nan

        rows.append(row)

    df = pd.DataFrame(rows)

    # Ensure some key numeric fields are numeric
    numeric_cols = [
        "SVLEN", "SOMATICSCORE",
        f"{normal_sample}_PR_ref", f"{normal_sample}_PR_alt",
        f"{normal_sample}_SR_ref", f"{normal_sample}_SR_alt",
        f"{tumor_sample}_PR_ref",  f"{tumor_sample}_PR_alt",
        f"{tumor_sample}_SR_ref",  f"{tumor_sample}_SR_alt",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Convenience: total alt counts
    df["normal_alt"] = (
        df.get(f"{normal_sample}_PR_alt", 0).fillna(0) +
        df.get(f"{normal_sample}_SR_alt", 0).fillna(0)
    )
    df["tumor_alt"] = (
        df.get(f"{tumor_sample}_PR_alt", 0).fillna(0) +
        df.get(f"{tumor_sample}_SR_alt", 0).fillna(0)
    )

    df["normal_sr_alt"] = df.get(f"{normal_sample}_SR_alt", 0).fillna(0)
    df["tumor_sr_alt"]  = df.get(f"{tumor_sample}_SR_alt", 0).fillna(0)
    df["tumor_pr_alt"]  = df.get(f"{tumor_sample}_PR_alt", 0).fillna(0)

    return df, normal_sample, tumor_sample


def get_strict_sv_set(df: pd.DataFrame) -> pd.DataFrame:
    allowed_types = {"DEL", "DUP", "INV", "INS", "BND"}

    # Handle missing columns safely
    has_svlen = "SVLEN" in df.columns
    has_score = "SOMATICSCORE" in df.columns

    q = (
        (df["filter"] == "PASS") &
        (df["SVTYPE"].isin(allowed_types)) &
        (df["tumor_sr_alt"] >= 2) &
        (df["tumor_alt"] >= 8) &
        (df["normal_alt"] <= 1)
    )

    if has_score:
        q &= df["SOMATICSCORE"].fillna(0) >= 25

    # if has_svlen:
    #     # SVLEN can be negative for some callers (orientation); take abs
    #     q &= df["SVLEN"].abs() >= 200
    return df[q].copy()


def get_lenient_sv_set(df: pd.DataFrame) -> pd.DataFrame:
    allowed_types = {"DEL", "DUP", "INV", "INS", "BND"}

    has_svlen = "SVLEN" in df.columns
    has_score = "SOMATICSCORE" in df.columns

    q = (
        (df["filter"] == "PASS") &
        (df["SVTYPE"].isin(allowed_types)) &
        ((df["tumor_sr_alt"] >= 1) | (df["tumor_pr_alt"] >= 5)) &
        (df["normal_alt"] <= 1)
    )

    if has_score:
        q &= df["SOMATICSCORE"].fillna(0) >= 15

    # if has_svlen:
    #     abs_len = df["SVLEN"].abs()
    #     q &= (abs_len >= 50) & (abs_len <= 50000)

    return df[q].copy()


def compute_signed_distance_with_overlap(del_start, del_end, gene_start, gene_end, strand):
    """
    Compute signed distance between deletion [del_start, del_end]
    and gene [gene_start, gene_end], and also return overlap positions.

    Returns:
        signed_distance, overlap_start, overlap_end
        - signed_distance = 0 if overlapping
        - overlap_start/end = genomic coordinates of overlap, else None
    """

    # --- Check overlap ---
    if del_end >= gene_start and del_start <= gene_end:
        overlap_start = max(del_start, gene_start)
        overlap_end = min(del_end, gene_end)
        return 0, overlap_start, overlap_end

    # --- No overlap ---
    overlap_start = None
    overlap_end = None

    # deletion completely before gene
    if del_end < gene_start:
        d = gene_start - del_end
        if strand == '+':
            return -d, overlap_start, overlap_end
        else:
            return +d, overlap_start, overlap_end

    # deletion completely after gene
    if del_start > gene_end:
        d = del_start - gene_end
        if strand == '+':
            return +d, overlap_start, overlap_end
        else:
            return -d, overlap_start, overlap_end

    # Should not reach here
    return 0, overlap_start, overlap_end



def classify_span_sv_gene_hit(del_start: int,
                          del_end: int,
                          gene_start: int,
                          gene_end: int,
                          strand: str,
                          signed_dist: int,
                          gene_features: pd.DataFrame,
                          promoter_up: int = 2000,
                          promoter_down: int = 500,
                          proximal_window: int = 5000) -> dict:
    """
    Return a dict with region-level annotations for a DEL–gene pair.
    gene_features: rows from `genes` for this gene (exon/CDS/UTR etc.).
    """

    # --- gene body overlap ---
    inter_body_start = max(del_start, gene_start)
    inter_body_end = min(del_end, gene_end)
    if inter_body_end >= inter_body_start:
        overlap_bp = inter_body_end - inter_body_start + 1
        overlap_percent = (overlap_bp / (gene_end - gene_start + 1)) * 100
        gene_body_flag = 1
    else:
        overlap_bp = 0
        overlap_percent = 0
        gene_body_flag = 0

    # --- promoter region (strand-aware) ---
    # TSS = gene_start on + strand, gene_end on - strand
    if strand == "+":
        tss = gene_start
        prom_start = tss - promoter_up
        prom_end = tss + promoter_down
    else:
        tss = gene_end
        prom_start = tss - promoter_down
        prom_end = tss + promoter_up

    inter_prom_start = max(del_start, prom_start)
    inter_prom_end = min(del_end, prom_end)
    promoter_flag = 1 if inter_prom_end >= inter_prom_start else 0

    # --- detailed gene-structure hits from `genes` ---
    hits_exon_like = False
    hits_utr = False
    hits_stop_codon = False
    hits_start_codon = False
    transcript_id = None
    transcript_type = None

    for _, feat in gene_features.iterrows():
        f_start = int(feat["start"])
        f_end = int(feat["end"])
        if del_end < f_start or del_start > f_end:
            continue  # no overlap
        ftype = feat["feature"]
        if ftype in ("exon", "CDS", "start_codon", "stop_codon", "transcript"):
            hits_exon_like = True
            if ftype == "stop_codon":
                hits_stop_codon = True
                hits_start_codon = True
            elif ftype == "transcript":
                 transcript_id = feat["transcript_id"]
                 transcript_type = feat["transcript_type"]
        elif ftype == "UTR":
            hits_utr = True


    exon_flag = 1 if (hits_exon_like or hits_utr) else 0
    intron_only_flag = 1 if (gene_body_flag and not exon_flag) else 0
    stop_codon_flag = 1 if hits_stop_codon else 0
    start_codon_flag = 1 if hits_start_codon else 0

    # --- classify region_hit (coarse) ---
    region_hit = ""
    if promoter_flag or overlap_bp > 0:
        if promoter_flag:
            region_hit += "promoter,"
        if overlap_bp > 0:
            if exon_flag:
                region_hit += "exon_or_UTR,"
            else:
                region_hit += "intron,"
    else:
        # no gene body overlap → maybe proximal upstream/downstream
        if abs(signed_dist) <= proximal_window:
            if signed_dist < 0:
                region_hit += "upstream_5kb,"
            else:
                region_hit += "downstream_5kb,"
        else:
            region_hit += "intergenic"

   
    upstream_5kb_flag = int("upstream_5kb" in region_hit)
    downstream_5kb_flag = int("downstream_5kb" in region_hit)

    # For deletions, the "hit_side" is basically the span
    hit_side = "span"

    return {
        "overlap_bp": overlap_bp,
        "overlap_percent": overlap_percent,
        "promoter_flag": promoter_flag,
        "gene_body_flag": gene_body_flag,
        "exon_flag": exon_flag,
        "intron_only_flag": intron_only_flag,
        "upstream_5kb_flag": upstream_5kb_flag,
        "downstream_5kb_flag": downstream_5kb_flag,
        "region_hit": region_hit,
        "hit_side": hit_side,
        "stop_codon_flag": stop_codon_flag,
        "start_codon_flag": start_codon_flag,
        "transcript_id" : transcript_id,
        "transcript_type": transcript_type
    }


def classify_span_sv_elem_hit(sv_start: int,
                              sv_end: int,
                              elem_start: int,
                              elem_end: int,
                              signed_dist: int,
                              proximal_window: int = 5000) -> dict:
    """
    Return a dict with region-level annotations for a span-SV–element pair.

    This is simpler than the gene version: no strand, no promoter/exon logic.
    """

    # --- element overlap ---
    inter_start = max(sv_start, elem_start)
    inter_end = min(sv_end, elem_end)
    if inter_end >= inter_start:
        overlap_bp = inter_end - inter_start + 1
        overlap_percent = (overlap_bp / (elem_end - elem_start + 1)) * 100
        overlaps_flag = 1
    else:
        overlap_bp = 0
        overlap_percent = 0
        overlaps_flag = 0

    # --- classify coarse region_hit ---
    if overlaps_flag:
        region_hit = "overlaps"
    else:
        if abs(signed_dist) <= proximal_window:
            # You can keep "upstream/downstream" semantics from signed_dist
            if signed_dist < 0:
                region_hit = "proximal_upstream"
            elif signed_dist > 0:
                region_hit = "proximal_downstream"
            else:
                region_hit = "proximal"  # degenerate case
        else:
            region_hit = "distal"

    proximal_flag = int(region_hit.startswith("proximal"))
    distal_flag = int(region_hit == "distal")

    # For DEL/DUP-like span events, hit_side is "span"
    hit_side = "span"

    return {
        "overlap_bp": overlap_bp,
        "overlap_percent": overlap_percent,
        "overlaps_flag": overlaps_flag,
        "region_hit": region_hit,
        "proximal_flag": proximal_flag,
        "distal_flag": distal_flag,
        "hit_side": hit_side,
    }



def map_dels_to_genes(strict_sv_set: pd.DataFrame,
                      genes: pd.DataFrame,
                      genes_only: pd.DataFrame, is_lncRNAs: bool,
                      window: int = 1_000_000) -> pd.DataFrame:
    """
    For DEL SVs in strict_sv_set, add a 'gene_hits' column:
    a list of dicts, one per gene hit.

    Each dict roughly:
    {
      "gene_name": ...,
      "gene_id": ...,
      "strand": ...,
      "signed_dist": ...,
      "overlap_bp": ...,
      "region_hit": "promoter" / "exon_or_UTR" / "intron" /
                    "upstream_5kb" / "downstream_5kb" / "intergenic",
      "promoter_flag": 0/1,
      "gene_body_flag": 0/1,
      "exon_flag": 0/1,
      "intron_only_flag": 0/1,
      "upstream_5kb_flag": 0/1,
      "downstream_5kb_flag": 0/1,
      "hit_side": "span",
    }
    """
    df = strict_sv_set.copy()

    # Only deletions
    dels = df[df["SVTYPE"] == "DEL"].copy()

    # new column as list-of-dicts
    gene_hits_col = []

    # Pre-slice features by chromosome once to avoid repeated full scans
    # (genes_only and genes are assumed to have 'chrom' and 'gene_name')
    for idx, sv in dels.iterrows():
        chrom = sv["chrom"]
        del_start = int(sv["pos"])
        del_end = int(sv["END"])

        # candidate genes: same chromosome, roughly within the window
        g_chr = genes_only[genes_only["chrom"] == chrom].copy()

        # coarse window filter
        g_chr = g_chr[
            (g_chr["end"] >= del_start - window) &
            (g_chr["start"] <= del_end + window)
        ]

        if g_chr.empty:
            gene_hits_col.append([])  # no genes nearby
            continue

        # full features for just those genes on this chromosome
        genes_features_chr = genes[genes["gene_name"].isin(g_chr["gene_name"])]

        connections = []

        for _, g in g_chr.iterrows():
            gene_start = int(g["start"])
            gene_end = int(g["end"])
            strand = g.get("strand", "+")
            gene_name = g["gene_name"]
            gene_id = g.get("gene_id", None)

            # signed distance, strand-aware
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                del_start, del_end, gene_start, gene_end, strand
            )

            # only keep genes within the window by absolute distance
            if abs(signed_dist) > window:
                continue

            # all features for this specific gene
            gene_feats = genes_features_chr[
                genes_features_chr["gene_name"] == gene_name
            ]

            region_info = classify_span_sv_gene_hit(
                del_start=del_start,
                del_end=del_end,
                gene_start=gene_start,
                gene_end=gene_end,
                strand=strand,
                signed_dist=signed_dist,
                gene_features=gene_feats,
            )

            gene_hit = {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(gene_hit)

        gene_hits_col.append(connections)
    if not is_lncRNAs:
        # attach the list-of-dicts back into the DEL subset
        dels["gene_hits"] = gene_hits_col

        # and write that column back into the full strict_sv_set
        df.loc[dels.index, "gene_hits"] = dels["gene_hits"]
    else:
        dels["lncRNA_hits"] = gene_hits_col
        df.loc[dels.index, "lncRNA_hits"] = dels["lncRNA_hits"]

    return df


def map_dels_to_elements(strict_sv_set: pd.DataFrame,
                         elements: pd.DataFrame,
                         window: int = 1_000_000) -> pd.DataFrame:
    """
    For DEL SVs in strict_sv_set, add an 'elem_hits' column:
    a list of dicts, one per element hit.

    Each dict roughly:
    {
      "elem_id": ...,
      "elem_type": ...,
      "chrom": ...,
      "elem_start": ...,
      "elem_end": ...,
      "signed_dist": ...,
      "overlap_bp": ...,
      "region_hit": "overlaps" / "proximal_upstream" /
                    "proximal_downstream" / "distal",
      "overlaps_flag": 0/1,
      "proximal_flag": 0/1,
      "distal_flag": 0/1,
      "hit_side": "span",
    }
    """
    df = strict_sv_set.copy()

    # Only deletions
    dels = df[df["SVTYPE"] == "DEL"].copy()

    elem_hits_col = []

    for idx, sv in dels.iterrows():
        chrom = sv["chrom"]
        del_start = int(sv["pos"])
        del_end = int(sv["END"])

        # candidate elements: same chromosome, roughly within the window
        e_chr = elements[elements["chrom"] == chrom].copy()
        e_chr = e_chr[
            (e_chr["end"]   >= del_start - window) &
            (e_chr["start"] <= del_end   + window)
        ]

        if e_chr.empty:
            elem_hits_col.append([])  # no elements nearby
            continue

        connections = []

        for _, e in e_chr.iterrows():
            elem_start = int(e["start"])
            elem_end   = int(e["end"])
            elem_id    = e["elem_id"]
            elem_type  = e.get("type", None)

            # signed distance, strand-agnostic: just pass "+"
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                del_start, del_end, elem_start, elem_end, strand="+"
            )

            # only keep elements within ±window bp by absolute distance
            if abs(signed_dist) > window:
                continue

            region_info = classify_span_sv_elem_hit(
                sv_start=del_start,
                sv_end=del_end,
                elem_start=elem_start,
                elem_end=elem_end,
                signed_dist=signed_dist,
            )

            elem_hit = {
                "elem_id": elem_id,
                "elem_type": elem_type,
                "chrom": chrom,
                "elem_start": elem_start,
                "elem_end": elem_end,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(elem_hit)

        elem_hits_col.append(connections)

    # attach the list-of-dicts back into the DEL subset
    dels["elem_hits"] = elem_hits_col

    # and write that column back into the full strict_sv_set
    df.loc[dels.index, "elem_hits"] = dels["elem_hits"]

    return df



def map_dups_to_genes(strict_sv_set: pd.DataFrame,
                      genes: pd.DataFrame,
                      genes_only: pd.DataFrame, is_lncRNAs: bool,
                      window: int = 1_000_000) -> pd.DataFrame:
    """
    For DUP SVs in strict_sv_set, add a 'gene_hits' column:
    a list of dicts, one per gene hit.

    Each dict roughly:
    {
      "gene_name": ...,
      "gene_id": ...,
      "strand": ...,
      "signed_dist": ...,
      "overlap_bp": ...,
      "region_hit": "promoter" / "exon_or_UTR" / "intron" /
                    "upstream_5kb" / "downstream_5kb" / "intergenic",
      "promoter_flag": 0/1,
      "gene_body_flag": 0/1,
      "exon_flag": 0/1,
      "intron_only_flag": 0/1,
      "upstream_5kb_flag": 0/1,
      "downstream_5kb_flag": 0/1,
      "hit_side": "span",
      # plus any extra flags you added in classify_del_gene_hit,
      # e.g. start_codon_flag, stop_codon_flag, etc.
    }
    """
    df = strict_sv_set.copy()

    # Only duplications
    dups = df[df["SVTYPE"] == "DUP"].copy()

    # new column as list-of-dicts
    gene_hits_col = []

    for idx, sv in dups.iterrows():
        chrom = sv["chrom"]
        dup_start = int(sv["pos"])

        # figure out end; fallback to SVLEN if END is missing
        if "END" in sv and not pd.isna(sv["END"]):
            dup_end = int(sv["END"])
        elif "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            dup_end = dup_start + abs(int(sv["SVLEN"]))
        else:
            # if we don't know the length, treat as point
            dup_end = dup_start

        # candidate genes: same chromosome, within coarse window
        g_chr = genes_only[genes_only["chrom"] == chrom].copy()
        g_chr = g_chr[
            (g_chr["end"]   >= dup_start - window) &
            (g_chr["start"] <= dup_end   + window)
        ]

        if g_chr.empty:
            gene_hits_col.append([])  # no genes nearby
            continue

        # full features for just those genes on this chromosome
        genes_features_chr = genes[genes["gene_name"].isin(g_chr["gene_name"])]

        connections = []

        for _, g in g_chr.iterrows():
            gene_start = int(g["start"])
            gene_end   = int(g["end"])
            strand     = g.get("strand", "+")
            gene_name  = g["gene_name"]
            gene_id    = g.get("gene_id", None)

            # signed distance, strand-aware
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                dup_start, dup_end, gene_start, gene_end, strand
            )

            # keep genes within ±window bp by absolute distance
            if abs(signed_dist) > window:
                continue

            # all features for this specific gene
            gene_feats = genes_features_chr[
                genes_features_chr["gene_name"] == gene_name
            ]

            # same classifier as for DELs; it just sees an interval [start, end]
            region_info = classify_span_sv_gene_hit(
                del_start=dup_start,
                del_end=dup_end,
                gene_start=gene_start,
                gene_end=gene_end,
                strand=strand,
                signed_dist=signed_dist,
                gene_features=gene_feats,
            )

            gene_hit = {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(gene_hit)

        gene_hits_col.append(connections)

    if not is_lncRNAs:
        # attach the list-of-dicts back into the DUP subset
        dups["gene_hits"] = gene_hits_col

        # and write that column back into the full strict_sv_set
        df.loc[dups.index, "gene_hits"] = dups["gene_hits"]
    else:
        dups["lncRNA_hits"] = gene_hits_col
        df.loc[dups.index, "lncRNA_hits"] = dups["lncRNA_hits"]

    return df


def map_dups_to_elements(strict_sv_set: pd.DataFrame,
                         elements: pd.DataFrame,
                         window: int = 1_000_000) -> pd.DataFrame:
    """
    For DUP SVs in strict_sv_set, add an 'elem_hits' column:
    a list of dicts, one per regulatory-element hit.

    Each dict roughly:
    {
      "elem_id": ...,
      "elem_type": ...,
      "chrom": ...,
      "elem_start": ...,
      "elem_end": ...,
      "signed_dist": ...,
      "overlap_bp": ...,
      "region_hit": "overlaps" / "proximal_upstream" /
                    "proximal_downstream" / "distal",
      "overlaps_flag": 0/1,
      "proximal_flag": 0/1,
      "distal_flag": 0/1,
      "hit_side": "span",
    }
    """
    df = strict_sv_set.copy()

    # Only duplications
    dups = df[df["SVTYPE"] == "DUP"].copy()

    # new column as list-of-dicts
    elem_hits_col = []

    for idx, sv in dups.iterrows():
        chrom = sv["chrom"]
        dup_start = int(sv["pos"])

        # figure out end; fallback to SVLEN if END is missing
        if "END" in sv and not pd.isna(sv["END"]):
            dup_end = int(sv["END"])
        elif "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            dup_end = dup_start + abs(int(sv["SVLEN"]))
        else:
            # if we don't know the length, treat as point
            dup_end = dup_start

        # candidate elements: same chromosome, within coarse window
        e_chr = elements[elements["chrom"] == chrom].copy()
        e_chr = e_chr[
            (e_chr["end"]   >= dup_start - window) &
            (e_chr["start"] <= dup_end   + window)
        ]

        if e_chr.empty:
            elem_hits_col.append([])  # no elements nearby
            continue

        connections = []

        for _, e in e_chr.iterrows():
            elem_start = int(e["start"])
            elem_end   = int(e["end"])
            elem_id    = e["elem_id"]
            elem_type  = e.get("type", None)

            # signed distance, strand-agnostic: just pass "+"
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                dup_start, dup_end, elem_start, elem_end, strand="+"
            )

            # keep elements within ±window bp by absolute distance
            if abs(signed_dist) > window:
                continue

            region_info = classify_span_sv_elem_hit(
                sv_start=dup_start,
                sv_end=dup_end,
                elem_start=elem_start,
                elem_end=elem_end,
                signed_dist=signed_dist,
            )

            elem_hit = {
                "elem_id": elem_id,
                "elem_type": elem_type,
                "chrom": chrom,
                "elem_start": elem_start,
                "elem_end": elem_end,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(elem_hit)

        elem_hits_col.append(connections)

    # attach the list-of-dicts back into the DUP subset
    dups["elem_hits"] = elem_hits_col

    # and write that column back into the full strict_sv_set
    df.loc[dups.index, "elem_hits"] = dups["elem_hits"]

    return df



def map_ins_to_genes(strict_sv_set: pd.DataFrame,
                     genes: pd.DataFrame,
                     genes_only: pd.DataFrame, is_lncRNAs: bool,
                     window: int = 1_000_000) -> pd.DataFrame:
    """
    For INS events in strict_sv_set, add a 'gene_hits' column:
    a list of dicts, one per gene hit.

    Each dict roughly:
    {
      "gene_name": ...,
      "gene_id": ...,
      "strand": ...,
      "signed_dist": ...,
      "overlap_bp": ...,
      "region_hit": "promoter" / "exon_or_UTR" / "intron" /
                    "upstream_5kb" / "downstream_5kb" / "intergenic",
      "promoter_flag": 0/1,
      "gene_body_flag": 0/1,
      "exon_flag": 0/1,
      "intron_only_flag": 0/1,
      "upstream_5kb_flag": 0/1,
      "downstream_5kb_flag": 0/1,
      "hit_side": "point"  # overridden from classifier
      # plus any extra flags you added (start/stop codon etc.)
    }
    """
    df = strict_sv_set.copy()

    # Only insertions
    ins_df = df[df["SVTYPE"] == "INS"].copy()

    gene_hits_col = []

    for idx, sv in ins_df.iterrows():
        chrom = sv["chrom"]
        ins_pos = int(sv["pos"])

        # SVLEN is not needed for coordinates, but we can read it if present
        ins_len = None
        if "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            ins_len = int(abs(sv["SVLEN"]))

        # Genes on same chromosome, coarse window
        g_chr = genes_only[genes_only["chrom"] == chrom].copy()
        g_chr = g_chr[
            (g_chr["start"] <= ins_pos + window) &
            (g_chr["end"]   >= ins_pos - window)
        ]

        if g_chr.empty:
            gene_hits_col.append([])  # no nearby genes
            continue

        # Full feature rows for these genes
        genes_features_chr = genes[genes["gene_name"].isin(g_chr["gene_name"])]

        connections = []

        for _, g in g_chr.iterrows():
            gene_start = int(g["start"])
            gene_end   = int(g["end"])
            strand     = g.get("strand", "+")
            gene_name  = g["gene_name"]
            gene_id    = g.get("gene_id", None)

            # Signed distance: point vs interval, treat INS as [ins_pos, ins_pos]
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                ins_pos, ins_pos, gene_start, gene_end, strand
            )

            # Keep genes within ±window bp
            if abs(signed_dist) > window:
                continue

            # All features for this gene
            gene_feats = genes_features_chr[
                genes_features_chr["gene_name"] == gene_name
            ]

            # Reuse the span classifier, with start=end=ins_pos
            region_info = classify_span_sv_gene_hit(
                del_start=ins_pos,
                del_end=ins_pos,
                gene_start=gene_start,
                gene_end=gene_end,
                strand=strand,
                signed_dist=signed_dist,
                gene_features=gene_feats,
            )

            # For INS, mark hit_side as "point"
            region_info["hit_side"] = "point"

            gene_hit = {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(gene_hit)

        gene_hits_col.append(connections)

    # Attach back
    if not is_lncRNAs:
        ins_df["gene_hits"] = gene_hits_col
        df.loc[ins_df.index, "gene_hits"] = ins_df["gene_hits"]
    else:
        ins_df["lncRNA_hits"] = gene_hits_col
        df.loc[ins_df.index, "lncRNA_hits"] = ins_df["lncRNA_hits"]

    return df


def map_ins_to_elements(strict_sv_set: pd.DataFrame,
                        elements: pd.DataFrame,
                        window: int = 1_000_000) -> pd.DataFrame:
    """
    For INS events in strict_sv_set, add an 'elem_hits' column:
    a list of dicts, one per regulatory-element hit.

    Each dict roughly:
    {
      "elem_id": ...,
      "elem_type": ...,
      "chrom": ...,
      "elem_start": ...,
      "elem_end": ...,
      "signed_dist": ...,
      "overlap_bp": ...,
      "region_hit": "overlaps" / "proximal_upstream" /
                    "proximal_downstream" / "distal",
      "overlaps_flag": 0/1,
      "proximal_flag": 0/1,
      "distal_flag": 0/1,
      "hit_side": "point",  # override from classifier for INS
    }
    """
    df = strict_sv_set.copy()

    # Only insertions
    ins_df = df[df["SVTYPE"] == "INS"].copy()

    elem_hits_col = []

    for idx, sv in ins_df.iterrows():
        chrom = sv["chrom"]
        ins_pos = int(sv["pos"])

        # SVLEN is not needed for coordinates, but we can read it if present
        ins_len = None
        if "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            ins_len = int(abs(sv["SVLEN"]))

        # Elements on same chromosome, coarse window
        e_chr = elements[elements["chrom"] == chrom].copy()
        e_chr = e_chr[
            (e_chr["start"] <= ins_pos + window) &
            (e_chr["end"]   >= ins_pos - window)
        ]

        if e_chr.empty:
            elem_hits_col.append([])  # no nearby elements
            continue

        connections = []

        for _, e in e_chr.iterrows():
            elem_start = int(e["start"])
            elem_end   = int(e["end"])
            elem_id    = e["elem_id"]
            elem_type  = e.get("type", None)

            # Signed distance: point vs interval, treat INS as [ins_pos, ins_pos]
            # strand is irrelevant for elements → pass "+"
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                ins_pos, ins_pos, elem_start, elem_end, strand="+"
            )

            # Keep elements within ±window bp
            if abs(signed_dist) > window:
                continue

            # Reuse the span classifier, with start=end=ins_pos
            region_info = classify_span_sv_elem_hit(
                sv_start=ins_pos,
                sv_end=ins_pos,
                elem_start=elem_start,
                elem_end=elem_end,
                signed_dist=signed_dist,
            )

            # For INS, mark hit_side as "point"
            region_info["hit_side"] = "point"

            elem_hit = {
                "elem_id": elem_id,
                "elem_type": elem_type,
                "chrom": chrom,
                "elem_start": elem_start,
                "elem_end": elem_end,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(elem_hit)

        elem_hits_col.append(connections)

    # Attach back
    ins_df["elem_hits"] = elem_hits_col
    df.loc[ins_df.index, "elem_hits"] = ins_df["elem_hits"]

    return df


def get_bnd_remote_coords(alt_value):
    """
    Extract remote chromosome + position from a Manta BND ALT.

    Assumes the repr looks like:
        BreakEnd('chr5', 16982822, '+', '+', 'T', True)

    Returns:
        (chrom2, pos2) or (None, None).
    """
    if alt_value is None or (isinstance(alt_value, float) and np.isnan(alt_value)):
        return None, None

    # If it's a list/tuple of ALTs, take the first one
    if isinstance(alt_value, (list, tuple)):
        if len(alt_value) == 0:
            return None, None
        alt_value = alt_value[0]

    s = str(alt_value)

    # Main pattern: BreakEnd('chr5', 16982822, ...
    m = re.search(r"BreakEnd\('([^']+)',\s*([0-9]+)", s)
    if m:
        chrom2 = m.group(1)
        pos2 = int(m.group(2))
        return chrom2, pos2

    # Fallback: anything like chrXXX:12345 or chrXXX, 12345
    m = re.search(r"(chr[0-9XYM]+)[^\d]+([0-9]+)", s)
    if m:
        chrom2 = m.group(1)
        pos2 = int(m.group(2))
        return chrom2, pos2

    return None, None


def map_bnds_to_genes(strict_sv_set: pd.DataFrame,
                      genes: pd.DataFrame,
                      genes_only: pd.DataFrame, is_lncRNAs: bool,
                      window: int = 1_000_000) -> pd.DataFrame:
    """
    For BND SVs in strict_sv_set, add a 'gene_hits' column.

    Each BND has two breakpoints:
        side1: (chrom, pos)
        side2: (bnd_remote_chrom, bnd_remote_pos)

    We produce one gene_hit dict per (SV, gene, side), e.g.:

      {
        "gene_name": ...,
        "gene_id": ...,
        "strand": ...,
        "signed_dist": ...,
        "overlap_bp": ...,
        "region_hit": "promoter" / "exon_or_UTR" / "intron" /
                      "upstream_5kb" / "downstream_5kb" / "intergenic",
        "promoter_flag": 0/1,
        "gene_body_flag": 0/1,
        "exon_flag": 0/1,
        "intron_only_flag": 0/1,
        "upstream_5kb_flag": 0/1,
        "downstream_5kb_flag": 0/1,

        "hit_side": "bp1" or "bp2",
        "bp_index": 1 or 2,
        "bp_chrom": ...,
        "bp_pos": ...,
        "mate_chrom": ...,
        "mate_pos": ...,
      }

    Later you can explode and aggregate per gene.
    """
    df = strict_sv_set.copy()
    bnds = df[df["SVTYPE"] == "BND"].copy()

    gene_hits_col = []

    for idx, sv in bnds.iterrows():
        chrom1 = sv["chrom"]
        pos1   = int(sv["pos"])

        # If you've already parsed remote coords into columns:
        chrom2 = sv.get("bnd_remote_chrom", None)
        pos2   = sv.get("bnd_remote_pos", None)

        if pd.isna(chrom2) or chrom2 is None or pd.isna(pos2) or pos2 is None:
            chrom2 = None
            pos2   = None
        else:
            chrom2 = str(chrom2)
            pos2   = int(pos2)

        # We’ll collect hits from both sides
        gene_hits = []

        # Small helper to avoid code duplication
        def process_side(bk_chrom, bk_pos, side_label, bp_index,
                         mate_chrom, mate_pos):
            nonlocal gene_hits

            if bk_chrom is None:
                return

            # Genes on same chromosome within coarse window
            g_side = genes_only[genes_only["chrom"] == bk_chrom].copy()
            g_side = g_side[
                (g_side["start"] <= bk_pos + window) &
                (g_side["end"]   >= bk_pos - window)
            ]

            if g_side.empty:
                return

            # Features for these genes
            genes_features_side = genes[genes["gene_name"].isin(g_side["gene_name"])]

            for _, g in g_side.iterrows():
                gene_start = int(g["start"])
                gene_end   = int(g["end"])
                strand     = g.get("strand", "+")
                gene_name  = g["gene_name"]
                gene_id    = g.get("gene_id", None)

                signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                    bk_pos, bk_pos, gene_start, gene_end, strand
                )

                if abs(signed_dist) > window:
                    continue

                # gene-specific features
                gene_feats = genes_features_side[
                    genes_features_side["gene_name"] == gene_name
                ]

                # Treat breakpoint as a 1bp "span"
                region_info = classify_span_sv_gene_hit(
                    del_start=bk_pos,
                    del_end=bk_pos,
                    gene_start=gene_start,
                    gene_end=gene_end,
                    strand=strand,
                    signed_dist=signed_dist,
                    gene_features=gene_feats,
                )

                # Override hit_side for clarity
                region_info["hit_side"] = side_label
                region_info["bp_index"] = bp_index
                region_info["bp_chrom"] = bk_chrom
                region_info["bp_pos"]   = bk_pos
                region_info["mate_chrom"] = mate_chrom
                region_info["mate_pos"]   = mate_pos

                gene_hit = {
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "strand": strand,
                    "signed_dist": signed_dist,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    **region_info,
                }

                gene_hits.append(gene_hit)

        # Process side 1
        process_side(
            bk_chrom=chrom1,
            bk_pos=pos1,
            side_label="bp1",
            bp_index=1,
            mate_chrom=chrom2,
            mate_pos=pos2,
        )

        # Process side 2
        process_side(
            bk_chrom=chrom2,
            bk_pos=pos2,
            side_label="bp2",
            bp_index=2,
            mate_chrom=chrom1,
            mate_pos=pos1,
        )

        gene_hits_col.append(gene_hits)
    if not is_lncRNAs:
        bnds["gene_hits"] = gene_hits_col
        df.loc[bnds.index, "gene_hits"] = bnds["gene_hits"]
    else:
        bnds["lncRNA_hits"] = gene_hits_col
        df.loc[bnds.index, "lncRNA_hits"] = bnds["lncRNA_hits"]

    return df


def map_bnds_to_elements(strict_sv_set: pd.DataFrame,
                         elements: pd.DataFrame,
                         window: int = 1_000_000) -> pd.DataFrame:
    """
    For BND SVs in strict_sv_set, add an 'elem_hits' column.

    Each BND has two breakpoints:
        side1: (chrom, pos)
        side2: (bnd_remote_chrom, bnd_remote_pos)

    We produce one elem_hit dict per (SV, element, side), e.g.:

      {
        "elem_id": ...,
        "elem_type": ...,
        "chrom": ...,
        "elem_start": ...,
        "elem_end": ...,
        "signed_dist": ...,
        "overlap_bp": ...,
        "region_hit": "overlaps" /
                      "proximal_upstream" /
                      "proximal_downstream" /
                      "distal",
        "overlaps_flag": 0/1,
        "proximal_flag": 0/1,
        "distal_flag": 0/1,

        "hit_side": "bp1" or "bp2",
        "bp_index": 1 or 2,
        "bp_chrom": ...,
        "bp_pos": ...,
        "mate_chrom": ...,
        "mate_pos": ...,
      }

    Later you can explode and aggregate per element.
    """
    df = strict_sv_set.copy()
    bnds = df[df["SVTYPE"] == "BND"].copy()

    elem_hits_col = []

    for idx, sv in bnds.iterrows():
        chrom1 = sv["chrom"]
        pos1   = int(sv["pos"])

        # If you've already parsed remote coords into columns:
        chrom2 = sv.get("bnd_remote_chrom", None)
        pos2   = sv.get("bnd_remote_pos", None)

        if pd.isna(chrom2) or chrom2 is None or pd.isna(pos2) or pos2 is None:
            chrom2 = None
            pos2   = None
        else:
            chrom2 = str(chrom2)
            pos2   = int(pos2)

        # Collect hits from both sides
        elem_hits = []

        def process_side(bk_chrom, bk_pos, side_label, bp_index,
                         mate_chrom, mate_pos):
            nonlocal elem_hits

            if bk_chrom is None:
                return

            # Elements on same chromosome within coarse window
            e_side = elements[elements["chrom"] == bk_chrom].copy()
            e_side = e_side[
                (e_side["start"] <= bk_pos + window) &
                (e_side["end"]   >= bk_pos - window)
            ]

            if e_side.empty:
                return

            for _, e in e_side.iterrows():
                elem_start = int(e["start"])
                elem_end   = int(e["end"])
                elem_id    = e["elem_id"]
                elem_type  = e.get("type", None)

                # signed distance: breakpoint (point) vs element interval
                # strand is irrelevant for elements → pass "+"
                signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                    bk_pos, bk_pos, elem_start, elem_end, strand="+"
                )

                if abs(signed_dist) > window:
                    continue

                # Treat breakpoint as 1bp span
                region_info = classify_span_sv_elem_hit(
                    sv_start=bk_pos,
                    sv_end=bk_pos,
                    elem_start=elem_start,
                    elem_end=elem_end,
                    signed_dist=signed_dist,
                )

                # Override hit_side / bp metadata
                region_info["hit_side"] = side_label
                region_info["bp_index"] = bp_index
                region_info["bp_chrom"] = bk_chrom
                region_info["bp_pos"]   = bk_pos
                region_info["mate_chrom"] = mate_chrom
                region_info["mate_pos"]   = mate_pos

                elem_hit = {
                    "elem_id": elem_id,
                    "elem_type": elem_type,
                    "chrom": bk_chrom,
                    "elem_start": elem_start,
                    "elem_end": elem_end,
                    "signed_dist": signed_dist,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    **region_info,
                }

                elem_hits.append(elem_hit)

        # Side 1
        process_side(
            bk_chrom=chrom1,
            bk_pos=pos1,
            side_label="bp1",
            bp_index=1,
            mate_chrom=chrom2,
            mate_pos=pos2,
        )

        # Side 2
        process_side(
            bk_chrom=chrom2,
            bk_pos=pos2,
            side_label="bp2",
            bp_index=2,
            mate_chrom=chrom1,
            mate_pos=pos1,
        )

        elem_hits_col.append(elem_hits)

    bnds["elem_hits"] = elem_hits_col
    df.loc[bnds.index, "elem_hits"] = bnds["elem_hits"]

    return df


def map_SVs_to_genes(strict_sv_set: pd.DataFrame, genes,
                     genes_only: pd.DataFrame,
                     window: int = 1_000_000) -> pd.DataFrame:
    strict_sv_set = map_dels_to_genes(strict_sv_set, genes, genes_only, is_lncRNAs=False, window=1_000_000)
    strict_sv_set = map_dups_to_genes(strict_sv_set, genes, genes_only, is_lncRNAs=False, window=1_000_000)
    strict_sv_set = map_ins_to_genes(strict_sv_set, genes, genes_only, is_lncRNAs=False, window=1_000_000)


    bnd_mask = strict_sv_set["SVTYPE"] == "BND"
    coords = strict_sv_set.loc[bnd_mask, "alt"].apply(get_bnd_remote_coords)

    coords_df = pd.DataFrame(
        coords.tolist(),          # list of (chrom2, pos2)
        index=coords.index,       # same index as BND rows
        columns=["bnd_remote_chrom", "bnd_remote_pos"]
    )

    # Make sure columns exist (will be object/int as needed)
    if "bnd_remote_chrom" not in strict_sv_set.columns:
        strict_sv_set["bnd_remote_chrom"] = np.nan
    if "bnd_remote_pos" not in strict_sv_set.columns:
        strict_sv_set["bnd_remote_pos"] = np.nan

    # Assign by matching index explicitly
    strict_sv_set.loc[coords_df.index, ["bnd_remote_chrom", "bnd_remote_pos"]] = coords_df

    strict_sv_set = map_bnds_to_genes(strict_sv_set, genes, genes_only, is_lncRNAs=False, window=1_000_000)
    return strict_sv_set

def map_SVs_to_elements(strict_sv_set: pd.DataFrame,
                        elements: pd.DataFrame,
                        window: int = 1_000_000) -> pd.DataFrame:
    """
    Attach elem_hits for DEL, DUP, INS, and BND SVs.

    Assumes the following helpers exist:
      - map_dels_to_elements
      - map_dups_to_elements
      - map_ins_to_elements
      - map_bnds_to_elements
      - get_bnd_remote_coords (parses remote coords from 'alt' for BNDs)
    """
    # Span-type + point-type handled first
    strict_sv_set = map_dels_to_elements(strict_sv_set, elements, window=window)
    strict_sv_set = map_dups_to_elements(strict_sv_set, elements, window=window)
    strict_sv_set = map_ins_to_elements(strict_sv_set, elements, window=window)
 
    # Map BND breakpoints to elements
    strict_sv_set = map_bnds_to_elements(strict_sv_set, elements, window=window)

    return strict_sv_set


def map_SVs_to_lncRNAs(strict_sv_set: pd.DataFrame, lncRNAs: pd.DataFrame, lncRNAs_only: pd.DataFrame, window: int = 1_000_000) -> pd.DataFrame:
    strict_sv_set = map_dels_to_genes(strict_sv_set, lncRNAs, lncRNAs_only, is_lncRNAs=True, window=1_000_000)
    strict_sv_set = map_dups_to_genes(strict_sv_set, lncRNAs, lncRNAs_only, is_lncRNAs=True, window=1_000_000)
    strict_sv_set = map_ins_to_genes(strict_sv_set, lncRNAs, lncRNAs_only, is_lncRNAs=True, window=1_000_000)
    strict_sv_set = map_bnds_to_genes(strict_sv_set, lncRNAs, lncRNAs_only, is_lncRNAs=True, window=1_000_000)
    return strict_sv_set


import re

def breakend_to_vcf_alt(be_str: str):
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


def parse_breakend_string(s: str):
    """
    Parses BreakEnd('chrX', 101138003, '-', '-', 'T', True)
    Returns (chrom2, pos2, orientation_self, orientation_remote)
    """
    if not isinstance(s, str):
        return None

    # REGEX: capture chrom, pos, orientation_self, orientation_remote
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



def add_vep_hits_columns(df: pd.DataFrame, csq_description: str, primary_genes: list[str], lncRNAs_names) -> pd.DataFrame:
    """
    From a df that contains a 'CSQ' column and one row per (variant, ALT),
    add the following columns:

        - gene_hits       : list[dict]  (for Feature_type == "Transcript")
        - regulatory_hits : list[dict]  (Feature_type == "RegulatoryFeature")
        - motif_hits      : list[dict]  (Feature_type == "MotifFeature")

    Summary columns:

        - gene_symbols               : list of unique gene SYMBOLs hit by this variant/ALT
        - hits_canonical             : bool, hits any transcript with CANONICAL == "1"

        Global impact flags (any transcript):
        - has_missense
        - has_nonsense
        - has_frameshift
        - has_splice_effect

        Canonical-specific impact flags (CANONICAL == "1"):
        - has_missense_canonical
        - has_nonsense_canonical
        - has_frameshift_canonical
        - has_splice_effect_canonical

        MANE_SELECT-specific impact flags (MANE_SELECT not None):
        - has_missense_mane
        - has_nonsense_mane
        - has_frameshift_mane
        - has_splice_effect_mane
    """

    # 1) Parse CSQ column names from description
    if "Format:" in csq_description:
        fmt = csq_description.split("Format:")[1].strip().strip('"').strip()
    else:
        fmt = csq_description.strip().strip('"').strip()

    csq_cols = fmt.split("|")
    csq_index = {name: i for i, name in enumerate(csq_cols)}

    def get_field(parts, name, default=None):
        """
        Safely pull a field by VEP name from the parts list.
        Handles missing columns or truncated entries gracefully.
        """
        idx = csq_index.get(name)
        if idx is None:
            return default
        if idx >= len(parts):
            return default
        val = parts[idx]
        if val == "":
            return default
        return val

    # Fields for inner dicts

    GENE_FIELDS = [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",      # transcript ID
        "BIOTYPE",
        "EXON",
        "INTRON",
        "HGVSc",
        "HGVSp",
        "cDNA_position",
        "CDS_position",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "Existing_variation",
        "VARIANT_CLASS",
        "CANONICAL",
        "MANE",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "TSL",
        "APPRIS",
        "CCDS",
        "ENSP",
        "SWISSPROT",
        "TREMBL",
        "UNIPARC",
        "UNIPROT_ISOFORM",
        "GENE_PHENO",
        "SIFT",
        "PolyPhen",
        "DOMAINS",
        "AF",
        "gnomADe_AF",
        "gnomADg_AF",
        "MAX_AF",
        "CLIN_SIG",
        "SOMATIC",
        "PHENO",
        "PUBMED",
    ]

    REGULATORY_FIELDS = [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",          # ENSR...
        "BIOTYPE",
        "Existing_variation",
        "DISTANCE",
        "VARIANT_CLASS",
        "AF",
        "gnomADe_AF",
        "gnomADg_AF",
        "MAX_AF",
        "CLIN_SIG",
        "SOMATIC",
        "PHENO",
        "PUBMED",
    ]

    MOTIF_FIELDS = [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",
        "Existing_variation",
        "VARIANT_CLASS",
        "MOTIF_NAME",
        "MOTIF_POS",
        "HIGH_INF_POS",
        "MOTIF_SCORE_CHANGE",
        "TRANSCRIPTION_FACTORS",
        "AF",
        "gnomADe_AF",
        "gnomADg_AF",
        "MAX_AF",
    ]

    def has_consequence(gene_hits, target_terms:
                        set[str]) -> bool:
        """
        Check if ANY transcript (anywhere in gene_hits) has a consequence in target_terms.
        VEP uses '&' to separate multiple SO terms in one 'Consequence' field.
        """
        for hit in gene_hits:
            cons = hit.get("Consequence") or ""
            for term in cons.split("&"):
                if term.strip() in target_terms:
                    return True
        return False

    def has_consequence_filtered(
        gene_hits,
        target_terms: set[str],
        filter_field: str,
        filter_predicate
    ) -> bool:
        """
        Check if ANY transcript that satisfies filter_predicate(hit[filter_field])
        has a consequence in target_terms.
        e.g. filter_field = "CANONICAL", filter_predicate = lambda v: v == "1"
             filter_field = "MANE_SELECT", filter_predicate = lambda v: v is not None
        """
        for hit in gene_hits:
            if not filter_predicate(hit.get(filter_field)):
                continue
            cons = hit.get("Consequence") or ""
            for term in cons.split("&"):
                if term.strip() in target_terms:
                    return True
        return False

    def parse_csq_for_row(row: pd.Series, primary_genes, lncRNAs_names) -> pd.Series:
        """
        For a single row (variant+ALT), build:
          - gene_hits, regulatory_hits, motif_hits
          - gene_symbols, hits_canonical
          - global and canonical/MANE impact flags
        """
        alt_allele = str(row["alt"])

        val = row.get("CSQ")
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return pd.Series(
                {
                    "gene_hits_vep": [],
                    "regulatory_hits_vep": [],
                    "motif_hits": [],
                    "gene_symbols": [],
                    "hits_canonical": False,
                    "has_missense": False,
                    "has_nonsense": False,
                    "has_frameshift": False,
                    "has_splice_effect": False,
                    "has_missense_canonical": False,
                    "has_nonsense_canonical": False,
                    "has_frameshift_canonical": False,
                    "has_splice_effect_canonical": False,
                    "has_missense_mane": False,
                    "has_nonsense_mane": False,
                    "has_frameshift_mane": False,
                    "has_splice_effect_mane": False,
                }
            )

        # CSQ might be a list or a comma-separated string
        if isinstance(val, list):
            entries = val
        else:
            entries = str(val).split(",")


        gene_hits = []
        regulatory_hits = []
        motif_hits = []

        for entry in entries:
            entry = entry.strip()
            if not entry:
                continue

            parts = entry.split("|")

            allele = get_field(parts, "Allele")
            # if allele is None:
            #     continue

            # # keep only lines corresponding to this ALT
            # if allele != alt_allele:
            #     continue
            
            feature_type = get_field(parts, "Feature_type")
            if feature_type is None:
                continue

            if feature_type == "Transcript":
                hit = {}
                biotype = get_field(parts, "BIOTYPE")
                gene__name = get_field(parts, "SYMBOL")
                if gene__name not in primary_genes and gene__name not in lncRNAs_names:
                    continue
                for fname in GENE_FIELDS:
                    hit[fname] = get_field(parts, fname)
                gene_hits.append(hit)

            elif feature_type == "RegulatoryFeature":
                hit = {}
                for fname in REGULATORY_FIELDS:
                    hit[fname] = get_field(parts, fname)
                regulatory_hits.append(hit)

            elif feature_type == "MotifFeature":
                hit = {}
                for fname in MOTIF_FIELDS:
                    hit[fname] = get_field(parts, fname)
                motif_hits.append(hit)

        # outer gene list: unique SYMBOLs from gene_hits
        gene_symbol_list  = sorted(
            {
                h.get("SYMBOL")
                for h in gene_hits
                if h.get("SYMBOL") is not None
            }
        )
        gene_symbols = ",".join(gene_symbol_list) if gene_symbol_list else ""


        # canonical flag: any transcript with CANONICAL == "1"
        hits_canonical = any(
            (h.get("CANONICAL") == "1")
            for h in gene_hits
        )

        # --- global consequence flags (any transcript) ---
        MISSENSE = {"missense_variant"}
        NONSENSE = {"stop_gained"}
        FRAMESHIFT = {"frameshift_variant"}
        SPLICE = {
            "splice_acceptor_variant",
            "splice_donor_variant",
            "splice_region_variant",
        }

        has_missense = has_consequence(gene_hits, MISSENSE)
        has_nonsense = has_consequence(gene_hits, NONSENSE)
        has_frameshift = has_consequence(gene_hits, FRAMESHIFT)
        has_splice = has_consequence(gene_hits, SPLICE)

        # --- canonical-specific flags (CANONICAL == "1") ---
        canonical_pred = lambda v: v == "1"

        has_missense_canonical = has_consequence_filtered(
            gene_hits, MISSENSE, "CANONICAL", canonical_pred
        )
        has_nonsense_canonical = has_consequence_filtered(
            gene_hits, NONSENSE, "CANONICAL", canonical_pred
        )
        has_frameshift_canonical = has_consequence_filtered(
            gene_hits, FRAMESHIFT, "CANONICAL", canonical_pred
        )
        has_splice_canonical = has_consequence_filtered(
            gene_hits, SPLICE, "CANONICAL", canonical_pred
        )

        # --- MANE_SELECT-specific flags (MANE_SELECT not None) ---
        mane_pred = lambda v: v is not None

        has_missense_mane = has_consequence_filtered(
            gene_hits, MISSENSE, "MANE_SELECT", mane_pred
        )
        has_nonsense_mane = has_consequence_filtered(
            gene_hits, NONSENSE, "MANE_SELECT", mane_pred
        )
        has_frameshift_mane = has_consequence_filtered(
            gene_hits, FRAMESHIFT, "MANE_SELECT", mane_pred
        )
        has_splice_mane = has_consequence_filtered(
            gene_hits, SPLICE, "MANE_SELECT", mane_pred
        )


        return pd.Series(
            {
                "gene_hits_vep": gene_hits,
                "regulatory_hits_vep": regulatory_hits,
                "motif_hits": motif_hits,
                "gene_symbols": gene_symbols,
                "hits_canonical": hits_canonical,
                "has_missense": has_missense,
                "has_nonsense": has_nonsense,
                "has_frameshift": has_frameshift,
                "has_splice_effect": has_splice,
                "has_missense_canonical": has_missense_canonical,
                "has_nonsense_canonical": has_nonsense_canonical,
                "has_frameshift_canonical": has_frameshift_canonical,
                "has_splice_effect_canonical": has_splice_canonical,
                "has_missense_mane": has_missense_mane,
                "has_nonsense_mane": has_nonsense_mane,
                "has_frameshift_mane": has_frameshift_mane,
                "has_splice_effect_mane": has_splice_mane,
            }
        )

    hits_df = df.apply(parse_csq_for_row, axis=1, primary_genes=primary_genes, lncRNAs_names=lncRNAs_names)
    df = pd.concat([df, hits_df], axis=1)
    return df








def match_sv_to_genes_from_vcfs(vcf_dir: str, genes_path: str, elements_path: str, lncRNAs_path: str, output_dir: str, sample_annotations_path: str) -> pd.DataFrame:
    annotations = pd.read_csv(sample_annotations_path, sep="\t")
    genes = pd.read_csv(genes_path)
    elements = pd.read_csv(elements_path)
    lncRNAs = pd.read_csv(lncRNAs_path)

    # genes = genes[genes["feature"] == "gene"]
    genes = genes[genes["frame"].isin(['0',"."])]
    genes_only = genes[genes["feature"] == "gene"]
    primary_genes = genes_only["gene_name"].tolist()
    lncRNAs = lncRNAs[lncRNAs["frame"].isin(['0',"."])]
    lncRNAs_only = lncRNAs[lncRNAs["feature"] == "gene"]
    lncRNAs_names = lncRNAs_only["gene_name"].tolist()

    for i,vcf_file in enumerate(os.listdir(vcf_dir)):
        if not vcf_file.endswith(".vcf"):
            continue
      
        vcf_path = os.path.join(vcf_dir, vcf_file)
        df_sv, normal_sample, tumor_sample = load_manta_sv_vcf(vcf_path)
        strict_sv_set = get_strict_sv_set(df_sv)
        lenient_sv_set = get_lenient_sv_set(df_sv)
        strict_sv_set = map_SVs_to_genes(strict_sv_set, genes, genes_only, window=1_000_000)
        strict_sv_set = map_SVs_to_elements(strict_sv_set, elements, window=1_000_000)
        strict_sv_set = map_SVs_to_lncRNAs(strict_sv_set, lncRNAs, lncRNAs_only, window=1_000_000)
        strict_sv_set.rename(columns={"alt":"orig_alt"}, inplace=True)
         
        strict_sv_set["alt"] = strict_sv_set["orig_alt"].apply(breakend_to_vcf_alt)
        strict_sv_set["alt"].fillna(strict_sv_set["orig_alt"], inplace=True)
        strict_sv_set["alt"].apply(lambda x: "<" + x + ">" if x == "DEL" else x)
        # lenient_sv_set = map_SVs_to_genes(lenient_sv_set, genes, genes_only, window=1_000_000)
        # lenient_sv_set = map_SVs_to_elements(lenient_sv_set, elements, window=1_000_000)
        # lenient_sv_set = map_SVs_to_lncRNAs(lenient_sv_set, lncRNAs, lncRNAs_only, window=1_000_000)
        # lenient_sv_set.rename(columns={"alt":"orig_alt"}, inplace=True)
        # lenient_sv_set["alt"] = lenient_sv_set["orig_alt"].apply(breakend_to_vcf_alt)
        # lenient_sv_set["alt"].fillna(lenient_sv_set["orig_alt"], inplace=True)
        # lenient_sv_set["alt"].apply(lambda x: "<" + x + ">" if x == "DEL" else x)

        if annotations[annotations["File Name"] == vcf_file + ".gz"]["Tumor Descriptor"].iloc[0].startswith("Not"):
            sample_id = annotations[annotations["File Name"] == vcf_file + ".gz"]["Sample ID"].iloc[0].split(",")[1].strip()
        else:
            sample_id = annotations[annotations["File Name"] == vcf_file + ".gz"]["Sample ID"].iloc[0].split(",")[0].strip()
        print(sample_id)
        # sample_id = "attempt2"

        csq_descr = (
        'Consequence annotations from Ensembl VEP. Format: '
        'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|ARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_MID_AF|gnomADe_NFE_AF|gnomADe_REMAINING_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_REMAINING_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS'
        )
        strict_sv_set = add_vep_hits_columns(strict_sv_set, csq_descr, primary_genes, lncRNAs_names )
        # strict_sv_set.drop_duplicates(subset=["chrom","pos","ref","alt"], inplace=True)


        os.makedirs(output_dir, exist_ok=True)
        strict_sv_set.to_csv(os.path.join(output_dir, f"{sample_id}_strict_sv_set.csv"), index=False)
        # lenient_sv_set.to_csv(os.path.join(output_dir, f"{sample_id}_lenient_sv_set.csv"), index=False)

if __name__ == "__main__":
    VCF_DIR = sys.argv[1]
    genes_path = sys.argv[2]
    annotations_path = sys.argv[3]
    output_dir = sys.argv[4]
    elements_path = sys.argv[5]
    lncRNAs_path = sys.argv[6]
    match_sv_to_genes_from_vcfs(VCF_DIR, genes_path, elements_path, lncRNAs_path, output_dir, annotations_path)

