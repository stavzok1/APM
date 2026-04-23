"""CNV pipeline orchestrator."""
import os
from typing import Optional

import pandas as pd
from .loader import load_cnv_file, extract_sample_id_from_annotations
from .features import add_basic_cnv_features
from .gene_hits import (annotate_cnv_with_gene_hits,
    annotate_cnv_with_lncrna_hits, annotate_cnv_with_mirna_hits)
from .elem_hits import annotate_cnv_with_elem_hits
from .gene_summary import GeneCnvRules, summarize_gene_cnv
from .gene_level_ascat import (
    build_ascat_gene_sample_table,
    is_ascat_gene_level_filename,
    load_ascat_gene_level_tsv,
)
from ..sample_ids import add_tcga_id_columns_inplace
from ..genes.symbol_normalization import normalize_annotation_gene_names

def annotate_single_sample(raw, ga, la, mi, el, pu=2000, pd_=500, pw=5000):
    c = add_basic_cnv_features(raw)
    c = annotate_cnv_with_gene_hits(c, ga, pu, pd_, pw)
    c = annotate_cnv_with_lncrna_hits(c, la, pu, pd_, pw)
    c = annotate_cnv_with_mirna_hits(c, mi, pu, pd_, pw)
    c = annotate_cnv_with_elem_hits(c, el, pw)
    # Atlas uses ``lncRNA_hits`` (no trailing "s") for the nested CNV column;
    # keep the legacy name and add the canonical alias so downstream users can
    # rely on either.
    if "lncRNAs_hits" in c.columns and "lncRNA_hits" not in c.columns:
        c["lncRNA_hits"] = c["lncRNAs_hits"]
    return c

def process_cnv_directory(
    cnv_dir,
    genes_path,
    lncrnas_path,
    elements_path,
    ann_path,
    out_dir,
    *,
    mirna_path=None,
    sep="\t",
    pu=2000,
    pd_=500,
    pw=5000,
    suf=(".txt", ".tsv", ".seg"),
    write_gene_tables=True,
    gene_tables_subdir="gene_tables",
    gene_tables_root=None,
    gene_min_body_overlap_pct=None,
    gene_tables_combined_name="cnv_gene_calls_all_samples.csv",
    process_ascat_gene_level: bool = True,
    ascat_gene_table_stem: Optional[str] = None,
    ascat_gene_tables_combined_name: str = "cnv_gene_calls_all_samples_ascat3.csv",
    ascat_regulatory_window_bp: Optional[int] = None,
):
    an = pd.read_csv(ann_path, sep="\t")
    ga_full = normalize_annotation_gene_names(
        pd.read_csv(genes_path),
        ("gene_name",),
    )
    ga = ga_full
    if "frame" in ga.columns:
        ga = ga[ga["frame"].isin(["0", "."])]
    la = normalize_annotation_gene_names(
        pd.read_csv(lncrnas_path, low_memory=False),
        ("gene_name",),
    )
    if "frame" in la.columns:
        la = la[la["frame"].isin(["0", "."])]
    # Prefer mature-arm loci table if available for arm-resolution by coordinates.
    if mirna_path is None:
        try:
            from ..config import PATHS
            mirna_path = getattr(PATHS, "mirna_mature_loci_csv", None) or getattr(PATHS, "mirna_path", None)
        except Exception:
            mirna_path = mirna_path
    mi = pd.read_csv(mirna_path)
    from ..regulatory_elements import load_regulatory_element_focus
    el = load_regulatory_element_focus(elements_path)
    os.makedirs(out_dir, exist_ok=True)
    gene_tables = []
    ascat_gene_tables = []
    # Segment-level annotated CSVs → out_dir; gene-centric tables → gene_tables_root (or out_dir/gene_tables_subdir)
    gdir_base = gene_tables_root if gene_tables_root is not None else os.path.join(out_dir, gene_tables_subdir)
    try:
        from ..config import THRESHOLDS as _TH
    except Exception:
        _TH = None
    if gene_min_body_overlap_pct is None:
        gene_min_body_overlap_pct = (
            float(getattr(_TH, "cnv_gene_min_body_overlap_pct", 30.0)) if _TH is not None else 30.0
        )
    ascat_stem = ascat_gene_table_stem
    if ascat_stem is None:
        ascat_stem = (
            str(getattr(_TH, "cnv_ascat3_gene_table_stem", "cnv_gene_calls_ascat3"))
            if _TH is not None
            else "cnv_gene_calls_ascat3"
        )
    rwin = ascat_regulatory_window_bp
    if rwin is None:
        rwin = int(getattr(_TH, "cnv_ascat3_regulatory_window_bp", 250_000)) if _TH is not None else 250_000
    gene_rules = GeneCnvRules(min_body_overlap_pct=float(gene_min_body_overlap_pct))
    for fn in sorted(os.listdir(cnv_dir)):
        path = os.path.join(cnv_dir, fn)
        if is_ascat_gene_level_filename(fn):
            if not process_ascat_gene_level:
                print(f"Skip ASCAT gene-level (disabled): {fn}")
                continue
            print(f"Processing ASCAT gene-level: {fn}")
            raw_gl = load_ascat_gene_level_tsv(path, sep=sep)
            sid = extract_sample_id_from_annotations(an, fn)
            print(f"  -> {sid}")
            gt = build_ascat_gene_sample_table(
                raw_gl,
                sid,
                genes_lookup=ga_full,
                regulatory_window_bp=int(rwin),
            )
            add_tcga_id_columns_inplace(gt, raw_id=sid)
            ascat_gene_tables.append(gt)
            if write_gene_tables:
                gdir = str(gdir_base)
                os.makedirs(gdir, exist_ok=True)
                gt.to_csv(os.path.join(gdir, f"{sid}_{ascat_stem}.csv"), index=False)
            continue

        if suf and not fn.endswith(suf):
            continue
        print(f"Processing: {fn}")
        raw = load_cnv_file(path, sep=sep)
        cnv = annotate_single_sample(raw, ga, la, mi, el, pu, pd_, pw)
        sid = extract_sample_id_from_annotations(an, fn)
        print(f"  -> {sid}")
        if "GDC_Aliquot" in cnv.columns:
            cnv.drop(columns=["GDC_Aliquot"], inplace=True)
        add_tcga_id_columns_inplace(cnv, raw_id=sid)
        # Atlas: one ``sample_id`` column carrying the TCGA sample id.
        if "sample_id" not in cnv.columns:
            cnv["sample_id"] = cnv.get("sample", sid)
        cnv.to_csv(os.path.join(out_dir, f"{sid}_cnv_annotated.csv"), index=False)
        if write_gene_tables:
            gt = summarize_gene_cnv(cnv, rules=gene_rules)
            gt.insert(0, "sample_id", sid)
            gene_tables.append(gt)
            gdir = str(gdir_base)
            os.makedirs(gdir, exist_ok=True)
            gt.to_csv(os.path.join(gdir, f"{sid}_cnv_gene_calls.csv"), index=False)
    gdir = str(gdir_base)
    os.makedirs(gdir, exist_ok=True)
    if write_gene_tables and gene_tables_combined_name:
        all_gt = pd.concat(gene_tables, ignore_index=True) if gene_tables else pd.DataFrame()
        all_gt.to_csv(os.path.join(gdir, gene_tables_combined_name), index=False)
    if write_gene_tables and ascat_gene_tables_combined_name and ascat_gene_tables:
        all_a = pd.concat(ascat_gene_tables, ignore_index=True)
        all_a.to_csv(os.path.join(gdir, ascat_gene_tables_combined_name), index=False)
    print(f"Done -> {out_dir}")
