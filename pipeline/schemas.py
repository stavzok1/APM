"""
Schema definitions and builders for nested data structures.

These helpers ensure consistent structure across all evidence types,
even when data is missing for certain biosamples or assays.

ATAC: the saved ``atac_peaks_annotated.parquet`` row schema is unchanged (full nested
columns). After TAD streaming, the in-memory table is only
``pipeline.atac_peaks.peak_table.ATAC_SLIM_MEMORY_COLUMNS``; use
``pipeline.atac_peaks.load_atac_peaks_annotated`` when code needs ``gene_links`` / ``TAD_*``.
"""

from typing import Dict, Any, List, Optional




# =============================================================================
# SCREEN EVIDENCE SCHEMAS
# =============================================================================

def empty_assay_entry() -> Dict[str, Any]:
    """Empty entry for a single assay type within a biosample."""
    return {
        "score": None,
        "p_value": None,
        "strength": "none",
    }


def empty_biosample_assays(assay_types: List[str]) -> Dict[str, Dict[str, Any]]:
    """Empty dict of assays for a single biosample."""
    return {assay: empty_assay_entry() for assay in assay_types}


def empty_conservation_entry() -> Dict[str, int]:
    """Empty conservation counts for a single assay type."""
    return {
        "n_biosamples": 0,
        "n_strong": 0,
        "n_weak": 0,
    }


def empty_conservation_block(assay_types: List[str]) -> Dict[str, Dict[str, int]]:
    """Empty conservation dict across all assay types."""
    return {assay: empty_conservation_entry() for assay in assay_types}


def empty_screen_block(biosamples: List[str], assay_types: List[str]) -> Dict[str, Any]:
    """
    Empty SCREEN evidence block structure.
    """
    return {
        "per_biosample": {bio: empty_biosample_assays(assay_types) for bio in biosamples},
        "conservation_global": empty_conservation_block(assay_types),
        "conservation_breast": empty_conservation_block(assay_types),
    }


def ensure_screen_block(
    val: Any,
    biosamples: List[str],
    assay_types: List[str],
) -> Dict[str, Any]:
    """Coerce any value into a valid SCREEN block, filling missing keys."""
    if not isinstance(val, dict):
        return empty_screen_block(biosamples, assay_types)

    val.setdefault("per_biosample", {})
    val.setdefault("conservation_global", {})
    val.setdefault("conservation_breast", {})

    for bio in biosamples:
        if bio not in val["per_biosample"] or not isinstance(val["per_biosample"][bio], dict):
            val["per_biosample"][bio] = empty_biosample_assays(assay_types)
        else:
            for assay in assay_types:
                val["per_biosample"][bio].setdefault(assay, empty_assay_entry())

    for cons_key in ["conservation_global", "conservation_breast"]:
        if not isinstance(val[cons_key], dict):
            val[cons_key] = empty_conservation_block(assay_types)
        else:
            for assay in assay_types:
                val[cons_key].setdefault(assay, empty_conservation_entry())

    return val


# =============================================================================
# ABC EVIDENCE SCHEMAS
# =============================================================================

def empty_abc_celltype_entry() -> Dict[str, Any]:
    """Empty ABC entry for a single cell type."""
    return {
        "ABC_score": None,
        "ABC_num": None,
        "activity": None,
        "distance": None,
        "element_class": None,
        "is_self_promoter": False,
        "hic_pl_scaled": None,
        "powerlaw_score": None,
        "gene_expr": None,
        "promoter_activity_q": None,
        "gene_is_expressed": False,
        "rank_within_gene": None,
        "is_present": False,
        "is_strong": False,
    }


def empty_abc_block(celltypes: List[str]) -> Dict[str, Dict[str, Any]]:
    """Empty ABC block across all cell types."""
    return {ct: empty_abc_celltype_entry() for ct in celltypes}


def empty_abc_enhancer_entry(celltypes: List[str]) -> Dict[str, Any]:
    """Empty entry for a single ABC enhancer mapped to a cCRE."""
    return {
        "start": None,
        "end": None,
        "ABC_full": empty_abc_block(celltypes),
    }


# =============================================================================
# HICHIP EVIDENCE SCHEMAS
# =============================================================================

def empty_hichip_loop_entry() -> Dict[str, Any]:
    """Empty entry for a single HiChIP loop."""
    return {
        "loop_id": None,
        "counts": None,
        "score": None,
        "n_reps": 0,
        "partner_genes": [],
        "partner_cCREs": [],
        "anchor_a": {"chr": None, "start": None, "end": None},
        "anchor_b": {"chr": None, "start": None, "end": None},
    }


def empty_hichip_celltype_entry() -> Dict[str, Any]:
    """Empty HiChIP summary for a single cell type."""
    return {
        "n_loops": 0,
        "max_counts": None,
        "max_score": None,
        "loops": [],
    }


def empty_hichip_block(celltypes: List[str]) -> Dict[str, Dict[str, Any]]:
    """Empty HiChIP block across all cell types."""
    return {ct: empty_hichip_celltype_entry() for ct in celltypes}


# =============================================================================
# GENE LINKS SCHEMA (for cCREs)
# =============================================================================

def empty_gene_link_entry(
    screen_exp_biosamples: List[str],
    screen_exp_assays: List[str],
    screen_comp_biosamples: List[str],
    screen_comp_assays: List[str],
    abc_celltypes: Optional[List[str]] = None,
    hichip_celltypes: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Complete empty gene_link entry with all evidence types."""
    entry = {
        "gene_id": None,
        "gene_type": None,
        "screen_exp": empty_screen_block(screen_exp_biosamples, screen_exp_assays),
        "screen_comp": empty_screen_block(screen_comp_biosamples, screen_comp_assays),
        "ABC_enhancers": [],
        "hichip": {},
    }
    
    if hichip_celltypes:
        entry["hichip"] = empty_hichip_block(hichip_celltypes)
    
    return entry


# =============================================================================
# CELL LINE SIGNAL SCHEMAS
# =============================================================================

def empty_cell_line_signals(signal_names: List[str]) -> Dict[str, Any]:
    """Empty signal dict for a single cell line."""
    return {sig: None for sig in signal_names}


def ensure_cell_line_signals(
    val: Any,
    signal_names: List[str],
    in_key: str,
) -> Dict[str, Any]:
    """Ensure cell line signal dict has all expected keys."""
    if not isinstance(val, dict):
        result = empty_cell_line_signals(signal_names)
        result[in_key] = False
        return result
    
    val.setdefault(in_key, False)
    for sig in signal_names:
        val.setdefault(sig, None)
    
    return val


# =============================================================================
# MIRNA SCHEMA
# =============================================================================

def empty_mirna_entry() -> Dict[str, Any]:
    """Empty miRNA target entry."""
    return {
        "miRNA": None,
        "weighted_context_score": None,
        "num_sites": 0,
        "UTR_starts": [],
        "UTR_ends": [],
    }


# =============================================================================
# ATAC PEAK SCHEMAS
# =============================================================================

def empty_atac_body_overlap_entry() -> Dict[str, Any]:
    """Empty gene body overlap entry for ATAC peaks."""
    return {
        "overlaps": False,
        "overlap_bp": 0,
        "overlap_interval": None,
        "overlap_frac_of_peak": 0.0,
        "overlap_frac_of_gene": 0.0,
        "overlap_type": None,  # "promoter", "gene_body", None
    }


def empty_atac_gene_link_entry() -> Dict[str, Any]:
    """
    Empty gene link entry for ATAC peaks.
    
    Structure for each gene linked to a peak:
    {
        "gene_id": str,
        "gene_type": str,
        "dist_to_tss": int,
        "tier": str,
        "tss_position": int,
        "strand": str,
        "body_overlap": {...}
    }
    """
    return {
        "gene_id": None,
        "gene_type": None,
        "dist_to_tss": None,
        "tier": None,
        "tss_position": None,
        "strand": None,
        "body_overlap": empty_atac_body_overlap_entry(),
    }


def empty_atac_ccre_link_entry() -> Dict[str, Any]:
    """
    Empty cCRE link entry for ATAC peaks.
    """
    return {
        "cCRE_id": None,
        "ENCODE_id": None,
        "raw_type": None,
        "distance": None,
        "overlap": {
            "overlaps": False,
            "overlap_bp": 0,
            "overlap_interval": None,
            "overlap_frac_of_peak": 0.0,
            "overlap_frac_of_ccre": 0.0,
        },
    }


def empty_tad_boundary_overlaps_biosample() -> Dict[str, Any]:
    """
    Empty TAD boundary overlaps for a single biosample.
    """
    return {
        "overlaps_boundary": False,
        "n_boundaries": 0,
        "boundaries": [],
    }



# =============================================================================
# SNV/VARIANT SCHEMAS
# =============================================================================

def empty_vep_gene_hit() -> Dict[str, Any]:
    """Empty VEP gene/transcript hit entry."""
    return {
        "Allele": None,
        "Consequence": None,
        "IMPACT": None,
        "SYMBOL": None,
        "Gene": None,
        "Feature_type": None,
        "Feature": None,
        "BIOTYPE": None,
        "EXON": None,
        "INTRON": None,
        "HGVSc": None,
        "HGVSp": None,
        "cDNA_position": None,
        "CDS_position": None,
        "Protein_position": None,
        "Amino_acids": None,
        "Codons": None,
        "CANONICAL": None,
        "MANE_SELECT": None,
        "SIFT": None,
        "PolyPhen": None,
        "gnomADe_AF": None,
        "gnomADg_AF": None,
        "MAX_AF": None,
        "CLIN_SIG": None,
    }


def empty_vep_regulatory_hit() -> Dict[str, Any]:
    """Empty VEP regulatory feature hit entry."""
    return {
        "Allele": None,
        "Consequence": None,
        "IMPACT": None,
        "Feature_type": None,
        "Feature": None,
        "BIOTYPE": None,
        "DISTANCE": None,
        "VARIANT_CLASS": None,
    }


def empty_vep_motif_hit() -> Dict[str, Any]:
    """Empty VEP motif feature hit entry."""
    return {
        "Allele": None,
        "Consequence": None,
        "Feature_type": None,
        "Feature": None,
        "MOTIF_NAME": None,
        "MOTIF_POS": None,
        "HIGH_INF_POS": None,
        "MOTIF_SCORE_CHANGE": None,
        "TRANSCRIPTION_FACTORS": None,
    }


def empty_ccre_hit() -> Dict[str, Any]:
    """Empty cCRE overlap hit entry."""
    return {
        "cCRE_id": None,
        "elem_type": None,
        "chrom": None,
        "elem_start": None,
        "elem_end": None,
        "genes_by_exact_dist": None,
    }


def empty_snv_entry() -> Dict[str, Any]:
    """
    Empty SNV/variant entry with all expected fields.
    
    Structure:
    {
        # Core variant info
        "chrom": str,
        "pos": int,
        "ref": str,
        "alt": str,
        "qual": float,
        "filter": str,
        
        # VAF info
        "tumor_vaf": float,
        "normal_vaf": float,
        
        # VEP parsed hits
        "gene_hits": List[Dict],
        "regulatory_hits": List[Dict],
        "motif_hits": List[Dict],
        
        # Summary flags
        "gene_symbols": str,
        "hits_canonical": bool,
        "has_missense": bool,
        "has_nonsense": bool,
        "has_frameshift": bool,
        "has_splice_effect": bool,
        
        # cCRE overlap
        "cCRE_hits": List[Dict],

        # Optional MEME FIMO on reference windows (``SNV.snv_fimo``)
        "fimo_hits": List[Dict],
        # ChIP unified overlap (``SNV.snv_chip``)
        "snv_chip_hits": List[Dict],
        "snv_chip_aggregate": Dict[str, Any],
    }
    """
    return {
        "chrom": None,
        "pos": None,
        "ref": None,
        "alt": None,
        "qual": None,
        "filter": None,
        "tumor_vaf": None,
        "normal_vaf": None,
        "gene_hits": [],
        "regulatory_hits": [],
        "motif_hits": [],
        "gene_symbols": "",
        "hits_canonical": False,
        "has_missense": False,
        "has_nonsense": False,
        "has_frameshift": False,
        "has_splice_effect": False,
        "has_missense_canonical": False,
        "has_nonsense_canonical": False,
        "has_frameshift_canonical": False,
        "has_splice_effect_canonical": False,
        "cCRE_hits": [],
        "fimo_hits": [],
        "snv_chip_hits": [],
        "snv_chip_aggregate": {"by_tf_source_stratum": []},
    }


def empty_sample_snv_summary() -> Dict[str, Any]:
    """
    Empty per-sample SNV summary.
    
    Used for aggregating SNV burden per sample.
    """
    return {
        "sample_id": None,
        "n_total_variants": 0,
        "n_coding_variants": 0,
        "n_missense": 0,
        "n_nonsense": 0,
        "n_frameshift": 0,
        "n_splice": 0,
        "n_regulatory": 0,
        "n_ccre_hits": 0,
        "genes_affected": [],
        "ccre_types_affected": {},
    }

# =============================================================================
# SV (STRUCTURAL VARIANT) SCHEMAS
# =============================================================================

def empty_sv_gene_hit() -> Dict[str, Any]:
    """
    Empty structure for an SV-to-gene hit.
    
    Captures spatial relationship, overlap details, and region classification.
    """
    return {
        "gene_name": None,
        "gene_id": None,
        "strand": None,
        "signed_dist": None,
        "overlap_start": None,
        "overlap_end": None,
        "overlap_bp": 0,
        "overlap_percent": 0.0,
        "promoter_flag": 0,
        "gene_body_flag": 0,
        "exon_flag": 0,
        "intron_only_flag": 0,
        "upstream_5kb_flag": 0,
        "downstream_5kb_flag": 0,
        "region_hit": "",
        "hit_side": "span",  # "span", "point", "bp1", "bp2"
        "stop_codon_flag": 0,
        "start_codon_flag": 0,
        "transcript_id": None,
        "transcript_type": None,
    }


def empty_sv_element_hit() -> Dict[str, Any]:
    """
    Empty structure for an SV-to-regulatory-element hit.
    
    Captures spatial relationship and overlap classification.
    """
    return {
        "elem_id": None,
        "elem_type": None,
        "chrom": None,
        "elem_start": None,
        "elem_end": None,
        "signed_dist": None,
        "overlap_start": None,
        "overlap_end": None,
        "overlap_bp": 0,
        "overlap_percent": 0.0,
        "overlaps_flag": 0,
        "region_hit": "",  # "overlaps", "proximal_upstream", "proximal_downstream", "distal"
        "proximal_flag": 0,
        "distal_flag": 0,
        "hit_side": "span",  # "span", "point", "bp1", "bp2"
    }


def empty_sv_motif_hit() -> Dict[str, Any]:
    """
    Empty structure for a motif hit from FIMO scanning.
    
    Used for both flank motifs and regulatory element motifs.
    """
    return {
        "start": None,
        "end": None,
        "TF": None,
        "motif_id": None,
        "score": None,
        "p_value": None,
        "q_value": None,
        "strand": None,
        "distance_to_pos": None,
    }


def empty_sv_flank_motif_hit() -> Dict[str, Any]:
    """
    Empty structure for a flank-region motif hit.
    
    Extends motif hit with flank side information.
    """
    base = empty_sv_motif_hit()
    base["flank_side"] = None  # "left" or "right"
    return base


def empty_sv_vep_gene_hit() -> Dict[str, Any]:
    """
    Empty structure for VEP gene annotation.
    
    Contains transcript-level VEP consequence information.
    """
    return {
        "Allele": None,
        "Consequence": None,
        "IMPACT": None,
        "SYMBOL": None,
        "Gene": None,
        "Feature_type": None,
        "Feature": None,
        "BIOTYPE": None,
        "EXON": None,
        "INTRON": None,
        "HGVSc": None,
        "HGVSp": None,
        "CANONICAL": None,
        "MANE_SELECT": None,
        "SIFT": None,
        "PolyPhen": None,
    }


def empty_sv_vep_regulatory_hit() -> Dict[str, Any]:
    """
    Empty structure for VEP regulatory feature annotation.
    """
    return {
        "Allele": None,
        "Consequence": None,
        "IMPACT": None,
        "Feature_type": None,
        "Feature": None,  # ENSR ID
        "BIOTYPE": None,
        "DISTANCE": None,
    }


def empty_sv_vep_motif_hit() -> Dict[str, Any]:
    """
    Empty structure for VEP motif feature annotation.
    """
    return {
        "Allele": None,
        "Consequence": None,
        "IMPACT": None,
        "Feature_type": None,
        "Feature": None,
        "MOTIF_NAME": None,
        "MOTIF_POS": None,
        "HIGH_INF_POS": None,
        "MOTIF_SCORE_CHANGE": None,
        "TRANSCRIPTION_FACTORS": None,
    }


def empty_sv_record() -> Dict[str, Any]:
    """
    Complete empty structure for a processed SV record.
    
    This is the main output structure containing:
    - Basic SV info (position, type, quality)
    - Gene hits from spatial mapping
    - Element hits from regulatory element mapping
    - lncRNA hits
    - VEP annotations
    - Motif hits (both flank and element-overlapping)
    """
    return {
        # Basic info
        "id": None,
        "chrom": None,
        "pos": None,
        "END": None,
        "SVTYPE": None,
        "SVLEN": None,
        "ref": None,
        "alt": None,
        "qual": None,
        "filter": None,
        
        # Somatic scores
        "SOMATICSCORE": None,
        "normal_alt": 0,
        "tumor_alt": 0,
        "tumor_sr_alt": 0,
        "tumor_pr_alt": 0,
        
        # BND-specific
        "bnd_remote_chrom": None,
        "bnd_remote_pos": None,
        
        # Spatial mappings (list of dicts)
        "gene_hits": [],
        "elem_hits": [],
        "lncRNA_hits": [],
        
        # VEP annotations (list of dicts)
        "gene_hits_vep": [],
        "regulatory_hits_vep": [],
        "motif_hits": [],
        
        # VEP summary flags
        "gene_symbols": "",
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
        
        # FIMO motif hits
        "flank_motif_hits": [],
    }


def empty_sv_bed_interval() -> Dict[str, Any]:
    """
    Empty structure for an SV BED interval (for FASTA extraction).
    """
    return {
        "chrom": None,
        "start": None,
        "end": None,
        "name": None,  # Format: "{sv_id}|{SVTYPE}_flank_left" or "{sv_id}|{SVTYPE}|elem:{elem_id}|..."
    }


# =============================================================================
# RPPA SCHEMAS
# =============================================================================

def empty_rppa_sample_entry() -> Dict[str, Any]:
    """
    Empty entry for RPPA data per sample.
    
    Structure mirrors protein-level evidence in gene_links.
    """
    return {
        "sample_id": None,
        "expression": {},  # target → value
        "activation_ratios": {},  # gene → ratio
        "panel_scores": {},  # panel_name → score
        "signaling_blocks": {},  # block_type → bool
        "visibility_state": None,
    }


def empty_rppa_target_entry() -> Dict[str, Any]:
    """
    Empty entry for a single RPPA target.
    """
    return {
        "target": None,
        "gene_name": None,
        "expression": None,
        "z_score": None,
        "percentile": None,
        "validation_status": None,
        "is_phospho": False,
        "phospho_site": None,
    }


def empty_rppa_gene_evidence() -> Dict[str, Any]:
    """
    Empty RPPA evidence block for a gene.
    
    Structure:
    {
        "total_protein": {target, expression, z_score},
        "phospho_forms": [{target, site, expression, z_score, activation_ratio}, ...],
        "validation_status": "Valid" | "Caution",
    }
    """
    return {
        "total_protein": {
            "target": None,
            "expression": None,
            "z_score": None,
        },
        "phospho_forms": [],
        "validation_status": None,
    }


def empty_rppa_panel_scores() -> Dict[str, Optional[float]]:
    """
    Empty panel scores dictionary.
    """
    return {
        "IFN_activated": None,
        "DDR_activation": None,
        "cGAS_STING": None,
        "DNA_repair": None,
        "checkpoint": None,
        "PI3K_AKT": None,
        "mTOR": None,
        "STAT3_suppressive": None,
        "lymphocyte_infiltration": None,
        "cytolytic": None,
        "proliferation": None,
        "apoptosis_balance": None,
    }


def empty_signaling_blocks() -> Dict[str, bool]:
    """
    Empty signaling blocks dictionary.
    """
    return {
        "ddr_sting_block": False,
        "sting_irf1_block": False,
        "stat3_override": False,
        "checkpoint_escape": False,
        "pi3k_activated": False,
        "ddr_ifn_chain_intact": False,
        "immune_desert": False,
    }


# =============================================================================
# VALIDATION HELPERS
# =============================================================================

def validate_gene_links_structure(gene_links: Dict[str, Dict]) -> List[str]:
    """Validate a gene_links dict and return list of issues."""
    issues = []
    
    if not isinstance(gene_links, dict):
        return ["gene_links is not a dict"]
    
    for gene_name, entry in gene_links.items():
        if not isinstance(entry, dict):
            issues.append(f"{gene_name}: entry is not a dict")
            continue
        
        for key in ["screen_exp", "screen_comp", "ABC_enhancers", "hichip"]:
            if key not in entry:
                issues.append(f"{gene_name}: missing key '{key}'")
        
        for screen_key in ["screen_exp", "screen_comp"]:
            if screen_key in entry and isinstance(entry[screen_key], dict):
                for sub_key in ["per_biosample", "conservation_global", "conservation_breast"]:
                    if sub_key not in entry[screen_key]:
                        issues.append(f"{gene_name}.{screen_key}: missing '{sub_key}'")
    
    return issues


def validate_atac_gene_links(gene_links: Dict[str, Dict]) -> List[str]:
    """Validate ATAC peak gene_links structure."""
    issues = []
    
    if not isinstance(gene_links, dict):
        return ["gene_links is not a dict"]
    
    required_keys = ["gene_id", "dist_to_tss", "tier", "body_overlap"]
    
    for gene_name, entry in gene_links.items():
        if not isinstance(entry, dict):
            issues.append(f"{gene_name}: entry is not a dict")
            continue
        
        for key in required_keys:
            if key not in entry:
                issues.append(f"{gene_name}: missing key '{key}'")
    
    return issues


def validate_atac_ccre_links(ccre_links: List[Dict]) -> List[str]:
    """Validate ATAC peak ccre_links structure."""
    issues = []
    
    if not isinstance(ccre_links, list):
        return ["ccre_links is not a list"]
    
    required_keys = ["cCRE_id", "ENCODE_id", "raw_type", "distance", "overlap"]
    
    for i, entry in enumerate(ccre_links):
        if not isinstance(entry, dict):
            issues.append(f"ccre_links[{i}]: entry is not a dict")
            continue
        
        for key in required_keys:
            if key not in entry:
                issues.append(f"ccre_links[{i}]: missing key '{key}'")
    
    return issues


def validate_snv_dataframe(df) -> List[str]:
    """
    Validate SNV DataFrame has expected columns.
    
    Returns list of issues (empty if valid).
    """
    issues = []
    
    required_cols = ["chrom", "pos", "ref", "alt"]
    for col in required_cols:
        if col not in df.columns:
            issues.append(f"Missing required column: {col}")
    
    recommended_cols = [
        "filter", "tumor_vaf", "normal_vaf",
        "gene_hits", "cCRE_hits",
    ]
    for col in recommended_cols:
        if col not in df.columns:
            issues.append(f"Missing recommended column: {col}")
    
    return issues



# =============================================================================
# SV VALIDATION HELPERS
# =============================================================================

def validate_sv_gene_hit(hit: Dict[str, Any]) -> List[str]:
    """
    Validate an SV gene hit dict and return list of issues.
    """
    issues = []
    required = ["gene_name", "signed_dist", "hit_side"]
    
    if not isinstance(hit, dict):
        return ["gene_hit is not a dict"]
    
    for key in required:
        if key not in hit:
            issues.append(f"missing required key '{key}'")
    
    if hit.get("hit_side") not in ["span", "point", "bp1", "bp2", None]:
        issues.append(f"invalid hit_side: {hit.get('hit_side')}")
    
    return issues


def validate_sv_element_hit(hit: Dict[str, Any]) -> List[str]:
    """
    Validate an SV element hit dict and return list of issues.
    """
    issues = []
    required = ["elem_id", "signed_dist", "hit_side"]
    
    if not isinstance(hit, dict):
        return ["elem_hit is not a dict"]
    
    for key in required:
        if key not in hit:
            issues.append(f"missing required key '{key}'")
    
    return issues


def validate_sv_record(record: Dict[str, Any]) -> List[str]:
    """
    Validate a complete SV record and return list of issues.
    """
    issues = []
    
    if not isinstance(record, dict):
        return ["sv_record is not a dict"]
    
    # Required fields
    required = ["id", "chrom", "pos", "SVTYPE"]
    for key in required:
        if key not in record:
            issues.append(f"missing required key '{key}'")
    
    # Validate nested gene_hits
    if "gene_hits" in record:
        if not isinstance(record["gene_hits"], list):
            issues.append("gene_hits should be a list")
        else:
            for i, hit in enumerate(record["gene_hits"]):
                hit_issues = validate_sv_gene_hit(hit)
                for issue in hit_issues:
                    issues.append(f"gene_hits[{i}]: {issue}")
    
    # Validate nested elem_hits
    if "elem_hits" in record:
        if not isinstance(record["elem_hits"], list):
            issues.append("elem_hits should be a list")
        else:
            for i, hit in enumerate(record["elem_hits"]):
                hit_issues = validate_sv_element_hit(hit)
                for issue in hit_issues:
                    issues.append(f"elem_hits[{i}]: {issue}")
    
    return issues


def ensure_sv_record_structure(record: Dict[str, Any]) -> Dict[str, Any]:
    """
    Ensure an SV record has all expected fields with defaults.
    """
    template = empty_sv_record()
    
    if not isinstance(record, dict):
        return template
    
    # Fill missing keys with defaults
    for key, default_value in template.items():
        if key not in record:
            record[key] = default_value
        elif isinstance(default_value, list) and not isinstance(record[key], list):
            record[key] = []
    
    return record


# =============================================================================
# RPPA VALIDATION HELPERS
# =============================================================================


def validate_rppa_expression_matrix(
    expr_matrix,  # pd.DataFrame
) -> List[str]:
    """
    Validate RPPA expression matrix structure.
    
    Returns:
        List of issues (empty if valid)
    """
    issues = []
    
    # Check type
    import pandas as pd
    if not isinstance(expr_matrix, pd.DataFrame):
        return ["Expression matrix is not a DataFrame"]
    
    # Check not empty
    if expr_matrix.empty:
        issues.append("Expression matrix is empty")
    
    # Check index (samples)
    if expr_matrix.index.duplicated().any():
        n_dup = expr_matrix.index.duplicated().sum()
        issues.append(f"Expression matrix has {n_dup} duplicate sample IDs")
    
    # Check columns (targets)
    if expr_matrix.columns.duplicated().any():
        n_dup = expr_matrix.columns.duplicated().sum()
        issues.append(f"Expression matrix has {n_dup} duplicate target names")
    
    # Check values are numeric
    non_numeric = []
    for col in expr_matrix.columns:
        if not pd.api.types.is_numeric_dtype(expr_matrix[col]):
            non_numeric.append(col)
    if non_numeric:
        issues.append(f"Non-numeric columns: {non_numeric[:5]}...")
    
    return issues


def validate_rppa_annotation(
    annotation,  # pd.DataFrame
) -> List[str]:
    """
    Validate RPPA annotation DataFrame.
    
    Returns:
        List of issues (empty if valid)
    """
    issues = []
    
    import pandas as pd
    if not isinstance(annotation, pd.DataFrame):
        return ["Annotation is not a DataFrame"]
    
    # Required columns
    required = ["AGID", "peptide_target", "gene_name", "validation_status"]
    missing = [c for c in required if c not in annotation.columns]
    if missing:
        issues.append(f"Missing required columns: {missing}")
    
    # Check for duplicates
    if "peptide_target" in annotation.columns:
        if annotation["peptide_target"].duplicated().any():
            n_dup = annotation["peptide_target"].duplicated().sum()
            issues.append(f"Annotation has {n_dup} duplicate peptide_targets")
    
    # Validate validation_status values
    if "validation_status" in annotation.columns:
        valid_statuses = {"Valid", "Caution", "Under Evaluation"}
        actual = set(annotation["validation_status"].str.strip().str.title().unique())
        unknown = actual - valid_statuses
        if unknown:
            issues.append(f"Unknown validation statuses: {unknown}")
    
    return issues


def validate_phospho_pairs(
    phospho_pairs: Dict[str, Dict],
) -> List[str]:
    """
    Validate phospho pairs structure.
    
    Returns:
        List of issues (empty if valid)
    """
    issues = []
    
    if not isinstance(phospho_pairs, dict):
        return ["phospho_pairs is not a dict"]
    
    for gene, pair in phospho_pairs.items():
        if not isinstance(pair, dict):
            issues.append(f"{gene}: pair entry is not a dict")
            continue
        
        if "total" not in pair:
            issues.append(f"{gene}: missing 'total' key")
        if "phospho" not in pair:
            issues.append(f"{gene}: missing 'phospho' key")
        elif not isinstance(pair["phospho"], list):
            issues.append(f"{gene}: 'phospho' is not a list")
    
    return issues


def build_rppa_gene_evidence(
    gene_name: str,
    expr_matrix,  # pd.DataFrame
    target_to_gene: Dict[str, str],
    phospho_pairs: Dict[str, Dict],
    sample_id: str,
) -> Dict[str, Any]:
    """
    Build RPPA evidence structure for a specific gene and sample.
    
    Args:
        gene_name: Gene of interest
        expr_matrix: Expression matrix (samples × targets)
        target_to_gene: Mapping from target to gene
        phospho_pairs: Phospho pair information
        sample_id: Sample to extract data for
    
    Returns:
        RPPA evidence dict for this gene-sample pair
    """
    evidence = empty_rppa_gene_evidence()
    
    if sample_id not in expr_matrix.index:
        return evidence
    
    # Find targets for this gene
    gene_targets = [t for t, g in target_to_gene.items() if g == gene_name]
    
    if not gene_targets:
        return evidence
    
    # Get phospho pair info
    pair_info = phospho_pairs.get(gene_name, {})
    total_target = pair_info.get("total")
    phospho_targets = pair_info.get("phospho", [])
    phospho_sites = pair_info.get("phospho_sites", {})
    
    # Total protein
    if total_target and total_target in expr_matrix.columns:
        val = expr_matrix.loc[sample_id, total_target]
        col_vals = expr_matrix[total_target]
        z = (val - col_vals.mean()) / col_vals.std() if col_vals.std() > 0 else 0
        
        evidence["total_protein"] = {
            "target": total_target,
            "expression": float(val) if not pd.isna(val) else None,
            "z_score": float(z) if not pd.isna(z) else None,
        }
    
    # Phospho forms
    for p_target in phospho_targets:
        if p_target not in expr_matrix.columns:
            continue
        
        val = expr_matrix.loc[sample_id, p_target]
        col_vals = expr_matrix[p_target]
        z = (val - col_vals.mean()) / col_vals.std() if col_vals.std() > 0 else 0
        
        # Activation ratio if total available
        activation_ratio = None
        if total_target and total_target in expr_matrix.columns:
            total_val = expr_matrix.loc[sample_id, total_target]
            if not pd.isna(val) and not pd.isna(total_val):
                activation_ratio = float(val - total_val)  # Log-space ratio
        
        evidence["phospho_forms"].append({
            "target": p_target,
            "site": phospho_sites.get(p_target, ""),
            "expression": float(val) if not pd.isna(val) else None,
            "z_score": float(z) if not pd.isna(z) else None,
            "activation_ratio": activation_ratio,
        })
    
    return evidence


def ensure_rppa_panel_scores(
    val: Any,
) -> Dict[str, Optional[float]]:
    """
    Ensure panel scores dict has all expected keys.
    """
    if not isinstance(val, dict):
        return empty_rppa_panel_scores()
    
    template = empty_rppa_panel_scores()
    for key in template:
        val.setdefault(key, None)
    
    return val


def ensure_signaling_blocks(
    val: Any,
) -> Dict[str, bool]:
    """
    Ensure signaling blocks dict has all expected keys.
    """
    if not isinstance(val, dict):
        return empty_signaling_blocks()
    
    template = empty_signaling_blocks()
    for key in template:
        val.setdefault(key, False)
    
    return val
