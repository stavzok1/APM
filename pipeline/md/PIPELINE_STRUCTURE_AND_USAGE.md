# Pipeline Structure & Usage Guide

## Directory Structure

```
pipeline/
├── __init__.py                 # Package entry point, exports main functions
├── config.py                   # Centralized paths, thresholds, biosamples
├── schemas.py                  # Nested dict structure helpers (empty_*, ensure_*)
├── utils.py                    # Shared utilities (distance calc, chr handling)
├── main.py                     # Main orchestrator with run_full_pipeline()
│
├── genes/                      # Gene & lncRNA processing
│   ├── __init__.py
│   ├── gene_loader.py          # Load GENCODE, filter, promoter coords
│   ├── lncrna_matching.py      # lncRNA ↔ gene proximity matching
│   └── mirna_targets.py        # TargetScan miRNA target processing
│
├── regulatory_elements/        # cCRE processing
│   ├── __init__.py
│   ├── ccre_loader.py          # Load cCREs, add cell-line signals
│   ├── distance_matching.py    # cCRE ↔ gene distance tiers
│   └── element_table.py        # Build element-focused output table
│
├── evidence/                   # Functional evidence sources
│   ├── __init__.py
│   ├── screen_links.py         # SCREEN experimental + computational
│   ├── abc_links.py            # ABC enhancer-gene predictions
│   ├── hichip_links.py         # HiChIP loop processing
│   └── evidence_merger.py      # Combine all evidence into gene_links
│
├── tad_annotation/             # TAD domain context
│   ├── __init__.py             # Clean exports
│   ├── relations.py            # Pure geometry (overlap, distance)
│   ├── annotator.py            # Feature → TAD (adds TAD_domains column)
│   ├── mirroring.py            # TAD → Feature (adds gene_hits to domains)
│   ├── tad_config.py           # Biosample registry, PAM50 metadata
│   └── loader.py               # Multi-biosample orchestration
│
└── atac_peaks/                 # ATAC peaks processing (NEW)
    ├── __init__.py             # Clean exports
    ├── peak_loader.py          # Load peaks, generate composite IDs
    ├── gene_matching.py        # TSS distance + gene body overlap
    ├── ccre_matching.py        # cCRE overlap + proximity
    ├── tad_annotation.py       # TAD domains + boundary overlaps
    ├── peak_table.py           # Main orchestrator + output saving
    └── annotate_df_with_peaks.py  # Annotate features with nearby peaks
```

---

## Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           REGULATORY ELEMENT PIPELINE                        │
└─────────────────────────────────────────────────────────────────────────────┘

STEP 1: Load & Filter Genes
    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
    │  GENCODE    │     │  lncRNA loci│     │    cCRE     │
    │  (genes +   │     │ (gene_type= │     │    BED      │
    │  panel)     │     │ lncRNA rows)│     │             │
    └──────┬──────┘     └──────┬──────┘     └──────┬──────┘
           │                   │                   │
           │     (lncRNA loci = gene_type lncRNA from same GENCODE table)   │
           │                   │                   │
           ▼                   ▼                   ▼
    ┌─────────────────────────────────────────────────────┐
    │              Harmonize chromosomes                   │
    │              Filter to gene panel                    │
    │              Add promoter coordinates                │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 2: lncRNA Matching
    ┌─────────────────────────────────────────────────────┐
    │     Match lncRNAs within 1Mb of panel genes         │
    │     Output: pairs_df, genes_with_lnc                │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 3: cCRE-Gene Distance Matching
    ┌─────────────────────────────────────────────────────┐
    │     Match cCREs to genes by distance tiers          │
    │     0-100kb │ 100-250kb │ 250-500kb │ 500kb-1Mb    │
    │     Add cell-line signals (H3K27ac, CTCF, etc.)    │
    │     Output: elem_focus table                        │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 4: TAD Annotation ──────────────────────────────────────────────┐
    ┌─────────────────────────────────────────────────────┐          │
    │     For each TAD biosample (24 available):          │          │
    │       • Annotate genes with TAD context             │          │
    │       • Annotate cCREs with TAD context             │          │
    │       • Annotate lncRNAs with TAD context           │          │
    │                                                     │          │
    │     Adds TAD_domains[biosample] = {                 │          │
    │       domains: {id: relation, ...},                 │          │
    │       primary: {domain_id, boundaries, normalized}  │          │
    │     }                                               │          │
    └─────────────────────────────────────────────────────┘          │
           │                                                         │
           │    ┌──────────────────────────────────────────────┐     │
           │    │  TAD Sources (from tads/ batch processing)   │     │
           │    │  ─────────────────────────────────────────   │     │
           │    │  Kim_T47D, Kim_HMEC, Kim_HCC70, Kim_HCC1954  │     │
           │    │  Kim_BT549, Kim_ZR7530, Kim_normal_tissue    │     │
           │    │  Kim_TNBC_tissue1/2/3                        │     │
           │    │  Rao_HMEC                                    │     │
           │    │  LeDily_BT474, LeDily_MCF10A, LeDily_SKBR3   │     │
           │    │  Golloshi_MDA231_control/top10/bottom10      │     │
           │    │  vandenBrand_HB1/HB2, vandenBrand_PB1-5      │     │
           │    └──────────────────────────────────────────────┘     │
           │                                                         │
           ▼                                                         │
STEP 5: SCREEN Evidence ◄────────────────────────────────────────────┘
    ┌─────────────────────────────────────────────────────┐
    │     Experimental links (Hi-C validated)             │
    │     Computational links (predicted)                 │
    │     Per-biosample scores + conservation             │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 6: ABC Evidence
    ┌─────────────────────────────────────────────────────┐
    │     Activity-by-Contact enhancer predictions        │
    │     Map to cCREs, collapse per gene-element link    │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 7: Merge Evidence
    ┌─────────────────────────────────────────────────────┐
    │     Combine SCREEN + ABC into gene_links dict       │
    │     Attach to elem_focus table                      │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 8: HiChIP Evidence
    ┌─────────────────────────────────────────────────────┐
    │     H3K27ac HiChIP loops                            │
    │     Integrate into elem_focus                       │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 9: miRNA Targets
    ┌─────────────────────────────────────────────────────┐
    │     TargetScan predictions for panel genes          │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 10: ATAC Peaks Processing
    ┌─────────────────────────────────────────────────────┐
    │     Load ATAC peaks, generate composite IDs         │
    │     Match to genes (TSS distance + body overlap)    │
    │     Match to cCREs (overlap + proximity)            │
    │     Add TAD domain + boundary overlap annotations   │
    │     Output: atac_peaks_annotated.parquet/csv        │
    └─────────────────────────────────────────────────────┘
           │
           ▼
STEP 11: Save Outputs
    ┌─────────────────────────────────────────────────────┐
    │     regulatory_element_focus_with_evidence.parquet  │
    │     (+ .csv when APM_WRITE_ELEM_FOCUS_CSV=1)         │
    │     Gene tables, BED files, lncRNA matching         │
    │     ATAC peaks annotated table                      │
    └─────────────────────────────────────────────────────┘
```

---

## Output Data Structure

### elem_focus DataFrame (main output)

```
┌──────────────┬────────────┬─────────────────┬────────────────┬──────────────┐
│  ENCODE_id   │   chrom    │  linked_genes   │  TAD_domains   │  gene_links  │
├──────────────┼────────────┼─────────────────┼────────────────┼──────────────┤
│ EH38E123456  │   chr6     │ [HLA-A, HLA-B]  │ {Kim_T47D:     │ {HLA-A:      │
│              │            │                 │   {primary:    │   {screen:   │
│              │            │                 │     {...},     │     {...},   │
│              │            │                 │    domains:    │    ABC:      │
│              │            │                 │     {...}},    │     {...}}}  │
│              │            │                 │  Kim_HCC70:    │              │
│              │            │                 │   {...}, ...}  │              │
└──────────────┴────────────┴─────────────────┴────────────────┴──────────────┘
```

### TAD_domains column structure (per row)

```python
{
    "Kim_T47D": {
        "domains": {"domain_001": "contains", "domain_002": "overlap_right"},
        "primary": {
            "domain_id": "domain_001",
            "rel": "contains",
            "domain": {"chrom": "chr6", "start": 1000000, "end": 2000000, "len": 1000000},
            "feature": {"chrom": "chr6", "start": 1500000, "end": 1502000, "mid": 1501000},
            "boundaries": {
                "left": {"boundary_id": "b_001", "dist_bp": 500000, "overlap": False},
                "right": {"boundary_id": "b_002", "dist_bp": 498000, "overlap": False},
                "nearest": {"side": "right", "dist_bp": 498000}
            },
            "normalized": {
                "frac_from_left": 0.501,
                "frac_from_right": 0.499,
                "frac_to_nearest_boundary": 0.499
            }
        }
    },
    "Kim_HCC70": {...},
    "Rao_HMEC": {...},
    # ... up to 24 biosamples
}
```

---

## ATAC Peaks Output Structure

### Main Table Columns

| Column | Type | Description |
|--------|------|-------------|
| `peak_id` | str | Composite ID: `{chrom}:{start}-{end}` |
| `chrom`, `start`, `end` | str, int, int | Peak coordinates |
| `center`, `length` | int, int | Peak center and length |
| `original_name` | str | Original peak name from input |
| `score`, `annotation`, `percentGC` | float, str, float | Original metadata |
| `linked_genes` | list[str] | All genes within window |
| `genes_by_tier` | dict | `{"0-100kb": [...], ...}` |
| `n_genes_total`, `n_genes_overlapping` | int, int | Gene counts |
| `gene_links` | dict | Per-gene detailed info |
| `n_ccres_total`, `n_ccres_overlapping` | int, int | cCRE counts |
| `ccre_types` | dict | `{"pELS": 3, "dELS": 1}` |
| `ccre_links` | list | Per-cCRE detailed info |
| `TAD_domains` | dict | Per-biosample TAD context |
| `TAD_boundary_overlaps` | dict | Per-biosample boundary overlaps |

### gene_links Structure (ATAC peaks)

```python
{
    "HLA-A": {
        "gene_id": "ENSG00000206503",
        "gene_type": "protein_coding",
        "dist_to_tss": 5000,
        "tier": "0-100kb",
        "tss_position": 29942470,
        "strand": "+",
        "body_overlap": {
            "overlaps": True,
            "overlap_bp": 800,
            "overlap_interval": [29942000, 29942800],
            "overlap_frac_of_peak": 0.12,
            "overlap_frac_of_gene": 0.001,
            "overlap_type": "promoter"  # or "gene_body"
        }
    },
}
```

### ccre_links Structure (ATAC peaks)

```python
[
    {
        "cCRE_id": "EH38E1234567.1",
        "ENCODE_id": "EH38E1234567",
        "raw_type": "pELS,CTCF-bound",
        "distance": 0,
        "overlap": {
            "overlaps": True,
            "overlap_bp": 500,
            "overlap_interval": [6012, 6512],
            "overlap_frac_of_peak": 0.08,
            "overlap_frac_of_ccre": 0.25,
        }
    },
]
```

### TAD_boundary_overlaps Structure (ATAC peaks)

```python
{
    "Kim_T47D": {
        "overlaps_boundary": True,
        "n_boundaries": 1,
        "boundaries": [
            {
                "boundary_id": "Kim_T47D_boundary_042",
                "overlap_bp": 200,
                "overlap_interval": [6012, 6212],
                "boundary_strength": "strong"
            }
        ]
    },
}
```

---

## Usage Examples

### Command Line

```bash
# From repo root on WSL: use the project venv so pandas and other deps resolve
# (bare `python3` may be system Python without those packages).
# .venv/bin/python3 -m pipeline.main

# Full pipeline with all steps
python -m pipeline.main

# Skip slow steps for testing
python -m pipeline.main --skip-hichip --skip-abc --skip-tads --skip-atac

# Skip only TADs
python -m pipeline.main --skip-tads

# Run only ATAC processing
python -m pipeline.main --atac-only --atac-peaks /path/to/peaks.csv

# Custom output directory
python -m pipeline.main --working-dir /path/to/output

# Gene processing only (fastest)
python -m pipeline.main --genes-only

# Evidence collection only
python -m pipeline.main --evidence-only
```

### Python Import - Full Pipeline

```python
from pipeline import run_full_pipeline
from pathlib import Path

# Run everything
elem_focus = run_full_pipeline()

# Custom configuration
elem_focus = run_full_pipeline(
    working_dir=Path("/my/output/dir"),
    gene_panel=["HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "B2M"],
    skip_hichip=True,      # Skip slow HiChIP step
    skip_abc=False,        # Include ABC
    skip_mirna=True,       # Skip miRNA
    skip_tads=False,       # Include TAD annotation
    skip_atac=False,       # Include ATAC processing
    tad_biosamples=["Kim_T47D", "Kim_HCC70", "Rao_HMEC"],  # Specific TADs only
)
```

### Python Import - Modular Functions

```python
# Just gene processing
from pipeline import run_genes_only
genes = run_genes_only()

# Just evidence collection
from pipeline import run_evidence_only
evidence_df = run_evidence_only(include_hichip=False, include_abc=True)

# Just TAD annotation on existing DataFrames
from pipeline.main import run_tad_annotation_only
genes, ccres, lncrnas = run_tad_annotation_only(
    genes_df,
    ccres_df,
    lncrnas_df,
    biosamples=["Kim_T47D", "Rao_HMEC"],
)
```

### Python Import - TAD Annotation Directly

```python
from pipeline.tad_annotation import (
    annotate_with_all_tad_sources,
    annotate_with_pam50_matched,
    discover_tad_sources,
    TAD_BIOSAMPLE_REGISTRY,
    get_biosamples_by_pam50,
)
from pathlib import Path

TAD_DIR = Path("/home/stavz/masters/gdc/TADs/processed")

# Upstream TAD processing details (how `processed/` is produced):
# see pipeline/md/module_specific_processing_md/tads_processing_original.md

# Annotate with ALL available TAD sources (24 biosamples)
genes_df, ccre_df, lncrnas_df = annotate_with_all_tad_sources(
    genes_df,
    ccre_df,
    lncrnas_df,
    processed_dir=TAD_DIR,
)

# Annotate with SPECIFIC biosamples only
genes_df, ccre_df, _ = annotate_with_all_tad_sources(
    genes_df,
    ccre_df,
    processed_dir=TAD_DIR,
    biosamples=["Kim_T47D", "Kim_HCC70", "Rao_HMEC"],
)

# Annotate with PAM50-MATCHED sources (for subtype-specific analysis)
genes_df, ccre_df, _ = annotate_with_pam50_matched(
    genes_df,
    sample_pam50="Basal",  # Uses Kim_HCC70, Kim_BT549, Golloshi_MDA231
    ccre_df=ccre_df,
    processed_dir=TAD_DIR,
)

# Query available sources
available = discover_tad_sources(TAD_DIR)
print(f"Found {len(available)} TAD sources")

# Query biosample metadata
meta = TAD_BIOSAMPLE_REGISTRY["Kim_T47D"]
print(f"PAM50: {meta.pam50}, Study: {meta.study}")

# Get all Basal cell lines
basal_sources = get_biosamples_by_pam50("Basal")
# ['Kim_HCC70', 'Kim_BT549', 'Golloshi_MDA231_control', ...]
```

### Python Import - Low-Level TAD Functions

```python
from pipeline.tad_annotation import (
    # Single-biosample annotation
    annotate_genes_with_tads,
    annotate_ccres_with_tads,
    annotate_svs_with_tads,
    
    # Mirroring (TAD → feature hits)
    mirror_genes_into_domains,
    mirror_all_features_into_domains,
    
    # Query helpers
    get_genes_in_domain,
    get_domains_containing_gene,
    get_domain_gene_count,
    
    # Loading
    load_tad_source,
)
from pathlib import Path

TAD_DIR = Path("/home/stavz/masters/gdc/TADs/processed")

# Load one biosample's TAD data
from pipeline.tad_annotation import TADSourcePaths
paths = TADSourcePaths(biosample="Kim_T47D", base_dir=TAD_DIR / "Kim_T47D")
domains, boundaries, flanks = load_tad_source(paths)

# Annotate genes with single biosample
genes_df = annotate_genes_with_tads(
    genes_df,
    tad_domains=domains,
    domain_flanks=flanks,
    boundaries=boundaries,
    biosample="Kim_T47D",
)

# Mirror back into domains
domains = mirror_genes_into_domains(domains, genes_df, biosample="Kim_T47D")

# Query: what genes are in domain X?
genes_in_domain = get_genes_in_domain(domains, "Kim_T47D_domain_042", "Kim_T47D")

# Query: which domains contain gene Y?
domains_with_hla = get_domains_containing_gene(domains, "HLA-A", "Kim_T47D")
```

### Python Import - SV Annotation (Sample-Specific)

```python
from pipeline.tad_annotation import annotate_svs_with_all_tad_sources
from pathlib import Path

TAD_DIR = Path("/home/stavz/masters/gdc/TADs/processed")

# Load your SV DataFrame (from MANTA output)
sv_df = pd.read_csv("sample_svs.vcf.csv")

# Annotate with TAD context
sv_df = annotate_svs_with_all_tad_sources(
    sv_df,
    processed_dir=TAD_DIR,
    biosamples=["Kim_T47D"],  # Use PAM50-matched for this sample
    ins_mode="point",         # Treat insertions as points
)

# Now sv_df["TAD_domains"]["Kim_T47D"] contains:
# - Which domain each SV falls in
# - Distance to nearest boundary
# - Whether SV overlaps a boundary (potential disruption!)
```

### Python Import - ATAC Peaks Processing

```python
from pipeline.main import run_atac_only
from pipeline.atac_peaks import (
    load_atac_peaks,
    match_peaks_to_genes,
    match_peaks_to_ccres,
    build_gene_links,
    build_ccre_links,
)
from pathlib import Path

# Quick ATAC processing with defaults
atac_table = run_atac_only(
    atac_peaks_path=Path("/path/to/peaks.csv"),
    include_tads=True,
    tad_biosamples=["Kim_T47D"],
)

# Or step-by-step for more control
peaks = load_atac_peaks("/path/to/peaks.csv")
peak_gene_pairs = match_peaks_to_genes(peaks, genes, window_bp=500_000)
peak_ccre_pairs = match_peaks_to_ccres(peaks, ccres, max_distance=0)
gene_links_df = build_gene_links(peak_gene_pairs)
ccre_links_df = build_ccre_links(peak_ccre_pairs)
```

### Jupyter Notebook Quick Start

```python
# Cell 1: Setup
import sys
sys.path.insert(0, '/home/stavz/masters/gdc')  # Parent of pipeline/

from pipeline import run_full_pipeline
from pipeline.config import PATHS, PRIMARY_GENES
from pathlib import Path

# Cell 2: Run pipeline
elem_focus = run_full_pipeline(
    skip_hichip=True,  # Skip for faster iteration
    skip_tads=False,
    tad_biosamples=["Kim_T47D", "Rao_HMEC"],
)

# Cell 3: Explore results
print(f"Total elements: {len(elem_focus)}")
print(f"Columns: {elem_focus.columns.tolist()}")

# Cell 4: Check TAD annotations
sample_row = elem_focus.iloc[0]
print(sample_row["TAD_domains"].keys())  # ['Kim_T47D', 'Rao_HMEC']
```

---

## Query Helpers for ATAC Peaks

```python
from pipeline.atac_peaks import (
    get_peaks_near_gene,
    get_peaks_overlapping_gene,
    get_peaks_at_promoters,
    get_peaks_at_boundaries,
)

# Get peaks within 100kb of HLA-A
hla_peaks = get_peaks_near_gene(atac_table, "HLA-A", max_distance=100_000)

# Get peaks overlapping HLA-A gene body
overlap_peaks = get_peaks_overlapping_gene(atac_table, "HLA-A")

# Get peaks at promoters of APM genes
promoter_peaks = get_peaks_at_promoters(
    atac_table,
    gene_names=["HLA-A", "HLA-B", "TAP1", "TAP2"]
)

# Get peaks at TAD boundaries
boundary_peaks = get_peaks_at_boundaries(atac_table, biosample="Kim_T47D")
```

---

## Configuration Reference

### config.py Additions

```python
@dataclass
class PathConfig:
    # ... existing paths ...
    
    # TADs processed directory (output from tads/ batch processing)
    tads_processed: Path = Path("/home/stavz/masters/gdc/TADs/processed")
    
    # ATAC peaks
    atac_peaks_csv: Path = Path("/path/to/tcga_atac_peaks.csv")

@dataclass  
class ThresholdConfig:
    # ... existing thresholds ...
    
    # ATAC peak thresholds
    atac_gene_window_bp: int = 1_000_000
    atac_ccre_max_distance: int = 0  # 0 = overlap only

OUTPUT_SUBDIRS = {
    # ... existing subdirs ...
    "atac_peaks": "atac_peaks",
}
```

### TAD Biosample Registry (in tad_config.py)

| Biosample | Study | Cell Line | PAM50 | Tissue Type |
|-----------|-------|-----------|-------|-------------|
| Kim_T47D | Kim | T47D | LumA | cell_line |
| Kim_HCC70 | Kim | HCC70 | Basal | cell_line |
| Kim_HCC1954 | Kim | HCC1954 | HER2 | cell_line |
| Kim_BT549 | Kim | BT549 | Basal | cell_line |
| Kim_ZR7530 | Kim | ZR75-30 | LumA | cell_line |
| Kim_HMEC | Kim | HMEC | - | normal_tissue |
| Kim_TNBC_tissue1/2/3 | Kim | - | Basal | tumor_tissue |
| Rao_HMEC | Rao | HMEC | - | normal_tissue |
| LeDily_BT474 | LeDily | BT474 | LumB | cell_line |
| LeDily_MCF10A | LeDily | MCF10A | - | cell_line |
| LeDily_SKBR3 | LeDily | SKBR3 | HER2 | cell_line |
| Golloshi_MDA231_* | Golloshi | MDA-MB-231 | Basal | cell_line |
| vandenBrand_HB1/2 | vandenBrand | - | - | healthy_breast |
| vandenBrand_PB1-5 | vandenBrand | - | - | tumor_tissue |

### PAM50 Reference Map (for matched analysis)

```python
PAM50_REFERENCE_MAP = {
    "Basal": ["Kim_HCC70", "Kim_BT549", "Golloshi_MDA231_control"],
    "LumA": ["Kim_T47D", "Kim_ZR7530"],
    "LumB": ["LeDily_BT474"],
    "HER2": ["Kim_HCC1954", "LeDily_SKBR3"],
    "Normal": ["Rao_HMEC", "Kim_HMEC", "LeDily_MCF10A"],
}
```

---

## Output Files

### From Full Pipeline

```
working_dir/
├── apm_genes.bed
├── primary_genes_all_features.csv   # PIPELINE_GENE_PANEL multifeature (legacy basename)
├── primary_genes_only.csv             # PIPELINE_GENE_PANEL gene rows + promoters
├── cnv_genes.csv                      # CNV_GENES union
├── lncRNAs_genes_all_features.csv
├── tier2_medium_genes_only.csv
├── tier2_medium_genes_all_features.csv
├── tier3_cnv_only_genes_only.csv
├── tier3_cnv_only_genes_all_features.csv
├── tier4_readout_genes_only.csv
├── tier4_readout_genes_all_features.csv
├── lncRNA_matching/
│   ├── lncrna_gene_pairs.csv
│   ├── genes_with_lncrnas.csv
│   └── lncrnas_with_genes.csv
├── regulatory_elements_matching/
│   ├── coding_element_focus.csv
│   ├── coding_gene_summary.csv
│   ├── coding_gene_to_elements.csv
│   ├── lncRNA_element_focus.csv
│   ├── lncRNA_gene_summary.csv
│   ├── regulatory_element_focus_with_evidence.parquet   # primary (JSON-encoded nested cols)
│   └── regulatory_element_focus_with_evidence.csv       # only when APM_WRITE_ELEM_FOCUS_CSV=1
├── miRNA/
│   ├── mirna_targets.csv
│   └── mirna_summary.csv
├── evidence/
│   ├── screen_exp_links.csv
│   ├── screen_comp_links.csv
│   ├── abc_links.csv
│   └── hichip_links.csv
└── atac_peaks/
    ├── atac_peaks_annotated.parquet
    ├── atac_peaks_annotated.csv
    ├── peak_id_mapping.csv
    └── summary.json
```
