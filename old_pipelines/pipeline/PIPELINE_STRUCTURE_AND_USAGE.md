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
└── tad_annotation/             # TAD domain context (NEW)
    ├── __init__.py             # Clean exports
    ├── relations.py            # Pure geometry (overlap, distance)
    ├── annotator.py            # Feature → TAD (adds TAD_domains column)
    ├── mirroring.py            # TAD → Feature (adds gene_hits to domains)
    ├── tad_config.py           # Biosample registry, PAM50 metadata
    └── loader.py               # Multi-biosample orchestration
```

---

## Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           REGULATORY ELEMENT PIPELINE                        │
└─────────────────────────────────────────────────────────────────────────────┘

STEP 1: Load & Filter Genes
    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
    │  GENCODE    │     │   lncRNA    │     │    cCRE     │
    │    GTF      │     │    CSV      │     │    BED      │
    └──────┬──────┘     └──────┬──────┘     └──────┬──────┘
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
STEP 10: Save Outputs
    ┌─────────────────────────────────────────────────────┐
    │     regulatory_element_focus_with_evidence.csv      │
    │     Gene tables, BED files, lncRNA matching         │
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

## Usage Examples

### Command Line

```bash
# Full pipeline with all steps
python -m pipeline.main

# Skip slow steps for testing
python -m pipeline.main --skip-hichip --skip-abc --skip-tads

# Skip only TADs
python -m pipeline.main --skip-tads

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

## Configuration Reference

### Add to config.py

```python
@dataclass
class PathConfig:
    # ... existing paths ...
    
    # TADs processed directory (output from tads/ batch processing)
    tads_processed: Path = Path("/home/stavz/masters/gdc/TADs/processed")
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
