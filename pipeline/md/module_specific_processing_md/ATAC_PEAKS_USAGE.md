# ATAC Peaks Module - Usage Guide

## Module Structure

```
atac_peaks/
├── __init__.py                  # Clean exports
├── peak_loader.py               # Load and normalize peaks, generate IDs
├── gene_matching.py             # TSS distance + gene body overlap
├── ccre_matching.py             # cCRE overlap + proximity
├── tad_annotation.py            # TAD domains + boundary overlaps
├── peak_table.py                # Main orchestrator + output saving
└── annotate_df_with_peaks.py    # Annotate features (genes/cCREs/SVs/SNVs) with peaks
```

## Quick Start

```python
from pathlib import Path
from atac_peaks import build_atac_peak_table, save_atac_outputs

# Load your existing gene and cCRE tables
genes = pd.read_csv("genes.csv")
ccres = pd.read_csv("ccres.csv")

# Build complete peak annotation table
peak_table = build_atac_peak_table(
    peaks_path="atac_peaks.csv",
    genes=genes,
    ccres=ccres,
    gene_panel=PRIMARY_GENES,          # Optional: filter to specific genes
    gene_window_bp=1_000_000,          # TSS distance window
    ccre_max_distance=0,               # 0 = overlap only, >0 includes proximity
    tad_processed_dir="/path/to/tads", # Optional: add TAD annotations
    tad_biosamples=["Kim_T47D"],       # Optional: specific cell lines
)

# Save outputs
save_atac_outputs(peak_table, output_dir="output/atac/")
```

## Output Table Structure

### Core Columns
| Column | Type | Description |
|--------|------|-------------|
| `peak_id` | str | Composite ID: `{chrom}:{start}-{end}` |
| `chrom` | str | Chromosome (with 'chr' prefix) |
| `start`, `end` | int | Peak coordinates |
| `center` | int | Peak center position |
| `length` | int | Peak length in bp |
| `original_name` | str | Original peak name from input file |
| `score` | float | Peak score (if present in input) |
| `annotation` | str | Original annotation (if present) |
| `percentGC` | float | GC content (if present) |

### Gene Link Columns
| Column | Type | Description |
|--------|------|-------------|
| `linked_genes` | list[str] | All genes linked to this peak |
| `genes_by_tier` | dict | `{"0-100kb": [...], "100-250kb": [...], ...}` |
| `n_genes_total` | int | Total linked genes |
| `n_genes_overlapping` | int | Genes with body overlap |
| `gene_links` | dict | Detailed per-gene info (see below) |

### gene_links Structure
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
            "overlap_type": "promoter"  # or "gene_body", None
        }
    },
    "HLA-B": {...},
}
```

### cCRE Link Columns
| Column | Type | Description |
|--------|------|-------------|
| `n_ccres_total` | int | Total linked cCREs |
| `n_ccres_overlapping` | int | cCREs with direct overlap |
| `ccre_types` | dict | `{"pELS": 3, "dELS": 1, ...}` |
| `ccre_links` | list | List of linked cCRE dicts (see below) |

### ccre_links Structure
```python
[
    {
        "cCRE_id": "EH38E1234567.1",
        "ENCODE_id": "EH38E1234567",
        "raw_type": "pELS,CTCF-bound",
        "distance": 0,  # 0 if overlap
        "overlap": {
            "overlaps": True,
            "overlap_bp": 500,
            "overlap_interval": [6012, 6512],
            "overlap_frac_of_peak": 0.08,
            "overlap_frac_of_ccre": 0.25,
        }
    },
    {...},
]
```

### TAD Annotation Columns
| Column | Type | Description |
|--------|------|-------------|
| `TAD_domains` | dict | Per-biosample TAD domain context |
| `TAD_boundary_overlaps` | dict | Per-biosample boundary overlaps |

### TAD_domains Structure
Same as your existing pipeline - keyed by biosample:
```python
{
    "Kim_T47D": {
        "domains": {"domain_001": "contains", ...},
        "primary": {
            "domain_id": "domain_001",
            "rel": "contains",
            "domain": {"chrom": "chr6", "start": 1000000, "end": 2000000, "len": 1000000},
            "feature": {"chrom": "chr6", "start": 1500000, "end": 1502000, "mid": 1501000},
            "boundaries": {"left": {...}, "right": {...}, "nearest": {...}},
            "normalized": {"frac_from_left": 0.501, ...}
        }
    },
    "Kim_HCC70": {...},
}
```

### TAD_boundary_overlaps Structure (NEW)
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
    "Kim_HCC70": {
        "overlaps_boundary": False,
        "n_boundaries": 0,
        "boundaries": []
    },
}
```

---

## Modular Usage

### Load Peaks Only
```python
from atac_peaks import load_atac_peaks, create_peak_id_mapping

peaks = load_atac_peaks("peaks.csv")
mapping = create_peak_id_mapping(peaks, output_path="peak_id_mapping.csv")
```

### Gene Matching Only
```python
from atac_peaks import match_peaks_to_genes, build_gene_links

# Get all peak-gene pairs
pairs = match_peaks_to_genes(
    peaks, genes,
    window_bp=500_000,
    promoter_upstream=2000,
    promoter_downstream=500,
)

# Build nested gene_links dict
gene_links_df = build_gene_links(pairs)
```

### cCRE Matching Only
```python
from atac_peaks import match_peaks_to_ccres, build_ccre_links

# Overlap only (max_distance=0)
pairs = match_peaks_to_ccres(peaks, ccres, max_distance=0)

# Include proximity (e.g., within 1kb)
pairs = match_peaks_to_ccres(peaks, ccres, max_distance=1000)

# Build nested ccre_links list
ccre_links_df = build_ccre_links(pairs)
```

### TAD Annotation Only
```python
from atac_peaks import (
    annotate_peaks_with_tads,
    annotate_peaks_with_boundary_overlaps,
)

# Single biosample
peaks = annotate_peaks_with_tads(
    peaks,
    tad_domains=domains_df,
    domain_flanks=flanks_df,
    boundaries=bounds_df,
    biosample="Kim_T47D",
)

peaks = annotate_peaks_with_boundary_overlaps(
    peaks,
    boundaries=bounds_df,
    biosample="Kim_T47D",
)

# Or use batch annotation for multiple biosamples
from atac_peaks.tad_annotation import annotate_peaks_with_all_tad_sources

peaks = annotate_peaks_with_all_tad_sources(
    peaks,
    processed_dir="/path/to/tads",
    biosamples=["Kim_T47D", "Kim_HCC70", "Rao_HMEC"],
)
```

---

## Query Helpers

### Find Peaks Near Specific Genes
```python
from atac_peaks import get_peaks_near_gene, get_peaks_overlapping_gene

# All peaks within 100kb of HLA-A
hla_a_peaks = get_peaks_near_gene(peak_table, "HLA-A", max_distance=100_000)

# Peaks overlapping HLA-A gene body
hla_a_overlap = get_peaks_overlapping_gene(peak_table, "HLA-A")
```

### Find Peaks at Promoters
```python
from atac_peaks import get_peaks_at_promoters

# All peaks at any gene's promoter
promoter_peaks = get_peaks_at_promoters(peak_table)

# Peaks at specific genes' promoters
apm_promoter_peaks = get_peaks_at_promoters(
    peak_table, 
    gene_names=["HLA-A", "HLA-B", "TAP1", "TAP2"]
)
```

### Find Peaks at TAD Boundaries
```python
from atac_peaks.tad_annotation import get_peaks_at_boundaries, summarize_boundary_overlaps

# Get peaks overlapping boundaries in any biosample
boundary_peaks = get_peaks_at_boundaries(peak_table)

# Get peaks at boundaries in specific cell line
boundary_peaks_t47d = get_peaks_at_boundaries(peak_table, biosample="Kim_T47D")

# Summary statistics
summary = summarize_boundary_overlaps(peak_table)
print(summary)
```

### Bidirectional cCRE Lookups
```python
from atac_peaks.ccre_matching import get_ccres_for_peak, get_peaks_for_ccre

# What cCREs does this peak overlap?
ccres = get_ccres_for_peak(peak_ccre_pairs, "chr6:29942000-29943000")

# What peaks overlap this cCRE?
peaks = get_peaks_for_ccre(peak_ccre_pairs, "EH38E1234567", id_col="ENCODE_id")
```

---

## Integration with Main Pipeline

To add this module to your existing pipeline, place the `atac_peaks/` folder in your pipeline directory:

```
pipeline/
├── __init__.py
├── config.py
├── main.py
├── ...
├── atac_peaks/           # <-- NEW
│   ├── __init__.py
│   ├── schemas.py
│   ├── peak_loader.py
│   ├── gene_matching.py
│   ├── ccre_matching.py
│   ├── tad_annotation.py
│   └── peak_table.py
```

Then import from the module:
```python
from pipeline.atac_peaks import build_atac_peak_table
```

### Add to config.py
```python
@dataclass
class PathConfig:
    # ... existing paths ...
    
    # ATAC peaks
    atac_peaks_csv: Path = Path("/path/to/atac_peaks.csv")
```

### Add to main.py (optional step)
```python
# After STEP 4 (TAD annotation) or wherever appropriate:

if not skip_atac:
    print("\n" + "-" * 40)
    print("STEP X: Processing ATAC peaks")
    print("-" * 40)
    
    from .atac_peaks import build_atac_peak_table, save_atac_outputs
    
    atac_table = build_atac_peak_table(
        peaks_path=PATHS.atac_peaks_csv,
        genes=genes,
        ccres=ccres,
        gene_panel=gene_panel,
        tad_processed_dir=PATHS.tads_processed,
        tad_biosamples=tad_biosamples,
    )
    
    save_atac_outputs(
        atac_table,
        output_dir=working_dir / "atac_peaks",
    )
```

---

## Output Files

When calling `save_atac_outputs()`:

| File | Format | Description |
|------|--------|-------------|
| `atac_peaks_annotated.parquet` | Parquet | Full table with nested dicts preserved |
| `atac_peaks_annotated.csv` | CSV | Full table with dicts as JSON strings |
| `peak_id_mapping.csv` | CSV | `peak_id ↔ original_name` lookup |
| `peak_gene_pairs.csv` | CSV | All peak-gene pairs (optional) |
| `peak_ccre_pairs.csv` | CSV | All peak-cCRE pairs (optional) |
| `summary.json` | JSON | Summary statistics |

---

## Annotating Features with Peaks (Reverse Direction)

The module also supports annotating other feature types (genes, cCREs, SVs, SNVs) with nearby ATAC peaks. This is the reverse of the main pipeline flow.

### Basic Usage

```python
from pipeline.atac_peaks import (
    load_atac_peaks,
    annotate_df_with_peaks,
    annotate_genes_with_peaks,
    annotate_ccres_with_peaks,
    annotate_svs_with_peaks,
    annotate_snvs_with_peaks,
)

# Load peaks
peaks = load_atac_peaks("/path/to/atac_peaks.csv")

# Annotate genes (uses strand for directional distance)
genes = annotate_genes_with_peaks(genes_df, peaks, window_bp=100_000)

# Annotate cCREs
ccres = annotate_df_with_peaks(ccres_df, peaks, kind="ccre", window_bp=50_000)

# Annotate SVs (handles DEL/DUP as intervals, INS/BND as points)
svs = annotate_svs_with_peaks(sv_df, peaks, window_bp=100_000)

# Annotate SNVs (point features using 'pos' column)
snvs = annotate_snvs_with_peaks(snv_df, peaks, window_bp=50_000)
```

### Feature Types Supported

| Kind | Input Columns | Behavior |
|------|--------------|----------|
| `"gene"`, `"ccre"`, `"lncrna"`, `"interval"` | `chrom`, `start`, `end` | Uses full interval |
| `"sv"` | `chrom`, `pos`, `END`, `SVTYPE`, `SVLEN` | DEL/DUP → interval; INS/BND → point at `pos` |
| `"snv"` | `chrom`, `pos` | Point feature at `pos` |
| `"point"` | `chrom`, `pos` or `position` or `start` | Point feature |

### Output: `atac_peak_links` Column

Each feature gets a list of dicts sorted by absolute distance:

```python
[
    {
        # Peak metadata
        "peak_id": "chr6:29000000-29001000",
        "chrom": "chr6",
        "start": 29000000,
        "end": 29001000,
        "center": 29000500,
        "score": 5.2,
        "percentGC": 0.45,
        
        # Spatial relationship
        "signed_dist": -5000,       # Negative = upstream (for + strand)
        "abs_dist": 5000,           # Absolute distance
        "overlap_bp": 0,            # Overlap in bp (0 if no overlap)
        "overlaps": False,          # Boolean overlap flag
        "upstream_5kb_flag": True,  # Peak within 5kb upstream
        "downstream_5kb_flag": False,
    },
    ...
]
```

### Distance Convention

- **`signed_dist`**: Distance from feature center to peak center
  - For **+ strand** features: negative = upstream (peak is 5' of feature)
  - For **- strand** features: signs are flipped
  - For features without strand: raw distance (positive = peak to the right)

- **`upstream_5kb_flag`** / **`downstream_5kb_flag`**: Configurable thresholds (default 5kb)
  - Only set for non-overlapping peaks
  - Respects strand orientation

### Parameters

```python
annotate_df_with_peaks(
    df,                            # Feature DataFrame
    peaks,                         # ATAC peaks DataFrame  
    kind="gene",                   # Feature type
    window_bp=100_000,             # Max distance from feature center to peak center
    out_col="atac_peak_links",     # Output column name
    strand_col="strand",           # Column for strand info (None to ignore)
    upstream_threshold=5000,       # Threshold for upstream flag
    downstream_threshold=5000,     # Threshold for downstream flag
)
```

### Query Helpers

```python
from pipeline.atac_peaks import (
    get_features_with_overlapping_peaks,
    get_features_with_upstream_peaks,
    get_features_with_downstream_peaks,
    get_closest_peak,
    summarize_peak_links,
    flatten_peak_links,
)

# Filter to features with overlapping peaks
overlapping = get_features_with_overlapping_peaks(genes)

# Filter to features with upstream peaks (within threshold)
upstream = get_features_with_upstream_peaks(genes)

# Get summary statistics per feature
summary = summarize_peak_links(genes)
# Returns: n_peaks, n_overlapping, n_upstream, n_downstream, 
#          min_dist, closest_peak_id, closest_peak_score

# Flatten to long-form DataFrame (one row per feature-peak pair)
pairs = flatten_peak_links(genes, id_col="gene_name")

# Get closest peak from a links list
closest = get_closest_peak(genes.iloc[0]["atac_peak_links"])
```

### Example: SVs Near Open Chromatin

```python
# Load SVs from MANTA output
svs = pd.read_csv("manta_svs.csv")

# Annotate with nearby peaks
svs = annotate_svs_with_peaks(svs, peaks, window_bp=50_000)

# Find SVs that overlap ATAC peaks (potential regulatory disruption)
sv_at_peaks = get_features_with_overlapping_peaks(svs)
print(f"SVs overlapping ATAC peaks: {len(sv_at_peaks)}")

# Get summary
summary = summarize_peak_links(svs)
svs = pd.concat([svs, summary], axis=1)
```

---

## Notes

1. **Peak ID**: Composite `{chrom}:{start}-{end}` ensures uniqueness and enables coordinate-based lookups

2. **Gene body overlap types**: Currently classifies as "promoter" or "gene_body". More detailed classification (exon/intron/UTR) would require transcript-level annotation

3. **cCRE matching**: Set `max_distance=0` for overlap-only, or increase for proximity matching

4. **TAD annotations**: Uses the same structure as your existing pipeline, just adds the new `TAD_boundary_overlaps` column

5. **Parquet recommended**: Use parquet format to preserve nested dict structures; CSV converts them to JSON strings
