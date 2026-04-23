# =============================================================================
# SNV MODULE DOCUMENTATION
# =============================================================================
# Add this section to your PIPELINE_STRUCTURE_AND_USAGE.md file


## SNV Processing Module

### Directory Structure

```
snv/
├── __init__.py             # Clean exports
├── vcf_loader.py           # Main VCF loading and orchestration
├── somatic_filter.py       # High-confidence somatic filtering
├── vep_parser.py           # VEP CSQ annotation parsing
└── ccre_matching.py        # cCRE overlap detection
```

### Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           SNV PROCESSING PIPELINE                            │
└─────────────────────────────────────────────────────────────────────────────┘

INPUT: VEP-annotated Mutect2 VCF files
       ├── Raw VCF from TCGA/GDC
       └── Annotated with VEP --everything

STEP 1: VCF Parsing (vcf_loader.py)
    ┌─────────────────────────────────────────────────────────────────────────┐
    │     • Parse VCF with vcfpy                                              │
    │     • Clean non-standard VEP header lines                               │
    │     • Expand multi-allelic sites (one row per ALT)                      │
    │     • Extract per-sample AD/DP values                                   │
    │     • Compute tumor_vaf and normal_vaf                                  │
    └─────────────────────────────────────────────────────────────────────────┘
           │
           ▼
STEP 2: Somatic Filtering (somatic_filter.py)
    ┌─────────────────────────────────────────────────────────────────────────┐
    │     Apply high-confidence somatic masks:                                │
    │       • FILTER == "PASS"                                                │
    │       • tumor_vaf >= 5%                                                 │
    │       • normal_vaf <= 2%                                                │
    │       • TLOD >= 6.0 (if present)                                        │
    │       • POPAF >= 3.0 (population AF ≤ 1e-3)                             │
    │       • Read depth thresholds                                           │
    └─────────────────────────────────────────────────────────────────────────┘
           │
           ▼
STEP 3: VEP Annotation Parsing (vep_parser.py)
    ┌─────────────────────────────────────────────────────────────────────────┐
    │     Parse CSQ field into structured hits:                               │
    │       • gene_hits: Transcript-level consequences                        │
    │       • regulatory_hits: RegulatoryFeature consequences                 │
    │       • motif_hits: MotifFeature (TF binding) consequences              │
    │                                                                         │
    │     Compute impact flags:                                               │
    │       • has_missense, has_nonsense, has_frameshift, has_splice_effect   │
    │       • *_canonical variants (CANONICAL == "1")                         │
    │       • *_mane variants (MANE_SELECT present)                           │
    └─────────────────────────────────────────────────────────────────────────┘
           │
           ▼
STEP 4: cCRE Matching (ccre_matching.py)
    ┌─────────────────────────────────────────────────────────────────────────┐
    │     Match variants to overlapping cCREs:                                │
    │       • Point overlap (SNV pos within cCRE interval)                    │
    │       • Captures linked genes via genes_by_exact_dist                   │
    │       • Enables regulatory impact analysis                              │
    └─────────────────────────────────────────────────────────────────────────┘
           │
           ▼
OUTPUT: DataFrame with annotated somatic variants
```

### Output DataFrame Columns

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `pos` | int | 1-based position |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele |
| `qual` | float | Variant quality |
| `filter` | str | FILTER status |
| `tumor_vaf` | float | Tumor variant allele frequency |
| `normal_vaf` | float | Normal VAF (should be ~0 for somatic) |
| `gene_hits` | list[dict] | VEP transcript-level hits |
| `regulatory_hits` | list[dict] | VEP regulatory feature hits |
| `motif_hits` | list[dict] | VEP motif feature hits |
| `gene_symbols` | str | Comma-separated affected genes |
| `hits_canonical` | bool | Hits any canonical transcript |
| `has_missense` | bool | Has missense consequence |
| `has_nonsense` | bool | Has stop-gained consequence |
| `has_frameshift` | bool | Has frameshift consequence |
| `has_splice_effect` | bool | Has splice region/acceptor/donor |
| `has_*_canonical` | bool | Impact on canonical transcripts |
| `has_*_mane` | bool | Impact on MANE SELECT transcripts |
| `cCRE_hits` | list[dict] | Overlapping regulatory elements |

### gene_hits Structure

```python
{
    "Allele": "T",
    "Consequence": "missense_variant",
    "IMPACT": "MODERATE",
    "SYMBOL": "HLA-A",
    "Gene": "ENSG00000206503",
    "Feature": "ENST00000376809",
    "BIOTYPE": "protein_coding",
    "EXON": "3/8",
    "HGVSc": "c.123A>T",
    "HGVSp": "p.Lys41Asn",
    "CANONICAL": "1",
    "MANE_SELECT": "NM_002116.8",
    "SIFT": "deleterious(0.01)",
    "PolyPhen": "probably_damaging(0.99)",
    "gnomADe_AF": "0.0001",
    "CLIN_SIG": "pathogenic",
}
```

### cCRE_hits Structure

```python
[
    {
        "cCRE_id": "EH38E1234567.1",
        "elem_type": "pELS,CTCF-bound",
        "chrom": "chr6",
        "elem_start": 29942000,
        "elem_end": 29943000,
        "genes_by_exact_dist": "HLA-A:5000,HLA-B:120000",
    },
]
```

### Usage Examples

#### Single VCF Loading

```python
from pipeline.snv import load_mutect_snv_vcf
from pipeline.config import PATHS, PRIMARY_GENES

# Load single VCF
df, normal, tumor = load_mutect_snv_vcf(
    vcf_path="/path/to/sample.vep.vcf",
    primary_genes=PRIMARY_GENES,
    elements_path=PATHS.regulatory_elements_table,
)

print(f"Loaded {len(df)} somatic variants")
print(f"Normal sample: {normal}, Tumor sample: {tumor}")
```

#### Batch Processing

```python
from pipeline.snv import load_mutect_snv_batch

# Load all VCFs in directory
combined_df = load_mutect_snv_batch(
    vcf_dir="/path/to/vep_vcfs/",
    primary_genes=PRIMARY_GENES,
    elements_path=PATHS.regulatory_elements_table,
    pattern="*.vep.vcf",
)

# Group by sample
by_sample = combined_df.groupby("source_file")
```

#### Custom Filtering

```python
from pipeline.snv import load_mutect_snv_vcf
from pipeline.config import THRESHOLDS

# Stricter filtering
df, _, _ = load_mutect_snv_vcf(
    vcf_path="/path/to/sample.vcf",
    primary_genes=PRIMARY_GENES,
    filter_kwargs={
        "min_tumor_vaf": 0.10,      # 10% VAF minimum
        "min_tumor_dp": 30,          # 30x depth
        "require_tlod": True,        # TLOD must be present
    },
)
```

#### Filtering Convenience Functions

```python
from pipeline.snv import (
    get_coding_variants,
    get_splice_variants,
    get_regulatory_variants,
    get_canonical_variants,
)

# Get different variant classes
coding = get_coding_variants(df)       # missense/nonsense/frameshift
splice = get_splice_variants(df)       # splice region/donor/acceptor
regulatory = get_regulatory_variants(df)  # cCRE or VEP regulatory hits
canonical = get_canonical_variants(df)    # affecting canonical transcripts
```

#### cCRE-Based Filtering

```python
from pipeline.snv.ccre_matching import (
    filter_snvs_by_ccre_type,
    get_snvs_in_enhancers,
    get_snvs_in_promoters,
    summarize_ccre_hits,
)

# Get variants in enhancers
enhancer_vars = get_snvs_in_enhancers(df)

# Get variants in promoters
promoter_vars = get_snvs_in_promoters(df)

# Custom filtering by cCRE type
ctcf_vars = filter_snvs_by_ccre_type(df, ccre_types=["CTCF-only", "CTCF-bound"])

# Add summary columns
df = summarize_ccre_hits(df)
```

### Integration with Main Pipeline

The SNV module is designed to work with the regulatory element table output:

```python
from pipeline import run_full_pipeline
from pipeline.snv import load_mutect_snv_batch
from pipeline.config import PATHS, PRIMARY_GENES

# 1. Run main pipeline to generate regulatory elements table
elem_focus = run_full_pipeline(skip_hichip=True)

# 2. Use that table for SNV annotation
snv_df = load_mutect_snv_batch(
    vcf_dir=PATHS.snv_vcf_dir,
    primary_genes=PRIMARY_GENES,
    elements_path=PATHS.working_dir / "regulatory_elements_matching/regulatory_element_focus_with_evidence.csv",
)

# 3. Now each variant has cCRE_hits linking back to the element table
```

### VEP Command Reference

For optimal results, run VEP with these options:

```bash
vep \
  --input_file sample.vcf \
  --output_file sample.vep.vcf \
  --offline \
  --cache \
  --dir_cache ~/.vep \
  --species homo_sapiens \
  --assembly GRCh38 \
  --fasta ~/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --format vcf \
  --vcf \
  --force_overwrite \
  --everything \           # Include all annotations
  --symbol \               # Gene symbols
  --variant_class \        # Variant class annotation
  --check_ref \            # Check reference alleles
  --fork 4
```

### Configuration Reference

Add to `config.py`:

```python
# In PathConfig:
snv_vcf_dir: Path = Path("/path/to/vep_vcfs")
regulatory_elements_table: Path = Path("/path/to/regulatory_element_focus_with_evidence.csv")

# In ThresholdConfig:
snv_min_tumor_vaf: float = 0.05
snv_max_normal_vaf: float = 0.02
snv_min_tlod: float = 6.0
snv_min_popaf: float = 3.0
snv_min_tumor_dp: int = 20
snv_min_normal_dp: int = 10

# In OUTPUT_SUBDIRS:
"snv": "snv_variants",
```
