# miRTarBase Module Documentation

## Overview

`pipeline/genes/mirtarbase.py` processes experimentally-validated miRNA–gene interactions from miRTarBase. It normalizes the raw download, assigns miRNA family labels (TargetScan seed families with a name-based fallback), and collapses multi-row records to study-level granularity. The output is five CSV tables at increasing levels of aggregation, all written to `{working_dir}/miRNA/mirtarbase/`.

This module complements the existing TargetScan-based `mirna_targets.py` (computational predictions) by adding the experimental evidence layer. Together they cover both arms of miRNA target annotation: predicted binding affinity and validated functional impact.

---

## Data Flow

```
miRTarBase CSV  +  TargetScan miR_Family_Info.txt
        │                       │
        ▼                       ▼
   Column detection       Family lookup (human only)
        │                       │
        └───────────┬───────────┘
                    ▼
          Normalize & filter to gene panel
          (support types, miRNA names, gene symbols)
                    │
                    ▼
          Explode study_ids (one row per PMID)
                    │
                    ▼
    ┌───────────────────────────────────────┐
    │  Collapse to (miRNA, gene, study)     │
    │  Binary flags: support × experiment   │
    └───────────────┬───────────────────────┘
                    │
        ┌───────────┼───────────┬───────────┐
        ▼           ▼           ▼           ▼
  Interaction    Gene       miRNA       Family
   Summary     Summary    Summary     Summary
  (miRNA×gene) (gene)    (miRNA)    (family)
```

---

## Input Files

| File | Config key | Description |
|------|-----------|-------------|
| `data/miRNA/mirtar.csv` | `PATHS.mirtarbase_csv` | miRTarBase full download (all species; filtered to gene panel at runtime) |
| `data/miRNA/miR_Family_Info.txt` | `PATHS.mir_family_info` | TargetScan family annotation (tab-separated); filtered to human (Species ID 9606). Optional — falls back to name-based proxy if absent |

---

## Output Files

All written to `{working_dir}/miRNA/mirtarbase/`:

```
mirtarbase/
├── mirtar_interaction_study_collapsed.csv
├── mirtar_interaction_summary.csv
├── mirtar_gene_summary.csv
├── mirtar_mirna_summary.csv
└── mirtar_family_summary.csv
```

---

## Normalization & Filtering

Before any aggregation the raw table passes through several normalization steps:

- **Support types** are mapped to four canonical strings: `Functional MTI`, `Functional MTI (Weak)`, `Non-Functional MTI`, `Non-Functional MTI (Weak)`. Rows with unrecognized support types are dropped.
- **miRNA names** have whitespace stripped (e.g. `hsa-miR- 21-5p` → `hsa-miR-21-5p`).
- **Gene symbols** are whitespace-trimmed.
- **Experiments** are split on `//`, `;`, `,`, `|` delimiters and classified into six broad classes: `reporter`, `binding`, `protein`, `rna`, `perturbation`, `other`.
- **Family assignment** uses TargetScan seed-family first (`miR_Family_Info.txt`, human only). For miRNAs absent from TargetScan the fallback strips the species prefix and arm suffix (e.g. `hsa-miR-155-5p` → `miR-155`).
- **Study IDs** are extracted from the references/PMID column as 5–9 digit numbers. When no numeric PMID is found, the raw reference string is used as a pseudo-ID (`REF::<text>`). If no references column exists at all, the miRTarBase ID is used (`MIRTAR::<id>`).
- **Gene panel filter**: Only rows matching `PRIMARY_GENES` (or a custom list) are retained.

---

## Table 1: `mirtar_interaction_study_collapsed.csv`

**Granularity:** One row per (miRNA, gene, study).

This is the foundational table from which all summaries are derived. Each row represents what a single publication reported about a specific miRNA–gene interaction.

### Columns

| Column | Type | Description |
|--------|------|-------------|
| `miRNA` | str | Normalized miRNA name (e.g. `hsa-miR-155-5p`) |
| `miRNA_family` | str | TargetScan family or name-based proxy (e.g. `miR-155`) |
| `gene` | str | Normalized gene symbol |
| `study_id` | str | PMID or pseudo-study identifier |
| `has_functional_mti` | 0/1 | Study reported Functional MTI |
| `has_functional_mti_weak` | 0/1 | Study reported Functional MTI (Weak) |
| `has_nonfunctional_mti` | 0/1 | Study reported Non-Functional MTI |
| `has_nonfunctional_mti_weak` | 0/1 | Study reported Non-Functional MTI (Weak) |
| `has_reporter` | 0/1 | Study used reporter assay (luciferase, GFP) |
| `has_binding` | 0/1 | Study used binding assay (CLIP, RIP, pull-down) |
| `has_protein` | 0/1 | Study used protein-level readout (western blot, ELISA, flow cytometry) |
| `has_rna` | 0/1 | Study used RNA-level readout (qRT-PCR, RNA-seq, microarray) |
| `has_perturbation` | 0/1 | Study used perturbation (overexpression, knockdown, mimic, inhibitor) |
| `has_other` | 0/1 | Study used experiment not matching above classes |
| `has_{experiment}__{support_slug}` | 0/1 | Cross flag: experiment class AND support type co-occurred in this study. 24 columns total (6 experiments × 4 support types). Example: `has_reporter__functional_mti` |

---

## Table 2: `mirtar_interaction_summary.csv`

**Granularity:** One row per (miRNA, gene) pair.

Aggregates across studies. The primary table for asking "how well-supported is the interaction between miR-X and gene Y?"

### Columns

| Column | Type | Description |
|--------|------|-------------|
| `miRNA` | str | Normalized miRNA name |
| `miRNA_family` | str | Family label |
| `gene` | str | Gene symbol |
| `n_studies` | int | Total unique studies reporting this interaction |
| `n_functional_mti_studies` | int | Studies with Functional MTI |
| `n_functional_mti_weak_studies` | int | Studies with Functional MTI (Weak) |
| `n_nonfunctional_mti_studies` | int | Studies with Non-Functional MTI |
| `n_nonfunctional_mti_weak_studies` | int | Studies with Non-Functional MTI (Weak) |
| `n_reporter_studies` | int | Studies using reporter assays |
| `n_binding_studies` | int | Studies using binding assays |
| `n_protein_studies` | int | Studies with protein readouts |
| `n_rna_studies` | int | Studies with RNA readouts |
| `n_perturbation_studies` | int | Studies using perturbation |
| `n_other_studies` | int | Studies with other experiment types |
| `n_{experiment}__{support_slug}_studies` | int | Cross-count columns (24 total). Example: `n_reporter__functional_mti_studies` = number of studies where reporter assay was used with Functional MTI support |
| `experiment_support_counts_json` | JSON str | Nested dict: `{experiment: {support_type: count}}` — same data as the cross-count columns in a single queryable field |
| `evidence_score` | int | Weighted composite score: `3×reporter + 3×binding + 2×protein + 1×rna + 1×perturbation` (summed across studies) |

---

## Table 3: `mirtar_gene_summary.csv`

**Granularity:** One row per gene.

The main gene-centric table. For each APM gene, summarizes all miRNAs that target it and the strength of evidence.

### Scalar Columns

| Column | Type | Description |
|--------|------|-------------|
| `gene` | str | Gene symbol |
| `n_unique_miRNAs` | int | Distinct miRNAs targeting this gene |
| `n_unique_families` | int | Distinct miRNA families |
| `n_total_studies` | int | Total unique studies across all miRNAs |
| `n_{support_slug}_studies` | int | Study counts per support type (4 columns) |
| `n_{experiment}_studies` | int | Study counts per experiment class (6 columns) |
| `n_{experiment}__{support_slug}_studies` | int | Cross-counts (24 columns) |

### JSON Dict Columns

These columns store JSON-serialized dicts for flexible downstream querying. Use `json.loads()` or `ast.literal_eval()` to recover.

| Column | Structure | Description |
|--------|-----------|-------------|
| `experiment_support_counts_json` | `{experiment: {support_type: int}}` | Cross-tabulated study counts |
| `mirna_study_counts_json` | `{miRNA: int}` | Per-miRNA study counts, sorted descending |
| `family_study_counts_json` | `{family: int}` | Per-family study counts, sorted descending |
| `support_to_mirnas_json` | `{support_type: [miRNA, ...]}` | Which miRNAs have each support type |
| `experiment_to_mirnas_json` | `{experiment: [miRNA, ...]}` | Which miRNAs have each experiment class |
| `experiment_support_to_mirnas_json` | `{experiment: {support_type: [miRNA, ...]}}` | Full cross-membership |
| `support_to_families_json` | `{support_type: [family, ...]}` | Family-level membership by support type |
| `experiment_to_families_json` | `{experiment: [family, ...]}` | Family-level membership by experiment |
| `experiment_support_to_families_json` | `{experiment: {support_type: [family, ...]}}` | Full cross-membership at family level |
| `mirna_detail_json` | `{miRNA: {n_studies, support_counts, experiment_support_counts}}` | Full per-miRNA breakdown (see structure below) |

**`mirna_detail_json` entry structure:**

```json
{
  "hsa-miR-155-5p": {
    "n_studies": 4,
    "support_counts": {
      "Functional MTI": 3,
      "Functional MTI (Weak)": 1,
      "Non-Functional MTI": 0,
      "Non-Functional MTI (Weak)": 0
    },
    "experiment_support_counts": {
      "reporter": {"Functional MTI": 2, "Functional MTI (Weak)": 0, ...},
      "binding":  {"Functional MTI": 1, ...},
      ...
    }
  }
}
```

---

## Table 4: `mirtar_mirna_summary.csv`

**Granularity:** One row per miRNA.

The miRNA-centric mirror of the gene summary. For each miRNA, summarizes which panel genes it targets and with what evidence.

### Scalar Columns

| Column | Type | Description |
|--------|------|-------------|
| `miRNA` | str | Normalized miRNA name |
| `miRNA_family` | str | Family label |
| `n_unique_targets` | int | Distinct genes targeted (within gene panel) |
| `n_total_studies` | int | Total unique studies |
| `n_{support_slug}_studies` | int | Study counts per support type (4 columns) |
| `n_{experiment}_studies` | int | Study counts per experiment class (6 columns) |
| `n_{experiment}__{support_slug}_studies` | int | Cross-counts (24 columns) |

### JSON Dict Columns

| Column | Structure | Description |
|--------|-----------|-------------|
| `experiment_support_counts_json` | `{experiment: {support_type: int}}` | Cross-tabulated study counts |
| `target_gene_study_counts_json` | `{gene: int}` | Per-gene study counts, sorted descending |
| `support_to_genes_json` | `{support_type: [gene, ...]}` | Which genes have each support type |
| `experiment_to_genes_json` | `{experiment: [gene, ...]}` | Which genes have each experiment class |
| `experiment_support_to_genes_json` | `{experiment: {support_type: [gene, ...]}}` | Full cross-membership |
| `target_gene_detail_json` | `{gene: {n_studies, support_counts, experiment_support_counts}}` | Full per-gene breakdown (same structure as `mirna_detail_json` in gene_summary) |

---

## Table 5: `mirtar_family_summary.csv`

**Granularity:** One row per miRNA family.

Lightweight summary for family-level analysis.

### Columns

| Column | Type | Description |
|--------|------|-------------|
| `miRNA_family` | str | Family label (TargetScan or name-proxy) |
| `n_unique_miRNAs` | int | Distinct miRNA members |
| `n_unique_targets` | int | Distinct gene targets |
| `n_total_studies` | int | Total unique studies |
| `n_{support_slug}_studies` | int | Study counts per support type (4 columns) |
| `n_{experiment}_studies` | int | Study counts per experiment class (6 columns) |

---

## Support Type Slugs Reference

The support type strings are converted to column-safe slugs used in all `has_*` and `n_*` column names:

| Support Type | Slug |
|-------------|------|
| `Functional MTI` | `functional_mti` |
| `Functional MTI (Weak)` | `functional_mti_weak` |
| `Non-Functional MTI` | `nonfunctional_mti` |
| `Non-Functional MTI (Weak)` | `nonfunctional_mti_weak` |

---

## Experiment Classification Reference

Raw miRTarBase experiment names are classified into six broad classes by keyword matching (evaluated in priority order):

| Class | Keywords (substring match, case-insensitive) |
|-------|----------------------------------------------|
| `reporter` | luciferase, reporter assay, reporter, gfp reporter |
| `binding` | clip, hits-clip, par-clip, clash, rip, rip-chip, pull-down, pulldown, immunoprecipitation, binding assay |
| `protein` | western blot, elisa, flow cytometry, immunoblot, protein assay, immunohistochemistry, ihc |
| `rna` | qrt-pcr, qpcr, rt-pcr, northern blot, microarray, rna-seq, sequencing, expression profiling, pcr |
| `perturbation` | overexpression, knockdown, transfection, sirna, shrna, mimic, inhibitor, antisense, gene silencing |
| `other` | anything not matching above |

---

## Evidence Score

The `evidence_score` column in `mirtar_interaction_summary.csv` is a simple weighted sum across studies for a given (miRNA, gene) pair:

```
evidence_score = 3 × Σ(reporter studies)
               + 3 × Σ(binding studies)
               + 2 × Σ(protein studies)
               + 1 × Σ(rna studies)
               + 1 × Σ(perturbation studies)
```

Reporter and binding assays receive the highest weight because they provide direct evidence of miRNA–target interaction (3'UTR binding + functional repression). Protein readouts get intermediate weight (confirms downstream effect). RNA-level and perturbation experiments get baseline weight (correlative or indirect).

---

## Usage

### From `main.py` (pipeline integration)

Called as step 9b inside `run_full_pipeline()` when `skip_mirna=False`:

```python
from .genes.mirtarbase import get_mirtarbase_targets

mirtarbase_tables = get_mirtarbase_targets(
    mirtarbase_csv=PATHS.mirtarbase_csv,
    family_info_tsv=PATHS.mir_family_info,
    gene_panel=gene_panel,
    output_dir=working_dir / "miRNA" / "mirtarbase",
)
```

### From Jupyter

```python
import sys
sys.path.insert(0, '/home/stavz/masters/gdc')

from pipeline.genes.mirtarbase import get_mirtarbase_targets
from pipeline.config import PATHS, PRIMARY_GENES

tables = get_mirtarbase_targets()

# Access individual tables
interaction_summary = tables["interaction_summary"]
gene_summary = tables["gene_summary"]

# Quick look: top miRNAs targeting TAP1 by evidence score
tap1 = interaction_summary[interaction_summary["gene"] == "TAP1"]
tap1.sort_values("evidence_score", ascending=False).head(10)
```

### Step-by-step (finer control)

```python
from pipeline.genes.mirtarbase import (
    load_mirtarbase,
    build_interaction_study_table,
    build_interaction_summary,
    build_gene_summary,
    build_mirna_summary,
    build_family_summary,
)

# Load with custom gene panel
df = load_mirtarbase(gene_panel=["HLA-A", "HLA-B", "TAP1", "B2M"])

# Build only the tables you need
interaction_study = build_interaction_study_table(df)
gene_summary = build_gene_summary(interaction_study)
```

### Recovering JSON columns

```python
import json

row = gene_summary[gene_summary["gene"] == "TAP1"].iloc[0]

# Per-miRNA study counts
mirna_counts = json.loads(row["mirna_study_counts_json"])
# → {"hsa-miR-155-5p": 4, "hsa-miR-21-5p": 2, ...}

# Full miRNA detail
detail = json.loads(row["mirna_detail_json"])
# → {"hsa-miR-155-5p": {"n_studies": 4, "support_counts": {...}, ...}}

# Membership: which miRNAs have reporter + Functional MTI?
memberships = json.loads(row["experiment_support_to_mirnas_json"])
reporter_functional = memberships["reporter"]["Functional MTI"]
# → ["hsa-miR-155-5p", "hsa-miR-21-5p"]
```

---

## Config Requirements

Two paths must be added to `PathConfig` in `config.py`:

```python
@dataclass
class PathConfig:
    # ... existing paths ...

    # miRTarBase
    mirtarbase_csv: Path = Path("/home/stavz/masters/gdc/APM/data/miRNA/mirtar.csv")
    mir_family_info: Path = Path("/home/stavz/masters/gdc/APM/data/miRNA/miR_Family_Info.txt")
```

No new thresholds are required. The module uses the existing `PRIMARY_GENES` and `OUTPUT_SUBDIRS["mirna"]` from config.
