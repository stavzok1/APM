## TAD processing (original `/home/stavz/masters/gdc/TADs`)

This repository contains a standalone “upstream” TAD-processing codebase under:

- `/home/stavz/masters/gdc/TADs/`

It ingests multiple study-specific input formats (Arrowhead, `.hic`, sparse matrices, `.mcool`, dense per-chrom matrices), normalizes them into a common representation, and writes standardized per-biosample outputs under:

- `/home/stavz/masters/gdc/TADs/processed/<biosample>/…`

Downstream (in this repo), `pipeline/tad_annotation/` consumes the **standardized** domain/boundary/flank tables and annotates pipeline features (genes, cCREs, lncRNAs, SVs) with TAD context.

### Common outputs (per biosample)

Across datasets, the processing code aims to produce (at minimum):

- **Domains** (intervals): TSV + BED
- **Boundaries** (intervals centered on insulation minima / inferred breakpoints): TSV + BED
- **Domain flanks**: TSV (left/right boundary metadata per domain)
- (Sometimes) intermediate contact containers (`.cool`) and insulation tracks (`*.insulation.tsv`)

Naming varies by dataset; examples seen in `processed/`:

- `Rao_HMEC_domains_hg38.tsv`, `Rao_HMEC_boundaries_hg38.tsv`, `Rao_HMEC_flanks_hg19.tsv`
- `Kim_T47D.domains_w200000.hg38.tsv`, `Kim_T47D.boundaries_w200000.hg38.tsv`, `Kim_T47D.domain_flanks_w200000.hg19.tsv`
- `vandenBrand_HB2.domains_w200000.hg38.tsv`, `vandenBrand_HB2.boundaries_w200000.hg38.tsv`, `vandenBrand_HB2.r20000.insulation.tsv`

### Dataset categories (pipelines)

The canonical implementations live in `/home/stavz/masters/gdc/TADs/datasets.py`, wired via `/home/stavz/masters/gdc/TADs/cli.py`.

#### Category A — Rao (GSE63525) Arrowhead domainlist (HMEC; hg19 → hg38)

**Entry point:** `process_rao_hmec(...)`

- Input: gzipped Arrowhead “domainlist” TSV (hg19)
- Convert to domains (BED-like)
- Derive boundaries from domain edges (with flank + merge)
- Write hg19 outputs
- Liftover domains + boundaries to hg38 (UCSC liftOver or pyliftover)
- Write hg38 outputs
- Also writes “unmapped” lists and a post-liftover “gaps” report

#### Category B — Kim (GSE167150) `.hic` (multiple cell lines; hg19 → hg38)

**Entry point:** `process_kim_hic(...)` and `process_kim_batch(...)`

Pipeline (per biosample):

- `.hic` → `.cool` (via `hic2cool convert`)
- Balance `.cool` (`cooler balance`)
- Insulation track (`cooltools insulation`) for configured windows (default includes 200kb)
- Convert insulation bins → **boundaries** (optionally “strong only”)
- Convert boundaries → **domains** (+ flanks) using chromosome sizes
- Write hg19 outputs
- Liftover domains + boundaries to hg38; write hg38 outputs

Outputs are named with:
- resolution: `r20000` (20kb bins)
- insulation window: `w200000` (200kb)
- genome build: `hg19` / `hg38`

**Default parameters (as implemented in `TADs/config.py` and CLI defaults):**
- **Resolution**: 20kb (`--resolution 20000`)
- **Insulation windows**: 200kb, 400kb (but **boundary calling** selects `pick_window_bp`, default 200kb)
- **Boundary region half-width**: 20kb (`flank_bp`)
- **Boundary merge distance**: 20kb (`merge_within_bp`)
- **Boundary selection**: `strong_only=True` (keeps `is_boundary_<window>` bins from cooltools)
- **Liftover**: hg19→hg38 using UCSC `liftOver` unless disabled

#### Category C — Le Dily (GSE109229) sparse contact matrix (hg38)

**Entry point:** `process_ledily_matrix(...)`

- Input: sparse triplet matrix (chrom, pos1, pos2, count), already hg38
- Build a `.cool` file via `cooler.create_cooler`
- Balance + insulation + boundary calling (same as Kim)
- Derive domains from boundaries + write outputs (hg38 only)

**Default parameters:** same insulation/boundary defaults as Category B (Kim), except the matrix is fixed at 100kb binning.

#### Category D — Golloshi (GSE143678) dense per-chrom matrices (hg19 → hg38)

**Entry point:** `process_golloshi_matrices(...)`

- Input: per-chrom dense matrices (`*.matrix.gz`) with ICE-corrected floats and bin labels encoding hg19 coordinates
- Parse matrices → global bins + upper-triangle pixels → `.cool` (hg19)
- Balance + insulation + boundary calling
- Derive domains from boundaries
- Write hg19 outputs
- Liftover domains + boundaries to hg38; write hg38 outputs + “unmapped” reports

**Default parameters:** similar to Category B (Kim), but resolution varies by the chosen matrix directory (`--resolution-kb`).

#### Category E — van den Brand (GSE273999) `.mcool` (hg38)

**Entry point:** `process_vandenbrand_mcool(...)`

- Input: multi-resolution cooler (`.mcool`)
- Select a resolution (default 20kb)
- Insulation + boundary calling + domains from boundaries
- Write hg38 outputs

**Default parameters:** same insulation/boundary defaults as Category B (Kim), but resolution is selected from the `.mcool` URI.

### Notes for downstream mirroring

The *processing* layer produces boundaries/domains as interval tables only. Any “boundary layers”
(e.g. overlaps with CTCF peaks, ATAC peaks, methylation probes, SV breakpoints, CNV deletions)
should be added in a **separate annotation/mirroring stage**, ideally in this repo’s
`pipeline/tad_annotation/` so it shares identifier normalization and output conventions.

