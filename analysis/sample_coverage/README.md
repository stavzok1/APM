## Sample coverage analyses

This folder contains **downstream analysis** scripts and their outputs for cohort-level
sample/module coverage summaries.

### Scripts

- `analysis/sample_module_coverage.py`: computes per-module sample counts, intersections, pairwise overlap, and writes figures.
- `analysis/clinical_omics_stratification.py`: stratifies omics coverage by PAM50 and collapsed stage using the unified clinical table.

### Output layout

Default run directory:

- `analysis/sample_coverage/output/run_<UTC timestamp>/`

Within each run:

- `tables/`
  - `omics/`: vial/participant omics tables (presence, counts, intersections, pairwise overlap)
  - `metadata/`: participant metadata-only tables (presence, counts, intersections, pairwise overlap)
  - `participant_presence.tsv`: full participant matrix (omics + metadata) for convenience
  - `module_counts.tsv`: combined legacy counts table for convenience
- `figures/`
  - `omics/`: figures for vial + participant omics bundles
  - `metadata/`: figures for participant metadata bundle
- `notes/`: non–patient-local reference module notes

