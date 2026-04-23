# SV pipeline (VEP, mapping, FIMO, neojunction, ChIP)

**Canonical SV documentation for APM** — narrative and operational detail for structural-variant processing live here only (not duplicated under `pipeline/SV/*.md`). Command cheat sheet: [`analysis/COMMANDS.md`](../../../analysis/COMMANDS.md).

**Code layout:** Python package [`pipeline/SV/`](../../SV/) (`pipeline.py`, `vcf_loader.py`, `sv_filtering.py`, `spatial_mapping.py`, `vep_annotation.py`, `bed_intervals.py`, `motif_scanning.py`). Shell helpers are checked in under [`pipeline/SV/shell/`](../../SV/shell/) and under [`scripts/SV_pipeline_scripts/`](../../../scripts/SV_pipeline_scripts/) (JSON `SCRIPTS_DIR` often points at `SV_scripts` there).

---

## How runs are orchestrated today

### Path A — Python (`pipeline.SV.pipeline`)

**`run_sv_pipeline(...)`** (see `pipeline/SV/pipeline.py`) is the integrated entry used by cohort drivers and smoke tools.

1. **VCF → strict/lenient tables (unless `resume_motifs_from_fimo`)**  
   - Optionally runs **VEP** on each raw VCF into `01_vep_vcfs/` (unless `skip_vep`).  
   - Then **`run_sv_vcf_processing`** reads VCFs from the active directory (raw or VEP-out) and writes `02_processed_sv_csv/<sample>_<filtering>_sv_set.csv`.  
   - Per-row work includes Manta parsing, **strict/lenient** filters, **gene / cCRE / lncRNA** windows, optional **miRNA mature-arm** proximity when a miRNA loci table is configured, BND alt harmonization, **VEP CSQ → nested columns** when CSQ exists (otherwise empty VEP-shaped columns so downstream schemas stay stable), and TCGA id harmonization.

2. **Motifs + downstream (unless `skip_motifs`) — `run_sv_motif_scanning(...)`**  
   Internal numbered logs inside that function:
   - **[STEP 1]** BED intervals → `03_sv_bed/`  
   - **[STEP 2]** reference FASTA → `04_sv_fasta/`  
   - **[STEP 3]** FIMO → `05_fimo_tsv/`  
   - **[STEP 4]** `bedtools intersect` → `06_sv_fimo_merged/`  
   - **[STEP 5]** Python **`recombine_all_sv_fimo`** → `07_final_sv_with_fimo/`  
   - **[STEP 6]** **BND neojunction** motifs (genome-wide cCRE FIMO lookup) → `08_neojunction_enriched/` when `PATHS.all_ccre_fimo_tsv` and cCRE table paths resolve.  
   - **[STEP 7]** **ChIP** peak-disruption annotation → `09_chip_enriched/` when ChIP inputs resolve and `skip_chip` is false.

**Deliverable preference:** if `09_chip_enriched` exists it is treated as the SV “final” CSV set; else `08_neojunction_enriched`; else `07_final_sv_with_fimo`.

**`resume_motifs_from_fimo=True`:** skips the VCF/VEP leg, reuses `02_processed_sv_csv`, and runs motif **steps 4–7** only using existing `03_sv_bed` / `05_fimo_tsv` (no new FASTA/FIMO in steps 1–3).

### Path B — Shell (`run_sv_motif_pipeline.sh`)

The bash driver ([`scripts/SV_pipeline_scripts/run_sv_motif_pipeline.sh`](../../../scripts/SV_pipeline_scripts/run_sv_motif_pipeline.sh) or the copy under `pipeline/SV/shell/`) stages: `vep` → `process` → `bed` → `fasta` → `fimo` → `intersect` → `recombine`, producing **`01`–`07` only**. It does **not** execute Python **STEP 6–7** (neojunction, ChIP). Use **`run_sv_motif_scanning`** / **`run_sv_pipeline`** for those.

---

## Overview (shell-aligned stages)

1. **VEP annotate** raw VCFs → annotated VCFs  
2. **Process** annotated VCFs → strict SV CSV tables  
3. **Create BED** SV intervals (with configurable flank)  
4. **Extract FASTA** sequences for BED intervals  
5. **Run FIMO** motif scanning on FASTA sequences  
6. **Intersect** FIMO hits with SV BED intervals  
7. **Recombine** motif hits back into processed SV CSV tables  

Steps **8–9** (neojunction, ChIP) are **Python-only** (see above).

---

## Configuration (`pipeline/config.py`)

### TF panel for MEME extraction (`SV_TARGET_TF_SYMBOLS`)

`pipeline/SV/motif_scanning.extract_selected_motifs()` defaults to the same list as **`SV_TARGET_TF_SYMBOLS`** in `pipeline/config.py`: case-insensitive **substring** match on each full HOCOMOCO `MOTIF` name (e.g. `STAT3.H14CORE.0.P.B`). The list is aligned with **HOCOMOCO v14 H14*** bundles for IFN / NF-κB / AP-1 / CTCF / MYC / **IRF9** / **NF-Y (NFYA–C)** / **RFX5** / **FOXA1, ESR1, GATA3** / **SMAD***.

Factors **not** named in standard H14 human bundles (extend the MEME library from JASPAR/CIS-BP/etc. if you need them in FIMO): **NLRC5**, **CIITA**, **EZH2**, **SUZ12**, **BRD4**, **RELA** (p65).

Multi-bundle name checks: [`scripts/motifs/hocomoco_multi_bundle_tf_report.py`](../../../scripts/motifs/hocomoco_multi_bundle_tf_report.py) (scan all `H14*_meme_format.meme` in a directory).

Motif **intersection** still uses FIMO’s `--thresh` via `THRESHOLDS.fimo_pvalue_threshold` when FIMO is invoked from Python. The **recombine** step (Python `recombine_sv_fimo` in `pipeline/SV/motif_scanning.py`) uses `ThresholdConfig` fields—**not environment variables**:

| Field | Role |
|------|------|
| `sv_fimo_recombine_chunk_rows` | Rows per chunk when reading each `06_sv_fimo_merged/*.bed` |
| `sv_fimo_max_flank_hits_per_tf` | Max flank motif hits kept per `(SV id, TF)` (best p-value within that TF) |
| `sv_fimo_max_elem_hits_per_tf` | Max motif hits kept per `(SV id, regulatory elem_id, TF)` |
| `sv_fimo_recombine_stream_if_merged_bed_bytes` | If the merged BED is at least this many bytes, STEP 5 reads/writes the SV CSV in chunks (avoids holding the full serialized CSV in RAM). `0` disables streaming. |
| `sv_fimo_recombine_sv_csv_chunk_rows` | SV rows per pandas chunk when streaming STEP 5 output |
| `sv_fimo_recombine_skip_csv_basenames` | Stems of `02_processed_sv_csv/*.csv` to skip for FASTA/FIMO/intersect/recombine; `07` gets a copy of that CSV without motif attachment. Empty list = skip nothing. |

Defaults live on `THRESHOLDS` in `pipeline/config.py`. Increase if you have RAM and want more hits stored; decrease chunk sizes or enable streaming at a lower byte threshold if the process is OOM-killed. A sample whose `06` merged BED is orders of magnitude larger than others usually reflects many more SV/flank/element intervals intersecting many FIMO hits upstream, not recombine alone.

### Curated MEME file (`PATHS.sv_meme_file`)

The default file (`data/SV/motifs/core_ifn_stat_ctcf_ap1_nlrc5.meme`; legacy basename) is a **subset** of motifs whose names match **`SV_TARGET_TF_SYMBOLS`**, extracted from a **HOCOMOCO v14** dump (typically **H14CORE**; each `MOTIF` line looks like `STAT3.H14CORE.0.P.B`). **SNV FIMO** reuses this same path.

**TF / factor tokens in `SV_TARGET_TF_SYMBOLS` (2026-04 refresh)** — substring match on full HOCOMOCO `MOTIF` names when running `extract_selected_motifs()`:

- **STAT / IRF:** `STAT1`, `STAT2`, `STAT3`, `IRF1`, `IRF2`, `IRF3`, `IRF9`
- **NF-κB:** `NFKB1`, `NFKB2`, `RELB`, `REL` (there is no separate `RELA` / `p65` motif name in H14 bundles; use JASPAR/CIS-BP if you need that label explicitly)
- **AP-1 / bZIP:** `FOS`, `FOSL1`, `FOSL2`, `JUN`, `JUNB`, `JUND`, `ATF3`
- **Insulator / MYC:** `CTCF`, `CTCFL`, `MYC`
- **NF-Y:** `NFYA`, `NFYB`, `NFYC`
- **Breast / lineage (H14):** `RFX5`, `FOXA1`, `ESR1`, `GATA3`
- **TGF-β:** `SMAD` (matches `SMAD1`, `SMAD2`, … in HOCOMOCO names)

**Not present** under standard **H14CORE / H14INVIVO / H14INVITRO / H14RSNP** human motif names (add from JASPAR or other sources if required): `NLRC5`, `CIITA`, `EZH2`, `SUZ12`, `BRD4`, `RELA` (gene product p65).

To rebuild after downloading a new master MEME database, use `extract_selected_motifs()` in `pipeline/SV/motif_scanning.py` (see **Curated SV/SNV MEME library** in `analysis/COMMANDS.md`). For a **composite** library, run extraction per source into separate MEME files, then **concatenate motif blocks** under one shared MEME v4 header (one `MEME version` / `ALPHABET` / background block; append all `MOTIF …` sections with **globally unique motif names** so FIMO can distinguish them).

### FIMO recombine (STEP 5) and memory

`bedtools intersect -wa -wb` between SV BED lines and FIMO hits can produce **many millions** of rows (every overlap × motif). The recombine step **streams** the merged BED in chunks of `THRESHOLDS.sv_fimo_recombine_chunk_rows` and keeps only the **best p-value hits per TF** separately for flanks (`sv_fimo_max_flank_hits_per_tf`) and for each regulatory element (`sv_fimo_max_elem_hits_per_tf`). TF here is the short prefix before the first `.` in FIMO’s `motif_id` (same as the `"TF"` field written into each hit dict).

When the merged `06_sv_fimo_merged` BED is larger than `sv_fimo_recombine_stream_if_merged_bed_bytes`, STEP 5 also **streams the SV CSV** in `sv_fimo_recombine_sv_csv_chunk_rows` chunks so the full table with serialized `elem_hits` / `flank_motif_hits` is never held in memory (a common cause of OOM spikes even when average RAM looks fine). Set `sv_fimo_recombine_stream_if_merged_bed_bytes` to `0` to force the legacy single-load path. Tune these thresholds in `pipeline/config.py` (`ThresholdConfig`) if jobs are killed for memory or if you need denser motif lists.

`sv_fimo_recombine_skip_csv_basenames` lists processed SV **file stems** (e.g. `TCGA-…_strict_sv_set`) to omit from FASTA extraction, FIMO, intersect, and recombine; `07_final` still receives a **copy** of that rowset from `02_processed_sv_csv` (no motif columns). Clear the list when you want those samples processed again.

### Duplicate-looking rows in `05_fimo_tsv` and `06_sv_fimo_merged`

**Not duplicates:** FIMO’s `fimo.tsv` normally has **many rows per locus** because each row is one **(sequence window, motif, strand, …)** hit. Different `motif_id` values at the same coordinates are expected.

**True duplicate lines (identical across all columns in the merged BED):** `bedtools intersect -wa -wb` emits one output row per **pair** of overlapping intervals from `-a` (SV flank / element-overlap BED) and `-b` (FIMO hits BED). Identical pairs only happen if the **same** line appears more than once in `-a` or `-b`. Typical causes:

- **`03_sv_bed`:** duplicate entries in `elem_hits` for the same SV (same `elem_id`, `overlap_start`, `overlap_end`) produced two identical BED lines. The builder now **drops duplicate BED rows** before writing.
- **`05` → FIMO hits BED:** rare duplicate rows in the parsed FIMO table; these are **deduped** when writing `*_fimo_hits.bed`.

**Same FIMO hit repeated with different BED `name` values:** that is **not** a duplicate row. Each overlapping **cCRE** becomes its own BED line in `build_sv_flanks_and_overlaps_bed`: coordinates may match another cCRE’s overlap, but the 4th column (`name`) includes `elem:{elem_id}`, so the full line differs. `bedtools` then reports one intersection row per cCRE interval; STEP 5 routes motifs into `elem_hits[].motif_hits` **by `elem_id`**, so each element gets its own motif list.

**Genome-wide cCRE FIMO (`all_ccre_fimo.tsv`, STEP 6 / neojunction):** `load_all_ccre_fimo` assigns each genomic motif hit to **every cCRE** whose coordinates it overlaps (see docstring: a hit may fall in zero or more cCREs). The lookup is **`cCRE_id → [motif dicts]`**; the same motif interval can appear in **multiple** cCRE buckets when several cCREs overlap that genomic span—by design, not an error.

### Output directories (`PATHS.output_subdirs` pattern)

Under `sv_output_root` / working tree:

- `01_vep_vcfs` — VEP-annotated VCFs (Python path; shell equivalent)  
- `02_processed_sv_csv` — processed SV tables  
- `03_sv_bed` — BED intervals  
- `04_sv_fasta` — FASTA per SV BED  
- `05_fimo_tsv` — FIMO outputs  
- `06_sv_fimo_merged` — intersected motif evidence  
- `07_final_sv_with_fimo` — SV CSVs with FIMO merged into rows  
- `08_neojunction_enriched` — post–STEP 5, BND neojunction motifs (Python)  
- `09_chip_enriched` — ChIP disruption columns (Python)  
- `logs` — where shell driver writes logs; Python steps log to stdout  

---

## Shell scripts (bash driver + per-step)

### 1) `run_sv_motif_pipeline.sh` (master)

**Purpose**  
Orchestrate the end-to-end SV→motif workflow with stage selection, resuming, and logging through **`07`** only.

**Stages**

- `vep`
- `process`
- `bed`
- `fasta`
- `fimo`
- `intersect`
- `recombine`

**Arguments (JSON mode — recommended)**

    ./run_sv_motif_pipeline.sh --config sv_paths.json [--threads N] [--force] [--from STEP] [--to STEP]

**Arguments (positional mode — backward compatible)**

    ./run_sv_motif_pipeline.sh [--threads N] [--force] [--from STEP] [--to STEP] \
        SCRIPTS_DIR RAW_VCF_DIR OUTPUT_ROOT GENE_FEATURES SAMPLES_TSV \
        REG_ELEM_FILE LNCRNA_FEATURES REF_FASTA MEME_FILE FLANK_SIZE

**Inputs**

- Raw VCF directory (`RAW_VCF_DIR`)
- Feature tables:
  - `GENE_FEATURES` (CSV)
  - `LNCRNA_FEATURES` (CSV)
  - `REG_ELEM_FILE` (CSV)
  - `SAMPLES_TSV` (TSV)
- Reference FASTA (`REF_FASTA`)
- MEME motif file (`MEME_FILE`)

**Outputs** (under `OUTPUT_ROOT`): `01`–`07` and `logs/` as above; `.done` markers for resume.

---

### 2) `run_vep_batch.sh`

**Purpose**  
Batch-run VEP on all `.vcf` files in a directory.

**Arguments**

    ./run_vep_batch.sh [--fork N] [--dir-cache PATH] [--fasta PATH] [--force] <input_vcf_dir> <output_dir>

**Outputs**

- `<output_dir>/*.vcf` — VEP-annotated VCFs (same basenames)

---

### 3) `SV_vcf_handling.py`

**Purpose**  
Convert VEP-annotated VCFs to strict SV tables and enrich them with gene/regulatory/lncRNA context (legacy CLI; **`run_sv_vcf_processing`** in Python is the maintained batch path).

**Arguments**

    python SV_vcf_handling.py <vep_vcf_dir> <gene_features.csv> <samples.tsv> <out_dir> <reg_elem.csv> <lncrna_features.csv>

---

### 4) `SV_create_bed.py`

**Purpose**  
Convert processed SV CSVs to BED intervals suitable for sequence extraction, with configurable flank.

**Arguments**

    python SV_create_bed.py <processed_sv_csv_dir> <out_bed_dir> <flank_size>

---

### 5) `create_fasta_for_sv_bed.sh`

**Purpose**  
Extract reference sequences for SV BED intervals from the reference FASTA.

**Arguments**

    ./create_fasta_for_sv_bed.sh <ref_fasta> <sv_bed_dir> <out_fasta_dir>

---

### 6) `run_fimo_on_fasta.sh`

**Purpose**  
Run motif scanning (FIMO) on FASTA sequences and collect `fimo.tsv` per input FASTA.

**Arguments**

    ./run_fimo_on_fasta.sh [--jobs N] [--thresh X] [--force] <meme_file> <fasta_dir> <out_dir>

---

### 7) `run_fimo_sv_intersect.sh`

**Purpose**  
Intersect motif hits (FIMO) back onto SV BED intervals.

**Arguments**

    ./run_fimo_sv_intersect.sh <fimo_tsv_dir> <sv_bed_dir> <out_dir>

---

### 8) `recombine_sv_fimo.py`

**Purpose**  
Legacy CLI wrapper; the maintained implementation is **`recombine_sv_fimo` / `recombine_all_sv_fimo`** in `pipeline/SV/motif_scanning.py` (STEP 5 above).

---

## Configuration: `sv_paths.json`

Example:

```json
{
  "SCRIPTS_DIR": "SV_scripts",
  "RAW_VCF_DIR": "SV_vcfs_fimo",
  "OUTPUT_ROOT": "SV_pipeline_out",

  "GENE_FEATURES": "SV_fimo_inputs/primary_genes_all_features.csv",
  "SAMPLES_TSV": "SV_fimo_inputs/samples.tsv",
  "REG_ELEM_FILE": "SV_fimo_inputs/regulatory_element_focus.csv",
  "LNCRNA_FEATURES": "SV_fimo_inputs/primary_lncRNAs_all_features.csv",

  "REF_FASTA": "/home/you/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  "MEME_FILE": "SV_fimo_inputs/core_ifn_stat_ctcf_ap1_nlrc5.meme",

  "FLANK_SIZE": 150,
  "THREADS": 8
}
```

---

## Typical usage

**Shell — full pipeline to `07`:**

```bash
./run_sv_motif_pipeline.sh --config sv_paths.json
```

Resume from a specific step (e.g., run only FASTA → recombine):

```bash
./run_sv_motif_pipeline.sh --config sv_paths.json --from fasta --to recombine
```

**Python — full stack including neojunction + ChIP:**

```python
from pathlib import Path
from pipeline.SV.pipeline import run_sv_pipeline

results = run_sv_pipeline(
    vcf_dir=Path("raw_vcfs/"),
    output_root=Path("sv_output/"),
    filtering="strict",
    skip_vep=False,
    skip_motifs=False,
    skip_chip=False,
    resume_motifs_from_fimo=False,
)
# Prefer results["final_csvs"] under 09 / 08 / 07 as documented above.
```

---

## Output columns (processed / final CSVs)

### Processed SV CSV columns (representative)

- Basic: `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `filter`
- SV info: `SVTYPE`, `SVLEN`, `END`, `SOMATICSCORE`
- Counts: `normal_alt`, `tumor_alt`, `normal_sr_alt`, `tumor_sr_alt`
- BND: `bnd_remote_chrom`, `bnd_remote_pos`
- Mappings:
  - `gene_hits`: list of gene hit dicts
  - `elem_hits`: list of element hit dicts
  - `lncRNA_hits`: list of lncRNA hit dicts
  - (optional) miRNA-related columns when miRNA table supplied
- VEP:
  - `gene_hits_vep`, `regulatory_hits_vep`, `motif_hits`
  - `gene_symbols`, `hits_canonical`
  - `has_missense`, `has_nonsense`, `has_frameshift`, `has_splice_effect`
  - `has_*_canonical`, `has_*_mane`
- Motifs (after FIMO):
  - `flank_motif_hits`: motifs in flanking regions
  - motifs also attached under `elem_hits[*].motif_hits`
- After neojunction (STEP 6): `neojunction_motif_hits` (serialized)
- After ChIP (STEP 7): `chip_hits` / aggregates per `DATA_STRUCTURE_ATLAS.md`

### `gene_hits` structure (conceptual)

```python
{
    "gene_name": str,
    "gene_id": str,
    "strand": str,
    "signed_dist": int,
    "overlap_start": int | None,
    "overlap_end": int | None,
    "overlap_bp": int,
    "overlap_percent": float,
    "promoter_flag": 0 | 1,
    "gene_body_flag": 0 | 1,
    "exon_flag": 0 | 1,
    "intron_only_flag": 0 | 1,
    "upstream_5kb_flag": 0 | 1,
    "downstream_5kb_flag": 0 | 1,
    "region_hit": str,
    "hit_side": "span" | "point" | "bp1" | "bp2",
    "stop_codon_flag": 0 | 1,
    "start_codon_flag": 0 | 1,
    "transcript_id": str | None,
    "transcript_type": str | None,
}
```

### `elem_hits` structure (conceptual)

```python
{
    "elem_id": str,
    "elem_type": str,
    "chrom": str,
    "elem_start": int,
    "elem_end": int,
    "signed_dist": int,
    "overlap_start": int | None,
    "overlap_end": int | None,
    "overlap_bp": int,
    "overlap_percent": float,
    "overlaps_flag": 0 | 1,
    "region_hit": str,
    "proximal_flag": 0 | 1,
    "distal_flag": 0 | 1,
    "hit_side": str,
    "motif_hits": [
        {
            "start": int,
            "end": int,
            "TF": str,
            "motif_id": str,
            "score": float,
            "p_value": float,
            "q_value": float,
            "strand": str,
            "distance_to_pos": int,
        }
    ],
}
```

---

## Dependencies

### Python packages

- pandas, numpy, vcfpy (and project `pipeline.*` imports)

### External tools

- VEP (Variant Effect Predictor)
- FIMO (MEME suite)
- bedtools

---

## Integration with the main regulatory pipeline

SV tables consume the same **gene** and **regulatory element** artifacts as the cCRE-centric pipeline (paths via `PATHS` / `working_dir`). After **`run_sv_pipeline`**, SV rows reference the same `elem_id` universe as `coding_element_focus` exports where inputs match.

---

## Notes

- The shell master script writes a pipeline log plus per-step logs under `OUTPUT_ROOT/logs/`.
- The `.done` marker mechanism supports robust resume after failures.
- For best performance, expose parallelism via:
  - VEP: `--fork <N>`
  - FIMO: `--jobs <N>` (parallelize across FASTA files)
