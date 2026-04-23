## lncRNA interactions module (miRNAs + RBPs)

This module builds a compact, deterministic set of interaction tables for the lncRNAs you are actively
modeling (Tier‑1 panel lncRNAs + a small “closest-to-primary-genes” extension).

It integrates:
- **miRNA→lncRNA binding evidence** from **ENCORI/starBase** (AGO CLIP–supported).
- **RBP→RNA binding peaks** from **POSTAR3** (CLIP/eCLIP peak catalogs), summarized for the selected lncRNA loci.

### 1) Inputs

#### Selected lncRNAs (deterministic)

We define the lncRNA target list as:
- the 8 Tier‑1 panel lncRNAs in `pipeline/config.py` (`TIER1_LNCRNA_GENES`)
- plus **N** additional lncRNAs that are **closest to any of the frozen 66 primary genes**
  using `data/lncRNA_matching/genes_lncRNAs_1000000bp_distances.csv` and sorting by `min_distance_bp`

Default is **N=20**.

#### ENCORI (miRNA-target)

ENCORI endpoint used:
- `miRNATarget` (`assembly=hg38`, `geneType=lncRNA`, `target=<lncRNA gene symbol>`, `miRNA=all`, `cellType=all`)

ENCORI returns **arm-level miRNAs** (e.g. `hsa-miR-339-5p`) and `miRNAid` (MiRBase MIMAT).

#### POSTAR3 (RBP peaks)

Expected local files:
- Raw POSTAR: `data/RBP-RNA/POSTAR3.txt` (BED-like 10-column TSV, very large)
- Converted parquet (recommended): `data/RBP-RNA/POSTAR3.parquet`

POSTAR3 columns (as used here):
`chrom`, `start`, `end`, `strand`, `rbp`, `assay`, `cell_tissue`, `source_accession`, `score`

Feature-aware overlap modes:
- `gene`: whole lncRNA gene locus
- `exons`: exon union only (splicing/stability RBPs tend to enrich here)
- `introns`: gene locus minus exon union
- `promoter`: strand-aware TSS window (default upstream/downstream 2000/500bp)

### 2) Outputs (canonical)

All outputs are written under:
- `data/lncRNA_interactions/`

Files:
- `selected_lncrnas.txt` — the exact lncRNA gene symbol list used (reproducibility)
- `encori_mirna_targets.parquet` — ENCORI miRNA→lncRNA rows for `selected_lncrnas`
- `encori_mirna_targets_diagnostics.csv` — per-lncRNA diagnostics for ENCORI fetching (symbol vs ENSG fallback)
- `encori_rbp_targets.parquet` — ENCORI RBPTarget rows for a **small curated RBP list** (top-N from POSTAR3 overlap summary) × `selected_lncrnas`
- `encori_rbp_targets_diagnostics.csv` — per (RBP, lncRNA) diagnostics for ENCORI RBPTarget fetching (symbol vs ENSG fallback)
- `postar3_overlaps.parquet` — POSTAR3 peaks overlapping any selected lncRNA **gene locus**
- `postar3_rbp_summary.parquet` — per-RBP summary across selected lncRNA loci
- `postar3_lncrna_summary.parquet` — per-lncRNA summary (how many POSTAR peaks overlap its locus)
- `recommended_rbps_from_postar3.csv` — convenience shortlist (RBP + counts) derived from POSTAR3

### 4) Optional: predicted miRNA targets (tool-based)

Because ENCORI miRNA→lncRNA coverage can be sparse for many lncRNAs, you can generate **predicted**
miRNA→lncRNA binding calls locally using `RNAhybrid` on **spliced exon sequences**.

Entry point:
- `python -m pipeline.lncRNA_interactions.predicted_targets_cli --mirna-fasta <mature_miRNAs.fa>`

Outputs are written under:
- `data/lncRNA_interactions/predicted_targets/rnahybrid/`

Practical notes:
- If you pass `--mirna-subset-from-encori`, the CLI will write `mirnas.subset.fa` from the set of miRNAs
  observed in `encori_mirna_targets.parquet` (avoids all-miRs × all-lncRNAs blowups).
- Long spliced lncRNAs can make RNAhybrid abort (“target too long”). Default mitigation is
  `--max-target-bp 8000` which truncates each target to the **last 8kb** (3′ end). Use `--max-target-bp 0`
  to disable truncation.

### 3) How to run

From repo root:

```bash
python -m pipeline.lncRNA_interactions.build_all \
  --n-extra-close-lncrnas 20 \
  --encori-clip-exp-num 2 \
  --encori-rbp-clip-exp-num 2 \
  --encori-rbp-top-n 12 \
  --postar3-region-mode exons \
  --postar3-chunk-rows 500000
```

Notes:
- ENCORI **RBP-Target** is *not* used here because ENCORI only supports RBPs present in its CLIP universe;
  many chromatin modifiers are not valid ENCORI RBPs. Instead, we use POSTAR3 to derive a realistic
  “RBP shortlist” for the lncRNAs of interest.

