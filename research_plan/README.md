# APM immune visibility — research plan

This folder is the **single source of truth** for the scientific scope, hypothesis catalog, analysis architecture, conservation-prior methodology, timeline, and DAGs for the APM immune-visibility thesis built on top of the `pipeline/` infrastructure.

It is deliberately kept outside `pipeline/md/` because it is cross-cutting (analysis + pipeline + external cohorts + thesis planning) and is not a module-level doc.

---

## Navigation

| File | Purpose |
|---|---|
| `01_gene_panel_extended.md` | **Canonical** tiered panel: every symbol in `pipeline/config.py`, justification, categories, connections, miRNA table. |
| `extended_genes.md` | Pointer to `01_gene_panel_extended.md` (legacy filename for links). |
| `02_hypothesis_catalog.md` | The full catalog **H1–H65**, with Atlas dependencies, composite features, automation tier, lead-set flag, and priority. Families M (HLA per-allele) and N (TNBC hot/cold) are the newest. |
| `03_architecture.md` | Layers L0–L5, per-sample × per-gene feature builders, mediation engine, scoring interface. |
| `04_conservation_prior.md` | The conservation-as-prior methodology. **Includes the 13 TCGA-BRCA per-tumor HiChIP samples as the calibration anchor** for the subtype-consensus prior. |
| `05_timeline.md` | PhD calendar, three paths (fast / realistic / conservative), gated by automation tier. |
| `06_automation_plan.md` | Hypothesis-level automation tiers, batch test design, ranking, human-in-the-loop steps. |
| `07_wet_lab.md` | Wet-lab validation roster (W1–W6). Recommended minimum package: **W1 (CRISPRi enhancer) + W3 (DNMT/EZH2 reversal)**. |
| `08_hla_allele_layer.md` | Per-allele HLA layer: allele-aware feature builders + Family M (H46, H63–H65). |
| `dag/` | Graphviz `.dot` sources for the five canonical diagrams + render instructions (`dag/README.md`). |

---

## Core scientific claim

Tumor-intrinsic immune visibility is a **decomposable, element-resolved, multi-omic property** whose dominant regulatory mechanism varies by gene, subtype, and patient, and can be measured with novel element-resolved quantities (CRML, enhancer redundancy, boundary-conservation-weighted disruption, evidence-weighted miR pressure, protein–RNA discordance, pathway coherence) that jointly predict immune infiltration, protein-level APM competence, and clinical outcome.

## Three aims, three papers

| Aim | Focus | Lead paper target |
|---|---|---|
| A1 Atlas | Per-gene variance decomposition + element-resolved mediation across the APM/NK/IFN circuit | *Cell Reports* / *Nat Commun* |
| A2 Scores | Define and clinically validate CRML, R_g, BCD_g, miR pressure, PRD, IVI; survival + external cohort | *Cancer Cell* / *Nat Cancer* |
| A3 Platform | Pathway-agnostic **mechanism fingerprint**, lead case APM + contrast pathways DDR & EMT | *Genome Biology* / *Nat Methods* |

## Key assets already in the repo (don't duplicate effort)

- `pipeline/` — regulatory-element–centric tables with SCREEN + ABC + HiChIP + ChIP + TAD (24 biosamples)
- `pipeline/config.py` — `PRIMARY_GENES` (66), `PATHS`, `THRESHOLDS`, `BIOSAMPLES`
- `pipeline/md/DATA_STRUCTURE_ATLAS.md` — per-table nested contracts
- Per-sample modules: RNA, SNV (with `motif_hits` Δscore + `cCRE_hits`), SV (with `flank_motif_hits`, `chip_hits`), CNV (`gene_hits`, `elem_hits`), Methylation (per-sample probe/gene/cCRE aggregation), ATAC case-level, RPPA (panel scores, `signaling_blocks`, `protein_rna_discordance`, `estimated_apm_capacity`), miRTarBase (evidence-weighted)
- Unified clinical + Thorsson immune + PAM50
- **13 per-tumor HiChIP samples** (`tcga_brca_hichip.csv`) — calibration anchor for conservation-as-prior

## How to use this folder

- If editing scope → update `02_hypothesis_catalog.md` and `01_gene_panel_extended.md` (keep `pipeline/config.py` tier lists in sync).
- If building analysis code → follow contracts in `03_architecture.md`.
- If deciding priorities → read `05_timeline.md` + `06_automation_plan.md`.
- If a future agent joins → start at `.cursor/skills/apm-research-plan/SKILL.md`, which points back here.

## Status

- Pipeline: mature, per-module Atlas contracts frozen.
- Analysis layer (`analysis/`): partially in place; feature builders contract not yet frozen.
- Hypothesis catalog: drafted (**H1–H65**), lead pre-registration set of **15** marked in `02_hypothesis_catalog.md`.
- External cohort access: METABRIC assumed; ICB cohort (IMvigor210 / IMpassion130 / Liu / Van Allen) pending.
- HiChIP tumor samples: **13 processed**; incorporated into conservation-prior calibration (see `04_conservation_prior.md`).
- HLA class-I typing: **available for the full cohort**; unlocks Family M (per-allele layer, `08_hla_allele_layer.md`).
- Wet-lab: available as a parallel track; roster and minimum package in `07_wet_lab.md`.
