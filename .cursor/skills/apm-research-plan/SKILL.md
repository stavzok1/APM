---
name: apm-research-plan
description: Guides any agent working on the APM immune-visibility PhD research plan. Use when the user asks about hypotheses, aims, scores (CRML, BCD, R_g, IVI, PRD, JSI, DAV, COMB_HLA, ARCH, BOTTLE, COHER), the conservation-prior methodology, the 13 per-tumor HiChIP samples, the PhD timeline, the gene panel (core 66 or extended tiers), DAG/diagram updates, or any task that touches research_plan/ in this repository.
---

# APM Research Plan

Cross-cutting research plan for a PhD thesis on tumor-intrinsic immune visibility built on top of the `pipeline/` infrastructure. **All substantive content lives in `research_plan/`.** This skill tells the agent how to operate on that folder.

## First action

Always read `research_plan/README.md` before anything else. It is the navigation for the whole plan.

## Folder map

```
research_plan/
├── README.md                         # nav + scope + status
├── 01_gene_panel_extended.md         # canonical tiered panel (all symbols + biology)
├── extended_genes.md                 # pointer → 01_gene_panel_extended.md
├── 02_hypothesis_catalog.md          # H1–H65 with schema; families A–N
├── 03_architecture.md                # L0–L5, feature builder contracts
├── 04_conservation_prior.md          # methodology (13 HiChIP anchor)
├── 05_timeline.md                    # PhD calendar, three paths, five gates
├── 06_automation_plan.md             # hypothesis tiers + ranked waves
├── 07_wet_lab.md                     # W1–W6 experiments, minimum package
├── 08_hla_allele_layer.md            # per-allele builders, Family M details
└── dag/                              # Graphviz .dot sources
    ├── README.md
    ├── 01_biological_circuit.dot
    ├── 02_pipeline_architecture.dot
    ├── 03_within_panel_relations.dot
    ├── 04_noncoding_space.dot
    └── 05_hypothesis_dependency.dot
```

## Common tasks

### Adding a new hypothesis

1. Append a new entry to the correct family in `02_hypothesis_catalog.md` using the standard schema: ID, Statement, Atlas dependencies, Derived features, Automation tier, Lead-set flag, Priority, Claim coverage.
2. Add its feature dependencies to `03_architecture.md` if any new builder is required.
3. Add the hypothesis node to `dag/05_hypothesis_dependency.dot` in the correct family cluster with edges to the feature builder nodes.
4. If it changes the execution order, update `06_automation_plan.md` ranking.

### Adding a new gene

Decide Tier (1: full integration, 2: medium, 3: CNV only, 4: readout only) per `research_plan/01_gene_panel_extended.md` criteria. Only promote to Tier 1 if it closes a circuit gap or mechanistically unblocks a lead hypothesis.

### Adding a new feature builder

1. Specify the builder in `03_architecture.md` under the correct layer (L1 raw, L1 novel, L1.5 conservation-normalized, L2 composite).
2. State the sample key (`sample`) and gene key (`gene_name`), unit, missingness policy, caching path.
3. If it adds a new score, register it in the L2 table.
4. Update the appropriate DAG (`02_pipeline_architecture.dot` or `05_hypothesis_dependency.dot`).

### Updating timelines or priorities

Edit `05_timeline.md` (calendar) and `06_automation_plan.md` (hypothesis rank). Keep the three paths (fast / realistic / conservative) coherent: if a gate slips, it slips in all three.

### Editing DAGs

Edit `.dot` files directly. Render locally or paste into https://dreampuf.github.io/GraphvizOnline. DAG conventions:

- `shape=box` for data/tables, `shape=ellipse` for biological entities, `shape=diamond` for composite scores.
- `solid` = direct mechanism, `dashed` = statistical/inferred, `dotted` = negative/inhibitory.
- Colors: `lightblue` data, `lightyellow` feature/score, `lightgreen` biological, `mistyrose` counter-regulator, `lightgrey` infra.

Do not replace DOT with Mermaid — DOT was chosen deliberately because the graphs are too dense for Mermaid's auto-layout.

## Five persistent claims to keep in mind

- **C1** — gene- and subtype-specific mechanism of visibility loss
- **C2** — element-resolved causal inference via CRML, BCD, conservation-weighted features
- **C3** — proteostatic visibility (RPPA discordance, signaling-refractory classes)
- **C4** — platform generalizes (pathway-agnostic fingerprint)
- **C5** — intrinsic visibility is prognostic (OS/PFS beyond stage + PAM50)

Every new hypothesis should identify which claim(s) it reinforces.

## Four hard-coded facts

1. **Per-tumor HiChIP count = 13** (`tcga_brca_hichip.csv`). They calibrate the subtype-consensus prior; they are not the cohort.
2. **PRIMARY_GENES = 66** (frozen core) in `pipeline/config.py`. Default pipeline work uses **`PIPELINE_GENE_PANEL`** (`FULL_INTEGRATION_GENES` when extended mode is on). Tier lists: `EXTENDED_PRIMARY_GENES`, `TIER2_MEDIUM_GENES`, `TIER3_CNV_ONLY_GENES`, `TIER4_READOUT_GENES`; **`CNV_GENES`** is the single deduped CNV union. Disable extension with env `APM_USE_EXTENDED_GENE_PANEL=0`.
3. **24 biosamples** in the Element Focus `TAD_domains` consensus. Do not re-derive; use the column.
4. **HLA class-I typing available for the full cohort**. Per-allele outcome variables (Family M) are in scope, not deferred.

## Pre-registered lead hypothesis set (15)

H1, H4, H9, H12, H14, H18, H19, H38, H41, H46, H51, H53, H56, **H59** (TNBC four-class), **H63** (regulatory ASE). These are the minimum the thesis must defend.

## When to update the plan

Update `research_plan/` whenever:

- New per-tumor data arrives (e.g., more HiChIP samples) → `04_conservation_prior.md` + `05_timeline.md`.
- A feature-builder contract change happens → `03_architecture.md`.
- A hypothesis finishes (success or null) → `02_hypothesis_catalog.md` plus a line in `negative_results.md` (create on first null).
- An external cohort is secured → `05_timeline.md` gate G5 + add replication rows to affected hypotheses.
- A wet-lab experiment completes → `07_wet_lab.md` status + figures; if a core claim (C1/C2/C3) gains causal evidence, cite the W-ID in the matching hypothesis entry.
- HLA typing, LOHHLA, or neoantigen-per-allele pipelines progress → `08_hla_allele_layer.md` and gate G4 in `05_timeline.md`.

## What NOT to do

- Do not duplicate pipeline-module docs here. For table contracts, read `pipeline/md/DATA_STRUCTURE_ATLAS.md`.
- Do not write code inside `research_plan/`. Analysis code belongs in `analysis/` and must match contracts in `03_architecture.md`.
- Do not collapse the three timeline paths. They intentionally differ in assumed slippage.
- Do not expand the lead-set beyond 13 without justification; it is an explicit commitment device.
