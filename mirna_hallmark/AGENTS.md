# Agent onboarding: `mirna_hallmark` subproject

> **This is a SEPARATE subproject** living inside the APM repository. It has its
> **own** documentation, analysis catalog, discovery registry, and run ledger
> under `mirna_hallmark/` — **do not** mix them with the main `analysis/` docs,
> `analysis/DISCOVERY_REGISTRY.md`, or `analysis/docs/ANALYSIS_RUN_LEDGER.md`.
> The root `AGENTS.md` only points here.

Read this before editing code or answering questions about `mirna_hallmark/`.

---

## 1. What this subproject is

A focused study of **miRNA regulation of MSigDB Hallmark programs** in TCGA-BRCA.
It reuses the parent APM data + loaders but writes everything under
`mirna_hallmark/output/`. Central questions:

1. Which Hallmark programs are most heavily targeted by **high-evidence**
   miRNAs (miRTarBase functional MTI + low-throughput validation)?
2. Does **miRNA pressure** (evidence-weighted repression) anti-correlate with
   Hallmark target-gene expression, **gated by AGO/RISC availability**?
3. How do the **miRNA universe** (loci) and Hallmark **target genes** vary by
   copy number across strata (PAM50 / TNBC / immune / stage)?

**Output lanes:** the core Hallmark spine (`run_all`) lives at `output/` top level.
**Tissue reference** (TCGA tumor / NAT / GTEx cross-state, UMAP, family depth) writes
to `output/tissue_reference/`. **Plasma EV** (circulating exosome GEO cohorts) is an
optional translational lane at `output/plasma_ev/` — not part of `run_all`.

## 2. Reuse, don't duplicate

The subproject imports from the parent repo rather than re-implementing:

| Need | Reused from |
|------|-------------|
| miRTarBase normalize | `pipeline.genes.mirtarbase.load_mirtarbase` (+ constants) |
| miRNA pressure | `mirna_hallmark.pressure_build` / `pressure_engine` — spine **softmax_z_logrpm + evidence_mass**; hybrid M7 **softmax_z + combined_mass** (parent `compute_pressure` is legacy) |
| miRNA-universe CNV | `analysis.cohort_landscapes.cnv.dosage_landscape_cnv` (catalog, ENSG map, cohort extraction, `summarize_by_strata`) |
| clinical / RNA / partial Spearman | `analysis.utils.common.loaders` |

**We do NOT call** the slow `pipeline.genes.mirtarbase.get_mirtarbase_targets`
for edge building — see `build_edges.compute_interaction_summary_fast` (vectorized).

## 3. Module spine (run order)

```
config.py              paths, strata (STRATUM_SPECS), AGO gate params, cutoffs
hallmark_sets.py       parse annotations/GSEA/h.all.v2026.1.Hs.txt (50 sets)
build_edges.py         vectorized Hallmark-scoped miRTarBase summary + edges
data_loaders.py        clinical strata, miRNA arms, RNA, AGO/RISC, target-gene CNV (cached)
stats.py               bh_fdr, kruskal_across_strata, hypergeom_enrichment, zscore_rows
ago_gate.py            per-sample RISC capacity + bounded saturating gate
mirna_locus_cnv.py     miRNA-universe CNV by stratum + CN->expression concordance
stratum_characterization.py  per-Hallmark target CNV + regulatory-miRNA expression by stratum
hallmark_interaction.py      AGO-gated pressure, anti-correlation, enrichment (CORE)
subtype_contrasts.py   luminal vs non-luminal (Normal aside) — cohort-wide lineage
pam50_gene_resolution.py  actual LumA/LumB/Her2/Basal per-gene pressure + drivers
visibility_archetype_contrasts.py  actual cold_Basal / hot_Luminal (APM archetypes)
robustness_checks.py   grant Aims 1-2: proliferation-adjusted coupling + curation-bias (evidence/binary/TargetScan) hub robustness
run_all.py             orchestrator -> output/run_manifest.json
```

Run anything as `.venv/bin/python3 -m mirna_hallmark.<module>` from the repo root.
Full pipeline: `.venv/bin/python3 -m mirna_hallmark.run_all` (use `--skip-mirna-cnv`
once the CNV cache exists; `--include-tnrc6` to add GW182 to the gate).

## 4. Key design decisions (read before changing behavior)

- **Gene universe** = union of all 50 Hallmark sets (~4,384 genes). A gene in N
  sets yields N `(miRNA, gene, hallmark_set)` edge rows; cohort-level counts
  dedup on `(miRNA, gene)`.
- **High-evidence edge** = ≥1 Functional MTI study AND a low-throughput
  validation (reporter or protein), else an `evidence_score` fallback. See
  `config.HIGH_EVIDENCE_*` and `build_edges._high_evidence_mask`.
- **AGO/RISC gate** is a **documented sensitivity layer, not a causal model**.
  Every interaction result is reported **gated AND ungated**. Gate ∈ [`gate_min`, 1].
- **Confounders**: inferential tests use partial Spearman adjusting CPE + HRD
  (and report the raw Spearman alongside), per APM analysis-conduct guardrails.
- **Caches** (rebuild with the relevant `--force*` flag):
  `output/matrices/cnv_target_genes.tsv.gz`, `output/matrices/cnv_ago_risc.tsv.gz`,
  `output/mirna_locus_cnv/tables/sample_entity_cnv.tsv.gz`.

## 5. Self-contained docs (this folder only)

| Doc | Purpose |
|-----|---------|
| `docs/INDEX.md` | **Map of every doc** — grouped by role (canonical / living / reports / parked / historical), with status + link graph. Start here to find a doc. |
| `README.md` | Quickstart + output map |
| `method_dev/` | **Method & conceptual development** subfolder, organized into 4 arcs (see `method_dev/README.md`): `aggregate_pressure/` (force-vs-abundance — `AGGREGATE_PRESSURE_FINDINGS.md` Q1–Q6: edge weight & promiscuity inert/harmful, **acquired abundance `max(a−h,0)` the one robust win**, tumor↔GTEx Δρ validated; + design log + followups), `arm_expression/` (canonical detectability floor docs), `site_ladder/` (seed-scan → filter ladder → experimental validation; code+data+doc), `landscape/` (`he_edge_arm_landscape.py` genome loci + `pten_competition.py`; code+data). Figures shared in `method_dev/figures/`. Run modules as `-m mirna_hallmark.method_dev.<subfolder>.<module>`. |
| `docs/HANDOFF_NORMAL_EXCLUDED_BATCH_RERUN.md` | **ACTIVE handoff** — Normal-like-discarded + batch-corrected rerun (running as systemd user service `mirna_rerun`); how to monitor/relaunch + the pending framework-doc exact-number refresh. Read if picking up that task. |
| `docs/MODELING_FRAMEWORK.md` | **Conceptual + implemented framework** — the connective layer (resolution hierarchy, edge universe, pressure, architecture, coupling/confounding, protein realization, healthy→tumor, validation). Start here for "how is this modeled and why" |
| `docs/BIOLOGICAL_SYNTHESIS.md` | **Biology-first findings surface + literature standing** (start here for "what did we learn") |
| `docs/MIRNA_CNV_DOSAGE_REPORT.md` | **Dedicated miRNA CNV dosage report** — segment-based locus CN, MIMAT paralog weights, PAM50 distributions, neighborhood co-amps, case interpretations |
| `docs/WHATS_NEXT.md` | Forward-looking extensions, discoveries, and caveat resolutions |
| `docs/REPORT.md` | Detailed results: tables, stratum differences, distinct miRNAs in context, literature comparison |
| `docs/METHODS.md` | Formal methods: pressure, gate, enrichment, concordance, stats |
| `docs/ANALYSES_CATALOG.md` | One row per analysis component (inputs/outputs/CLI) |
| `docs/DISCOVERY_REGISTRY.md` | Claims + strength tags, anchored to outputs |
| `docs/ANALYSIS_RUN_LEDGER.md` | Per-component last-run timestamp + status |

Update `DISCOVERY_REGISTRY.md` and `ANALYSIS_RUN_LEDGER.md` after any run that
changes a finding or refreshes an output.

## 6. Conventions inherited from the parent repo

- Cohort: TCGA-BRCA primary tumors; participant key = 12-char TCGA id; multi-vial
  collapse by mean (expression) / median (CNV). See parent
  `analysis/docs/COHORT_AND_JOIN_CONVENTIONS.md`.
- Running under WSL-from-Windows: follow root `.cursor/rules/wsl-shell-bridge.mdc`.
