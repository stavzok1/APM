# miRNA × Hallmark (TCGA-BRCA)

A self-contained subproject of the APM repository analyzing **miRNA regulation
of MSigDB Hallmark gene programs**, with evidence-weighted miRNA pressure,
**AGO/RISC availability gating**, copy-number characterization of both the
**miRNA universe** and Hallmark **target genes**, and stratification by PAM50 /
TNBC / Thornsson immune / pathologic stage.

> Onboarding for agents: see [`AGENTS.md`](AGENTS.md). Methods, catalog,
> discoveries, and the run ledger live under [`docs/`](docs/). These are
> **separate** from the main `analysis/` documentation.

## Quickstart

From the repository root (WSL), with the project venv:

```bash
# 1) build Hallmark-scoped miRTarBase edges (vectorized; ~30 s)
.venv/bin/python3 -m mirna_hallmark.build_edges

# 2) full pipeline (build -> gate -> miRNA-universe CNV -> stratum -> interaction)
.venv/bin/python3 -m mirna_hallmark.run_all --include-tnrc6

# faster reruns once the CNV cache exists:
.venv/bin/python3 -m mirna_hallmark.run_all --skip-mirna-cnv
```

Each step is also runnable standalone, e.g.
`.venv/bin/python3 -m mirna_hallmark.hallmark_interaction`.

## Inputs (all from the parent repo)

- `annotations/GSEA/h.all.v2026.1.Hs.txt` — 50 Hallmark gene sets
- miRTarBase raw (`data/miRNA/mirtar.csv`) + family info
- miRNA-seq log2(RPM+1), RNA log2(TPM+1), ASCAT3 gene-level CNV
- `annotations/BRCA_clinical_immune_unified.tsv` + Lehmann TNBC subtypes

## Output map (`output/`)

| Path | Contents |
|------|----------|
| `edges/mirna_hallmark_edges.tsv.gz` | `(miRNA, gene, hallmark_set)` edges + `high_evidence` |
| `edges/mirtarbase/mirtar_interaction_summary.csv` | Hallmark-scoped miRTarBase summary (vectorized) |
| `ago_gate/` | per-sample RISC capacity + gate; AGO expression/CNV by stratum |
| `mirna_locus_cnv/` | miRNA-universe CNV by stratum + CN→expression concordance; `maps/` contains segment-based locus/MIMAT genome maps and subtype-depth tables |
| `stratum_characterization/` | per-Hallmark target CNV + regulatory-miRNA expression by stratum |
| `hallmark_interaction/` | AGO-gated Hallmark pressure, anti-correlation, high-evidence target enrichment |
| `cptac_validation/` | CPTAC **TCGA-105** proteome validation (same-patient iTRAQ; gap/protein_resid/protein-anticorr × 4 pressure variants); see MH-33 |
| `cptac_validation/prospective/` | CPTAC **independent** prospective-breast validation (CPTAC-2; pressure self-contained from CPTAC-2 miRNA-seq) — ZEB1/EMT protein anti-corr validates; see MH-34 |
| `matrices/` | cached CNV matrices |
| `run_manifest.json` | last orchestrated run: per-step status + timing |

### Tissue reference (`output/tissue_reference/`)

Cross-state TCGA tumor / NAT / GTEx analyses, expression panels, UMAP, and
miR-301 family depth on **tissue** miRNA-seq. Not plasma EV.

### Plasma EV (`output/plasma_ev/`) — optional, not in `run_all`

Circulating exosome miRNA from GEO (`data/miRNA/ev/`). Modules: `ev_mirna_screening`,
`ev_mirna_replication`, `ev_mirna_multi_cohort`, `ev_mirna_followup`, `ev_mirna_deep`,
`ev_mirna_gse255660_screen`, `ev_mirna_gse270497_ts`. Kept in this subproject for
shared pressure/TargetScan helpers; candidate for a future standalone package if the
lane grows.

Every output directory carries a `method_manifest.json` describing how it was built.

## Headline results (current run)

- **Actual PAM50 per-gene resolution:** Basal shows strong pressure on EMT,
  apoptosis/tumor-suppressor and immune-priming nodes (PTEN, BCL2L11/BIM, FOXO1,
  CDKN1A, IRF1; drivers include miR-17/20a/106b/23a/130b). Her2 has a weaker
  TLR3/WNT5A/CENPF/IL15 pattern; LumB leans E2F/MYC/WNT/NF-kB; LumA carries
  TFAP2C/WNT2/luminal-differentiation pressure.
- **Cold/hot archetype checks are secondary:** cold_Basal vs other_Basal did not
  show a strong Hallmark-pressure split; most useful signal is intrinsic PAM50.
- Canonical miRNA-regulated programs (TGF-β, WNT/β-catenin, apoptosis,
  TNFα/NF-κB, EMT) are the most enriched for high-evidence miRNA targets;
  program-specific leaders (miR-29→EMT, miR-34a→WNT, miR-21/155→NF-κB) match the
  mechanism literature.
- **23/50** Hallmarks show a significant **negative** pressure↔expression
  coupling (q<0.05, CPE+HRD-adjusted); proliferation/immune positives are
  co-expression/purity confounds.
- **237/756** mature miRNA arms show positive locus-CN→expression concordance
  after the ASCAT3 segment refresh; 8q/chr1/chr17 miRNA amplicon neighborhoods
  are subtype-structured (LumB/Her2 high, Basal miR-200c gain).

**Full doc map:** [`docs/INDEX.md`](docs/INDEX.md) — every doc grouped by role, with status.

**Start here for the biology:** [`docs/BIOLOGICAL_SYNTHESIS.md`](docs/BIOLOGICAL_SYNTHESIS.md)
(findings surface + literature standing). Full tables and stratum detail:
[`docs/REPORT.md`](docs/REPORT.md). Dedicated miRNA CNV dosage report:
[`docs/MIRNA_CNV_DOSAGE_REPORT.md`](docs/MIRNA_CNV_DOSAGE_REPORT.md).
Next steps: [`docs/WHATS_NEXT.md`](docs/WHATS_NEXT.md). Tracked claims:
[`docs/DISCOVERY_REGISTRY.md`](docs/DISCOVERY_REGISTRY.md).
