# Handoff: CPTAC orphan-edge discovery → next analysis

**For an agent picking up the orphan-edge thread.** Read `mirna_hallmark/AGENTS.md` first, then
this. Full session context: `docs/CPTAC_VALIDATION_SESSION_SUMMARY.md`. Registry: MH-36.

---

## Where things stand (concise session recap)

The miRNA-Hallmark pressure spine was validated against CPTAC proteome:
- **MH-33/34:** spine validates at protein in the **independent CPTAC-2 prospective** cohort
  (ZEB1/EMT; 53/1143 genes FDR) but not the same-patient TCGA-105 iTRAQ. Repression is
  **mostly mRNA-mediated** (a slope-free `protein_resid` layer isolates the small translational part).
- **MH-35:** the residual survivors are driven by **known** subtype-structured oncomiR arms
  (PTEN←miR-17~92 hub at protein).
- **MH-36 (this thread):** an **orphan**-edge scan (`cptac_orphan_discovery.py`) tested 20,663
  miRNA→gene edges with a TargetScan/ENCORI prior but **no miRTarBase functional-MTI curation**
  against CPTAC protein, and nominated novel candidate edges. Method self-validated by recovering
  the miR-29→collagen axis de novo.

## The orphan output you're inheriting

`output/cptac_validation/orphan_discovery/`
- `orphan_edges_tested.tsv.gz` — all 15,512 tested edges, both-cohort stats, priors, tags.
- `orphan_candidates.tsv` — the 594 FDR-sig (protein or `protein_resid`) candidates.
- `method_manifest.json`.

**Candidate tiers** (column `strength_tag`):
| tag | meaning | n |
|-----|---------|---|
| `S_translational_replicated` | `protein_resid` FDR-neg (prospective) **and** TCGA-105 sign-concordant | **93** |
| `S_protein_replicated` | protein FDR-neg + TCGA-replicated | 237 |
| `P_translational_prospective` | `protein_resid` FDR-neg, not replicated | 89 |
| `P_protein_prospective` | protein FDR-neg only | 175 |

**Highest-confidence shortlist** = `translational_candidate & tcga_replicated & ~mirtar_any &
prior=="targetscan+encori"` → **15 edges** (fully uncurated, sequence **and** CLIP prior, translational, two-cohort).
Leaders: miR-30e→MINPP1/VIM, miR-29a→COL11A1/COL6A2/BMP1, miR-17→NEDD4L, miR-15b→AP2B1, miR-16→FKBP1A.

**Hard caveat to carry forward:** these are **wet-lab NOMINATIONS**, not validated edges. Protein
anti-correlation ≠ direct edge. The dominant confounds to break next are **(a) co-expression**
(ECM/collagen genes co-vary with stromal content despite purity+CIN adjustment) and **(b)
seed-family ambiguity** (miR-29a/b/c all hit; cannot attribute to one arm).

## Suggested next steps (pick up here)

1. **Breast-PMID + fresh-literature triage of the 15-edge shortlist.** Run the MH-31
   `edge_breast_context` fetch over just these edges (and a live PubMed pass) to separate
   "uncurated-but-known-elsewhere" from "genuinely uncharted" — mirror MH-26's `literature_status`
   triage. Expectation from MH-26: most resolve to under-curation, a small tail is truly novel.
2. **Co-expression / indirect-path control.** For each candidate, re-fit the `protein_resid`
   coupling adding the target's **own mRNA-state neighbours** or a stromal/EMT composite as a
   covariate; an edge that survives is less likely pure co-expression. ECM candidates especially
   need a stromal-fraction (xCell/ESTIMATE StromalScore from the CPTAC CLI) adjustment beyond CIN.
3. **Seed-family disambiguation.** Where multiple family arms hit one gene (miR-29a/b/c→COL11A1),
   test which arm's *abundance* best explains the residual and whether the TargetScan site is
   shared (seed) or arm-specific — you usually can only claim the family, not the arm.
4. **Single-edge pressure (the open WHATS_NEXT 3.3 refinement).** The scan used arm *expression*;
   build true per-(arm,gene) pressure via `pressure_build.compute_gene_pressure_contributions`
   restricted to the candidate edge and re-test — sharper than expression alone.
5. **CLIP cross-check.** For the ENCORI-prior candidates, pull the actual CLIP `clipExpNum` /
   binding-site evidence from `encori_edges` and rank by physical-binding depth — a candidate with
   strong CLIP + sequence + two-cohort protein residual is the strongest wet-lab nomination.
6. **Register** any refined shortlist as a new MH-# and update the ledger/catalog (don't overwrite MH-36).

## Key code / data pointers

- Module: `mirna_hallmark/cptac_orphan_discovery.py` (orphan def in `build_orphan_edges`; screen in `screen_cohort`).
- Orphan priors: `targetscan_orphan_coupling.build_orphan_edge_table` (TargetScan), `encori_edges.load_collapsed_encori_pairs` (CLIP).
- CPTAC layers: `cptac_validation.load_cptac_layers("prospective"|"tcga105")` → `protein_z`, `protein_resid`, `gap`.
- Prospective miRNA arms: `cptac_validation.load_prospective_mirna_arms()` (built from `data/CPTAC/cptac2_brca_mirna_isoform_log2rpm.tsv`).
- Literature layer: `output/tissue_reference/edge_breast_context/edge_breast_context_scores.tsv` (MH-31).
- Do **not** treat protein anti-correlation as edge validation in any registered claim — keep the P/H tag.
