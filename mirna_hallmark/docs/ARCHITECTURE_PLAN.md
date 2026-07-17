> **Goal:** the plan for ONE canonical, non-redundant, maintainable architecture map that joins the four
> dimensions of the work — **axes × models × analyses × results** — as a materialized VIEW over the existing
> one-home docs (never a content-duplicating megadoc).
> **What belongs here:** the plan, the axis taxonomy (the join key), the map schema, the generator design, and
> the tagging convention. NOT the map itself (that is generated → `ARCHITECTURE.md` + artifact) and NOT any
> finding content (that stays in its one home).
> **Update trigger:** when the axis taxonomy changes, or a phase completes.
> **Sync-partner:** `STATE_OF_PLAY.md` (per-axis verdicts — the map links to it), `DISCOVERY_REGISTRY.md`
> (results — the map links to MH-##), `ANALYSES_CATALOG.md` (analyses — the map links to modules).

# Canonical architecture — the plan (2026-07-18)

## The diagnosis
Each axis is discussed in **6–20 docs** ("attribution" 20, "discovery"/"CPTAC"/"subtype" 19). This is the
signature of a schema **normalized by document TYPE** (registry=results, catalog=analyses, framework=models,
state=axis-verdicts). That normalization is *deliberately good* — it enforces the one-home rule that stops
retraction-rot. Its known cost: **there is no JOIN VIEW.** Understanding one axis end-to-end means reconstructing
the join across ~20 docs by hand. That join is the missing artifact.

## The design principle — a materialized view, not a megadoc
Keep the normalized tables (the one-home docs). Add a **VIEW that joins them BY REFERENCE.** Every cell of the
map is a **link/tag into a normalized home, never copied content.** When a finding updates in the registry, the
map's link stays valid because the map never held the content. Because it is a view, it is **generated**, not
hand-maintained — refreshed from its sources like a DB view refreshes from tables.

## Decisions (2026-07-18, user)
1. **Form = BOTH** — the generator emits the canonical markdown `ARCHITECTURE.md` (in-repo, source of the view)
   AND a rendered interactive artifact (click an axis → its models/analyses/results/status). The markdown is
   canonical; the artifact is a rendered view of it.
2. **Map first, converge later** — build the join view over the docs AS THEY ARE. Defer merging overlapping
   type-docs (RATIONALE / VALIDATION / FRAMEWORK / METHODS / FORMULAS) to a later pass; the map alone removes
   most cross-referencing pain.
3. **Tagging = frontmatter / registry-column** — each `MH-##` row gets an **axis** column; each module gets an
   `# axis: <name>` comment. Enforced at creation like the one-home rule. Self-describing, travels with the
   artifact. (NOT a central axis_map.tsv; NOT auto-infer.)

## The map schema (one row per AXIS; every cell a link/tag, no prose)
| axis | scientific question | model(s) | key analyses (→catalog/code) | current verdict (→ MH-##) | status | open (→board) |

## Phase 0 — the axis taxonomy (the JOIN KEY) — ✅ FINALIZED 2026-07-18 (user)
The closed list of **12 axes**. Every module, MH-row, and model tags to exactly ONE **primary** axis (it may
touch others, but has one home-axis). A finding that fits none is a signal to consciously add an axis, not to
sprawl. **Resolutions (user):** `coupling` MERGED into `edge-existence` (one question: do these edges exist and
act); `cn-causal` and `discovery` are their OWN axes — both **extensible beyond edge-existence** (CN dose informs
dosage/causal questions broadly; discovery is a full lane — site-free null, convergent evidence, subtype-
stratified — not just edge-finding).

| axis-tag | scientific question | primary home (STATE_OF_PLAY §) |
|---|---|---|
| `model` | what the learned estimator IS (dense Gibbs posterior, 2 readouts, gauge, site-free null) | Axis 1 |
| `edge-existence` | do the curated/predicted edges exist AND act (coupling \| C, exogenous validation)? | Axis 2 |
| `cn-causal` | CN-dose causal identification (instrument, exclusion) — **its own axis, extensible to dosage** | CN_INSTRUMENT |
| `attribution` | WHO owns a gene's regulation — identity vs magnitude (Shapley)? | Axis 3 |
| `decoy` | does curation beat an abundance/variance-matched fake (specificity)? | Axis 4 |
| `discovery` | novel edges beyond curation — **its own lane** (site-free null, convergent evidence, MH-155/156) | Axis 2/E |
| `protein` | does it hold at the protein layer (CPTAC)? | Axis 5 |
| `progression` | GTEx→NAT→tumor trajectory / state | Axis 6 |
| `subtype` | PAM50-stratified coupling / who-is-pressured | (new; §F) |
| `outcome` | does it predict survival (prognostic)? | Axis 7 |
| `external` | independent-cohort replication (Buffa / METABRIC / SCAN-B) | Axis 7 |
| `dcis-ev` | DCIS / pre-malignant / extracellular-vesicle lane | §H |

## Phase 1 — the linkage layer (mechanical, mostly auto-extractable)
Add the axis tag at the three sources:
- **registry:** an `axis` value per `MH-##` row (parse the existing table; assign from its content).
- **modules:** an `# axis: <name>` line per module (top-level + learned + eval + analyses + method_dev).
- **models:** MODELING_FRAMEWORK sections tagged.
Auto-extract the rest: module→outputs (catalog), MH→evidence (registry `evidence` column), import graph (the
session's tooling). The ONLY semantic step is the axis assignment; everything else joins automatically.

## Phase 2 — the generator (the view)
A script (`docs/gen_architecture.py` or `analyses/ops/`) that JOINS the linkage layer against
registry + catalog + code and emits: (a) `ARCHITECTURE.md` (the canonical matrix, per-axis, all links);
(b) an HTML artifact (rendered, navigable). Regenerable on demand ⇒ cannot drift.

## Phase 3 — convergence (DEFERRED per decision 2)
Once the view exists, audit the type-docs for role overlap (RATIONALE vs VALIDATION vs FRAMEWORK; METHODS vs
FORMULAS) and collapse to a minimal non-redundant set. Safe only AFTER the map exists (docs that only provided
cross-axis context can then collapse into pointers).

## Phase 4 — maintenance
The map regenerates from the Phase-1 tags. Discipline: every new MH-row and module gets its `axis:` tag at
creation (enforced like the one-home rule). The map is NEVER hand-edited — edit the tags, regenerate.

## Guardrails (so the map does not become the problem it solves)
- The map LINKS, never restates. A cell with finding-content is a bug — it must be an MH-## link.
- The registry stays the source of truth; the map is downstream.
- Generate the mechanical parts; hand-curate ONLY the axis taxonomy + assignment.
- Closed axis list — a misfit finding means consciously adding an axis, not sprawling the taxonomy.
