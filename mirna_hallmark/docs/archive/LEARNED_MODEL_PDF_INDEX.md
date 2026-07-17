# Learned model — PDF bundle reading guide (what's what)

This bundle collects the four core documents of the **learned miRNA-arm → Hallmark-gene
regulatory model** (`mirna_hallmark/learned/`) as standalone PDFs, plus this guide. Read them in
the order below; each answers a different question and they are written to be held side by side.

### 1. LEARNED_MODEL_DESIGN_RESPONSE — *why this design at all?*

The decision log. The full reasoning for building a program-wise hierarchical Bayesian
occupancy-regression, with every major choice carrying an *Alternatives kept in play* block. Start
here for the big picture and the rationale behind the architecture.

### 2. LEARNED_MODEL_METHODS — *what is the actual formula?*

The clean formula-level spec of every estimator: notation, objective, formula, one line of what it
computes, section by section (§0–§17). The reference sheet, deliberately kept free of prose.

### 3. LEARNED_MODEL_RATIONALE — *why this construction and not the naive one?*

The reasoning companion to METHODS; its section numbers mirror METHODS exactly. Why each
construction, what breaks otherwise, what the estimator really identifies, and how to read its
output honestly.

### 4. LEARNED_MODEL_VALIDATION — *does it actually work?*

The living evidence dossier: real numbers, recovered positive controls, sharpening/novelty
examples, and the honest negatives (§4), each with the command to reproduce it.

## How they fit together

- **DESIGN_RESPONSE** is the *why-at-all* — the argument for the whole approach, written before/around the build.
- **METHODS** is the *spec* — the frozen math. It is the single source of truth for notation and formulas.
- **RATIONALE** is the *why-this-formula* — it shadows METHODS section-for-section with the reasoning a formula sheet cannot carry.
- **VALIDATION** is the *does-it-work* — the empirical record that each estimator earns its place.

A typical path: skim **DESIGN_RESPONSE §1** (executive summary) → pull up **METHODS** and **RATIONALE**
together for any estimator you care about → check **VALIDATION** for whether that estimator is verified,
carried from a prior run, or still to be (re)generated.

## Not in this bundle (referenced by the four docs)

- `LEARNED_MODEL_ESTIMATOR_MAP.md` — which estimator does which job, and chosen over what.
- `LEARNED_MODEL_BUILD_PLAN.md`, `LEARNED_MODEL_WHATS_NEXT.md` — build sequencing and open items.
- `LEARNED_MODEL_DISCUSSION_PROMPT.md` — the §9 prompt DESIGN_RESPONSE answers.
- `FORMULAS.md`, `MODELING_FRAMEWORK.md`, `EDGE_QUESTION_TAXONOMY.md`, `GENE_QUESTION_TAXONOMY.md` — surrounding heuristic/framework context.

*Status glyphs in VALIDATION are rendered as text tokens in these PDFs: `[OK]` = verified this arc,
`[FILE]` = from a prior run (cited), `[TODO]` = to (re)generate under canonical M.*
