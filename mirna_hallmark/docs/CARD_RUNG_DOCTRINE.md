# Card & rung doctrine — what unit does this number live on?

> **Scope: every column on every card.** The canonical source is the `learned/card_rungs.py` module
> docstring and `CARDS` dict; this doc is its doc-level home and the entry point for "which rung am I on".
> Verdicts live in `DISCOVERY_REGISTRY.md` (MH-179, 187, 188, 191, 214, 222, 226, 227).

## 0. The one rule

> **⛔ THE RUNG IS A PROPERTY OF `(CARD, COLUMN)` — NEVER OF THE COLUMN NAME.**

`beta` means a different thing on each card:

| card | what `beta` is |
|---|---|
| `edge_card.tsv` | **EDGE** rung — `readouts.run(level="arm")`, the SS8 collapse **REMOVED** |
| `gene_family_card.tsv` | **FAMILY** rung — `readouts.run(level="family")`, the SS8 collapse **APPLIED** |

Same name, same estimator, different unit. A single shared prefix map cannot express that, which is why
each card carries its own map and why `domain_of(col, card)` is **card-scoped** — prefixes were once
matched globally and silently attached one card's caveat to another's columns (three bugs; `sub_` on the
arm card vs the edge card, fixed by renaming to `esub_`).

**Four defects in one day came from not knowing the unit:** MH-179 (a FAMILY-estimated β applied to RAW
ARM abundance) · MH-187 (a family weight beside an arm correlation, unmarked) · MH-188 (a within-cell
Shapley compared against a GENE-level OOF statistic) · MH-191 (β labelled `family` when the fit was
per-arm). All the same error.

## 1. The five cards

| rung | card | key | expresses |
|---|---|---|---|
| **edge** | `realization/edge_card.tsv` | `[gene, arm]` | one (miRNA, gene) pair |
| **family_edge** | `gene_family_card.tsv` | `[gene, family]` | one (seed family, gene) pair |
| **gene** | `realization/gene_card.tsv` | `[gene]` | the gene's total incoming regulation |
| **arm** | `arm_card.tsv` | `[arm]` | the arm itself, gene-free |
| **seed_family** | `seed_family_card.tsv` | `[seed_family]` | the family itself, gene-free |

⚠ **`family_card` is keyed `[gene, family]` — it is the GENE×FAMILY card and CANNOT express a property of
the family itself.** That is exactly why `seed_family_card` was split off as the 5th rung. If your
quantity is a property of the family alone (size, seed heterogeneity, member composition), it belongs on
`seed_family_card`, not `family_card`.

## 2. Two orthogonal facts per column — do not conflate them

- **`rung`** — the unit the value is **defined on**: `key` / `gene` / `family` / `arm` / `edge` /
  `arm-in-family`.
- **`agg_of`** — if the value **summarises a lower unit**, which one. E.g. on the gene card
  `cptac_prosp_agg_rho_prot` is `rung=gene, agg_of=arm`: one row per gene, but computed by summing
  `β·X` over the gene's ARMS, so **it inherits the arm rung's caveats**.

⚠ **A `domain` entry is not a rung.** Domain says *where a column is defined*; rung says *what unit it
lives on*. `--check` needs both (MH-214).

⚠⚠ **AND `domain` IS NOT A DESCRIPTION OF WHAT THE COLUMN MEANS.** It is a **row-applicability**
statement — *which rows the column is defined on* — e.g. *"multi-arm cells only (n_arm_in_cell > 1) —
20.3% of edges"*, *"edges with a matched NAT leg (n~104 paired participants)"*. ⇒ **a column with NO
domain entry is defined on EVERY row, which is the correct default, not a gap.** Measured 2026-08-04:
no-domain columns are **91–100%** populated (arm card 1.000, edge 0.995, gene_family 0.985) while
with-domain columns are **44–77%** — exactly what the definition predicts.
⛔ **Do not "fill in the missing domains".** Reading the field name rather than its values made me
report a nonexistent "179 undescribed columns" gap twice; the check that refuted it was one query.
✅ **THE MEANING GAP IS CLOSED (2026-08-18).** This block used to end *"there is separately NO per-column
MEANING documentation anywhere ... it lives only in the producing modules' docstrings."* `learned/card_glossary.py`
now carries one for **all 702 columns (100%)**, and it is **CARD-SCOPED for the same reason `domain_of` is**:
`describe("beta", card="edge")` and `describe("beta", card="gene_family")` return *different* text, because
they are different quantities. Resolution order: exact `(card, column)` override → the card's own prefix
block → the global prefix block → **`None`, counted as missing and never guessed**.
```
python3 -m mirna_hallmark.learned.card_glossary            # coverage report
python3 -m mirna_hallmark.learned.card_glossary --emit     # -> output/learned/card_glossary.tsv
python3 -m mirna_hallmark.learned.card_glossary --col beta # one column, every card it lives on
```
⚠ Still true, and the reason the glossary is a THIRD field rather than a rename: `domain` is
row-applicability, `rung` is the unit, `description` is the meaning. Three orthogonal facts.

## 3. ⭐ The labels are TESTED, not asserted

`card_rungs.verify(card)` checks each declared rung against the data's actual invariance **given that
card's key**:

| card key | gene | family | arm | edge |
|---|---|---|---|---|
| `[gene, arm]` | constant within gene | constant within (gene, seed_family) | constant within arm **across genes** | free |
| `[gene, family]` | constant within gene | **free — it IS the grain** | — | — |
| `[gene]` | one row per gene ⇒ **nothing checkable**; `agg_of` carries the meaning | | | |

**A declared rung the data contradicts is a mislabel** — that is how MH-191's wrong labels and the
`dose_rank_*` / `family_role` errors were caught. Run it:

```
.venv/bin/python3 -m mirna_hallmark.learned.card_rungs           # full report
.venv/bin/python3 -m mirna_hallmark.learned.card_rungs --check   # non-zero exit on any gap
```

## 4. Adding a column — the checklist

1. **Name the rung and the `agg_of`** before writing the column.
2. **Register it** in `CARDS[<card>]["explicit"]` (or a card-scoped prefix). Unregistered columns are
   flagged unassigned by `--check`.
3. **Run `--check`.** If the invariance fails, the label is wrong — fix the label, not the test.
4. **Both annotation passes run.** `realization.edge_card()`/`gene_card()` call `card_context.annotate()`
   **and** `card_ladders.annotate()` by default (MH-222/227). Writing a card without them silently drops
   ~57 columns.
   ⛔⛔ **AND THE ARM CARD IS WORSE, BECAUSE IT HAS NO SUCH DEFAULT — added 2026-08-19 after hitting it.**
   `python -m mirna_hallmark.learned.arm_card` does **NOT** call either annotator. `card_ladders.annotate()`
   attaches the arm block at its own call site (`_annotate(OUT/"arm_card.tsv", [(arm_admissibility_rollup(),
   ["arm"])])`), so a bare arm-card rebuild **silently drops `adm_n_edges` / `adm_n_with_site` /
   `adm_n_admissible` / `adm_frac_with_site`** — the arm rollup of `adm_has_site`, which MH-216 calls *the
   project's single most load-bearing conditioning variable*. It reappears only after:
   ```
   .venv/bin/python3 -m mirna_hallmark.learned.arm_card            # 287 -> 283 cols  ⚠ adm_* GONE
   .venv/bin/python3 -m mirna_hallmark.learned.card_ladders --annotate   # 283 -> 287, "bit-identical"
   ```
   **RULE: after ANY card rebuild — arm included — re-run `card_ladders --annotate`, then
   `card_rungs --check`, then `gen_cards --build`.** The annotator prints
   *"pre-existing columns bit-identical"*, which is the confirmation to look for.
5. ⚠ **`gene_family_card.tsv` has no `_finish_card` call site** — adding columns there means extending
   `readouts.py`'s promotion or adding one, **plus** registering in `CARDS["family"]["explicit"]`.

## 5. Traps with a recorded cost

- **`w_max` is the max curated EVIDENCE weight**, not a β or dose share (`top_beta_frac` /
  `concentration` are the shares). It is gene-rung on the gene card and family-rung on the family card.
- **The gene card is `realization/gene_card.tsv`.** `learned/gene_card.tsv` does not exist; looking there
  and concluding the atlas block is missing has already happened once.
- **`retention` names two unrelated quantities** — see `PATIENT_QUESTION_TAXONOMY.md` §5.
- **An arm-rung column must be constant within arm ACROSS genes.** If it varies by gene it is not an arm
  property, whatever its name suggests (MH-214).

## 6. Which rung is my question on? — full coverage

| rung | question taxonomy | status |
|---|---|---|
| edge | `EDGE_QUESTION_TAXONOMY.md` | ✅ + skill `apm-edge-question` |
| gene | `GENE_QUESTION_TAXONOMY.md` | ✅ + skill `apm-gene-question` |
| patient *(orthogonal axis)* | `PATIENT_QUESTION_TAXONOMY.md` | ✅ + skill `apm-patient-question` |
| **arm** | `ARM_QUESTION_TAXONOMY.md` | ✅ + skill `apm-arm-question` |
| **seed_family** | `ARM_QUESTION_TAXONOMY.md` §6 | ✅ covered as rung **C** there |
| family_edge | — | ⬜ none, and probably not needed |

**ARM / SEED_FAMILY — now covered** by `ARM_QUESTION_TAXONOMY.md`, as one doc rather than two: its
organising lemma is that **the arm is a measured but not an estimated unit** (`M` is fit per seed family
and only broadcast down), which partitions the rung into **A** arm-alone · **B** arm-in-family · **C** the
family itself. They were kept together because the `arm-in-family` resolvability machinery belongs to
neither pole alone.

**FAMILY_EDGE** — measured (MH-241) to carry *less* per-patient structure than the edge rung and to be
largely the edge question at a coarser grain. **Not obviously worth its own taxonomy.**
