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

Same name, same estimator, different unit.

⭐⭐ **AND THE TWO `beta`s DIFFER EVEN WHERE THE COLLAPSE IS A NO-OP — measured 2026-08-19.** The natural
reading of the table above is *"they coincide for singleton families, since there is nothing to collapse"*.
**They do not.** On **SINGLETON** cells `beta_edge == beta_family` on only **21.6%** (median |diff| 0.0004,
max **0.2444**); on multi-arm cells, **0.2%** (median 0.0100, max **0.9680**). The reason is that the
family-grain fit **re-designs the WHOLE gene**: collapsing a gene's *other* families changes this family's
competitors, so its β moves even though the family itself was never collapsed. ⇒ **the rung difference is
not local to the collapsed cells — never mix the two cards' `beta`, not even for singleton families.**
✅ And the edge card's `beta` is genuinely EDGE-rung, not a broadcast: it varies across arms in
**467 of 467** multi-arm cells. A single shared prefix map cannot express that, which is why
each card carries its own map and why `domain_of(col, card)` is **card-scoped** — prefixes were once
matched globally and silently attached one card's caveat to another's columns (three bugs; `sub_` on the
arm card vs the edge card, fixed by renaming to `esub_`).

⛔⛔ **AND CARD-SCOPING A PREFIX IS NOT SUFFICIENT — A PREFIX CAN BE AMBIGUOUS *WITHIN ONE CARD*.**
*(Added 2026-08-19, column-review unit 12, after causing a 4th instance of the bug this section describes.)*
The edge card carries **two unrelated families under `arm_`**: the arm-resolution FIT quantities
(`arm_sep_z`, `arm_dbeta`, `arm_resolvable`, `arm_credit_share`, `arm_id_status`) **and** seven lifted
ABUNDANCE/TRAJECTORY columns (`arm_med_rpm`, `arm_pct_floor`, `arm_iqr`, `arm_lfc_*`). A global `arm_`
glossary block added during the rename unit silently re-described all seven as *"arm-resolved fit
quantities … 20.3% of edges"* — **wrong meaning, wrong domain, wrong rung, on 7 columns × 2 cards.**
⇒ **RULE: a prefix block is a DEFAULT, never a guarantee. After adding one, print what it now captures and
read the list.** Columns a prefix mis-serves get an exact `(card, column)` entry, which outranks every
prefix on every card.

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

### ⛔⛔ THE COST OF A BROADCAST COLUMN — measured 2026-08-19, and it is large

A column whose rung is **coarser than the card's key** repeats down the rows. Everyone knows that. What
was never quantified is what it does to a statistic computed over those rows:

> **A row-wise statistic over a broadcast column weights each parent unit by HOW MANY ROWS IT HAS.**

On the edge card that is the gene's **design width** — 5,649 rows over 1,420 genes, up to **91 edges on a
single gene**. And ceiling, dose, abundance and the decoy gap all TRACK width, so the weighting is not a
harmless imprecision: **it is the confound itself**. Same column, two rungs, two answers:

| column | edge-row median | gene-row median | |
|---|---|---|---|
| `ctx_ceiling` | **0.0809** | **0.0138** | **~6×** |
| `ctx_n_abund` | 4.00 | 1.00 | 4× |
| `ctx_dose_max` | 13.06 | 8.97 | |
| `ctx_gap_deconv` | −0.0196 | −0.0085 | 2.3× |
| `ctx_gap_core` | −0.0186 | −0.0148 | |

**RULE: compute a broadcast column's statistic on the card that OWNS its rung, or de-duplicate to one row
per parent first.** Two of the columns above are load-bearing — the decoy gap and the measurability
ceiling — and a 6× shift in the ceiling would silently move any claim conditioned on measurability.
*(Checked: the MH-248 decoy stratification used the gene card and reproduces exactly at −0.0148; on edge
rows it would have read −0.0186.)*

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

## 2b. ⭐ WHEN IS A BLOCK "THE SAME BLOCK AT ANOTHER RUNG"? — THE NO-OP TEST

*(Added 2026-08-19, user-asked: "you say cptac can be inherited from edge to gene? in what way?")*

Two blocks sharing a prefix on two cards can stand in three different relations, and **the review shortcut
of closing one by inheritance is only valid for the first**:

| relation | test | example |
|---|---|---|
| **SUBSET** — same columns, same rungs | name-set inclusion + rung equality | `ctx_` on gene/gene_family ⊂ edge (units 6/7, closed by inheritance) |
| **AGGREGATE** — same estimand, coarser unit | ⭐ **the NO-OP TEST below** | `cptac_` edge → gene |
| **UNRELATED** — one prefix, two quantities | read the values | `arm_` on the edge card (unit 12) |

⛔ **`cptac_` shares ZERO column names between the edge and gene cards, so it cannot be closed by
inheritance — but it is NOT unrelated either.** The gene card inserts an infix naming the aggregation:

```
edge   cptac_{cohort}_rho_{layer}[_raw]        ← ONE arm vs the gene's protein/RNA
gene   cptac_{cohort}_agg_rho_{layer}[_raw]    ← the β-WEIGHTED SUM over the gene's arms
gene   cptac_{cohort}_abund_rho_{layer}        ← the UNWEIGHTED abundance-sum reference
```

⇒ **the gene value is NOT derivable from the edge values**, because `ρ(Σβx, y) ≠ mean ρ(xᵢ, y)`. Measured:
spearman **+0.870** (prospective) / **+0.782** (t105) against the per-gene mean, and the aggregate is
**~40% STRONGER** (median |ρ| 0.0838 vs 0.0588) — which is the point of aggregating, not a discrepancy.

### ⭐ THE NO-OP TEST — the check that makes an aggregate relation falsifiable

> **On units where the aggregation is a NO-OP, the two rungs MUST agree exactly.** For `cptac_` that is
> genes with exactly one arm: a weighted sum of one term IS that term.

It is the only cheap check that distinguishes *"a coarser view of the same estimand"* from *"a different
quantity that happens to share a prefix"* — and it immediately found a defect:

| cohort | true no-op genes | disagree | median \|diff\| | max |
|---|---|---|---|---|
| **prospective** | 466 | **0 (0.0%)** | — | — |
| **tcga105** | 389 | **112 (28.8%)** | 0.0587 | **0.6446** |

⛔⛔ **The prospective cohort proves the identity SHOULD hold exactly. TCGA-105 breaks it on 29% of genes
where the aggregate is mathematically the single edge, by up to 0.64.** Ruled out: it is not a sample-count
difference (edge-side n = 105 on both the matching and mismatching sets) and not a mis-specified no-op set
(all 112 have `n_arms == 1` *and* exactly one scored edge). ⬜ **UNDIAGNOSED — the cause is somewhere in
`cptac_card`'s t105 path.** Until it is resolved, **do not cross-reference edge-rung and gene-rung
`cptac_*` values on the tcga105 cohort.** (MH-262, `PROGRAM_FORWARD_BOARD.md`.)

**RULE: whenever you add a block at a second rung, write down its no-op set and assert the identity there.**
An aggregate that cannot reproduce itself on a one-element set is not an aggregate.

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
5b. ⛔⛔ **A PRUNE MUST RUN AFTER THE LAST THING THAT CAN RE-INTRODUCE THE COLUMN — and "it printed
   DROPPED" is not evidence it is gone.** *(Added 2026-08-19, unit 25, after two silent failures.)*
   `_annotate` **preserves every pre-existing column** by design. So a prune placed in a normaliser (which
   runs FIRST) is undone for any column that lives on the BASE card: the run prints `✅ DROPPED`, the
   annotator re-reads the file, and the column is still delivered. `adm_expressed` survived two such
   attempts before the prune was moved to a **post-annotation** step.
   **RULE — a column can enter the card from three places, and a prune has to cover the one it uses:**
   | where it comes from | where to prune |
   |---|---|
   | an annotation block (`card_ladders`) | the block's own `keep`/`return` |
   | the BASE card (`canonical_card` / `realization`) | the **post-annotation** step, and the base builder |
   | a normaliser-owned rename | the normaliser |
   ⇒ **verify by re-reading the delivered file, never by trusting the print.** Both are now in place for
   `adm_expressed`: `edge_admissibility()` no longer emits it, and the post-prune covers the stale base
   until the rebuild.

5. ⚠ **`gene_family_card.tsv` has no `_finish_card` call site** — ⛔⛔ **AND THIS BIT, EXACTLY AS
   WRITTEN (2026-08-19).** Column-review unit C pruned `p_fam` and renamed `n`→`n_samples`; that landed via
   `realization._normalise_edge_names`, which funnels the **edge card only**, so the family card kept BOTH
   for two weeks while the change was recorded as applied. ⇒ the card now has a funnel of its own,
   `card_ladders.normalise_family_card()`. **A schema change is not applied until you have checked it on
   every card that carries the column** (axiom 2, downstream ripple).
   The original note stands too: — adding columns there means extending
   `readouts.py`'s promotion or adding one, **plus** registering in `CARDS["family"]["explicit"]`.

## 5. Traps with a recorded cost

- **`w_max` is the max curated EVIDENCE weight**, not a β or dose share (`top_beta_frac` /
  `concentration` are the shares). It is gene-rung on the gene card and family-rung on the family card.
- **The gene card is `realization/gene_card.tsv`.** `learned/gene_card.tsv` does not exist; looking there
  and concluding the atlas block is missing has already happened once.
- ⛔⛔ **GATING THE INPUTS DOES NOT GATE A DERIVED COLUMN — cost: two blow-ups, one found 11 units after
  the other.** `identity` is a SIGNED share, so any max/sum over it is unbounded. The `identity_reliable`
  gate shipped to the edge and gene_family cards in the unit-1/2 column review; **11 units later the gene
  card's `top_identity` still read +740.007 on 80 genes**, because `readouts` had already computed
  `float(d.identity.max())` from the pre-gate values and nothing recomputed it. Same shape as
  `arb_max_identity` on the arm card. ⇒ **when you add a gate, grep for every column DERIVED from the
  gated one and recompute it** — `top_identity_gated` / `arb_max_identity` now do. And ship the gated
  version BESIDE the raw column, never over it: silently rewriting a value another module wrote is how
  provenance rots.
- ⛔⛔ **`retention` names FIVE unrelated quantities, not two.** MH-258 found four estimands sharing
  the name (`β_deconv/β_core` · `ρ_adj/ρ_raw` · `gap_deconv/gap_core` · `state.py`'s `ρ_H/ρ_Tsub`), and
  unit 25 found a fifth: **`field_retention` is not a ratio at all** — it equals `field_r_own −
  field_r_perm`, a permutation-corrected EXCESS (verified exact within rounding on 571 arms). ⚠ Because
  it is bounded like a correlation it never trips `ratio_blowup_audit`, so **the tooling cannot find this
  class — only reading the estimator can.** `learned/retention.py` is the home for the four true ratios;
  this one is deliberately NOT in it.
- **`retention` (original note) names two unrelated quantities** — see `PATIENT_QUESTION_TAXONOMY.md` §5.
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
