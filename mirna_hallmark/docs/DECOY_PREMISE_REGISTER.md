# The decoy bench — premise register

> **Goal:** enumerate every ASSUMPTION the site-free decoy control rests on, and for each one name the
> test that checked it and its current verdict — so "is the decoy sound?" is a single read.
> **What belongs here:** the PREMISE LIST and each premise's status + pointer. **NOT** the numbers
> (they live in their `MH-##` row), **NOT** the construction recipe (that is the `eval/decoy_bench.py`
> docstring, which sits next to the code that implements it), **NOT** the headline verdict
> (→ `STATE_OF_PLAY.md` Axis 4), **NOT** open work (→ `PROGRAM_FORWARD_BOARD.md` §B).
> **Update trigger:** a new premise is discovered, or an existing premise's verdict changes.
> **Sync-partner:** none. This doc must stay a pointer table — **if you find yourself copying a number
> in here, stop: that is how MH-127/128/130 rotted.**

**Why this exists (2026-08-01).** Claim 1 of the program — *curated edges anti-correlate more than a
matched fake* — is carried entirely by this control, and its soundness was distributed across a module
docstring, ~12 registry rows, one STATE_OF_PLAY block and a board section. Nobody could answer *"which
assumptions has this design actually earned?"* without reassembling it. The archived `MH127_*` /
`MH128_*` docs tried to be that and instead became the doc-trap they are now listed as.

**How to read the verdict column.** **HOLDS** = tested and passed · **INERT** = the assumption is
violated but measured not to matter · **BOUND** = holds only as an upper/lower bound · **OPEN** = not
tested · **CLOSED-BY-DATA** = untestable with what is in this repo.

---

## A. What makes the FAKE a valid negative control

| # | Premise | Verdict | Where it was tested |
|---|---|---|---|
| A1 | A decoy arm **cannot repress** the gene — no context++ site, no scanMiR duplex, no Poisson-significant 6mer, no evidence-source hit | **HOLDS** | 4-layer rule in `decoy_bench._site_maps`; the 6mer layer was added after 14.6% of "site-free" decoys were found to carry one |
| A2 | **No site ⇒ no gap.** If the rule works, genes whose *real* edges are all seedless must show zero gap | **HOLDS** ⭐ | **MH-136** — the arc's strongest result and its only a-priori, sequence-only, non-circular positive control; correct on *both* sides |
| A3 | The evidence exclusion covers **all** evidence, not just high-confidence curation | **HOLDS** (was a 126× hole) | **MH-137**; Manakov chimeric is in the cache — verify with the `[sites]` log line, never by assumption |
| A4 | Remaining uncovered evidence would not overturn it | **CLOSED-BY-DATA** | POSTAR3 miRNA-target is **not in this repo** (only its lncRNA/RBP tables) ⇒ acquisition, not a re-run. Board §B |
| A5 | A decoy is not a disguised real regulator via family or genomic cluster | **HOLDS** | seed-family + cluster-mate exclusion in `build_decoys` |

## B. What makes the fake MATCHED (so the contrast is sites, not confounds)

| # | Premise | Verdict | Where it was tested |
|---|---|---|---|
| B1 | Matched on **signed** global loadings, not \|loading\| | **HOLDS** (was violated) | mPC2 is the purity axis; 39% of v4 decoys were anti-matched before the fix |
| B2 | Decorrelated from real regulators **in the sparse residual**, not raw | **HOLDS** (was violated) | a raw cap systematically rejects high-loading arms — the v3 error that ~2× overstated the gap |
| B3 | Matched on **dose** | **BOUND** — a residual survives and is removed post hoc | **MH-137** (caliper) / **MH-139** (caliper-off + `b·Δ`, b re-derived per decoy). ⛔ never import `b` |
| B4 | Matched on **variance / SD** | **HOLDS** — the mismatch runs the wrong way to manufacture the result | **MH-129(a)** |
| B5 | Matched on **design WIDTH** (families per design) | **VIOLATED but INERT** | **MH-167** — real 3.096 vs fake 3.510 families; gap contamination bounded near zero |
| B6 | Matched on **within-design collinearity** | **VIOLATED but INERT** | **MH-167** — same bound |
| B7 | Any single construction is an **upper bound** — every mismatch inflates | **BOUND, standing** | `decoy_bench` header; the trajectory −0.045 → −0.0306 → −0.0147 → −0.0129/−0.0124 has only ever shrunk |

## C. What makes the SCORING honest

| # | Premise | Verdict | Where it was tested |
|---|---|---|---|
| C1 | **Strictly out-of-fold** — β, z-params, variance floor and C's coefficients all train-only | **HOLDS** | `oof_budget` |
| C2 | Pooling fold predictions is **scale-safe** | **VIOLATED, FIXED, effect measured at ~4%** | **MH-181** — β is refit per fold, so pooling an arbitrary-scale index is not scale-invariant. Both arms shifted; the *gap* barely moved. ⚠ single-arm users were affected |
| C3 | The **internal null** (1 family, β uniform) must be ≈0 | **HOLDS** | **MH-137/139**; exact-zero form proven in `weight_gain` |
| C4 | Both **C blocks** are run and retention reported — never silently conditioned | **HOLDS** | axiom 2a; `_C_blocks` always runs both |
| C5 | The family collapse is a **true RPM pool**, never a mean | **HOLDS** | `pool_family` |

## D. What the result may and may not be used for

| # | Premise | Verdict | Where it was tested |
|---|---|---|---|
| D1 | The gap is a **SET-level** quantity, not per-edge significance | **STANDING CONSTRAINT** | **MH-123/124** — the per-edge null is 3–4× too narrow. ⛔ never call one gene "validated" from `ctx_gap_*` |
| D2 | Competence is a **gradient in design width**, not a partition | **HOLDS** | **MH-147**, reproduced **MH-169** |
| D3 | Width and abundance are **separable** axes | **FALSE** | **MH-169** — they are one axis here; the gap cannot be attributed to either |
| D4 | The evidence axis (`w_max`) adds something given width | **UNDERPOWERED, not "no"** | **MH-147** as recorded; softened by **MH-169** |
| D5 | An **abundance baseline is not a control** — benchmark against a *fitted* fake | **STANDING RULE** | **MH-115/127**; a fitted matched decoy also beats unfitted abundance |
| D6 | The decoy gap validates the **protein** transfer too | **FALSE** | **MH-130d** — a site-free fitted fake reaches ~99% of it; only the CPTAC-mRNA cell is decoy-controlled |

---

## The one-line summary

**Every premise that makes the fake a valid negative control (A) and honestly scored (C) currently HOLDS
or has been measured INERT; the matching premises (B) hold as an UPPER BOUND with dose handled post hoc;
and the live constraints are all in (D) — the result is set-level, width-and-abundance are inseparable,
and abundance is not a control.** The single unclosable gap is **A4 (POSTAR3)**.

⚠ **This table is a map, not a source.** Cite the `MH-##` row, never this file.
