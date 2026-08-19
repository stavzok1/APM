# gene_landscape — the four gene-rung characterisations (2026-08-18)

Producers for the dossier's **Biology** section. Run from the repo root with
`PYTHONPATH=. .venv/bin/python3 <script>`; they read the cards and print, they do not write cards.

| script | what it answers |
|---|---|
| `axes_build.py` | Builds `output/learned/gene_axes_matrix.tsv` (52 per-gene axes: `card_*`, `ident_*`, `reg_*`, arm-identifiability rollups) and runs an unconditioned FDR scan of 8 outcomes. |
| `conditional.py` | ⭐ **The one that matters.** Re-scans every outcome PARTIALLED on `n_fam` + `n_arms`. Unconditioned, three of the outcomes are mechanically monotone in design width and every scan just rediscovers `n_fam`. |
| `three.py` | Arm-vs-family variance capture (`oof_drho` by cell size), the CPTAC cohort inventory, and the first cut at handoff. |
| `handoff.py` | Splits `regulatory_handoff` into its two legs — it fires on HLY≠TUM **or** NAT≠TUM, and 27 of 328 genes fire only via the NAT leg and are not healthy→tumour switches. |

⚠ **These are DESCRIPTIVE.** No decoy, no null beyond the FDR, conditioning on design width only.
They are registry *candidates*, not registry rows — send them through `rigor-auditor` before promoting.

⚠ `axes_build.py` assembles the axes by hand because **`gene_axes.build_axes()` does not exist** — the
module docstring's usage example names it, but the module ships the pieces (`regulator_axes`, `self_axes`,
`weight_axes`, `hhi`, `scan`, `contrast`, `mask_degenerate`) and expects the caller to assemble. Deviation
recorded: the `reg_var_*` axes use the card's `arm_iqr` as the dynamic-range measure where `regulator_axes`
would use a per-sample SD.

⚠ A BH bug lived in `conditional.py` for one run: `cummin()` on a p-ascending series propagates the
minimum from the top, which reported **36/36 axes significant for every outcome**. BH step-up needs the
REVERSE cumulative minimum. Fixed; the corrected counts are 2–14 of 36.

## Round 2 (2026-08-19) — strata, extremes, discovery, ladder, audit

| script | what it answers |
|---|---|
| `decoy_strat.py` | The decoy gene STRATIFICATION (competence class / width / measurability) and both TAILS with named genes. The internal null (`n_fam==1`) lands at −0.0024, p=0.49. |
| `discovery.py` | The 157-edge gold set's real composition (11 families, 61% one family) **and the convergent evidence ladder** — seedless → site → +1 chimeric → +2 chimeric, ρ=−0.241 (p=5.4e−55), abundance-matched 5/5 quintiles. |
| `extremes.py` | Named extremes for attribution, arm-in-family, protein and state. ⚠ Its attribution ranking is SUPERSEDED by `attr_gated.py` — it ranked on an ungated `top_identity`. |
| `topident.py` | ⛔ Diagnoses the `top_identity` blow-up: edge `identity` is a SIGNED share summing to 1, 10.5% negative (suppressor effects), so `max()` is unbounded — 740.0 at worst, concentrated on unmeasurable genes. |
| `attr_gated.py` | The GATED attribution disagreement list (`top_identity ≤ 1` ∧ `ceiling > 0.02`): 19.8% vs 31.1% on unmeasurable genes. |
| `magn_vs_ident.py` | Do magnitude and identity pick different KINDS of family? Paired within gene: magnitude's pick has 17 PMIDs vs identity's 73 (p=3.1e−12). Control on agreeing genes fires correctly. |
| `col_audit.py` | Block-by-block column audit: 23 constant, 9 near-empty, 15 duplicate pairs, 0 sparse-without-domain. |

### Defects these found
1. **`p_fam` is NOT a p-value** — it is the design dimension (count of family predictors), bit-identical to `n_fam` on the gene_family card.
2. **`arb_n_edges` == `arb_n_genes` == `arb_n_identity_reliable`** — one vector under three names.
3. **`top_identity` is unbounded** and needs the `add_reliability` gate treatment.
4. **Five `fame_assay_perturbation*` columns are all-zero**; the `fame_` block is 43 of the arm card's 297 columns.
5. **`echim_any` is never False** — only True or NaN; test `!= 1`, never `== 0`.
