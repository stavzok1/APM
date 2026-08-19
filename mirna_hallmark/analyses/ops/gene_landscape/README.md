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
