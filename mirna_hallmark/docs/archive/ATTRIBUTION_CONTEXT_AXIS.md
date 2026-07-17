# Attribution-context axis — forward design (documented 2026-07-05, build later)

Parked per user (#3, 2026-07-05): the bridge from **edge-discovery / coupling** (built) to an
**annotation taxonomy** — for every anti-coupling, *what regime is it in, how much of the gene it owns,
and how does that change between states*. Grounded in the precursor conceptual files
**`EDGE_QUESTION_TAXONOMY.md` (E1–E10)** and **`GENE_QUESTION_TAXONOMY.md` (G1–G10)** and
`ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md` (share = identity ⊥ magnitude, FORMULAS §5a). This doc adds the
**annotation/context layer** those taxonomies name but don't yet render per-edge.

**What the learned model already built of the precursor's "open" estimands** (so the annotation has data):
- **E5** (paired-delta `ρ(Δx,Δy)`, same-patient NAT — "acquired change tracks acquired change") = now
  `states.realization` (was "open"). 15/15 hub genes realized, mean ρ=−0.30 (2026-07-05).
- **E6 / G7** (where the coupling lives — subtype/state) = `states.realization_by_subtype`, `cross_state_coupling`.
- **E7 / G4** (share of the gene's budget) = the identity axis; the "place in the regulatory budget" the
  user means.
- composition Z (E-Group-II confounder) = the **deconv-retention tag** (#0): miR-29→collagen 6–44% (stroma),
  miR-135b→GATA3 86% (cell-intrinsic).

---

## 1. Per-edge annotation card (the immediate add)

Each anti-coupling (m→g) carries, beyond `ρ`:

**(a) Expression-regime context** — so a coupling is read in the right regime:
- arm: RPM range, % samples above the functional floor (RPM≥10), median in tumour / NAT / GTEx.
- gene: expression range + variance; whether the anti-coupling sits in the arm's high- or low-abundance
  regime (occupancy/sub-saturation context — coupling is cross-sample variance, and needs *both* sides to vary).

**(b) Budget share (E7/G4)** — the arm's share of the gene's total realized pressure among its co-regulators
(identity, abundance-removed): concentrated vs diffuse, the arm's rank, promiscuity/co-expression split.

**(c) Composition tag (mandatory, #0)** — deconv-retention = `ρ_deconv/ρ_raw`: **cell-intrinsic (≥0.7) /
partial (0.4–0.7) / composition-explained (<0.4)**. A "composition-explained decoupling" is **not discarded**
— it is a real *tissue-architecture* signal (stroma fraction co-varies), tagged as such and kept as its own
class (user: "a composition-explained decoupling is possibly interesting as well").

**(d) Shift-status across states — TWO orthogonal axes, and the canonical taxonomy already exists:
`mirna_state_class.joint_edge_class`.** Do NOT reinvent it. It joins:

- **LEVEL / rank axis (potential — "how much, vs which reference"):** the miRNA's cross-state pressure
  trajectory via within-state percentile **rank** deltas (GTEx is cross-platform → rank/QN, not raw): `dHN`
  (healthy→NAT = field effect), `dNT` (NAT→tumour = tumour-specific increment), `dHT` (healthy→tumour = the
  **acquired** headline; healthy-anchored QN-logFC). A *gainer* has dHT>0 — the user's "actually represses in
  tumour AND is much higher than healthy". Also the arm's **rank-within-the-gene's-budget** shift (E7/G4 across
  states) — the same level axis at the attribution level.
- **REALIZATION / coupling axis (realised):** composition+proliferation-adjusted partial-Spearman **per
  state**, on two ladders (abundance→target; realized-pressure→target) — does the edge actually couple in
  that state (has variance AND anti-correlates).

Their **join = `joint_edge_class` ∈ {acquired_realized, field_established_realized, acquired_unrealized,
constitutive, lost, nat_decoupled, stable, non_monotonic, uncoupled}** — for both the healthy (dHT) and NAT
(dHN/dNT) references. `acquired_realized` = the case the user named (elevated over healthy **and** realized).

**What the learned model ADDS to `joint_edge_class` (the reason to extend, not replace):**
1. **Paired-delta realization (E5, `states.realization`)** — within-*patient* (Δx vs Δy, same person),
   strictly stronger than the precursor's per-state *cohort* partial-Spearman: it removes the patient
   baseline, so the realization axis is baseline-free. Feed it as a second, sharper realization ladder.
2. **composition-explained as an explicit KEPT class (not just an adjustment).** The precursor *adjusts for*
   composition inside the coupling axis; the user wants a realized edge whose coupling **vanishes under
   deconv** tagged `composition_explained` and retained as a tissue-architecture signal (deconv-retention tag
   c). So: `realized` splits into `realized_cell_intrinsic` (retention≥0.7) vs `realized_composition` (<0.4).
3. **low-variance / undetectable** as an explicit detectability state distinct from `uncoupled`: the arm OR
   gene has ~no within-state variance → **no coupling is measurable ≠ no regulation** (range-restriction
   honesty guard). The precursor's `uncoupled` conflates "tested, flat" with "couldn't test"; split them.
4. the **learned cell-intrinsic M** (deconv) as the edge weight, vs the precursor's static curated `w(m,g)`.

---

## 2. Per-gene aggregate card (the same picture, one rung up — G-series)

The gene's *whole* incoming-regulation situation and how it shifts (mirrors §1 at `P_agg`, the G-taxonomy):
- **net-repression status (G1)**: is `g` net-repressed by its stack? `partial-ρ(P_agg, y_g | Z)`.
- **stack coherence + role (G8)**: signed vs abs share — a canceling stack (~35% of net-repressed genes)
  is not real net pressure; overlay TSG/oncogene role (`gene_roles`).
- **budget concentration (G4)**: is the pressure held by one arm or spread (Herfindahl of the shares).
- **composition fraction of the aggregate**: how much of `P_agg`'s coupling is stroma vs cell-intrinsic
  (aggregate analog of tag c).
- **shift class NAT→tumour (G10, per subtype G7)**: the gene-level version of §1(d) — does the gene
  *acquire* / *lose* net repression, and is its **budget redistributed** among arms (the attribution
  hand-off / recruitment taxonomy) or does it just scale with abundance? Uses `states.realization` (E5) at
  the aggregate + the per-subtype split (E6/G7).

---

## 3. Build order (when resumed; infra now mostly present)
1. `w_eff` (realized, per-sample, floor/occupancy-aware) → shares (E7/G4). **Build on the deconv /
   cell-intrinsic M only** — never the stroma-contaminated M, or shares are composition shadows.
2. Per-edge card = join {coupling, deconv-retention, expr-range, share, shift-class}. Table, one row/edge.
3. Per-gene card = the G-series roll-up + shift class.
4. Validate a known **hand-off** (a gene whose dominant regulator changes NAT→tumour or by subtype) before
   scaling — the analog of validating miR-135b→GATA3 for the discovery lane.

**First concrete deliverable:** a per-gene "attribution card" — regulators ranked by share, each with its
expression-range + composition tag + NAT→tumour shift-class — starting from the cell-intrinsic identified
edges (`programs --deconv` output) and the paired Δ (`states`).
