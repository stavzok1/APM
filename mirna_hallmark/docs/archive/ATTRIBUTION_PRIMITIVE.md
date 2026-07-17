# Attribution primitive — the shared `exonerate → distribute → annotate` grammar

Every "who owns g" question in the learned model (identity, budget, supply-shift, discovery, the CN portion)
faces the **same** collinearity problem and — right now — solves it **ad-hoc or not at all**. This doc pins the
**one primitive** they should all call, so attribution is collinearity-honest everywhere at once. Status:
**active** (design). Companions: `COLLINEARITY_AND_IDENTIFIABILITY.md` (the two kinds), `CN_INSTRUMENT.md`
(where this was first worked out), `ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md` (identity vs magnitude).

---

## 0. The problem

Shapley/LMG **divides variance fairly but cannot exonerate a passenger**: a co-varying non-regulator B
(ρ≈0.8 with a real driver A) still gets ≈½·R²(B) from the *shared* variance. So any lens that runs raw Shapley
**over-credits collinear riders**. Today only `cluster_resolution` (retrofitted this arc) handles it; the core
lenses (`identity_vs_magnitude`, `budget_shift`, `state_pressure_attribution`) run raw Shapley/concentration —
their numbers currently credit passengers. One primitive fixes all of them.

## 1. The grammar — one primitive, two levels

```
attribute(g, units, level):                # units = families (level=between) or members of a family (level=within)
    1. EXONERATE  (between-family only)   conditioned partial anti-corr | (mates + C); rider (no unique signal) → 0
    2. DISTRIBUTE — TWO first-class shares (the identity / magnitude pair — NOT one number):
         (a) identity = Shapley/LMG of R²(X·M, Y)     # who OWNS the coupling  (variation, abundance-removed)
         (b) budget   = M·X̄   (X̄ = level abundance)    # who DELIVERS the dose  (level)
       [within-family, |ρ|≥0.95: only (b) is resolvable — variation can't separate, level still can]
    3. ANNOTATE   driver-class  +  the identity−budget GAP (loud-passenger / quiet-owner)  +  identity/dose reading
```

The steps are **separable and always in this order**: exoneration is a *gate* (real-regulator test), distribution
is the *quantitative split* — and it is **two shares, not one**: `identity` (Shapley, variation-based, who owns
the coupling) and `budget` (`M·X̄`, level-based, who delivers the dose). Abundance lives in `budget` as a
**first-class axis**, never a fallback (except literally at |ρ|≥0.95 where variation is exhausted). Their
**divergence is itself the signal** (the loud-passenger / quiet-owner gap, `ATTRIBUTION_IDENTITY_VS_MAGNITUDE`).
Only step 1 differs by level.

## 2. Between families — `exonerate → Shapley → class`

Families have **different seeds ⇒ different mechanisms ⇒ a passenger can be a genuine non-regulator** (a
co-transcribed family with no site). So:
- **exonerate:** conditioned partial anti-corr of each family | (other cluster-mates + C); a family whose *unique*
  residual variation doesn't survive (≥ −0.1) is a rider → **share 0** (this is the step Shapley cannot do).
- **distribute:** Shapley among survivors (`attribution.shapley_identity`) → collinearity-fair shares summing to
  the whole.
- **class:** `owner` (site ∧ couples) · `rider` (no site ∧ no coupling) · `designed-not-coupling` (site, no
  residual coupling) · `couples-no-site` (reference-blind — SNV/editing/APA candidate). *(from `cluster_resolution`.)*

## 3. Within a family — `Shapley-by-variation → class`, NO exoneration

Members share the seed ⇒ **identical sites ⇒ every member is a real regulator** ⇒ **never exonerate**. But they
are usually **resolvable** (real data: median member-abundance |ρ|=**0.44**, 62% of pairs <0.8, only 8% ≥0.95),
so:
- **distribute (two shares, both first-class):** `identity` = Shapley-by-variation (whose fluctuations carried
  the coupling) AND `budget` = level-abundance share (`2^a/Σ2^a` — who delivers the most dose). At |ρ|≥0.95 only
  the **level** share is resolvable (variation exhausted); everywhere else both are, and their gap is the
  within-family loud/quiet signal.
- **class:** `single-driver / co-drivers / family-only` — the conditioned partial anti-corr run *for annotation,
  not pruning* (`resolution_report`): which members *independently* drive vs which are passive co-contributors.
- Shapley itself never zeros a member (a low-variance member → small share, not 0); the conditioned partial's
  *binary pruning* is what's banned within a family, its *classing* is kept.

## 4. The interface (what to build)

One reusable function, called by every lens:

```
attribute(gene, level="between"|"within", *, on="coupling"|"supply") -> DataFrame
    index    = family (between) or member (within)
    columns  = identity  (Shapley/variation share, normalised)
             · budget    (M·X̄ / level-abundance share, normalised)     # first-class, not a fallback
             · gap       (identity − budget: loud-passenger / quiet-owner)
             · class     (owner/rider/… between · single/co/family-only within)
             · exonerated (bool; between only)
```

`on` picks what is divided (coupling R², or the supply-shift magnitude). The grammar is identical; only the
quantity fed in changes. Both `identity` and `budget` are always returned — attribution is the *pair*, never one.

## 5. Where each lens plugs in (the underexploitation)

| lens | today | with the primitive |
|---|---|---|
| `attribution.identity_vs_magnitude` | **raw Shapley, no exoneration** → over-credits collinear riders (live correctness issue on ESR1/PTEN) | `attribute(g, "between")` → exonerated riders drop to 0 before the split |
| `states.budget_shift` | within-family split by **abundance only** | add the within-family Shapley-by-variation + `single/co/family-only` class |
| `parallel_view.state_pressure_attribution` + `change_trajectory` | supply-shift (`dForce`/`dBud`) ranked by **raw concentration** → the global-shifter confound (miR-200 over-nominated) | exonerate against the other families' shifts *before* concentration |
| `program_network` | family **engagement / hub roll-up** sums M with no collinearity check → **co-transcribed hub families over-credited** | exonerate before the roll-up → hubs are real co-targeters, not co-transcribed riders |
| `card` (per-edge attribution card) | reports coupling/budget/identity per edge with **no class/exoneration flag** | stamp each edge with its `attribute()` class (owner/rider/…) + exonerated bool |
| `subtype` (`subtype_wiring`) + `state.cross_state_transfer` (`top_delta`) | per-family **ΔM reads** with no collinearity handling → collinear families' wiring-change confounded | run ΔM through `attribute()` so a rider's ΔM isn't read as rewiring |
| `discovery` | ad-hoc | the **sequence-prune = between-family exoneration** for the orphan set |
| `cluster_resolution` / CN portion | already retrofitted | *is* the reference implementation of the grammar |

**Out of scope (deliberately):** `structural_identity` — it's a *design/sequence* score (potential × specificity ×
confidence), expression-free, so the coupling-based exoneration doesn't apply (there's no coupling to condition
on). `hierarchical` — cross-*gene* pooling, an orthogonal axis. So the primitive is for the **coupling/dose**
lenses, not the design or cross-gene ones.

## 6. The class typology = the regulation fingerprint

Steps 1+3 produce, per (gene, family/member), a **class label**. Propagated genome-wide and joined onto the
card / discovery / registry, the cross-product of the between-class {owner/rider/designed-not-coupling/
couples-no-site} and within-class {single-driver/co-drivers/family-only} **is** the "regulation-pattern
fingerprint" parked in `LEARNED_MODEL_WHATS_NEXT §PARALLEL-VIEW` — this primitive is its engine.

## 7. Identity vs dose (the interpretation layer)

The same Shapley number reads differently by level, and the label must carry it:
- **between-family share** = attribution to a *distinct repressive mechanism* (different sites).
- **within-family share** = attribution to a *dose-deliverer of one mechanism* (identical sites) — never a claim
  that one member out-represses its mate.
So `attribute(..., level=)` must stamp the interpretation, not just the number.

## 8. Built vs the gap

**Built:** `shapley_identity` (the split), `resolution_report` (within-family class + conditioned coupling),
`cluster_resolution` (between-family exonerate→Shapley→class + CN).

**There is NO shared primitive today — each lens is inline and independent** (grep-verified): `shapley_identity`
is called in exactly ONE place (`identity_vs_magnitude`), the exoneration in ONE (`within_gene_ranks`), and every
other lens is **`shapley:0`** (rolls its own from budget/abundance/ΔM/roll-up). So a change does **not**
auto-propagate — there's nothing shared to change. The build is **factor once, then route each**:

1. **build `attribute()` once** — the grammar (§1) as one function; and
2. **route each consumer** through it — per-lens, because each divides a **different** quantity (`on=`: coupling
   R² / supply-shift / ΔM / roll-up). Cost per lens: **thin** (`identity_vs_magnitude` already Shapleys → just add
   the gate; `card` → surface the class/flag), **additive** (`budget_shift` already has the level share → gains
   Shapley+gate+class), **real** (`program_network`/`subtype`/`state` don't Shapley or exonerate at all → adopt
   the gate before their roll-up/ΔM). It is **no new statistics** (a unification of existing estimators), but it
   is **not free** and **not a reimplementation** — a one-time routing cost per lens.
3. **then it propagates** — after routing, changing the exoneration/Shapley logic hits every routed lens at once.

**Order:** do `identity_vs_magnitude` first — it's both the live-bug lens *and* the one already Shapleying, so
`attribute()` is born there with least friction and the rest become adapters against it.

## 9. Code map

| step | function |
|------|----------|
| exonerate (conditioned partial anti-corr) | `families.resolution_report`, `families.cluster_resolution` (cross-seed) |
| distribute (Shapley/LMG) | `attribution.shapley_identity` |
| annotate — within-family class | `families.resolution_report` (`single/co/family-only`) |
| annotate — between-family class | `families.cluster_resolution` (`owner/rider/designed-not-coupling/couples-no-site`) |
| abundance (level dose) fallback | `states.budget_shift` (`2^a/Σ2^a`) |
| **target primitive (to build)** | `attribution.attribute(gene, level, on)` |
| consumers to route through it | `attribution.identity_vs_magnitude`, `states.budget_shift`, `parallel_view.state_pressure_attribution` |
