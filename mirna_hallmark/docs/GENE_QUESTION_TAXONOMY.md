# Gene-level question taxonomy — the gene's total incoming miRNA regulation

> **Scope: the GENE.** A target gene `g`, its **aggregate incoming pressure**
> `P_agg(g,s) = Σ_{m∈R(g)} c(m,g,s)`, vs its own expression `y_g`, across samples. Question:
> *"is `g` net-repressed by its miRNA stack, and which construction best detects it?"* Companion to
> `EDGE_QUESTION_TAXONOMY.md` (single edges) and one rung below the **program** level
> (`hallmark_interaction`, out of scope here). Numbers/claims live in modules + `DISCOVERY_REGISTRY.md`.

## 1. The aggregation lemma — *the edge lemma flips here*

At the edge, a single-arm rank coupling is **invariant** to per-arm-monotone transforms (reference/
scale inert; only the softmax bites; abundance wins — MH-78). **At the gene aggregate the opposite
holds**, because summing breaks monotonicity: `Σ_m f(x_m)` is **not** a monotone transform of `Σ_m x_m`.

Three consequences:
1. **Per-arm construction is consequential** — z-centering/referencing each arm before summing
   reweights which arms dominate the sum, changing `P_agg`'s sample-ranks and thus the coupling.
   Empirically `softmax_z`/`softmax_z_logrpm` (411/407 neg-sig) **beat naive `abundance_sum` (398)**.
2. **Bare share-sum is degenerate** — `Σ_{m∈R(g)} share_m ≡ 1` (constant per sample) → ~no coupling
   (143 neg-sig). A share predictor must carry a per-arm **magnitude/dynamics** term to be non-constant.
3. **Aggregation function + coherence matter** — sum vs mean, and signed vs positive-part vs abs
   (a canceling/incoherent stack is not real net pressure).

**So the gene level is where construction can genuinely improve detection** — the current best beats the
naive baseline by only ~2%, leaving headroom. The grid (§ below) probes it.

## 2. Axis Group I — per-arm construction (now consequential)

`c(m,s) = comp(transform(x_m; ref, scale); promisc, aff)` — same Axis-Group-I as the edge
(`pressure_engine.ExprMode`): reference ∈ {cohort, healthy(GTEx), NAT}, scale ∈ {level, z, ratio},
comp ∈ {standalone, softmax-share}, promiscuity, affinity. Unlike the edge these move the coupling.

## 3. Axis Group II — aggregation

| Axis | Values | Engine knob |
|---|---|---|
| **function** | sum · mean · max · top-k · share-sum (degenerate) · mass-action-occupancy-sum | `compute_gene_pressure(aggregate=…)` |
| **contrib transform** | signed · positive-part · abs | `contrib_transform=…` (coherence) |
| **promiscuity** | none · evidence_mass · degree · ts/combined | `target_norm=…` |
| **prune method** | all-HE · evidence-tertile · evidence-decile · abundance-floor · top-k/gene | filter `R(g)` before aggregating |

## 4. Axis Group III — estimator

| Axis | Values | Note |
|---|---|---|
| confounders Z₀ | CPE · HRD · target-CN · proliferation · composition | as edge |
| gene-own conditioning | target-CN · methylation · ATAC | the gene's *own* genomic state — analog of the edge `P_others` (decoupling does this) |
| **AGO gating** | ungated · full · partial-λ · purified · ±TNRC6 co-limit | per-sample multiplier on `P_agg` (`ago_gate`); the NATURAL home of the gate (vs the edge where it injected a proliferation confound). Purified = residualize AGO capacity on proliferation+CPE+HRD before gating |
| state / reference | cohort · subtype · acquired (healthy/NAT) | as edge |

## 5. Axis Group IV — object

target layer (mRNA · protein L2/L1b · isoform); gene **role** (TSG/oncogene + `malignancy_sign`,
`gene_roles`); stack **coherence** (`global_signed_share` vs `global_abs_share` — canceling ⇒ |signed|≪abs).

## 6. Named questions (G-series)

| ID | Question | Estimator | Status / home |
|---|---|---|---|
| G1 | Is `g` net-repressed by its aggregate pressure? | `partial-ρ(P_agg, y_g \| Z)` | **built** — `mirna_comovement.gene_corepression` (282/1424 tumor) |
| **G2** | Which **aggregation** (function/transform/promiscuity/prune) best detects net repression? | construction grid vs `abundance_sum` | **run (MH-79)**: function barely matters (linear≈mass-action≈temp≈sparsemax≈cohort-softmax); **AGO-gating HURTS monotonically** (ungated best); **pruning neutral-to-harmful** (all best; evidence-pruning collapses on sparsity) |
| **G3** | Which **per-arm construction** (ref×scale) maximizes it? | construction grid | **run (MH-79)**: **NAT-(acquired-)referenced softmax share = 448 vs abundance-sum 398 (+12.6%, 96 net-new genes)** — the **reference** is the winning axis (NAT>cohort>rawx>healthy); z-sum ≈ baseline |
| G4 | Which arm **owns** the gene's pressure? | `gene_struct_share` | built (⊥ coupling, §5a) |
| G5 | Net of the gene's **own genomic state**? | +target-CN/meth/ATAC | built (decoupling) |
| G6 | Realization onto **protein/isoform**? | L2/L1b | built (CPTAC) |
| G7 | Differs by **state/subtype**? | cohort/subtype/acquired | built (pam50_gene_resolution, cross_state) |
| **G8** | Stack **coherent**? does **role** flip the meaning? | `global_signed_share` + role overlay | **run (MH-79)**: net-repression concentrates on known cancer genes (TSG 52%, onco 42% vs dual/unknown 27%, ~2× bg); coherence median 0.62; **~35% of net-repressed genes have an incoherent (canceling) stack** |
| **G8b** | Does **AGO gating** the aggregate help (and which gating)? | gated vs ungated `P_agg` | **run (MH-79): NO** — gating hurts monotonically (proliferation confound generalizes from edge MH-78); purified mitigates vs full but < ungated |
| G9 | Program/gene-set aggregate (next rung) | per-Hallmark pressure vs mean target | built (`hallmark_interaction`) |
| G10 | **Acquired/lost** net repression | escape_scope | built (decoupling §14a) |

## 7. The construction/prune/gating grid (G2+G3+G8+G8b — `gene_construction_grid`)

`coupling_predictor_comparison --gene-grid`. Sweep, scored as `partial-ρ(P_agg, y_g | CPE+HRD+target-CN)`
per gene (BH-FDR per config, neg-sig = ρ<0 & q<FDR), vs the `abundance_sum` baseline:
**construction** (expr_mode) × **prune** (all/evidence-tier/abundance-floor/top-k) × **AGO gating**
(ungated/full/partial/purified) × aggregation (sum/mean) × contrib (signed/pos/abs) × promiscuity.
Gating is a cheap post-multiply on `P_agg`, so pressure is computed once per (construction, prune) and
each gate applied. **Role/coherence overlay**: net-repression rate split by `malignancy_sign`; per-gene
coherence from `compute_gene_pressure_contributions`.

**Open questions the grid answers:** (a) can any construction beat `abundance_sum` by **> the current
~2%**, and which axis drives it (z-centering / healthy ref / coherence / mean-vs-sum)? (b) does
**pruning** to high-confidence/expressed regulators sharpen the aggregate? (c) does **AGO gating** the
aggregate help here (unlike the edge), and which gating (purified expected best)? (d) does net
repression concentrate on **TSGs** (pro-tumor), and on **coherent** (non-canceling) stacks?

Outputs: `output/coupling_predictor_comparison/gene_construction_grid_{comparison.tsv.gz,summary.tsv}`,
`gene_construction_role_coherence.tsv`.
