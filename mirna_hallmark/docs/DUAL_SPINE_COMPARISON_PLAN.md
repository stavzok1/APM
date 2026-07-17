# Dual-spine comparison plan — miRTar M0 vs ENCORI M0′

**Status:** planned (2026-06-22)  
**Prerequisite:** `mirna_hallmark/output/encori/encori_mirna_mrna_targets.parquet` (21,346 rows; 1,951 genes with hits)

**Goal:** Run two edge backbones in parallel with the **same pressure engine** and **same coupling benchmark**, then decide empirically whether miRTar-first or ENCORI-first is the better headline spine.

---

## 1. Research questions

| ID | Question |
|----|----------|
| R1 | Does **S1** (miRTar + ENCORI boost) hold **Basal ≥ 42/50, key 8/8** while improving hub attribution? |
| R2 | Does **S2** (ENCORI M0′ + miRTar boost) tell a **more coherent story on cliff genes** (IRF1) without destroying global coupling? |
| R3 | Where do the spines **disagree** (edge set, rank, pressure) — hub panel + disagreement export? |
| R4 | Is **ENCORI depth in share** (Lane 2) worth testing *after* edge-backbone winner is clear? |

**Not in scope for v1:** changing registered spine, DISCOVERY_REGISTRY, or default `run_all` until review of this package.

---

## 2. Shared pressure layer (fixed across spines)

Same as current L1 coupling benchmark:

```
expr_mode     = softmax_z_logrpm
target_norm   = evidence_mass
aggregate     = sum
ago_gate      = default (TNRC6 off)
min_abundance = config PRESSURE_ABUNDANCE_FLOOR
```

Per-edge contribution (unchanged):

```
c(m,g,s) = share(m,g,s) × z(m,s) × log2(RPM+1) × edge_w(m,g)
edge_w   = log1p(evidence_score_eff) / log1p(Σ_g log1p(evidence_score_eff))
```

Only **`evidence_score_eff`** and the **edge table** differ between S1 and S2.

---

## 3. Spine definitions

### S0 — Baseline (reference)

| Field | Value |
|-------|--------|
| Edge set | **M0** — tiered permissive, `min_evidence ≥ 2`, arm-resolved |
| `evidence_score_eff` | tiered permissive (no ENCORI) |
| Pairs | **5,315** (current) |
| Role | Reproduce registered benchmark |

### S1 — miRTar backbone + ENCORI depth boost (Lane 1)

| Field | Value |
|-------|--------|
| Edge set | **Same as M0** (fixed pairs) |
| Base score | `confidence_logclass` (preferred) or `tiered_permissive` if grid shows tie |
| ENCORI modifier | only where `(m,g)` exists in ENCORI parquet (arm-resolved) |

```python
enc_depth(m,g) = 2.0·log1p(clipExpNum)
               + 1.0·log1p(degraExpNum)
               + 0.5·I(TargetScan ∧ miRanda ∧ PITA)
               + 0.25·log1p(pancancerNum)

log_evidence = log1p(base_score(m,g)) + α · enc_depth(m,g)   # α=0 if no ENCORI row
evidence_score_eff = exp(log_evidence) - 1   # or store log_evidence directly in engine
```

**Sensitivity grid (S1 only):** `α ∈ {0, 0.25, 0.5, 1.0}` on fixed M0 pairs.

**Expected:** ~25% of M0 pairs get non-zero boost (1,352 / 5,315); coupling stable; PTEN/IRF1 rank-ρ may improve.

### S2 — ENCORI backbone M0′ + miRTar functional boost

| Field | Value |
|-------|--------|
| Edge set | **Built from ENCORI**, not miRTar M0 |
| Gate | `clipExpNum ≥ 2` (primary); sensitivity `≥ 1` as row |
| Arm filter | `hsa-` only; resolve to expression matrix |
| Abundance | `cohort median logRPM ≥ PRESSURE_ABUNDANCE_FLOOR` |
| Per-gene cap | top **K** arms by `enc_depth` (default **K=25**; hub genes **K=50** via `HUB_ROUTES`) |
| Base score | `enc_depth` (raw, pre-log) |

**miRTar modifier** (inverse of S1 — functional validation on CLIP rows):

```python
mirtar_boost(m,g) = log1p(tiered_permissive(m,g))   # 0 if pair not in summary
evidence_score_eff = enc_depth(m,g) + β · mirtar_boost(m,g)
```

**Sensitivity grid (S2 only):** `clipExpNum_min ∈ {1, 2, 3}`, `β ∈ {0, 0.5, 1.0}`, `K ∈ {15, 25, 50}`.

**Expected:** denser on IRF1-like genes; sparser on TP53-like; **different** total edge count (~8–15k before cap); coupling **unknown** — must measure.

---

## 4. ENCORI table prep (one-time)

Module: `mirna_hallmark/encori_edges.py`

1. Load parquet `output/encori/encori_mirna_mrna_targets.parquet`
2. Collapse to **one row per (miRNAname, geneName)** — keep row with max `clipExpNum` (tie-break `degraExpNum`, then `pancancerNum`)
3. Harmonize `miRNAname` → expression matrix arms via `resolve_edges_mirna`
4. Compute `enc_depth` column (formula above)
5. Optional flags: `program_consensus`, `has_degra`, `cellline_breast` (descriptive; not in v1 gate)

Outputs:

- `output/encori/encori_pair_table.tsv.gz` — collapsed pairs + depth
- `output/encori/encori_arm_resolve_audit.tsv` — dropped arms

---

## 5. Comparison module

Module: `mirna_hallmark/dual_spine_comparison.py`  
CLI: `.venv/bin/python3 -m mirna_hallmark.dual_spine_comparison`

### 5.1 Runs

| run_id | spine | α / β / clip / K |
|--------|-------|------------------|
| S0_baseline | M0, no ENCORI | — |
| S1_a0 … S1_a1 | M0 + ENCORI boost | α ∈ {0, 0.25, 0.5, 1.0} |
| S2_grid | ENCORI M0′ + miRTar boost | small factorial (see §3) |

Reuse: `hallmark_coupling_by_subtype`, `compute_gene_pressure`, `compute_gene_pressure_contributions`, hub rank-ρ from `evidence_scoring_sensitivity`.

### 5.2 Metrics (per run)

**Coupling (primary):**

- Basal `n_neg_sig / 50`, `n_key_neg_sig / 8`
- Median `rho_prolif_cn_wsd_adj` (Basal)
- Key hallmark pivot (same as `pressure_evidence_sensitivity`)

**Attribution (secondary):**

- Median within-gene `edge_w` max/min ratio
- Hub panel: PTEN, CDKN1A, E2F1, IRF1, TP53, RB1, BRCA1
  - `spearman(edge_w, mean_abs_contribution)`
  - top-5 arms by structural vs realized rank
- IRF1-specific: list M0 vs M0′ arms, ENCORI-only arms

**Coverage:**

- `n_edges`, `n_genes`, `n_pairs_with_encori` (S1), `n_pairs_with_mirtar` (S2)
- Jaccard edge sets: S0 vs S2 best

### 5.3 Outputs

```
mirna_hallmark/output/dual_spine_comparison/
  method_manifest.json
  spine_spec_manifest.tsv
  coupling_summary_by_spine.tsv
  basal_key_pivot.tsv
  edge_hierarchy_by_spine.tsv
  hub_rank_correlation.tsv
  hub_arm_attribution/          # per hub × spine top arms
  edge_set_jaccard.tsv          # S0 vs S2
  disagreement/IRF1_edge_diff.tsv
  disagreement/PTEN_edge_diff.tsv
```

---

## 6. Decision rules (after grid)

Promote a spine variant only if **all** hold:

1. **Basal key 8/8** neg-sig on key hallmarks (same grid as S0)
2. **Basal breadth:** `n_neg_sig ≥ 40/50` (allow small drop from 42 if attribution gain is large)
3. **Hub attribution:** median hub rank-ρ **≥ S0 + 0.05** OR IRF1 rank-ρ **≥ 0.25**
4. **Biology sanity:** top arms on IRF1 / CDKN1A not dominated by binding-only miRTar tiers

**If S1 wins:** update attribution lane + optional `α=0.5` in sensitivity docs; spine coupling stays M0.  
**If S2 wins:** headline “ENCORI-first” spine with documented M0′ gates; miRTar as functional boost.  
**If neither wins:** keep S0; ENCORI only in hub review tables (Lane 3).

---

## 7. Implementation order

| Step | Task | Est. |
|------|------|------|
| 1 | `encori_edges.py` — collapse parquet, `enc_depth`, arm resolve | 1–2 h |
| 2 | Extend `evidence_scoring.py` — `apply_encori_boost(edges, α)`, `apply_mirtar_boost(enc_edges, β)` | 1 h |
| 3 | `build_encori_m0_edges(spec)` — gates, cap, hub K | 1–2 h |
| 4 | `dual_spine_comparison.py` — orchestrate S0/S1/S2 grids | 2–3 h |
| 5 | Fix `hub_gene_comparison.tsv` — `overlap_m0` vs `overlap_mirtar_all` columns | 15 min |
| 6 | Run full grid (~30–60 min wall) | — |
| 7 | Short results memo in `output/dual_spine_comparison/README.md` | 30 min |

**Phase 2 (only if S1 or S2 wins):** Lane 2 — `enc_depth` in softmax logits (`β_share` grid) on winning edge set.

---

## 8. Risks and mitigations

| Risk | Mitigation |
|------|------------|
| S2 edge flood | strict `clipExpNum ≥ 2`, per-gene cap K, abundance floor |
| Arm name mismatch ENCORI ↔ GDC | `resolve_edges_mirna` + audit table |
| S1 α too high → coupling drop | grid α; never promote without key 8/8 |
| IRF1 wins S2 but global coupling fails | report S2 as **gene-specific attribution lane**, not cohort spine |
| Double-counting in future Lane 2 | edge boost OR share boost, not both at full strength |

---

## 9. Commands (after implementation)

```bash
# Full comparison
.venv/bin/python3 -m mirna_hallmark.dual_spine_comparison

# S1 only (quick)
.venv/bin/python3 -m mirna_hallmark.dual_spine_comparison --spines S0 S1 --alpha 0 0.5

# S2 only
.venv/bin/python3 -m mirna_hallmark.dual_spine_comparison --spines S2 --clip-min 2 --beta 0.5
```

---

## 10. Anchors from ENCORI ↔ miRTar comparison

| Metric | Value |
|--------|------:|
| M0 pairs | 5,315 |
| ENCORI collapsed pairs | 19,701 |
| ENCORI ∩ M0 | 1,352 (25.4% of M0) |
| IRF1 M0 / ENCORI∩M0 | 4 / 2 |
| IRF1 ENCORI∩miRTar (not M0) | 6 extra arms |

These numbers set expectations for S1 boost coverage and S2 IRF1 expansion.
