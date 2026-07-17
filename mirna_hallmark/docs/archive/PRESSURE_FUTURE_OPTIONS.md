# Pressure formula — deferred options (A–E)

Companion to the evidence-scoring work in `evidence_scoring.py` / `evidence_scoring_sensitivity.py`.
These are **not** in the default spine until sensitivity + biology review say otherwise.

Context (2026-06): on M0 edges, **within-gene edge weights are often flat** (median max/min ≈ 1.8×;
57% of genes <2×), while **realized pressure ranks follow expression**, not `edge_w` (PTEN
rank-ρ ≈ 0.09). The upstream evidence table is the first lever to pull **before** changing
share/z/logRPM.

---

## D — Finer miRTarBase evidence scoring (ACTIVE TRACK)

**What:** Recompute `evidence_score` from the full interaction summary (study × experiment ×
support cross-counts), not only coarse `n_reporter_studies + n_binding_studies + …`.

**Note:** The registered spine already applies **tiered permissive** weights via
`load_mirtar_edges` (functional-MTI cross-count tiers), *not* the raw CSV `evidence_score`
from `build_edges`. Sensitivity grid includes `tiered_permissive` as the baseline row.

**Rationale:**
- Current score counts **studies per experiment class**, ignoring **Functional vs Weak MTI**
  and **reporter/protein vs binding-only** quality.
- 44% of M0 edges share `evidence_score = 8.5` → near-categorical weights.
- Cross-columns (`n_reporter__functional_mti_studies`, etc.) are already built in
  `build_edges.compute_interaction_summary_fast` but under-used.

**Expectations:**
- **Best case:** wider `edge_w` spread within genes → structural rank aligns with realized
  contribution; Basal coupling **≥ 42/50** with key 8/8; hub routes (PTEN, CDKN1A) rank miR-17
  cluster by evidence not only abundance.
- **Neutral:** coupling unchanged but **attribution tables** become trustworthy.
- **Failure mode:** over-weighting reporter/protein collapses edge count → weaker Hallmark coverage.

**Recommendation:** Run `evidence_scoring_sensitivity` grid; promote a scorer only if Basal key 8/8
holds **and** PTEN/IRF1 rank-ρ(edge_w, |contrib|) improves materially (target > 0.25).

**ENCORI extension (Phase D2):** APM’s ENCORI client (`pipeline/lncRNA_interactions/encori.py`) is
**lncRNA-first** today. For protein-coding Hallmark genes, query `miRNATarget` with
`geneType=mRNA` (or `protein_coding` per ENCORI API) per gene symbol, cache under
`data/external_cache/encori/miRNATarget_mRNA/`. Use ENCORI as:
- **boost** on existing miRTar rows (CLIP + target prediction support), not orphan flood;
- **tie-breaker** when miRTar `evidence_score` ties at 8.5;
- optional **minimum CLIP exp** filter (`clipExpNum ≥ 2`).

Do **not** merge ENCORI-only pairs into M0 spine without abundance floor + per-gene cap (M7-style).

---

## A — Gene-specific specificity in `edge_w`

**Formula:**
```
spec(m,g) = log1p(ev_m,g) / Σ_{g' ∈ targets(m)} log1p(ev_m,g')
edge_w(m,g) = spec(m,g) / log1p(degree(m))   [optional mild hub cap]
```

**Rationale:** Global `evidence_mass` penalizes hub miRNAs across the universe but leaves
**within-gene** weights flat. `spec` is “what fraction of miRNA m’s curated budget lands on g?”

**Expectations:**
- Stronger hierarchy on crowded genes (PTEN 81 arms).
- Basal coupling: tested as `inner_sum_log` → **37/50** but **key 8/8** (slightly weaker breadth).
- Realized rank may still follow abundance unless combined with sharper upstream scoring (D).

**Recommendation:** Sensitivity-only (`target_norm=evidence_sum_log`) until D improves raw scores;
consider **spec + D2** combined before spine change.

---

## B — Specificity in softmax logits (evidence-in-share)

**Formula:**
```
logit(m,g,s) = (x_m,s − median_m) + log1p(ev_m,g)
share = softmax_m∈R(g)(logit)
edge_w = 1   (or degree_only)
c = share × z × logRPM × edge_w
```

**Rationale:** Fuses “who is abundant” and “who has evidence on g” in one competition step.

**Expectations (tested):** Basal **41/50**, key **8/8** — safe but **not better** than spine.
Harder to interpret; does not fix E2F co-regulation.

**Recommendation:** Footnote / sensitivity row only. Keep evidence in `edge_w` for default reporting.

**Lane 2 run (2026-06, `encori_share_sensitivity`):** ENCORI `enc_depth` in softmax logits on
fixed M0 (`β_share` grid), `edge_w` unchanged (`evidence_mass`). β=0 reproduces S0 coupling
(42/50, key 8/8). β≥0.25 **degrades** hub median rank-ρ; IRF1 ρ → −0.40 at β≥0.5.
**S1 edge-weight boost (α·enc_depth) remains the better ENCORI lane** for attribution.

---

## C — Decouple structural vs dynamic pressure  ✅ IMPLEMENTED (resolution axis, 2026-06-27)

**Status:** the structural/abundance split is **live for the within-gene edge-resolution &
decoupling axis**. `gene_struct_share_*` (mode `softmax`, `c_struct = edge_w·sm`, no `logRPM`
double-count) ships alongside the abundance-driven `gene_share_*`; the §9 quadrant and §14
competitor logic key off the structural rank. Coupling spine and cross-state mass untouched. See
FORMULAS §5a/§9/§11c and the ledger entry (2026-06-27). Below is the original design note.



**Definitions:**
```
P_struct(g)   = Σ_m spec(m,g)                         static structural load
P_dyn(g,s)    = Σ_m spec(m,g) × share(m,s) × z(m,s)   sample dynamics (no logRPM)
P_dose(g,s)   = P_dyn × anchored abundance            optional logRPM or absratio
```

**Rationale:** Single product `c = expr × edge_w` conflates “who should regulate g” with “who is
up today.” Coupling should use **P_dyn**; attribution papers should cite **P_struct** or
`softmax_logrpm` shares.

**Expectations:**
- Clearer prose in methods; may allow dropping logRPM from coupling without losing key 8/8.
- Does not require new edge data — reorganizes existing terms.

**Recommendation:** Report **dual tracks** in hub-gene tables before changing spine. Coupling spine
candidate: `softmax_z + spec edge_w` (no logRPM).

---

## E — TargetScan as tie-breaker only (not M7 fusion)

**Formula (miRTar rows only):**
```
edge_w = spec_miRTar(m,g) × (1 + α · ts_spec(m,g))     α ≈ 0.25–0.5
ts_spec(m,g) = ts_weight(m,g) / Σ_g' ts_weight(m,g')
```

**Rationale:** TS adds continuous specificity where miRTar ties at 8.5; avoids M11/M7 orphan
dilution (Basal 33–41/50 in hybrid grid).

**Expectations:**
- Modest rank reordering on tied edges; IRF1 alt arms (miR-130b) may rise on sequence support.
- Risk: sequence prior disagrees with functional validation (PTEN → miR-19).

**Recommendation:** Only after D (+ optional ENCORI boost). Compare α on hub panel; never replace
miRTar in M0 without coupling grid sign-off.

---

## Decision order (agreed 2026-06)

1. **D** upstream miRTar (+ ENCORI mRNA fetch) — current priority  
2. **C** dual-track reporting (low code cost, clarity)  
3. **A** spec normalization if D alone insufficient  
4. **E** TS tie-break on miRTar rows  
5. **B** evidence-in-softmax — last resort (interpretability cost)

---

## Related outputs

| Path | Content |
|------|---------|
| `output/pressure_evidence_sensitivity/` | edge mass log variants (inner/outer log) |
| `output/pressure_layer_comparison/` | L1–L4 layers, hub cancellation |
| `output/evidence_scoring_sensitivity/` | upstream scorer grid (this track) |
