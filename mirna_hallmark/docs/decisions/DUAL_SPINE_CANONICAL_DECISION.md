# Canonical spine decision — coupling vs attribution (2026-06-22)

This memo locks the **registered** miRNA→Hallmark pressure spine for **coupling claims**
and documents **sensitivity lanes** for attribution. It supersedes the automated
`README.md` verdict (`S2_promoted`), which applied hub-rank rules without edge-set
comparability checks.

**Representative runs compared below**

| Label | run_id | Role |
|-------|--------|------|
| **S0** | `S0` | miRTar M0, `tiered_permissive`, 5,315 pairs |
| **S1** | `S1__alpha=0.5__base=confidence_logclass` | Same 5,315 pairs; `evidence_score += 0.5·enc_depth` |
| **S2** | `S2__beta=0.0__clip_min=2__top_k=15` | ENCORI-built graph, 10,652 pairs (β=0 = no miRTar boost) |

Shared pressure layer (all three): `expr_mode=softmax_z_logrpm`, `target_norm=evidence_mass`,
`aggregate=sum`, AGO gate on.

---

## 1. Primary coupling benchmark

**Metric:** partial Spearman ρ between Hallmark pressure and Hallmark target expression,
| CPE + HRD + **e2f_g2m** + mean member CN + within-subtype z-scoring (`rho_prolif_cn_wsd_adj`).
**Significance:** ρ < 0 and BH-FDR q < 0.10 within PAM50 subtype.

### 1.1 Breadth by subtype (neg-sig hallmarks / 50)

| Subtype | S0 | S1 α=0.5 | S2 enc | S2/S0 ratio |
|---------|----|----------|--------|-------------|
| **Basal** | **42** | **41** | 46 | 1.10× |
| LumA | 17 | 16 | 31 | **1.82×** |
| LumB | 19 | 15 | 36 | **1.89×** |
| Her2 | 0 | 0 | 0 | — |

**Key 8 hallmarks (Basal):** S0 **8/8**, S1 **8/8**, S2 **8/8** — no sign flips on any key program.

### 1.2 Effect size (Basal median ρ, primary metric)

| Spine | Median ρ | Mean ρ among neg-sig |
|-------|----------|----------------------|
| S0 | −0.288 | — |
| S1 α=0.5 | −0.278 | — |
| S2 enc | −0.369 | — |

S2 shows **stronger** (more negative) median ρ, but on a **2× denser graph** — see §4.

### 1.3 Basal hallmark-set stability (Jaccard of neg-sig sets)

| Pair | Shared | Jaccard | S0-only lost | S2-only gained |
|------|--------|---------|--------------|----------------|
| S0 vs S1 | 41/42 | **0.976** | HEME_METABOLISM | — |
| S0 vs S2 | 42/42 | 0.913 | — | 4 programs* |
| S1 vs S2 | 41/41 | 0.891 | — | 5 programs |

\*S2-only Basal neg-sig gains (not in S0): `FATTY_ACID_METABOLISM`, `KRAS_SIGNALING_DN`,
`PROTEIN_SECRETION`, `SPERMATOGENESIS`.

**Hallmark-level ρ correlation** (41 hallmarks neg-sig in both S0 and S1): Spearman **0.965**
(p ≈ 2.6×10⁻²⁴) — near-perfect rank stability on the fair head-to-head lane.

### 1.4 Key 8 ρ deltas (Basal, primary metric)

| Program | S0 ρ | Δ(S1−S0) | Δ(S2−S0) |
|---------|------|----------|----------|
| APOPTOSIS | −0.541 | +0.000 | +0.004 |
| P53_PATHWAY | −0.477 | +0.020 | −0.080 |
| IFN-γ | −0.332 | +0.006 | −0.051 |
| EMT | −0.446 | +0.002 | +0.007 |
| PI3K/AKT/mTOR | −0.305 | −0.025 | −0.029 |
| NOTCH | −0.176 | +0.022 | −0.144 |
| TGF-β | −0.245 | +0.003 | −0.159 |
| TNFα/NFκB | −0.496 | +0.015 | −0.011 |

No key hallmark **sign flips**. S1 perturbations are ≤0.025 |ρ|; S2 shifts some programs
(NOTCH, TGF-β) toward stronger negative coupling on the ENCORI graph.

---

## 2. Proliferation robustness ladder (Basal)

Partial ρ | CPE + HRD + proxy (no CN/WSD in this table — cohort/Basal scope from
`hallmark_coupling_prolif`). **Survives** = ρ < 0 and p < 0.05.

| Proxy | S0 basal | S0 key 8 | S1 α=0.5 basal | S1 key 8 | S2 basal | S2 key 8 |
|-------|----------|----------|----------------|----------|----------|----------|
| e2f_g2m (primary) | 40/50 | 8/8 | 39/50 | 8/8 | 46/50 | 8/8 |
| mki67 | 30/50 | 6/8 | 30/50 | **7/8** | 39/50 | 8/8 |
| ortho_noE2F_MYC | 36/50 | 8/8 | 36/50 | 8/8 | 40/50 | 8/8 |
| pc1 | 40/50 | 8/8 | 39/50 | 8/8 | 46/50 | 8/8 |
| **ALL four** | **28/50** | **6/8** | **29/50** | **7/8** | 38/50 | 8/8 |

**Interpretation**

- **S1 vs S0:** −1 hallmark on e2f_g2m/pc1 (HEME loss); **+1 key** on mki67 (PI3K/AKT/mTOR
  borderline under MKI67 alone on S0, survives on S1). ALL4 ladder: 28→29 basal, 6→7 key.
- **S2 vs S0:** Higher survive counts everywhere, but **non-Basal inflation is ~1.8×** (§1.1)
  and the graph adds ~5k ENCORI-only edges — breadth gains are **not** a fair upgrade test.

Cohort scope (e2f_g2m): S0 **41/50** (key 7/8), S1 **41/50** (key 7/8), S2 **44/50** (key 8/8).

---

## 3. Edge-set comparability (why S2 ≠ fair coupling benchmark)

| | S0 / S1 | S2 enc |
|--|---------|--------|
| Pairs | 5,315 | 10,652 |
| Jaccard vs S0 | 1.0 / 1.0 | **0.060** |
| Shared pairs | — | 898 / 5,315 (16.9%) |
| ENCORI∩M0 (hub table) | — | 1,352 / 5,315 (25.4%) |
| Median within-gene max/min edge_w ratio | ~1.86–1.97 | **1.34** (flatter weights) |
| frac genes with <2× spread | ~51% | **~86%** |

S2 **compresses** within-gene edge weight spread while **doubling** pair count → more
Hallmarks receive pressure mass → **artifactually** more neg-sig calls outside Basal.

---

## 4. Attribution QA (separate from coupling registration)

Hub **rank-ρ**: Spearman(edge_w, |mean arm contribution|) per hub gene — tests whether
high-evidence arms drive pressure, not whether coupling holds.

| | S0 | S1 α=0.5 | S2 enc |
|--|----|----------|--------|
| Median hub rank-ρ | 0.113 | **0.321** (+0.21) | 0.040 |
| PTEN | 0.076 | **0.153** | 0.275 |
| IRF1 (n=4 arms) | **0.775** | 0.400 | **−0.402** |

S1 improves attribution coherence on fixed M0 without breaking coupling.
S2 **inverts IRF1** rank-ρ despite higher PTEN — cliff-gene story degrades on ENCORI-first graph.

---

## 5. Decision rules (applied)

| Criterion | S0 | S1 α=0.5 | S2 enc |
|-----------|----|----------|--------|
| Fair edge set (M0 5,315) | ✓ | ✓ | ✗ |
| Basal key 8/8 (primary) | ✓ | ✓ | ✓* |
| Basal breadth ≥40/50 | ✓ (42) | ✓ (41) | ✓ (46)* |
| Basal ρ rank stability vs S0 | — | ✓ (ρ≈0.97) | partial |
| Non-Basal inflation ≤1.2× | ✓ | ✓ | ✗ (1.8×) |
| IRF1 attribution sane | ✓ | ✓ | ✗ |
| ALL4 prolif ladder (key 8) | 6/8 | **7/8** | 8/8* |

\*Passes on raw counts but **fails comparability** — different graph, inflated Lum strata.

---

## 6. Canonical assignment (effective immediately for new work)

| Layer | Canonical choice | Notes |
|-------|------------------|-------|
| **Registered production spine** | **S1, α=0.5** | M0 miRTar pairs; `confidence_logclass` + ENCORI depth in `evidence_score`; wired in `config.PRESSURE_*` + `pressure_build` (2026-06-23) |
| **Legacy S0 reference** | S0 | tiered_permissive only; kept for sensitivity / dual_spine grid |
| **ENCORI review / disagreement** | `encori/mirtar_comparison/`, Lane-2 share grid | Not registered coupling |
| **S2 ENCORI-first graph** | Exploratory only | Do not cite Basal 46/50 vs registered 41–42/50 |
| **S1+share (softmax bonus)** | Rejected | Coupling 41/50; PTEN flat; see `encori_share_sensitivity/VALIDATION.md` |

### What to cite in papers / registry

- **Coupling headline:** S1 Basal **41/50** neg-sig (prolif+CN+WSD), key **8/8**; hybrid M0 Basal **42/50**; cohort gated **26/50**, partial CPE+HRD **21/50**.
- **Robust coupling:** require ALL4 ladder or minimum ortho_noE2F_MYC + pc1 survival for strong claims.
- **Attribution:** S1 α=0.5 is now the **registered production spine** (`config.PRESSURE_EVIDENCE_SCORER` + `PRESSURE_ENCORI_ALPHA`).

### Outputs referenced

| File | Content |
|------|---------|
| `coupling_summary_by_spine.tsv` | Subtype breadth |
| `coupling_detail_by_spine.tsv` | Per-hallmark ρ/q |
| `basal_key_pivot.tsv` | Key 8 wide |
| `canonical_decision_coupling.tsv` | Side-by-side subtype stats |
| `canonical_basal_hallmark_by_spine.tsv` | Per-hallmark Basal flags |
| `prolif_robustness_summary.tsv` | Four-proxy ladder |
| `edge_set_jaccard.tsv` | S2 vs S0 overlap |
| `hub_rank_correlation.tsv` | Attribution appendix |

---

## 7. Not yet wired (await explicit sign-off)

- `mirna_hallmark/docs/DISCOVERY_REGISTRY.md` — spine paragraph
- `pressure_build.py` / `config.PRESSURE_*` — still S0 defaults (correct for coupling)
- Optional: `--attribution-spine S1 --alpha 0.5` flag in downstream runners

---

## 8. Claim audit (S0 vs S1 α=0.5, 2026-06-22)

Module: `mirna_hallmark.spine_claim_audit` → `claim_audit/`

| Check | Result |
|-------|--------|
| Hallmark pressure Spearman (median across programs) | **0.994** (key-8 min **0.983**) |
| Registry claims unchanged | **8/9** (only MH-30 headline 42→41) |
| Gene pressure Basal survival flips (CDKN1A/PTEN/TGFBR2/VIM/IRF1) | **0** |
| Hub route Basal survival flips (e.g. PTEN–miR-106b, CDKN1A–miR-17) | **0** |
| Bootstrap Basal neg-sig (n=500, primary metric) | S0 mean **39.7**, S1 **39.3**; median **40** both |
| P(S0 > S1) neg-sig count | **0.42** (not significant — 40% ties) |
| HEME neg-sig flip rate (S0+ but S1−) | **12%** of bootstraps |
| Key 8 identical S0 vs S1 | **422/500** bootstraps |

**Conclusion:** No substantive registry biology changes under S1. The 42 vs 41 split is **bootstrap-unstable** for both spines (median 40/50); HEME is a borderline call. Unified S1 is **safe for gene/hub claims**; only the breadth headline moves.

---

*Generated from grid under `mirna_hallmark/output/dual_spine_comparison/`; prolif ladder
re-run 2026-06-22; claim audit `spine_claim_audit` 2026-06-22.*
