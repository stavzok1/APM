# S1 + share validation memo

Generated after independent recomputation vs dual-spine baselines.

## 1. Coupling — PASS

| Spec | expr_mode | Basal neg-sig | Key 8/8 | dual-spine match? |
|------|-----------|---------------|---------|-------------------|
| S0 | softmax_z_logrpm | **42/50** | 8/8 | Yes (S0) |
| S1 α=0.5 | softmax_z_logrpm | **41/50** | 8/8 | Yes (S1 α=0.5) |
| Lane2 β=0 | softmax_z_logrpm_enc | **42/50** | 8/8 | Yes (M0 ref) |
| S1 α=0.5, β=0.25 | softmax_z_logrpm_enc | **41/50** | 8/8 | Stable |
| S1 α=0.5, β=0.5 | softmax_z_logrpm_enc | **41/50** | 8/8 | Stable |

Coupling is **stable** under S1+share; no gain vs S1 alone; still −1 vs S0 (HEME_METABOLISM).

## 2. Hub rank-ρ — mixed (do not promote combined lane yet)

Fair recompute (same code path):

| Spec | Hub median ρ | PTEN (n=81) | CDKN1A | IRF1 (n=4) | BCL2L11 |
|------|--------------|-------------|--------|------------|---------|
| S0 | +0.08 | +0.09 | **−0.15** | +0.20 | +0.08 |
| S1 α=0.5 (dual expr) | +0.09 | +0.09 | **−0.15** | +0.20 | +0.31 |
| S1 α=0.5, β=0 | +0.09 | +0.09 | **−0.15** | +0.20 | +0.31 |
| **S1 α=0.5, β=0.25** | **+0.27** | **+0.06** | **−0.20** | **+0.80** | +0.68 |
| S1 α=0.25, β=0.25 | +0.23 | +0.01 | **−0.24** | +0.80 | +0.59 |
| Lane2 β=0.5 (α=0) | **−0.21** | **−0.21** | **−0.31** | **−0.40** | +0.17 |

### PTEN — FAIL for attribution claim

PTEN ρ stays **~0.05–0.09** with S1+share (81 arms). No material gain vs S0 (+0.09) or dual-spine S1 (+0.15 in original grid). **Combined lane does not fix the main crowded-hub problem.**

### CDKN1A — FAIL

Negative rank-ρ across **all** specs including S1+share.

### IRF1 — FRAGILE (n=4)

At S1 α=0.5, β=0.25:

| miRNA | edge_w | share bonus | mean |contrib| |
|-------|--------|-------------|----------------|
| miR-106b-5p | 0.51 | 0.96 | 0.92 |
| miR-130b-3p | 0.55 | 1.50 | 0.87 |
| miR-23a-3p | 0.44 | 0.00 | 0.49 |
| miR-383-5p.1 | ~0 | 0.00 | ~0 |

ρ=0.80 is driven by **2 ENCORI-supported arms** out of 4; not robust for cohort-level claims.

**Permutation test** (shuffle `share_logit_bonus` across IRF1 arms, n=500):  
observed ρ=0.80, perm mean=0.47, **p=0.31** — not significant at n=4.

### BCL2L11 — only consistent winner

ρ rises from ~0.08 (S0) to ~0.59–0.72 under S1+share.

## 3. Synergy interpretation (revised)

- **Share without S1 edge boost (α=0)** hurts attribution — confirmed.
- **Share with S1 edge boost** raises hub **median** ρ mainly via BCL2L11/TGFBR2/VIM; **not PTEN/CDKN1A**.
- IRF1 improvement is **real but small-n** and should not drive spine promotion.

## 4. Verdict

| Criterion | S1+share (α=0.5, β=0.25) |
|-----------|---------------------------|
| Coupling ≥ S0 breadth | No (41 vs 42) |
| Key 8/8 | Yes |
| PTEN rank-ρ gain | **No** |
| CDKN1A rank-ρ | **No** |
| IRF1 rank-ρ | Fragile (n=4) |
| **Promote to registered spine?** | **No** |

**Keep S0 for coupling registry.** Use **S1 α-only** (dual-spine, softmax_z_logrpm) for attribution sensitivity. S1+share may be cited for **gene-specific** panels (BCL2L11, IRF1 exploratory) only.
