# Why we remove "silent" miRNA arms (and keep the 111)

Canonical: `mirna_hallmark/arm_expression.py`. Full reference: `ARM_EXPRESSION_FLOOR.md`. Concise rationale below.

## The cut (detectability, max-based)
- **REMOVE = never reaches 10 RPM in any tumor** (truly silent): **307 of 789** HE-edge arms
  (234 never-detected in the TCGA matrix + 73 absent from TCGA miRNA-seq entirely).
- **KEEP = detected at ≥10 RPM in ≥1 tumor**: **482 arms** = 220 robust (median≥10) + 262 conditional.

## Why remove the silent ones — from the data
1. **Can't act.** A removed arm never reaches 10 RPM in any of ~1,000 tumors → never loaded into RISC in functional amounts.
2. **Can't co-occupy RISC even as a group.** All removed arms **combined are ≤2% of the total miRNA pool even in the single most extreme tumor** → cannot compete with the high-abundance hubs for AGO.
3. **Can't show coupling.** Near-flat across tumors → ≈0 cross-sample variance → invisible to anti-correlation (the only cohort-level impact test we have).
4. **Only add noise** to regulator counts (`n_eff`) and any construction that isn't abundance-weighted.

## Why we KEEP the 111
111 arms are median-silent but **spike >10 RPM in a few tumors** (<1% of the cohort). We keep them ("conditional"):
- Induction in a tumor *subset* can be **real and functional there** — the cohort median hides it (this is the whole reason the cut is max-based, not median-based).
- They **self-attenuate** in the abundance-weighted aggregate wherever they're low, so keeping them is low-risk.
- We remove only the *never-detected*; `frac_expressed` flags the weakest conditional arms (seen in very few tumors).

## Edge-count impact
| set | edges | arms | genes |
|-----|-------|------|-------|
| HE all | 5,219 | 789 | 1,424 |
| **kept (expressed)** | **4,628** | 482 | 1,319 |
| removed (silent) | 591 | 307 | — |

Removing ~**half the arms drops only 11% of edges** (silent arms are low-degree); keeping the 111 spikers adds back **~443 edges** vs a strict cohort-fraction cut (1% → 4,185).

## It changes no conclusion
Acquired-coupling headline holds on the kept universe (acq wins, promiscuity harmful, weighting inert; base tightens) — the removal is **pure noise reduction**, not a functional judgement.
