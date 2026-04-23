# 04 — Conservation-as-prior methodology

Every small-sample module (HiChIP, cell-line ATAC, ChIP biosamples) becomes usable for per-tumor inference if you treat subtype-conservation as a prior rather than ground truth. The 13 TCGA-BRCA HiChIP tumor samples calibrate that prior.

---

## The core idea

Let \( E_{s,g} \) be a latent regulatory feature (e.g., enhancer-gene edge strength, ATAC accessibility, boundary strength) for sample `s` at gene/element `g`. For most BRCA tumors \( E_{s,g} \) is unobserved.

Define the subtype-consensus prior:

\[
\mu_{\text{sub}(s), g} = \text{aggregate}\big( E_{b, g} : b \in \text{biosamples with subtype} = \text{sub}(s) \big)
\]

computed over reference biosamples (cell lines, primary ENCODE/BioGPS, etc.).

We then use \( \mu_{\text{sub}(s), g} \) as a prior for \( E_{s, g} \), and build features on **(i) the prior itself** and **(ii) the residual** once per-sample evidence becomes available.

---

## The 13 HiChIP tumor samples — calibration anchor

You have 13 BRCA tumor HiChIP samples (`tcga_brca_hichip.csv`). These are invaluable not because they cover the cohort (they do not) but because they let you **calibrate the subtype-consensus prior to real tumor ground truth** before extrapolating.

### Calibration workflow

1. For each loop `e → g` in the 13 tumor HiChIP, compute a conserved prior from same-subtype reference biosamples: `prior_loop[e→g]`.
2. Fit a calibration function on the 13 tumors:

   \[
   P(\text{loop observed} \mid \text{prior}) = \sigma\big(\alpha + \beta \cdot \text{prior}\big)
   \]

   plus a continuous version for loop strength.
3. Report calibration quality (Brier, calibration curve, Spearman).
4. For the remaining ~1000+ BRCA tumors without per-tumor HiChIP, apply the calibrated prior:

   - Every enhancer–gene edge gets a probability + expected strength.
   - Feature builders (`build_enhancer_redundancy`, `build_crml_*`, `build_mirna_pressure` when it uses enhancer context) use **calibrated expected strength**, not raw prior.

This upgrades every hypothesis using `gene_links` or enhancer-gene inference in the entire catalog. It is the single most consequential methodological move and should be done first in the feature-builder phase.

### Secondary use of the 13: subtype-specific calibration split

Split the 13 by PAM50 (basal vs luminal vs HER2; numbers permitting). Fit subtype-specific calibration to check that the prior does not systematically over- or under-predict in one subtype. Fallback to a shared prior with a subtype offset if per-subtype calibration is under-powered.

### Tertiary use: held-out validation of mediation claims

Any hypothesis whose claim runs through enhancer–gene edges can be cross-validated within the 13 (using per-tumor HiChIP as ground truth) before being extrapolated cohort-wide.

---

## Where conservation-as-prior applies across the catalog

| Data layer | Subtype-consensus source | Calibrated via |
|---|---|---|
| HiChIP loops | ENCODE/GEO BRCA cell lines (MCF7, HCC1954, MDA-MB-231, etc.) | **13 tumor HiChIP** (primary anchor) |
| ATAC peaks | 74 TCGA-BRCA peak set | Per-tumor ATAC when the sample is in the ~74 cohort |
| ChIP binding | Same-subtype ENCODE ChIP biosamples | None; treat prior as best available |
| TAD boundaries | 24-biosample consensus (already in Element Focus) | 13 tumor HiChIP for boundary loops |
| ABC scores | Breast cell-line ABC | 13 HiChIP + consensus ATAC |
| SCREEN cCREs | Pan-tissue consensus | None; treat as element catalog |

---

## Conservation-normalized feature spec

Every enhancer-dependent feature gets three variants:

| Variant | Semantics |
|---|---|
| `<feat>_raw` | Observed per-sample value if present; else missing |
| `<feat>_prior` | Subtype-consensus value |
| `<feat>_resid` | `raw − prior` (or log-ratio / scaled residual) |

The hypothesis tests decide per-claim whether to use raw, prior, or residual — typically residual (H53). Using raw alone is reserved for samples inside a per-tumor module (e.g., 13 HiChIP or ~74 ATAC).

---

## Conservation levels to track (per enhancer / boundary)

For each enhancer \( e \) to APM gene \( g \):

| Level | Definition |
|---|---|
| L_core | enhancer present at ≥k reference biosamples across all subtypes, all methods agree |
| L_subtype | enhancer present only within one subtype |
| L_opportunistic | enhancer in ≤1 biosample (questionable) |

`build_enhancer_redundancy` returns a breakdown by level, so H54 ("core enhancers carry most of the effect") can be tested directly.

Similarly for TAD boundaries across 24 biosamples (H9 weights), and for ChIP binding (consensus strength).

---

## Subtype-specific vs cohort-wide effects

Conservation also lets you state the effect explicitly by subtype.

- Step 1: fit subtype-specific effects of each feature on expression (interactions with PAM50).
- Step 2: in the conservation-prior variant, also fit subtype-specific priors (e.g., luminal-specific HiChIP consensus).
- Step 3: compare the variance explained by (a) shared prior + interaction, vs (b) subtype-specific prior + main effect.

Winning (b) supports claim C4 (platform is subtype-aware). Claim C2 (element-resolved causal) uses (a) where subtypes share the mechanism.

---

## Why this is defensible

1. Ground truth HiChIP does not exist for every tumor in any large cohort. Conservation-prior is not a workaround; it is the standard way to make enhancer-gene inference scalable.
2. The 13 tumor HiChIP samples ensure the prior is not biased by cell-line artifacts, which is a documented problem in enhancer catalogs.
3. Residualization (`<feat>_resid`) preserves a clean causal interpretation: the "prior" is the part of regulatory architecture shared with the tissue class; the "residual" is the part unique to the tumor.
4. Element-resolved claims (CRML, BCD) become honest because they are weighted by calibrated confidence in the enhancer-gene edge rather than by a raw cell-line edge.

---

## Open experimental leverage

- Additional tumor HiChIP runs would directly lower calibration variance. Each new tumor HiChIP sample **compounds** by tightening calibration for the entire cohort.
- If experimental budget allows, prioritize HiChIP in subtypes currently underrepresented in the 13 (check subtype distribution of the 13 first).
