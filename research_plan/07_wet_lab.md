# 07 — Wet-lab validation plan

Wet-lab work is optional for thesis defensibility (the computational story stands on its own) but **high-leverage** because the platform produces per-element, per-variant, per-class predictions that most computational frameworks cannot. The experiments here are chosen to defend the *distinctive* claims of the platform, not to reconfirm textbook APM biology.

## Selection principles

1. **Every experiment must validate a specific numbered hypothesis in `02_hypothesis_catalog.md`**, not a general biological fact.
2. **Prefer experiments that only the platform can even design** — i.e., that require per-element or per-variant resolution.
3. **Two primary lines of defense**: methodological (is the attribution causal at element/variant resolution?) and translational (are the mechanism-attributed classes therapeutically actionable?).
4. Cell lines first; organoids / in vivo only if an experiment produces a headline-grade result that needs native-context confirmation.

---

## Experiment roster

| ID | Title | Defends | Tier | Time | Cost |
|---|---|---|---|---|---|
| **W1** | CRISPRi of a CRML-predicted enhancer | C2 element resolution; H1, H7 | Methodological core | 4–6 mo | mid |
| **W2** | Allele-specific reporter for high-ΔFIMO CRML variants | C2 variant resolution; H1 | Methodological core | 2–3 mo | low |
| **W3** | DNMT / EZH2 inhibitor reversal in cold-intrinsic class | C3, C5 therapeutic stratification; H19, H21, H22, H59 | Translational | 4–6 mo | mid |
| **W4** | STING agonist response in cold-primed class | H48, H56, H59 (primed-not-triggered) | Translational | 4–6 mo | mid |
| **W5** | NK killing assay across MHCI ± / NKL ± quadrants | H14 dual-axis escape | Hypothesis-specific | 3 mo | mid |
| **W6** | Allele-specific reporter for regulatory ASE loss | H46, H63 per-allele attribution | Methodological | 2–3 mo | low |

## Recommended minimum package

**W1 + W3**. One methodological, one translational. Together these answer the two most common reviewer objections in one PhD-scale wet-lab commitment.

If budget and time allow:
- Add **W2** — it is cheap (a weekend of cloning + reporter reads) and directly defends the variant-Δ claim that underpins CRML.
- Add **W6** — trivially extends W2 to per-allele contexts; nearly free once W2 is in place.

Skip **W4** and **W5** unless those hypotheses become the paper's centerpiece.

---

## W1 — CRISPRi of a CRML-predicted enhancer

### Design

1. Rank all cCREs in the Element Focus table by `CRML_enhancer_score` × `ABC_calibrated × chip_support`. Top 10 across the 66-gene panel.
2. Pick 2 to test, constrained by:
   - cCRE is distal (≥ 10 kb from TSS) to isolate enhancer contribution from promoter
   - High ABC edge to a single panel gene (clean target)
   - Target gene is well-expressed in the chosen cell line (MDA-MB-231 basal or MCF7 luminal)
   - ChIP support for IRF1 or STAT1 in a matched-subtype biosample
3. CRISPRi: dCas9-KRAB stable line; 3 gRNAs per element; non-targeting control + a safe-harbor control.
4. Readouts:
   - Target APM mRNA (qPCR) at baseline and under IFN-γ stimulation (24 h, 100 ng/mL)
   - MHC-I surface by flow (anti-HLA-A/B/C W6/32)
   - Optional: bulk RNA-seq to catch trans effects.

### Success criteria

- At least one of the two elements shows ≥ 30 % reduction in target expression at baseline or under IFN-γ vs. non-targeting control.
- Effect direction matches CRML's sign prediction.
- Bulk RNA-seq shows effect localised to cis-target, not pathway-wide.

### Failure mode and interpretation

- If neither element shows an effect → CRML is overcalling. Response: examine whether the hit element has ABC-edge evidence in the chosen cell line's subtype specifically (not generic). Tighten the scoring.
- If effect exists but is smaller than CRML predicts → quantitative calibration issue; useful data anyway.

---

## W2 — Allele-specific reporter for high-ΔFIMO CRML variants

### Design

1. Take top 5 SNVs by `|ΔFIMO score| × ChIP-support × ABC-edge-to-panel-gene` observed in the cohort.
2. Clone ±200 bp around each SNV (WT and variant haplotypes) into a minimal-promoter luciferase reporter.
3. Transfect into an appropriate cell line; treat ± IFN-γ.
4. Measure WT vs variant luciferase output.

### Success criteria

- At least 3 of 5 variants show WT > variant (or vice versa matching ΔFIMO sign) at p < 0.05 in biological triplicate.
- Direction of effect matches ΔFIMO sign.

### Why this is powerful

Luciferase allele comparisons are the gold-standard falsification for single-base enhancer claims. A clean 3/5 hit rate defends CRML as a real variant-resolved quantity in one figure.

---

## W3 — DNMT / EZH2 inhibitor reversal in cold-intrinsic class

### Design

1. Score every BRCA cell line in CCLE against the four-class hot/cold taxonomy (H59) using the same L2 composite scores used on TCGA tumors.
2. Pick 3 predicted **cold-intrinsic** TNBC lines (high `COMB_HLA`, low `CRML`, low CNV damage on APM, high `EZH2/SUZ12/EED` expression) and 2 predicted **cold-primed** or **hot-competent** controls.
3. Treat with 5-azacytidine (1–10 μM, 72 h) or tazemetostat (EZH2i, 5 μM, 5–7 d) vs vehicle.
4. Readouts:
   - MHC-I surface by flow
   - APM panel expression (Nanostring or targeted qPCR of 20 APM genes)
   - Optional: co-culture with NK92 or primary NK for killing assay.

### Success criteria

- Cold-intrinsic lines show ≥ 2× MHC-I surface induction under 5-aza / EZH2i; control classes show < 1.3×.
- APM panel re-activation correlates with predicted reversible loci (those with high β, low CRML).

### Why it earns its keep

This turns the four-class hot/cold decomposition from descriptive into a **biomarker-for-therapy** story. If cold-intrinsic class predicts response to demethylating agents and cold-primed does not, that is a direct translational implication publishable on its own.

---

## W4 — STING agonist response in cold-primed class

### Design

1. Pick 2 predicted **cold-primed** TNBC lines (high `JSI` + high `COHER` + low `CRML`) and 2 cold-intrinsic controls.
2. Treat with ADU-S100 or diABZI (STING agonist) at standard doses, 24–48 h.
3. Measure IFN-β secretion (ELISA), STAT1 phosphorylation (Western), APM panel induction (qPCR), MHC-I surface (flow).

### Success criteria

- Cold-primed lines show robust STING-induced APM induction (≥ 2× MHC-I, p-STAT1 up); cold-intrinsic lines show blunted response despite similar STING pathway integrity.

### Coverage

Defends H48 (cGAS–STING parallel) and H56 (dual-axis intactness) therapeutically.

---

## W5 — NK killing across quadrants

### Design

1. Pick 4 BRCA lines, one per MHCI ± / NKL ± quadrant from H14.
2. Co-culture with NK92 (or primary NK from healthy donor) at 5:1, 10:1 E:T.
3. Measure killing (flow, CellTox, or impedance-based).

### Success criteria

- MHCI− / NKL+ quadrant shows highest killing; MHCI+ / NKL− lowest; intermediate cases between.

### Coverage

Defends H14 as causal rather than correlational.

---

## W6 — Per-allele reporter for regulatory ASE loss

### Design

1. Identify 3 cohort samples with **regulatory ASE** (allele-1 expression drop without copy loss) at HLA-A, B, or C based on H46/H63 classification.
2. Identify the haplotype-defining SNVs in the promoter/enhancer region and clone allele-1 vs allele-2 haplotypes into reporter.
3. Compare.

### Success criteria

- Observed reporter ratio matches the per-allele RNA-seq ratio in the sample ± 30 %.

### Coverage

Defends H46 / H63 as a measurement of regulatory, not genomic, allele loss.

---

## Practical sequencing (6-month wet-lab window)

| Month | Activity |
|---|---|
| 0–1 | Cell-line scoring against cohort taxonomy; CCLE matching; pick lines; order reagents |
| 1–2 | W2 / W6 cloning and pilot luciferase |
| 2–4 | W1 CRISPRi stable line generation; W3 inhibitor treatment pilot |
| 4–5 | W1 full readout; W3 full readout with NK killing if available |
| 5–6 | Analysis, figure drafting; decision on W4 / W5 based on how W1 + W3 land |

Runs in parallel with Aim-2 computational work. Does not block the timeline.

---

## Mapping back to the catalog

| Hypothesis | Wet-lab backing |
|---|---|
| H1 (CRML beyond CNV/meth) | W1, W2 |
| H7 (R_g × burden) | W1 (secondary) |
| H14 (dual-axis escape) | W5 |
| H19/H21/H22 (methylation axis) | W3 |
| H46/H63 (per-allele ASE) | W6 |
| H48/H56 (cGAS–STING) | W4 |
| H59 (four-class hot/cold decomposition) | W3, W4 |
