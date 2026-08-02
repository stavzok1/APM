# Germline cis-eQTL MR — restore an exogenous edge-existence handle (MH-159, 2026-07-18)

**Goal:** the program's top open item — an EXOGENOUS validation that curated miRNA→target edges are real,
since both CN instruments are dead (MH-124r/126/133) and guide/passenger is foreclosed (MH-157). Germline
genotype is genuinely exogenous (randomized at meiosis; immune to reverse causation and tumour-state
confounding), so a germline cis-eQTL for a miRNA is a Mendelian-randomization instrument for its dose.

## What was built (persists in this dir)
- **Cohort:** 1,148 TCGA-BRCA blood-normal Affy SNP6 birdseed files (`data/SNV/germline_array/raw/`,
  genotypes PRE-CALLED 0/1/2 + confidence). Overlap with tumour miRNA+mRNA cohort = **1,075** (near-complete).
- **Annotation:** GDC hg38-remapped SNP6 probeset (`snp6.na35.remap.hg38.subset.txt.gz`, MD5-verified),
  788k usable SNP probes; **497/506 miRNA loci (98%) have cis-SNPs within ±1Mb, median 519/locus.**
- **Genotype matrix:** `geno_cis.npy` (1075 × 128,263 cis-SNPs, int8, conf>0.1→missing), `geno_participants.txt`,
  `geno_cis_annot.tsv`. QC-pass (MAF≥0.05, call≥0.95): 98,027 SNPs. Scripts: `extract_geno.py`.

## First stage (single best cis-SNP; `firststage.py`)
Winner's-curse-guarded by a per-arm permutation null (max-F over ~500 cis-SNPs). **38/334 arms perm-sig at
F>10, only 2 at F>20, 1 at F>30.** ⚠ Instrument signal DOMINATED by the imprinted DLK1-DIO3 polycistron
(chr14q32: miR-370/299/485/433/134/487b/409) — where the exclusion restriction is compromised by
co-transcription (same problem as CN clusters).

## Reduced-form MR (`mr.py`) + target specificity (`specificity.py`)
Sign-aligned reduced form (SNP↑→miRNA↑→target↓ ⇒ repression-consistent = positive), C=[PCs,CPE,target_cn].
- Edge-level (128 curated edges, 14 arms): mean +0.0074, p=0.003. **But the unit is the ARM** (per memory
  `exogenous-identification-or-nothing` axiom 2): arm-level sign-test 82–86% arms repression-dir,
  p=0.006–0.033; arm-level t marginal (p=0.04–0.095). Biggest contributor miR-26a-5p is the textbook
  host-pleiotropy case (CTDSP2/CTDSPL host).
- ⛔ **DECISIVE — TARGET SPECIFICITY FAILS:** curated aligned-π **+0.0027** vs matched non-target genes
  **+0.0049** (paired Δ=−0.0022, **p=0.66**, sign p=0.89). Non-targets respond AS MUCH as curated targets ⇒
  the signal is GENERIC (pleiotropy / residual structure / weak-instrument noise), **not miRNA-mediated.**

## Multi-SNP cross-fitted instrument (`multisnp.py`) — the "stronger instrument" test (user-asked, in lieu of imputation)
5-fold cross-fitted Ridge(miRNA|PCs ~ all cis-SNPs) → OOF allelic score (unbiased, no winner's curse).
- **Honest OOF first stage: only 10 arms F>10, 1 F>20, 0 F>30** — the single-SNP "38" was winner's-curse
  inflation; true typed-SNP cis-instrument strength is LOW.
- **Specificity STILL null at strong instruments** (4 testable arms): curated −0.010 vs control +0.002,
  Δ=−0.012, p=0.42, repression-dir in 1/4. ⇒ Stronger instruments do NOT rescue specificity.

## VERDICT
Germline cis-eQTL is genuinely exogenous but at achievable power (directly-typed SNP6, n=1,075) the
instruments are too WEAK (10 arms honest-F>10) and NON-SPECIFIC to validate the curated edges. Aggregating
SNPs (multi-SNP score) does not help ⇒ the ceiling is weak miRNA cis-heritability + exclusion/pleiotropy
(cis-SNPs still hit host genes / clusters / paralogs — germline fixes EXOGENEITY, not EXCLUSION), not
instrument aggregation. **Edge existence still rests on ONE observational line.**
⬜ ONE untested extension: full genotype IMPUTATION (denser SNPs → possibly stronger instruments). Unlikely
to change the verdict given the low OOF-F AND the specificity null, but not run (needs a phasing+imputation
pipeline: tools + reference panel + allele annotation, none on disk). Infrastructure here is imputation-ready.
