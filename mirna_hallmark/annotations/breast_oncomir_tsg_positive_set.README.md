# Breast oncomiR–TSG positive set (Δρ validation)

Hand-curated positive control for the aggregate force-vs-abundance **tumor↔GTEx Δρ** test
(`method_dev/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md` §4, §9.3): tumor-suppressor genes that
**independent functional literature** establishes as directly repressed by oncomiRs in
breast (or pan-cancer with breast relevance) — luciferase + phenotype-rescue grade, not
mere database membership.

## Non-circularity

The **label** ("known oncomiR-repressed TSG") comes from the functional literature, **not**
from any pressure/coupling/edge quantity we measure. These genes are (expectedly) also in
miRTarBase, but the force-vs-abundance comparison is run on the *same* edge set for both
predictors — so the test is whether the **force** Δρ recovers these genes *more evidently
than abundance* does, a within-construction contrast that the label cannot bias. Do **not**
re-derive the positive set from `gene_force_coupling.tsv` or any Δρ output.

## Composition (checked 2026-06-29)

- 22 curated; **18 scored** (have HE edges + coupling in the Hallmark universe).
- **Dropped (not in Hallmark universe / no HE edge):** HOXD10, DICER1, RECK, SMAD4 —
  curated and real, but unusable here; revisit if the gene universe is expanded (TarBase v9).
- Positives sit in the **multi-regulator regime** (PTEN n_reg 87, CDKN1A 43, TGFBR2 28,
  FOXO1 27, BCL2L11 25 …) — i.e. the n_eff≥2 tail where the construction can act. One
  single-driver exception: **TPM1** (n_reg 1, miR-21 only) — force ≡ abundance for it by
  Spearman scale-invariance; keep but expect no force/abundance difference.

## Negative set

**miRNA-poor genes** = scored genes with `n_regulators ≤ 1`, not in the positive set: **648
genes** (median regulation across scored genes is n_reg = 2, so ≤1 is a clean floor). The
claim to test is that the force Δρ **separates** the positive set from this negative set —
stronger than "the positive set shows an effect."

`landmark` column = first-author/year anchors for the well-established relationships; entries
marked "(verify)" need PMID confirmation before any external write-up.
