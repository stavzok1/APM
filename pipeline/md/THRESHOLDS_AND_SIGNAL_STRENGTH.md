# Thresholds and “weak / strong” signal rules

This document explains **where numeric thresholds live** (`pipeline/config.py` → `ThresholdConfig`)
and how **weak vs strong** (and related) classifications are applied in the pipeline.

All values below refer to the **defaults** in `ThresholdConfig` unless you override them when calling APIs.

---

## 1. SCREEN experimental & computational links (`screen_links.py`)

**Config**: `THRESHOLDS.screen_weak_quantile`, `THRESHOLDS.screen_strong_quantile`, `THRESHOLDS.intact_hic_pvalue_threshold`

After gene filtering, each row’s `score` is compared to empirical quantiles computed **per `assay_type` only** (`groupby("assay_type")["score"].quantile(...)`). All biosamples that share an assay therefore use the **same** `q_weak` / `q_strong` cutoffs for that assay in a given run.

- **Strong**: `score >= q_strong` where `q_strong` is the **`screen_strong_quantile`** (default **0.95**) quantile of scores for that assay.
- **Weak**: `score >= q_weak` **and not strong**, where `q_weak` is the **`screen_weak_quantile`** (default **0.90**) quantile.

Row-level flags are then aggregated to one row per `(cCRE, gene, biosample, assay)`; `any_strong` / `any_weak` use **max** across duplicate rows, and `classify_strength` in `pipeline/utils.py` maps flags to the string **`"strong"`**, **`"weak"`**, or **`"none"`** (strong wins over weak).

**Hi-C “intact” links** (when Hi-C p-values exist): for `assay_type == "Intact-HiC"`, **strong** additionally requires  
`p_value <= intact_hic_pvalue_threshold` (default **1e-6**).

**Interpretation**: “weak/strong” is **relative ranking within each assay’s score distribution** in the extracted subset, not an absolute physical unit, and is **not** recomputed separately per biosample.

---

## 2. ABC enhancer–gene links (`abc_links.py`)

**Config**: `THRESHOLDS.abc_present_threshold`, `THRESHOLDS.abc_strong_threshold`

For aggregated ABC rows (`ABC_score` is the primary scalar):

- **Present**: `ABC_score >= abc_present_threshold` (default **0.015**).
- **Strong**: `ABC_score >= abc_strong_threshold` (default **0.05**) **and** `element_class != "promoter"`  
  (promoters are treated as present but not “strong” enhancer evidence).

**Rank within gene**: `rank_within_gene = ABC_score / max(ABC_score)` within `(gene_name, cell_type)`.

---

## 3. TargetScan miRNA weights (`mirna_targets` integration)

**Config**: `THRESHOLDS.mirna_weight_threshold`, `THRESHOLDS.mirna_top_n`

- **`mirna_weight_threshold`** (default **-0.2**): context-score cutoff used when filtering / prioritizing predicted targets (more negative = stricter / fewer retained targets in typical usage).
- **`mirna_top_n`** (default **200**): cap on how many top targets are retained per miRNA when building summaries.

These are **not** “weak/strong” labels; they are **cutoffs on model scores**.

---

## 4. Methylation hyper / hypo thresholds

**Config**: `THRESHOLDS.meth_hypermeth_threshold`, `THRESHOLDS.meth_hypometh_threshold`,  
`meth_intermediate_low`, `meth_intermediate_high`

Beta values are on **[0, 1]**:

- **Hypermethylated**: beta **>** `meth_hypermeth_threshold` (default **0.7**).
- **Hypomethylated**: beta **<** `meth_hypometh_threshold` (default **0.3**).
- **Intermediate band**: between **0.3** and **0.7** (inclusive boundaries per column usage in code).

**QC**: `meth_min_valid_probe_pct`, `meth_detection_pval_threshold` gate low-quality samples/probes.

---

## 5. SNV somatic confidence (`SNV` filters)

**Config**: `snv_min_tumor_vaf`, `snv_max_normal_vaf`, `snv_min_tlod`, `snv_min_popaf`, depths, `snv_require_pass`, etc.

These define **high-confidence somatic** variants (VAF / depth / population AF / TLOD). They are **binary pass/fail**, not graded weak/strong.

---

## 6. SV strict vs lenient filters

**Config**: `sv_min_tumor_alt`, `sv_max_normal_alt`, `sv_min_somatic_score`, lenient counterparts, windows (`sv_gene_window_bp`, …).

Strict sets require **higher tumor support** and **lower normal support**; lenient thresholds relax these. Again, **binary** inclusion, not graded signal strength.

---

## 7. FIMO motif scanning (SV)

**Config**: `fimo_pvalue_threshold` (default **1e-4**)

Motif hits with p-value **above** this threshold are typically discarded downstream. This is a **hard cutoff**, not a weak/strong tier.

---

## 8. ChIP-seq overlap summaries (`chip_hits` on cCREs)

ChIP peaks carry `score_norm` (ENCODE signalValue normalized in loader; may be `None` for CHIP_ATLAS).

There is **no global weak/strong quantile** applied in `chip_hits` today: the nested dict stores **per-TF / cell line / source** aggregates (`n_peaks`, `max_overlap_bp`, `max_score_norm`, `mean_score_norm`). Interpret **strength** by comparing `score_norm` or overlap mass **relative to other cCREs** or **within a TF**.

---

## 9. Distance tiers (cCRE ↔ gene)

**Config**: `THRESHOLDS.tier_edges_bp`, `tier_labels`

Genes are binned into distance bands (e.g. 0–100kb, …) relative to cCREs. These are **geometry buckets**, not signal strength.

---

## Quick reference: where to change defaults

Edit `ThresholdConfig` in `pipeline/config.py`, or pass explicit parameters into the relevant builder functions in:

- `pipeline/evidence/screen_links.py`
- `pipeline/evidence/abc_links.py`
- `pipeline/Methylation/*`
- `pipeline/SNV/*`, `pipeline/SV/*`
