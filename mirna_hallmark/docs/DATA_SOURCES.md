# Data sources — `mirna_hallmark`

Every dataset the subproject reads, **what it is**, and **what its preprocessing is based on**.
Grouped by lane. Paths are relative to repo root (`/sci/labs/michall/stavzok/APM`). The single
source of truth for paths is `mirna_hallmark/config.py` (+ parent `pipeline/config.py:PATHS`);
this doc is the human-readable provenance map. Newly **ingested 2026-06-26** sets are flagged ⊕.
**2026-07-04 provenance sweep (flagged ⊗):** added TarBase v9, Manakov chimeric eCLIP, SCAN-B and
full METABRIC-mRNA cohorts, arm-resolved isomiR, and CPTAC phospho/acetyl layers; **corrected the
miRGeneDB entry** (it is a guide/passenger *sequence* catalog only — no experimental context).

> Cross-platform rule (applies throughout): external cohorts are integrated with the TCGA spine
> by **rank / direction**, never absolute scale. miRNA arm names are bridged to the modern
> `-5p/-3p` convention (miRBase-v16 array names and bare names resolved to guide/passenger);
> mature↔MIMAT via `data/miRNA/locus_mimat_weight_reference.tsv` + `mimat_to_arm`.

## 0. Technology, annotation build & how the cross-platform/nomenclature gap was overcome

One row per dataset: the exact assay, the annotation build it ships, the integration challenge,
and ✓ how it was resolved (or ✗ if it remains a caveat).

| Dataset | Assay / platform | Annotation build | Challenge | Overcome |
|---|---|---|---|---|
| TCGA-BRCA miRNA | Illumina miRNA-seq (BCGSC) | mature arms by **MIMAT** (miRBase) | — (reference scale) | ✓ spine reference; log2(RPM+1) |
| **TCGA-BRCA isomiR** ⊗ | BCGSC **isoform** expression quant (raw `.isoforms.quantification.txt`) | per-isoform `isoform_coords` (hg38 5′-3′) + `miRNA_region` (mature,MIMAT) | 5′-isomiR seed shifts vs paralog loci (resolved ±6nt window) | ✓ **full cohort fetched** via GDC API (1207-file manifest `data/miRNA/gdc_mirna_isoform_download_manifest.tsv` → `data/miRNA/GDCdata_isoform/`, `scripts/analysis/fetch_tcga_isomir.py`); powers `learned/within_family.seed_shift` (§8 collapse QC, MH-96); + CPTAC `data/CPTAC/cptac2_mirna_isoform/` cross-cohort. NB the *processed* `TCGA-BRCA.mirna_isoform_rpm_gdc.tsv` is MIMAT-aggregated (no 5′ coord) → the RAW files are required for seed QC |
| TCGA-BRCA mRNA | Illumina RNA-seq | gene symbol / Ensembl | — | ✓ spine reference |
| TCGA-BRCA CNV | Affy SNP6 → **ASCAT3** segments | gene/locus overlap | broad segments | ✓ overlap→per-gene/locus CN; focal claims limited (✗ caveat) |
| GTEx breast | Illumina TPM (miRNA + gene) v10 | MIMAT / Ensembl, genome **GRCh38** | cross-platform scale vs TCGA | ✓ **quantile-normalised onto TCGA scale**; arms absent→floor-0 |
| miRTarBase | curated DB | mature-name (current) | mixed evidence tiers | ✓ evidence_score tiers; Hallmark-scoped |
| TargetScan / ENCORI / miRGeneDB | seq prediction / CLIP / curation | seed-family / mature | name drift | ✓ seed-family map; guide/passenger 3-tier resolution |
| **TarBase v9** ⊗ | curated + high-throughput (CLIP/CLASH/degradome) DB | mature-name + **MIMAT**; Ensembl gene/transcript; **3′UTR coords** | huge (4.72M human rows); mixed methods; directional | ✓ MIMAT/arm bridge; **carries `experimental_method`, `regulation` direction, tissue, cell_line, confidence, PubMed id, site start/end** — ingested (`data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz`), **not yet wired into the spine** |
| **Manakov chimeric eCLIP** ⊗ | AGO chimeric eCLIP (direct miRNA–target duplex reads) | mature arm + target coords | per-run SRR tables; HEK293T context | ✓ direct-duplex physical support for orphan edges (`data/external_cache/manakov_chimeric/`); **used ad hoc for orphan validation, not canonicalised** |
| CPTAC-2 prospective | TMT-11 MS proteome + GDC miRNA-seq | gene symbol; **MIMAT** | MS batch | ✓ self-contained pressure; **literal TMT plex (PDC000120)** covariate |
| CPTAC TCGA-105 | iTRAQ-4 MS proteome | gene symbol | MS batch | ✓ "Experiment Order" plex covariate |
| Buffa (GSE22216/22219) | Illumina bead miRNA + Affy mRNA arrays | **bare** miRNA names; Affy probes | bare-name drift; pairing | ✓ bare-stem→guide map; **paired by clinical fingerprint** |
| Plasma/serum EV cohorts (GSE270497/197020/255660/301416/241784/785) | EV-seq / exosome-seq / qPCR / NanoString | mature names (mixed) | heterogeneous platforms | ✓ per-cohort normalize; rank-bridge; binary contrasts |
| **GSE59247** ⊕ miRNA | Agilent **GPL15019** "V16.0" array | array labelled V16.0 but ships **miRBase r14 IDs** (`hsa-let-7a`, `hsa-miR-29a*`, `_v16.0` tags) | old/bare nomenclature; star arms | ✓ **v16→modern guide/passenger bridge** (exact→abundance-resolution; `bridge_method` per arm; 663 arms) |
| **GSE59246** ⊕ mRNA | Agilent **GPL13607** 8×60K *Feature-Number* | numeric feature id (no symbols in matrix) | feature id ≠ gene | ✓ persisted **`GPL13607_feature_to_gene.tsv.gz`** (GEO platform table; 58,717→34,809 genes); probes→gene mean |
| **GSE297448** ⊕ EV miRNA | small-RNA-seq (sEV), MCF10A isogenic | modern `-5p/-3p` names | **global EV-cargo increase** (95% of miRNA rise) | ✓ **CLR compositional EV-share** (removes per-sample loading); ✓ exact cell-line-token stage map |
| **GSE93740** ⊕ cellular miRNA | array, MCF10A isogenic (cellular) | modern names | linking to EV lineage | ✓ same MCF10A lineage as GSE297448 → **direct** cellular-vs-EV export |
| **GSE162670** ⊕ LCM miRNA | Affymetrix RMA, FFPE **LCM** TNBC | modern names; morphology labels in **sample name suffix** | messy GEO-template matrix; morphology not compartment | ✓ skip `#`-comment header → `ID_REF` row; suffix→morphology (N/cis/T); ✗ no epi/stroma compartment (used as composition-reduced progression) |
| GSE297447 ⊕ EV mRNA | RNA-seq (sEV) | gene symbol | (staged) | — staged for recipient-delivery follow-up |
| GSE22216 ⊕ METABRIC miRNA | Illumina bead array | bare names | bare-name drift | ✓ via `buffa_validation` bridge |
| **Risom 2022** ⊕ MIBI-TOF | multiplexed ion-beam imaging, 37 protein markers | single-cell tables, cell lineage/phenotype + tissue type | **protein, no miRNA**; cross-modality | ✓ effector readout (miR-29 targets COLI/SMA/FAP + myoep CK5/P63) + CAF/myoepithelial compartments; ✗ cannot give miR-29c (modality gap → §5.8) |
| **GSE37527** ⊕ CAF/NF miRNA | Agilent **GPL14767** miRNA V16 (Feature-Number) | numeric feature id; miRBase-v16 names (resolves a/b/c) | feature id ≠ miRNA name | ✓ persisted `GPL14767_feature_to_mirna.tsv.gz` (GEO platform table; 15,024→939 miRNAs); the **only non-circular** compartment-resolved breast miRNA measurement (CAF vs NF) |

---

## 1. Core TCGA-BRCA spine (the model is built on this)

| Dataset | What | File(s) | Preprocessing basis |
|---|---|---|---|
| TCGA-BRCA miRNA | mature-arm expression, **2,236 arms × 1,207 samples (1,096 tumour, 104 NAT, 7 met)** | `data/miRNA/…` (RPM/MIMAT) | TCGA miRNA-seq RPM; collapsed to mature arms by MIMAT; log2(RPM+1) for coupling; detectable subset far smaller (floor ~RPM≥10); Normal-like excluded for clinical cohort |
| TCGA-BRCA isomiR ⊗ | arm-resolved **isomiR** counts/RPM (5′/3′ end heterogeneity) | `data/miRNA/TCGA-BRCA.mirna_isoform_*` | GDC isoform quantification → MIMAT; sub-arm resolution available but collapsed to arm for the spine |
| TCGA-BRCA mRNA | gene expression (the target side) | parent `PATHS.rna_matrix` | TCGA RNA-seq; log2; gene-symbol keyed |
| TCGA-BRCA CNV | gene + miRNA-locus copy number | `data/CNV_TCGA/…`, `mirna_locus_cnv` cache | ASCAT3 segments → gene/locus overlap → per-sample CN; classified gain/loss/amp/deep_del |
| TCGA clinical / strata | PAM50, ER/PR/HER2, CPE purity, HRD, aneuploidy, immune subtype, batch | `annotations/BRCA_clinical_immune_unified.tsv`, Thorsson/Taylor | unified clinical table; partial-correlation confounders = **CPE + thornsson_hrd_score** (+ TCGA plate batch `plate_both`) |
| GTEx breast (true-healthy) | small-RNA + gene TPM, normal breast | `data/GTEx/{miRNA_TPM_matrix_v10,smallRNA_annotation_v10,gene_tpm_v10_breast,GTEx_v10_SampleAttributesDS}` | GTEx v10 TPM; **quantile-normalised onto TCGA scale** for the healthy anchor (arms absent in GTEx floored to 0 = tumour-acquired); GTEx is the healthy anchor, **never NAT** (NAT is field-effect) |
| **CIBERSORTx cell-type fractions** (the composition block `C`) | per-sample cell-type fractions for **all three states** — tumour, NAT, GTEx | `data/external/cibersortx/` (**see its `README.md` — full inventory + usage rules**); re-exported to `mirna_hallmark/output/brca_deconvolution/tcga_cibersortx_fractions.tsv` | CIBERSORTx Fractions, S-mode batch correction, **Wu-major** scRNA ref (GSE176078, 9 lineages). Wired via `learned/data._DECONV_COLS` = the **8 non-malignant** columns (`Cancer Epithelial` held out — it feeds `_malignant_prolif`). **Use:** tumour=`tcga_tumor_nat_pooled` (1172) · NAT=`tcga_nat_alone` (113) · GTEx=`gtex_wu_major` (514). **Established 2026-07-12:** β is panel-**invariant** (ρ=0.94 across entirely different atlases) but *no-C* halves \|β\| (ρ=0.58) ⇒ having a C matters ~10× more than which; **zero new runs needed**. HBCA runs = second-reference sensitivity; LM22/Thorsson = tumour-only **validators**, never covariates. **(→ `LEARNED_MODEL_STATE_CHANNEL_PLAN.md` §1/§3/§4)** |

States used across the spine: **GTEx (true-healthy) → NAT (field-effect) → tumour**; the DCIS
work below adds an **in-situ / invasive** axis on top (see lane 4).

---

## 2. Edge universe + non-circular priors (evidence, not expression)

| Source | What | File / cache | Preprocessing basis |
|---|---|---|---|
| miRTarBase | curated miRNA→gene interactions → Hallmark-scoped edges | `data/miRNA/mirtar.csv` → `output/edges/mirna_hallmark_edges.tsv.gz` | filtered to MSigDB-Hallmark target genes; evidence tiers (functional-MTI / reporter / protein); `evidence_score`, `high_evidence` |
| MSigDB Hallmark | 50 hallmark gene sets | `annotations/GSEA/h.all.v2026.1.Hs.txt` | parsed by `hallmark_sets.py` |
| TargetScan | seed-family map + predicted (orphan) edges | `MIR_FAMILY_INFO` | seed-family collapse for dependence-aware FDR; orphan = TargetScan∪ENCORI minus miRTarBase-functional |
| ENCORI / starBase | CLIP-supported edges | `data/external_cache/encori/` | evidence overlay + orphan universe; depth boost (α=0.5) |
| miRGeneDB ⊗ | **guide/passenger mature *sequences* only** — 559 guide (`hsa_guide_mat.fas`) + 468 star (`hsa_star.fas`) FASTA | `data/external_cache/mirgenedb/` | curated arm-assignment / validated-arm universe. **Correction:** a sequence catalog with **no experimental / assay context** — use for guide-vs-passenger resolution and arm breadth, not as an evidence-context source |
| **TarBase v9** ⊗ | curated + HT (CLIP/CLASH/degradome) validated interactions, method- and direction-labelled | `data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz` (4.72M human rows) | columns: `experimental_method`, `regulation` (Positive/Negative), tissue, cell_line, confidence, PubMed id, **3′UTR site start/end**, microT score. A method-graded, directional, partly site-resolved prior source. **Ingested, not yet wired** |
| **Manakov chimeric eCLIP** ⊗ | direct AGO **chimeric duplex** reads (physical miRNA–target pairing) | `data/external_cache/manakov_chimeric/*.tsv` (per-SRR) | HEK293T; validated ~60% of high-value orphan edges. **Used ad hoc for orphan validation, not canonicalised** |
| PolyASite + APAatlas + Xia 2014 | APA / 3′UTR diagnostics | `data/external_cache/{polyasite,apaatlas}/`, `data/miRNA/XIA_APA.xlsx` | gene-level APA overlay (diagnostic; cancels in re-ranking); APAatlas = GTEx-normal breast, Xia = TCGA tumour ΔPDUI |
| PubMed breast-context | literature prior per edge (non-circular) | `data/external_cache/pubmed/` | breast-context PMID counts per edge |

---

## 3. Independent validation cohorts (replication, not training)

| Cohort | What | Path / loader | Preprocessing basis |
|---|---|---|---|
| CPTAC-2 prospective | independent breast proteome (**+ phosphoproteome + acetylome** ⊗) + miRNA | `data/CPTAC/…`, `scripts/cptac/*` | self-contained pressure from CPTAC-2 GDC miRNA-seq; protein anti-corr / residual; confounders purity+CIN; **batch = literal TMT plex (PDC000120)**; phospho/acetyl layers available for effector/PTM readouts |
| CPTAC TCGA-105 | iTRAQ proteome (same-patient) | `data/CPTAC/…` | iTRAQ; batch = "Experiment Order" plex from Sample-Mapping xlsx |
| Buffa 2011 (GSE22216 miRNA + GSE22219 mRNA) | 207 paired primary tumours, no TCGA overlap | `data/external_cache/geo/`, `buffa_validation.py` | microarray; **paired by clinical fingerprint** (GEO title index is NOT the pairing); adj prolif (E2F/G2M) + ER |
| Plasma/serum EV miRNA cohorts | GSE270497, GSE197020, GSE255660, GSE301416, GSE241784/785 | `config.EV_MIRNA_*` paths; `ev_mirna_multi_cohort.py`, `ev_mirna_followup.py`, `ev_mirna_gse255660_screen.py` | per-cohort normalized (log2 TPM/CPM/expr/dcq); binary contrasts (cancer vs normal/benign, carrier vs wt, Basal vs LumA); rank-bridged to TCGA pressure |
| **SCAN-B** ⊗ | large independent breast **mRNA** (RNA-seq) — GSE96058 (3,273 samples + 136 reps) + GSE60789 | `data/SCAN-B/` | external target-side replication; rank/direction-bridged, never absolute scale |
| **METABRIC (full mRNA)** ⊗ | ~1,900-tumour breast **mRNA microarray** + CNA + clinical (distinct from the GSE22216 miRNA subset) | `data/METABRIC/{data_mrna_illumina_microarray,data_cna,data_clinical_*}.txt` | independent target-side + outcome cohort (outcome modules); Illumina array |

---

## 4. Pre-malignant / DCIS / progression corpus ⊕ (ingested 2026-06-26)

Staged under `data/external/dcis_geo/` (manifest: `data/external/dcis_geo/README.md`).
Downloaded from NCBI GEO FTP; all public. These add the **in-situ → invasive** axis and the
**EV × progression** lane. See discovery registry **MH-48/49/50**.

| GSE | What | n | usable file | loader / module | preprocessing basis |
|---|---|--:|---|---|---|
| **GSE59247** ⊕ | DCIS+IBC **miRNA** (Agilent **GPL15019**, miRBase-v16) | 48→**40** | `GSE59247_series_matrix.txt.gz` | `dcis_geo_loader.load_gse59247` | hsa mature arms; **state from curated `tissue type`, NOT title** (GSM1431408 mislabelled); 4 technical-replicate pairs collapsed; Ambion controls dropped; **no matched synchronous DCIS+IBC pairs** |
| **GSE59246** ⊕ | matched **mRNA**, same lesions (Agilent **GPL13607**, *Feature-Number*) | 105→**102** | `GSE59246_series_matrix.txt.gz` | `dcis_geo_loader.load_gse59246` | numeric feature id → `GeneName` via persisted **`GPL13607_feature_to_gene.tsv.gz`** (fetched from GEO platform table, 58,717 probes→34,809 genes); probes collapsed per gene by mean; pairs to GSE59247 by `rep_group` (**40 paired samples**) |
| **GSE297448** ⊕ | **EV (exosomal) miRNA-seq**, MCF10A **isogenic** series | 12 | `GSE297448_Normalized_microRNA_readcounts_All_samples.txt.gz` | `dcis_ev_progression.py` | **isogenic** MCF10A→AT1→DCIS→CA1a (Normal→benign→DCIS→invasive), 3 reps/stage; modern arm names; normalized readcounts carry a **global EV-cargo increase** → analysis uses **CLR (compositional EV-share)**; ordered trajectory legitimate *only* because isogenic |
| GSE297447 ⊕ | EV mRNA-seq companion to 297448 | 12 | `GSE297447_normalized_mRNA_gene_counts_all_samples.txt.gz` | (staged) | EV target side; normalized gene counts |
| GSE93740 ⊕ | **cellular** MCF10A isogenic progression (P→AT1→NeoT→DCIS→CA1) | 18 | `GSE93740_bedrosian-mirna.exact4.txt.gz` | `dcis_mir29c_gse93740.py` | cellular counterpart of GSE297448 EVs → **direct export test** + miR-29c validation (MH-51) |
| GSE255660 ⊕ | paired-miRNA mammogram cohort (also in lane 3) | 319 | `GSE255660_series_matrix.txt.gz` | `ev_mirna_gse255660_screen.py` | benign/malignant discrimination |
| GSE162670 ⊕ | LCM miRNA by tissue morphology (TNBC) | 74 | `GSE162670_RMA_normalized_matrix.txt.gz` | `dcis_lcm_compartment_gse162670.py` | composition-reduced (microdissected) progression tiebreaker (MH-51) |
| GSE22216 ⊕ | METABRIC invasive miRNA (= Buffa miRNA arm) | 210 | `GSE22216_series_matrix.txt.gz` | `buffa_validation.py` | independent invasive cohort |
| **Risom 2022** ⊕ | **single-cell spatial PROTEOMICS** (MIBI-TOF, 37-plex), matched normal/DCIS/IBC | 69,672 cells | `risom2022_dcis_mibi/risom2022_single_cell.csv` (+ `tissue_features/Table_S*.csv`) | `dcis_htan_compartment.py` | **effector/compartment proof** of the invasion programs (myoepithelial loss; CAF activation + ECM); MH-52. Mendeley `d87vg86zd8` (placed manually). NB protein, **not miRNA** |
| **GSE37527** ⊕ | **paired CAF vs NF miRNA** (breast, primary cultures), the non-circular compartment test | 6 pairs (12) | `caf_nf_gse37527/GSE37527_series_matrix.txt.gz` (+ `GPL14767_feature_to_mirna.tsv.gz`) | `dcis_caf_nf_gse37527.py` | **directly measures miR-29c per compartment** (MH-54): miR-29c NOT down in CAFs (overturns the CAF hypothesis); miR-143/145/126/205 are CAF-lost. GEO FTP. **Cultured CAF → culture-drift caveat** |
| **GSE196354** ⊕ | **paired CAF vs NF mRNA** (breast, RNA-seq), direct TGFβ test | 7 pairs (14) | `caf_nf_mrna_gse196354/GSE196354_NFs_CAFs_RNAseq_FPKM.txt.gz` | `dcis_caf_tgfb_gse196354.py` | TGFβ/ECM mRNA test (MH-54 addendum): cultured CAFs **culture-drifted** (weak myCAF program) → TGFβ-desmoplasia rests on in-vivo Risom. GEO FTP |

**Caveats specific to this corpus**
- Modest n (12–48) for several → used to **slot in a state / read direction**, not for within-stage power.
- EV ≠ intracellular pressure: exported miRNA can be inversely related to intracellular repression → EV sets are a **biomarker/export** question kept separate from the tissue-pressure mechanism, and require compositional (CLR) analysis (GSE297448).
- `GSE297448` is an **isogenic cell-line** series (corrects an earlier note that called it patient DCIS→IDC); this is *why* its ordered trajectory is usable.

**Not downloaded** (controlled access, dbGaP/EGA DAR): HTAN MCL Pre-Cancer Atlas (`phs002225`),
HTAN DCIS (TBCRC 038 + RAHBT, with progression outcomes). Upgrades, not prerequisites.

---

## Reproducing the ingests
GEO FTP pulls are documented in `data/external/dcis_geo/README.md`. The GPL13607 feature→gene
map was built once from the GEO platform table and persisted to
`data/external/dcis_geo/GPL13607_feature_to_gene.tsv.gz` (no live fetch at run time).

---

## miRNA–target INTERACTION EVIDENCE — by assay class × tissue × method (2026-07-09)

The miRNA→gene edge evidence, organized as requested: **per data category (assay class) → subcategorized by
breast / non-breast and chimeric / non-chimeric**, with specific (per-source PMID) and pooled counts. Counts are
**Hallmark-gene-scoped** (the model's universe) unless noted. Assay→class mapping and weights are the single source
of truth in `learned/evidence/methods.py` (`classify`, `CLASS_WEIGHT`); ingestion in `learned/evidence/ledger.py`.

### How membership vs weight differ (READ FIRST — I kept conflating these)
- **HE MEMBERSHIP** (`ledger.pooled_he_edges`, ≈5,940 edges) = miRTarBase-HE ∪ TarBase **{reporter, western,
  proteomics}** (low-throughput FUNCTIONAL), non-weak. **Excludes ago_clip, qpcr_rna, chimeric, prediction.**
- **EVIDENCE WEIGHT** `w` (`ledger.edge_weights`) = Σ_class `CLASS_WEIGHT[class]`·log1p(#distinct PMIDs), PMID-
  deduped across miRTarBase+TarBase. **ago_clip DOES contribute to w (0.5)**, just not to HE membership.
- CLASS_WEIGHT: reporter 3.0 · western 2.5 · proteomics 2.5 · qpcr_rna 1.5 · degradome 1.0 · chimeric 1.0 ·
  ago_clip 0.5 · other 0.0.

### Inventory (TarBase v9, Hallmark genes; edges = distinct (arm,gene))

| assay class | tissue | chimeric? | edges | PMIDs | methods | in HE? | weight |
|---|---|---|---|---|---|---|---|
| **ago_clip** | breast | no | 127,721 | 7 | AGO-IP, Biotin-µarray, HITS-CLIP, PAR-CLIP | ✗ | 0.5 |
| **ago_clip** | non-breast | no | 340,166 | 47 | +Biotin-Seq/qPCR, RIP-Seq | ✗ | 0.5 |
| **chimeric** | breast | yes | 40 | 1 | IMPACT-Seq | ✗ | 1.0 |
| **chimeric** | non-breast | yes | 19,260 | 5 | CLASH, Chimeric fragments, qCLASH | ✗ | 1.0 |
| reporter | breast | — | 252 | 97 | Luciferase | ✓ | 3.0 |
| reporter | non-breast | — | 1,232 | 531 | Luciferase | ✓ | 3.0 |
| western | breast | — | 70 | 34 | Western blot | ✓ | 2.5 |
| western | non-breast | — | 655 | 277 | Western blot | ✓ | 2.5 |
| proteomics | non-breast | — | 1,327 | 7 | (p)SILAC, 2D-DIGE | ✓ | 2.5 |
| qpcr_rna | breast | — | 3,690 | 33 | µarray, RNA-Seq, qPCR, sRNA-Seq | ✗ | 1.5 |
| qpcr_rna | non-breast | — | 33,838 | 286 | +Northern, RPF-Seq | ✗ | 1.5 |

**Named breast AGO CLIP sources** (all in TarBase; NOT in HE — CLIP excluded): MCF7 PAR-CLIP **PMID 24398324**
(Farazi; also in miRTarBase-weak; local copy `data/external/breast_parclip_farazi2014/`, 11,516 sites),
MCF7 HITS-CLIP **27150721** (154k pairs), MDA-MB-231 HITS-CLIP **26061048** (144k), BT-474/MCF7 HITS-CLIP **24906430**.
**No breast CHIMERIC exists** (only IMPACT-Seq×40; true CLASH/chimeric-eCLIP breast = none) → L2/arm-level breast
validation impossible; only L1/family via seed-inferred CLIP (partial seed-circularity).

### Local out-of-ledger physical caches (used for L2 arm nomination / loading, NOT for w)
| dataset | assay | tissue | chimeric | pairs | usage |
|---|---|---|---|---|---|
| Manakov chimeric eCLIP | chimeric | HEK293T (non-breast) | yes | ~620k | `chimeric_evidence` (L2), `ago_loading` |
| TarBase qCLASH/CLASH/Chimeric | chimeric | non-breast | yes | 56k/1.8k/11k | `chimeric_evidence` (L2, overlap-resolved) |
| GTEx miRNA | small-RNA | non-breast (pan-GTEx) | — | — | `ago_loading` arm-dominance fallback |
| TargetScan context++ | prediction | — | — | 1.4M | `seq_specificity` (identity percentile, L1) |
| scanMiR K_D (HE scan) | prediction | HE genes only | — | — | `kd.py` legacy HE-restricted scan (occupancy κ) |
| **scanMiR K_D (genome-wide)** ⊗ | prediction | **18,852 MANE genes × 746 detectable arms** | — | ~10.4M | `kd.genome_affinity{,_pct}`; built by `scanmir_genomewide_par.R` (BiocParallel, ~4× faster) → `data/external_cache/scanmir/genomewide_kd{,_new,_disc}.tsv.gz`. **MH-87: OVERTURNS the MH-86 discard** — all-site genome-wide beats context++ at specialist recovery (0.89 vs 0.79 fair, `kd_fair_bench.tsv`); strong-site thresholding HURTS. Consumers: `discovery.py` (biochemical support), `seq_specificity.affinity_percentile_kd` (per-arm specificity), `instrument.family_causal_attribution` |
| **GSE178127 AGO2-RIP (HeLa)** ⊗ | small-RNA loaded (RIP) + total | HeLa (non-breast) | — | 171 arms | AGO-loading cross-context test (MH-87); `data/external_cache/ago_loading_gse178127/hela_loading_enrichment.tsv` |
| **GSE106810 AGO2-RIP (cortex/GBM)** ⊗ | small-RNA loaded (RIP) + input | brain/GBM (non-breast) | — | 642 miRNA (6 pairs) | AGO-loading SAME-ASSAY test (MH-87) → cross-arm loading does NOT transfer (incl. vs HeLa) ⇒ **INERT, dropped**; only within-hairpin passenger-exclusion robust. `data/external_cache/ago_loading_gse106810/` |
| Farazi 2014 MCF7 AGO2-PAR-CLIP | small-RNA loaded + total | **breast MCF7** | — | 158 miRNA | breast loading context (MH-87, `mcf7_loading_enrichment.tsv`); also L1 validation (MH-86). `data/external/breast_parclip_farazi2014/` |

### Which functions dictate USAGE of this data
| function | what it decides | uses |
|---|---|---|
| `learned/evidence/ledger.pooled_he_edges` | **HE membership** (the model's edge universe) | miRTarBase-HE + TarBase{reporter,western,proteomics} |
| `learned/evidence/ledger.edge_weights` | **evidence weight w** (adaptive prior in π/slab) | PMID-deduped ledger × CLASS_WEIGHT |
| `learned/evidence/methods.{classify,CLASS_WEIGHT}` | assay→class + per-class weight | the taxonomy |
| `data_loaders.high_evidence_edges` | miRTarBase `high_evidence` flag | miRTarBase Functional MTI + LT |
| `learned/chimeric_evidence.{evidence,evidence_matrix}` | **L2 arm nomination** (which mature) | Manakov + TarBase chimeric, per-source/pooled |
| `learned/ago_loading.{loading,empirical_arm}` | **guide/passenger + L2** | Manakov AGO-dominance + GTEx |
| `learned/seq_specificity.{affinity_percentile,seq_spec}` | **L1 identity ownership** | TargetScan context++ |
| `learned/genomic_context.classify_he_arms` | **miRNA locus GENOMIC CONTEXT** (sense-coding-host / sense-lncRNA-host / antisense-overlap / intergenic; strand-aware) | GENCODE miRNA loci × host gene bodies |

**miRNA arm→locus localization — CANONICALIZED (2026-07-11).** Two files, ORTHOGONAL roles — do not conflate:
one is the coordinate authority, the other is a quality flag. (Earlier drafts wrongly used the CN file as a
coordinate *fallback* and called it "correctly narrower" — both corrected below.)
- **COORDINATES = `data/gencode.v49.annotation.gtf.csv` (`gene_type=="miRNA"`, 1879 loci, WITH strand) — the SOLE
  authoritative source** for `learned/genomic_context`. Standard, versioned, exhaustive, and the *same coordinate
  system as the host genes* (no cross-source build/strand mismatch). Arm→locus via `_norm_name` + `_arm_gencode_loci`,
  which handles ALL of GENCODE's copy conventions — concatenated (`MIRLET7A1/2/3`), hyphenated (`MIR1-1`), copy-specific
  arm names (`miR-125b-1`→`MIR125B1`), family letters (`MIR517`→`MIR517A/B/C`) — **without** the `MIR93`→`MIR935`
  collision. **Coverage = 714/714 HE arms (100%), single provenance, NO fallback.**
- **QUALITY FLAG = `data/miRNA/cnv_miRNA.csv` (`PATHS.mirna_precursor_loci_csv`, 506 **MirGeneDB** precursors) →
  `genomic_context.mirgenedb` (526/714 HE arms validated).** MirGeneDB is a *curated bona-fide* subset (it excludes
  ~1400 dubious/mis-annotated miRBase entries GENCODE keeps) — so it is a **confidence layer, NOT a coordinate
  authority**. It remains, separately, the **CN-pipeline's** precursor-loci file (`mirna_locus_cnv`/`instrument.locus_cn`;
  the CN universe = these 506 = `_entity_long` locus entities); the 188 non-MirGeneDB HE arms are absent from the CN
  system by that curation choice (a confidence scope, not a bug).
- **MULTI-LOCUS arms** (72/714; 20 context-MIXED): every GENCODE copy is classified; `n_loci` = copy count,
  `context_mixed` = loci disagree on class, `mir_class` = confound-aware representative (majority, ties→coding-host,
  since abundance is the SUM over loci so one host-coupled copy already injects the host program). Dose-weighted
  host-coupled fraction (per-locus expression / `w_norm`) is the Phase-1 refinement.

**⚠ ARM-NAME `.N` NORMALIZATION BUG → UNIVERSE REDEFINITION (GLOBAL, pending next refresh).** The miRNA matrix `X`
stores some arms with a trailing `.N` (`hsa-miR-101-3p.1`, a duplicate-suffix artifact) while `pooled_he_edges` uses
the plain name — so `assemble_gene`'s old `regs = [m for m in ed.miRNA if m in Xall.index]` **silently dropped** them.
**✅ FIXED (`data._arm_name_map`, live in `assemble_gene`)** — transforms: exact / `.N`-strip / case / **suffixless→GUIDE
arm** (old-style `hsa-miR-181a` → most-abundant arm miR-181a-5p, NOT the -3p passenger; abundance-resolved). Recovers
**~18 HE arms (16 expressed) + ~300 edges, incl. canonical miR-101-3p / miR-124-3p / miR-126-3p** ⇒ corrected universe
**HE arms in X 714→730, expressed 480→494, in-X edges 5499→5799**. LIVE for future runs; **persisted outputs stay
old-universe until a refresh regenerates them**. The 90 missing = ~18 recovered + **6 sister-arm-present** (only the
other arm measured — miR-211-3p/216b-3p/…; not substitutable, different seed) + **64 entirely absent** (of which only
miR-21a-5p is important+rescuable = miR-21-5p already in the universe, a mislabeled-edge merge; the rest genuinely
unmeasured — dubious/non-BRCA-expressed). Memory `universe-redefinition-pending-refresh`.

**miRNA PROMOTER architecture (FANTOM5 CAGE, de Rie 2017):** `data/external_cache/fantom5_human_miRNA_promoters.tsv`
(fetched from `fantom.gsc.riken.jp/5/suppl/De_Rie_et_al_2017/vis_viewer/data/human/human.promoters.tsv`). The `promoter`
field = the promoter each miRNA actually uses (`p1@HOSTGENE` = host promoter). **DEFINITIVE own-vs-host-promoter call** for
intragenic miRNAs (211 sense_coding_host arms: **74% host-shared/co-transcribed, 26% own-promoter**) — the gold-standard
substitute for the mature↔host expression correlation, which is confounded by post-transcriptional processing (host-shared
ρ 0.29 ≈ own 0.28). **✅ FOLDED INTO `genomic_context.tsv` (2026-07-11) as the `promoter`+`promoter_class` columns**
(`_fantom5_promoters`/`_classify_promoter`): every HE arm carries its co-transcription call {host_shared / shared_other /
own / unknown}; sense_coding_host = 271 host_shared (80%) / 33 shared_other / 32 own / 1 unknown. **de Rie's LABEL is the
discriminator + a two-part ACCESSION RESOLVER**: (a) `pN@ENST…` → gene via the **GENCODE transcript_id→gene_name map**
(exact); (b) a non-ENST / v49-retired accession (`AK…`) confirmed as the host when its FANTOM5 **hg19 TSS lifted to hg38**
(`data/external_cache/hg19ToHg38.over.chain.gz`, UCSC, via `pyliftover`) lands inside the host body — applied ONLY to
accession labels (miR-21→VMP1, miR-101→JAK1 recovered), NOT coordinate labels (intronic-own peaks stay `own`). **Plus an INDEPENDENT coordinate annotation** (`_promoter_coord`/`_locus_at`): lift TSS→hg38, read GENCODE `promoter_gene`/`promoter_gene_type`/`promoter_feature`(exon/intron)/`promoter_coord_class` — a label-free cross-check (label↔coord agree 78%; disagreements flag positionally-intronic-but-independent promoters); + `promoter_host_tss_kb` = distance to the host's EXACT hg38 TSS (no liftover) — host_shared-label median 0.1kb vs own 2.3kb, exposes alt-promoter co-transcription miR-21@VMP1 130kb). Liftover artifact-free here (0% unmapped); de Rie promoters are CAGE-experimental not distance.

**miRNA-LOCUS CN (extended instrument):** `output/matrices/mirna_locus_cn_gencode.tsv.gz` (`scripts/analysis/build_mirna_locus_cn.py`)
— ASCAT3 allele-specific segments × GENCODE miRNA loci → **1862 loci × 1060 participants** (vs the 506 MirGeneDB loci
of `mirna_locus_cnv`/`instrument.locus_cn`). Gives CN instruments for the newly-localized arms. Host pleiotropy on
these instruments (axis BP) is real only where the host relates to the target — `genomic_context.host_target_relatedness`
(partial Spearman(host,target|CPE,**miRNA**)) → `output/learned/host_target_relatedness.tsv` — now **PER-LOCUS complete**
(2026-07-11): enumerates every distinct coding host across an arm's loci + a `host_locus` column (403→**444 related**; the 6
multi-coding-host arms miR-26a→CTDSP2+CTDSPL, miR-218→SLIT2+SLIT3, … score both).

**PER-LOCUS host map — `genomic_context.locus_host_map()` → `output/learned/locus_context.tsv` (2026-07-11).** Classifies
each **CN `locus_id` (`MI*`)** at its OWN hg38 coordinate, so the CN-exclusion host gate conditions each ACTIVE locus on the
host it *actually* sits in (not the arm's coding-first representative, which for a multi-locus arm can be a different or
SILENT locus). Bridge: **`_precursor_loci`** (`PATHS.mirna_precursor_loci_csv`) maps `MI*` → (chrom,start,end,strand) — the
same ids as `instrument.locus_cn` columns / `arm_loci_map` (verified 506/506) — then `_classify_locus` (GENCODE host bodies).
Output: 506 CN loci = 211 coding-host / 235 lncRNA-host / 32 intergenic / 28 antisense. Consumed by `instrument._locus_coding_host`
→ `exclusion`/`host_pleiotropy`/`dose_attribute_host_downweights`.

### CORUM 5.3 — protein complexes (added 2026-07-12, MH-106)
`data/CORUM/humanComplexes.txt` (+ `humanComplexes.manifest.json`). **5,628 human complexes / 5,150 subunit genes.**
Release 5.3 (2026-04-14), fetched from the CORUM FastAPI backend
`https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=human&file_format=txt`.
⚠ **The host's TLS chain is broken** (cert verification fails in Python *and* WebFetch), so the fetch relaxed
verification and the payload was **VALIDATED post-hoc** (5,628 complexes; ribosome/proteasome/spliceosome present;
RPL5/PSMA1/TP53 present; not an HTML error page — an earlier naive fetch returned a 423-byte SPA shell).
**Used by:** the protein-buffering test on the mRNA→protein propagation slope `a_g` (MH-105/106) — buffering is a
**large-obligate-complex** effect (size ≥31 → `a_g` 0.261 vs 2–3 → 0.407; Spearman −0.079, p=2.6e-6), grounding
CHANNEL_FUSION axis **BM** (protein-complex buffering) empirically.

---

## TF-regulation databases — the shared-inducer test (MH-131, 2026-07-16)

**WHY:** the positive-coupling curated genes (101/482 = 21%, whose curated regulators go **UP** with their
target) survive composition and every proliferation axis. The leading hypothesis was a **SHARED INDUCER** — a TF
driving BOTH the miRNA and its target. Testing that needs **both halves of the circuit**; a co-expression proxy
is NOT a substitute (it failed twice: a sign artifact, and the target's neighbours are **11–14× enriched** for
being targets of the same regulators, so the "module" is not a control).

| source | what it is | coverage | provenance |
|---|---|---|---|
| **CollecTRI** (TF→**gene**) | the successor to DoRothEA (Müller-Dott et al., *NAR* 2023); TRRUST ∪ DoRothEA ∪ ExTRI ∪ HTRIdb ∪ … with per-edge PMIDs | **62,411** interactions, **1,201** TFs × **6,628** genes | `data/external_cache/collectri/collectri_human.tsv` — ⚠ **was already in the repo but UNDOCUMENTED until now**; fetch date/hash unknown, no manifest. **Carries ~NO TF→miRNA edges** (1 apparent row is a UniProt id, `Q3MIR4`, not a miRNA). |
| **TransmiR v2.0** (TF→**miRNA**) ⊗ | literature-curated TF→miRNA regulation with effect direction + PMIDs | **2,678** interactions, **369** TFs × **402** miRNAs; effects: Activation 1,536 / Repression 897 / Regulation 235 | `data/external_cache/transmir/literature_hsa.tsv.gz` + **`manifest.json`** (url, sha256 `9b12aa4a…`, date, row/TF/miRNA counts). **TLS VERIFIED** — standard curl, **no `--insecure`** (contrast: the CORUM 5.3 fetch above, which had to relax verification). Cite: Tong Z, Cui Q, Wang J, Zhou Y. *NAR* 2019;47(D1):D253-D258. |

**Used by:** `scratchpad/shared_inducer.py` → `output/learned/shared_inducer.tsv`.

### ⚠ THE RESULT, AND THE COVERAGE BIAS THAT LIMITS IT

**The shared-inducer hypothesis is NOT CONFIRMED — and it fails in the OPPOSITE direction to the prediction.**
Testable subset **436/482 genes (90%)** with both halves annotated:

| | POSITIVE (n=90) | NEGATIVE (n=346) | |
|---|---|---|---|
| P(≥1 shared TF) | **60.0%** | **74.9%** | OR 0.50, Fisher p=0.998 |
| n shared TFs | 2.23 | 3.31 | MWU p=0.997 |

**⚠ CONFOUNDED BY ANNOTATION COVERAGE, NOT BIOLOGY:** the positive genes' miRNAs carry **37% FEWER annotated
TFs** (23.1 vs 36.8, **p<0.001**) while the TF count on the *gene* side is balanced (24.7 vs 23.3, p=0.53). Fewer
annotated TFs ⇒ fewer chances of a coincidental overlap ⇒ the 60%-vs-75% gap is a **database artifact**.
⇒ **The shared-inducer route is NEITHER confirmed NOR cleanly refuted: TransmiR's 402-miRNA coverage cannot
separate it from its own annotation bias.**

**⛔ A permutation null reported during this run (POSITIVE 60% vs null 39%, p<1e-4) is VOID for this question**:
it permutes the *gene's* TF set while holding the miRNA side fixed, so it asks "is the overlap above chance for
these miRNAs" — not "do POSITIVE genes share more TFs than NEGATIVE ones". The correct null must hold the
**overlap opportunity** fixed on **both** sides.

**⭐ THE INCIDENTAL FINDING WORTH KEEPING:** **60–75% of curated miRNA→target edges have a TF driving BOTH the
miRNA and its target** — and it is *higher* for the normally anti-coupled genes. **Shared-TF feed-forward
co-regulation is the DEFAULT architecture of this edge set, not a pathology of the 21%.** A curated
miRNA→target pair is almost never a clean two-node system.
