#############################
## PAM50 recomputation for TCGA-BRCA
## Matrix input: genes = HGNC symbols
#############################

## 0. Load packages and PAM50 model ----
if (!requireNamespace("genefu", quietly = TRUE)) {
  stop("Package 'genefu' is not installed. In R: BiocManager::install('genefu')")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Package 'org.Hs.eg.db' is not installed. In R: BiocManager::install('org.Hs.eg.db')")
}

library(genefu)
library(org.Hs.eg.db)

data(pam50.robust)  # PAM50 model used by molecular.subtyping

## 1. File paths (edit if needed) ----
expr_file     <- "/home/stavz/masters/gdc/RNAexp/BRCA_log2_expression"
clinical_file <- "/home/stavz/masters/gdc/RNAexp/BRCA_clinical"
output_file   <- "/home/stavz/masters/gdc/RNAexp/BRCA_PAM50_merged.tsv"

## 2. Load expression matrix (genes x samples, rownames = symbols) ----
cat("Reading expression matrix from:", expr_file, "\n")
expr_df <- read.delim(expr_file,
                      header = TRUE,
                      check.names = FALSE,
                      stringsAsFactors = FALSE)

## Assume first column is gene symbol:
gene_col <- colnames(expr_df)[1]
cat("Assuming first column is gene symbol:", gene_col, "\n")

rownames(expr_df) <- expr_df[[gene_col]]
expr_df[[gene_col]] <- NULL

expr_mat <- as.matrix(expr_df)

cat("Expression matrix dimensions (genes x samples):",
    paste(dim(expr_mat), collapse = " x "), "\n")

## 3. Map SYMBOL -> Entrez IDs ----
gene_symbols <- rownames(expr_mat)

cat("Mapping SYMBOL -> Entrez IDs using org.Hs.eg.db ...\n")
entrez <- mapIds(org.Hs.eg.db,
                 keys      = gene_symbols,
                 keytype   = "SYMBOL",
                 column    = "ENTREZID",
                 multiVals = "first")

keep <- !is.na(entrez)
cat("Genes with Entrez ID:", sum(keep), "out of", length(entrez), "\n")

expr_mat <- expr_mat[keep, , drop = FALSE]
entrez   <- entrez[keep]

## Collapse duplicates by Entrez (mean expression) ----
cat("Collapsing duplicated Entrez IDs by mean expression...\n")
f <- factor(entrez)
expr_entrez <- rowsum(expr_mat, group = f) / as.vector(table(f))

cat("After collapsing, matrix dimensions (genes x samples):",
    paste(dim(expr_entrez), collapse = " x "), "\n")

## 4. Prepare data / annot for molecular.subtyping ----
## molecular.subtyping expects:
##  - data: samples x "probes"
##  - annot: data.frame with at least EntrezGene.ID column

ddata <- t(expr_entrez)  # samples x genes (Entrez IDs)
col_ids <- colnames(ddata)

dannot <- data.frame(
  probe         = col_ids,     # arbitrary, but consistent
  EntrezGene.ID = col_ids,     # TRUE Entrez gene IDs
  row.names     = col_ids,
  stringsAsFactors = FALSE
)

## 5. Run PAM50 classification ----
cat("centering genes around median...\n")

expr_entrez_centered <- t(apply(expr_entrez, 1, function(x) x - median(x, na.rm=TRUE)))
ddata <- t(expr_entrez_centered)

cat("Running PAM50 molecular subtyping with genefu...\n")

pam50_res <- molecular.subtyping(
  sbt.model  = "pam50",
  data       = ddata,
  annot      = dannot,
  do.mapping = TRUE,
  verbose    = TRUE
)

cat("PAM50 subtype counts (recomputed):\n")
print(table(pam50_res$subtype))

pam50_df <- data.frame(
  sample_id        = rownames(ddata),
  PAM50_recomputed = as.character(pam50_res$subtype),
  stringsAsFactors = FALSE
)

## 6. Load Xena clinical matrix + PAM50Call_RNAseq ----
cat("Reading Xena clinical matrix from:", clinical_file, "\n")
clin_df <- read.delim(clinical_file,
                      header = TRUE,
                      check.names = FALSE,
                      stringsAsFactors = FALSE)

if (!"sampleID" %in% colnames(clin_df)) {
  stop("Clinical file does not contain a 'sampleID' column. Check column names.")
}
if (!"PAM50Call_RNAseq" %in% colnames(clin_df)) {
  warning("Clinical file does not contain 'PAM50Call_RNAseq'. These values will be NA.")
}

clin_pam <- clin_df[, c("sampleID", intersect("PAM50Call_RNAseq", colnames(clin_df)))]

## 7. Harmonize sample IDs (cut to 16 chars) ----
short_id <- function(x) substr(x, 1, 16)

pam50_df$sample_short <- short_id(pam50_df$sample_id)
clin_pam$sample_short <- short_id(clin_pam$sampleID)

## 8. Merge recomputed PAM50 with Xena PAM50Call_RNAseq ----
merged_pam <- merge(
  pam50_df,
  clin_pam,
  by = "sample_short",
  all.x = TRUE,
  suffixes = c("", "_clin")
)

if (!"PAM50Call_RNAseq" %in% colnames(merged_pam)) {
  merged_pam$PAM50Call_RNAseq <- NA
}

merged_pam$PAM50_final <- ifelse(
  is.na(merged_pam$PAM50Call_RNAseq),
  merged_pam$PAM50_recomputed,
  merged_pam$PAM50Call_RNAseq
)

cat("Final PAM50 subtype counts (PAM50_final):\n")
print(table(merged_pam$PAM50_final))

## 9. Save result ----
cat("Writing merged PAM50 table to:", output_file, "\n")
write.table(
  merged_pam,
  file      = output_file,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

cat("Done.\n")
