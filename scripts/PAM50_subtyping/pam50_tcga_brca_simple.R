#############################
## Simple PAM50 recomputation for TCGA-BRCA
#############################

library(genefu)
library(readr)
library(dplyr)
library(stringr)

expr_file   <- "RNAexp/BRCA_log2_expression"
clinical_file <- "RNAexp/BRCA_clinical"
output_file  <- "RNAexp/BRCA_PAM50_merged"

#############################
## PAM50 recomputation for TCGA-BRCA (base R only)
#############################
data(pam50.robust)

## 2. Load expression matrix ----
cat("Reading expression matrix from:", expr_file, "\n")
expr_df <- read.delim(expr_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

## If your first column is not literally called 'gene_symbol',
## you can either rename it in the file, or change this line:
if (!"gene_symbol" %in% colnames(expr_df)) {
  stop("Could not find a 'gene_symbol' column in expression file. Rename your first column to 'gene_symbol' or change the script.")
}

# Put gene symbols into rownames and drop the column
rownames(expr_df) <- expr_df$gene_symbol
expr_df$gene_symbol <- NULL

expr_mat <- as.matrix(expr_df)

cat("Expression matrix dimensions (genes x samples):",
    paste(dim(expr_mat), collapse = " x "), "\n")

## 3. Prepare data for genefu::molecular.subtyping ----
## molecular.subtyping expects:
##   data  = samples x genes
##   annot = data.frame with 'Gene.Symbol' and rownames == colnames(data)

ddata <- t(expr_mat)   # samples x genes

dannot <- data.frame(
  probe         = colnames(ddata),   # <- what genefu expects
  EntrezGene.ID = NA,                # dummy, OK
  Gene.Symbol   = colnames(ddata),   # assume your columns are HGNC symbols
  row.names     = colnames(ddata),
  stringsAsFactors = FALSE
)
## 4. Run PAM50 classification ----
cat("Running PAM50 molecular subtyping with genefu...\n")

pam50_res <- molecular.subtyping(
  sbt.model  = "pam50",
  data       = ddata,
  annot      = dannot,
  do.mapping = TRUE
)

cat("PAM50 subtype counts (recomputed):\n")
print(table(pam50_res$subtype))

pam50_df <- data.frame(
  sample_id        = rownames(ddata),
  PAM50_recomputed = as.character(pam50_res$subtype),
  stringsAsFactors = FALSE
)

## 5. Load Xena clinical matrix ----
cat("Reading Xena clinical matrix from:", clinical_file, "\n")
clin_df <- read.delim(clinical_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"sample" %in% colnames(clin_df)) {
  stop("Clinical file does not contain a 'sample' column. Check its column names.")
}
if (!"PAM50Call_RNAseq" %in% colnames(clin_df)) {
  warning("Clinical file does not contain 'PAM50Call_RNAseq'. PAM50Call_RNAseq will be NA.")
}

clin_pam <- clin_df[, c("sample", intersect("PAM50Call_RNAseq", colnames(clin_df)))]

## 6. Harmonize sample IDs (cut to first 16 characters) ----
short_id <- function(x) substr(x, 1, 16)

pam50_df$sample_short <- short_id(pam50_df$sample_id)
clin_pam$sample_short <- short_id(clin_pam$sample)

## 7. Merge recomputed PAM50 with Xena PAM50Call_RNAseq ----
merged_pam <- merge(
  pam50_df,
  clin_pam,
  by = "sample_short",
  all.x = TRUE,
  suffixes = c("", "_clin")
)

# Prefer Xena PAM50Call_RNAseq when available, otherwise use recomputed
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

## 8. Save results ----
cat("Writing merged PAM50 table to:", output_file, "\n")
write.table(
  merged_pam,
  file      = output_file,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

cat("Done.\n")
