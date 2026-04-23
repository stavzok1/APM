#############################
## Recompute PAM50 for TCGA-BRCA and merge with Xena
## Save as: recompute_pam50_tcga_brca.R
#############################

## 0. Packages ----

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("genefu", "dplyr", "readr", "tibble", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "genefu") {
      BiocManager::install("genefu")
    } else {
      install.packages(pkg)
    }
  }
}

library(genefu)
library(dplyr)
library(readr)
library(tibble)
library(stringr)

## 1. File paths (edit if needed) ----

expr_file   <- "RNAexp/BRCA_log2_expression"
clinical_file <- "RNAexp/BRCA_clinical"
output_file  <- "RNAexp/BRCA_PAM50_merged"

## 2. Load expression matrix ----
## Assumes:
##   - first column = gene_symbol
##   - columns 2..n = samples (TCGA barcodes)
##   - log2-normalized values (e.g., log2(RSEM+1))

cat("Reading expression matrix from:", expr_file, "\n")
expr_df <- read_tsv(expr_file, col_types = cols())

# Check first few columns
cat("Expression columns:\n")
print(head(colnames(expr_df)))

# Set gene symbols as rownames
if (!"gene_symbol" %in% colnames(expr_df)) {
  stop("Could not find 'gene_symbol' column in expression file. Edit the script to match your column names.")
}

expr_mat <- expr_df %>%
  rename(gene_symbol = gene_symbol) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

cat("Expression matrix dimensions (genes x samples):", 
    paste(dim(expr_mat), collapse = " x "), "\n")

## 3. Load PAM50 data from genefu ----

data(pam50)  # loads pam50 object

# Inspect mapping
# str(pam50$centroids.map)

pam50_symbols <- unique(pam50$centroids.map$gene.symbol)
cat("Number of PAM50 genes in genefu object:", length(pam50_symbols), "\n")

## 4. Keep only genes present in both expr and PAM50 ----

common_genes <- intersect(rownames(expr_mat), pam50_symbols)
cat("Number of PAM50 genes found in your matrix:", length(common_genes), "\n")

if (length(common_genes) < 40) {
  warning("Fewer than 40 PAM50 genes found. Check that your rownames are HGNC symbols.")
}

expr_pam50 <- expr_mat[common_genes, , drop = FALSE]

## 5. Median-center genes (Perou-style) ----
## Center each gene across all samples: x_ij <- x_ij - median_i

cat("Median-centering PAM50 genes across samples...\n")
gene_medians <- apply(expr_pam50, 1, median, na.rm = TRUE)
expr_centered <- sweep(expr_pam50, 1, gene_medians, FUN = "-")

## 6. Run PAM50 classification with genefu::molecular.subtyping ----

cat("Running PAM50 molecular subtyping with genefu...\n")

res_pam50 <- molecular.subtyping(
  sdata      = NULL,
  geneid     = rownames(expr_centered),
  data       = expr_centered,
  annot      = pam50$centroids.map,
  do.mapping = TRUE,      # map your gene ids to PAM50
  std        = "none",    # we already centered; no additional std
  model      = "pam50",
  verbose    = TRUE
)

# Extract subtypes and correlations
pam50_subtypes <- res_pam50$subtype
pam50_cor      <- res_pam50$cor   # correlations to centroids

cat("PAM50 subtype counts (recomputed):\n")
print(table(pam50_subtypes))

## 7. Build PAM50 results data frame ----

pam50_df <- data.frame(
  sample_id        = colnames(expr_centered),
  PAM50_recomputed = as.character(pam50_subtypes),
  stringsAsFactors = FALSE
)

# Add centroid correlations if you like
pam50_cor_df <- as.data.frame(t(pam50_cor))
pam50_cor_df <- pam50_cor_df %>%
  mutate(sample_id = rownames(pam50_cor_df))

pam50_df <- pam50_df %>%
  left_join(pam50_cor_df, by = "sample_id")

## 8. Load Xena clinical matrix with PAM50Call_RNAseq ----

cat("Reading Xena clinical matrix from:", clinical_file, "\n")
clin_df <- read_tsv(clinical_file, col_types = cols())

# Check columns
cat("Clinical columns:\n")
print(head(colnames(clin_df)))

if (!"sample" %in% colnames(clin_df)) {
  stop("Clinical file does not contain 'sample' column. Edit the script to match your clinical file format.")
}
if (!"PAM50Call_RNAseq" %in% colnames(clin_df)) {
  warning("Clinical file does not contain 'PAM50Call_RNAseq'. You may be using a different Xena file.")
}

clin_pam <- clin_df %>%
  select(sample, PAM50Call_RNAseq = PAM50Call_RNAseq)

## 9. Harmonize sample IDs (expression vs clinical) ----
## Often:
##   - Expression samples look like: TCGA-XX-XXXX-01A-11R-2001-01
##   - Clinical 'sample' is:        TCGA-XX-XXXX-01
## Below we cut both to the first 16 characters (patient+sample),
## adjust here if your IDs look different.

cat("Harmonizing sample IDs (using first 16 characters)...\n")

pam50_df <- pam50_df %>%
  mutate(
    sample_short = str_sub(sample_id, 1, 16)
  )

clin_pam <- clin_pam %>%
  mutate(
    sample_short = str_sub(sample, 1, 16)
  )

## 10. Merge recomputed PAM50 with Xena PAM50Call_RNAseq ----

merged_pam <- pam50_df %>%
  left_join(clin_pam, by = "sample_short") %>%
  # Keep some clean columns at front
  select(
    sample_id,
    sample_short,
    PAM50_recomputed,
    PAM50Call_RNAseq,
    everything()
  )

cat("Merged table dimensions:", paste(dim(merged_pam), collapse = " x "), "\n")

# How many samples had Xena PAM50?
cat("Samples with Xena PAM50Call_RNAseq (non-NA):\n")
print(sum(!is.na(merged_pam$PAM50Call_RNAseq)))

# Define final PAM50 label: prefer Xena when available, otherwise recomputed
merged_pam <- merged_pam %>%
  mutate(
    PAM50_final = ifelse(
      !is.na(PAM50Call_RNAseq),
      PAM50Call_RNAseq,
      PAM50_recomputed
    )
  )

cat("Final PAM50 subtype counts (PAM50_final):\n")
print(table(merged_pam$PAM50_final))

## 11. Save result ----

cat("Writing merged PAM50 table to:", output_file, "\n")
write.table(
  merged_pam,
  file      = output_file,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

cat("Done.\n")
