#!/usr/bin/env Rscript
# estimate_purity_from_rnaseq.R (robust to duplicate column names)

suppressPackageStartupMessages({
  if (!requireNamespace("tidyestimate", quietly = TRUE)) {
    install.packages("tidyestimate", repos = "https://cloud.r-project.org")
  }
  library(tidyestimate)
})

# =========================
# USER SETTINGS (EDIT ME)
# =========================
in_path  <- "/home/stavz/masters/gdc/RNAexp/BRCA_log2_expression"        # <-- your input TSV
out_path <- "estimate_scores.tsv"    # <-- your output TSV

gene_col_name <- NULL
# If your gene symbols are in the first column, leave as NULL.
# If not, set to the exact column name as it appears in the file header.

input_is_logged <- TRUE
# TRUE if already log2(TPM+1) / log2(FPKM+1)
# FALSE if TPM/FPKM-like numeric values

# =========================
# Helpers
# =========================
collapse_duplicates_mean <- function(df, gene_col) {
  gene <- trimws(as.character(df[[gene_col]]))
  df[[gene_col]] <- gene
  df <- df[!is.na(df[[gene_col]]) & df[[gene_col]] != "", , drop = FALSE]

  sample_cols <- setdiff(colnames(df), gene_col)

  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  keep <- rowSums(!is.na(mat)) > 0
  mat <- mat[keep, , drop = FALSE]
  gene <- df[[gene_col]][keep]

  agg <- aggregate(mat, by = list(gene = gene), FUN = function(x) mean(x, na.rm = TRUE))
  colnames(agg)[1] <- gene_col
  agg
}

# =========================
# Read input (and sanitize colnames)
# =========================
expr <- read.table(
  in_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Make column names unique (THIS fixes your error)
colnames(expr) <- make.unique(colnames(expr), sep = ".")

# Pick gene column
gene_col <- if (!is.null(gene_col_name)) gene_col_name else colnames(expr)[1]
if (!(gene_col %in% colnames(expr))) {
  stop("Gene column '", gene_col, "' not found after sanitizing names.\n",
       "Available columns: ", paste(colnames(expr), collapse = ", "), "\n",
       "Tip: set gene_col_name to the correct header name.")
}

# Identify sample columns (everything except gene column)
sample_cols <- setdiff(colnames(expr), gene_col)
if (length(sample_cols) < 2) {
  stop("Not enough sample columns detected. Check that your TSV has samples as columns.")
}

# Convert expression columns to numeric
for (cn in sample_cols) {
  expr[[cn]] <- suppressWarnings(as.numeric(expr[[cn]]))
}

# Optional log2 transform
if (!input_is_logged) {
  expr[, sample_cols] <- log2(expr[, sample_cols] + 1)
}

# Collapse duplicate genes (mean)
expr_f <- collapse_duplicates_mean(expr, gene_col)

# # =========================
# # Filter to ESTIMATE genes
# # =========================
# expr_f <- tidyestimate::filter_common_genes(
#   expr,
#   id = gene_col,
#   tidy = FALSE,
#   tell_missing = TRUE,
#   find_alias = TRUE
# )

# =========================
# Compute scores (RNA-seq)
# =========================
scores <- tidyestimate::estimate_score(expr_f, is_affymetrix = FALSE)
names(scores) <- gsub("\\.", "_", names(scores))

# Optional: purity heuristic (RNA-seq = relative proxy)
estimate_col <- NULL
for (cand in c("estimate_score", "ESTIMATEScore", "estimatescore", "estimateScore")) {
  if (cand %in% names(scores)) { estimate_col <- cand; break }
}
if (!is.null(estimate_col)) {
  scores$purity_heuristic <- cos(0.6049872018 + 0.0001467884 * scores[[estimate_col]])
}

# Write output
write.table(scores, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\nInput:  ", in_path, "\nOutput: ", out_path, "\n", sep = "")
cat("Gene column used: ", gene_col, "\n", sep = "")
