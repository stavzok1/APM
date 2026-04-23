#!/usr/bin/env Rscript

library(edgeR)
library(preprocessCore)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: normalize_atac_edger.R <input> <output> <first_sample_col>")
}

input_file <- args[1]
output_file <- args[2]
first_sample_col <- as.integer(args[3])

# 1. Read the data
df <- read.delim(input_file, check.names = FALSE)  # keep column names as-is

meta <- df[, 1:first_sample_col]
counts <- df[, (first_sample_col+1):ncol(df)]

# 2. edgeR logCPM with prior.count=5
# counts is peaks x samples (as in your Python)
logCPM <- cpm(as.matrix(counts), log = TRUE, prior.count = 5)

# 3. Quantile normalization across samples
logCPM_qn <- normalize.quantiles(logCPM)
colnames(logCPM_qn) <- colnames(counts)
rownames(logCPM_qn) <- rownames(counts)

# 4. Combine meta + normalized matrix, write out
out <- cbind(meta, as.data.frame(logCPM_qn, check.names = FALSE))
write.table(out, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
