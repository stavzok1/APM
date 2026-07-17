#!/usr/bin/env Rscript
# Run the actual ESTIMATE algorithm (Yoshihara et al. 2013) on a GCT expression
# matrix: official 141-gene stromal + 141-gene immune signatures scored by the
# package's ssGSEA, producing StromalScore / ImmuneScore / ESTIMATEScore (and
# TumorPurity only for platform="affymetrix").
#
# Usage: Rscript run_estimate.R <input.gct> <output_scores.gct> [platform]
#   platform in {affymetrix, agilent, illumina}; default illumina (RNA-seq).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run_estimate.R <input.gct> <output_scores.gct> [platform]")
input_gct <- args[1]
output_scores <- args[2]
platform <- if (length(args) >= 3) args[3] else "illumina"

suppressMessages(library(estimate))
common <- tempfile(fileext = ".gct")
filterCommonGenes(input.f = input_gct, output.f = common, id = "GeneSymbol")
estimateScore(input.ds = common, output.ds = output_scores, platform = platform)
cat("ESTIMATE_DONE platform=", platform, "\n", sep = "")
