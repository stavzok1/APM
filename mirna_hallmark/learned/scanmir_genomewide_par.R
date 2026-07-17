#!/usr/bin/env Rscript
# Parallel genome-wide scanMiR K_D scan (BiocParallel MulticoreParam). Args: arms_file out_gz [ncores]
.libPaths(c(file.path(Sys.getenv("HOME"), "R_scanmir"), .libPaths()))
suppressMessages({library(scanMiR); library(scanMiRData); library(Biostrings); library(BiocParallel)})
a <- commandArgs(trailingOnly = TRUE)
arms_file <- a[1]; out_gz <- a[2]; nc <- if (length(a) >= 3) as.integer(a[3]) else 7L
BP <- MulticoreParam(nc, progressbar = FALSE)
mods <- getKdModels("hsa")
arms <- intersect(readLines(arms_file), names(mods))
seqs <- readDNAStringSet("data/external_cache/scanmir/genomewide_3utr.fa")
cat("[par] arms:", length(arms), " UTRs:", length(seqs), " cores:", nc, "\n"); flush.console()
chunks <- split(arms, ceiling(seq_along(arms) / 40)); resA <- list(); resS <- list(); t0 <- Sys.time()
for (i in seq_along(chunks)) {
  mm <- findSeedMatches(seqs, mods[chunks[[i]]], BP = BP, verbose = FALSE, ret = "GRanges")
  resA[[i]] <- aggregateMatches(mm)
  mm_s <- mm[grepl("8mer|7mer", as.character(mm$type))]
  resS[[i]] <- if (length(mm_s)) aggregateMatches(mm_s) else resA[[i]][0, ]
  cat("[par] chunk", i, "/", length(chunks), " elapsed=",
      round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n"); flush.console()
}
agA <- do.call(rbind, resA); agS <- do.call(rbind, resS)
oA <- agA[, c("miRNA", "transcript", "repression")]; colnames(oA) <- c("arm", "gene", "repression_all")
oS <- agS[, c("miRNA", "transcript", "repression")]; colnames(oS) <- c("arm", "gene", "repression_strong")
m <- merge(oA, oS, by = c("arm", "gene"), all = TRUE); gz <- gzfile(out_gz, "w")
write.table(m, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
cat("[par] DONE", nrow(m), "rows,", length(unique(m$arm)), "arms\n")
