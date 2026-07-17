#!/usr/bin/env Rscript
# scanMiR biochemical K_D scan → per-(arm, gene) predicted repression for the occupancy κ (learned/kd.py).
# Inputs:  data/external_cache/scanmir/{hallmark_3utr.fa, model_arms.txt}   (written by kd.build)
# Output:  data/external_cache/scanmir/kd_repression.tsv.gz  (arm, gene, repression, site counts)
# scanMiR McGeary-Bartel KdModels (getKdModels("hsa")); aggregateMatches → per-transcript predicted logFC.
.libPaths(c(file.path(Sys.getenv("HOME"), "R_scanmir"), .libPaths()))
suppressMessages({library(scanMiR); library(scanMiRData); library(Biostrings)})
mods <- getKdModels("hsa")
arms <- intersect(readLines("data/external_cache/scanmir/model_arms.txt"), names(mods))
seqs <- readDNAStringSet("data/external_cache/scanmir/hallmark_3utr.fa")
cat("[scan] arms:", length(arms), " UTRs:", length(seqs), "\n"); flush.console()
chunks <- split(arms, ceiling(seq_along(arms) / 60))
res <- list(); t0 <- Sys.time()
for (i in seq_along(chunks)) {
  mm <- findSeedMatches(seqs, mods[chunks[[i]]], verbose = FALSE, ret = "GRanges")
  res[[i]] <- aggregateMatches(mm)
  cat("[scan] chunk", i, "/", length(chunks), " cum_rows=", sum(sapply(res, nrow)),
      " elapsed=", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n"); flush.console()
}
ag <- do.call(rbind, res)
ag <- ag[, c("miRNA", "transcript", "repression", "8mer", "7mer", "6mer", "non-canonical")]
colnames(ag) <- c("arm", "gene", "repression", "n_8mer", "n_7mer", "n_6mer", "n_noncanon")
gz <- gzfile("data/external_cache/scanmir/kd_repression.tsv.gz", "w")
write.table(ag, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
cat("[scan] DONE wrote", nrow(ag), "rows,", length(unique(ag$arm)), "arms x",
    length(unique(ag$gene)), "genes\n")
