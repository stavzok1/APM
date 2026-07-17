#!/usr/bin/env Rscript
# Genome-wide scanMiR K_D scan for the strong-site-vs-context++ specificity bench (per-arm, global denominator).
# Outputs BOTH all-site and strong-site (8mer+7mer) aggregated repression per (arm, gene).
.libPaths(c(file.path(Sys.getenv("HOME"), "R_scanmir"), .libPaths()))
suppressMessages({library(scanMiR); library(scanMiRData); library(Biostrings)})
mods <- getKdModels("hsa")
arms <- intersect(readLines("data/external_cache/scanmir/genomewide_arms.txt"), names(mods))
seqs <- readDNAStringSet("data/external_cache/scanmir/genomewide_3utr.fa")
cat("[gw] arms:", length(arms), " UTRs:", length(seqs), "\n"); flush.console()
chunks <- split(arms, ceiling(seq_along(arms) / 20))
resA <- list(); resS <- list(); t0 <- Sys.time()
for (i in seq_along(chunks)) {
  mm <- findSeedMatches(seqs, mods[chunks[[i]]], verbose = FALSE, ret = "GRanges")
  resA[[i]] <- aggregateMatches(mm)
  ty <- as.character(mm$type)
  mm_s <- mm[grepl("8mer|7mer", ty)]
  resS[[i]] <- if (length(mm_s)) aggregateMatches(mm_s) else resA[[i]][0,]
  cat("[gw] chunk", i, "/", length(chunks), " elapsed=",
      round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n"); flush.console()
}
agA <- do.call(rbind, resA); agS <- do.call(rbind, resS)
oA <- agA[, c("miRNA","transcript","repression")]; colnames(oA) <- c("arm","gene","repression_all")
oS <- agS[, c("miRNA","transcript","repression")]; colnames(oS) <- c("arm","gene","repression_strong")
m <- merge(oA, oS, by=c("arm","gene"), all=TRUE)
gz <- gzfile("data/external_cache/scanmir/genomewide_kd.tsv.gz","w")
write.table(m, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
cat("[gw] DONE", nrow(m), "rows,", length(unique(m$arm)),"arms x", length(unique(m$gene)),"genes\n")
