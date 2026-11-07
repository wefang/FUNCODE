# GH-CO-B bigWigs for every element with a gene-homology partner, not only the significant ones.
#
# A pair's CO-B in the extended set is pmin(human hval, mouse hval) (02_compute_cob_significant_only.R),
# so an element's best CO-B over all its gene-homology partners is pmin(own hval, the largest partner
# hval) within each linked homology-set pair, maximised over those pairs. That needs the element hvals
# and the homology sets only, not the 1.3 billion pairs. Track value = best CO-B^5, as for the cob and
# the earlier gh_cob tracks.
#
# Check: the p-value is 1 - ecdf(CO-B) from one unstratified null, so significance (FDR < 0.25) is a
# CO-B cutoff; every element in the significant-only tracks must keep exactly its old value.
suppressPackageStartupMessages({ library(data.table); library(rtracklayer); library(GenomicRanges); library(readr) })
setwd("F:/0_Projects/funcode_revision")

load("funcCons/output/all_dhs_pairs.rda")          # hidset_pairs, human_hid_set, mouse_hid_set
long <- function(hs, col) {
  n <- vapply(hs[[col]], NROW, integer(1))
  data.table(set = rep(hs$hid_set, n), id = as.character(unlist(lapply(hs[[col]], function(x) x$id))))
}
hl <- unique(long(human_hid_set, "human_dhs")); ml <- unique(long(mouse_hid_set, "mouse_dhs"))
pairs <- unique(data.table(hs = hidset_pairs$human_hid_set, ms = hidset_pairs$mouse_hid_set))
cat("homology-set pairs", nrow(pairs), " human memberships", nrow(hl), " mouse memberships", nrow(ml), "\n")

hd <- read_tsv("funcCons/metadata/ENCODE_DHS/ENCFF503GCK.tsv", col_types = cols(identifier = "c", .default = "?"), progress = FALSE)
hd <- data.table(id = hd$identifier, chr = as.character(hd$seqname), start = as.integer(hd$start), end = as.integer(hd$end))
md <- fread("funcCons/metadata/ENCODE_DHS/ENCFF910SRW.tsv")
md <- data.table(id = as.character(md$index_num), chr = as.character(md$chr), start = as.integer(md$start), end = as.integer(md$end))

best_cob <- function(own, own_col, other, other_col, hv_own, hv_other) {
  o <- own[id %in% names(hv_own)][, v := hv_own[id]]
  x <- other[id %in% names(hv_other)][, v := hv_other[id]][, .(pmax_other = max(v)), by = set]
  p <- merge(pairs[, .(own_set = get(own_col), other_set = get(other_col))], x, by.x = "other_set", by.y = "set")
  p <- p[, .(pmax_other = max(pmax_other)), by = own_set]     # best partner over all linked sets
  r <- merge(o, p, by.x = "set", by.y = "own_set")
  r[, .(cob = max(pmin(v, pmax_other))), by = id]
}

dirs <- c(chromatin_accessibility = "output/MTD-3_extended_cob/", H3K27ac = "output/MTD-3_extended_H3K27ac/",
          H3K4me1 = "output/MTD-3_extended_H3K4me1/", H3K4me3 = "output/MTD-3_extended_H3K4me3/")
sigf <- c(chromatin_accessibility = "cob_significant_fdr25.rds", H3K27ac = "significant_fdr25.rds",
          H3K4me1 = "significant_fdr25.rds", H3K4me3 = "significant_fdr25.rds")
summ <- list()
for (assay in names(dirs)) {
  hv <- readRDS(paste0(dirs[[assay]], "human_element_hval.rds")); mv <- readRDS(paste0(dirs[[assay]], "mouse_element_hval.rds"))
  hb <- best_cob(hl, "hs", ml, "ms", hv, mv); mb <- best_cob(ml, "ms", hl, "hs", mv, hv)
  # check against the significant pairs
  s <- as.data.table(readRDS(paste0(dirs[[assay]], sigf[[assay]])))
  if (!"cob" %in% names(s)) s[, cob := score]
  cut <- min(s$cob)
  for (side in c("human", "mouse")) {
    b <- if (side == "human") hb else mb
    old <- s[, .(old = max(cob)), by = .(id = as.character(get(paste0(side, "_id"))))]
    j <- merge(old, b, by = "id", all.x = TRUE)
    stopifnot(!anyNA(j$cob), all(abs(j$cob - j$old) < 1e-9), sum(b$cob >= cut - 1e-9) == nrow(old))
  }
  cat(assay, ": human", nrow(hb), "elements, mouse", nrow(mb), "; significant elements unchanged (",
      length(unique(s$human_id)), "human,", length(unique(s$mouse_id)), "mouse, CO-B >=", cut, ")\n")
  for (g in c("hg38", "mm10")) {
    b <- merge(if (g == "hg38") hb else mb, if (g == "hg38") hd else md, by = "id")
    lost <- nrow(if (g == "hg38") hb else mb) - nrow(b)
    si <- seqinfo(BigWigFile(paste0("funcode/track_hubs/data/gh_cov_", assay, "_", g, ".bw")))
    b <- b[chr %in% seqlevels(si)]
    gr <- GRanges(b$chr, IRanges(b$start, b$end), score = b$cob^5, seqinfo = si)
    dj <- disjoin(gr); ol <- findOverlaps(gr, dj)
    score(dj) <- aggregate(score(gr)[queryHits(ol)], list(subjectHits(ol)), max)$x
    f <- paste0("funcode/track_hubs/data/gh_cob_", assay, "_", g, ".bw")
    export(sort(dj), BigWigFile(f))
    back <- import(f)
    summ[[paste(assay, g)]] <- data.frame(track = basename(f), elements = nrow(b), no_coordinates = lost,
      intervals = length(back), median = round(median(back$score), 4), share_ge_0.859 = round(mean(back$score >= 0.97^5 - 1e-9), 4),
      size_mb = round(file.size(f) / 1e6, 1))
  }
}
print(do.call(rbind, summ), row.names = FALSE)
