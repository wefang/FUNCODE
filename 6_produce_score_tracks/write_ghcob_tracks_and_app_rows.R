# GH-CO-B for the track hub and the shiny app, built from the significant pairs (FDR < 0.25).
# Hub: 8 bigWigs, value = (max CO-B over an element's significant pairs)^5, as for the cob tracks.
# App: one row per significant pair and modality, in shinyapp/data/unalign_ghcob.rds. app.R appends
# these to the unalign table, so the existing data files stay as they are.
# GH-CO-V tracks cover every scored element. These cover significant elements only, because the
# unfiltered histone scores are not on this machine.
suppressPackageStartupMessages({library(rtracklayer); library(GenomicRanges); library(dplyr); library(readr)})
setwd("F:/0_Projects/funcode_revision")

files <- c("chromatin accessibility" = "output/MTD-3_extended_cob/cob_significant_fdr25.rds",
           "H3K27ac" = "output/MTD-3_extended_H3K27ac/significant_fdr25.rds",
           "H3K4me1" = "output/MTD-3_extended_H3K4me1/significant_fdr25.rds",
           "H3K4me3" = "output/MTD-3_extended_H3K4me3/significant_fdr25.rds")
sig <- bind_rows(lapply(names(files), function(m) {
  s <- as.data.frame(readRDS(files[[m]]))
  data.frame(human_dhs_id = as.character(s$human_id), mouse_dhs_id = as.character(s$mouse_id), modality = m,
             cob = if ("cob" %in% names(s)) s$cob else s$score, cob_fdr = s$fdr)
}))
cat("significant rows", nrow(sig), "\n")

# element coordinates
hd <- read_tsv("funcCons/metadata/ENCODE_DHS/ENCFF503GCK.tsv", col_types = cols(identifier = "c", .default = "?"), progress = FALSE)
hd <- data.frame(human_dhs_id = hd$identifier, human_chromosome = as.character(hd$seqname), human_start = as.integer(hd$start), human_end = as.integer(hd$end))
md <- read.delim("funcCons/metadata/ENCODE_DHS/ENCFF910SRW.tsv")
md <- data.frame(mouse_dhs_id = as.character(md$index_num), mouse_chromosome = as.character(md$chr), mouse_start = as.integer(md$start), mouse_end = as.integer(md$end))

# element genes (an element in a multi-gene homology set can have several)
collapse_sym <- function(tb, id, sym) {
  x <- as.data.frame(tb)[, c(id, sym)]; names(x) <- c("id", "sym"); x$id <- as.character(x$id)
  x <- distinct(x[x$id %in% unique(sig[[ifelse(id == "human_id", "human_dhs_id", "mouse_dhs_id")]]), ])
  x %>% group_by(id) %>% summarize(sym = paste(sym, collapse = ","), .groups = "drop")
}
hs <- collapse_sym(readRDS("funcCons/non_alignable_env/annot/human_dhs_info.rds"), "human_id", "human_symbol")
ms <- collapse_sym(readRDS("funcCons/non_alignable_env/annot/mouse_dhs_info.rds"), "mouse_id", "mouse_symbol")

# shared motif, per pair and modality
mf <- readRDS("funcCons/output/nonalign_ghcob_motif_shared.rds")$pair_modality_tb
mf$modality[mf$modality == "CA"] <- "chromatin accessibility"
mf$mouse_dhs_id <- sub("^MouseDHS-", "", mf$mouse_id)
mf <- transmute(mf, human_dhs_id = as.character(human_id), mouse_dhs_id, modality, shared_motif = n_motif_shared > 0)

u <- readRDS("funcode/shinyapp/data/unalign.rds")
gh <- sig %>% left_join(hd, by = "human_dhs_id") %>% left_join(md, by = "mouse_dhs_id") %>%
  left_join(rename(hs, human_dhs_id = id, human_symbol = sym), by = "human_dhs_id") %>%
  left_join(rename(ms, mouse_dhs_id = id, mouse_symbol = sym), by = "mouse_dhs_id") %>%
  left_join(mf, by = c("human_dhs_id", "mouse_dhs_id", "modality")) %>%
  mutate(cov = NA_real_, cov_fdr = NA_real_, index_id = max(u$index_id) + row_number())
stopifnot(nrow(gh) == nrow(sig), !anyNA(gh$human_start), !anyNA(gh$mouse_start), !anyNA(gh$shared_motif))
cat("missing symbols: human", sum(is.na(gh$human_symbol)), " mouse", sum(is.na(gh$mouse_symbol)), "\n")

# the coordinate tables must agree with the app's rows for the same elements
h1 <- distinct(as.data.frame(u)[, c("human_dhs_id", "human_chromosome", "human_start", "human_end")]) %>% inner_join(hd, by = "human_dhs_id", suffix = c("", ".n"))
m1 <- distinct(as.data.frame(u)[, c("mouse_dhs_id", "mouse_chromosome", "mouse_start", "mouse_end")]) %>% mutate(mouse_dhs_id = as.character(mouse_dhs_id)) %>% inner_join(md, by = "mouse_dhs_id", suffix = c("", ".n"))
cat("human coord agreement", nrow(h1), mean(h1$human_chromosome == h1$human_chromosome.n & h1$human_start == h1$human_start.n & h1$human_end == h1$human_end.n), "\n")
cat("mouse coord agreement", nrow(m1), mean(m1$mouse_chromosome == m1$mouse_chromosome.n & m1$mouse_start == m1$mouse_start.n & m1$mouse_end == m1$mouse_end.n), "\n")
key <- function(x) paste(x$human_dhs_id, x$mouse_dhs_id, x$modality)
cat("pairs already in the app with a GH-CO-V row", sum(key(gh) %in% key(as.data.frame(u)[, c("human_dhs_id", "mouse_dhs_id", "modality")])), "\n")

saveRDS(as.data.frame(gh), "funcode/shinyapp/data/unalign_ghcob.rds")

# bigWigs
assay_file <- c("chromatin accessibility" = "chromatin_accessibility", H3K27ac = "H3K27ac", H3K4me1 = "H3K4me1", H3K4me3 = "H3K4me3")
write_track <- function(tb, chr, start, end, id, genome, assay) {
  per <- tb %>% group_by(.data[[id]]) %>% summarize(chr = first(.data[[chr]]), start = first(.data[[start]]), end = first(.data[[end]]), score = max(cob)^5, .groups = "drop")
  ref <- BigWigFile(paste0("funcode/track_hubs/data/gh_cov_", assay, "_", genome, ".bw"))
  si <- seqinfo(ref)
  per <- per[per$chr %in% seqlevels(si), ]
  g <- GRanges(per$chr, IRanges(per$start, per$end), score = per$score, seqinfo = si)
  dj <- disjoin(g)
  ol <- findOverlaps(g, dj)
  score(dj) <- aggregate(score(g)[queryHits(ol)], list(subjectHits(ol)), max)$x
  f <- paste0("funcode/track_hubs/data/gh_cob_", assay, "_", genome, ".bw")
  export(dj, BigWigFile(f))
  back <- import(f)
  cat(basename(f), "elements", nrow(per), "intervals", length(back), "range", round(range(back$score), 3), "size_kb", round(file.size(f) / 1024), "\n")
}
for (m in names(files)) {
  x <- gh[gh$modality == m, ]
  write_track(x, "human_chromosome", "human_start", "human_end", "human_dhs_id", "hg38", assay_file[[m]])
  write_track(x, "mouse_chromosome", "mouse_start", "mouse_end", "mouse_dhs_id", "mm10", assay_file[[m]])
}
