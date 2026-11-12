library(tidyverse)
library(GenomicRanges)
script_output_dir = "./output/reciprocal_mapping/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
# forward (hg38 to mm10) and reverse (mm10 to hg38) mapped DHSs
regions_fwd = readRDS("./output/indexed_dhs_mapped_regions_v1.rds")
regions_rev = readRDS("./output/indexed_mouse_dhs_mapped_regions.rds")
parse_region <- function(region_str) {
        parts = strsplit(region_str, "[:-]")
        GRanges(seqnames = sapply(parts, `[`, 1),
                ranges = IRanges(start = as.numeric(sapply(parts, `[`, 2)),
                                 end = as.numeric(sapply(parts, `[`, 3))))
}
fwd_human_gr = parse_region(regions_fwd$hg38_region)
fwd_mouse_gr = parse_region(regions_fwd$mm10_region)
rev_mouse_gr = parse_region(regions_rev$mm10_region)
rev_human_gr = parse_region(regions_rev$hg38_region)

# match forward and reverse mappings on the mouse side, reciprocal if the human sides also overlap
mouse_overlaps = findOverlaps(fwd_mouse_gr, rev_mouse_gr, minoverlap = 50)
human_overlap_widths = width(pintersect(fwd_human_gr[queryHits(mouse_overlaps)],
                                        rev_human_gr[subjectHits(mouse_overlaps)]))
is_pair_reciprocal = human_overlap_widths > 50
fwd_tb = tibble(has_mouse_match = seq_along(fwd_human_gr) %in% queryHits(mouse_overlaps),
                is_reciprocal = seq_along(fwd_human_gr) %in% queryHits(mouse_overlaps)[is_pair_reciprocal]) %>%
        mutate(category = case_when(is_reciprocal ~ "Reciprocal",
                                    !has_mouse_match ~ "Forward-only",
                                    has_mouse_match & !is_reciprocal ~ "Non-reciprocal"))
rev_tb = tibble(has_mouse_match = seq_along(rev_mouse_gr) %in% subjectHits(mouse_overlaps),
                is_reciprocal = seq_along(rev_mouse_gr) %in% subjectHits(mouse_overlaps)[is_pair_reciprocal]) %>%
        mutate(category = case_when(is_reciprocal ~ "Reciprocal",
                                    !has_mouse_match ~ "Reverse-only",
                                    has_mouse_match & !is_reciprocal ~ "Non-reciprocal"))
# four categories: total REPs = forward mappings + reverse-only mappings
category_tb = tibble(category = c("Reciprocal", "Forward-only", "Reverse-only", "Non-reciprocal"),
                     n = c(sum(fwd_tb$category == "Reciprocal"),
                           sum(fwd_tb$category == "Forward-only"),
                           sum(rev_tb$category == "Reverse-only"),
                           sum(fwd_tb$category == "Non-reciprocal")))
category_tb$pct = 100 * category_tb$n / (nrow(regions_fwd) + sum(rev_tb$category == "Reverse-only"))

# overlap of the human side with segmental duplications (UCSC hg38 genomicSuperDups)
segdup_raw = read_tsv("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
                      col_names = F, col_types = cols(.default = "c"))
segdups = GRanges(seqnames = segdup_raw$X2,
                  ranges = IRanges(start = as.numeric(segdup_raw$X3) + 1,
                                   end = as.numeric(segdup_raw$X4)))
fwd_tb$overlaps_segdup = overlapsAny(fwd_human_gr, segdups)
rev_tb$overlaps_segdup = overlapsAny(rev_human_gr, segdups)
segdup_tb = bind_rows(fwd_tb, filter(rev_tb, category == "Reverse-only")) %>%
        group_by(category) %>%
        summarise(n = n(),
                  n_segdup = sum(overlaps_segdup),
                  pct_segdup = 100 * mean(overlaps_segdup)) %>%
        arrange(desc(pct_segdup))
segdup_tb$enrichment_vs_reciprocal = segdup_tb$pct_segdup / segdup_tb$pct_segdup[segdup_tb$category == "Reciprocal"]
# baseline: the three other categories pooled
other_tb = filter(segdup_tb, category != "Non-reciprocal")
pct_segdup_other = 100 * sum(other_tb$n_segdup) / sum(other_tb$n)
# non-reciprocal vs reciprocal forward mappings
nonrecip = fwd_tb$overlaps_segdup[fwd_tb$category == "Non-reciprocal"]
recip = fwd_tb$overlaps_segdup[fwd_tb$category == "Reciprocal"]
fisher_res = fisher.test(matrix(c(sum(nonrecip), sum(!nonrecip), sum(recip), sum(!recip)), nrow = 2))
write_tsv(category_tb, paste0(script_output_dir, "four_category_summary.tsv"))
write_tsv(segdup_tb, paste0(script_output_dir, "four_category_segdup_summary.tsv"))
write_tsv(tibble(pct_segdup_other = pct_segdup_other,
                 odds_ratio = fisher_res$estimate,
                 conf_low = fisher_res$conf.int[1],
                 conf_high = fisher_res$conf.int[2],
                 p_value = fisher_res$p.value),
          paste0(script_output_dir, "segdup_tests.tsv"))
