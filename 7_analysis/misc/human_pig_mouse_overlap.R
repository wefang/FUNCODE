library(tidyverse)
script_output_dir = "./output/human_pig_conservation/"
scores = readRDS(paste0(script_output_dir, "conservation_scores.rds"))
regions = readRDS("./output/indexed_dhs_mapped_regions_v1.rds")
# human elements in both the human-pig and the human-mouse core sets
common = intersect(scores$identifier, regions$id)
scores = scores[match(common, scores$identifier), ]
regions = regions[match(common, regions$id), ]
# high in human-pig at FDR < 0.25, high in human-mouse at FDR < 0.1, NA excluded
overlap_test = function(hp_fdr, hm_fdr) {
        ok = !is.na(hp_fdr) & !is.na(hm_fdr)
        x = hp_fdr[ok] < 0.25
        y = hm_fdr[ok] < 0.1
        fisher_res = fisher.test(table(factor(x, c(F, T)), factor(y, c(F, T))))
        tibble(hp_thr = 0.25,
               hm_thr = 0.1,
               n_evaluable = sum(ok),
               n_hp_high = sum(x),
               n_hm_high = sum(y),
               n_both_high = sum(x & y),
               n_hp_only = sum(x & !y),
               n_hm_only = sum(!x & y),
               odds_ratio = unname(fisher_res$estimate),
               fisher_p = fisher_res$p.value)
}
overlap_tb = bind_rows(
        mutate(overlap_test(scores$cov_fdr, regions$cov_chromatin_accessibility_fdr), score_type = "CO-V", .before = 1),
        mutate(overlap_test(scores$cob_fdr, regions$cob_chromatin_accessibility_fdr), score_type = "CO-B", .before = 1))
write_tsv(overlap_tb, paste0(script_output_dir, "human_mouse_overlap.tsv"))
