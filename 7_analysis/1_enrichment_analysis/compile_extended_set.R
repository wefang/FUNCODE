library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
source("./helper/helper.R")

sign_tb_list = c("./non_alignable_env/chromatin_total_ws_sign_fdr0.25.tsv",
                 "./non_alignable_env/H3K27ac_total_ws_sign_fdr0.25.tsv",
                 "./non_alignable_env/H3K4me3_total_ws_sign_fdr0.25.tsv",
                 "./non_alignable_env/H3K4me1_total_ws_sign_fdr0.25.tsv") %>%
        map(., readr::read_tsv, col_types = c("human_id" = "c", "mouse_id" = "c"))
sign_tb = bind_rows(list(mutate(sign_tb_list[[1]], modality = "chromatin"),
                         mutate(sign_tb_list[[2]], modality = "H3K27ac"),
                         mutate(sign_tb_list[[3]], modality = "H3K4me3"),
                         mutate(sign_tb_list[[4]], modality = "H3K4me1")))
table(is.na(sign_tb$weighted_spearman))

sign_tb$mouse_id = paste0("MouseDHS-", sign_tb$mouse_id)
all(sign_tb$mouse_id %in% rownames(mouse_motif_red))
sign_tb$n_motif_shared = compute_motif_shared(human_motif_red, mouse_motif_red, sign_tb$human_id, sign_tb$mouse_id)
table(sign_tb$n_motif_shared)
saveRDS(sign_tb, file = "./output/nonalign_sign_combined.rds")

# Figure 2d: plot fraction with motif analysis
source("./helper/def_color.R")
sign_tb = readRDS("./output/nonalign_sign_combined.rds")
sign_tb$motif_sign = sign_tb$n_motif_shared >= 1.
sign_tb$weighted_spearman_bin = cut(sign_tb$weighted_spearman, breaks = c(-Inf, seq(0, 0.75, by = 0.05), 1))
g_frac_motif = sign_tb %>% group_by(weighted_spearman_bin, modality) %>% filter(!is.na(weighted_spearman_bin)) %>%
        summarise(mean_motif = mean(motif_sign),
                  se_motif = sqrt(mean_motif * (1 - mean_motif) / n()),
                  n = n()) %>%
        ggplot(data = .) + geom_line(aes(x = as.numeric(weighted_spearman_bin), y = mean_motif, color = modality)) +
        geom_errorbar(aes(x = as.numeric(weighted_spearman_bin),
                          ymin = mean_motif - se_motif, ymax = mean_motif + se_motif,
                          color = modality)) +
        scale_color_manual(values = data_mod_col1) +
        scale_x_continuous(breaks = unique(as.numeric(sign_tb$weighted_spearman_bin)), labels = unique(sign_tb$weighted_spearman_bin)) +
        labs(x = "GH-CO-V", y = "Fraction with shared TFBS motif")+
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.7),
              legend.position = "right")
push_pdf(g_frac_motif, file = "frac_shared_motif", w = 4, h = 2.5)

# GH-CO-B: shared motifs of significant pairs
ghcob_files = c(CA = "./output/MTD-3_extended_cob/cob_significant_fdr25.rds",
                H3K27ac = "./output/MTD-3_extended_H3K27ac/significant_fdr25.rds",
                H3K4me1 = "./output/MTD-3_extended_H3K4me1/significant_fdr25.rds",
                H3K4me3 = "./output/MTD-3_extended_H3K4me3/significant_fdr25.rds")
ghcob_tb = map2_dfr(ghcob_files, names(ghcob_files), function(f, m) {
        s = readRDS(f)
        tibble(human_id = as.character(s$human_id), mouse_id = paste0("MouseDHS-", s$mouse_id), modality = m)
})
ghcob_pair_tb = distinct(ghcob_tb, human_id, mouse_id)
ghcob_pair_tb$n_motif_shared = compute_motif_shared(human_motif_red, mouse_motif_red,
                                                    ghcob_pair_tb$human_id, ghcob_pair_tb$mouse_id)
human_motif_ids = sort(unique(ghcob_pair_tb$human_id[ghcob_pair_tb$n_motif_shared >= 1]))
mouse_motif_ids = sort(unique(ghcob_pair_tb$mouse_id[ghcob_pair_tb$n_motif_shared >= 1]))
saveRDS(list(pair_tb = ghcob_pair_tb,
             pair_modality_tb = left_join(ghcob_tb, ghcob_pair_tb, by = c("human_id", "mouse_id")),
             human_motif_ids = human_motif_ids,
             human_no_motif_ids = setdiff(unique(ghcob_pair_tb$human_id), human_motif_ids),
             mouse_motif_ids = mouse_motif_ids,
             mouse_no_motif_ids = setdiff(unique(ghcob_pair_tb$mouse_id), mouse_motif_ids)),
        file = "./output/nonalign_ghcob_motif_shared.rds")

# extended set sizes: GH-CO-V, GH-CO-B, their overlap and union, GH-CO-B pairs sharing motifs
ghcov_pair_tb = distinct(tibble(human_id = sign_tb$human_id, mouse_id = sign_tb$mouse_id))
n_overlap = nrow(inner_join(ghcov_pair_tb, ghcob_pair_tb, by = c("human_id", "mouse_id")))
c(gh_cov = nrow(ghcov_pair_tb), gh_cob = nrow(ghcob_pair_tb), overlap = n_overlap,
  union = nrow(ghcov_pair_tb) + nrow(ghcob_pair_tb) - n_overlap,
  gh_cob_motif = sum(ghcob_pair_tb$n_motif_shared >= 1))

# Figure S7e: fraction with shared motif by GH-CO-B bin, pairs with p < 0.25
pval_files = c(CA = "./output/MTD-3_extended_cob/cob_pval_filtered.rds",
               H3K27ac = "./output/MTD-3_extended_H3K27ac/pval_filtered.rds",
               H3K4me1 = "./output/MTD-3_extended_H3K4me1/pval_filtered.rds",
               H3K4me3 = "./output/MTD-3_extended_H3K4me3/pval_filtered.rds")
ghcob_frac_motif = map2_dfr(pval_files, names(pval_files), function(f, m) {
        s = readRDS(f)
        tibble(cob_bin = cut(s$cob, breaks = c(-Inf, seq(0.5, 1.0, by = 0.1), Inf), include.lowest = T, dig.lab = 3),
               motif_sign = compute_motif_shared(human_motif_red, mouse_motif_red,
                                                 s$human_id, paste0("MouseDHS-", s$mouse_id)) >= 1) %>%
                group_by(cob_bin) %>%
                summarise(n_motif = sum(motif_sign), n = n()) %>%
                mutate(modality = m, mean_motif = n_motif / n, se_motif = sqrt(mean_motif * (1 - mean_motif) / n))
})
ghcob_frac_motif
