library(tidyverse)
source("./helper/helper.R")
script_output_dir = "./output/novel_element/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
regions = readRDS("./output/indexed_dhs_mapped_regions_v1.rds")
regions$human_cre = ifelse(regions$human_cre_class_prio_v4 != "non-cCRE", "cCRE", "non-cCRE")
regions$mouse_cre = ifelse(regions$mouse_cre_class_prio_v4 != "non-cCRE", "cCRE", "non-cCRE")

# mean, max and standardized variance of signal over anchor samples, per mark and species
mark_tb = tibble(mark = c("CA", "H3K27ac", "H3K4me1", "H3K4me3"),
                 score_mark = c("chromatin_accessibility", "H3K27ac", "H3K4me1", "H3K4me3"),
                 human_file = c("./intermediate_data/human_chromatin_correct.rds",
                                paste0("./processed_data/histone_norm_correct/Homo_sapiens_", c("H3K27ac", "H3K4me1", "H3K4me3"), "_norm_correct.rds")),
                 mouse_file = c("./intermediate_data/mouse_chromatin_correct.rds",
                                paste0("./processed_data/histone_norm_correct/Mus_musculus_", c("H3K27ac", "H3K4me1", "H3K4me3"), "_norm_correct.rds")),
                 anchor_file = c("./output/chromatin_conservation/anchors.tsv",
                                 paste0("./output/histone_conservation/anchors_", c("H3K27ac", "H3K4me1", "H3K4me3"), "_v1.rds")))
for (k in 1:nrow(mark_tb)) {
        mark = mark_tb$mark[k]
        if (mark == "CA") {
                anchors = read_tsv(mark_tb$anchor_file[k])
        } else {
                anchors = readRDS(mark_tb$anchor_file[k])
        }
        anchors = anchors[anchors$score > 0.5, ]
        human_mat = readRDS(mark_tb$human_file[k])
        mouse_mat = readRDS(mark_tb$mouse_file[k])
        anchors = anchors[anchors$accession1 %in% colnames(human_mat) & anchors$accession2 %in% colnames(mouse_mat), ]
        human_mat = as.matrix(human_mat[, anchors$accession1, drop = F])
        mouse_mat = as.matrix(mouse_mat[, anchors$accession2, drop = F])
        human_log = log2(human_mat + 1)
        mouse_log = log2(mouse_mat + 1)
        regions[[paste0("human_", mark, "_mean")]] = rowMeans(human_log)
        regions[[paste0("human_", mark, "_max")]] = apply(human_log, 1, max)
        regions[[paste0("human_", mark, "_var_std")]] = calc_var_std_ver1(human_mat)
        regions[[paste0("mouse_", mark, "_mean")]] = rowMeans(mouse_log)
        regions[[paste0("mouse_", mark, "_max")]] = apply(mouse_log, 1, max)
        regions[[paste0("mouse_", mark, "_var_std")]] = calc_var_std_ver1(mouse_mat)
        rm(human_mat, mouse_mat, human_log, mouse_log)
        gc()
}

# categories by cCRE annotation and mark-specific conservation (CO-V or CO-B significant)
compare_categories = c("cCRE_conserved", "cCRE_nonconserved", "novel_conserved", "DHS_nonconserved")
category_list = map(c("human", "mouse"), function(species) {
        if (species == "human") {
                regions_data = regions
        } else {
                regions_data = filter(regions, !is.na(mouse_id))
        }
        cre = regions_data[[paste0(species, "_cre")]]
        map(1:nrow(mark_tb), function(k) {
                mean_col = paste0(species, "_", mark_tb$mark[k], "_mean")
                max_col = paste0(species, "_", mark_tb$mark[k], "_max")
                var_col = paste0(species, "_", mark_tb$mark[k], "_var_std")
                mark_conserved = ifelse(regions_data[[paste0("cov_", mark_tb$score_mark[k], "_sign")]] == 1 |
                                                regions_data[[paste0("cob_", mark_tb$score_mark[k], "_sign")]] == 1,
                                        "Conserved", "non-Conserved")
                regions_data$category = case_when(cre == "cCRE" & mark_conserved == "Conserved" ~ "cCRE_conserved",
                                                  cre == "cCRE" & mark_conserved != "Conserved" ~ "cCRE_nonconserved",
                                                  cre == "non-cCRE" & mark_conserved == "Conserved" ~ "novel_conserved",
                                                  cre == "non-cCRE" & mark_conserved != "Conserved" ~ "DHS_nonconserved",
                                                  TRUE ~ "other")
                # 5th percentile of mean signal among cCREs
                ccre_threshold = quantile(regions_data[[mean_col]][cre == "cCRE"], probs = 0.05, na.rm = T)
                summary_tb = regions_data %>%
                        filter(category %in% compare_categories) %>%
                        group_by(mark_novel_category = category) %>%
                        summarise(n = n(),
                                  mean_signal_mean = mean(.data[[mean_col]], na.rm = T),
                                  median_signal_mean = median(.data[[mean_col]], na.rm = T),
                                  mean_signal_max = mean(.data[[max_col]], na.rm = T),
                                  median_signal_max = median(.data[[max_col]], na.rm = T),
                                  mean_var_std = mean(.data[[var_col]], na.rm = T),
                                  median_var_std = median(.data[[var_col]], na.rm = T),
                                  mean_cov = mean(.data[[paste0("cov_", mark_tb$score_mark[k])]], na.rm = T),
                                  frac_above_threshold = mean(.data[[mean_col]] >= ccre_threshold, na.rm = T),
                                  .groups = "drop") %>%
                        mutate(mark = mark_tb$mark[k], species = species)
                wilcox_res = wilcox.test(regions_data[[mean_col]][regions_data$category == "novel_conserved"],
                                         regions_data[[mean_col]][regions_data$category == "cCRE_conserved"])
                list(summary_tb = summary_tb,
                     wilcox_tb = tibble(mark = mark_tb$mark[k], species = species,
                                        wilcox_pval_novel_vs_ccre_conserved = wilcox_res$p.value))
        })
}) %>% unlist(recursive = F)
write_tsv(bind_rows(map(category_list, "summary_tb")), paste0(script_output_dir, "category_summary_all_marks.tsv"))
write_tsv(bind_rows(map(category_list, "wilcox_tb")), paste0(script_output_dir, "category_wilcox_tests.tsv"))
# per-element mean signal for the density panels
saveRDS(regions[c("id", "mouse_id", "human_cre_class_prio_v4", "mouse_cre_class_prio_v4",
                  paste0(c("cov_", "cob_"), rep(mark_tb$score_mark, each = 2), "_sign"),
                  "human_cre", "mouse_cre",
                  paste0(c("human_", "mouse_"), rep(mark_tb$mark, each = 2), "_mean"))],
        paste0(script_output_dir, "signal_mean.rds"))
