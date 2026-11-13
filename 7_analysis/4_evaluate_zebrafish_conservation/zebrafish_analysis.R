library(tidyverse)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
source("./helper/helper.R")
source("./helper/def_color.R")

script_output_dir = "./output/zebrafish_peak/"
script_plot_dir = "./plots_v1/zebrafish_peak/"
regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")

seq_score_chr = c("phastCons4way",
                  "phastCons7way",
                  "phastCons17way",
                  "phastCons20way",
                  "phastCons30way",
                  "phastCons100way",
                  "phyloP4way",
                  "phyloP7way",
                  "phyloP17way",
                  "phyloP20way",
                  "phyloP30way",
                  "phyloP100way",
                  "percentage",
                  "lecif")
grey_ramp <- colorRampPalette(c("#525252", "#d9d9d9"))
seq_score_col = grey_ramp(length(seq_score_chr))
names(seq_score_col) = seq_score_chr

all_score_chr = c("cob_chromatin_accessibility",
                  "cob_H3K4me3",
                  "cob_H3K27ac",
                  "cob_H3K4me1",
                  "cov_chromatin_accessibility",
                  "cov_H3K4me3",
                  "cov_H3K27ac",
                  "cov_H3K4me1",
                  "cov_manual_chromatin_accessiblity",
                  "cov_manual_H3K4me3",
                  "cov_manual_H3K27ac",
                  "cov_manual_H3K4me1",
                  seq_score_chr)

score_col = c("cov_chromatin_accessibility" = "#d73027",
              "cov_H3K27ac" = "#f46d43",
              "cov_H3K4me3" = "#fdae61",
              "cov_H3K4me1" = "#fee08b",
              # "cob_mean" = "#006837",
              "cob_chromatin_accessibility" = "#1a9850",
              "cob_H3K4me3" = "#66bd63",
              "cob_H3K27ac" = "#a6d96a",
              "cob_H3K4me1" = "#d9ef8b",
              "cov_manual_chromatin_accessiblity" = "#0c2c84",
              "cov_manual_H3K27ac" = "#225ea8",
              "cov_manual_H3K4me3" = "#1d91c0",
              "cov_manual_H3K4me1" = "#41b6c4",
              seq_score_col)
regions_subset = regions[c("id", all_score_chr)]

# add within-species baselines (human mean and variability) and AdaLiftOver scores
baseline_tb = readRDS("./output/chromatin_conservation_v1/CA_species_baselines.rds")
adaliftover_tb = readRDS("./output/adaliftover/adaliftover_scores.rds")
regions_combined = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
ada_match = match(regions_subset$id, regions_combined$human_id)
regions_subset$adaliftover_h2m = adaliftover_tb$AdaLiftOver_hg38_to_mm10[ada_match]
regions_subset$adaliftover_m2h = adaliftover_tb$AdaLiftOver_mm10_to_hg38[ada_match]
baseline_match = match(regions_subset$id, baseline_tb$id)
regions_subset$human_CA_mean = baseline_tb$human_CA_mean[baseline_match]
regions_subset$human_CA_var_std = baseline_tb$human_CA_var_std[baseline_match]
new_score_chr = c("adaliftover_h2m", "adaliftover_m2h",
                  "human_CA_mean", "human_CA_var_std")


# compile zebrafish peaks
# loop over all tissue collect and construct peak indicator matrix for one modality at a time
# peak_names = paste0(c("ATAC", "H3K27ac", "H3K4me3"), ".peakSignal")
peak_names = paste0(c("ATAC", "H3K27ac", "H3K4me3"), ".peak")
all_peak_mat = map(peak_names, function(assay_name) {
        peak_mat = do.call(cbind, map(human_zebrafish_tissue_names$...1, function(zebrafish_tissue_name) {
                zebrafish_data = readRDS(paste0("./processed_data/Chaoran_Zebrafish_YueLab/Peak.check/", zebrafish_tissue_name, " (1).rds"))
                # zebrafish_data = readRDS(paste0("./processed_data/Zebrafish_YueLab_signal/", zebrafish_tissue_name, ".rds"))
                zebrafish_data[[assay_name]]
        }))
        colnames(peak_mat) = human_zebrafish_tissue_names$...1
        peak_mat
})
names(all_peak_mat) = peak_names
human_zebrafish_tissue_names = readxl::read_excel("./custom_metadata/zebrafish_human_tissue_names.xlsx", col_names = F)

eval_all_mod = bind_rows(map(peak_names, function(cur_peak_name) {
        message(cur_peak_name)
        if (cur_peak_name == "ATAC.peak") {
                human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
                human_meta = xlsx::read.xlsx("table_s2_temp1.xlsx", sheetIndex = 1)
                human_meta = human_meta[human_meta$Organism == "Homo sapiens", ]
                
                tissue_ind_list = map(1:nrow(human_zebrafish_tissue_names), function(i) {
                        human_tissue_name = strsplit(human_zebrafish_tissue_names$...2[i], ";")[[1]]
                        human_meta_sel = human_meta[human_meta$Life.stage %in% c("adult") & human_meta$tissue %in% human_tissue_name & human_meta$type %in% c("tissue"), ]
                        # 
                        human_col_ind = which(colnames(human_correct) %in% human_meta_sel$Accession)
                        human_col_ind
                })
                tissue_common_ind = which(map_dbl(tissue_ind_list, length) > 0)
                print(tissue_common_ind)
                human_zebrafish_tissue_sel = human_zebrafish_tissue_names[tissue_common_ind, ]
                human_correct_tissue = do.call(cbind, map(tissue_ind_list[tissue_common_ind], function(human_col_ind) {
                        rowMeans(human_correct[, human_col_ind, drop = F])
                }))
                colnames(human_correct_tissue) = human_zebrafish_tissue_sel$...2
        }
        if (cur_peak_name == "H3K4me3.peak") {
                human_correct = readRDS(paste0("./processed_data/histone_norm_correct/Homo_sapiens_", "H3K4me3", "_norm_correct.rds"))
                rownames(human_correct) = regions$id
                human_meta = xlsx::read.xlsx("table_s2_temp1.xlsx", sheetIndex = 2)
                human_meta = human_meta[human_meta$Organism == "Homo sapiens", ]
                tissue_ind_list = map(1:nrow(human_zebrafish_tissue_names), function(i) {
                        human_tissue_name = strsplit(human_zebrafish_tissue_names$...2[i], ";")[[1]]
                        human_meta_sel = human_meta[human_meta$Life.stage %in% c("adult") & human_meta$tissue %in% human_tissue_name & human_meta$type %in% c("tissue"), ]
                        # 
                        human_col_ind = which(colnames(human_correct) %in% human_meta_sel$Accession)
                        human_col_ind
                })
                tissue_common_ind = which(map_dbl(tissue_ind_list, length) > 0)
                human_zebrafish_tissue_sel = human_zebrafish_tissue_names[tissue_common_ind, ]
                human_correct_tissue = do.call(cbind, map(tissue_ind_list[tissue_common_ind], function(human_col_ind) {
                        rowMeans(human_correct[, human_col_ind, drop = F])
                }))
                colnames(human_correct_tissue) = human_zebrafish_tissue_sel$...2
        }
        if (cur_peak_name == "H3K27ac.peak") {
                human_correct = readRDS(paste0("./processed_data/histone_norm_correct/Homo_sapiens_", "H3K27ac", "_norm_correct.rds"))
                rownames(human_correct) = regions$id
                human_meta = xlsx::read.xlsx("table_s2_temp1.xlsx", sheetIndex = 3)
                human_meta = human_meta[human_meta$Organism == "Homo sapiens", ]
                tissue_ind_list = map(1:nrow(human_zebrafish_tissue_names), function(i) {
                        human_tissue_name = strsplit(human_zebrafish_tissue_names$...2[i], ";")[[1]]
                        human_meta_sel = human_meta[human_meta$Life.stage %in% c("adult") & human_meta$tissue %in% human_tissue_name & human_meta$type %in% c("tissue"), ]
                        human_meta_sel = human_meta[human_meta$tissue %in% human_tissue_name & human_meta$type %in% c("tissue"), ]
                        # 
                        human_col_ind = which(colnames(human_correct) %in% human_meta_sel$Accession)
                        human_col_ind
                })
                tissue_common_ind = which(map_dbl(tissue_ind_list, length) > 0)
                human_zebrafish_tissue_sel = human_zebrafish_tissue_names[tissue_common_ind, ]
                human_correct_tissue = do.call(cbind, map(tissue_ind_list[tissue_common_ind], function(human_col_ind) {
                        rowMeans(human_correct[, human_col_ind, drop = F])
                }))
                colnames(human_correct_tissue) = human_zebrafish_tissue_sel$...2
        }
        
        collect_list = list()
        for (i in 1:nrow(human_zebrafish_tissue_sel)) {
                # change for each tissue:
                zebrafish_tissue_name = human_zebrafish_tissue_sel$...1[i]
                human_tissue_name = strsplit(human_zebrafish_tissue_sel$...2[i], ";")[[1]]
                message(zebrafish_tissue_name)
                
                zebrafish_data = readRDS(paste0("./processed_data/Zebrafish_YueLab_signal/", zebrafish_tissue_name, ".rds"))
                zebrafish_data = left_join(zebrafish_data, regions_subset, by = c("dhs.id" = "id"))
                
                # zebrafish_data$human_chromatin_var_std_global_anchors = regions$human_chromatin_var_std_global_anchors[match(zebrafish_data$dhs.id, regions$id)]
                # 
                # zebrafish_data$cov_chromatin_accessibility = regions$cov_chromatin_accessibility[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$cov_manual_chromatin_accessiblity = regions$cov_manual_chromatin_accessiblity[match(zebrafish_data$dhs.id, regions$id)]
                # 
                # zebrafish_data$cov_H3K4me3 = regions$cov_H3K4me3[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$cov_manual_H3K4me3 = regions$cov_manual_H3K4me3[match(zebrafish_data$dhs.id, regions$id)]
                # 
                # zebrafish_data$cov_H3K27ac = regions$cov_H3K27ac[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$cov_manual_H3K27ac = regions$cov_manual_H3K27ac[match(zebrafish_data$dhs.id, regions$id)]
                # 
                # zebrafish_data$cov_H3K4me1 = regions$cov_H3K4me1[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$cov_manual_H3K4me1 = regions$cov_manual_H3K4me1[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$cov_mean3 = (zebrafish_data$cov_chromatin_accessibility + zebrafish_data$cov_H3K4me3 + zebrafish_data$cov_H3K27ac)/3
                # 
                # zebrafish_data$phastCons4way = regions$phastCons4way[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$phastCons100way = regions$phastCons100way[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$phyloP4way = regions$phyloP4way[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$phyloP100way = regions$phyloP100way[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$percentage = regions$percentage[match(zebrafish_data$dhs.id, regions$id)]
                # zebrafish_data$lecif = regions$lecif[match(zebrafish_data$dhs.id, regions$id)]
                
                sim_df = as_tibble(expand.grid(
                        score = all_score_chr,
                        min_rank = c(1e3, 2.5e3, 5e3, 10e3, 1.5e4, 2e4, 2.5e4, 5e4)))
                sim_df = arrange(sim_df, score)
                
                # fold change of peak coverage in zebrafish
                peak_mat = all_peak_mat[[cur_peak_name]]
                col_ind = which(colnames(peak_mat) %in% zebrafish_tissue_name)
                assertthat::assert_that(length(col_ind) > 0)
                zebrafish_log2fc = rowMeans(log2(1 + peak_mat[, col_ind, drop = F])) - rowMeans(log2(1 + peak_mat[, -col_ind, drop = F]))
                
                # human_tissue_name = strsplit(human_zebrafish_tissue_names$...2[i], ";")[[1]]
                # human_meta_sel = human_meta[human_meta$`Life stage` %in% c("adult", "child") & human_meta$`Tissue/cell types` %in% human_tissue_name & human_meta$`Biosample type` %in% c("tissue", "primary cell"), ]
                # tissue_indices = which(colnames(human_correct) %in% human_meta_sel$Accession)
                # 
                # human_log2fc = rowMeans(log2(1+human_correct[match(zebrafish_data$dhs.id, rownames(human_correct)), tissue_indices, drop = F])) -
                #         rowMeans(log2(1+human_correct[match(zebrafish_data$dhs.id, rownames(human_correct)), -tissue_indices, drop = F]))
                
                assertthat::assert_that(colnames(human_correct_tissue)[i] == human_tissue_name)
                human_log2fc = log2(1+human_correct_tissue[match(zebrafish_data$dhs.id, rownames(human_correct)), i]) -
                        rowMeans(log2(1+human_correct_tissue[match(zebrafish_data$dhs.id, rownames(human_correct)), -i]))
                
                zebrafish_data$human_log2fc = human_log2fc
                zebrafish_data$zebrafish_log2fc = zebrafish_log2fc
                
                sim_df$value = map_dbl(1:nrow(sim_df), function(i_sim) {
                        score_chr = as.character(sim_df$score[i_sim])
                        zebrafish_fil = zebrafish_data %>%
                                arrange(desc(!! rlang::sym(score_chr))) %>%
                                filter(!is.na(cov_chromatin_accessibility)) %>%
                                filter(row_number() <= sim_df$min_rank[i_sim])
                        mean(map_dbl(1:100, function(i) {
                                sample_indices = sample(nrow(zebrafish_fil), size = 1e3, replace = F)
                                cor(zebrafish_fil$zebrafish_log2fc[sample_indices],
                                    zebrafish_fil$human_log2fc[sample_indices], method = "spearman")
                        }))
                })
                sim_df$tissue = zebrafish_tissue_name
                collect_list = append(collect_list, list(sim_df))
        }
        eval_all = bind_rows(collect_list)
        eval_all$modality = cur_peak_name
        eval_all
}))
saveRDS(eval_all_mod, file = "./output/zebrafish_cov_res1.rds")

# CO-V evaluation of within-species baselines and AdaLiftOver
human_data_tb = tibble(peak_name = peak_names,
                       file = c("./intermediate_data/encodev4_human_chromatin_correct.rds",
                                "./processed_data/histone_norm_correct/Homo_sapiens_H3K27ac_norm_correct.rds",
                                "./processed_data/histone_norm_correct/Homo_sapiens_H3K4me3_norm_correct.rds"),
                       sheet = c(1, 3, 2),
                       cov_score = c("cov_chromatin_accessibility", "cov_H3K27ac", "cov_H3K4me3"))
eval_new_mod = bind_rows(map(1:nrow(human_data_tb), function(k) {
        cur_peak_name = human_data_tb$peak_name[k]
        cov_score = human_data_tb$cov_score[k]
        human_correct = readRDS(human_data_tb$file[k])
        if (cur_peak_name != "ATAC.peak") rownames(human_correct) = regions$id
        human_meta = readxl::read_excel("table_s2_temp1.xlsx", sheet = human_data_tb$sheet[k])
        human_meta = human_meta[human_meta$Organism == "Homo sapiens", ]
        tissue_ind_list = map(1:nrow(human_zebrafish_tissue_names), function(i) {
                human_tissue_name = strsplit(human_zebrafish_tissue_names$...2[i], ";")[[1]]
                if (cur_peak_name == "ATAC.peak") {
                        human_meta_sel = human_meta[human_meta$Life.stage %in% c("adult") & human_meta$Tissue.cell.type %in% human_tissue_name & human_meta$Biosample.type %in% c("tissue"), ]
                } else {
                        human_meta_sel = human_meta[human_meta$Tissue.cell.type %in% human_tissue_name & human_meta$Biosample.type %in% c("tissue"), ]
                }
                which(colnames(human_correct) %in% human_meta_sel$Accession)
        })
        tissue_common_ind = which(map_dbl(tissue_ind_list, length) > 0)
        human_zebrafish_tissue_sel = human_zebrafish_tissue_names[tissue_common_ind, ]
        human_correct_tissue = do.call(cbind, map(tissue_ind_list[tissue_common_ind], function(human_col_ind) {
                rowMeans(human_correct[, human_col_ind, drop = F])
        }))
        colnames(human_correct_tissue) = human_zebrafish_tissue_sel$...2
        peak_mat = all_peak_mat[[cur_peak_name]]
        bind_rows(map(1:nrow(human_zebrafish_tissue_sel), function(i) {
                zebrafish_tissue_name = human_zebrafish_tissue_sel$...1[i]
                zebrafish_file = paste0("./processed_data/Zebrafish_YueLab_signal/", zebrafish_tissue_name, ".rds")
                if (!file.exists(zebrafish_file)) return(NULL)
                zebrafish_data = left_join(readRDS(zebrafish_file), regions_subset, by = c("dhs.id" = "id"))
                col_ind = which(colnames(peak_mat) == zebrafish_tissue_name)
                zebrafish_data$zebrafish_log2fc = rowMeans(log2(1 + peak_mat[, col_ind, drop = F])) - rowMeans(log2(1 + peak_mat[, -col_ind, drop = F]))
                row_ind = match(zebrafish_data$dhs.id, rownames(human_correct))
                zebrafish_data$human_log2fc = log2(1 + human_correct_tissue[row_ind, i]) -
                        rowMeans(log2(1 + human_correct_tissue[row_ind, -i, drop = F]))
                sim_df = as_tibble(expand.grid(score = new_score_chr,
                                               min_rank = c(5e3, 10e3, 1.5e4, 2e4, 2.5e4, 5e4)))
                sim_df = arrange(sim_df, score)
                sim_df = sim_df[map_dbl(as.character(sim_df$score), function(s) sum(!is.na(zebrafish_data[[s]]))) >= 1000, ]
                sim_df$value = map_dbl(seq_len(nrow(sim_df)), function(i_sim) {
                        score_chr = as.character(sim_df$score[i_sim])
                        zebrafish_fil = zebrafish_data %>%
                                arrange(desc(!! rlang::sym(score_chr))) %>%
                                filter(!is.na(!! rlang::sym(cov_score))) %>%
                                filter(row_number() <= sim_df$min_rank[i_sim])
                        if (nrow(zebrafish_fil) < 1e3) return(NA)
                        mean(map_dbl(1:100, function(i) {
                                sample_indices = sample(nrow(zebrafish_fil), size = 1e3, replace = F)
                                cor(zebrafish_fil$zebrafish_log2fc[sample_indices],
                                    zebrafish_fil$human_log2fc[sample_indices], method = "spearman", use = "complete.obs")
                        }), na.rm = T)
                })
                sim_df$tissue = zebrafish_tissue_name
                filter(sim_df, !is.na(value))
        })) %>%
                mutate(modality = cur_peak_name)
}))
saveRDS(eval_new_mod, file = "./output/zebrafish_cov_res_new_scores.rds")
# eval_all_mod = readRDS("./output/zebrafish_cov_res.rds")

filter(eval_all_mod, min_rank == 1000,
       modality == "H3K27ac.peak")  %>% ggboxplot(x = "tissue", y = "value")

# begin complete version
g1 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "ATAC.peak") %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("Chromatin Accessibility") +
        theme(legend.position = "right")

g2 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "H3K4me3.peak") %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("H3K4me3") +
        theme(legend.position = "right")

g3 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "H3K27ac.peak") %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("H3K27ac") +
        theme(legend.position = "right")
(g1 + g2 + g3 + plot_layout(guides = "collect")) %>%
        push_pdf(file_name = "zebrafish_cov_eval_complete_new", w = 12.5, h = 3.)
# end complete version

# main figure version
score_col = c(# "cov_mean3" = "#377eb8",
        "cov_chromatin_accessibility" = "#D63127",
        "cov_manual_chromatin_accessiblity" = "#EE634C",
        "cov_H3K4me3" = "#D63127",
        "cov_manual_H3K4me3" = "#EE634C",
        "cov_H3K27ac" = "#D63127",
        "cov_manual_H3K27ac" = "#EE634C",
        "cov_H3K4me1" = "#D63127",
        "cov_manual_H3K4me1" = "#EE634C",
        "lecif" = "#D8D8D7",
        "phastCons4way" = "#BCBCBC",
        "phastCons100way" = "#BCBCBC",
        "phyloP4way" = "#969696",
        "phyloP100way" = "#969696",
        "percentage" = "#737374")
g1 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "ATAC.peak") %>%
        filter(score %in% c("cov_chromatin_accessibility",
                            "cov_manual_chromatin_accessiblity",
                            "lecif",
                            "phastCons100way",
                            "phyloP100way",
                            "percentage")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("Chromatin Accessibility") +
        theme(legend.position = "right")

g2 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "H3K4me3.peak") %>%
        filter(score %in% c("cov_H3K4me3",
                            "cov_manual_H3K4me3",
                            "lecif",
                            "phastCons100way",
                            "phyloP100way",
                            "percentage")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("H3K4me3") +
        theme(legend.position = "right")
(g1 + g2 + plot_layout(guides = "collect")) %>%
        push_pdf(file_name = "zebrafish_cov_eval_main_new", w = 7., h = 2.6)

# Figure 4c
score_col = c(# "cov_mean3" = "#377eb8",
        "cov_chromatin_accessibility" = "#D63127",
        "cov_manual_chromatin_accessiblity" = "#EE634C",
        "cov_H3K4me3" = "#D63127",
        "cov_manual_H3K4me3" = "#EE634C",
        "cov_H3K27ac" = "#D63127",
        "cov_manual_H3K27ac" = "#EE634C",
        "cov_H3K4me1" = "#D63127",
        "cov_manual_H3K4me1" = "#EE634C",
        "lecif" = "#D8D8D7",
        "phastCons4way" = "#BCBCBC",
        "phastCons100way" = "#BCBCBC",
        "phyloP4way" = "#969696",
        "phyloP100way" = "#969696",
        "percentage" = "#737374")

g1 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(modality == "ATAC.peak",
               score %in% c(#"cov_mean3",
                       "cov_chromatin_accessibility",
                       "cov_manual_chromatin_accessiblity",
                       "phastCons100way",
                       "phastCons100way",
                       "phyloP100way",
                       "percentage",
                       "lecif")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05) +
        scale_color_manual(values = score_col) +
        theme(legend.position = "bottom")

g2 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(modality == "H3K4me3.peak",
               score %in% c(#"cov_mean3",
                       "cov_H3K4me3",
                       "cov_manual_H3K4me3",
                       "phastCons100way",
                       "phyloP100way",
                       "percentage",
                       "lecif")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05) +
        scale_color_manual(values = score_col) +
        theme(legend.position = "bottom")

g3 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(modality == "H3K4me3.peak",
               score %in% c(#"cov_mean3",
                       "cov_H3K27ac",
                       "cov_manual_H3K27ac",
                       "phastCons100way",
                       "phyloP100way",
                       "percentage",
                       "lecif")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05) +
        # ggplot(aes(x = min_rank, y = value, color = score, group = score),
        #        data = .) + geom_line() + geom_point() +
        scale_color_manual(values = score_col) +
        theme(legend.position = "bottom")

(g1 + g2 + g3 + plot_layout(guides = "collect")) %>%
        push_pdf(file_name = "zebrafish_cov_eval", w = 12, h = 2.)
# end main figure version

# Figure 4d: CO-B evaluation
# zebrafish_data_all = bind_rows(map(1:nrow(human_zebrafish_tissue_names), function(i) {
#         # change for each tissue:
#         zebrafish_tissue_name = human_zebrafish_tissue_names$...1[i]
#         zebrafish_data = readRDS(paste0("./processed_data/Chaoran_Zebrafish_YueLab/Peak.check/", zebrafish_tissue_name, " (1).rds"))
#         zebrafish_data$tissue = zebrafish_tissue_name
#         zebrafish_data
# }))
i = 1
zebrafish_tissue_name = human_zebrafish_tissue_names$...1[i]
zebrafish_data = readRDS(paste0("./processed_data/Chaoran_Zebrafish_YueLab/Peak.check/", zebrafish_tissue_name, " (1).rds"))
zebrafish_data = left_join(zebrafish_data, regions_subset, by = c("dhs.id" = "id"))

eval_all_mod = bind_rows(map(peak_names, function(cur_peak_name) {
        message(cur_peak_name)
        sim_df = as_tibble(expand.grid(
                score = c(all_score_chr, new_score_chr),
                min_rank = c(1e3, 2.5e3, 5e3, 10e3, 1.5e4, 2e4, 2.5e4, 5e4)))
        sim_df = arrange(sim_df, score)

        peak_mat = all_peak_mat[[cur_peak_name]]
        zebrafish_data$zebrafish_peak_frac = rowMeans(peak_mat > 0)
        
        sim_df$value = map_dbl(1:nrow(sim_df), function(i_sim) {
                score_chr = as.character(sim_df$score[i_sim])
                zebrafish_fil = zebrafish_data %>%
                        arrange(desc(!! rlang::sym(score_chr))) %>%
                        filter(!is.na(cob_chromatin_accessibility)) %>%
                        filter(row_number() <= sim_df$min_rank[i_sim])
                mean(map_dbl(1:100, function(i) {
                        sample_indices = sample(nrow(zebrafish_fil), size = 1e3, replace = F)
                        mean(zebrafish_fil$zebrafish_peak_frac)
                }))
        })
        sim_df$modality = cur_peak_name
        sim_df
}))
score_col = c(
        # "cob_mean3" = "#bf812d",
        "cob_chromatin_accessibility" = "#72C064",
        "cob_H3K4me3" = "#72C064",
        "cob_H3K27ac" = "#72C064",
        "lecif" = "#D8D8D7",
        "phastCons4way" = "#BCBCBC",
        "phastCons100way" = "#BCBCBC",
        "phyloP4way" = "#969696",
        "phyloP100way" = "#969696",
        "percentage" = "#737374")

saveRDS(eval_all_mod, file = "./output/zebrafish_cob_res.rds")

eval_all_mod = readRDS("./output/zebrafish_cob_res.rds")

# Figure S8c: complete version with more scores
g1 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "ATAC.peak") %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("Chromatin Accessibility") +
        theme(legend.position = "right")

g2 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "H3K4me3.peak") %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("H3K4me3") +
        theme(legend.position = "right")
g3 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(min_rank > 1000) %>%
        filter(modality == "H3K27ac.peak") %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05,
               position = position_dodge()) +
        scale_x_discrete(breaks = sort(unique(eval_all_mod$min_rank)),
                         labels = paste0(sort(unique(eval_all_mod$min_rank))/1e3, "k")) +
        scale_color_manual(values = score_col) +
        ggtitle("H3K27ac") +
        theme(legend.position = "right")
(g1 + g2 + g3 + plot_layout(guides = "collect")) %>%
        push_pdf(file_name = "zebrafish_cob_eval_complete_v1", w = 12.5, h = 3.)
# end complete version

# Figure S8d: complete version with more scores
g1 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(modality == "ATAC.peak",
               score %in% c(#"cob_mean3",
                       "cob_chromatin_accessibility",
                       "phastCons100way",
                       "phastCons100way",
                       "phyloP100way",
                       "percentage",
                       "lecif")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05) +
        scale_color_manual(values = score_col)
g2 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(modality == "H3K4me3.peak",
               score %in% c(#"cob_mean3",
                       "cob_H3K4me3",
                       "phastCons100way",
                       "phastCons100way",
                       "phyloP100way",
                       "percentage",
                       "lecif")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05) +
        scale_color_manual(values = score_col)
g3 = eval_all_mod %>% group_by(min_rank, score, modality) %>%
        filter(modality == "H3K27ac.peak",
               score %in% c(#"cob_mean3",
                       "cob_H3K27ac",
                       "phastCons100way",
                       "phastCons100way",
                       "phyloP100way",
                       "percentage",
                       "lecif")) %>%
        summarise(value = mean(value)) %>%
        mutate(score = factor(score, levels = names(score_col))) %>%
        arrange(desc(score)) %>%
        ggline(x = "min_rank", y = "value", color = "score", point.size = 0.05) +
        scale_color_manual(values = score_col)

(g1 + g2 + g3 + plot_layout(guides = "collect")) # %>%
# push_pdf(file_name = "zebrafish_eval_cob1", w = 13, h = 2.)
