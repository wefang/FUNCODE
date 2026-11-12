library(tidyverse)
library(GenomicRanges)
library(furrr)
source("./helper/helper.R")
script_output_dir = "./output/human_pig_cross_validation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
plan(multisession, workers = 6)
options(future.globals.maxSize = 16 * 1024^3)
regions = readRDS("./intermediate_data/human_pig_mapped_regions.rds")
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds")
regions = regions[regions$identifier %in% rownames(human_correct), ]
human_correct = human_correct[regions$identifier, ]
pig_correct = readRDS("./intermediate_data/pig_chromatin_correct.rds")[regions$identifier, ]
n_regions = nrow(regions)
exp_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)

# load anchors with score >= 0.7
anchors_all = distinct(read_tsv("./output/human_pig_conservation/anchors.tsv"))
anchors_all = filter(anchors_all, score >= 0.7,
                     accession1 %in% colnames(human_correct),
                     accession2 %in% colnames(pig_correct))

# five tissues, pig samples from anchors, human samples by biosample name
pig_tissue_tb = anchors_all %>%
        distinct(accession = accession2, tissue = tolower(pig_tissue)) %>%
        mutate(tissue = case_when(tissue %in% c("fat", "adipose") ~ "adipose",
                                  tissue %in% c("heart") ~ "heart",
                                  tissue %in% c("liver") ~ "liver",
                                  tissue %in% c("muscle", "skeletal muscle") ~ "muscle",
                                  tissue %in% c("spleen") ~ "spleen"))
human_tissue_tb = exp_meta %>%
        mutate(b = tolower(`Biosample term name`),
               tissue = case_when(grepl("adipose|fat|adipocyte", b) ~ "adipose",
                                  grepl("heart|cardiac|atri|ventric|aort", b) ~ "heart",
                                  grepl("liver|hepatocyte|hep g2|hepg2", b) ~ "liver",
                                  grepl("muscle|myoblast|myocyte|myotube|psoas|tongue", b) ~ "muscle",
                                  grepl("spleen", b) ~ "spleen")) %>%
        filter(Accession %in% colnames(human_correct), !is.na(tissue)) %>%
        distinct(accession = Accession, tissue)
tissue_tb = distinct(bind_rows(mutate(filter(pig_tissue_tb, !is.na(tissue)), species = "pig"),
                               mutate(human_tissue_tb, species = "human")))
samples_in = function(species_name, tissue_name) {
        tissue_tb$accession[tissue_tb$species == species_name & tissue_tb$tissue == tissue_name]
}

# human-mouse scores and sequence baselines, joined by largest hg38 overlap
regions_hm = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds") %>%
        filter(!is.na(hg38_region)) %>%
        separate(hg38_region, into = c("b_chr", "b_coords"), sep = ":", remove = F) %>%
        separate(b_coords, into = c("b_start", "b_end"), sep = "-") %>%
        mutate(b_start = as.numeric(b_start), b_end = as.numeric(b_end)) %>%
        filter(!is.na(b_chr), !is.na(b_start), !is.na(b_end))
hm_gr = GRanges(regions_hm$b_chr, IRanges(regions_hm$b_start, regions_hm$b_end))
core_coords = str_split_fixed(regions$hg38_extend_region, "[:-]", 3)
core_gr = GRanges(core_coords[, 1], IRanges(as.integer(core_coords[, 2]), as.integer(core_coords[, 3])))
ol = as.data.frame(findOverlaps(core_gr, hm_gr))
ol$width = width(pintersect(core_gr[ol$queryHits], hm_gr[ol$subjectHits]))
ol = ol %>% arrange(queryHits, desc(width)) %>% distinct(queryHits, .keep_all = T)
baseline_cols = c("phyloP4way", "phyloP7way", "phyloP17way", "phyloP20way", "phyloP30way", "phyloP100way",
                  "phastCons4way", "phastCons7way", "phastCons17way", "phastCons20way", "phastCons30way", "phastCons100way",
                  "cov_chromatin_accessibility", "cob_chromatin_accessibility", "cov_manual_chromatin_accessiblity", "lecif")
baseline = tibble(id = regions$identifier)
for (cc in baseline_cols) {
        baseline[[cc]] = NA_real_
        baseline[[cc]][ol$queryHits] = regions_hm[[cc]][ol$subjectHits]
}
rm(regions_hm)
# human-pig sequence homology
psh = read_tsv("./intermediate_data/hg38_susScr11_percentage.txt", col_types = cols_only(identifier = "c", percentage = "d")) %>%
        group_by(id = identifier) %>%
        summarize(percentage = suppressWarnings(max(percentage, na.rm = T))) %>%
        mutate(percentage = ifelse(is.finite(percentage), percentage, NA))
baseline = left_join(baseline, psh, by = "id")
# AdaLiftOver, both directions
load("./metadata_processed/indexed_aligned_combined.rda")
ada = readRDS("./output/adaliftover/adaliftover_scores.rds") %>%
        mutate(id = as.character(hg38_regions[human_idx]$identifier)) %>%
        group_by(id) %>%
        summarize(adaliftover_h2m = suppressWarnings(max(AdaLiftOver_hg38_to_mm10, na.rm = T)),
                  adaliftover_m2h = suppressWarnings(max(AdaLiftOver_mm10_to_hg38, na.rm = T))) %>%
        mutate(adaliftover_h2m = ifelse(is.finite(adaliftover_h2m), adaliftover_h2m, NA),
               adaliftover_m2h = ifelse(is.finite(adaliftover_m2h), adaliftover_m2h, NA))
baseline = left_join(baseline, ada, by = "id")
seq_methods = c("phyloP4way", "phyloP7way", "phyloP17way", "phyloP20way", "phyloP30way", "phyloP100way",
                "phastCons4way", "phastCons7way", "phastCons17way", "phastCons20way", "phastCons30way", "phastCons100way",
                "percentage", "lecif", "adaliftover_h2m", "adaliftover_m2h")
num_sim = 20
loci_sample_size = 1000

# CO-V validation, leave two tissues out
tissue_pairs = tribble(~pair_id, ~tissue1, ~tissue2,
                       1, "adipose", "liver",
                       2, "heart", "muscle",
                       3, "liver", "spleen",
                       4, "adipose", "muscle",
                       5, "heart", "liver",
                       6, "spleen", "adipose",
                       7, "adipose", "heart",
                       8, "heart", "spleen",
                       9, "liver", "muscle",
                       10, "muscle", "spleen")
for (i_valid in 1:nrow(tissue_pairs)) {
        message(paste0("running valid ", i_valid))
        human_acc1 = samples_in("human", tissue_pairs$tissue1[i_valid])
        human_acc2 = samples_in("human", tissue_pairs$tissue2[i_valid])
        pig_acc1 = samples_in("pig", tissue_pairs$tissue1[i_valid])
        pig_acc2 = samples_in("pig", tissue_pairs$tissue2[i_valid])
        exclude_indices = which(anchors_all$accession1 %in% c(human_acc1, human_acc2) |
                                        anchors_all$accession2 %in% c(pig_acc1, pig_acc2))
        if (length(exclude_indices) == 0) {
                anchors = anchors_all
        } else {
                anchors = anchors_all[-exclude_indices, ]
        }
        if (nrow(anchors) < 5 || length(c(human_acc1, human_acc2)) == 0 || length(c(pig_acc1, pig_acc2)) == 0) next
        set.seed(37 + tissue_pairs$pair_id[i_valid])
        human_anchor = as.matrix(human_correct[, anchors$accession1])
        pig_anchor = as.matrix(pig_correct[, anchors$accession2])
        loci_cor = future_map_dbl(1:n_regions, function(j) {
                x = human_anchor[j, ]
                y = pig_anchor[j, ]
                if (var(x) > 0 && var(y) > 0) {
                        wCorr::weightedCorr(x, y, weights = anchors$score, method = "spearman")
                } else {
                        NA_real_
                }
        }, .options = furrr_options(seed = T))
        human_var_std = calc_var_std_ver1(human_anchor)
        pig_var_std = calc_var_std_ver1(pig_anchor)
        human_diff = rowMeans(log2(1 + human_correct[, human_acc1, drop = F])) -
                rowMeans(log2(1 + human_correct[, human_acc2, drop = F]))
        pig_diff = rowMeans(log2(1 + pig_correct[, pig_acc1, drop = F])) -
                rowMeans(log2(1 + pig_correct[, pig_acc2, drop = F]))
        score_tb = bind_cols(tibble(`CO-V` = loci_cor,
                                    `CO-V-humanmouse` = baseline$cov_chromatin_accessibility,
                                    `CO-V-humanmouse-manual` = baseline$cov_manual_chromatin_accessiblity,
                                    human_var_std = human_var_std,
                                    pig_var_std = pig_var_std),
                             baseline[seq_methods])
        bg_set_var = order(rank(human_var_std) + rank(pig_var_std), decreasing = T)[1:5e5]
        sim_df = expand.grid(method = colnames(score_tb),
                             min_rank = c(15e3, 25e3, 50e3, 75e3, 100e3),
                             sim_id = 1:num_sim,
                             KEEP.OUT.ATTRS = F, stringsAsFactors = F) %>%
                as_tibble() %>%
                arrange(method, min_rank, sim_id)
        sim_df$cor = future_map_dbl(1:nrow(sim_df), function(k) {
                method_name = sim_df$method[k]
                min_rank = sim_df$min_rank[k]
                if (method_name %in% c("CO-V", "CO-V-humanmouse", "CO-V-humanmouse-manual", "human_var_std", "pig_var_std")) {
                        bg_set = bg_set_var
                } else {
                        bg_set = 1:n_regions
                }
                bg_score = score_tb[[method_name]][bg_set]
                if (sum(!is.na(bg_score)) < min_rank) return(NA_real_)
                fg_set = sample(bg_set[order(bg_score, decreasing = T, na.last = NA)][1:min_rank], loci_sample_size)
                suppressWarnings(cor(human_diff[fg_set], pig_diff[fg_set], method = "spearman"))
        }, .options = furrr_options(seed = T))
        sim_df = mutate(sim_df, pair_id = tissue_pairs$pair_id[i_valid],
                        T1 = tissue_pairs$tissue1[i_valid], T2 = tissue_pairs$tissue2[i_valid])
        saveRDS(sim_df, output_rds_path(paste0("res_cov_", i_valid)))
}
valid_res_tb = bind_rows(map(list.files(script_output_dir, pattern = "^res_cov_", full.names = T), readRDS))
write_tsv(valid_res_tb, paste0(script_output_dir, "cov_loto.tsv"))

# CO-B validation, leave one tissue out
tissues = c("adipose", "heart", "liver", "muscle", "spleen")
for (i_valid in 1:length(tissues)) {
        message(paste0("running valid ", i_valid))
        human_acc = samples_in("human", tissues[i_valid])
        pig_acc = samples_in("pig", tissues[i_valid])
        exclude_indices = which(anchors_all$accession1 %in% human_acc | anchors_all$accession2 %in% pig_acc)
        if (length(exclude_indices) == 0) {
                anchors = anchors_all
        } else {
                anchors = anchors_all[-exclude_indices, ]
        }
        if (nrow(anchors) < 5 || length(human_acc) == 0 || length(pig_acc) == 0) next
        set.seed(37 + i_valid)
        prob_vec = seq(0.1, 1.0, by = 0.01)
        human_anchor = as.matrix(human_correct[, anchors$accession1])
        pig_anchor = as.matrix(pig_correct[, anchors$accession2])
        sample_indices = sample(n_regions, 1e5)
        human_quantile_vec = quantile(human_anchor[sample_indices, ], prob_vec, na.rm = T)
        pig_quantile_vec = quantile(pig_anchor[sample_indices, ], prob_vec, na.rm = T)
        loci_hval = compute_hval(human_anchor, pig_anchor, human_quantile_vec, pig_quantile_vec, prob_vec)
        human_mean = rowMeans(log2(1 + human_correct[, human_acc, drop = F]))
        pig_mean = rowMeans(log2(1 + pig_correct[, pig_acc, drop = F]))
        human_cutoff = quantile(human_mean, 0.95, na.rm = T)
        pig_cutoff = quantile(pig_mean, 0.95, na.rm = T)
        score_tb = bind_cols(tibble(`CO-B` = loci_hval,
                                    `CO-B-humanmouse` = baseline$cob_chromatin_accessibility,
                                    human_mean = rowMeans(human_anchor),
                                    pig_mean = rowMeans(pig_anchor)),
                             baseline[seq_methods])
        sim_df = expand.grid(method = colnames(score_tb),
                             min_rank = c(1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2e5),
                             sim_id = 1:num_sim,
                             KEEP.OUT.ATTRS = F, stringsAsFactors = F) %>%
                as_tibble() %>%
                arrange(method, min_rank, sim_id)
        sim_df$jaccard = future_map_dbl(1:nrow(sim_df), function(k) {
                bg_score = score_tb[[sim_df$method[k]]]
                min_rank = sim_df$min_rank[k]
                if (sum(!is.na(bg_score)) < min_rank) return(NA_real_)
                fg_set = sample(order(bg_score, decreasing = T, na.last = NA)[1:min_rank], loci_sample_size)
                union_n = sum(human_mean[fg_set] > human_cutoff | pig_mean[fg_set] > pig_cutoff)
                if (union_n == 0) return(NA_real_)
                sum(human_mean[fg_set] > human_cutoff & pig_mean[fg_set] > pig_cutoff) / union_n
        }, .options = furrr_options(seed = T))
        sim_df = mutate(sim_df, tissue_id = i_valid, tissue = tissues[i_valid])
        saveRDS(sim_df, output_rds_path(paste0("res_cob_", i_valid)))
}
valid_res_tb = bind_rows(map(list.files(script_output_dir, pattern = "^res_cob_", full.names = T), readRDS))
write_tsv(valid_res_tb, paste0(script_output_dir, "cob_loto.tsv"))
