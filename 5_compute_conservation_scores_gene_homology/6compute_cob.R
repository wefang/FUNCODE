## Calculate GH-CO-B for human-mouse DHS pairs near homologous genes and call significant pairs
library(tidyverse)
library(data.table)

modality = commandArgs(trailingOnly = TRUE)[1] # chromatin, H3K27ac, H3K4me1, H3K4me3
script_output_dir = paste0("./output/MTD-3_extended_", ifelse(modality == "chromatin", "cob", modality), "/")
if (!dir.exists(script_output_dir)) dir.create(script_output_dir)
file_prefix = ifelse(modality == "chromatin", "cob_", "")

# h-value of each element, separately for each species
compute_hval_single <- function(mat, quantile_vec, prob_vec) {
        hval = numeric(nrow(mat))
        names(hval) = rownames(mat)
        for (i in seq_along(prob_vec)) {
                hval[rowMeans(mat > quantile_vec[i]) >= prob_vec[i]] = prob_vec[i]
        }
        hval
}
prob_vec = seq(0.01, 1.0, by = 0.01)

# load matched samples
if (modality == "chromatin") {
        anchor_df = read_tsv("./output/chromatin_conservation/anchors.tsv")
        human_file = "./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds"
        mouse_file = "./intermediate_data/encodev4_all_dhs_mouse_chromatin_correct.rds"
} else {
        anchor_df = read_tsv(paste0("./output/histone_conservation_v1/anchors_", modality, "_v1.tsv"))
        human_file = paste0("./non_alignable_env/proc_data/human_", modality, ".rds")
        mouse_file = paste0("./non_alignable_env/proc_data/mouse_", modality, ".rds")
}

# human h-value
human_mat = readRDS(human_file)
anchor_df = filter(anchor_df, accession1 %in% colnames(human_mat))
set.seed(42)
human_quantile_vec = quantile(as.vector(human_mat[sample(nrow(human_mat), 1e5), anchor_df$accession1]), prob_vec)
human_hval = compute_hval_single(human_mat[, anchor_df$accession1, drop = F], human_quantile_vec, prob_vec)
saveRDS(human_hval, paste0(script_output_dir, "human_element_hval.rds"))
rm(human_mat)

# mouse h-value
mouse_mat = readRDS(mouse_file)
anchor_df = filter(anchor_df, accession2 %in% colnames(mouse_mat))
if (modality == "chromatin") set.seed(42)
mouse_quantile_vec = quantile(as.vector(mouse_mat[sample(nrow(mouse_mat), 1e5), anchor_df$accession2]), prob_vec)
mouse_hval = compute_hval_single(mouse_mat[, anchor_df$accession2, drop = F], mouse_quantile_vec, prob_vec)
saveRDS(mouse_hval, paste0(script_output_dir, "mouse_element_hval.rds"))
rm(mouse_mat)

# single null of randomly paired elements
set.seed(123)
null_human = sample(names(human_hval), 1e6, replace = T)
null_mouse = sample(names(mouse_hval), 1e6, replace = T)
null_cob = pmin(human_hval[null_human], mouse_hval[null_mouse])
null_ecdf = ecdf(null_cob)
cob_threshold = quantile(null_cob, 1 - 0.25)

# score pairs within each homologous gene set, keep pairs with p < 0.25 and count all pairs
load("./output/all_dhs_pairs.rda")
human_dhs_info = readRDS("./non_alignable_env/annot/human_dhs_info.rds")
mouse_dhs_info = readRDS("./non_alignable_env/annot/mouse_dhs_info.rds")
human_dist_lookup = setNames(human_dhs_info$human_dist_group, human_dhs_info$human_id)
mouse_dist_lookup = setNames(mouse_dhs_info$mouse_dist_group, mouse_dhs_info$mouse_id)
total_pairs = 0
pval_list = list()
for (i in 1:nrow(hidset_pairs)) {
        if (i %% 1000 == 0) message(i)
        human_idx = which(human_hid_set$hid_set == hidset_pairs$human_hid_set[i])
        mouse_idx = which(mouse_hid_set$hid_set == hidset_pairs$mouse_hid_set[i])
        if (length(human_idx) == 0 || length(mouse_idx) == 0) next
        human_ids = human_hid_set$human_dhs[[human_idx]]$id
        mouse_ids = mouse_hid_set$mouse_dhs[[mouse_idx]]$id
        human_ids = human_ids[human_ids %in% names(human_hval)]
        mouse_ids = mouse_ids[mouse_ids %in% names(mouse_hval)]
        total_pairs = total_pairs + length(human_ids) * length(mouse_ids)
        human_ids = human_ids[human_hval[human_ids] >= cob_threshold]
        mouse_ids = mouse_ids[mouse_hval[mouse_ids] >= cob_threshold]
        if (length(human_ids) == 0 || length(mouse_ids) == 0) next
        pairs = as.data.table(expand.grid(human_id = human_ids, mouse_id = mouse_ids, stringsAsFactors = F))
        pairs[, cob := pmin(human_hval[human_id], mouse_hval[mouse_id])]
        pairs = pairs[cob >= cob_threshold]
        if (nrow(pairs) > 0) {
                pairs[, pval := 1 - null_ecdf(cob)]
                pairs[, `:=`(human_hid_set = hidset_pairs$human_hid_set[i],
                             mouse_hid_set = hidset_pairs$mouse_hid_set[i])]
                pairs[, human_dist_group := human_dist_lookup[human_id]]
                pairs[, mouse_dist_group := mouse_dist_lookup[mouse_id]]
                pairs[, group := paste0(human_dist_group, mouse_dist_group)]
                pval_list[[length(pval_list) + 1]] = pairs
        }
}
pval_tb = rbindlist(pval_list, fill = T)
saveRDS(pval_tb, paste0(script_output_dir, file_prefix, "pval_filtered.rds"))

# BH over all pairs, including the pairs not kept
setorder(pval_tb, pval)
pval_tb[, rank := as.numeric(1:nrow(pval_tb))]
pval_tb[, fdr := pmin(1, pval * total_pairs / rank)]
pval_tb[, fdr := rev(cummin(rev(fdr)))]
saveRDS(pval_tb[fdr <= 0.25], paste0(script_output_dir, file_prefix, "significant_fdr25.rds"))
saveRDS(pval_tb[fdr <= 0.10], paste0(script_output_dir, file_prefix, "significant_fdr10.rds"))
