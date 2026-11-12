library(tidyverse)
library(qvalue)
source("./helper/helper.R")
script_output_dir = "./output/human_pig_conservation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
# core set pairs, batch corrected human and normalized pig data
regions = readRDS("./intermediate_data/human_pig_mapped_regions.rds")
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds")
regions = regions[regions$identifier %in% rownames(human_correct), ]
human_correct = human_correct[regions$identifier, ]
pig_correct = readRDS("./intermediate_data/pig_chromatin_correct.rds")[regions$identifier, ]

# anchors with score >= 0.5 for CO-V and >= 0.7 for CO-B
anchors_all = distinct(read_tsv(paste0(script_output_dir, "anchors.tsv")))
anchors_all = filter(anchors_all, accession1 %in% colnames(human_correct), accession2 %in% colnames(pig_correct))
anchor_df = filter(anchors_all, score >= 0.5)
anchor_df_cob = filter(anchors_all, score >= 0.7)

#### CO-V ####
# weighted pearson correlation of within-region ranks
weighted_spearman = function(x, y, w) {
        if (var(x) > 0 & var(y) > 0) {
                rx = rank(x)
                ry = rank(y)
                w = w / sum(w)
                rx = rx - sum(w * rx)
                ry = ry - sum(w * ry)
                den = sqrt(sum(w * rx^2)) * sqrt(sum(w * ry^2))
                if (den > 0) {
                        return(sum(w * rx * ry) / den)
                }
        }
        NA_real_
}
human_anchor = as.matrix(human_correct[, anchor_df$accession1])
pig_anchor = as.matrix(pig_correct[, anchor_df$accession2])
loci_cor = map_dbl(1:nrow(human_anchor), function(j) {
        weighted_spearman(human_anchor[j, ], pig_anchor[j, ], anchor_df$score)
})
loci_cor_unweighted = map_dbl(1:nrow(human_anchor), function(j) {
        weighted_spearman(human_anchor[j, ], pig_anchor[j, ], rep(1, nrow(anchor_df)))
})
# pooled null from random region pairs
set.seed(42)
null_size = 1e5
sample_indices1 = sample(nrow(human_anchor), null_size, replace = T)
sample_indices2 = sample(nrow(pig_anchor), null_size, replace = T)
cor_null = map_dbl(1:null_size, function(j) {
        weighted_spearman(human_anchor[sample_indices1[j], ], pig_anchor[sample_indices2[j], ], anchor_df$score)
})
spearman_pval = rep(NA_real_, length(loci_cor))
spearman_pval[!is.na(loci_cor)] = empPvals(loci_cor[!is.na(loci_cor)], cor_null[!is.na(cor_null)])
spearman_pval_adj = p.adjust(spearman_pval, method = "fdr")

#### CO-B ####
prob_vec = seq(0.01, 1.0, by = 0.01)
compute_hval_single = function(mat, quantile_vec) {
        hval = numeric(nrow(mat))
        for (i in seq_along(prob_vec)) {
                hval[rowMeans(mat > quantile_vec[i]) >= prob_vec[i]] = prob_vec[i]
        }
        hval
}
human_anchor = as.matrix(human_correct[, anchor_df_cob$accession1])
pig_anchor = as.matrix(pig_correct[, anchor_df_cob$accession2])
set.seed(42)
human_quantile_vec = quantile(human_anchor[sample(nrow(human_anchor), 1e5), ], prob_vec)
set.seed(42)
pig_quantile_vec = quantile(pig_anchor[sample(nrow(pig_anchor), 1e5), ], prob_vec)
human_hval = compute_hval_single(human_anchor, human_quantile_vec)
pig_hval = compute_hval_single(pig_anchor, pig_quantile_vec)
loci_hval = pmin(human_hval, pig_hval)
# activity filter, >= 50 anchors above the 80th percentile of anchor means in each species
n_h_active = as.integer(rowSums(human_anchor > quantile(rowMeans(human_anchor), 0.8, na.rm = T)))
n_p_active = as.integer(rowSums(pig_anchor > quantile(rowMeans(pig_anchor), 0.8, na.rm = T)))
active = n_h_active >= 50 & n_p_active >= 50
# pooled null from random region pairs, same filter
set.seed(42)
null_size = 1e6
sample_indices1 = sample(nrow(human_anchor), null_size, replace = T)
sample_indices2 = sample(nrow(pig_anchor), null_size, replace = T)
hval_null = pmin(human_hval[sample_indices1], pig_hval[sample_indices2])
null_active = n_h_active[sample_indices1] >= 50 & n_p_active[sample_indices2] >= 50
hval_pval = rep(NA_real_, length(loci_hval))
hval_pval[active] = empPvals(loci_hval[active], hval_null[null_active])
hval_pval[active & loci_hval == 0] = 1
hval_pval_adj = p.adjust(hval_pval, method = "fdr")

scores = tibble(identifier = regions$identifier,
                source = regions$source,
                cov_all_weighted = loci_cor,
                cov_all_unweighted = loci_cor_unweighted,
                human_hval = human_hval,
                pig_hval = pig_hval,
                cob_all = ifelse(active, loci_hval, NA),
                n_h_active = n_h_active,
                n_p_active = n_p_active,
                cov_pval = spearman_pval,
                cov_fdr = spearman_pval_adj,
                cov_sign = as.integer(spearman_pval_adj < 0.1),
                cob_pval = hval_pval,
                cob_fdr = hval_pval_adj,
                cob_sign = as.integer(hval_pval_adj < 0.1),
                hg38_summit = regions$hg38_summit,
                susScr11_summit = regions$susScr11_summit)
saveRDS(scores, paste0(script_output_dir, "conservation_scores.rds"))
write_tsv(scores, paste0(script_output_dir, "conservation_scores.tsv"))
save(cor_null, hval_null, null_active, file = paste0(script_output_dir, "null.rda"))
