library(tidyverse)
library(Seurat)
library(ggpubr)
library(ComplexHeatmap)
region_file = "./output/indexed_dhs_mapped_regions.rds"
source("./helper/helper.R")
source("./helper/def_color.R")
script_output_dir = "./output/chromatin_conservation_v1/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
script_plot_dir = "./plots/conservation_chromatin_v1/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}
# batch corrected data
human_correct = readRDS("./intermediate_data/human_chromatin_correct.rds")
mouse_correct = readRDS("./intermediate_data/mouse_chromatin_correct.rds")
regions = readRDS(region_file)

anchor_df = read_tsv(paste0(script_output_dir, "anchors.tsv"))
script_output_dir = "./output/chromatin_conservation_v1/"
regions$human_dist_cat = factor(cut(regions$dist_tss_hg38, c(-Inf, 200, 2000, Inf)), labels = c("promoter", "proximal", "distal"))
regions$mouse_dist_cat = cut(regions$dist_tss_mm10, c(-Inf, 200, 2000, Inf), labels = c("promoter", "proximal", "distal"))
regions$dist_pair_cat = paste0(regions$human_dist_cat, "_", regions$mouse_dist_cat)
# assign category of DHS based on distance
null_tb = as_tibble(expand.grid(human_dist_cat = c("promoter", "proximal", "distal"),
                                mouse_dist_cat = c("promoter", "proximal", "distal")))
null_tb$dist_pair_cat = paste0(null_tb$human_dist_cat, "_", null_tb$mouse_dist_cat)
# computing weighted spearman correlation
loci_cor = sapply(1:nrow(human_correct), function(j) {
        if (j %% 1e5 == 0) message(j)
        x = human_correct[j, anchor_df$accession1]
        y = mouse_correct[j, anchor_df$accession2]
        if (var(x) > 0 & var(y) > 0) {
                wCorr::weightedCorr(x,
                                    y,
                                    weights = anchor_df$score,
                                    method = "spearman")
        } else {
                return(NA)
        }
})
null_tb$spearman_null = map(1:nrow(null_tb), function(i) {
        message(i)
        human_cat = null_tb$human_dist_cat[i]
        mouse_cat = null_tb$mouse_dist_cat[i]
        null_size = 1e5
        sample_indices1 = sample(which(regions$human_dist_cat == human_cat), null_size, replace = T)
        sample_indices2 = sample(which(regions$mouse_dist_cat == mouse_cat), null_size, replace = T)
        human_null = human_correct[sample_indices1, ]
        mouse_null = mouse_correct[sample_indices2, ]
        loci_cor_null = map_dbl(1:null_size, function(j) {
                x = human_null[j, anchor_df$accession1]
                y = mouse_null[j, anchor_df$accession2]
                if (var(x) > 0 & var(y) > 0) {
                        wCorr::weightedCorr(x,
                                            y,
                                            weights = anchor_df$score,
                                            method = "spearman")
                } else {
                        return(NA)
                }
        })
        loci_cor_null
})
prob_vec = seq(0.01, 1.0, by = 0.01)
human_quantile_vec = quantile(human_correct[sample(nrow(human_correct), 1e5), anchor_df$accession1], prob_vec)
mouse_quantile_vec = quantile(mouse_correct[sample(nrow(mouse_correct), 1e5), anchor_df$accession2], prob_vec)
loci_hval = compute_hval(human_correct[, anchor_df$accession1],
                         mouse_correct[, anchor_df$accession2],
                         human_quantile_vec, mouse_quantile_vec, prob_vec)
loci_hval
null_tb$hval_null = map(1:nrow(null_tb), function(i) {
        message(i)
        human_cat = null_tb$human_dist_cat[i]
        mouse_cat = null_tb$mouse_dist_cat[i]
        null_size = 1e5
        sample_indices1 = sample(which(regions$human_dist_cat == human_cat), null_size, replace = T)
        sample_indices2 = sample(which(regions$mouse_dist_cat == mouse_cat), null_size, replace = T)
        human_null = human_correct[sample_indices1, ]
        mouse_null = mouse_correct[sample_indices2, ]
        
        hval = compute_hval(human_null[, anchor_df$accession1],
                            mouse_null[, anchor_df$accession2],
                            human_quantile_vec, mouse_quantile_vec, prob_vec)
        hval
})
cor_null_all_cat = bind_rows(map(1:nrow(null_tb), function(i) {
        cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
        tibble(score = c(loci_cor[cat_indices],
                         null_tb$spearman_null[[i]]),
               source = c(rep("ZAlt", length(cat_indices)),
                          rep("Null", length(null_tb$spearman_null[[i]]))),
               category = null_tb$dist_pair_cat[i])
}))
null_tb$spearman_cutoff = map_dbl(1:nrow(null_tb), function(i) {
        cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
        min(loci_cor[cat_indices][spearman_pval_adj[cat_indices] < 0.1])
})

((cor_null_all_cat %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "Variable Conservation Score") +
                geom_vline(aes(xintercept = spearman_cutoff), color = "#e66101",
                           data = transmute(null_tb, category = dist_pair_cat, spearman_cutoff))) %>%
                facet(facet.by = "category") +
                scale_color_manual(values = c("#868686", "#e66101"))+
                scale_fill_manual(values = c("#868686", "#e66101")) + theme(text = element_text(size = 10), legend.position = "none")) %>%
        push_pdf(file_name = "hist_spearman", w = 6., h = 6.)

(cor_null_all_cat %>%
                ggdensity(x = "score", color = "category", fill = NA, xlab = "Variable Conservation Score",
                          facet.by = "source") + 
                theme(text = element_text(size = 10), legend.position = "top"))  %>%
        push_pdf(file_name = "hist_spearman_cat", w = 6., h = 6.)

hval_null_all_cat = bind_rows(map(1:nrow(null_tb), function(i) {
        cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
        tibble(score = c(loci_hval[cat_indices],
                         null_tb$hval_null[[i]]),
               source = c(rep("ZAlt", length(cat_indices)),
                          rep("Null", length(null_tb$hval_null[[i]]))),
               category = null_tb$dist_pair_cat[i])
}))
((hval_null_all_cat %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "Baseline Conservation Score") +
                geom_vline(aes(xintercept = hval_cutoff), color = "#5e3c99",
                           data = transmute(null_tb, category = dist_pair_cat, hval_cutoff))) %>%
        facet(facet.by = "category") +
                scale_color_manual(values = c("#868686", "#5e3c99"))+
                scale_fill_manual(values = c("#868686", "#5e3c99")) +
                theme(text = element_text(size = 10), legend.position = "none")) %>%
push_pdf(file_name = "hist_hval", w = 6.0, h = 6.0)

(hval_null_all_cat %>%
                ggdensity(x = "score", color = "category", fill = NA, xlab = "Baseline Conservation Score",
                          facet.by = "source") + 
                theme(text = element_text(size = 10), legend.position = "top")) %>%
        push_pdf(file_name = "hist_hval_cat", w = 6.0, h = 6.0)

# pool and adjust
spearman_pval = numeric(nrow(regions))
walk(1:nrow(null_tb), function(i) {
        cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
        spearman_pval[cat_indices] <<- qvalue::empPvals(loci_cor[cat_indices], null_tb$spearman_null[[i]])
})
spearman_pval_adj = p.adjust(spearman_pval, method = "fdr")
sum(spearman_pval_adj < 0.1)
table(spearman_pval_adj < 0.1, regions$human_cre_class_prio_v4)

hval_pval = numeric(nrow(regions))
walk(1:nrow(null_tb), function(i) {
        cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
        hval_pval[cat_indices] <<- qvalue::empPvals(loci_hval[cat_indices], null_tb$hval_null[[i]])
})
hval_pval_adj = p.adjust(hval_pval, method = "fdr")
sum(hval_pval_adj < 0.1)
table(hval_pval_adj < 0.1, regions$human_cre_class_prio_v4)

null_tb$hval_cutoff = map_dbl(1:nrow(null_tb), function(i) {
        cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
        min(loci_hval[cat_indices][hval_pval_adj[cat_indices] < 0.1])
})
save(loci_cor, loci_hval, null_tb,
     spearman_pval, spearman_pval_adj,
     hval_pval, hval_pval_adj, file = paste0(script_output_dir, "scores.rda"))

#### within-species baselines: mean and standardized variance over anchor samples ####
anchors_baseline = read_tsv("./output/chromatin_conservation/anchors.tsv")
anchors_baseline = anchors_baseline[anchors_baseline$score > 0.5, ]
anchors_baseline = anchors_baseline[anchors_baseline$accession1 %in% colnames(human_correct) &
                                            anchors_baseline$accession2 %in% colnames(mouse_correct), ]
human_anchor_mat = as.matrix(human_correct[, anchors_baseline$accession1, drop = F])
mouse_anchor_mat = as.matrix(mouse_correct[, anchors_baseline$accession2, drop = F])
baseline_tb = tibble(id = regions$id,
                     human_CA_mean = rowMeans(log2(human_anchor_mat + 1)),
                     human_CA_var_std = calc_var_std_ver1(human_anchor_mat),
                     mouse_CA_mean = rowMeans(log2(mouse_anchor_mat + 1)),
                     mouse_CA_var_std = calc_var_std_ver1(mouse_anchor_mat))
saveRDS(baseline_tb, paste0(script_output_dir, "CA_species_baselines.rds"))
