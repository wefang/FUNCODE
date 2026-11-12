library(tidyverse)
library(Seurat)
library(ggpubr)
library(ComplexHeatmap)
region_file = "./output/indexed_dhs_mapped_regions.rds"
source("./helper/helper.R")
source("./helper/def_color.R")
script_output_dir = "./output/histone_conservation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
script_plot_dir = "./plots/histone_chromatin/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}
regions = readRDS(region_file)
regions$human_dist_cat = factor(cut(regions$dist_tss_hg38, c(-Inf, 200, 2000, Inf)), labels = c("promoter", "proximal", "distal"))
regions$mouse_dist_cat = cut(regions$dist_tss_mm10, c(-Inf, 200, 2000, Inf), labels = c("promoter", "proximal", "distal"))
regions$dist_pair_cat = paste0(regions$human_dist_cat, "_", regions$mouse_dist_cat)
null_tb = as_tibble(expand.grid(human_dist_cat = c("promoter", "proximal", "distal"),
                                mouse_dist_cat = c("promoter", "proximal", "distal")))
null_tb$dist_pair_cat = paste0(null_tb$human_dist_cat, "_", null_tb$mouse_dist_cat)
set.seed(37)
null_size = 5e5
for (assay in c("H3K4me1", "H3K4me3", "H3K27ac")) {
        message(assay)
        human_histone = readRDS(paste0("E:\\GlobusDownload/histone_norm_correct/Homo_sapiens_", assay, "_norm_correct.rds"))
        mouse_histone = readRDS(paste0("E:\\GlobusDownload/histone_norm_correct/Mus_musculus_", assay, "_norm_correct.rds"))
        anchors = readRDS(paste0("./output/histone_conservation_v1/", "anchors_", assay, "_v1.rds"))
        all(anchors$accession1 %in% colnames(human_histone))
        all(anchors$accession2 %in% colnames(mouse_histone))
        # anchors = anchors[anchors$score > 0.5, ]
        # calling conserved elements
        loci_cor = map_dbl(1:nrow(human_histone), function(j) {
                if (j %% 1e5 == 0) message(j)
                x = human_histone[j, anchors$accession1]
                y = mouse_histone[j, anchors$accession2]
                if (var(x) > 0 & var(y) > 0) {
                        wCorr::weightedCorr(x,
                                            y,
                                            weights = anchors$score,
                                            method = "spearman")
                } else {
                        return(NA)
                }
        })
        null_tb$spearman_null = map(1:nrow(null_tb), function(i) {
                message(i)
                human_cat = null_tb$human_dist_cat[i]
                mouse_cat = null_tb$mouse_dist_cat[i]
                sample_indices1 = sample(which(regions$human_dist_cat == human_cat), null_size, replace = T)
                sample_indices2 = sample(which(regions$mouse_dist_cat == mouse_cat), null_size, replace = T)
                human_null = human_histone[sample_indices1, ]
                mouse_null = mouse_histone[sample_indices2, ]
                loci_cor_null = map_dbl(1:null_size, function(j) {
                        x = human_null[j, anchors$accession1]
                        y = mouse_null[j, anchors$accession2]
                        if (var(x) > 0 & var(y) > 0) {
                                wCorr::weightedCorr(x,
                                                    y,
                                                    weights = anchors$score,
                                                    method = "spearman")
                        } else {
                                return(NA)
                        }
                })
                loci_cor_null
        })
        null_tb$spearman = map(1:nrow(null_tb), function(i) {
                cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
                out = loci_cor[cat_indices]
                names(out) = regions$id[cat_indices]
                out
        })
        null_tb$spearman_pval = map(1:nrow(null_tb), function(i) {
                cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
                spearman_pval_cat = qvalue::empPvals(loci_cor[cat_indices], null_tb$spearman_null[[i]])
                names(spearman_pval_cat) = regions$id[cat_indices]
                spearman_pval_cat
        })
        null_tb$spearman_pval_adj = map(null_tb$spearman_pval, function(x) {
                p.adjust(x, method = "fdr")
        })
        prob_vec = seq(0.01, 1.0, by = 0.01)
        human_quantile_vec = quantile(human_histone[sample(nrow(human_histone), 1e5), anchors$accession1], prob_vec)
        mouse_quantile_vec = quantile(mouse_histone[sample(nrow(mouse_histone), 1e5), anchors$accession2], prob_vec)
        loci_hval = compute_hval(human_histone[, anchors$accession1],
                                 mouse_histone[, anchors$accession2],
                                 human_quantile_vec, mouse_quantile_vec, prob_vec)
        null_tb$hval_null = map(1:nrow(null_tb), function(i) {
                message(i)
                human_cat = null_tb$human_dist_cat[i]
                mouse_cat = null_tb$mouse_dist_cat[i]
                null_size = 1e5
                sample_indices1 = sample(which(regions$human_dist_cat == human_cat), null_size, replace = T)
                sample_indices2 = sample(which(regions$mouse_dist_cat == mouse_cat), null_size, replace = T)
                human_null = human_histone[sample_indices1, ]
                mouse_null = mouse_histone[sample_indices2, ]
                
                hval = compute_hval(human_null[, anchors$accession1],
                                    mouse_null[, anchors$accession2],
                                    human_quantile_vec, mouse_quantile_vec, prob_vec)
                hval
        })
        null_tb$hval = map(1:nrow(null_tb), function(i) {
                cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
                out = loci_hval[cat_indices]
                names(out) = regions$id[cat_indices]
                out
        })
        null_tb$hval_pval = map(1:nrow(null_tb), function(i) {
                cat_indices = which(regions$dist_pair_cat == null_tb$dist_pair_cat[i])
                hval_pval_cat = qvalue::empPvals(loci_hval[cat_indices], null_tb$hval_null[[i]])
                names(hval_pval_cat) = regions$id[cat_indices]
                hval_pval_cat
        })
        null_tb$hval_pval_adj = map(null_tb$hval_pval, function(x) {
                p.adjust(x, method = "fdr")
        })
        sum(map_dbl(null_tb$hval_pval_adj, function(x) sum(x < 0.1)))
        save(loci_cor, loci_hval, null_tb, file = paste0(script_output_dir, assay, "_scores.rda"))
}
print(load(paste0("./output/chromatin_conservation/scores.rda")))
cor_null_all_cat = bind_rows(map(1:nrow(null_tb), function(i) {
        cat_indices = which(regions_all$dist_pair_cat == null_tb$dist_pair_cat[i])
        cat_indices = cat_indices[cat_indices <= length(loci_cor)]
        tibble(score = c(loci_cor[cat_indices],
                         null_tb$spearman_null[[i]]),
               source = c(rep("ZAlt", length(cat_indices)),
                          rep("Null", length(null_tb$spearman_null[[i]]))),
               category = null_tb$dist_pair_cat[i])
}))
null_tb$spearman_cutoff = map_dbl(1:nrow(null_tb), function(i) {
        cat_indices = which(regions_all$dist_pair_cat == null_tb$dist_pair_cat[i])
        cat_indices = cat_indices[cat_indices <= length(loci_cor)]
        min(loci_cor[cat_indices][spearman_pval_adj[cat_indices] < 0.1])
})
((cor_null_all_cat %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "CO-V (CA)") +
                geom_vline(aes(xintercept = spearman_cutoff), color = "#99154e",
                           data = transmute(null_tb, category = dist_pair_cat, spearman_cutoff))) %>%
                facet(facet.by = "category") +
                scale_color_manual(values = c("#868686", "#99154e"))+
                scale_fill_manual(values = c("#868686", "#99154e")) + theme(text = element_text(size = 10), legend.position = "none")) %>%
        push_pdf(file_name = "hist_cov_ca", w = 5., h = 3.)

hval_null_all_cat = bind_rows(map(1:nrow(null_tb), function(i) {
        cat_indices = which(regions_all$dist_pair_cat == null_tb$dist_pair_cat[i])
        cat_indices = cat_indices[cat_indices <= length(loci_cor)]
        tibble(score = c(loci_hval[cat_indices],
                         null_tb$hval_null[[i]]),
               source = c(rep("ZAlt", length(cat_indices)),
                          rep("Null", length(null_tb$hval_null[[i]]))),
               category = null_tb$dist_pair_cat[i])
}))
null_tb$hval_cutoff = map_dbl(1:nrow(null_tb), function(i) {
        cat_indices = which(regions_all$dist_pair_cat == null_tb$dist_pair_cat[i])
        cat_indices = cat_indices[cat_indices <= length(loci_cor)]
        min(loci_hval[cat_indices][hval_pval_adj[cat_indices] < 0.1])
})
((hval_null_all_cat %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "CO-B (CA)") +
                geom_vline(aes(xintercept = hval_cutoff), color = "#ffd56b",
                           data = transmute(null_tb, category = dist_pair_cat, hval_cutoff))) %>%
        facet(facet.by = "category") +
                scale_color_manual(values = c("#868686", "#ffd56b"))+
                scale_fill_manual(values = c("#868686", "#ffd56b")) +
                theme(text = element_text(size = 10), legend.position = "none")) %>%
push_pdf(file_name = "hist_cob_ca", w = 5.0, h = 3.0)

assay = "H3K4me3"
print(load(paste0("./output/histone_conservation_v2/", assay, "_scores.rda")))
cor_null_all_cat = bind_rows(map(1:nrow(null_tb), function(i) {
        tibble(score = c(null_tb$spearman[[i]],
                         null_tb$spearman_null[[i]]),
               source = c(rep("ZAlt", length(null_tb$spearman[[i]])),
                          rep("Null", length(null_tb$spearman_null[[i]]))),
               category = null_tb$dist_pair_cat[i])
}))
null_tb$spearman_cutoff = map_dbl(1:nrow(null_tb), function(i) {
        min(null_tb$spearman[[i]][null_tb$spearman_pval_adj[[i]] < 0.1])
})
((cor_null_all_cat %>% 
                ggdensity(x = "score", color = "source", fill = "source", xlab = paste0("CO-V (", assay, ")")) +
                geom_vline(aes(xintercept = spearman_cutoff), color = "#99154e",
                           data = transmute(null_tb, category = dist_pair_cat, spearman_cutoff))) %>%
                facet(facet.by = "category") +
                scale_color_manual(values = c("#868686", "#99154e"))+
                scale_fill_manual(values = c("#868686", "#99154e")) + theme(text = element_text(size = 10), legend.position = "none")) %>%
       push_pdf(file_name = paste0("hist_cov_", assay), w = 5., h = 3.)

hval_null_all_cat = bind_rows(map(1:nrow(null_tb), function(i) {
        tibble(score = c(null_tb$hval[[i]],
                         null_tb$hval_null[[i]]),
               source = c(rep("ZAlt", length(null_tb$hval[[i]])),
                          rep("Null", length(null_tb$hval_null[[i]]))),
               category = null_tb$dist_pair_cat[i])
}))
null_tb$hval_cutoff = map_dbl(1:nrow(null_tb), function(i) {
        min(null_tb$hval[[i]][null_tb$hval_pval_adj[[i]] < 0.1])
})
((hval_null_all_cat %>% 
          ggdensity(x = "score", color = "source", fill = "source", xlab = paste0("CO-V (", assay, ")")) +
          geom_vline(aes(xintercept = hval_cutoff), color = "#ffd56b",
                     data = transmute(null_tb, category = dist_pair_cat, hval_cutoff))) %>%
                facet(facet.by = "category") +
                scale_color_manual(values = c("#868686", "#ffd56b"))+
                scale_fill_manual(values = c("#868686", "#ffd56b")) + theme(text = element_text(size = 10), legend.position = "none")) %>%
        push_pdf(file_name = paste0("hist_cob_", assay), w = 5., h = 3.)

#### within-species baselines: mean and standardized variance over anchor samples ####
for (assay in c("H3K27ac", "H3K4me1", "H3K4me3")) {
        human_histone = readRDS(paste0("./processed_data/histone_norm_correct/Homo_sapiens_", assay, "_norm_correct.rds"))
        mouse_histone = readRDS(paste0("./processed_data/histone_norm_correct/Mus_musculus_", assay, "_norm_correct.rds"))
        anchors_baseline = readRDS(paste0(script_output_dir, "anchors_", assay, "_v1.rds"))
        anchors_baseline = anchors_baseline[anchors_baseline$score > 0.5, ]
        anchors_baseline = anchors_baseline[anchors_baseline$accession1 %in% colnames(human_histone) &
                                                    anchors_baseline$accession2 %in% colnames(mouse_histone), ]
        human_anchor_mat = as.matrix(human_histone[, anchors_baseline$accession1, drop = F])
        mouse_anchor_mat = as.matrix(mouse_histone[, anchors_baseline$accession2, drop = F])
        baseline_tb = tibble(id = regions$id,
                             human_mean = rowMeans(log2(human_anchor_mat + 1)),
                             human_var_std = calc_var_std_ver1(human_anchor_mat),
                             mouse_mean = rowMeans(log2(mouse_anchor_mat + 1)),
                             mouse_var_std = calc_var_std_ver1(mouse_anchor_mat))
        names(baseline_tb)[-1] = paste0(c("human_", "human_", "mouse_", "mouse_"), assay,
                                        c("_mean", "_var_std", "_mean", "_var_std"))
        saveRDS(baseline_tb, paste0(script_output_dir, assay, "_species_baselines.rds"))
}






















