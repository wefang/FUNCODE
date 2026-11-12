# Rscript downsample_samples.R <CA|H3K27ac|H3K4me1|H3K4me3> <fraction> <seed>, then Rscript downsample_samples.R summarize
library(tidyverse)
options(Seurat.object.assay.version = "v3")
library(Seurat)
source("./helper/helper.R")
script_output_dir = "./output/downsample_samples/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
modalities = c("CA", "H3K27ac", "H3K4me1", "H3K4me3")
args = commandArgs(trailingOnly = T)
compute_cov = function(human_mat, mouse_mat, anchor_df) {
        sapply(1:nrow(human_mat), function(j) {
                x = human_mat[j, anchor_df$accession1]
                y = mouse_mat[j, anchor_df$accession2]
                if (var(x) > 0 & var(y) > 0) {
                        wCorr::weightedCorr(x, y, weights = anchor_df$score, method = "spearman")
                } else {
                        return(NA)
                }
        })
}

#### begin downsampled integration ####
if (args[1] %in% modalities) {
        modality = args[1]
        fraction = as.numeric(args[2])
        seed = as.integer(args[3])
        regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
        exp_meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
        human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
        mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
        rownames(human_correct) = rownames(mouse_correct) = regions$id
        data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
        if (modality == "CA") {
                set.seed(seed)
                human_keep = sort(sample(colnames(human_correct), max(5, round(ncol(human_correct) * fraction))))
                mouse_keep = sort(sample(colnames(mouse_correct), max(5, round(ncol(mouse_correct) * fraction))))
                human_obj = CreateSeuratObject(counts = human_correct[data_cross_var_feature, human_keep], assay = "OpenChromatin")
                mouse_obj = CreateSeuratObject(counts = mouse_correct[data_cross_var_feature, mouse_keep], assay = "OpenChromatin")
                human_obj$Accession = colnames(human_obj)
                mouse_obj$Accession = colnames(mouse_obj)
                VariableFeatures(human_obj) = rownames(human_obj)
                VariableFeatures(mouse_obj) = rownames(mouse_obj)
                anchors = FindIntegrationAnchors(list(human = human_obj, mouse = mouse_obj),
                                                 anchor.features = 10000,
                                                 k.anchor = min(5, length(mouse_keep) - 1),
                                                 k.filter = min(10, length(mouse_keep) - 1),
                                                 k.score = min(30, length(mouse_keep) - 1),
                                                 dims = 1:min(30, length(human_keep) - 1, length(mouse_keep) - 1))
                anchor_df = make_anchor_tb(anchors, human_obj, mouse_obj)
                co_v = compute_cov(human_correct, mouse_correct, anchor_df)
        } else {
                # bridge integration with the histone samples of the mark downsampled
                human_im_mat = readRDS(paste0("./intermediate_data/histone_rev_transfer_v1/", "human_rev_transfer_", modality, ".rds"))
                mouse_im_mat = readRDS(paste0("./intermediate_data/histone_rev_transfer_v1/", "mouse_rev_transfer_", modality, ".rds"))
                if (is.null(rownames(human_im_mat)) || all(rownames(human_im_mat) == "")) {
                        rownames(human_im_mat) = rownames(mouse_im_mat) = data_cross_var_feature
                }
                human_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Homo_sapiens_", modality, "_norm_correct.rds"))
                mouse_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Mus_musculus_", modality, "_norm_correct.rds"))
                rownames(human_histone) = rownames(mouse_histone) = regions$id
                set.seed(seed)
                human_keep = sort(sample(colnames(human_im_mat), max(5, round(ncol(human_im_mat) * fraction))))
                mouse_keep = sort(sample(colnames(mouse_im_mat), max(5, round(ncol(mouse_im_mat) * fraction))))
                human_obj = CreateSeuratObject(counts = cbind(human_correct[data_cross_var_feature, ], human_im_mat[, human_keep]),
                                               assay = "OpenChromatin")
                mouse_obj = CreateSeuratObject(counts = cbind(mouse_correct[data_cross_var_feature, ], mouse_im_mat[, mouse_keep]),
                                               assay = "OpenChromatin")
                rm(human_correct, mouse_correct, human_im_mat, mouse_im_mat)
                gc()
                for (col in c("Biosample term name", "Target of assay")) {
                        human_obj[[col]] = exp_meta[[col]][match(colnames(human_obj), exp_meta$Accession)]
                        mouse_obj[[col]] = exp_meta[[col]][match(colnames(mouse_obj), exp_meta$Accession)]
                }
                # removing over-represented samples for human H3K27ac
                if (modality == "H3K27ac") {
                        human_obj = human_obj[, !human_obj$`Biosample term name` %in% c("skin epidermis", "middle frontal area 46")]
                }
                human_obj = FindVariableFeatures(human_obj, nfeatures = 100000)
                mouse_obj = FindVariableFeatures(mouse_obj, nfeatures = 100000)
                human_obj$organism = "human"
                mouse_obj$organism = "mouse"
                anchors = FindIntegrationAnchors(list(human = human_obj, mouse = mouse_obj),
                                                 anchor.features = 10000,
                                                 k.anchor = 5,
                                                 k.filter = min(20, ncol(mouse_obj) - 1),
                                                 k.score = 30,
                                                 dims = 1:min(30, ncol(human_obj) - 1, ncol(mouse_obj) - 1))
                integrated = IntegrateData(anchors, dims = 1:15, k.weight = 20)
                DefaultAssay(integrated) = "integrated"
                VariableFeatures(integrated) = anchors@anchor.features
                integrated = ScaleData(integrated, features = anchors@anchor.features)
                integrated = RunPCA(integrated, features = anchors@anchor.features, npcs = 30)
                integrated_list = SplitObject(integrated[, which(integrated$`Target of assay` == modality)], split.by = "organism")
                n_min = min(map_dbl(integrated_list, ncol))
                assay_anchors = FindIntegrationAnchors(integrated_list,
                                                       assay = c("integrated", "integrated"),
                                                       anchor.features = 10000,
                                                       reduction = "rpca",
                                                       k.anchor = min(3, n_min - 1),
                                                       k.filter = min(25, n_min - 1),
                                                       k.score = min(30, n_min - 1),
                                                       dims = 1:min(15, n_min - 1))
                integrated_list[[1]]$Accession = colnames(integrated_list[[1]])
                integrated_list[[2]]$Accession = colnames(integrated_list[[2]])
                anchor_df = make_anchor_tb(assay_anchors, data_obj1 = integrated_list[[1]], data_obj2 = integrated_list[[2]])
                anchor_df = anchor_df[anchor_df$accession1 %in% colnames(human_histone) &
                                              anchor_df$accession2 %in% colnames(mouse_histone), ]
                if (nrow(anchor_df) < 3) {
                        co_v = rep(NA_real_, nrow(human_histone))
                } else {
                        co_v = compute_cov(human_histone, mouse_histone, anchor_df)
                }
        }
        names(co_v) = regions$id
        saveRDS(list(modality = modality, fraction = fraction, seed = seed,
                     n_human = length(human_keep), n_mouse = length(mouse_keep),
                     n_anchors = nrow(anchor_df), co_v = co_v, anchor_df = anchor_df),
                paste0(script_output_dir, sprintf("scores_%s_frac%.2f_seed%d.rds", modality, fraction, seed)))
}

#### begin summary over fractions and seeds ####
if (args[1] == "summarize") {
        exp_meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
        tissue_map = read_tsv("./custom_metadata/tissue_system_mapping.tsv")
        acc_tb = exp_meta %>%
                select(Accession, biosample_term = `Biosample term name`, species = Organism) %>%
                distinct(Accession, .keep_all = T) %>%
                mutate(species = if_else(grepl("sapiens", species, ignore.case = T), "human", "mouse")) %>%
                left_join(tissue_map, by = c("biosample_term", "species"))
        human_lookup = acc_tb %>% filter(species == "human") %>% select(Accession, system_human = tissue_system)
        mouse_lookup = acc_tb %>% filter(species == "mouse") %>% select(Accession, system_mouse = tissue_system)
        score_files = list.files(script_output_dir, pattern = "^scores_(CA|H3K27ac|H3K4me1|H3K4me3)_frac.*_seed.*\\.rds$", full.names = T)
        results = map(score_files, readRDS)
        # reference CO-V, mean over seeds with all samples
        refs = map(set_names(modalities), function(mod) {
                full_results = keep(results, ~ .x$modality == mod && abs(.x$fraction - 1.0) < 1e-6)
                rowMeans(do.call(cbind, map(full_results, "co_v")), na.rm = T)
        })
        per_task = map_dfr(results, function(r) {
                valid = !is.na(refs[[r$modality]]) & !is.na(r$co_v)
                evaluable = as_tibble(r$anchor_df) %>%
                        left_join(human_lookup, by = c("accession1" = "Accession")) %>%
                        left_join(mouse_lookup, by = c("accession2" = "Accession")) %>%
                        filter(!is.na(system_human), !is.na(system_mouse))
                tibble(modality = r$modality, fraction = r$fraction, seed = r$seed,
                       n_human = r$n_human, n_mouse = r$n_mouse, n_anchors = r$n_anchors,
                       spearman_r = if (sum(valid) < 100) NA_real_ else cor(refs[[r$modality]][valid], r$co_v[valid], method = "spearman"),
                       n_evaluable = nrow(evaluable),
                       precision = if (nrow(evaluable) == 0) NA_real_ else mean(evaluable$system_human == evaluable$system_mouse))
        })
        write_csv(per_task, paste0(script_output_dir, "downsample_per_task.csv"))
        cor_summary = per_task %>%
                filter(fraction < 1.0) %>%
                group_by(modality, fraction) %>%
                summarise(mean_r = mean(spearman_r, na.rm = T), sd_r = sd(spearman_r, na.rm = T),
                          n_seeds = sum(!is.na(spearman_r)), mean_anchors = mean(n_anchors, na.rm = T), .groups = "drop") %>%
                bind_rows(per_task %>%
                                  filter(abs(fraction - 1) < 1e-6) %>%
                                  group_by(modality) %>%
                                  summarise(fraction = 1.0, mean_r = 1.0, sd_r = 0.0, n_seeds = NA_integer_,
                                            mean_anchors = mean(n_anchors), .groups = "drop")) %>%
                arrange(factor(modality, levels = modalities), fraction)
        prec_summary = per_task %>%
                group_by(modality, fraction) %>%
                summarise(mean_precision = mean(precision, na.rm = T), sd_precision = sd(precision, na.rm = T),
                          n_seeds = sum(!is.na(precision)), mean_anchors = mean(n_anchors, na.rm = T), .groups = "drop") %>%
                arrange(factor(modality, levels = modalities), fraction)
        anchor_summary = per_task %>%
                group_by(modality, fraction) %>%
                summarise(mean = mean(n_anchors, na.rm = T), sd = sd(n_anchors, na.rm = T), .groups = "drop") %>%
                arrange(factor(modality, levels = modalities), fraction)
        write_csv(cor_summary, paste0(script_output_dir, "downsample_correlation_summary.csv"))
        write_csv(prec_summary, paste0(script_output_dir, "downsample_precision_summary.csv"))
        write_csv(anchor_summary, paste0(script_output_dir, "downsample_anchor_counts.csv"))
        # precision with all samples as plotted in Fig. S2g, taken from bridge_vs_direct_matching.R (cell lines excluded)
        canonical_prec = read_tsv("./output/sample_matching_evaluation/precision_summary_with_celllines.tsv") %>%
                filter(method == "CA" | (method == "Bridge" & mark %in% modalities)) %>%
                transmute(modality = if_else(method == "CA", "CA", mark), fraction = 1,
                          mean_precision = precision_nocl, sd_precision = se_nocl, n_seeds = 3, mean_anchors = NA_real_)
        write_csv(rows_update(prec_summary, canonical_prec, by = c("modality", "fraction"), unmatched = "ignore"),
                  paste0(script_output_dir, "downsample_precision_summary_fig.csv"))
        # per-region CO-V SD across seeds
        unc_summary = map_dfr(modalities, function(mod) {
                map_dfr((1:10) / 10, function(frac) {
                        if (abs(frac - 1.0) < 1e-6) {
                                return(tibble(modality = mod, fraction = 1.0, n_regions_evaluable = NA_integer_,
                                              median_sd_cov = 0, q25_sd_cov = 0, q75_sd_cov = 0, mean_sd_cov = 0,
                                              sd_aggregate_spearman_r = 0))
                        }
                        triplet = keep(results, ~ .x$modality == mod && abs(.x$fraction - frac) < 1e-6)
                        sd_vec = matrixStats::rowSds(do.call(cbind, map(triplet, "co_v")), na.rm = T)
                        sd_valid = sd_vec[!is.na(sd_vec)]
                        qs = quantile(sd_valid, c(0.25, 0.5, 0.75), names = F)
                        tibble(modality = mod, fraction = frac, n_regions_evaluable = sum(!is.na(sd_vec)),
                               median_sd_cov = qs[2], q25_sd_cov = qs[1], q75_sd_cov = qs[3], mean_sd_cov = mean(sd_valid),
                               sd_aggregate_spearman_r = sd(per_task$spearman_r[per_task$modality == mod & abs(per_task$fraction - frac) < 1e-6], na.rm = T))
                })
        })
        write_tsv(unc_summary, paste0(script_output_dir, "cov_uncertainty_summary.tsv"))
}
