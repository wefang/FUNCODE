# Rscript feature_noise_robustness.R prep, then run <none|unselected|low_var|gene_hom|random|shuffled> <fraction> <seed> per task, then summarize
library(tidyverse)
library(data.table)
options(Seurat.object.assay.version = "v3")
library(Seurat)
source("./helper/helper.R")
script_output_dir = "./output/feature_noise_robustness/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
args = commandArgs(trailingOnly = T)
compute_cov = function(human_mat, mouse_mat, anchor_df) {
        sapply(1:nrow(human_mat), function(j) {
                x = human_mat[j, anchor_df$accession1]
                y = mouse_mat[j, anchor_df$accession2]
                if (var(x) > 0 & var(y) > 0) {
                        wCorr::weightedCorr(x, y, weights = anchor_df$score, method = "spearman")
                } else {
                        return(NA_real_)
                }
        })
}

#### begin replacement feature pools ####
if (args[1] == "prep") {
        regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
        human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
        mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
        rownames(human_correct) = rownames(mouse_correct) = regions$id
        data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
        # the 10,000 matching features selected as in integrate_chromatin.R
        human_obj = FindVariableFeatures(CreateSeuratObject(counts = human_correct[data_cross_var_feature, ], assay = "OpenChromatin"),
                                         nfeatures = 40000)
        mouse_obj = FindVariableFeatures(CreateSeuratObject(counts = mouse_correct[data_cross_var_feature, ], assay = "OpenChromatin"),
                                         nfeatures = 40000)
        canonical_features = SelectIntegrationFeatures(list(human_obj, mouse_obj), nfeatures = 10000)
        saveRDS(canonical_features, paste0(script_output_dir, "canonical_10k_features.rds"))
        rm(human_obj, mouse_obj)
        # least variable aligned pairs above the signal floor of integrate_chromatin.R
        human_var_std = calc_var_std_ver1(human_correct)
        mouse_var_std = calc_var_std_ver1(mouse_correct)
        region_filter = matrixStats::rowMaxs(human_correct) > 20 & matrixStats::rowMaxs(mouse_correct) > 20
        low_var_pool = tibble(id = regions$id, comb_rank = (rank(human_var_std) + rank(mouse_var_std)) / 2)[region_filter, ] %>%
                arrange(comb_rank) %>%
                slice_head(n = 100000) %>%
                pull(id)
        saveRDS(low_var_pool, paste0(script_output_dir, "low_var_pool.rds"))
        rm(human_correct, mouse_correct)
        gc()

        # gene-homology mispairs, human element of a pair with the mouse element of another pair near the same homologous genes
        regions_homolog = readRDS("./output/indexed_dhs_mapped_regions_homolog_genes.rds")
        ids = as.character(regions_homolog$id)
        keep = lengths(regions_homolog$homolog_gene_filter_tb) > 0
        flat = rbindlist(regions_homolog$homolog_gene_filter_tb[keep], idcol = "row", fill = T)
        flat[, region_id := ids[keep][row]]
        setorder(flat, row, dist_to_tss_rank)
        best = flat[!is.na(hid) & hid != "", .SD[1], by = row][, .(region_id, hid, gene_human, gene_mouse)]
        setorder(best, hid, region_id)
        best[, n_in_group := .N, by = hid]
        grp = best[n_in_group >= 2]
        grp[, idx := seq_len(.N), by = hid]
        set.seed(42)
        grp[, pick := ((idx - 1L + sample.int(.N - 1L, .N, replace = T)) %% .N) + 1L, by = hid]
        grp[, mouse_id := region_id[pick], by = hid]
        gh_mispairs = grp[region_id != mouse_id, .(human_id = region_id, mouse_id = mouse_id, hid = hid,
                                                    gene_human = gene_human, gene_mouse = gene_mouse)]
        if (nrow(gh_mispairs) > 200000L) gh_mispairs = gh_mispairs[sample(.N, 200000L)]
        saveRDS(as.data.frame(gh_mispairs), paste0(script_output_dir, "core_gene_homology_mispair_index.rds"))

        # 100,000 extended set element pairs, sampled uniformly over unique DHS pairs
        load("./output/all_dhs_pairs.rda")
        set.seed(20260825)
        h_ids = lapply(human_hid_set$human_dhs, `[[`, "id")
        m_ids = lapply(mouse_hid_set$mouse_dhs, `[[`, "id")
        hi = match(hidset_pairs$human_hid_set, human_hid_set$hid_set)
        mi = match(hidset_pairs$mouse_hid_set, mouse_hid_set$hid_set)
        set_size = as.numeric(lengths(h_ids)[hi]) * lengths(m_ids)[mi]
        h_sets = split(rep(human_hid_set$hid_set, lengths(h_ids)), unlist(h_ids))
        m_sets = split(rep(mouse_hid_set$hid_set, lengths(m_ids)), unlist(m_ids))
        pair_keys = paste(hidset_pairs$human_hid_set, hidset_pairs$mouse_hid_set)
        multiplicity = function(h, m) {
                sum(outer(h_sets[[h]], m_sets[[m]], paste) %in% pair_keys)
        }
        sampled = character(0)
        while (length(sampled) < 100000L) {
                draw = sample.int(length(set_size), 2L * 100000L, replace = T, prob = set_size)
                h = vapply(h_ids[hi[draw]], function(v) v[sample.int(length(v), 1L)], "")
                m = vapply(m_ids[mi[draw]], function(v) v[sample.int(length(v), 1L)], "")
                k = mapply(multiplicity, h, m)
                accept = runif(length(k)) < 1 / k
                sampled = unique(c(sampled, paste(h[accept], m[accept], sep = "\t")))
        }
        sampled = str_split_fixed(sampled[1:100000L], "\t", 2)
        ext_pairs = tibble(human_id = sampled[, 1], mouse_id = sampled[, 2])
        write_tsv(ext_pairs, paste0(script_output_dir, "extended_pairs_sample_100k.tsv"))
        human_all = readRDS("./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds")
        human_ext = as.matrix(human_all[intersect(unique(ext_pairs$human_id), rownames(human_all)), , drop = F])
        rm(human_all)
        gc()
        mouse_all = readRDS("./intermediate_data/encodev4_all_dhs_mouse_chromatin_correct.rds")
        mouse_ext = as.matrix(mouse_all[intersect(unique(ext_pairs$mouse_id), rownames(mouse_all)), , drop = F])
        saveRDS(list(pairs = ext_pairs, human = human_ext, mouse = mouse_ext), paste0(script_output_dir, "extended_sample_matrices.rds"))
}

#### begin integration with replaced matching features ####
if (args[1] == "run") {
        noise_type = args[2]
        fraction = as.numeric(args[3])
        seed = as.integer(args[4])
        regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
        human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
        mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
        rownames(human_correct) = rownames(mouse_correct) = regions$id
        data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
        exp_meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
        tissue_map = read_tsv("./custom_metadata/tissue_system_mapping_v3.tsv")
        set.seed(seed)
        n_noise = as.integer(round(10000L * fraction))
        n_real = 10000L - n_noise
        canonical_features = readRDS(paste0(script_output_dir, "canonical_10k_features.rds"))
        feat_human = feat_mouse = sample(canonical_features, n_real)
        if (n_noise > 0) {
                if (noise_type == "unselected") {
                        pool = setdiff(regions$id, data_cross_var_feature)
                        pool = pool[matrixStats::rowMaxs(human_correct[pool, , drop = F]) > 20 &
                                            matrixStats::rowMaxs(mouse_correct[pool, , drop = F]) > 20]
                        feat_human = feat_mouse = c(feat_human, sample(pool, n_noise))
                } else if (noise_type == "low_var") {
                        pool = readRDS(paste0(script_output_dir, "low_var_pool.rds"))
                        feat_human = feat_mouse = c(feat_human, sample(pool, n_noise))
                } else if (noise_type == "gene_hom") {
                        gh_mispairs = readRDS(paste0(script_output_dir, "core_gene_homology_mispair_index.rds"))
                        gh_mispairs = gh_mispairs[gh_mispairs$human_id %in% rownames(human_correct) &
                                                          gh_mispairs$mouse_id %in% rownames(mouse_correct), ]
                        idx = sample(nrow(gh_mispairs), n_noise)
                        feat_human = c(feat_human, gh_mispairs$human_id[idx])
                        feat_mouse = c(feat_mouse, gh_mispairs$mouse_id[idx])
                } else if (noise_type == "random") {
                        feat_human = c(feat_human, sample(regions$id, n_noise))
                        feat_mouse = c(feat_mouse, sample(regions$id, n_noise))
                } else if (noise_type == "shuffled") {
                        feat_human = feat_mouse = c(feat_human, sample(regions$id, n_noise))
                }
        }
        # features are paired by row name, so both matrices get shared synthetic names
        synthetic = sprintf("f%06d", 1:10000)
        human_mat = human_correct[feat_human, , drop = F]
        mouse_mat = mouse_correct[feat_mouse, , drop = F]
        rownames(human_mat) = rownames(mouse_mat) = synthetic
        if (noise_type == "shuffled" && n_noise > 0) {
                for (i in (n_real + 1L):10000L) {
                        human_mat[i, ] = human_mat[i, sample.int(ncol(human_mat))]
                        mouse_mat[i, ] = mouse_mat[i, sample.int(ncol(mouse_mat))]
                }
        }
        human_obj = CreateSeuratObject(counts = human_mat, assay = "OpenChromatin")
        mouse_obj = CreateSeuratObject(counts = mouse_mat, assay = "OpenChromatin")
        human_obj$Accession = colnames(human_obj)
        mouse_obj$Accession = colnames(mouse_obj)
        VariableFeatures(human_obj) = rownames(human_obj)
        VariableFeatures(mouse_obj) = rownames(mouse_obj)
        anchors = FindIntegrationAnchors(list(human = human_obj, mouse = mouse_obj),
                                         anchor.features = synthetic,
                                         reduction = "cca",
                                         k.anchor = 5,
                                         k.filter = 20,
                                         k.score = 30,
                                         dims = 1:30)
        anchor_df = make_anchor_tb(anchors, human_obj, mouse_obj)
        # CO-V on a fixed subset of 100,000 core regions
        set.seed(20260825)
        score_idx = sort(sample(nrow(human_correct), 100000L))
        co_v = compute_cov(human_correct[score_idx, ], mouse_correct[score_idx, ], anchor_df)
        # tissue-system precision, embryonic liver and limb reassigned
        assign_system = function(acc, sp) {
                species_map = tissue_map[tissue_map$species == sp, ]
                system_base = species_map$tissue_system[match(exp_meta$`Biosample term name`[match(acc, exp_meta$Accession)],
                                                              species_map$biosample_term)]
                life_stage = exp_meta$`Life stage`[match(acc, exp_meta$Accession)]
                case_when(system_base == "Liver" & life_stage == "embryonic" ~ "Immune/Blood",
                          system_base == "Limb/Facial" & life_stage == "embryonic" ~ "Muscle",
                          TRUE ~ system_base)
        }
        anchor_tb = anchor_df %>%
                mutate(system_human = assign_system(accession1, "human"),
                       system_mouse = assign_system(accession2, "mouse")) %>%
                filter(!is.na(system_human), !is.na(system_mouse))
        saveRDS(list(noise_type = noise_type, fraction = fraction, seed = seed,
                     n_real = n_real, n_noise = n_noise,
                     n_anchors = nrow(anchor_df), n_evaluable = nrow(anchor_tb),
                     precision = mean(anchor_tb$system_human == anchor_tb$system_mouse),
                     w_precision = sum((anchor_tb$system_human == anchor_tb$system_mouse) * anchor_tb$score) / sum(anchor_tb$score),
                     anchor_df = anchor_df, score_ids = rownames(human_correct)[score_idx], co_v = co_v),
                paste0(script_output_dir, sprintf("match_%s_frac%.2f_seed%d.rds", noise_type, fraction, seed)))
}

#### begin summary over runs ####
if (args[1] == "summarize") {
        res = map(list.files(script_output_dir, pattern = "^match_.*\\.rds$", full.names = T), readRDS)
        ref = Filter(function(r) r$fraction == 0, res)[[1]]
        noise_summary = map_dfr(res, function(r) {
                ok = !is.na(r$co_v) & !is.na(ref$co_v)
                tibble(noise_type = r$noise_type, fraction = r$fraction, seed = r$seed,
                       n_anchors = r$n_anchors, precision = r$precision, w_precision = r$w_precision,
                       n_evaluable = r$n_evaluable,
                       cov_spearman = if (sum(ok) > 100) cor(r$co_v[ok], ref$co_v[ok], method = "spearman") else NA_real_,
                       n_scored = sum(ok))
        })
        # the undiluted run starts every replacement class
        noise_summary = bind_rows(noise_summary %>% filter(fraction > 0),
                                  crossing(noise_type = c("low_var", "unselected", "gene_hom", "random", "shuffled"),
                                           noise_summary %>% filter(fraction == 0) %>% select(-noise_type))) %>%
                arrange(noise_type, fraction, seed)
        write_tsv(noise_summary, paste0(script_output_dir, "feature_noise_summary.tsv"))

        # GH-CO-V on the extended set pairs from the anchors of each run
        ext = readRDS(paste0(script_output_dir, "extended_sample_matrices.rds"))
        ext_pairs = ext$pairs[ext$pairs$human_id %in% rownames(ext$human) & ext$pairs$mouse_id %in% rownames(ext$mouse), ]
        human_ext = ext$human[ext_pairs$human_id, , drop = F]
        mouse_ext = ext$mouse[ext_pairs$mouse_id, , drop = F]
        ext_cov = map(res, function(r) compute_cov(human_ext, mouse_ext, r$anchor_df))
        ext_ref = ext_cov[[which(map_lgl(res, function(r) r$noise_type == "none" && r$fraction == 0))]]
        ext_concordance = map2_dfr(res, ext_cov, function(r, co_v) {
                ok = !is.na(co_v) & !is.na(ext_ref)
                tibble(noise_type = r$noise_type, fraction = r$fraction, seed = r$seed, n_scored = sum(ok),
                       ghcov_spearman = if (sum(ok) > 100) cor(co_v[ok], ext_ref[ok], method = "spearman") else NA_real_)
        }) %>%
                arrange(noise_type, fraction)
        write_tsv(ext_concordance, paste0(script_output_dir, "extended_cov_concordance.tsv"))
}
