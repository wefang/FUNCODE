library(tidyverse)
library(Seurat)
source("./helper/helper.R")
script_output_dir = "./output/sample_matching_evaluation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
exp_meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
exp_meta$`Tissue/cell types` = NA
tissue_map = read_tsv("./custom_metadata/tissue_system_mapping.tsv")
cell_line_terms = readxl::read_excel("./metadata/table_s2_experiment_metadata_annotated_final.xlsx") %>%
        filter(Biosample.type == "cell line") %>%
        distinct(Biosample.term.name) %>%
        pull(Biosample.term.name)

# tissue system of a sample, embryonic liver and limb reassigned
assign_system = function(biosample, life_stage, sp) {
        species_map = tissue_map[tissue_map$species == sp, ]
        system_base = species_map$tissue_system[match(biosample, species_map$biosample_term)]
        case_when(system_base == "Liver" & life_stage == "embryonic" ~ "Immune/Blood",
                  system_base == "Limb/Facial" & life_stage == "embryonic" ~ "Muscle",
                  TRUE ~ system_base)
}
classify_anchors = function(anchor_df) {
        acc_meta = function(acc, col) exp_meta[[col]][match(acc, exp_meta$Accession)]
        tibble(human_accession = anchor_df$accession1,
               mouse_accession = anchor_df$accession2,
               human_biosample = acc_meta(human_accession, "Biosample term name"),
               mouse_biosample = acc_meta(mouse_accession, "Biosample term name"),
               human_assay = acc_meta(human_accession, "Assay name"),
               mouse_assay = acc_meta(mouse_accession, "Assay name"),
               human_target = acc_meta(human_accession, "Target of assay"),
               mouse_target = acc_meta(mouse_accession, "Target of assay"),
               human_life_stage = acc_meta(human_accession, "Life stage"),
               mouse_life_stage = acc_meta(mouse_accession, "Life stage"),
               human_tissue_system = assign_system(human_biosample, human_life_stage, "human"),
               mouse_tissue_system = assign_system(mouse_biosample, mouse_life_stage, "mouse"),
               either_cell_line = human_biosample %in% cell_line_terms | mouse_biosample %in% cell_line_terms,
               tissue_match = ifelse(is.na(human_tissue_system) | is.na(mouse_tissue_system), NA,
                                     human_tissue_system == mouse_tissue_system),
               anchor_score = anchor_df$score)
}

#### begin published chromatin accessibility anchors ####
anchors_all = read_tsv("./output/ENCODE4_conservation/anchors.tsv")
anchor_tb = classify_anchors(anchors_all) %>% filter(!is.na(tissue_match))
p_human = table(anchor_tb$human_tissue_system) / nrow(anchor_tb)
p_mouse = table(anchor_tb$mouse_tissue_system) / nrow(anchor_tb)
shared_systems = intersect(names(p_human), names(p_mouse))
p_random = sum(p_human[shared_systems] * p_mouse[shared_systems])
n_match = sum(anchor_tb$tissue_match)
expected_matches = nrow(anchor_tb) * p_random
published_summary = tibble(n_pairs = nrow(anchors_all),
                           n_valid = nrow(anchor_tb),
                           n_match = n_match,
                           precision = n_match / nrow(anchor_tb),
                           weighted_precision = sum(anchor_tb$tissue_match * anchor_tb$anchor_score) / sum(anchor_tb$anchor_score),
                           p_random = p_random,
                           expected_matches = expected_matches,
                           odds_ratio = (n_match / (nrow(anchor_tb) - n_match)) / (expected_matches / (nrow(anchor_tb) - expected_matches)),
                           chisq_pvalue = chisq.test(c(n_match, nrow(anchor_tb) - n_match), p = c(p_random, 1 - p_random))$p.value)
write_tsv(published_summary, paste0(script_output_dir, "published_anchors_precision.tsv"))

#### begin re-integration with tuned parameters ####
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
histone_meta = filter(exp_meta, `Target of assay` %in% c("H3K4me1", "H3K4me3", "H3K27ac"), `Assay name` != "Mint-ChIP-seq")
best_params = list(H3K27ac = list(k_anchor = 20, k_filter = 20, anchor_features = 20000),
                   H3K4me3 = list(k_anchor = 20, k_filter = 50, anchor_features = 10000),
                   H3K4me1 = list(k_anchor = 5, k_filter = 50, anchor_features = 20000))
# add_metadata of helper.R resets cell names under Seurat v5, metadata is matched instead
add_sample_meta = function(obj, meta) {
        obj$Accession = colnames(obj)
        for (col in c("Biosample term name", "Target of assay")) {
                obj[[col]] = meta[[col]][match(colnames(obj), meta$Accession)]
        }
        obj
}
find_anchors = function(human_mat, mouse_mat, meta, mark, nfeatures, k_anchor, k_filter, anchor_features) {
        human_obj = add_sample_meta(CreateSeuratObject(counts = human_mat, assay = "OpenChromatin"), meta)
        mouse_obj = add_sample_meta(CreateSeuratObject(counts = mouse_mat, assay = "OpenChromatin"), meta)
        # removing over-represented samples for human H3K27ac
        if (mark == "H3K27ac") {
                human_obj = human_obj[, !human_obj$`Biosample term name` %in% c("skin epidermis", "middle frontal area 46")]
        }
        human_obj = FindVariableFeatures(NormalizeData(human_obj), nfeatures = nfeatures)
        mouse_obj = FindVariableFeatures(NormalizeData(mouse_obj), nfeatures = nfeatures)
        anchors = FindIntegrationAnchors(list(human = human_obj, mouse = mouse_obj),
                                         anchor.features = anchor_features,
                                         k.anchor = k_anchor,
                                         k.filter = k_filter,
                                         k.score = 30,
                                         dims = 1:30)
        make_anchor_tb(anchors, human_obj, mouse_obj)
}

# chromatin accessibility
anchor_df = find_anchors(human_correct[data_cross_var_feature, ], mouse_correct[data_cross_var_feature, ], exp_meta, "CA",
                         nfeatures = 40000, k_anchor = 10, k_filter = 20, anchor_features = 20000)
classified_all = classify_anchors(anchor_df) %>% mutate(method = "CA", mark = NA_character_)

# bridge integration through chromatin accessibility, histone-histone anchors of the mark
for (mark in names(best_params)) {
        human_im_mat = readRDS(paste0("./intermediate_data/histone_rev_transfer_v1/", "human_rev_transfer_", mark, ".rds"))
        mouse_im_mat = readRDS(paste0("./intermediate_data/histone_rev_transfer_v1/", "mouse_rev_transfer_", mark, ".rds"))
        anchor_df = exec(find_anchors, cbind(human_correct[data_cross_var_feature, ], human_im_mat),
                         cbind(mouse_correct[data_cross_var_feature, ], mouse_im_mat), histone_meta, mark,
                         nfeatures = 100000, !!!best_params[[mark]])
        anchor_df = anchor_df[anchor_df$`Target of assay_1` %in% mark & anchor_df$`Target of assay_2` %in% mark, ]
        classified_all = bind_rows(classified_all, classify_anchors(anchor_df) %>% mutate(method = "Bridge", mark = mark))
}
region_names = rownames(human_correct)
rm(human_correct, mouse_correct)
gc()

# direct integration of each mark
for (mark in names(best_params)) {
        human_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Homo_sapiens_", mark, "_norm_correct.rds"))
        mouse_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Mus_musculus_", mark, "_norm_correct.rds"))
        if (is.null(rownames(human_histone)) || all(rownames(human_histone) == "")) {
                rownames(human_histone) = region_names[1:nrow(human_histone)]
                rownames(mouse_histone) = region_names[1:nrow(mouse_histone)]
        }
        shared_features = intersect(intersect(data_cross_var_feature, rownames(human_histone)), rownames(mouse_histone))
        mark_meta = filter(exp_meta, `Target of assay` == mark, `Assay name` != "Mint-ChIP-seq")
        anchor_df = exec(find_anchors, human_histone[shared_features, ], mouse_histone[shared_features, ], mark_meta, mark,
                         nfeatures = 100000, !!!best_params[[mark]])
        classified_all = bind_rows(classified_all, classify_anchors(anchor_df) %>% mutate(method = "Direct", mark = mark))
}
write_tsv(classified_all, paste0(script_output_dir, "anchor_pair_classifications.tsv"))

# precision with and without cell lines, bootstrap SE
set.seed(42)
n_boot = 1000
precision_summary = classified_all %>%
        group_by(method, mark) %>%
        summarise(n_total_anchors = n(),
                  n_match_all = sum(tissue_match, na.rm = T),
                  precision_all = n_match_all / n_total_anchors,
                  n_nocl_anchors = sum(!either_cell_line),
                  n_match_nocl = sum(tissue_match[!either_cell_line], na.rm = T),
                  precision_nocl = n_match_nocl / n_nocl_anchors,
                  n_cl_anchors = sum(either_cell_line),
                  n_match_cl = sum(tissue_match[either_cell_line], na.rm = T),
                  precision_cl = ifelse(n_cl_anchors > 0, n_match_cl / n_cl_anchors, NA_real_),
                  .groups = "drop") %>%
        left_join(classified_all %>%
                          group_by(method, mark) %>%
                          summarise(se_all = sd(replicate(n_boot, mean(sample(tissue_match, replace = T), na.rm = T))),
                                    se_nocl = sd(replicate(n_boot, mean(sample(tissue_match[!either_cell_line], replace = T), na.rm = T))),
                                    .groups = "drop"),
                  by = c("method", "mark"))
write_tsv(precision_summary, paste0(script_output_dir, "precision_summary_with_celllines.tsv"))
