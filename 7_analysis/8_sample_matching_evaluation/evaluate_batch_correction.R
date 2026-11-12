library(tidyverse)
library(Seurat)
source("./helper/helper.R")
script_output_dir = "./output/sample_matching_evaluation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
exp_meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
tissue_map = read_tsv("./custom_metadata/tissue_system_mapping.tsv")
data_norm = readRDS("./intermediate_data/normalied_chromatin_matrices.rds")
names(data_norm) = c("human_dnase", "human_atac", "mouse_dnase", "mouse_atac")

# anchor features for batch correction, as in remove_chromatin_batch_effect.R
region_filter_human = which(matrixStats::rowMaxs(data_norm$human_dnase) > 20 |
                                    matrixStats::rowMaxs(data_norm$human_atac) > 20)
data_vars_rank = map(data_norm[c("human_dnase", "human_atac")], function(x) rank(calc_var_std_ver2(x)))
data_human_var_feature = rownames(data_norm$human_dnase)[order(data_vars_rank[[1]] + data_vars_rank[[2]], decreasing = T)[1:25000]]
data_human_var_feature = data_human_var_feature[data_human_var_feature %in% rownames(data_norm$human_dnase)[region_filter_human]]

region_filter_mouse = which(matrixStats::rowMaxs(data_norm$mouse_dnase) > 20 |
                                    matrixStats::rowMaxs(data_norm$mouse_atac) > 20)
data_vars_std_mouse_dnase = calc_var_std_ver1(data_norm$mouse_dnase)
data_vars_std_mouse_atac = calc_var_std_ver1(data_norm$mouse_atac)
data_mouse_var_feature = rownames(data_norm$mouse_dnase)[order(rank(data_vars_std_mouse_dnase) + rank(data_vars_std_mouse_atac), decreasing = T)[1:25000]]
data_mouse_var_feature = data_mouse_var_feature[data_mouse_var_feature %in% rownames(data_norm$mouse_dnase)[region_filter_mouse]]

# PCA and UMAP before and after DNase-ATAC correction
embed_before_after = function(dnase, atac, var_feature, species, k_anchor, k_filter, k_weight, umap_dims) {
        data_seurat = map2(list(dnase, atac), paste(str_to_title(species), c("Dnase", "ATAC")), function(x, source) {
                y = CreateSeuratObject(counts = x[var_feature, ], assay = "OpenChromatin")
                y = SetAssayData(y, layer = "data", new.data = GetAssayData(y, layer = "counts"))
                y$source = source
                y
        })
        merged = JoinLayers(merge(data_seurat[[1]], data_seurat[[2]]))
        merged = ScaleData(merged, features = var_feature)
        merged = RunPCA(merged, features = var_feature, npcs = 30)
        merged = RunUMAP(merged, dims = 1:30)
        anchors = FindIntegrationAnchors(data_seurat,
                                         anchor.features = var_feature,
                                         reduction = "cca",
                                         k.anchor = k_anchor,
                                         k.filter = k_filter,
                                         dims = 1:30)
        integrated = IntegrateData(anchors, features.to.integrate = anchors@anchor.features,
                                   dims = 1:15,
                                   k.weight = k_weight)
        integrated = quick_process(integrated, var_feature = var_feature, dims = umap_dims)
        map2_dfr(list(merged, integrated), c("Before", "After"), function(obj, condition) {
                bind_cols(tibble(species = species, condition = condition,
                                 Accession = colnames(obj), source = obj$source),
                          as_tibble(Embeddings(obj, "pca")[, 1:20]),
                          as_tibble(Embeddings(obj, "umap")))
        })
}
embed_tb = bind_rows(embed_before_after(data_norm$human_dnase, data_norm$human_atac, data_human_var_feature, "human",
                                        k_anchor = 8, k_filter = 20, k_weight = 10, umap_dims = 1:30),
                     embed_before_after(data_norm$mouse_dnase, data_norm$mouse_atac, data_mouse_var_feature, "mouse",
                                        k_anchor = 3, k_filter = 10, k_weight = 20, umap_dims = 1:15))
embed_tb$biosample = exp_meta$`Biosample term name`[match(embed_tb$Accession, exp_meta$Accession)]
embed_tb$tissue_system = tissue_map$tissue_system[match(paste(embed_tb$species, embed_tb$biosample),
                                                        paste(tissue_map$species, tissue_map$biosample_term))]
write_tsv(embed_tb, paste0(script_output_dir, "batch_correction_embeddings.tsv"))

# ENCODE biosample classification of tissue-labeled samples, to exclude cell lines
labeled_acc = unique(embed_tb$Accession[!is.na(embed_tb$tissue_system)])
biosample_class = map_dfr(split(labeled_acc, ceiling(seq_along(labeled_acc) / 80)), function(acc) {
        res = jsonlite::fromJSON(paste0("https://www.encodeproject.org/search/?type=Experiment&limit=all&format=json&",
                                        paste0("accession=", acc, collapse = "&"),
                                        "&field=accession&field=biosample_ontology.classification"))$`@graph`
        tibble(Accession = res$accession, classification = res$biosample_ontology$classification)
})
stopifnot(all(labeled_acc %in% biosample_class$Accession))
write_tsv(biosample_class, paste0(script_output_dir, "batch_correction_biosample_class.tsv"))
cell_lines = biosample_class$Accession[biosample_class$classification == "cell line"]

# assay offset and kNN tissue purity per tissue system, in 20 PCs
k_nn = 10
per_tissue = embed_tb %>%
        filter(!is.na(tissue_system), !Accession %in% cell_lines) %>%
        mutate(assay = ifelse(grepl("dnase", source, ignore.case = T), "DNase", "ATAC")) %>%
        group_by(species, condition) %>%
        group_modify(function(d, key) {
                pc_mat = as.matrix(select(d, PC_1:PC_20))
                cent = t(sapply(sort(unique(d$tissue_system)), function(x) colMeans(pc_mat[d$tissue_system == x, , drop = F])))
                scale_med = median(as.numeric(dist(cent)))
                nn = FNN::get.knn(pc_mat, k = k_nn)$nn.index
                knn_purity = rowMeans(matrix(d$tissue_system[nn], ncol = k_nn) == d$tissue_system)
                keep = d %>% count(tissue_system, assay) %>%
                        pivot_wider(names_from = assay, values_from = n, values_fill = 0) %>%
                        filter(DNase >= 2, ATAC >= 2) %>%
                        pull(tissue_system)
                map_dfr(keep, function(x) {
                        dnase_pc = pc_mat[d$tissue_system == x & d$assay == "DNase", , drop = F]
                        atac_pc = pc_mat[d$tissue_system == x & d$assay == "ATAC", , drop = F]
                        tibble(tissue_system = x, n_dnase = nrow(dnase_pc), n_atac = nrow(atac_pc),
                               offset_raw = sqrt(sum((colMeans(dnase_pc) - colMeans(atac_pc))^2)),
                               offset_norm = offset_raw / scale_med,
                               purity = mean(knn_purity[d$tissue_system == x]),
                               purity_all = mean(knn_purity))
                })
        }) %>%
        ungroup()
write_tsv(select(per_tissue, -purity_all), paste0(script_output_dir, "batch_correction_per_tissue.tsv"))
batch_summary = per_tissue %>%
        select(species, tissue_system, condition, offset_norm, purity, purity_all) %>%
        pivot_wider(names_from = condition, values_from = c(offset_norm, purity, purity_all)) %>%
        group_by(species) %>%
        summarize(n_tissues = n(),
                  offset_before_mean = mean(offset_norm_Before),
                  offset_after_mean = mean(offset_norm_After),
                  offset_n_shrink = sum(offset_norm_After < offset_norm_Before),
                  offset_p = wilcox.test(offset_norm_Before, offset_norm_After, paired = T, alternative = "greater")$p.value,
                  purity_before_mean = mean(purity_Before),
                  purity_after_mean = mean(purity_After),
                  purity_n_up = sum(purity_After > purity_Before),
                  purity_p = wilcox.test(purity_Before, purity_After, paired = T, exact = F)$p.value,
                  purity_all_before = purity_all_Before[1],
                  purity_all_after = purity_all_After[1])
write_tsv(batch_summary, paste0(script_output_dir, "batch_correction_summary.tsv"))
