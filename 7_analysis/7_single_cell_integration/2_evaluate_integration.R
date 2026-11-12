# Evaluations on single cell integration
library(Seurat)
library(tidyverse)
source("./scatac_helpers.R")
source("./helper/helper.R")
source("./helper/def_color.R")

script_plot_dir =  "./plots/single_cell/ren_sciatac/"
script_output_dir = "./output/single_cell_ren/"

# ndim = 15
ari_collect = list()
ndim = 15
kval = 40
num_features = 48160
# num_features = 25000
for (method in c("method1", "method2", "method3", "method4_phyloP4way", "method5_spearman_manual",
                 "method6_lecif", "method7_percentage", "method8_human_var", "method9_mouse_var",
                 "method10_adaliftover_h2m", "method11_adaliftover_m2h")) {
        if (method %in% c("method1", "method2", "method3", "method6_lecif", "method7_percentage")) {
                data_dir = "E:\\large_data/sciatac_seurat_ren/"
        }
        if (method %in% c("method4_phyloP4way", "method5_spearman_manual")) {
                data_dir = "E:\\GlobusDownload/Ruzhang_sciatac_integrated_Feb7/"
        }
        if (method %in% c("method8_human_var", "method9_mouse_var", "method10_adaliftover_h2m", "method11_adaliftover_m2h")) {
                data_dir = "/dcs04/hongkai/data/rzhao/sciATAC/mid_data/"
        }
        if (num_features == 25000) {
                int_obj = readRDS(paste0(data_dir, "integrated_human_mouse_example1_", method, ".rds"))
        }
        if (num_features == 48160) {
                int_obj = readRDS(paste0(data_dir, "integrated_human_mouse_example1_", method, "_48.rds"))
        }
        int_obj = process_scatac(int_obj, features = rownames(int_obj),
                                 dims = 2:ndim,
                                 run_umap = T, normalize = F)
        all_cell_types = sort(unique(int_obj$cell_type))
        matched_label_col_all = matched_label_col[all_cell_types]
        g_umap_spearman = UMAPPlot(int_obj, group.by = "cell_type", split.by = "species", pt.size = 0.5) +
                scale_color_manual(values = matched_label_col_all) +
                theme(text = element_text(size = 8),
                      # legend.position = "none",
                      legend.text = element_text(size = 10), legend.position = "none")
        g_umap_spearman
        push_png(g_umap_spearman, paste0("umap_example1_", method, "_nfeature", num_features, "_ndim", ndim , "_k", kval), w = 4, h = 2.5)
        # HoverLocator(g_umap_spearman, information = FetchData(int_obj, "cell_type"))
        
        u_embed = int_obj@reductions$lsi@cell.embeddings
        dataset = int_obj$species
        data1 <- u_embed[which(dataset == "human"), 2:ndim]
        data2 <- u_embed[which(dataset == "mouse"), 2:ndim]
        k1 <- kval
        k2 <- kval
        n1 <- nrow(data1)
        n2 <- nrow(data2)
        k1<-min(k1,n1)
        k2<-min(k2,n2)
        n.total <- n1 + n2
        
        W21 <- FNN::get.knnx(data2, query=data1, k=k1)
        js1 <- matrix(seq_len(n1), n1, k1)
        is1 <- n1 + W21$nn.index
        indices1 <- cbind(as.vector(js1), as.vector(is1))
        
        W12 <- FNN::get.knnx(data1, query=data2, k=k2)
        js2 <- matrix(n1 + seq_len(n2), n2, k2)
        is2 <- W12$nn.index
        indices2 <- cbind(as.vector(js2), as.vector(is2))
        
        W = Matrix::sparseMatrix(i = c(indices1[, 1], indices2[, 1]),
                                 j = c(indices1[, 2], indices2[, 2]),
                                 x = rep(1, nrow(indices1) + nrow(indices2)),
                                 dims = c(n.total, n.total))
        W <- W * Matrix::t(W) # elementwise multiplication to keep mutual nns only
        Wsp = as(W, "dgTMatrix")
        
        mnn_indices = cbind(Wsp@i, Wsp@j)
        mnn_indices = mnn_indices[mnn_indices[, 1] < n1, ]
        mnn_indices = mnn_indices+1
        mnn_indices[, 2] = mnn_indices[, 2] - n1
        A1 = mnn_indices[, 1]
        A2 = mnn_indices[, 2]
        
        # calculate weight vector for label transfer
        sigma<-20
        kk <- min(length(A2),100)
        W <- FNN::get.knnx(data2[A2, ], query=data2, k=kk)
        
        G = Matrix::sparseMatrix(i = c(row(W$nn.index)), j = c(W$nn.index), x = exp(-(c(W$nn.dist))^2/sigma))
        G_norm <- Matrix::Diagonal(x = 1 / Matrix::rowSums(G)) %*% G
        
        L_A = Matrix::sparse.model.matrix(object = ~ cell_type - 1,
                                          data = tibble(cell_type = factor(int_obj$cell_type[int_obj$species == "human"][A1])))
        L_P = G_norm %*% L_A
        label_pred = sort(unique(int_obj$cell_type[int_obj$species == "human"][A1]))[apply(as.matrix(L_P), 1, which.max)]
        saveRDS(label_pred, file = output_rds_path(paste0("label_pred_", method, "_nfeature", num_features, "_ndim", ndim , "_k", kval)))
        label_true = int_obj$cell_type[int_obj$species == "mouse"]

        ari_val = mclust::adjustedRandIndex(matched_label_mapper[label_true],
                                    matched_label_mapper[label_pred])
        print(ari_val)
        ari_tb = tibble(method = method,
                        num_features = num_features,
                        ari = ari_val,
                        ndim = ndim,
                        kval = kval)
        ari_collect = append(ari_collect, list(ari_tb))

        # visualize transferred labels
        # int_obj$cell_type_pred = rlang::duplicate(int_obj$cell_type)
        # int_obj$cell_type_pred[int_obj$species == "mouse"] = label_pred
        # g_pred = UMAPPlot(int_obj, group.by = "cell_type_pred", split.by = "species", pt.size = 0.5) +
        #         scale_color_manual(values = matched_label_col_all) +
        #         theme(# legend.position = "none",
        #                 legend.text = element_text(size = 10),
        #                 legend.position = "none", text = element_text(size = 8))
        # push_png(g_pred, paste0("umap_example1_", method, "_nfeature", num_features, "_ndim", ndim , "_k", kval, "_pred"), w = 4, h = 2.5)
}
ari_tb = bind_rows(ari_collect)

method_map = c("method1" = "CACO-V",
               "method5_spearman_manual" = "CACO-V_Manual",
               "method4_phyloP4way" = "PhyloP4Way",
               "method3" = "PhastCons4Way",
               "method2" = "Shared",
               "method6_lecif" = "LECIF",
               "method7_percentage" = "BasePercent")

(ari_tb %>% mutate(method = method_map[method]) %>%
        ggline(x = "num_features",y  ="ari", col = "method",
                  ylab = "Adjusted Rand Index (ARI)", xlab = "Number of features") +
        scale_color_manual(values = method_col) +
        theme(legend.position = "right", text = element_text(size = 12))) %>%
        push_pdf("air_eval", width = 4., h = 3., ps = 12)

#### Evaluation of calling conserved differential signals ####
ndim = 20
kval = 40
num_features = 48160
celltype1 =  "Hepatocytes"
celltype2 = "Endothelial cells"
# int_obj is loaded for label information
data_dir = "E:\\large_data/sciatac_seurat_ren/"
method = "method1"
int_obj = readRDS(paste0(data_dir, "integrated_human_mouse_example1_", method, "_48.rds"))
int_obj$matched_cell_type = matched_label_mapper[int_obj$cell_type]
UMAPPlot(int_obj, group.by = "cell_type") + scale_color_manual(values = matched_label_col)
features_eval = dhs_common_all

res_list = map(c("truth", names(method_map)), function(x) {
        print(x)
        # if truth, use true mouse labels, otherwise use transferred labels
        if (x == "truth") {
                x_label = label_true
        } else {
                x_label = readRDS(output_rds_path(paste0("label_pred_", x, "_nfeature", num_features, "_ndim", ndim , "_k", kval)))
        }
        transfer_tb = tibble(cell = colnames(int_obj)[int_obj$species == "mouse"],
                             label = matched_label_mapper[x_label],
                             label_true = matched_label_mapper[label_true])

        cells = transfer_tb$cell[transfer_tb$label == celltype1]
        mouse_psb1 = Matrix::rowSums(mouse_all@assays$DHS@counts[features_eval, cells, drop = F])
        mouse_psb1 = log2(1 + mouse_psb1 / sum(mouse_psb1) * 1e6)

        cells = transfer_tb$cell[transfer_tb$label == celltype2]
        mouse_psb2 = Matrix::rowSums(mouse_all@assays$DHS@counts[features_eval, cells, drop = F])
        mouse_psb2 = log2(1 + mouse_psb2 / sum(mouse_psb2) * 1e6)
        # report differential signal between cell types
        list(mouse_diff = mouse_psb1 - mouse_psb2)
})
# compute fraction of differential signals recovered for different ranks of human and mouse differential signals
diff_cutoff = 3
human_diff = human_psb1 - human_psb2
# true mouse difference
true_diff = res_list[[1]][[1]]
true_indices = which((abs(true_diff) > diff_cutoff) & (abs(human_diff) > diff_cutoff) & ((human_diff * true_diff) > 0))
length(true_indices)
res_tb = bind_rows(map(seq(1500, 3000, by = 250), function(n_rank) {
        bind_rows(map(2:8, function(i) {
                obs_diff = res_list[[i]][[1]]
                obs_direct = as.numeric(obs_diff * human_diff)
                # obs_diff_bg = abs(obs_diff) * bg_set
                # test_indices = order(obs_diff_bg, decreasing = T)[1:n_rank]
                diff_rank = rank(abs(obs_diff * obs_direct)) + rank(abs(human_diff * obs_direct))
                test_indices = order(diff_rank, decreasing = T)[1:n_rank]
                # test_indices = which((abs(obs_diff) > diff_cutoff) & (abs(human_diff) > diff_cutoff) & ((human_diff * obs_diff) > 0))
                tibble(method = method_map[names(method_map)[i - 1]],
                       n_rank = n_rank,
                       recall = length(intersect(true_indices, test_indices)) / n_rank,
                       precision = length(intersect(true_indices, test_indices)) / length(true_indices),
                )
        }))
}))
# Figure 6f:
(res_tb %>% ggline(x = "n_rank", y = "precision", color = "method", numeric.x.axis = T,
                  ylab = "Precision", xlab = "DHS Rank", size = 0.25) +
        scale_color_manual(values = method_col) +
        theme(text = element_text(size = 10),
              legend.position = "none")) %>%
        push_pdf("example_precision_rank", w = 2, h = 2.75)

# Figure 6e: PIE CHART 
pie_list = map(names(method_map), function(x) {
        print(x)
        if (x == "truth") {
                x_label = label_true
        } else {
                x_label = readRDS(output_rds_path(paste0("label_pred_", x, "_nfeature", num_features, "_ndim", ndim , "_k", kval)))
        }
        transfer_tb = tibble(cell = colnames(int_obj)[int_obj$species == "mouse"],
                             label = matched_label_mapper[x_label],
                             label_true = matched_label_mapper[label_true])
        g1 = transfer_tb %>%
                group_by(label) %>%
                count(label_true) %>%
                filter(label == celltype1) %>%
                ggpie(x = "n", label = "label_true", fill = "label_true") +
                scale_fill_manual(values = matched_col[1:8]) +
                theme(legend.position = "none")
        g2 = transfer_tb %>%
                group_by(label) %>%
                count(label_true) %>%
                filter(label == celltype2) %>%
                ggpie(x = "n", label = "label_true", fill = "label_true") +
                scale_fill_manual(values = matched_col[1:8]) +
                theme(legend.position = "none")
        g1 + g2
})
purrr::reduce(pie_list, `/`) %>% push_pdf("pie_chart_example", height = 10)

