#' Add experiment metadata to Seurat object
#' @export
add_metadata <- function(data_obj, exp_tb) {
        data_obj$Accession = colnames(data_obj)
        data_obj@meta.data = left_join(data_obj@meta.data,
                                       dplyr::select(exp_tb,
                                              Accession,
                                              `Biosample term name`,
                                              `Biosample treatment`,
                                              `Biosample age`,
                                              `Life stage`,
                                              Lab,
                                              `Assay title`,
                                              `Assay title`,
                                              `Target of assay`,
                                              `Tissue/cell types`))
        rownames(data_obj@meta.data) = colnames(data_obj)
        data_obj
}
#' Generate Hover Plot that shows metadata
#' #' @export
hover_metadata <- function(data_obj, group.by, ...) {
        HoverLocator(plot = UMAPPlot(data_obj, pt.size = 3, group.by = group.by, ...),
                     information = FetchData(data_obj,
                                             vars = c("Accession", "Biosample term name", "Biosample treatment", "Biosample age", "Life stage", "Lab", "Assay title", "Assay title", "Target of assay")))
}
#' Compute standardized variance based from data matrix
#' TODO: refactor this function to a better name
calc_var_std_ver1 <- function(mat, f = 1/10, plot_sample_size = 5000) {
        means = rowMeans(mat)
        means = pmax(means, 0)
        vars = matrixStats::rowVars(mat)
        sel = means !=0 & vars != 0
        
        means = log10(means[sel])
        vars = log10(vars[sel])
        lowess_fit = lowess(means, vars, f = f)
        vars_exp = approxfun(lowess_fit)(means)
        # negative vars_exp may be produced
        # mat_std = (mat[sel, ] - 10^means)/sqrt(10^vars_exp)
        # mat_std = pmin(mat_std, sqrt(ncol(mat_std)))
        # vars_std1 = matrixStats::rowVars(mat_std)
        num_batches = 50
        batch_size = ceiling(length(means)/num_batches)
        vars_std = map(1:num_batches, function(j) {
                message(j)
                indices = (1:length(means))
                indices_batch = which(ceiling(indices / batch_size) == j)
                mat_std = (mat[which(sel)[indices_batch], ] - 10^means[indices_batch])/sqrt(10^vars_exp[indices_batch])
                mat_std = pmin(mat_std, sqrt(ncol(mat_std)))
                vars_std = matrixStats::rowVars(mat_std)
        }) %>% purrr::reduce(c)
        # sample_idx = sample(nrow(mat), size = plot_sample_size)
        # smoothScatter(means[sample_idx], vars[sample_idx])
        # lines(lowess_fit)
        vars_std_all = numeric(nrow(mat))
        vars_std_all[sel] = vars_std
        vars_std_all[!sel] = NA
        vars_std_all
}
calc_var_std_ver2 <- function(mat, f = 1/10, plot_sample_size = 5000) {
        num_batches = 50
        batch_size = ceiling(nrow(mat)/num_batches)
        means = map(1:num_batches, function(j) {
                message(j)
                indices = (1:nrow(mat))
                indices_batch = which(ceiling(indices / batch_size) == j)
                out = rowMeans(mat[indices_batch, ])
        }) %>% purrr::reduce(c)
        vars = map(1:num_batches, function(j) {
                message(j)
                indices = (1:length(means))
                indices_batch = which(ceiling(indices / batch_size) == j)
                out = matrixStats::rowVars(mat[indices_batch, ])
        }) %>% purrr::reduce(c)
        sel = means !=0 & vars != 0
        
        means = log10(means[sel])
        vars = log10(vars[sel])
        lowess_fit = lowess(means, vars, f = f)
        vars_exp = approxfun(lowess_fit)(means)
        # negative vars_exp may be produced
        # mat_std = (mat[sel, ] - 10^means)/sqrt(10^vars_exp)
        # mat_std = pmin(mat_std, sqrt(ncol(mat_std)))
        # vars_std1 = matrixStats::rowVars(mat_std)
        vars_std = map(1:num_batches, function(j) {
                message(j)
                indices = (1:length(means))
                indices_batch = which(ceiling(indices / batch_size) == j)
                mat_std = (mat[which(sel)[indices_batch], ] - 10^means[indices_batch])/sqrt(10^vars_exp[indices_batch])
                mat_std = pmin(mat_std, sqrt(ncol(mat_std)))
                vars_std = matrixStats::rowVars(mat_std)
        }) %>% purrr::reduce(c)
        # sample_idx = sample(nrow(mat), size = plot_sample_size)
        # smoothScatter(means[sample_idx], vars[sample_idx])
        # lines(lowess_fit)
        vars_std_all = numeric(nrow(mat))
        vars_std_all[sel] = vars_std
        vars_std_all[!sel] = NA
        vars_std_all
}
#' Quickly process Seurat object using standard pipeline
quick_process <- function(integrated, dims = 1:15, var_feature = NULL) {
        if (is.null(var_feature)) {
                message("Setting variable features to all features.")
                var_feature = rownames(integrated)
        }
        VariableFeatures(integrated) = var_feature
        integrated <-
                ScaleData(integrated,
                          features = var_feature,
                          verbose = FALSE)
        integrated <-
                RunPCA(integrated, features = var_feature, npcs = 30)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = dims)
}
make_heatmap_col <- function(max_val, mid_val = max_val/2, min_val = 0, col = c("blue", "white", "red")) {
        circlize::colorRamp2(breaks = c(min_val, mid_val, max_val), colors = col)
}
#' push pdf to plot directory
push_pdf <- function(g, file_name, width = 4, height = 4, ps = 10, open_file = T) {
        # get directory to output plots
        if (exists("script_plot_dir")) {
                dir = script_plot_dir
        } else {
                dir = './plots/'
        }
        pdf(paste0(dir, file_name, ".pdf"), width = width, height = height, pointsize = ps)
        tryCatch({
                print(g)
        },
        error = function(e){ 
                print(e)
        })
        dev.off()
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".pdf"))
                shell.exec(file_path)
        }
}
#' push png to plot directory
push_png <- function(g, file_name, width = 4, height = 4, ps = 10, res = 300, open_file = T) {
        if (exists("script_plot_dir")) {
                dir = script_plot_dir
        } else {
                dir = './plots/'
        }
        png(paste0(dir, file_name, ".png"), res = res, units = "in", width = width, height = height, pointsize = ps)
        tryCatch({
                print(g)
        },
        error = function(e){ 
                print(e)
        })
        dev.off()
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".png"))
                shell.exec(file_path)
        }
}
push_widget <- function(w, file_name, open_file = T) {
        if (exists("script_plot_dir")) {
                dir = script_plot_dir
        } else {
                dir = './plots/'
        }
        htmlwidgets::saveWidget(w, file= paste0(dir, file_name, ".html"))
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".html"))
                shell.exec(file_path)
        }
}
push_tsv <- function(tb, file_name, open_file = T) {
        if (exists("script_output_dir")) {
                dir = script_output_dir
        } else {
                dir = './output/'
        }
        write_tsv(tb, file = paste0(dir, file_name, ".tsv"))
        if (open_file) {
                file_path = gsub("/", "\\\\", paste0(dir, file_name, ".tsv"))
                shell.exec(file_path)
        }
}
output_rds_path <- function(file_name) {
        if (exists("script_output_dir")) {
                dir = script_output_dir
        } else {
                dir = './output/'
        }
        paste0(dir, file_name, ".rds")
}

#' check vector is all the same element
check_equal <- function(x) {
        assertthat::assert_that(all(x == x[1]))
        x[1]
}
#' extract anchor tibble from anchor object integrating two objects
make_anchor_tb <- function(anchors, data_obj1, data_obj2) {
        anchor_df = as_tibble(anchors@anchors)
        if (all(c("dataset1", "dataset2") %in% colnames(anchor_df))) {
                anchor_df = anchor_df[anchor_df$dataset1 == 1 & anchor_df$dataset2 == 2, ]
        }
        anchor_df = anchor_df %>% mutate(accession1 = data_obj1$Accession[cell1],
                                         accession2 = data_obj2$Accession[cell2])
        anchor_df = anchor_df %>%
                left_join(data_obj1@meta.data, by = c("accession1" = "Accession")) %>%
                left_join(data_obj2@meta.data, by = c("accession2" = "Accession"), suffix = c("_1", "_2"))
        anchor_df
}
#' strip gene id version
strip_version <- function(x) {
        sapply(strsplit(x, "\\."), "[[", 1)
}
#' subset gtf object
make_genes_gr <- function(gtf, protein_coding = F) {
        out = subset(gtf, type == "gene")
        if (protein_coding) {
                out = subset(out, gene_type == "protein_coding")        
        }
        out
}
#' convert gene ids to gene names
gene_id_to_names <- function(gene_ids, gtf) {
        gtf$gene_name[match(strip_version(gene_ids), strip_version(gtf$gene_id))]
}
#' as_strata
as_strata <- function(x, breaks = c(-Inf, seq(-4.5, 4.5, by = 0.5), Inf)) {
        cut(x, breaks)
}
#' compute hvalue
compute_hval <- function(human_mat, mouse_mat, human_quantile_vec, mouse_quantile_vec, prob_vec) {
        human_hval_ind = do.call(cbind, map(1:length(prob_vec), function(i) {
                rowMeans(human_mat > human_quantile_vec[i]) > prob_vec[i]
        }))
        human_hval = apply(human_hval_ind, 1, function(x) {
                if (!any(x)) {
                        return(0)
                }
                prob_vec[tail(which(x), 1)]
        })
        mouse_hval_ind = do.call(cbind, map(1:length(prob_vec), function(i) {
                rowMeans(mouse_mat > mouse_quantile_vec[i]) > prob_vec[i]
        }))
        mouse_hval = apply(mouse_hval_ind, 1, function(x) {
                if (!any(x)) {
                        return(0)
                }
                prob_vec[tail(which(x), 1)]
        })
        hval = pmin(human_hval, mouse_hval)
        hval
}
#' number of TFBS motifs shared by each human-mouse pair, a motif is present when its row-scaled score > 2
compute_motif_shared <- function(human_motif, mouse_motif, human_id, mouse_id) {
        human_ind = match(human_id, rownames(human_motif))
        mouse_ind = match(mouse_id, rownames(mouse_motif))
        human_uniq = unique(human_ind)
        mouse_uniq = unique(mouse_ind)
        human_bin = t(scale(t(as.matrix(human_motif[human_uniq, , drop = F])))) > 2
        mouse_bin = t(scale(t(as.matrix(mouse_motif[mouse_uniq, , drop = F])))) > 2
        human_bin[is.na(human_bin)] = F
        mouse_bin[is.na(mouse_bin)] = F
        human_ind = match(human_ind, human_uniq)
        mouse_ind = match(mouse_ind, mouse_uniq)
        n_shared = integer(length(human_ind))
        for (k in 1:ncol(human_bin)) {
                n_shared = n_shared + (human_bin[human_ind, k] & mouse_bin[mouse_ind, k])
        }
        n_shared
}
heatmap_wrap <- function(sample_indices) {
        human_mat = t(scale(t(log2(1+human_correct[sample_indices, anchor_tb_sel$accession1]))))
        mouse_mat = t(scale(t(log2(1+mouse_correct[sample_indices, anchor_tb_sel$accession2]))))
        col_ord = hclust(dist(t(rbind(human_mat, mouse_mat))))
        
        Heatmap(human_mat,
                cluster_columns = col_ord,
                show_column_names = F,
                show_row_names = F,
                col = val_col,
                name = "Scaled Count",
                show_row_dend = F,
                show_column_dend = F,
                bottom_annotation = HeatmapAnnotation(Assay = anchor_tb_sel$`Assay title_1`,
                                                      col = list(Assay = assay_col_hu),
                                                      annotation_label = "Human Assay",
                                                      show_annotation_name = F)
        ) + 
                Heatmap(mouse_mat,
                        cluster_columns = col_ord,
                        show_column_names = F,
                        show_row_names = F,
                        col = val_col,
                        name = "Scaled Count",
                        show_row_dend = F,
                        show_column_dend = F,
                        bottom_annotation = HeatmapAnnotation(Assay =  anchor_tb_sel$`Assay title_2`,
                                                              col = list(Assay = assay_col_mm),
                                                              annotation_label = "Mouse_Assay",
                                                              show_annotation_name = F))#,
        # right_annotation = rowAnnotation(Component = regions$component[sample_indices0],
        #                                  col = list(Component = comp_col)))
}

