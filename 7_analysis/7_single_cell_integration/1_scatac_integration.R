# single cell integration
library(tidyverse)
library(Seurat)
library(Signac)
library(plotly)
save_path = "/dcs04/hongkai/data/rzhao/sciATAC/"
source("/dcs04/hongkai/data/rzhao/sciATAC/sciATAC_human_mouse_immune/load_mouse_sciatac2.R")
source("./scatac_helper.R")

regions = readRDS("/dcs04/hongkai/data/wfang/indexed_dhs_mapped_regions.rds")
path_data = "/dcl01/hongkai/data/wfang/human_sciatac_ren/"
meta = read_tsv(paste0(path_data,"meta.tsv"))

data_list = map(c("lung_SM-A62E9", "lung_SM-A8WNH",
                  "liver_SM-A8WNZ", "liver_SM-A8WNZ",
                  "heart_lv_SM-IOBHO", "heart_lv_SM-JF1NY", "heart_atrial_appendage_SM-JF1NX", "heart_atrial_appendage_SM-IOBHN",
                  "colon_transverse_SM-A9HOW", "colon_transverse_SM-A9VP4",
                  "colon_transverse_SM-ACCQ1", "colon_transverse_SM-BZ2ZS", "colon_transverse_SM-CSSDA",
                  "colon_sigmoid_SM-AZPYO", "colon_sigmoid_SM-JF1O8"),
                load_ren_sample)

human_all = SeuratObject:::merge.Seurat(data_list[[1]],
                                        data_list[2:length(data_list)])

list_mouse_samples()
mouse_data_list = map(c("Lung1_62216", "Lung2_62216", "Liver_62016", "HeartA_62816",
                        "LargeIntestineA_62816", "LargeIntestineB_62816"),
                      load_mouse_sample, min.cells.pct = 0.00)
mouse_all = SeuratObject:::merge.Seurat(mouse_data_list[[1]],
                                        mouse_data_list[2:length(mouse_data_list)])

mouse_sample_obj = mouse_all[, mouse_all$cell_type %in% unlist(matched_label)]
mouse_sample_obj$species = "mouse"
human_sample_obj = human_all[, human_all$cell_type %in% unlist(matched_label)]
human_sample_obj$species = "human"
mouse_sample_obj@assays$DHS@counts@x[mouse_sample_obj@assays$DHS@counts@x > 0] = 1
human_sample_obj@assays$DHS@counts@x[human_sample_obj@assays$DHS@counts@x > 0] = 1

mouse_sample_obj = mouse_sample_obj[, colSums(mouse_sample_obj@assays$DHS@counts) >= 900]
human_sample_obj = human_sample_obj[, colSums(human_sample_obj@assays$DHS@counts) >= 900]

common_features = intersect(rownames(human_sample_obj),
                            rownames(mouse_sample_obj))
# overall non-zero feature ranks
feature_rank = rank(rowSums(human_sample_obj@assays$DHS@counts[common_features, ])) + 
        rank(rowSums(mouse_sample_obj@assays$DHS@counts[common_features, ]))
# intitial filtering to get candidate features
dhs_common = common_features[order(feature_rank, decreasing = T)[1:300000]]

# Integration
# set different feature numbers and feature filtering
num_features = 48160
# num_features = 25000
# method 1: CO-V filtering
method = "method1"
dhs_sel = dhs_common[order(regions$spearman_global[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# method 2: raw feature rank
method = "method2"
dhs_sel = common_features[order(feature_rank, decreasing = T)[1:num_features]]
# method 3: PhastCons
method = "method3"
dhs_sel = dhs_common[order(regions$phastCons4way[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# method 4: PhastCons
method = "method4_phyloP4way"
dhs_sel = dhs_common[order(regions$phyloP4way[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# method 5: CO-V Manual
method = "method5_spearman_manual"
dhs_sel = dhs_common[order(regions$spearman_manual[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# method 6: LECIF
method = "method6_lecif"
dhs_sel = dhs_common[order(regions$lecif[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# method 7: PSH
method = "method7_percentage"
dhs_sel = dhs_common[order(regions$percentage[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# methods 8-11: within-species variability and AdaLiftOver
baseline_tb = readRDS("./output/chromatin_conservation_v1/CA_species_baselines.rds")
adaliftover_tb = readRDS("./output/adaliftover/adaliftover_scores.rds")
regions$human_CA_var_std = baseline_tb$human_CA_var_std[match(regions$id, baseline_tb$id)]
regions$mouse_CA_var_std = baseline_tb$mouse_CA_var_std[match(regions$id, baseline_tb$id)]
regions$adaliftover_h2m = adaliftover_tb$AdaLiftOver_hg38_to_mm10[match(1:nrow(regions), adaliftover_tb$human_idx)]
regions$adaliftover_m2h = adaliftover_tb$AdaLiftOver_mm10_to_hg38[match(1:nrow(regions), adaliftover_tb$human_idx)]
method = "method8_human_var"
dhs_sel = dhs_common[order(regions$human_CA_var_std[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
method = "method9_mouse_var"
dhs_sel = dhs_common[order(regions$mouse_CA_var_std[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
method = "method10_adaliftover_h2m"
dhs_sel = dhs_common[order(regions$adaliftover_h2m[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
method = "method11_adaliftover_m2h"
dhs_sel = dhs_common[order(regions$adaliftover_m2h[match(dhs_common, regions$id)], decreasing = T)[1:num_features]]
# end specifying features for methods

human_mouse_obj = SeuratObject:::merge.Seurat(human_sample_obj[common_features, ], mouse_sample_obj[common_features, ])
# there are too many hepatocytes in mouse data, subsampling here
cells_sel = colnames(mouse_sample_obj)[mouse_sample_obj$cell_type == "Hepatocytes"]
# cells_sel = cells_sel[order(mouse_sample_obj$nCount_DHS[cells_sel], decreasing = F)[1:3500]]
set.seed(37)
cells_sel = sample(cells_sel, 2500, replace = F)
human_mouse_obj = human_mouse_obj[, !colnames(human_mouse_obj) %in% cells_sel]

# begin example
#x = integrate_harmony(dhs_sel, nclust = 25, max.iter.harmony = 8, max.iter.cluster = 15)
#x <- RunUMAP(x, dims = 2:30, reduction = 'harmony')
human_mouse_obj_list <- SplitObject(human_mouse_obj, split.by = "species")

anchors_human_mouse = FindIntegrationAnchors(
        human_mouse_obj_list,
        anchor.features = dhs_sel,
        reduction = "cca",
        k.anchor = 40, # sc need to be larger # changedhere
        k.filter = 100, #100-200
        dims = 1:30)

save_path = "/dcs04/hongkai/data/rzhao/sciATAC/"
x = IntegrateData(anchors_human_mouse,
                  features.to.integrate = dhs_sel,
                  dims = 1:30,
                  k.weight = 100)
x@assays$DHS<-NULL
#x = quick_process(x, var_feature = dhs_sel, dims = 1:30)
x = process_scatac(x, features = rownames(x), run_umap = T, normalize = F)
### newly added 
#int_obj = process_scatac(int_obj, features = rownames(int_obj), run_umap = T, normalize = F)

save_path = "/dcs04/hongkai/data/rzhao/sciATAC/"
saveRDS(x, paste0(save_path,'mid_data/integrated_human_mouse_','example1_',method,'.rds'))
