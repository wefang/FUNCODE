library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
source("./helper/helper.R")
source("./hic_helper_func.R")

print(load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda"))
print(load("metadata_processed/indexed_aligned_combined.rda"))
regions_all = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
regions_noalign = readRDS("./output/nonalign_sign_combined.rds")
regions_noalign$motif_ind = regions_noalign$n_motif_shared >= 1

mm10_region_indexed = rlang::duplicate(mm10_regions)
mm10_region_indexed$identifier = paste0("align_", 1:length(mm10_region_indexed))
mouse_regions$identifier = paste0("MouseDHS-", mouse_regions$identifier)
mm10_region_indexed = c(mm10_region_indexed, mouse_regions)

eval_name = "Tcell"
dir.create(paste0("./processed_data/intact_hic_eval/", eval_name, "/"))
human_t_cell_types = c("CD4-positive, alpha-beta T cell",
                       "activated CD4-positive, alpha-beta T cell")
# insitu Hi-C loops
mouse_data_dir = "./processed_data/mouse_insitu_hic_loop_cell_types/Tcells/"
mouse_data_files = list.files(mouse_data_dir, full = T)
mouse_loop_list = map(mouse_data_files, rtracklayer::import)
mouse_loop_gr = do.call(c, mouse_loop_list)
mouse_tcells_loop_gr = sort(mouse_loop_gr)

eval_name = "Bcell"
dir.create(paste0("./processed_data/intact_hic_eval/", eval_name, "/"))
human_b_cell_types = c("B cell")
mouse_data_dir = "./processed_data/mouse_insitu_hic_loop_cell_types/Bcells/"
mouse_data_files = list.files(mouse_data_dir, full = T)
mouse_loop_list = map(mouse_data_files, rtracklayer::import)
mouse_loop_gr = do.call(c, mouse_loop_list)
mouse_bcells_loop_gr = sort(mouse_loop_gr)

meta_tb = read_tsv("./processed_data/intact_hic_loops/intact_hic_file_metadata.txt")
meta_tb = filter(meta_tb, `Output type` == "long range chromatin interactions")
download_dir = "./processed_data/intact_hic_loops/"

all_human_cell_types = sort(unique(meta_tb$`Biosample term name`))
# Load CTCF binding motifs
ctcf_tb = read_tsv("processed_data/mapeed_motif/map_MA0139.1_CTCF.txt")
ctcf_tb = read_tsv("processed_data/mapeed_motif/map_MA0139.1_CTCF.txt")
ctcf_gr <- GRanges(
        seqnames = Rle(ctcf_tb$chromosome),
        ranges = IRanges(start = ctcf_tb$start, end = ctcf_tb$end),
        strand = Rle(ctcf_tb$strand))
# Predictions
for (cell_type in all_human_cell_types) {
        message(cell_type)
        human_meta_tb = filter(meta_tb, `Biosample term name` %in% cell_type)
        human_loop_files = paste0(download_dir, basename(human_meta_tb$`File download URL`))
        human_loop_gr = map(human_loop_files, rtracklayer::import)
        human_loop_gr = do.call(c, human_loop_gr)
        human_loop_gr = sort(human_loop_gr)

        # human generate DHS pairs based on loops
        high_res_start = as.numeric(human_loop_gr@elementMetadata[["NA..14"]])
        high_res_end = as.numeric(human_loop_gr@elementMetadata[["NA..15"]])

        high_res_start1 = as.numeric(human_loop_gr@elementMetadata[["NA..16"]])
        high_res_end1 = as.numeric(human_loop_gr@elementMetadata[["NA..17"]])

        high_res_indices = which(!is.na(high_res_start))
        high_res_loop_first = GRanges(seqnames(human_loop_gr@first)[high_res_indices],
                                      IRanges(start = high_res_start[high_res_indices],
                                              end = high_res_end[high_res_indices]))
        high_res_loop_second = GRanges(seqnames(human_loop_gr@second)[high_res_indices],
                                       IRanges(start = high_res_start1[high_res_indices],
                                               end = high_res_end1[high_res_indices]))

        first_region_ol = findOverlaps(high_res_loop_first, subject = human_regions)
        second_region_ol = findOverlaps(high_res_loop_second, subject = human_regions)

        first_region_sp = split(subjectHits(first_region_ol),
                                queryHits(first_region_ol))
        second_region_sp = split(subjectHits(second_region_ol),
                                 queryHits(second_region_ol))
        human_loop_dhs_pairs = bind_rows(map(1:length(high_res_loop_first), function(i) {
                if (length(first_region_sp[[as.character(i)]]) > 0 &
                    length(second_region_sp[[as.character(i)]]) > 0) {
                        return(
                                as_tibble(expand.grid(dhs1 = first_region_sp[[as.character(i)]],
                                                      dhs2 = second_region_sp[[as.character(i)]]))
                        )
                } else {
                        return(NULL)
                }
        }))
        human_loop_dhs_pairs = distinct(human_loop_dhs_pairs)
        saveRDS(human_loop_dhs_pairs, file = paste0("./processed_data/human_intact_hic_dhs_pairs/", cell_type, ".rds"))
        # end human generate DHS pairs based on loops
        human_loop_dhs_pairs$dhs1_id = human_regions$identifier[human_loop_dhs_pairs$dhs1]
        human_loop_dhs_pairs$dhs2_id = human_regions$identifier[human_loop_dhs_pairs$dhs2]
        human_loop_dhs_pairs$dhs1_cat = decide_dhs_category(human_loop_dhs_pairs$dhs1_id)
        human_loop_dhs_pairs$dhs2_cat = decide_dhs_category(human_loop_dhs_pairs$dhs2_id)

        human_loop_dhs_pairs$dhs1_index = map_pair_mouse_index(human_loop_dhs_pairs$dhs1_id, human_loop_dhs_pairs$dhs1_cat)
        human_loop_dhs_pairs$dhs2_index = map_pair_mouse_index(human_loop_dhs_pairs$dhs2_id, human_loop_dhs_pairs$dhs2_cat)

        human_loop_dhs_mapped = human_loop_dhs_pairs[human_loop_dhs_pairs$dhs1_cat != "other" & human_loop_dhs_pairs$dhs2_cat != "other", ]

        batch_ind = ceiling(1:nrow(human_loop_dhs_mapped) / 1e4)
        all_loop_ind = do.call(c, map(1:max(batch_ind), function(i) {
                message(i)
                loop_ind = check_loop(resize(mm10_region_indexed[match(human_loop_dhs_mapped$dhs1_index[batch_ind == i], mm10_region_indexed$identifier)],
                                             width = 5000, fix = "center"),
                                      resize(mm10_region_indexed[match(human_loop_dhs_mapped$dhs2_index[batch_ind == i], mm10_region_indexed$identifier)],
                                             width = 5000, fix = "center"),
                                      mouse_tcells_loop_gr@first,
                                      mouse_tcells_loop_gr@second)
        }))
        human_loop_dhs_mapped$tcell_loop_ind = all_loop_ind
        
        batch_ind = ceiling(1:nrow(human_loop_dhs_mapped) / 1e4)
        all_loop_ind = do.call(c, map(1:max(batch_ind), function(i) {
                message(i)
                loop_ind = check_loop(resize(mm10_region_indexed[match(human_loop_dhs_mapped$dhs1_index[batch_ind == i], mm10_region_indexed$identifier)],
                                             width = 5000, fix = "center"),
                                      resize(mm10_region_indexed[match(human_loop_dhs_mapped$dhs2_index[batch_ind == i], mm10_region_indexed$identifier)],
                                             width = 5000, fix = "center"),
                                      mouse_bcells_loop_gr@first,
                                      mouse_bcells_loop_gr@second)
        }))
        human_loop_dhs_mapped$bcell_loop_ind = all_loop_ind
        
        human_loop_dhs_mapped$ctcf1 = countOverlaps(mm10_region_indexed[match(human_loop_dhs_mapped$dhs1_index, mm10_region_indexed$identifier)],
                                      ctcf_gr)
        human_loop_dhs_mapped$ctcf2 = countOverlaps(mm10_region_indexed[match(human_loop_dhs_mapped$dhs2_index, mm10_region_indexed$identifier)],
                                      ctcf_gr)
        human_loop_dhs_mapped$dhs1_motif = human_loop_dhs_mapped$ctcf1 > 0
        human_loop_dhs_mapped$dhs2_motif = human_loop_dhs_mapped$ctcf2 > 0
        
        human_loop_dhs_mapped$dist = distance(mm10_region_indexed[match(human_loop_dhs_mapped$dhs1_index, mm10_region_indexed$identifier)],
                                              mm10_region_indexed[match(human_loop_dhs_mapped$dhs2_index, mm10_region_indexed$identifier)])
        saveRDS(human_loop_dhs_mapped,
                file = paste0("./processed_data/intact_hic_eval/final/human_", cell_type, ".rds"))
}

# Collect results for all cell types
cell_type_files = list.files("./processed_data/intact_hic_eval/final/")
tb_comb = bind_rows(map(cell_type_files, function(x) {
        tb = readRDS(paste0("./processed_data/intact_hic_eval/final/", x))
        tb$cell_type = x
        tb
}))
tb_comb$ctcf1 = countOverlaps(mm10_region_indexed[match(tb_comb$dhs1_index, mm10_region_indexed$identifier)],
                              ctcf_gr)
tb_comb$ctcf2 = countOverlaps(mm10_region_indexed[match(tb_comb$dhs2_index, mm10_region_indexed$identifier)],
                              ctcf_gr)
tb_comb$dhs1_motif = tb_comb$ctcf1 > 0
tb_comb$dhs2_motif = tb_comb$ctcf2 > 0
tb_comb$dist = distance(mm10_region_indexed[match(tb_comb$dhs1_index, mm10_region_indexed$identifier)],
                        mm10_region_indexed[match(tb_comb$dhs2_index, mm10_region_indexed$identifier)])

lvls = c("align_sign", "align_sign_cov", "align", "noalign_sign_motif",  "noalign_sign")
tb_comb = tb_comb %>% 
        mutate(dhs1_cat = factor(dhs1_cat, levels = lvls),
               dhs2_cat = factor(dhs2_cat, levels = lvls)) %>%
        mutate(dhs1_cat_int = as.integer(dhs1_cat),
               dhs2_cat_int = as.integer(dhs2_cat),
               dhs1_cat_ord_int = pmin(dhs1_cat_int, dhs2_cat_int),
               dhs2_cat_ord_int = pmax(dhs1_cat_int, dhs2_cat_int)) %>%
        mutate(dhs1_cat_ord = factor(dhs1_cat_ord_int, levels = 1:length(lvls), labels = lvls),
               dhs2_cat_ord = factor(dhs2_cat_ord_int, levels = 1:length(lvls), labels = lvls)) %>%
        mutate(dhs1_id_ord = map2_chr(dhs1_id, dhs2_id, function(a, b) {
                ifelse(a > b, b, a)
        }),
        dhs2_id_ord = map2_chr(dhs1_id, dhs2_id, function(a, b) {
                ifelse(a > b, a, b)
        }))

# Summary statistics: All predictions, no distance constraint
tb_comb %>% filter(dhs1_motif | dhs2_motif) %>%
        filter(dhs2_cat_ord %in% c("align_sign", "align_sign_cov")) %>%
        select(dhs1_id_ord, dhs2_id_ord) %>% distinct()
# Compile predicted loops
pred_all = tb_comb %>% filter(dhs1_motif | dhs2_motif) %>%
        filter(dhs2_cat_ord %in% c("align_sign", "align_sign_cov"))
pred_compiled = select(pred_all, cell_type, dhs1_id, dhs2_id, dhs1_index, dhs2_index, dhs1_motif, dhs2_motif)
pred_compiled$cell_type = str_match(pred_compiled$cell_type, "human_(.*?)\\.rds")[, 2]
pred_compiled$dhs1_index = as.numeric(map_chr(strsplit(pred_compiled$dhs1_index, "_"), 2))
pred_compiled$dhs2_index = as.numeric(map_chr(strsplit(pred_compiled$dhs2_index, "_"), 2))
assay_name_mapping = c("dnase_atac" = "chromatin_accessibility",
                       "H3K27ac" ="H3K27ac",
                       "H3K4me1" = "H3K4me1",
                       "H3K4me3" = "H3K4me3")
cols_append = c("mouse_id", "mm10_region", "mouse_cre_class_prio_v4",
                "human_id", "hg38_region",
                paste0("cov_", assay_name_mapping, "_sign"),
                paste0("cob_", assay_name_mapping, "_sign"))
for (x in cols_append) {
        pred_compiled[[paste0("anchor1_", x)]] = regions_all[[x]][pred_compiled$dhs1_index]
        pred_compiled[[paste0("anchor2_", x)]] = regions_all[[x]][pred_compiled$dhs2_index]
}
# drop self-loops, anchors without a conserved pair, and mouse anchors on different chromosomes
anchor1_sign = rowSums(as.matrix(pred_compiled[paste0("anchor1_", cols_append[6:13])]))
anchor2_sign = rowSums(as.matrix(pred_compiled[paste0("anchor2_", cols_append[6:13])]))
same_chr = str_extract(pred_compiled$anchor1_mm10_region, "^[^:]+") == str_extract(pred_compiled$anchor2_mm10_region, "^[^:]+")
pred_compiled = pred_compiled[pred_compiled$dhs1_id != pred_compiled$dhs2_id & anchor1_sign > 0 & anchor2_sign > 0 & same_chr, ]
readr::write_csv(pred_compiled, "table_s10_intact_pred_compiled.csv")

# Filter loops by distance and Motifs for evaluations
tb_comb1 =  tb_comb %>% filter(dist > 65000) %>%
        filter(dhs1_motif | dhs2_motif)

# Cell type specific evaluations
library(ComplexHeatmap)
# T cell specific evaluations
all_loop_summarize = tb_comb1 %>%
        filter(cell_type %in% paste0("human_", human_t_cell_types, ".rds")) %>%
        group_by(dhs1_id_ord, dhs2_id_ord) %>%
        filter(row_number() <= 1) %>%
        group_by(dhs1_cat_ord, dhs2_cat_ord) %>%
        summarise(mean_loop = mean(tcell_loop_ind),
                  total = n(),
                  .groups = 'drop')
arrange(all_loop_summarize, desc(mean_loop))
# Null set statistics for estimating the average FDR
tb_comb1 %>%
        filter(dhs2_cat_ord %in% c("align_sign", "align_sign_cov")) %>%
        filter(cell_type %in% paste0("human_", human_t_cell_types, ".rds")) %>%
        group_by(dhs1_id_ord, dhs2_id_ord) %>%
        filter(row_number() <= 1) %>%
        ungroup() %>%
        summarise(mean_loop = mean(tcell_loop_ind),
                  total = n(),
                  .groups = 'drop')
null_loop_summarize =  tb_comb %>% filter(dist > 65000)  %>%
        # filter(dhs1_motif | dhs2_motif) %>%
        filter(cell_type %in% paste0("human_", human_t_cell_types, ".rds")) %>%
        filter(dhs1_cat_ord %in% c("noalign_sign", "noalign_sign_motif") & dhs2_cat_ord %in% c("noalign_sign", "noalign_sign_motif")) %>% 
        group_by(dhs1_id_ord, dhs2_id_ord) %>%
        filter(row_number() <= 1) %>% 
        ungroup() %>%
        summarise(mean_loop = mean(tcell_loop_ind),
                  total = n(),
                  .groups = 'drop')
null_loop_summarize
# End estimating average FDR

# Heatmap of evaluation accuracies for T cell
interact_mat = matrix(NA, 5, 5)
rownames(interact_mat) = colnames(interact_mat) = levels(all_loop_summarize$dhs1_cat_ord)
interact_mat[cbind(all_loop_summarize$dhs2_cat_ord, all_loop_summarize$dhs1_cat_ord)] = all_loop_summarize$mean_loop
interact_mat[cbind(all_loop_summarize$dhs2_cat_ord, all_loop_summarize$dhs1_cat_ord)[all_loop_summarize$total <= 100, ]] = NA
h = Heatmap(interact_mat, cluster_rows = F, cluster_columns = F,
            name = "Validation rate",
            column_title = "CD4+ T cells")
push_pdf(h,
         w = 5, h = 4,
         file = "Tcells_eval")
# Summarize B cell evaluations
all_loop_summarize = tb_comb1 %>%
        filter(cell_type %in% paste0("human_", human_b_cell_types, ".rds")) %>%
        group_by(dhs1_id_ord, dhs2_id_ord) %>%
        filter(row_number() <= 1) %>%
        group_by(dhs1_cat_ord, dhs2_cat_ord) %>%
        summarise(mean_loop = mean(bcell_loop_ind),
                  total = n(),
                  .groups = 'drop')
arrange(all_loop_summarize, desc(mean_loop))
# B cell FDR
tb_comb1 %>%
        filter(dhs2_cat_ord %in% c("align_sign", "align_sign_cov")) %>%
        filter(cell_type %in% paste0("human_", human_b_cell_types, ".rds")) %>%
        group_by(dhs1_id_ord, dhs2_id_ord) %>%
        filter(row_number() <= 1) %>%
        ungroup() %>%
        summarise(mean_loop = mean(tcell_loop_ind),
                  total = n(),
                  .groups = 'drop')
null_loop_summarize =  tb_comb %>% filter(dist > 65000)  %>%
        # filter(dhs1_motif | dhs2_motif) %>%
        filter(cell_type %in% paste0("human_", human_b_cell_types, ".rds")) %>%
        filter(dhs1_cat_ord %in% c("noalign_sign", "noalign_sign_motif") & dhs2_cat_ord %in% c("noalign_sign", "noalign_sign_motif")) %>% 
        group_by(dhs1_id_ord, dhs2_id_ord) %>%
        filter(row_number() <= 1) %>% 
        ungroup() %>%
        summarise(mean_loop = mean(tcell_loop_ind),
                  total = n(),
                  .groups = 'drop')
null_loop_summarize
# End B cell FDR 

# Heatmap of evaluation accuracies for B cell
interact_mat = matrix(NA, 5, 5)
rownames(interact_mat) = colnames(interact_mat) = levels(all_loop_summarize$dhs1_cat_ord)
interact_mat[cbind(all_loop_summarize$dhs2_cat_ord, all_loop_summarize$dhs1_cat_ord)] = all_loop_summarize$mean_loop
interact_mat[cbind(all_loop_summarize$dhs2_cat_ord, all_loop_summarize$dhs1_cat_ord)[all_loop_summarize$total <= 50, ]] = NA
h = Heatmap(interact_mat, cluster_rows = F, cluster_columns = F,
            name = "Validation rate",
            column_title = "B cells")
push_pdf(h,
         w = 5, h = 4,
         file = "Bcells_eval")

#### Visualize loops on genomic tracks ####
mouse_range_gr <- GRanges(
        seqnames = Rle("chr7"),
        ranges = IRanges(start = 99000001, end = 100000000)
)
# Loops with both anchors conserved
tb_pred = tb_comb1 %>%
        filter(cell_type %in% paste0("human_", human_t_cell_types, ".rds")) %>%
        select(-cell_type) %>%
        distinct() %>%
        filter(dhs1_cat_ord %in% c("align_sign", "align_sign_cov") &
                       dhs2_cat_ord %in% c("align_sign", "align_sign_cov"))
tb_pred$gr_ind1 = countOverlaps(mm10_region_indexed[match(tb_pred$dhs1_index, mm10_region_indexed$identifier)], mouse_range_gr)
tb_pred$gr_ind2 = countOverlaps(mm10_region_indexed[match(tb_pred$dhs2_index, mm10_region_indexed$identifier)], mouse_range_gr)
# Loops with zero or one anchor conserved
tb_pred1 = tb_comb %>%
        filter(cell_type %in% paste0("human_", human_t_cell_types, ".rds")) %>%
        select(-cell_type) %>%
        distinct() %>%
        filter((!(dhs1_cat_ord %in% c("align_sign", "align_sign_cov") &
                         dhs2_cat_ord %in% c("align_sign", "align_sign_cov"))))
tb_pred1$gr_ind1 = countOverlaps(mm10_region_indexed[match(tb_pred1$dhs1_index, mm10_region_indexed$identifier)], mouse_range_gr)
tb_pred1$gr_ind2 = countOverlaps(mm10_region_indexed[match(tb_pred1$dhs2_index, mm10_region_indexed$identifier)], mouse_range_gr)

# Plot in situ Hi-C loops
all_loop_gr = mouse_tcells_loop_gr
loop_gr_choice = intersect(queryHits(findOverlaps(all_loop_gr@first, mm10_region_indexed)),
                           queryHits(findOverlaps(all_loop_gr@second, mm10_region_indexed)))
all_loop_gr = all_loop_gr[loop_gr_choice]
all_loop_dist = calculateDistances(all_loop_gr@first, all_loop_gr@second)
all_loop_gr = all_loop_gr[all_loop_dist < 4e5 & all_loop_dist > 6e4]
first_ol = queryHits(findOverlaps(all_loop_gr@first, mouse_range_gr))
second_ol = queryHits(findOverlaps(all_loop_gr@second, mouse_range_gr))
hic_tb = as_tibble(all_loop_gr[intersect(first_ol, second_ol)])
g = hic_tb[c("first.start", "second.start")] %>%
        graph_from_data_frame(directed = F)

layout_mod = cbind(as.numeric(V(g)$name), 1)
rownames(layout_mod) = V(g)$name
g_hic_insitu = ggraph(g, layout = layout_mod) +
        # geom_node_point() +
        geom_edge_arc2(fold = T, color = "#3690c0") +
        xlim(c(start(mouse_range_gr)- 10000, end(mouse_range_gr)+ 10000)) +
        theme(panel.background = element_blank(),
              axis.ticks.x = element_line(),
              axis.line.x = element_line(),
              axis.text.x = element_text())
g_hic_insitu
# Plot loops with both anchors conserved
all_loop_dhs_pairs = tb_pred
all_loop_dhs_pairs$dhs1_int = map_dbl(strsplit(all_loop_dhs_pairs$dhs1_index, "_"), function(x) as.numeric(x[2]))
all_loop_dhs_pairs$dhs2_int = map_dbl(strsplit(all_loop_dhs_pairs$dhs2_index, "_"), function(x) as.numeric(x[2]))
dhs1_coord = as_tibble(mm10_region_indexed[match(all_loop_dhs_pairs$dhs1_index, mm10_region_indexed$identifier)]) %>%
        transmute(dhs1_seqnames = seqnames,
                  dhs1_start = start,
                  dhs1_end = end)
dhs2_coord = as_tibble(mm10_region_indexed[match(all_loop_dhs_pairs$dhs2_index, mm10_region_indexed$identifier)]) %>%
        transmute(dhs2_seqnames = seqnames,
                  dhs2_start = start,
                  dhs2_end = end)
all_loop_dhs_pairs = bind_cols(all_loop_dhs_pairs, dhs1_coord, dhs2_coord)
gene_loop_sel = intersect(queryHits(findOverlaps(mm10_region_indexed[match(all_loop_dhs_pairs$dhs1_index, mm10_region_indexed$identifier)],
                                                 mouse_range_gr)),
                          queryHits(findOverlaps(mm10_region_indexed[match(all_loop_dhs_pairs$dhs2_index, mm10_region_indexed$identifier)],
                                                 mouse_range_gr)))
loop_dhs_sel = all_loop_dhs_pairs[gene_loop_sel, ]
g = loop_dhs_sel %>%
        dplyr::select(dhs1_id, dhs2_id, tcell_loop_ind) %>%
        graph_from_data_frame(directed = F)
layout_mod = cbind(c(loop_dhs_sel$dhs1_start, loop_dhs_sel$dhs2_start)[match(V(g)$name, c(loop_dhs_sel$dhs1_id, loop_dhs_sel$dhs2_id))], 1)
rownames(layout_mod) = V(g)$name
g_hic = ggraph(g, layout = layout_mod) +
        # geom_node_point() +
        geom_edge_arc2(aes(color = tcell_loop_ind), fold = T) +
        scale_edge_color_manual(values = c("#8B0000", "#228B22")) +
        xlim(c(start(mouse_range_gr)- 10000, end(mouse_range_gr)+ 10000)) +
        theme(panel.background = element_blank(),
              axis.line.x = element_line(),
              axis.text.x = element_text(),
              axis.ticks.x = element_line())
g_hic
# Plot loops with zero or one anchor conserved
all_loop_dhs_pairs = tb_pred1
all_loop_dhs_pairs$dhs1_int = map_dbl(strsplit(all_loop_dhs_pairs$dhs1_index, "_"), function(x) as.numeric(x[2]))
all_loop_dhs_pairs$dhs2_int = map_dbl(strsplit(all_loop_dhs_pairs$dhs2_index, "_"), function(x) as.numeric(x[2]))
dhs1_coord = as_tibble(mm10_region_indexed[match(all_loop_dhs_pairs$dhs1_index, mm10_region_indexed$identifier)]) %>%
        transmute(dhs1_seqnames = seqnames,
                  dhs1_start = start,
                  dhs1_end = end)
dhs2_coord = as_tibble(mm10_region_indexed[match(all_loop_dhs_pairs$dhs2_index, mm10_region_indexed$identifier)]) %>%
        transmute(dhs2_seqnames = seqnames,
                  dhs2_start = start,
                  dhs2_end = end)
all_loop_dhs_pairs = bind_cols(all_loop_dhs_pairs, dhs1_coord, dhs2_coord)
gene_loop_sel = intersect(queryHits(findOverlaps(mm10_region_indexed[match(all_loop_dhs_pairs$dhs1_index, mm10_region_indexed$identifier)],
                                                 mouse_range_gr)),
                          queryHits(findOverlaps(mm10_region_indexed[match(all_loop_dhs_pairs$dhs2_index, mm10_region_indexed$identifier)],
                                                 mouse_range_gr)))
loop_dhs_sel = all_loop_dhs_pairs[gene_loop_sel, ]
g = loop_dhs_sel %>%
        dplyr::select(dhs1_id, dhs2_id, tcell_loop_ind) %>%
        graph_from_data_frame(directed = F)
layout_mod = cbind(c(loop_dhs_sel$dhs1_start, loop_dhs_sel$dhs2_start)[match(V(g)$name, c(loop_dhs_sel$dhs1_id, loop_dhs_sel$dhs2_id))],                   1)
rownames(layout_mod) = V(g)$name
g_hic1 = ggraph(g, layout = layout_mod) +
        # geom_node_point() +
        geom_edge_arc2(aes(color = tcell_loop_ind), fold = T) +
        scale_edge_color_manual(values = c("#8B0000", "#228B22")) +
        xlim(c(start(mouse_range_gr)- 10000, end(mouse_range_gr)+ 10000)) +
        theme(panel.background = element_blank(),
              axis.line.x = element_line(),
              axis.text.x = element_text(),
              axis.ticks.x = element_line())
g_hic1
# Combining all three examples
g_hic1 / g_hic / g_hic_insitu
#### end plot ####