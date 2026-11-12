library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)

source("./helper/helper.R")
source("./helper/def_color.R")

# load result files
sign_tb = readRDS("./output/nonalign_sign_combined.rds")
ghcob = readRDS("./output/nonalign_ghcob_motif_shared.rds")
regions_all = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
print(load("metadata_processed/DHS/human_mouse_dhs_raw.rda"))

human_dhs_tb = tibble(id = human_regions$identifier)
human_cob_id = regions_all$human_id[regions_all$cob_any]
human_cov_id = regions_all$human_id[regions_all$cov_any]
human_dhs_tb$sign_type = "non-alignable non-conserved"
human_dhs_tb$sign_type[human_dhs_tb$id %in% regions_all$human_id] = "Alignable non-conserved"
human_dhs_tb$sign_type[human_dhs_tb$id %in% unique(sign_tb$human_id)] = "non-alignable CO-V"
human_dhs_tb$sign_type[human_dhs_tb$id %in% unique(sign_tb$human_id[sign_tb$n_motif_shared >= 1])] = "non-alignable CO-V motif"
human_dhs_tb$sign_type[human_dhs_tb$id %in% ghcob$human_no_motif_ids] = "non-alignable CO-B"
human_dhs_tb$sign_type[human_dhs_tb$id %in% ghcob$human_motif_ids] = "non-alignable CO-B motif"
human_dhs_tb$sign_type[human_dhs_tb$id %in% human_cob_id] = "Alignable CO-B"
human_dhs_tb$sign_type[human_dhs_tb$id %in% human_cov_id] = "Alignable CO-V"
sign_type_levels = c("non-alignable non-conserved", "Alignable non-conserved",
                     "non-alignable CO-V", "non-alignable CO-V motif",
                     "non-alignable CO-B", "non-alignable CO-B motif",
                     "Alignable CO-B", "Alignable CO-V")
human_dhs_tb$sign_type = factor(human_dhs_tb$sign_type, levels = sign_type_levels)
# Figure 2e: Pie chart breaking down conservation types
human_dhs_tb %>% count(sign_type) %>%
        mutate(pct = n / sum(n) * 100) %>%
        ggpie(x = "n", fill = "sign_type") + 
        scale_fill_manual(values = sign_type_col) +
        theme(legend.position = "right")

# load cCRE annotations
print(load("./metadata_processed/DHS_cCRE/all_DHS_cCRE_v4.rda"))
human_dhs_tb = left_join(human_dhs_tb, select(human_cre_tb, id, human_cre_class_prio))
human_dhs_tb$human_cre_class_prio[is.na(human_dhs_tb$human_cre_class_prio)] = "non-cCRE"
saveRDS(human_dhs_tb, file = "./output/dhs_level_tb.rds")

human_dhs_tb = readRDS("./output/dhs_level_tb.rds")
# Figure 2g
cre_prob_tb = human_dhs_tb %>% group_by(sign_type) %>%
        count(human_cre_class_prio) %>% mutate(frac = n / sum(n))
cre_prob_tb %>% ggbarplot(x = "sign_type", y = "frac", position = position_stack(), fill = "human_cre_class_prio", color = NA) +
        scale_fill_manual(values = cCRE_v4_col) + coord_flip()
sign_prob_tb = human_dhs_tb %>% group_by(human_cre_class_prio) %>%
        count(sign_type) %>% mutate(frac = n / sum(n))
# Figure S6d
sign_prob_tb %>% ggbarplot(x = "human_cre_class_prio", y = "frac", position = position_stack(), fill = "sign_type", color = NA) +
        scale_fill_manual(values = sign_type_col) + coord_flip()
# GWAS overlap
gwas_prob_tb = human_dhs_tb %>% group_by(sign_type) %>%
        summarise(mean_gwas = mean(gwas_snp == "TRUE"))
# Figure 2h
gwas_prob_tb %>% ggbarplot(x = "sign_type", y = "mean_gwas", fill = "sign_type", color = NA) +
        scale_fill_manual(values = sign_type_col) +
        coord_flip() + theme(legend.position = "none")

# Figure S6: Repeat above for mouse DHS
load("metadata_processed/DHS/human_mouse_dhs_raw.rda")
# mouse DHS level analysis:
mouse_dhs_tb = tibble(id = paste0("MouseDHS-", mouse_regions$identifier))
mouse_cob_id = paste0("MouseDHS-", regions_all$mouse_id[regions_all$cob_any])
mouse_cov_id = paste0("MouseDHS-", regions_all$mouse_id[regions_all$cov_any])

sign_type_col = c("non-alignable non-conserved" = "#a3a3a3",
                  "Alignable non-conserved" = "#868586",
                  "non-alignable CO-V" = "#a8cee1",
                  "non-alignable CO-V motif" = "#2178b4",
                  "Alignable CO-B" = "#ffce0b",
                  "Alignable CO-V" = "#9d1d50")

sign_tb = readRDS("./output/nonalign_sign_combined.rds")
mouse_dhs_tb$sign_type = "non-alignable non-conserved"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% paste0("MouseDHS-", regions_all$mouse_id)] = "Alignable non-conserved"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% unique(sign_tb$mouse_id)] = "non-alignable CO-V"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% unique(sign_tb$mouse_id[sign_tb$n_motif_shared >= 1])] = "non-alignable CO-V motif"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% ghcob$mouse_no_motif_ids] = "non-alignable CO-B"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% ghcob$mouse_motif_ids] = "non-alignable CO-B motif"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% mouse_cob_id] = "Alignable CO-B"
mouse_dhs_tb$sign_type[mouse_dhs_tb$id %in% mouse_cov_id] = "Alignable CO-V"
mouse_dhs_tb$sign_type = factor(mouse_dhs_tb$sign_type, levels = sign_type_levels)
# pie chart for number of conserved called
(mouse_dhs_tb %>% count(sign_type) %>%
                mutate(pct = n / sum(n) * 100) %>%
                ggpie(x = "n", fill = "sign_type") + 
                scale_fill_manual(values = sign_type_col) +
                theme(legend.position = "right")) %>% 
        push_pdf("mouse_sign_pie", width = 8, height = 4)

print(load("./metadata_processed/DHS_cCRE/all_DHS_cCRE_v4.rda"))
mouse_cre_tb$id = paste0("MouseDHS-", mouse_cre_tb$id)
mouse_dhs_tb = left_join(mouse_dhs_tb, select(mouse_cre_tb, id, mouse_cre_class_prio))
mouse_dhs_tb$mouse_cre_class_prio[is.na(mouse_dhs_tb$mouse_cre_class_prio)] = "non-cCRE"
saveRDS(mouse_dhs_tb, file = "./output/dhs_level_tb_mouse.rds")

cre_prob_tb = mouse_dhs_tb %>% group_by(sign_type) %>%
        count(mouse_cre_class_prio) %>% mutate(frac = n / sum(n))
(cre_prob_tb %>% ggbarplot(x = "sign_type", y = "frac", position = position_stack(), fill = "mouse_cre_class_prio", color = NA) +
                scale_fill_manual(values = cCRE_v4_col) + coord_flip()) %>% 
        push_pdf("mouse_cre_prob", width = 8, height = 5)

sign_prob_tb = mouse_dhs_tb %>% group_by(mouse_cre_class_prio) %>%
        count(sign_type) %>% mutate(frac = n / sum(n))

(sign_prob_tb %>% ggbarplot(x = "mouse_cre_class_prio", y = "frac", position = position_stack(), fill = "sign_type", color = NA) +
                scale_fill_manual(values = sign_type_col) + coord_flip()) %>% 
        push_pdf("mouse_sign_prob", width = 8, height = 5)


