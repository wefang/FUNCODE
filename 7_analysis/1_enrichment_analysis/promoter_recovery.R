library(tidyverse)
library(data.table)
source("./helper/helper.R")

# promoter-like (PLS) elements
print(load("./metadata_processed/DHS_cCRE/all_DHS_cCRE_v4.rda"))
human_pls = human_cre_tb$id[human_cre_tb$human_cre_class_prio == "PLS"]
mouse_pls = mouse_cre_tb$id[mouse_cre_tb$mouse_cre_class_prio == "PLS"]

# ortholog genes, and the genes linked to each core-set human element
sh_homo = readRDS("./output/indexed_dhs_mapped_regions_homolog_genes.rds")
hid_lookup = bind_rows(sh_homo$homolog_gene_tb) %>%
        distinct(hid, gene_human, gene_mouse) %>%
        filter(!is.na(hid)) %>%
        mutate(hid = as.character(hid))
sh_hids = tibble(human_id = sh_homo$id,
                 hids = map_chr(sh_homo$homolog_gene_filter_tb, function(tb) {
                         h = unique(tb$hid)
                         h = h[!is.na(h)]
                         if (length(h) == 0) return(NA_character_)
                         paste(h, collapse = ",")
                 })) %>%
        filter(!is.na(hids))

# number of distinct PLS-PLS pairs per gene
count_pls_pairs <- function(pair_tb, col_name) {
        pair_tb %>%
                filter(human_id %in% human_pls, mouse_id %in% mouse_pls) %>%
                separate_rows(hids, sep = ",") %>%
                distinct(hid = hids, human_id, mouse_id) %>%
                count(hid, name = col_name)
}

# core set: all aligned pairs, high CO-B pairs and high CO-V pairs
regions_all = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
sh_pairs = regions_all %>% filter(!is.na(mouse_id)) %>%
        select(human_id, mouse_id, cob_any, cov_any) %>%
        inner_join(sh_hids, by = "human_id")
sh_total = count_pls_pairs(sh_pairs, "n_SH_PLS_total")
sh_cob = count_pls_pairs(filter(sh_pairs, !is.na(cob_any) & cob_any), "n_SH_PLS_COB")
sh_cov = count_pls_pairs(filter(sh_pairs, !is.na(cov_any) & cov_any), "n_SH_PLS_COV")

# extended set GH-CO-B: union over modalities, genes shared by the two homologous gene sets of each pair
ghcob_files = c(CA = "./output/MTD-3_extended_cob/cob_significant_fdr25.rds",
                H3K27ac = "./output/MTD-3_extended_H3K27ac/significant_fdr25.rds",
                H3K4me1 = "./output/MTD-3_extended_H3K4me1/significant_fdr25.rds",
                H3K4me3 = "./output/MTD-3_extended_H3K4me3/significant_fdr25.rds")
gh_cob_pairs = map_dfr(ghcob_files, function(f) {
        as_tibble(readRDS(f)) %>% select(human_id, mouse_id, human_hid_set, mouse_hid_set)
}) %>%
        distinct(human_id, mouse_id, .keep_all = T) %>%
        mutate(hids = map2_chr(strsplit(human_hid_set, "_"), strsplit(mouse_hid_set, "_"), function(h, m) {
                hids = intersect(h, m)
                if (length(hids) == 0) return(NA_character_)
                paste(hids, collapse = ",")
        })) %>%
        filter(!is.na(hids))
gh_cob = count_pls_pairs(gh_cob_pairs, "n_GH_PLS_COB")

# extended set GH-CO-V: union over modalities
gh_cov_pairs = readRDS("./output/nonalign_sign_combined.rds") %>%
        transmute(human_id, mouse_id = sub("^MouseDHS-", "", mouse_id), hids = as.character(hid)) %>%
        distinct()
gh_cov = count_pls_pairs(gh_cov_pairs, "n_GH_PLS_COV")

# extended set, all PLS-PLS pairs of homologous gene sets sharing a gene
load("./output/all_dhs_pairs.rda")
human_pls_list = setNames(map(human_hid_set$human_dhs, ~ as.character(.x$id)[as.character(.x$id) %in% human_pls]),
                          human_hid_set$hid_set)
mouse_pls_list = setNames(map(mouse_hid_set$mouse_dhs, ~ as.character(.x$id)[as.character(.x$id) %in% mouse_pls]),
                          mouse_hid_set$hid_set)
human_hid_list = setNames(human_hid_set$hid_list, human_hid_set$hid_set)
mouse_hid_list = setNames(mouse_hid_set$hid_list, mouse_hid_set$hid_set)
gh_total = map_dfr(1:nrow(hidset_pairs), function(i) {
        human_ids = human_pls_list[[hidset_pairs$human_hid_set[i]]]
        mouse_ids = mouse_pls_list[[hidset_pairs$mouse_hid_set[i]]]
        hids = intersect(human_hid_list[[hidset_pairs$human_hid_set[i]]], mouse_hid_list[[hidset_pairs$mouse_hid_set[i]]])
        if (length(human_ids) == 0 || length(mouse_ids) == 0 || length(hids) == 0) return(NULL)
        expand_grid(hid = hids, human_id = human_ids, mouse_id = mouse_ids)
}) %>%
        distinct() %>%
        count(hid, name = "n_GH_PLS_total")

# ortholog genes in the extended set universe
gh_hids = intersect(unique(unlist(strsplit(unique(hidset_pairs$human_hid_set), "_"))),
                    unique(unlist(strsplit(unique(hidset_pairs$mouse_hid_set), "_"))))
per_gene = hid_lookup %>% filter(hid %in% gh_hids) %>%
        left_join(sh_cob, by = "hid") %>%
        left_join(gh_cob, by = "hid") %>%
        left_join(sh_total, by = "hid") %>%
        left_join(gh_total, by = "hid") %>%
        left_join(sh_cov, by = "hid") %>%
        left_join(gh_cov, by = "hid") %>%
        mutate(across(starts_with("n_"), ~ replace_na(.x, 0L))) %>%
        rename(human_gene = gene_human, mouse_gene = gene_mouse) %>%
        arrange(desc(n_GH_PLS_COB - n_SH_PLS_COB), desc(n_GH_PLS_COB))
write_tsv(per_gene, "./output/gh_promoter_recovery_per_gene.tsv")

# Figure S7f: promoter conservation category of each ortholog gene pair
per_gene %>%
        mutate(category = case_when(n_SH_PLS_COB > 0 & n_GH_PLS_COB == 0 ~ "High CO-B, core only",
                                    n_SH_PLS_COB == 0 & n_GH_PLS_COB > 0 ~ "High CO-B, extended only",
                                    n_SH_PLS_COB > 0 & n_GH_PLS_COB > 0 ~ "High CO-B, both core & extended",
                                    n_SH_PLS_COV > 0 & n_GH_PLS_COV == 0 ~ "High CO-V, core only",
                                    n_SH_PLS_COV == 0 & n_GH_PLS_COV > 0 ~ "High CO-V, extended only",
                                    n_SH_PLS_COV > 0 & n_GH_PLS_COV > 0 ~ "High CO-V, both core & extended",
                                    n_SH_PLS_total > 0 ~ "Aligned promoters, not conserved",
                                    n_GH_PLS_total > 0 ~ "Non-aligned promoters only, not conserved",
                                    TRUE ~ "No promoter pair")) %>%
        count(category) %>%
        mutate(pct = n / sum(n) * 100)
