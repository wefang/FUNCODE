library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(patchwork)
source("./helper/helper.R")
# load human and mouse gene regulatory domains
human_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_gene_domains.rds")
human_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_genes.rds")
mouse_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_gene_domains.rds")
mouse_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_genes.rds")

assign_gene <- function(regions, genes_domain, genes_gr) {
        gene_promoters = promoters(genes_gr, 0, 1)
        domain_ol = findOverlaps(regions, genes_domain)
        tb = data.frame(id = regions$identifier[queryHits(domain_ol)],
                        gene = as.character(genes_domain$gene_name[subjectHits(domain_ol)]),
                        dist_to_tss = distance(regions[queryHits(domain_ol)],
                                               gene_promoters[subjectHits(domain_ol)],
                                               ignore.strand = T),
                        stringsAsFactors = F)
        tb
}
print(load("metadata_processed/indexed_aligned_combined.rda"))
mouse_tb = assign_gene(mm10_regions, mouse_genes_domain, mouse_genes_gr)
human_tb = assign_gene(hg38_regions, human_genes_domain, human_genes_gr)

# homologous genes
source("./load_homolog.R")
common_indices = which((human_symbols %in% human_genes_gr$gene_name) &
                               (mouse_symbols %in% mouse_genes_gr$gene_name))
homolog_id = homolog_id[common_indices]
human_symbols = human_symbols[common_indices]
mouse_symbols = mouse_symbols[common_indices]

# below if only unique gene names/id
# human_gene_name_count = table(human_genes_gr$gene_name)
# human_gene_name_unqiue = names(human_gene_name_count)[human_gene_name_count == 1]
# mouse_gene_name_count = table(mouse_genes_gr$gene_name)
# mouse_gene_name_unqiue = names(mouse_gene_name_count)[mouse_gene_name_count == 1]
# end homologous genes

human_tb$hid = homolog_id[match(human_tb$gene, human_symbols)]
mouse_tb$hid = homolog_id[match(mouse_tb$gene, mouse_symbols)]
human_tb_nest = nest(human_tb, human_gene_tb = -id)
mouse_tb_nest = nest(mouse_tb, mouse_gene_tb = -id)

regions = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
regions$id = hg38_regions$identifier
# need to remove older results
regions = regions[!names(regions) %in% c("human_gene_tb", "mouse_gene_tb")]
regions = left_join(regions, human_tb_nest)
regions = left_join(regions, mouse_tb_nest)

# mouse id need to be matched manually
regions$n_mouse = map_dbl(regions$mouse_gene_tb, function(x) {
        if (is.null(x)) {
                return(0)
        } else {
                return(nrow(x))
        }
})
regions$n_human = map_dbl(regions$human_gene_tb, function(x) {
        if (is.null(x)) {
                return(0)
        } else {
                return(nrow(x))
        }
})
# table(regions$n_human > 0, regions$n_mouse > 0)
regions$hid_count = map2_dbl(regions$human_gene_tb, regions$mouse_gene_tb, function(h_tb, m_tb) {
        if (!is.null(h_tb) & !is.null(m_tb)) {
                return(sum(h_tb$hid %in% m_tb$hid))
        } else {
                return(0)
        }
})
regions$consist_target = regions$hid_count > 0
g1 = regions %>% ggviolin(x = "consist_target", y = "cov_chromatin_accessibility",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-CA")
        # stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_chromatin_accessibility",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-CA")
        # stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_ca", width = 4., h = 3.5)

g1 = regions %>% ggviolin(x = "consist_target", y = "cov_H3K4me3",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-H3K4me3")
        # stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_H3K4me3",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-H3K4me3")
        # stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_H3K4me3", width = 4., h = 3.5)

g1 = regions %>% ggviolin(x = "consist_target", y = "cov_H3K27ac",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-H3K27ac")
# stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_H3K27ac",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-H3K27ac")
# stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_H3K27ac", width = 4., h = 3.5)

g1 = regions %>% ggviolin(x = "consist_target", y = "cov_H3K4me1",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-H3K4me1")
# stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_H3K4me1",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-H3K4me1")
# stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_H3K4me1", width = 4., h = 3.5)


regions_gene = regions[which(regions$hid_count >=1 & regions$hid_count <= 2), ]
j = 0
regions_gene$homolog_gene_tb = map2(regions_gene$human_gene_tb, regions_gene$mouse_gene_tb, function(h_tb, m_tb) {
        j <<- j + 1
        if (j %% 100000 == 0) {
                message(j)
        }
        if (!is.null(h_tb) & !is.null(m_tb)) {
                if (sum(h_tb$hid %in% m_tb$hid) > 0) {
                        h_tb = filter(h_tb, !is.na(hid))
                        b_tb = filter(m_tb, !is.na(hid))
                        hid_tb = inner_join(h_tb, m_tb, by = "hid", suffix = c("_human", "_mouse"))
                        return(hid_tb)
                } else {
                        return(NULL)
                }
        } else {
                return(NULL)
        }
})
saveRDS(regions_gene, file = "./intermediate_data/regions_gene_v1.rds")


regions_gene = readRDS("./intermediate_data/regions_gene_v1.rds")
# Method 2: select overall closest gene
j = 0
regions_gene$homolog_gene_filter_tb = map(regions_gene$homolog_gene_tb, function(tb) {
        j <<- j + 1
        if (j %% 1e4 == 0) {
                message(j)
        }
        tb = mutate(tb,
                    dist_to_tss_human_rank = rank(dist_to_tss_human),
                    dist_to_tss_mouse_rank = rank(dist_to_tss_mouse),
                    dist_to_tss_rank = dist_to_tss_human_rank + dist_to_tss_mouse_rank)
        tb = filter(tb, dist_to_tss_rank == min(dist_to_tss_rank))
        tb  
})
regions_gene_ext = select(regions_gene, id,
                          human_cre_class_prio_v4,
                          mouse_cre_class_prio_v4,
                          cov_chromatin_accessibility,
                          cov_H3K4me1,
                          cov_H3K4me3,
                          cov_H3K27ac,
                          cob_chromatin_accessibility,
                          cob_H3K4me1,
                          cob_H3K4me3,
                          cob_H3K27ac,
                          percentage, phastCons4way, phyloP4way,
                          homolog_gene_filter_tb) %>%
        unnest(cols = "homolog_gene_filter_tb")

gene_tb = readRDS("./output/rna_conservation/rna_conservation_scores.rds")
regions_gene_ext = left_join(regions_gene_ext, gene_tb)

# add within-species baselines, LECIF and AdaLiftOver scores
baseline_tb = purrr::reduce(list(readRDS("./output/chromatin_conservation_v1/CA_species_baselines.rds"),
                          readRDS("./output/histone_conservation/H3K27ac_species_baselines.rds"),
                          readRDS("./output/histone_conservation/H3K4me1_species_baselines.rds"),
                          readRDS("./output/histone_conservation/H3K4me3_species_baselines.rds")),
                     left_join, by = "id")
regions_v1 = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
adaliftover_tb = readRDS("./output/adaliftover/adaliftover_scores.rds")
adaliftover_tb = tibble(id = regions_v1$id[adaliftover_tb$human_idx],
                        adaliftover_h2m = adaliftover_tb$AdaLiftOver_hg38_to_mm10,
                        adaliftover_m2h = adaliftover_tb$AdaLiftOver_mm10_to_hg38)
regions_gene_ext = regions_gene_ext %>%
        left_join(baseline_tb, by = "id") %>%
        left_join(select(regions_v1, id, lecif), by = "id") %>%
        left_join(adaliftover_tb, by = "id")
score_chr = c("cov_chromatin_accessibility", "cov_H3K27ac", "cov_H3K4me1", "cov_H3K4me3",
              "cob_chromatin_accessibility", "cob_H3K27ac", "cob_H3K4me1", "cob_H3K4me3",
              "percentage", "phastCons4way", "phyloP4way", "lecif",
              setdiff(names(baseline_tb), "id"),
              "adaliftover_h2m", "adaliftover_m2h")

# correlation of RE scores with GE CO-V (gene_spearman) and GE CO-B (gene_hval) by human and mouse cCRE category
produce_gene_cor <- function(score, ge_score, cre_col) {
        regions_gene_ext %>%
                filter(!is.na(.data[[score]]) & !is.na(.data[[ge_score]])) %>%
                group_by(cre_category = .data[[cre_col]]) %>%
                summarise(dhs_gene_cor = cor(.data[[score]], .data[[ge_score]], method = "spearman"),
                          n = n(),
                          .groups = "drop") %>%
                mutate(Method = score, ge_score = ge_score, cre_col = cre_col)
}
gene_cor_tb = bind_rows(map(c("human_cre_class_prio_v4", "mouse_cre_class_prio_v4"), function(cre_col) {
        bind_rows(map(c("gene_spearman", "gene_hval"), function(ge_score) {
                bind_rows(map(score_chr, produce_gene_cor, ge_score = ge_score, cre_col = cre_col))
        }))
}))
write_tsv(gene_cor_tb, "./output/rna_conservation/re_ge_cor_by_cre.tsv")

(gene_cor_tb %>%
        filter(cre_category %in% c("PLS", "pELS", "dELS", "non-cCRE")) %>%
        mutate(cre_category = factor(cre_category, levels = c("PLS", "pELS", "dELS", "non-cCRE")),
               Method = factor(Method, levels = score_chr)) %>%
        ggbarplot(x = "cre_category",
                  y = "dhs_gene_cor",
                  fill = "Method",
                  color = NA,
                  position = position_dodge(),
                  facet.by = c("cre_col", "ge_score")) +
        theme(legend.position = "right",
              text = element_text(size = 8),
              legend.key.size = unit(0.25, 'cm'))) %>%
        push_pdf(file_name = "domain_dhs_gene_cor_by_cre", width = 12, h = 6.)


