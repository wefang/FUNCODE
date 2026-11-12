## human-pig extended set: DHS pairs near 1:1 orthologous genes, CO-V and CO-B with FDR
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(qvalue)
source("./helper/helper.R")
script_output_dir = "./output/human_pig_extended/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
pig_chr = c(paste0("chr", 1:18), "chrX")
bam_files = list.files("./processed_data/pig_atac_zhao/aligned/", pattern = "\\.sorted\\.bam$", full.names = T)
pig_samples = sub("\\.sorted\\.bam$", "", basename(bam_files))

#### pig DHS: union of MACS2 peaks ####
peak_dir = paste0(script_output_dir, "peaks/")
if (!dir.exists(peak_dir)) {
        dir.create(peak_dir)
}
walk2(bam_files, pig_samples, function(f, s) {
        system(paste("macs2 callpeak -t", f, "-f BAMPE -g 2.5e9 --nomodel --shift -75 --extsize 150 -q 0.01 -n", s, "--outdir", peak_dir))
})
pig_dhs = map(paste0(peak_dir, pig_samples, "_peaks.narrowPeak"), function(f) {
        peaks = read.table(f, sep = "\t")
        peaks = peaks[peaks$V1 %in% pig_chr, ]
        GRanges(peaks$V1, IRanges(peaks$V2 + 1, peaks$V3))
})
pig_dhs = reduce(unlist(GRangesList(pig_dhs)))
pig_dhs = pig_dhs[order(as.character(seqnames(pig_dhs)), start(pig_dhs))]
pig_dhs$identifier = paste0("pig_dhs_", seq_along(pig_dhs))

#### pig signal: RPM on the DHS union ####
pig_counts = matrix(0L, length(pig_dhs), length(bam_files), dimnames = list(pig_dhs$identifier, pig_samples))
lib_size = numeric(length(bam_files))
for (i in seq_along(bam_files)) {
        aln = readGAlignments(bam_files[i], param = ScanBamParam(what = c("rname", "pos"),
                                                                 flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F)))
        aln = aln[as.character(seqnames(aln)) %in% pig_chr]
        lib_size[i] = length(aln)
        pig_counts[, i] = countOverlaps(pig_dhs, aln)
}
pig_union_correct = t(t(pig_counts) / lib_size * 1e6)

#### 1:1 protein-coding orthologs (Ensembl 115 BioMart) ####
mart_query = function(attrs) {
        xml = paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
                     '<Query virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6">',
                     '<Dataset name = "hsapiens_gene_ensembl" interface = "default">',
                     paste0('<Attribute name = "', attrs, '" />', collapse = ""), '</Dataset></Query>')
        read_tsv(paste0("https://sep2025.archive.ensembl.org/biomart/martservice?query=", URLencode(xml, reserved = T)),
                 col_types = cols(.default = "c"), quote = "", na = "")
}
orth = merge(mart_query(c("ensembl_gene_id", "sscrofa_homolog_ensembl_gene", "sscrofa_homolog_associated_gene_name",
                          "sscrofa_homolog_orthology_type", "sscrofa_homolog_orthology_confidence",
                          "sscrofa_homolog_perc_id", "sscrofa_homolog_perc_id_r1")),
             mart_query(c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "chromosome_name")),
             by = "Gene stable ID", all.x = T)
colnames(orth) = c("ensembl_gene_id", "sscrofa_homolog_ensembl_gene", "sscrofa_homolog_associated_gene_name",
                   "sscrofa_homolog_orthology_type", "sscrofa_homolog_orthology_confidence",
                   "sscrofa_homolog_perc_id", "sscrofa_homolog_perc_id_r1", "hgnc_symbol", "gene_biotype", "chromosome_name")
write_tsv(orth, paste0(script_output_dir, "human_pig_orthologs.tsv"))
orth = filter(orth, !is.na(sscrofa_homolog_ensembl_gene),
              sscrofa_homolog_orthology_type == "ortholog_one2one",
              gene_biotype == "protein_coding")
orth$hid = as.character(1:nrow(orth))

#### candidate pairs: DHSs in regulatory domains of the same ortholog ####
make_domain = function(genes, basal_up = 5000, basal_down = 1000, domain_max = 1e6) {
        genes_basal = restrict(promoters(genes, upstream = basal_up, downstream = basal_down), start = 0)
        pre = precede(genes, genes_basal, ignore.strand = T)
        fol = follow(genes, genes_basal, ignore.strand = T)
        pre_dist = rep(domain_max, length(pre))
        pre_dist[!is.na(pre)] = distance(genes[!is.na(pre)], genes_basal[pre[!is.na(pre)]], ignore.strand = T)
        pre_dist = pmin(pre_dist, domain_max)
        fol_dist = rep(domain_max, length(fol))
        fol_dist[!is.na(fol)] = distance(genes[!is.na(fol)], genes_basal[fol[!is.na(fol)]], ignore.strand = T)
        fol_dist = pmin(fol_dist, domain_max)
        genes_domain = resize(genes, width(genes) + pre_dist, fix = "start", ignore.strand = T)
        genes_domain = resize(genes_domain, width(genes_domain) + fol_dist, fix = "end", ignore.strand = T)
        genes_domain = punion(restrict(genes_domain, start = 0), genes_basal)
        genes_domain$gene_id = genes$gene_id
        genes_domain
}
pig_gtf = import("./metadata/Sus_scrofa.Sscrofa11.1.115.gtf.gz")
pig_genes_domain = make_domain(subset(pig_gtf, type == "gene" & gene_biotype == "protein_coding"))
human_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_genes.rds")
human_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_gene_domains.rds")
human_genes_domain$gene_id = setNames(strip_version(human_genes_gr$gene_id), human_genes_gr$gene_name)[human_genes_domain$gene_name]
load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda")
pig_dhs_ensembl = renameSeqlevels(pig_dhs, sub("^chr", "", seqlevels(pig_dhs)))

human_genes_domain = human_genes_domain[!is.na(human_genes_domain$gene_id) & human_genes_domain$gene_id %in% orth$ensembl_gene_id]
pig_genes_domain = pig_genes_domain[pig_genes_domain$gene_id %in% orth$sscrofa_homolog_ensembl_gene]
human_ol = findOverlaps(human_regions, human_genes_domain)
pig_ol = findOverlaps(pig_dhs_ensembl, pig_genes_domain)
human_by_hid = split(human_regions$identifier[queryHits(human_ol)],
                     orth$hid[match(human_genes_domain$gene_id[subjectHits(human_ol)], orth$ensembl_gene_id)])
pig_by_hid = split(pig_dhs$identifier[queryHits(pig_ol)],
                   orth$hid[match(pig_genes_domain$gene_id[subjectHits(pig_ol)], orth$sscrofa_homolog_ensembl_gene)])
pairs = rbindlist(map(intersect(names(human_by_hid), names(pig_by_hid)), function(h) {
        data.table(human_dhs_id = rep(human_by_hid[[h]], each = length(pig_by_hid[[h]])),
                   pig_dhs_id = rep(pig_by_hid[[h]], times = length(human_by_hid[[h]])),
                   hid = h)
}))
pairs = merge(pairs, as.data.table(orth)[, .(hid, human_gene_ensembl = ensembl_gene_id, pig_gene_ensembl = sscrofa_homolog_ensembl_gene)],
              by = "hid", all.x = T)
pairs = merge(pairs, data.table(human_dhs_id = human_regions$identifier, human_start = start(human_regions), human_end = end(human_regions)),
              by = "human_dhs_id", all.x = T)
pairs = merge(pairs, data.table(pig_dhs_id = pig_dhs$identifier, pig_start = start(pig_dhs), pig_end = end(pig_dhs)),
              by = "pig_dhs_id", all.x = T)
pairs = unique(pairs, by = c("human_dhs_id", "pig_dhs_id"))

#### distance group of each pair, to the TSS of its ortholog ####
gene_tss = function(gtf) {
        genes = gtf[gtf$type == "gene"]
        unique(data.table(gene_id = strip_version(genes$gene_id),
                          tss = ifelse(as.character(strand(genes)) == "-", end(genes), start(genes))), by = "gene_id")
}
human_tss = gene_tss(import("./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"))
pig_tss = gene_tss(pig_gtf)
dist_group = function(d) ifelse(d <= 200, 1L, ifelse(d <= 2000, 2L, 3L))
pairs[, human_dist := abs((human_start + human_end) %/% 2L - human_tss[.(human_gene_ensembl), on = "gene_id", x.tss])]
pairs[, pig_dist := abs((pig_start + pig_end) %/% 2L - pig_tss[.(pig_gene_ensembl), on = "gene_id", x.tss])]
pairs[, `:=`(human_dist_group = dist_group(human_dist), pig_dist_group = dist_group(pig_dist))]
pairs[, cov_group := 10L * human_dist_group + pig_dist_group]
# DHS level group from the closest TSS among its pairs, for the null
human_dhs_group = pairs[!is.na(human_dist), .(min_dist = min(human_dist)), by = human_dhs_id]
human_dhs_group = setNames(dist_group(human_dhs_group$min_dist), human_dhs_group$human_dhs_id)
pig_dhs_group = pairs[!is.na(pig_dist), .(min_dist = min(pig_dist)), by = pig_dhs_id]
pig_dhs_group = setNames(dist_group(pig_dhs_group$min_dist), pig_dhs_group$pig_dhs_id)
n_orthologs = uniqueN(pairs, by = c("human_gene_ensembl", "pig_gene_ensembl"))

#### per-DHS h-value, activity and scaled ranks over anchors with score >= 0.7 ####
anchors = distinct(read_tsv("./output/human_pig_conservation/anchors.tsv")) %>% filter(score >= 0.7)
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds")
anchors = filter(anchors, accession1 %in% colnames(human_correct), accession2 %in% colnames(pig_union_correct))
human_anchor = human_correct[intersect(unique(pairs$human_dhs_id), rownames(human_correct)), anchors$accession1]
rm(human_correct)
pig_anchor = as.matrix(pig_union_correct[intersect(unique(pairs$pig_dhs_id), rownames(pig_union_correct)), anchors$accession2])
prob_vec = seq(0.01, 1.0, by = 0.01)
compute_hval_single = function(mat, quantile_vec) {
        hval = numeric(nrow(mat))
        for (i in seq_along(prob_vec)) {
                hval[rowMeans(mat > quantile_vec[i]) >= prob_vec[i]] = prob_vec[i]
        }
        hval
}
# centered unit-norm ranks, so Spearman is a dot product; NA for constant rows
rank_scale = function(mat) {
        mat_rank = matrixStats::rowRanks(mat, ties.method = "average")
        mat_rank = mat_rank - rowMeans(mat_rank)
        ss = rowSums(mat_rank^2)
        mat_rank = mat_rank / ifelse(ss == 0, 1, sqrt(ss))
        mat_rank[ss == 0, ] = NA
        mat_rank
}
set.seed(42)
human_quantile_vec = quantile(human_anchor[sample(nrow(human_anchor), 1e5), ], prob_vec)
human_hval = compute_hval_single(human_anchor, human_quantile_vec)
n_h_active = as.integer(rowSums(human_anchor > human_quantile_vec["80%"]))
human_rank = rank_scale(human_anchor)
set.seed(42)
pig_quantile_vec = quantile(pig_anchor[sample(nrow(pig_anchor), 1e5), ], prob_vec)
pig_hval = compute_hval_single(pig_anchor, pig_quantile_vec)
n_p_active = as.integer(rowSums(pig_anchor > pig_quantile_vec["80%"]))
pig_rank = rank_scale(pig_anchor)
human_ids = rownames(human_anchor)
pig_ids = rownames(pig_anchor)
human_group = human_dhs_group[human_ids]
pig_group = pig_dhs_group[pig_ids]
rm(human_anchor, pig_anchor)

#### pair scores ####
human_index = match(pairs$human_dhs_id, human_ids)
pig_index = match(pairs$pig_dhs_id, pig_ids)
pair_cor = function(i, j) {
        out = map(seq(1, length(i), by = 5e5), function(s) {
                b = s:min(s + 5e5 - 1, length(i))
                rowSums(human_rank[i[b], , drop = F] * pig_rank[j[b], , drop = F])
        })
        pmin(pmax(unlist(out), -1), 1)
}
scores = pairs[, .(human_dhs_id, pig_dhs_id,
                   n_h_active = n_h_active[human_index],
                   n_p_active = n_p_active[pig_index],
                   cov_group, human_dist_group, pig_dist_group)]
rm(pairs)
scores[, cov_eligible := !is.na(n_h_active) & !is.na(n_p_active) & n_h_active >= 5 & n_p_active >= 5]
scores[, eligible := !is.na(n_h_active) & !is.na(n_p_active) & n_h_active >= 50 & n_p_active >= 50]
scores[, cov := pair_cor(human_index, pig_index)]
scores[, cob := pmin(human_hval[human_index], pig_hval[pig_index])]

#### CO-V: null in each of 9 distance categories, joint FDR ####
scores[, cov_pval := NA_real_]
for (h_g in 1:3) {
        for (p_g in 1:3) {
                set.seed(42 + 100 * h_g + p_g)
                null_indices1 = sample(which(human_group == h_g & n_h_active >= 5), 1e5, replace = T)
                null_indices2 = sample(which(pig_group == p_g & n_p_active >= 5), 1e5, replace = T)
                cor_null = pair_cor(null_indices1, null_indices2)
                cat_indices = scores[, which(cov_group == 10 * h_g + p_g & cov_eligible & !is.na(cov))]
                set(scores, cat_indices, "cov_pval", empPvals(scores$cov[cat_indices], cor_null[!is.na(cor_null)]))
        }
}
scores[, cov_fdr := p.adjust(cov_pval, method = "fdr")]

#### CO-B: pooled null, same activity filter, zero-clip ####
set.seed(42)
null_indices1 = sample.int(length(human_hval), 1e6, replace = T)
null_indices2 = sample.int(length(pig_hval), 1e6, replace = T)
hval_null = pmin(human_hval[null_indices1], pig_hval[null_indices2])
hval_null = hval_null[n_h_active[null_indices1] >= 50 & n_p_active[null_indices2] >= 50]
scores[, cob_pval := NA_real_]
scores[eligible == T, cob_pval := empPvals(cob, hval_null)]
scores[eligible == T & cob == 0, cob_pval := 1]
scores[, cob_fdr := p.adjust(cob_pval, method = "fdr")]

setcolorder(scores, c("human_dhs_id", "pig_dhs_id", "n_h_active", "n_p_active", "cov_eligible", "eligible",
                      "cov_group", "human_dist_group", "pig_dist_group", "cov", "cov_pval", "cov_fdr", "cob", "cob_pval", "cob_fdr"))
saveRDS(scores, paste0(script_output_dir, "hp_extended_scores_with_fdr.rds"))
fwrite(scores[eligible == T | (!is.na(cov_fdr) & cov_fdr <= 0.25)], paste0(script_output_dir, "hp_extended_scores_with_fdr.tsv"), sep = "\t", na = "NA")
write_tsv(tibble(metric = c("Total pairs", "Orthologs", "CO-V FDR <= 0.10", "CO-V FDR <= 0.25", "CO-B FDR <= 0.10", "CO-B FDR <= 0.25"),
                 value = c(nrow(scores), n_orthologs,
                           sum(scores$cov_fdr <= 0.1, na.rm = T), sum(scores$cov_fdr <= 0.25, na.rm = T),
                           sum(scores$cob_fdr <= 0.1, na.rm = T), sum(scores$cob_fdr <= 0.25, na.rm = T))),
          paste0(script_output_dir, "hp_extended_sig_counts.tsv"))
