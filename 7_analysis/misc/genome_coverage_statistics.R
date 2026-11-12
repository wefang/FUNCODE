library(tidyverse)
library(GenomicRanges)
script_output_dir = "./output/genome_coverage/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
# all human and mouse DHSs
load("./metadata_processed/DHS/human_mouse_dhs_raw.rda")
n_human_total = length(human_regions)
n_mouse_total = length(mouse_regions)

# core set: DHSs overlapping aligned regions
load("./metadata_processed/indexed_aligned_combined.rda")
human_dhs_in_core = unique(queryHits(findOverlaps(human_regions, hg38_regions)))
mouse_dhs_in_core = unique(queryHits(findOverlaps(mouse_regions, mm10_regions)))

# extended set: DHSs matched within regulatory domains of homologous genes
matched_dhs = readRDS("./metadata_processed/human_mouse_all_DHS_matched_protein_coding.rds")
human_dhs_extended_ids = unique(unlist(map(matched_dhs$human_dhs, function(x) {
        if (!is.null(x)) as.character(x$id) else NULL
})))
mouse_dhs_extended_ids = unique(unlist(map(matched_dhs$mouse_dhs, function(x) {
        if (!is.null(x)) as.character(x$id) else NULL
})))
human_dhs_in_extended = which(as.character(human_regions$identifier) %in% human_dhs_extended_ids)
mouse_dhs_in_extended = which(as.character(mouse_regions$identifier) %in% mouse_dhs_extended_ids)

human_covered = unique(c(human_dhs_in_core, human_dhs_in_extended))
mouse_covered = unique(c(mouse_dhs_in_core, mouse_dhs_in_extended))
n_human = c(n_human_total,
            length(human_dhs_in_core),
            length(human_dhs_in_extended),
            length(intersect(human_dhs_in_core, human_dhs_in_extended)),
            length(human_covered),
            length(setdiff(1:n_human_total, human_covered)))
n_mouse = c(n_mouse_total,
            length(mouse_dhs_in_core),
            length(mouse_dhs_in_extended),
            length(intersect(mouse_dhs_in_core, mouse_dhs_in_extended)),
            length(mouse_covered),
            length(setdiff(1:n_mouse_total, mouse_covered)))
format_count = function(n, total) {
        c(format(n[1], big.mark = ","),
          sprintf("%s (%.1f%%)", sapply(n[-1], format, big.mark = ","), n[-1] / total * 100),
          sprintf("%.1f%%", n[6] / total * 100))
}
stats_table = tibble(Metric = c("Total DHSs",
                                "DHSs in core set (sequence-aligned)",
                                "DHSs in extended set (gene homology)",
                                "DHSs in both core and extended",
                                "DHSs covered (core OR extended)",
                                "DHSs NOT covered by FUNCODE",
                                "% NOT covered"),
                     Human = format_count(n_human, n_human_total),
                     Mouse = format_count(n_mouse, n_mouse_total))

# genomic span of merged DHSs relative to genome size (chr1-22/19 and chrX, excluding chrY)
hg38_chrom_sizes = read_tsv("./1_map_regions/chrom_size/hg38.chrom.sizes.txt", col_names = c("chr", "size"))
mm10_chrom_sizes = read_tsv("./1_map_regions/chrom_size/mm10.chrom.sizes.txt", col_names = c("chr", "size"))
genome_size_human_mb = sum(hg38_chrom_sizes$size[hg38_chrom_sizes$chr %in% paste0("chr", c(1:22, "X"))]) / 1e6
genome_size_mouse_mb = sum(mm10_chrom_sizes$size[mm10_chrom_sizes$chr %in% paste0("chr", c(1:19, "X"))]) / 1e6
human_dhs_span_mb = sum(width(reduce(human_regions))) / 1e6
mouse_dhs_span_mb = sum(width(reduce(mouse_regions))) / 1e6
human_covered_span_mb = sum(width(reduce(human_regions[human_covered]))) / 1e6
mouse_covered_span_mb = sum(width(reduce(mouse_regions[mouse_covered]))) / 1e6
genome_table = tibble(Metric = c("Total DHS span (Mb)",
                                 "Covered DHS span (Mb)",
                                 "Genome size (Mb)",
                                 "% genome covered by DHSs",
                                 "% genome covered by FUNCODE"),
                      Human = c(sprintf("%.1f", human_dhs_span_mb),
                                sprintf("%.1f", human_covered_span_mb),
                                format(genome_size_human_mb, big.mark = ","),
                                sprintf("%.2f%%", human_dhs_span_mb / genome_size_human_mb * 100),
                                sprintf("%.2f%%", human_covered_span_mb / genome_size_human_mb * 100)),
                      Mouse = c(sprintf("%.1f", mouse_dhs_span_mb),
                                sprintf("%.1f", mouse_covered_span_mb),
                                format(genome_size_mouse_mb, big.mark = ","),
                                sprintf("%.2f%%", mouse_dhs_span_mb / genome_size_mouse_mb * 100),
                                sprintf("%.2f%%", mouse_covered_span_mb / genome_size_mouse_mb * 100)))
write_tsv(stats_table, paste0(script_output_dir, "dhs_coverage_statistics.tsv"))
write_tsv(genome_table, paste0(script_output_dir, "genome_coverage_statistics.tsv"))
