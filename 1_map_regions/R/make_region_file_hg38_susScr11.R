# human-pig core set: forward mapped pairs plus reverse mapped pairs not duplicated in both species
library(tidyverse)
library(GenomicRanges)
map_dir = "./processed_data/indexed_dhs_hg38_susScr11_mapped/"
forward_map = read_tsv(paste0(map_dir, "hg38tosusScr11_base_alignment.txt"), col_types = cols(identifier = "c"))
reverse_map = read_tsv(paste0(map_dir, "susScr11tohg38_base_alignment.txt"), col_types = cols(identifier = "c"))

# forward_reverse identifier pairs whose regions overlap by more than 50 bp
overlap_pairs = function(forward_region, reverse_region) {
        forward_gr = GRanges(forward_region)
        reverse_gr = GRanges(reverse_region)
        hits = findOverlaps(forward_gr, reverse_gr)
        hits = hits[width(pintersect(forward_gr[queryHits(hits)], reverse_gr[subjectHits(hits)])) > 50]
        paste(forward_map$identifier[queryHits(hits)], reverse_map$identifier[subjectHits(hits)], sep = "_")
}
duplicate_pairs = intersect(overlap_pairs(forward_map$hg38_extend_region, reverse_map$hg38_extend_region),
                            overlap_pairs(forward_map$susScr11_extend_region, reverse_map$susScr11_extend_region))
duplicate_ids = unique(sapply(strsplit(duplicate_pairs, "_"), `[`, 2))

core_set = bind_rows(mutate(forward_map, source = "forward"),
                     mutate(filter(reverse_map, !identifier %in% duplicate_ids), source = "reverse")) %>%
        select(identifier, hg38_summit, susScr11_summit, hg38_extend_region, susScr11_extend_region,
               axtnet_id, axtnet_direction, gap_status, source)
write_tsv(core_set, paste0(map_dir, "hg38_susScr11_core_set.txt"))
saveRDS(core_set, file = "./intermediate_data/human_pig_mapped_regions.rds")
