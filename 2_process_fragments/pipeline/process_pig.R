library(tidyverse)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)

#### Zhao et al. 2021 pig ATAC-seq sample metadata (GSE143288, PRJNA597497), run before align_pig.sh ####
ena_file = "./raw_data/pig_atac_zhao/ENA_PRJNA597497_full.tsv"
download.file(paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?",
                     "accession=PRJNA597497&result=read_run&",
                     "fields=run_accession,experiment_title,library_strategy,",
                     "sample_title,scientific_name,fastq_ftp,fastq_bytes&",
                     "format=tsv&limit=0"), ena_file, mode = "w")
# original 2021 samples are SRR107646xx, later SRR293687xx runs are from a 2024 study
sample_mapping = read_tsv(ena_file) %>%
        filter(library_strategy == "ATAC-seq", grepl("^SRR107646", run_accession)) %>%
        transmute(SRR = run_accession,
                  tissue = str_replace(str_extract(experiment_title, "(?<=: )[a-z]+ tissue(?= of)"), " tissue", ""),
                  breed = trimws(str_extract(experiment_title, "(?<=tissue of )[A-Za-z ]+(?= at age)")),
                  experiment_title = experiment_title,
                  fastq_bytes_total = sapply(strsplit(fastq_bytes, ";"), function(x) sum(as.numeric(x))),
                  ftp_r1 = sapply(strsplit(fastq_ftp, ";"), `[`, 1),
                  ftp_r2 = sapply(strsplit(fastq_ftp, ";"), `[`, 2)) %>%
        mutate(size_gb = round(fastq_bytes_total / 1e9, 2)) %>%
        arrange(tissue, breed, SRR)
write_tsv(sample_mapping, "./custom_metadata/zhao_2021_sample_mapping.tsv")

#### count reads in human-pig core set regions, run after align_pig.sh ####
regions = readRDS("./intermediate_data/human_pig_mapped_regions.rds")
susScr11_regions = GRanges(regions$susScr11_extend_region)
susScr11_regions$identifier = regions$identifier
chrlen_df = read_tsv("./1_map_regions/chrom_size/susScr11.chrom.sizes.txt", col_names = c("chr", "size"))
susScr11_regions = susScr11_regions[seqnames(susScr11_regions) %in% chrlen_df$chr]
sample_mapping = read_tsv("./custom_metadata/zhao_2021_sample_mapping.tsv")
bam_files = list.files("./processed_data/pig_atac_zhao/aligned/", pattern = "\\.sorted\\.bam$", full.names = TRUE)
sample_ids = gsub("\\.sorted\\.bam$", "", basename(bam_files))
param = ScanBamParam(what = c("rname", "pos"),
                     flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE))
count_matrix = matrix(0, nrow = length(susScr11_regions), ncol = length(bam_files),
                      dimnames = list(susScr11_regions$identifier, sample_ids))
libsizes = setNames(numeric(length(bam_files)), sample_ids)
for (i in seq_along(bam_files)) {
        message(sample_ids[i])
        aln = readGAlignments(bam_files[i], param = param)
        aln = aln[seqnames(aln) %in% chrlen_df$chr]
        libsizes[i] = length(aln)
        count_matrix[, i] = countOverlaps(susScr11_regions, aln)
}
sample_info = data.frame(SRR = sample_ids) %>%
        left_join(sample_mapping[, c("SRR", "tissue", "breed")], by = "SRR")
saveRDS(count_matrix, "./intermediate_data/pig_chromatin_counts.rds")
saveRDS(libsizes, "./intermediate_data/pig_chromatin_libsize.rds")
saveRDS(sample_info, "./intermediate_data/pig_sample_info.rds")

# normalize by library size (reads per million)
saveRDS(t(t(count_matrix) / libsizes * 1e6), "./intermediate_data/pig_chromatin_correct.rds")
