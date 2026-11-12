library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(AdaLiftOver)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
# AdaLiftOver and motifmatchr are installed from GitHub (ThomasDCY/AdaLiftOver, GreenleafLab/motifmatchr)
script_output_dir = "./output/adaliftover/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
load("./metadata_processed/indexed_aligned_combined.rda")
data("epigenome_hg38")
data("epigenome_mm10")
data("jaspar_pfm_list")
chunk_size = 1e4

compute_score <- function(epi, gram) {
        1 / (1 + exp(-(-3 + 4 * epi + 5 * gram)))
}
# map one chunk of query regions and keep the best scoring candidate per query
run_chunk <- function(query_gr, chain, epigenome_query, epigenome_target, genome_query, genome_target) {
        gr_target_list = adaptive_liftover(query_gr, chain, window = 2000, step_size = 200)
        gr_target_list = compute_similarity_epigenome(query_gr, gr_target_list, epigenome_query, epigenome_target)
        grammar_ok = FALSE
        tryCatch({
                gr_target_list = compute_similarity_grammar(query_gr, gr_target_list,
                                                            genome_query, genome_target, jaspar_pfm_list)
                grammar_ok = TRUE
        }, error = function(e) {
                message("grammar similarity failed: ", conditionMessage(e))
        })
        bind_rows(map(seq_along(gr_target_list), function(i) {
                candidates = gr_target_list[[i]]
                if (length(candidates) == 0) return(NULL)
                has_epigenome = "epigenome" %in% names(mcols(candidates))
                has_grammar = grammar_ok && "grammar" %in% names(mcols(candidates))
                if (!has_epigenome && !has_grammar) return(NULL)
                if (has_epigenome && has_grammar) {
                        all_scores = compute_score(mcols(candidates)$epigenome, mcols(candidates)$grammar)
                } else if (has_epigenome) {
                        all_scores = compute_score(mcols(candidates)$epigenome, 0)
                } else {
                        all_scores = compute_score(0, mcols(candidates)$grammar)
                }
                best_idx = which.max(all_scores)
                if (length(best_idx) == 0 || is.na(all_scores[best_idx])) return(NULL)
                tibble(query_idx = i,
                       target_chr = as.character(seqnames(candidates[best_idx])),
                       target_start = start(candidates[best_idx]),
                       target_end = end(candidates[best_idx]),
                       score = all_scores[best_idx])
        }))
}

direction_tb = tibble(direction = c("hg38_to_mm10", "mm10_to_hg38"),
                      chain_file = c("hg38.mm10.rbest.chain", "mm10.hg38.rbest.chain"))
complete_list = map(1:nrow(direction_tb), function(k) {
        if (direction_tb$direction[k] == "hg38_to_mm10") {
                query_regions = hg38_regions
                target_regions = mm10_regions
                epigenome_query = epigenome_hg38
                epigenome_target = epigenome_mm10
                genome_query = BSgenome.Hsapiens.UCSC.hg38
                genome_target = BSgenome.Mmusculus.UCSC.mm10
        } else {
                query_regions = mm10_regions
                target_regions = hg38_regions
                epigenome_query = epigenome_mm10
                epigenome_target = epigenome_hg38
                genome_query = BSgenome.Mmusculus.UCSC.mm10
                genome_target = BSgenome.Hsapiens.UCSC.hg38
        }
        chain = import.chain(paste0("./metadata_processed/chain_files/", direction_tb$chain_file[k]))
        chunk_starts = seq(1, length(query_regions), by = chunk_size)
        results_tb = bind_rows(map(chunk_starts, function(start_idx) {
                message(direction_tb$direction[k], " ", start_idx)
                end_idx = min(start_idx + chunk_size - 1, length(query_regions))
                run_chunk(query_regions[start_idx:end_idx], chain,
                          epigenome_query, epigenome_target, genome_query, genome_target) %>%
                        mutate(query_idx = start_idx - 1 + query_idx)
        })) %>%
                arrange(query_idx)
        # assign each mapped target to the nearest region of the target species
        target_centers = GRanges(results_tb$target_chr,
                                 IRanges(start = (results_tb$target_start + results_tb$target_end) %/% 2, width = 1))
        results_tb$target_idx = nearest(target_centers, target_regions)
        # regions without an AdaLiftOver mapping get score 0
        complete_tb = tibble(query_idx = 1:length(query_regions),
                             target_idx = 1:length(query_regions),
                             score = 0)
        complete_tb$target_idx[results_tb$query_idx] = results_tb$target_idx
        complete_tb$score[results_tb$query_idx] = results_tb$score
        complete_tb
})
names(complete_list) = direction_tb$direction

score_h2m = complete_list$hg38_to_mm10$score
score_m2h = complete_list$mm10_to_hg38$score
adaliftover_tb = tibble(human_idx = complete_list$hg38_to_mm10$query_idx,
                        mouse_idx = complete_list$hg38_to_mm10$target_idx,
                        AdaLiftOver_hg38_to_mm10 = score_h2m,
                        AdaLiftOver_mm10_to_hg38 = score_m2h,
                        AdaLiftOver_bidirectional = ifelse(score_h2m > 0 & score_m2h > 0,
                                                           (score_h2m + score_m2h) / 2,
                                                           pmax(score_h2m, score_m2h)))
saveRDS(adaliftover_tb, paste0(script_output_dir, "adaliftover_scores.rds"))
