This step computes CO-V scores for candidate human/mouse DHS pairs near homologous genes and calls significant pairs by comparing to a null background.

- 1-1annotate_dhs_pairs.R: compile and annotate human/mouse DHS pairs
- 1-2preprocess_data.R: retain matched samples with anchor score > 0.5 and obtain the human/mouse normalized matrix for conservation score calculation
- 1-3split_annot.R: split the annotation files into smaller chunks
- 1-4get_coordinate.R: get coordinates
- 2compute_stat.R: calculate conservation scores for human-mouse DHS pairs
- 3null.R: construct a null distribution for each DHS pair category (9)
- 4-1compute_pvalue.R: calculate p-value
- 4-2compute_fdr_joint.R: combine all 9 categories and calculate FDR
- 5result_bygene.R: split the results by gene
- 6compute_cob.R: calculate CO-B scores for human/mouse DHS pairs of one data modality
- cnsv_score.R: core functions for calculating conservation scores

