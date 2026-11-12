library(tidyverse)
library(Seurat)
source("./helper/helper.R")
script_output_dir = "./output/human_pig_conservation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
# core set pairs present in the human batch corrected data, normalized pig data
regions = readRDS("./intermediate_data/human_pig_mapped_regions.rds")
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds")
regions = regions[regions$identifier %in% rownames(human_correct), ]
human_correct = human_correct[regions$identifier, ]
pig_correct = readRDS("./intermediate_data/pig_chromatin_correct.rds")[regions$identifier, ]
pig_sample_info = readRDS("./intermediate_data/pig_sample_info.rds")
pig_tissue = setNames(pig_sample_info$tissue, pig_sample_info$SRR)
meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
human_tissue = setNames(meta$`Biosample term name`, meta$Accession)

# top 2500 features by combined rank of the coefficient of variation in both species
calc_cv = function(mat) {
        sqrt(apply(mat, 1, var)) / (Matrix::rowMeans(mat) + 0.01)
}
region_filter = matrixStats::rowMaxs(human_correct) > 10 & matrixStats::rowMaxs(pig_correct) > 5
combined_rank = rank(calc_cv(human_correct)) + rank(calc_cv(pig_correct))
data_cross_var_feature = rownames(human_correct)[region_filter][order(-combined_rank[region_filter])][1:2500]
saveRDS(data_cross_var_feature, file = "./intermediate_data/human_pig_chromatin_var_features.rds")

human_correct_obj = CreateSeuratObject(counts = human_correct[data_cross_var_feature, ],
                                       assay = "OpenChromatin")
human_correct_obj$organism = "human"
human_correct_obj$Accession = colnames(human_correct)
human_correct_obj = NormalizeData(human_correct_obj)
pig_correct_obj = CreateSeuratObject(counts = pig_correct[data_cross_var_feature, ],
                                     assay = "OpenChromatin")
pig_correct_obj$organism = "pig"
pig_correct_obj$tissue = pig_tissue[colnames(pig_correct_obj)]
pig_correct_obj$Accession = colnames(pig_correct_obj)
pig_correct_obj = NormalizeData(pig_correct_obj)
anchors <-
        FindIntegrationAnchors(
                list(human = human_correct_obj,
                     pig = pig_correct_obj),
                anchor.features = data_cross_var_feature,
                reduction = "cca",
                k.anchor = 5,
                k.filter = NA,
                k.score = 10,
                dims = 1:10)
integrated = IntegrateData(anchors,
                           dims = 1:10,
                           k.weight = 10)
integrated$organism = ifelse(colnames(integrated) %in% colnames(human_correct_obj), "human", "pig")
integrated$tissue = unname(c(pig_tissue, human_tissue)[colnames(integrated)])
DefaultAssay(integrated) = "integrated"
integrated = ScaleData(integrated, features = data_cross_var_feature)
integrated = RunPCA(integrated, features = data_cross_var_feature)
integrated = RunUMAP(integrated, dims = 1:10)
saveRDS(integrated, file = "./intermediate_data/human_pig_integrated.rds")

# in silico matched pairs, those with score >= 0.5 are used for CO-V and >= 0.7 for CO-B
anchor_df = make_anchor_tb(anchors, human_correct_obj, pig_correct_obj) %>%
        mutate(pig_tissue = unname(pig_tissue[accession2]),
               human_tissue = unname(human_tissue[accession1]))
write_tsv(anchor_df, paste0(script_output_dir, "anchors.tsv"))

#### tissue-system precision of matched pairs ####
# human tissue system of samples in the batch correction embedding, else from biosample keywords; cell lines are not evaluable
human_emb = read_tsv("./output/sample_matching_evaluation/batch_correction_embeddings.tsv") %>%
        filter(species == "human") %>%
        distinct(Accession, .keep_all = TRUE)
human_meta = filter(meta, Organism == "Homo sapiens") %>%
        distinct(Accession, .keep_all = TRUE)
cell_lines = c("HepG2", "K562", "GM12878", "MCF-7", "MCF7", "HuH-7", "HuH-7.5", "HCT116", "HCT-116", "A549",
               "Caco-2", "Caco2", "H1-hESC", "H1 ES", "WA09", "H9-hESC", "iPSC", "Hep G2", "T47D", "T-47D",
               "IMR-90", "IMR90", "SH-SY5Y", "PC-3", "PC3", "LNCaP", "U2OS", "HeLa", "HEK293", "OCI-LY1",
               "OCI-LY7", "Panc1", "PANC-1", "MCF 10A", "MCF10A", "Karpas-422", "DOHH2", "Jurkat", "BJ ", "BJ",
               "GM23338", "GM12891", "GM12892", "SK-N-SH", "neuroblastoma", "myeloid leukemia", "K-562", "MEL",
               "C803", "RKO", "G401", "AG09309", "AG04449", "AG09319", "AG10803", "AG04450")
is_cell_line = function(x) {
        any(sapply(tolower(cell_lines), function(cl) grepl(cl, tolower(as.character(x)), fixed = TRUE)))
}
# first matching pattern wins
keyword_system = c("muscle|myoblast|myotube|myocyte|myofibroblast" = "Muscle",
                   "fibroblast|stromal|mesenchymal|connective" = "Connective/Fibroblast",
                   "adipo|fat " = "Adipose",
                   "hepatocyte|liver" = "Liver",
                   "cardiac|heart|aorta|aortic|coronary|atrium|ventric|left atrium" = "Heart",
                   "vascular|endothelial|umbilical vein|tibial artery|aortic|pulmonary artery|hmvec" = "Cardiovascular",
                   "brain|cerebellum|cortex|astrocyte|neuron|spinal|substantia|hippocamp|caudate|putamen|midbrain|cerebrum|choroid plexus|ciliary|retin|iris|pineal|nervous system|frontal|parietal|temporal|occipital|brodmann|amygdala|thalam|nucleus accumbens" = "Brain/CNS/Eye",
                   "lung|bronchial|tracheal|alveolar|airway" = "Lung",
                   "kidney|renal|nephron|glomerul" = "Kidney",
                   "intestine|colon|stomach|duodenum|jejunum|ileum|gastric|esophagus|pancreatic duct|sigmoid|rectum" = "Digestive",
                   "pancreas|pancreatic|islet|beta cell|alpha cell" = "Pancreas",
                   "thyroid" = "Thyroid",
                   "adrenal" = "Adrenal",
                   "ovary|testis|uterus|prostate|placenta|amniotic|trophoblast|cervix|endometri|vagin|sperm|chorion|oocyte|spermatogenic|extra-embryonic|extraembryonic" = "Reproductive",
                   "erythro|leukocyte|monocyte|macrophage|t cell|t-cell|b cell|b-cell|nk cell|cd4|cd8|cd14|cd34|natural killer|dendritic|thymus|spleen|lymph|blood|granulocyte|hematopoietic|hsc|naive|memory|peyer|plasma cell" = "Immune/Blood",
                   "bone|osteo|chondrocyte|cartilage" = "Bone",
                   "epider|melanocyte|keratinocyte|skin|epithel" = "Epithelial",
                   "limb|facial|jaw|gingiva|tongue|tooth|dental" = "Limb/Facial",
                   "stem cell|esc|ipsc|induced pluripotent|h1-hesc|h9|embryonic stem" = "Embryonic/Stem")
keyword_to_system = function(term) {
        hit = which(map_lgl(names(keyword_system), grepl, x = tolower(term)))
        if (is.na(term) | length(hit) == 0) return(NA_character_)
        unname(keyword_system[hit[1]])
}
# related systems are merged before matching
group_system = function(x) {
        case_when(x %in% c("Adipose", "Connective/Fibroblast") ~ "Connective/Adipose",
                  x %in% c("Heart", "Cardiovascular") ~ "Cardiovascular",
                  x %in% c("Muscle", "Limb/Facial") ~ "Musculoskeletal",
                  x %in% c("Liver", "Pancreas", "Digestive") ~ "Digestive",
                  x %in% c("Brain/CNS/Eye", "Nervous") ~ "Nervous",
                  TRUE ~ x)
}
pig_system = c(fat = "Adipose", heart = "Heart", liver = "Liver", muscle = "Muscle", spleen = "Immune/Blood")
anchor_precision = anchor_df %>%
        mutate(human_biosample = human_meta$`Biosample term name`[match(accession1, human_meta$Accession)],
               human_summary = human_meta$`Biosample summary`[match(accession1, human_meta$Accession)],
               human_system = coalesce(human_emb$tissue_system[match(accession1, human_emb$Accession)],
                                       map_chr(human_biosample, keyword_to_system)),
               pig_system = unname(pig_system[pig_tissue]),
               is_cell_line = map_lgl(human_tissue, is_cell_line) | map_lgl(human_biosample, is_cell_line) |
                       map_lgl(human_summary, is_cell_line),
               evaluable = !is_cell_line & !is.na(human_system) & !is.na(pig_system),
               match = group_system(human_system) == group_system(pig_system))
precision_tb = map_dfr(c(0, 0.5, 0.7, 1), function(cutoff) {
        filter(anchor_precision, score >= cutoff, evaluable) %>%
                summarise(score_cutoff = cutoff, n_evaluable = n(), n_match = sum(match), precision = mean(match))
})
write_tsv(precision_tb, paste0(script_output_dir, "anchor_tissue_system_precision.tsv"))
