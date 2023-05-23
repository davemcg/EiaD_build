#remotes::install_github("davemcg/eigenProjectR")
library(eigenProjectR)
library(tidyverse)
library(data.table)
library(matrixStats)

# load metadata and filter to eye only samples
meta <- read_csv("data/eyeIntegration22_meta_2023_03_03.csv.gz") %>% filter(Cohort == "Eye") %>% select(-run_accession) %>% unique()

# load count data for eyeIntegration
gene_counts <- fread("counts/gene_counts.csv.gz") %>%
  filter(!grepl("ENST", Gene)) %>% 
  mutate(Gene = gsub(" \\(.*","",Gene)) %>% 
  group_by(Gene) %>% 
  summarise_all(sum) %>% 
  column_to_rownames("Gene")

# create list of age groups
#age_groups <- c("Adult", "Fetal", "Unspecified")

# create list to store label guesses for each age group
label_guesses_list <- list()
label_val_list <- list()

# # loop through each age group
# for (age_group in age_groups) {
#   
#   # Which samples are we working with?
#   #cat("We are now processing", age_group, "samples", "\n")
#   
#   # filter metadata to the current age group
#   # meta_filtered <- if (age_group == "Unspecified") {meta %>% filter(is.na(Age))
#   # } else {meta %>% filter(Age == age_group)
#   # }
meta_filtered <- meta %>% filter(!Tissue %in% c('EyeLid', 'iPSC','ESC', 'Optic Nerve'))
gene_counts <- gene_counts[,meta_filtered$sample_accession]


# loop through excluded studies
for (study in unique(meta_filtered$study_accession)) {
  print(study)
  # filter to exclude one study
  meta_filtered_excluded <- meta_filtered %>% filter(study_accession != study)
  
  # line up gene counts to metadata
  gene_counts_filtered <- gene_counts[, meta_filtered_excluded$sample_accession]
  
  # run pca
  pca_output <- run_pca(gene_counts_filtered %>% as.matrix(), meta_filtered_excluded)
  
  # new counts to be labeled
  new_counts <- gene_counts[pca_output$PCA$rotation %>% row.names(),
                            meta_filtered %>% filter(study_accession == study) %>% 
                              pull(sample_accession), drop = FALSE] %>% as.matrix()
  
  # apply to new data (excluded study)
  row.names(new_counts) <- gsub(' \\(.*','',row.names(new_counts))
  projected_new_data <- eigenProjectR(new_counts, 
                                      pca_output$PCA$rotation)
  
  trained_model <- model_build(pca_output$PCA$x, pca_output$meta$Tissue)
  
  label_guesses <- model_apply(trained_model, projected_new_data, meta_filtered %>% filter(study_accession == study) %>% pull(Tissue))
  label_mat <- model_apply(trained_model, projected_new_data, meta_filtered %>% filter(study_accession == study) %>% pull(Tissue),
                           return_predictions = TRUE)
  
  #label_guesses_list[[paste0(age_group, "_", study)]] <- label_guesses
  label_guesses_list[[study]] <- label_guesses
  label_val_list[[study]] <- label_mat
  #cat("***", "Study", study, "has been processed for the", age_group, "age group", "***", "\n")
  
}

# }

# combine label guesses into a single data frame
predictions <- bind_rows(label_guesses_list, .id = "study_accession")
pred_values <- bind_rows(label_val_list, .id = "study_accession")

write_tsv(predictions, 'data/2023_05_18_ML_tissue_predictions.tsv.gz')
write_tsv(pred_values, 'data/2023_05_18_ML_tissue_values.tsv.gz')
