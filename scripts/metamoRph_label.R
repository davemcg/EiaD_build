#remotes::install_github("davemcg/metamoRph")
library(metamoRph)
library(tidyverse)
library(data.table)
library(matrixStats)

# load metadata and filter to common eye  samplese
meta <- data.table::fread('data/eyeIntegration23_meta_2023_09_01.csv.gz') %>% 
  filter(Cohort == "Eye",
         Tissue %in% c("Conjunctiva", "Cornea","Lens", "RPE","Retina", "Trabecular Meshwork")) %>% 
  filter(!grepl("AMD", Perturbation)) %>% 
  select(sample_accession:study_title,Source_details ) %>% 
  mutate(sample_accession = gsub("^1","X1",sample_accession),
         sample_accession = gsub("-",".",sample_accession)) %>% 
  
  unique()
# load count data for eyeIntegration

gene_counts <- vroom::vroom("counts/gene_counts.csv.gz") %>% data.frame() %>% 
  column_to_rownames("Gene")

sample_frac <- 0.1
sample_min <- 5
rounds <- 20
label_guesses_list <- list()
label_val_list <- list()

training_mat <- gene_counts
training_metadata <- meta

projection_mat <- gene_counts
projection_metadata <- meta

# e.g. "CellType" or "Tissue"
sample_label <- 'Tissue'

training_sample_id <- 'sample_accession'
# column that contains the UNIQUE sample identifier in
# like SRS12341234 or "BrainSample1"
projection_sample_id <- 'sample_accession'


na_if_sample_in_training <- FALSE
for (bin in seq(1:rounds)) {
  #set.seed(sample(10))
  print(paste(bin, ' of ', rounds))
  training_metadata$sample_label <- training_metadata %>% pull(sample_label)
  # run twice, first by sample_frac
  meta_for_training_by_frac <- training_metadata %>%
    group_by(sample_label) %>%
    slice_sample(prop = sample_frac) %>%
    mutate(SCount = n()) %>% 
    unique()
  # second time by sample_n
  meta_for_training_by_n <- training_metadata %>%
    group_by(sample_label) %>%
    slice_sample(n = sample_min) %>%
    mutate(SCount = n()) %>% 
    unique()
  # combine, remove any meta_for_training_by_frac samples with 
  # a Scount < sample_min and replace with the meta_for_training_by_n
  meta_for_training <- meta_for_training_by_frac %>% filter(SCount > sample_min) 
  meta_for_training <- bind_rows(meta_for_training,
                                 meta_for_training_by_n %>% 
                                   filter(!sample_label %in% (meta_for_training %>% pull(sample_label))))
  # line up gene counts to metadata
  gene_counts_train <- training_mat[,meta_for_training %>% pull(training_sample_id)]
  
  # run pca
  pca_output <- run_pca(gene_counts_train, meta_for_training)
  # train the model
  trained_model <- model_build(pca_output$PCA$x,
                               pca_output$meta %>% pull(sample_label),
                               model = 'lm', verbose = FALSE, num_PCs = 30)
  
  # project data from the pca
  rownames(projection_mat) <- gsub(' \\(.*','',rownames(projection_mat))
  projected_data <- metamoRph((projection_mat[,projection_metadata %>% pull(projection_sample_id)]),
                                  pca_output$PCA$rotation,
                                  pca_output$center_scale)
  # apply model
  label_guesses <- model_apply(trained_model,
                               projected_data,
                               projection_metadata %>% pull(sample_label)
                               
  )
  label_mat <- model_apply(trained_model,
                           projected_data,
                           projection_metadata %>% pull(sample_label),
                           return_predictions = TRUE)
  if (na_if_sample_in_training){
    label_mat[(projection_metadata[training_sample_id] %>% pull(1)) %in% (meta_for_training[training_sample_id] %>% pull(1)),] <- NA
  }
  
  label_guesses_list[[bin]] <- label_guesses
  label_val_list[[bin]] <- label_mat
}

# combine label guesses into a single data frame
predictions <- label_guesses_list %>% bind_rows() %>% 
  group_by(sample_id) %>% summarise(miscall_count = sum(sample_label != predict),
                                    sample_label = unique(sample_label), 
                                    predict = paste(unique(predict), collapse = ', '), 
                                    mean_max_score = mean(max_score),
                                    min = min(max_score),
                                    max = max(max_score))
#pred_values <- bind_rows(label_val_list, .id = "study_accession")

# create 3D array to stack all of the model scores
model_score_3D <- array(unlist(label_val_list),
                        c(nrow(label_val_list[[1]]),
                          ncol(label_val_list[[1]]),
                          length(label_val_list)))
model_scores <- rowMeans(model_score_3D, dims = 2, na.rm = TRUE)
colnames(model_scores) <- colnames(label_val_list[[1]])
calls <- cbind(label_guesses_list[[1]] %>% pull(sample_id),
               label_guesses_list[[1]] %>% pull(sample_label),
               colnames(model_scores)[apply(model_scores, 1, which.max)],
               apply(model_scores, 1, max),
               model_scores) %>%
  data.frame()
colnames(calls)[1:4] <- c("sample_id","sample_label","predict","mean_max_score")
calls$mean_max_score <- as.numeric(calls$mean_max_score)
calls %>% filter(sample_label == predict) %>% dim()
calls %>% filter(sample_label != predict) %>% dim()

predictions <- predictions %>% select(-predict) %>% 
  left_join(calls %>% select(sample_id, predict))



write_tsv(predictions, 'data/2023_09_01_ML_tissue_predictions.tsv.gz')
write_tsv(calls, 'data/2023_09_01_ML_tissue_values.tsv.gz')
