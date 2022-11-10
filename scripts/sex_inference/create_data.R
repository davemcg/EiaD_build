### All data files necessary for analysis can be found here (files must be uncompressed prior to use):
# https://hpc.nih.gov/~parikhpp/sex_inference/insert_file_name_here.csv.gz

# Set working directory
setwd("/Users/parikhpp/git/sex_inference")
# Import libraries
library(tidyverse)
# Import Data
eyeIntegration22 <-
  read.csv("data/eyeIntegration22_meta_2022_10_25.csv")
gender_data <- read.csv("data/gender_metadata.csv") # Downloaded from the SRA
old_ei <- read.csv("https://hpc.nih.gov/~parikhpp/EiaD/2022_metadata.csv") # Contained in the parikhpp EiaD datashare directory

# Sex inference metadata -----

# Only keep samples in EyeIntegration
eyeIntegration22_gtex <- eyeIntegration22 %>% filter(study_accession == "SRP012682")
gender_data_gtex <- gender_data %>% filter(SRA.Study == "SRP012682")

eyeIntegration22_not_gtex <- eyeIntegration22 %>% filter(study_accession != "SRP012682")
gender_data_not_gtex <- gender_data %>% filter(SRA.Study != "SRP012682")

gtex_keep <- intersect(eyeIntegration22_gtex$gtex_accession, gender_data_gtex$Sample.Name)
gtex_data <- gender_data_gtex %>% filter(Sample.Name %in% gtex_keep)

not_gtex_keep <- intersect(eyeIntegration22_not_gtex$run_accession, gender_data_not_gtex$Run)
not_gtex_data <- gender_data_not_gtex %>% filter(Run %in% not_gtex_keep)

gender_metadata <- bind_rows(not_gtex_data, gtex_data)

# Subset for ocular samples where gender is given

sex_inf <- gender_metadata %>% 
  filter(SRA.Study!="SRP012682") %>% 
  filter(sex != "") %>% 
  filter(sex != "not applicable") %>% 
  filter(sex != "not determined")
sex_inf <- rename(sex_inf, run_accession = Run)

# Combine data

sex_inf_meta <- old_ei %>% left_join(sex_inf, by = "run_accession")
sex_inf_meta <- sex_inf_meta %>% filter(!(is.na(sex)))

# Create inference metadata

write.csv(sex_inf_meta, "data/sex_inference_metadata.csv")

# Regression data -----

# XIST ENSG00000229807 -- y chromosome
# EIF1AY ENSG00000198692 -- y chromosome
# KDM5D ENSG00000012817 -- y chromosome
# UTY ENSG00000183878 -- y chromosome
# DDX3Y ENSG00000067048 -- y chromosome
# RPS4Y1 ENSG00000129824 -- y chromosome

#Imported from the create_sex_inference_metadata.R script output
sex_inf_metadata <- 
  read.csv("data/sex_inference_metadata.csv")

# Download gene counts files from the gene_counts folder in the EiaD_build repository
gene_counts_recount3 <- 
  vroom::vroom("data/recount3_transformed_counts.csv")
gene_counts_local_additions <- 
  vroom::vroom("data/local_data_additions_transformed_counts.csv")

# Combine recount3 and local samples for sex inference
combined_gene_counts <- gene_counts_local_additions %>% left_join(gene_counts_recount3, by = "gene_id")

genes <- c("ENSG00000229807.11", "ENSG00000198692.9", "ENSG00000012817.15", "ENSG00000183878.15",
           "ENSG00000067048.16", "ENSG00000129824.15")

combined_gene_counts_subset <- combined_gene_counts %>% filter(gene_id %in% genes)

# Log normalization
log_normalized_counts <- combined_gene_counts_subset %>% select(-gene_id) %>% +1 %>%  log()

# Only use recount3 samples included in eyeIntegration
log_normalized_counts_keep <- intersect(names(log_normalized_counts), sex_inf_metadata$run_accession)
combined_gene_counts_subset <- combined_gene_counts_subset %>% 
  select(one_of(c("gene_id", log_normalized_counts_keep)))
#Only use samples which pass quality check
qc_keep <- old_ei %>% filter(Kept == "Kept") %>% pull(run_accession)
combined_gene_counts_subset_quality <- combined_gene_counts_subset %>% select(one_of(c("gene_id", qc_keep)))

#Transform data long-ways

gene_counts_long <- t(combined_gene_counts_subset_quality) %>% as.data.frame()
gene_counts_long <- gene_counts_long %>% row_to_names(row_number = 1)
gene_counts_long <- tibble::rownames_to_column(gene_counts_long, "run_accession")

# Add sex data

regression_data <- 
  gene_counts_long %>% left_join(sex_inf_metadata %>% select(run_accession, sex), by = "run_accession")
colnames(regression_data) = c("run_accession", "XIST","EIF1AY","KDM5D","UTY","DDX3Y","RPS4Y1","sex")

# Add mapping data and x/y alignment percentages

recount3_mapping_information <- vroom::vroom("data/recount3_mapping_information.csv")
local_data_additions_mapping_information <- vroom::vroom("data/local_data_additions_mapping_information.csv")

mapping_data <- bind_rows(recount3_mapping_information, local_data_additions_mapping_information)

# Create regression data

regression_data_final <- regression_data %>% 
  left_join(mapping_data %>% 
              select(recount_qc.aligned_reads..chrx, recount_qc.aligned_reads..chry, external_id), 
            by = c("run_accession" = "external_id"))

write.csv(regression_data_final, "data/inference_regression_data.csv", row.names = FALSE)
