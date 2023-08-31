library(RSQLite)
library(tidyverse)
library(dbplyr)
library(readr)
library(rtracklayer)
library(dtplyr)
library(Biostrings)
library(metamoRph)
gene_pool_2023 <- dbConnect(RSQLite::SQLite(), dbname = "~/git/eyeIntegration_app/inst/app/www/2023/eyeIntegration_2023_human_counts.sqlite")

#Load necessary files
eyeIntegration23 <- data.table::fread("~/git/EiaD_build/data/eyeIntegration23_meta_2023_08_28.csv.gz") %>% 
  as_tibble() %>%
  dplyr::rename(Age_Years = age_number) %>% 
  mutate(Sex = case_when(study_accession == 'SRP012682' & grepl("female", sample_title) ~ 'female',
                         study_accession == 'SRP012682' & grepl(" male", sample_title) ~ 'male',
                         TRUE ~ Sex)) %>% 
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details %>% gsub('Adult Tissue ', '', .))) 

gene_counts <- vroom::vroom("counts/gene_counts.csv.gz") %>% data.table::as.data.table()
tx_counts <- vroom::vroom("counts/tx_counts.csv.gz") %>% data.table::as.data.table()
####################
# apply sex ML via metamoRph
load('~/data/metamoRph_models/sex.Rdata')
gene_mat <- gene_counts[,2:ncol(gene_counts)] %>% as.matrix()
row.names(gene_mat) <- gene_counts$Gene %>% gsub('\\.\\d+','',.) #str_extract(gene_counts$Gene, 'ENSG\\d+')
mm_pca <- metamoRph::metamoRph(gene_mat, gtex_sex$PCA$rotation, gtex_sex$center_scale)
sex_labels <- metamoRph::model_apply(list_of_models = sex_model, 
                                     experiment_data = mm_pca)
eyeIntegration23 <- eyeIntegration23 %>% 
  left_join(sex_labels %>% dplyr::select(sample_accession = sample_id, Sex_ML = predict, Sex_Score = max_score))
##################### 
# output new meta
write_csv(eyeIntegration23, file= 'data/eyeIntegration23_meta_2023_08_28.built.csv.gz')

#####################################################
mapping_data <- vroom::vroom("run_meta.csv.gz") %>% mutate(sample_accession = gsub('salmon_quant\\/|\\/aux_info\\/meta_info.json','',log)) %>% relocate(sample_accession) %>% select(-log)

full_annotations <- data.table::fread('counts/gene_anno.csv.gz') %>% select(-Name, -length, -empty)
full_annotations <- full_annotations %>%
  mutate(Gene = paste0(gene_name, ' (', gene_id, ')'),
         Transcript = paste0(gene_name, ' (', transcript_id, ')'))

gene_annotation <- full_annotations %>% select(Gene, gene_name, gene_id, type) %>% unique()
gene_names <- gene_annotation %>% select(Gene) %>% unique()

tx_annotation <-  full_annotations %>% select(Transcript, gene_name, transcript_id, type) %>% unique()
tx_names <- tx_annotation %>% select(Transcript)

# pull full name and location
conv_table <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
                                    keys=gsub('\\.\\d+','',gene_annotation$gene_name),
                                    columns=c("SYMBOL", "MAP","GENENAME"), keytype="SYMBOL") %>% 
  group_by(SYMBOL) %>% 
  summarise(MAP = paste0(unique(MAP), collapse = ', '), 
            GENENAME = paste0(unique(GENENAME), collapse = ', ')) %>% 
  dplyr::rename(gene_name = SYMBOL, Location = MAP, Description = GENENAME)
gene_annotation <- gene_annotation %>% left_join(conv_table, by = 'gene_name')
tx_annotation <- tx_annotation %>% left_join(conv_table, by = 'gene_name')

############ load count matrix to ID zero count genes#
# transform gene ids into the gene (ens) format
gene_counts$Gene <- gene_counts %>% select(Gene) %>% 
  left_join(gene_annotation %>% select(nID = Gene, gene_id) %>% unique(), by = c('Gene' = 'gene_id')) %>% 
  pull(nID)
genes_above_zero <- gene_counts$Gene[rowSums(gene_counts[,2:ncol(gene_counts)]) > gene_counts %>% ncol - 1]
# transform tx ids into the gene (ens tx) format
tx_counts$Transcript <- tx_counts %>% select(Transcript) %>% 
  mutate(Transcript = str_extract(Transcript, 'ENST\\d+\\.\\d+')) %>% 
  left_join(tx_annotation, by = c("Transcript" = 'transcript_id')) %>% as_tibble() %>% pull(Transcript.y)
tx_above_zero <- tx_counts$Transcript[rowSums(tx_counts[,2:ncol(tx_counts)]) > gene_counts %>% ncol - 1]

# scale data
gene_counts_m <- gene_counts[,2:ncol(gene_counts)] %>% as.matrix()
row.names(gene_counts_m) <- gene_counts$Gene
tx_counts_m <- tx_counts[,2:ncol(tx_counts)] %>% as.matrix()
row.names(tx_counts_m) <- tx_counts$Transcript

s_gene_counts <- metamoRph::normalize_data(gene_counts_m, log1p = FALSE)
s_tx_counts <- metamoRph::normalize_data(tx_counts_m, log1p = FALSE)

# # make RSE
# coldata <- eyeIntegration23 %>%
#  dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession, study_accession, BioSample) %>%
#  unique()
# s_gene_counts_SE <- s_gene_counts[,coldata$sample_accession]
# rse <- SummarizedExperiment(assays = list(zcounts = s_gene_counts_SE), colData = coldata)

# make counts long 
TPM_gene <- s_gene_counts %>%
  as_tibble(rownames = 'ID') %>% 
  filter(ID %in% genes_above_zero)  %>%
  pivot_longer(-ID, names_to = 'sample_accession') %>%
  collect() 

TPM_tx <- s_tx_counts %>%
  as_tibble(rownames = 'ID') %>% 
  filter(ID %in% tx_above_zero)  %>%
  pivot_longer(-ID, names_to = 'sample_accession') %>%
  collect() 

#Calculating Mean Rank Decile
mean_rank_decile_gene <- TPM_gene %>% 
  #mutate(value = log2(value+1)) %>% 
  left_join(., eyeIntegration23 %>% 
              dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession) %>%
              unique(), 
            by = "sample_accession") %>% 
  dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession, ID, value ) %>% 
  group_by(Tissue, Sub_Tissue, Source, Perturbation, Age, ID) %>% summarise(meanlsTPM = mean(value, na.rm = TRUE)) %>% 
  mutate(Rank = rank(-meanlsTPM, ties.method='first'), Decile = ntile(meanlsTPM, n = 10)) %>% 
  arrange(Tissue, Rank) %>% 
  collect() 

mean_rank_decile_tx <- TPM_tx %>%
  #mutate(value = log2(value+1)) %>% 
  left_join(., eyeIntegration23 %>%
              dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession) %>%
              unique(),
            by = "sample_accession") %>%
  dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession, ID, value ) %>%
  group_by(Tissue, Sub_Tissue, Source, Perturbation, Age,ID) %>% summarise(meanlsTPM = mean(value)) %>%
  mutate(Rank = rank(-meanlsTPM, ties.method='first'), Decile = ntile(meanlsTPM, n = 10)) %>%
  arrange(Tissue, Rank) %>%
  collect() 

#Writing the tables
dbWriteTable(gene_pool_2023, 'mean_rank_decile_gene', mean_rank_decile_gene, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2023,'mean_rank_decile_gene', 'ID')

dbWriteTable(gene_pool_2023, 'mean_rank_decile_tx', mean_rank_decile_tx, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2023,'mean_rank_decile_tx', 'ID')

dbWriteTable(gene_pool_2023, 'mapping_information', mapping_data, row.names = FALSE, overwrite = TRUE)

dbWriteTable(gene_pool_2023, 'metadata', eyeIntegration23 %>%
               dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, 
                             Age, Age_Days, Age_Years, 
                             Sex, Sex_ML,
                             sample_accession, study_accession, Cohort,region, study_title, 
                             study_abstract, sample_attribute, BioSample, 
                             Source_details, Library_Notes, Comment, Origin, 
                             Sample_comment, geo) %>% unique(), 
             row.names = FALSE, overwrite = TRUE)

dbWriteTable(gene_pool_2023, 'lsTPM_gene', TPM_gene, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2023, 'lsTPM_gene', 'ID')

dbWriteTable(gene_pool_2023, 'lsTPM_tx', TPM_tx, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2023, 'lsTPM_tx', 'ID')

dbWriteTable(gene_pool_2023, 'gene_IDs', gene_annotation %>% filter(Gene%in% genes_above_zero) %>% rename(Gene = 'ID'), row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2023, 'gene_IDs', 'ID')

dbWriteTable(gene_pool_2023, 'tx_IDs', tx_annotation %>% filter(Transcript %in% tx_above_zero) %>% rename(Transcript = 'ID'), row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2023, 'tx_IDs', 'ID')

excluded_samples <- scan("data/excluded_samples.txt", what = "character")
dbWriteTable(gene_pool_2023, 'sample_outliers', excluded_samples %>% enframe() %>% select(outlier = value))

dbWriteTable(gene_pool_2023, 'Date_DB_Created', Sys.Date() %>% as.character() %>% enframe(name = NULL) %>% 
               select(DB_Created = value), row.names = FALSE, overwrite=TRUE)

system("mv eyeIntegration_2023_human_counts.sqlite ~/git/eyeIntegration_app/inst/app/www/2023/")
