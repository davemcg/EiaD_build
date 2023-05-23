library(RSQLite)
library(tidyverse)
library(dbplyr)
library(readr)
library(rtracklayer)
library(dtplyr)
library(Biostrings)

gene_pool_2023 <- dbConnect(RSQLite::SQLite(), dbname = "eyeIntegration_2023_human_counts.sqlite")

#Load necessary files
eyeIntegration23 <- data.table::fread("~/git/EiaD_build/data/eyeIntegration22_meta_2023_03_03.csv.gz") %>% 
  as_tibble() %>% 
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details %>% gsub('Adult Tissue ', '', .))) 
gene_counts <- vroom::vroom("counts/gene_counts.csv.gz") %>% data.table::as.data.table()
tx_counts <- vroom::vroom("counts/tx_counts.csv.gz") %>% data.table::as.data.table()
#####################################################
mapping_data <- vroom::vroom("run_meta.csv.gz") %>% mutate(sample_accession = gsub('salmon_quant\\/|\\/aux_info\\/meta_info.json','',log)) %>% relocate(sample_accession) %>% select(-log)

full_annotations <- readDNAStringSet("default/gentrome.fa") %>% names() %>% enframe()
full_annotations <- full_annotations %>%
  mutate(gene_symbol = str_extract(value, 'gene_symbol:\\w+') %>% gsub(".*:","",.),
         transcript_id = str_extract(value, 'ENST\\d+\\.\\d+'),
         gene_id = str_extract(value, 'ENSG\\d+\\.\\d+'),
         gene_biotype = str_extract(value, 'gene_biotype:\\w+') %>% gsub(".*:","",.),
         transcript_biotype = str_extract(value, 'transcript_biotype:\\w+') %>% gsub(".*:","",.),
         description = str_extract(value, "description:.*")) %>%
  mutate(Gene = paste0(gene_symbol, ' (', gene_id, ')'),
         Transcript = paste0(gene_symbol, ' (', transcript_id, ')'))

gene_annotation <- full_annotations %>% select(Gene, gene_symbol, gene_id, gene_biotype, description) %>% unique()
gene_names <- gene_annotation %>% select(Gene) %>% unique()

tx_annotation <-  full_annotations %>% select(Transcript, gene_symbol, transcript_id, transcript_biotype, description) %>% unique()
tx_names <- tx_annotation %>% select(Transcript)


############ load count matrix to ID zero count genes#
count_matrix <- vroom::vroom("counts/gene_counts.csv.gz") %>% data.table::as.data.table()
genes_above_zero <- count_matrix$Gene[rowSums(count_matrix[,2:ncol(count_matrix)]) > 0]
tx_matrix <- vroom::vroom("counts/tx_counts.csv.gz") %>% data.table::as.data.table()
# transform tx ids into the gene (ens tx) format
tx_counts$Transcript <- tx_counts %>% select(Transcript) %>% left_join(tx_annotation, by = c("Transcript" = 'transcript_id')) %>% as_tibble() %>% pull(Transcript.y)
tx_matrix$Transcript <- tx_matrix %>% select(Transcript) %>% left_join(tx_annotation, by = c("Transcript" = 'transcript_id')) %>% as_tibble() %>% pull(Transcript.y)
tx_above_zero <- tx_matrix$Transcript[rowSums(tx_matrix[,2:ncol(tx_matrix)]) > 0]

# scale data
gene_counts_m <- gene_counts[,2:ncol(gene_counts)] %>% as.matrix()
row.names(gene_counts_m) <- gene_counts$Gene
tx_counts_m <- tx_counts[,2:ncol(tx_counts)] %>% as.matrix()
row.names(tx_counts_m) <- tx_counts$Transcript

s_gene_counts <- gene_counts_m %>% scale(., center = FALSE)
s_tx_counts <- tx_counts_m %>% scale(., center = FALSE)
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
  mutate(value = log1p(value)) %>% 
  left_join(., eyeIntegration23 %>% 
              dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession) %>%
              unique(), 
            by = "sample_accession") %>% 
  dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession, ID, value ) %>% 
  group_by(Tissue, Sub_Tissue, Source, Perturbation, Age, ID) %>% summarise(meanlsTPM = mean(value)) %>% 
  mutate(Rank = rank(-meanlsTPM, ties.method='first'), Decile = ntile(meanlsTPM, n = 10)) %>% 
  arrange(Tissue, Rank) %>% 
  collect() 

mean_rank_decile_tx <- TPM_tx %>%
  mutate(value = log1p(value)) %>% 
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
                             Age, sample_accession, Cohort,region, study_title, 
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

dbWriteTable(gene_pool_2023, 'Date_DB_Created', Sys.Date() %>% as.character() %>% enframe(name = NULL) %>% 
               select(DB_Created = value), row.names = FALSE, overwrite=TRUE)

