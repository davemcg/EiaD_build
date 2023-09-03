library(RSQLite)
library(tidyverse)
library(pool)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])
gene_tpm <- read_csv(args[2]) # smoothed_filtered_tpms_transcript.csv
tx_tpm <- read_csv(args[3]) # smoothed_filtered_tpms_gene.csv
load(args[4]) # tx_names.Rdata
load(args[5]) # limma_DE_listDF_gene.Rdata
limma_lists_gene <- limma_de_data
load(args[6]) # limma_DE_listDF_transcript.Rdata
limma_lists_tx <- limma_de_data 
load(args[7]) # de_comparison_name_list_gene.Rdata
load(args[8]) # all_vs_all_GO.Rdata
mrd_gene <- read_tsv(args[9]) # mean_rank_decile_gene.tsv
mrd_tx <- read_tsv(args[10]) # mean_rank_decile_transcript.tsv
load(args[11]) # core_tight.Rdata
load(args[12]) # gene_tx_gtf_info.Rdata
load(args[13]) # tSNE_coords.Rdata
sqlite_file <- args[14] # eyeIntegration_human_expression_2019_v100.sqlite

# create connection to sqlite file
expression_pool <- dbPool(SQLite(), dbname = sqlite_file)

# make data long, if not already
## gather gene and tx data from wide to long
long_gene <- gene_tpm %>% gather(key = 'sample_accession', value='value', -ID)
long_tx <- tx_tpm %>% gather(key = 'sample_accession', value='value', -ID)
## make gene ID df
gene_IDs <- gene_tpm %>% dplyr::select(ID) %>% arrange(ID) %>% 
              left_join(., anno %>% select(gene_name, gene_type), by = c('ID' = 'gene_name')) %>%
              unique()
## change tx ID from ENST to Gene (ENST)
tx_names <- geneTX_names %>% enframe() %>% select(ID = name, value)
tx_tpm <- left_join(tx_tpm, tx_names, by = c("ID")) %>% mutate(ENST = ID) %>% mutate(ID = value) %>% select(-value)
long_tx <- tx_tpm %>% select(-ENST) %>% gather(key = 'sample_accession', value='value', -ID)
## make tx ID df
tx_IDs <- tx_tpm %>% 
            left_join(., anno %>% select(transcript_id, transcript_type), by = c('ENST' = 'transcript_id')) %>%
            dplyr::select(ID, transcript_type) %>% arrange(ID) %>% unique()
## convert limma DE data from list to one DF
limma_DE_gene <- bind_rows(.id = 'Comparison', lapply(limma_lists_gene, rownames_to_column)) %>% 
					rename('ID' = rowname)
limma_DE_tx <- bind_rows(.id = 'Comparison', lapply(limma_lists_tx, rownames_to_column)) %>% 
					rename('ID' = rowname) %>% 
					left_join(., tx_names, by = c("ID")) %>% # add gene name to tx ID here, same as above
					mutate(ID = value) %>%
					select(-value)
## turn limma comparisons into a DF
de_tests <- de_comparison_contrast_names %>% enframe('Name') %>% dplyr::rename(Comparison = value)
# convert mean_rank_decile ID from ENST to ENGS (ENST)
tx_db <- geneTX_names %>% enframe() %>% select(ID = name, value)
mrd_tx <- left_join(mrd_tx, tx_db, by = 'ID') %>% mutate(ID = value) %>% dplyr::select(-value)

# write tables!
dbWriteTable(expression_pool, 'lsTPM_gene', long_gene, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool, 'lsTPM_gene', 'ID')
dbWriteTable(expression_pool, 'lsTPM_tx', long_tx, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool, 'lsTPM_tx', 'ID')
dbWriteTable(expression_pool, 'all_vs_all_go', all_vs_all_go, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool, 'all_vs_all_go', 'Set')
dbWriteTable(expression_pool, 'limma_DE_gene', limma_DE_gene, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool, 'limma_DE_gene', 'Comparison')
dbWriteTable(expression_pool, 'limma_DE_tx', limma_DE_tx, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool, 'limma_DE_tx', 'Comparison')
dbWriteTable(expression_pool, 'limma_DE_tests', de_tests, row.names = FALSE, overwrite = TRUE)
dbWriteTable(expression_pool, 'mean_rank_decile_gene', mrd_gene, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool,'mean_rank_decile_gene', 'ID')
dbWriteTable(expression_pool, 'mean_rank_decile_tx', mrd_tx, row.names = FALSE, overwrite = TRUE)
db_create_index(expression_pool,'mean_rank_decile_tx', 'ID')
dbWriteTable(expression_pool, 'metadata', core_tight, row.names = FALSE, overwrite = TRUE)
dbWriteTable(expression_pool, 'tx_IDs', tx_IDs, row.names = FALSE, overwrite = TRUE)
dbWriteTable(expression_pool, 'gene_IDs', gene_IDs, row.names = FALSE, overwrite = TRUE)
dbWriteTable(expression_pool, 'tSNE_bulk_RNA', all_tsne_plot_prepped, row.names = FALSE, overwrite = TRUE)
dbWriteTable(expression_pool, 'Date_DB_Created', Sys.Date() %>% as.character() %>% enframe(name = NULL) %>% select(DB_Created = value), 
    row.names = FALSE, overwrite=TRUE)
# close pool, disconnect
poolClose(expression_pool)
