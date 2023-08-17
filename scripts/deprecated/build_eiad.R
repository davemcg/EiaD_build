library(RSQLite)
library(tidyverse)
library(dbplyr)
library(readr)
library(rtracklayer)
library(dtplyr)
library(Biostrings)

gene_pool_2023 <- dbConnect(RSQLite::SQLite(), dbname = "eyeIntegration_2023_human.sqlite")

#Load necessary files
eyeIntegration23 <- data.table::fread("~/git/EiaD_build/data/eyeIntegration22_meta_2023_03_03.csv.gz") %>% 
  as_tibble() %>% 
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details %>% gsub('Adult Tissue ', '', .))) 
gene_counts <- vroom::vroom("counts/gene_tpm.csv.gz") %>% data.table::as.data.table()
tx_counts <- vroom::vroom("counts/tx_tpm.csv.gz") %>% data.table::as.data.table()
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


# make counts long 
TPM_gene <- gene_counts %>%
  filter(Gene %in% genes_above_zero)  %>%
  pivot_longer(-Gene, names_to = 'sample_accession') %>%
  collect() %>%
  rename(Gene = 'ID')
TPM_tx <- tx_counts %>%
  filter(Transcript %in% tx_above_zero)  %>%
  pivot_longer(-Transcript, names_to = 'sample_accession') %>%
  collect() %>% 
  rename(Transcript = 'ID')
#Calculating Mean Rank Decile
mean_rank_decile_gene <- TPM_gene %>% 
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



###########################################################################
# build sc database
#############################################################################
# load in scRNA pseudobulk counts data from plae/scEiaD
pb <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType_predict-Homo_sapiens.Staged.pseudoCounts.tsv.gz') 

######################################
# update gene names
gtf35 <- rtracklayer::readGFF("gencode.v35.annotation.gtf.gz")
gene_names_tib <- gtf35 %>% dplyr::select(gene_name, gene_id) %>% unique() %>% mutate(gene_id = gsub('\\.\\d+','', gene_id))
###########################################

# only keep pb genes that also exist gtf (because our PB matrix has mouse/macaque gene names)
pbH <- pb %>% filter(Gene %in% gene_names_tib$gene_id)
genes <- pbH$Gene
pbH <- pbH %>% dplyr::select(-Gene)

# scale norm (z scaling)
scalePB <- scale(as.matrix(pbH), center = FALSE)
# make gene id / ensg names
nnames <- genes %>% enframe() %>% left_join(gene_names_tib, by = c('value'='gene_id'))
if (nrow(nnames) != nrow(pbH)){
  print("OH NO")
  stop()
}

row.names(scalePB) <- paste0(nnames$gene_name, ' (', nnames$value, ')')
# make long
scalePB_out <-  scalePB %>% as_tibble(rownames = 'Gene')  %>% pivot_longer(-Gene) %>% separate(name, c("study_accession","CellType_predict","Stage"), "__")


# load in table info from scEiaD
# make with EiaD_build/scripts/pull_scEiaD.R
load("data/scEiaD_CT_table_2023_03_01.Rdata")
meta_filter <- fst::read_fst('~/data/scEiaD_2022_02//meta_filter.fst') 

# get gene names synced up
genes <- scEiaD_CT_table %>% ungroup() %>% dplyr::select(Gene) %>% arrange(Gene) %>% unique()
scalePB_out_filter <- scalePB_out %>% filter(Gene %in% genes$Gene, value > 0)


# get cell counts by the groupings
meta_counts <- meta_filter %>% 
  filter(organism == 'Homo sapiens', Source == 'Tissue') %>% 
  mutate(Organ = case_when(Organ == 'Eye' ~ 'Eye', TRUE ~ 'Body'))  %>% 
  group_by(study_accession, Organ,  CellType_predict, Stage) %>% 
  summarise(Count = n()) %>% 
  filter(Count > 50) 

# build some rough anatomical sections for the cell types for plot splitting
ct_site <- meta_filter %>% filter(Organ == 'Eye', !is.na(CellType)) %>% group_by(Tissue, CellType) %>% summarise(Count = n()) %>% filter(Count > 10) %>% group_by(CellType) %>% summarise(Tissue = paste0(Tissue, collapse = ', '))

ct_site <- ct_site %>% mutate(Site = case_when(CellType %in% c('B-Cell', 'Blood Vessel', 'Endothelial','Fibroblast','Macrophage', 'Mast', 'Melanocyte', 'Monocyte', 'Pericyte', 'Red Blood Cell', 'Schwann', 'Smooth Muscle', 'T/NK-Cell') ~ 'Eye',
                                               grepl("Ciliary", Tissue) ~ 'Front Eye',
                                               grepl("Retina|Choroid|RPE", Tissue) ~ 'Back Eye',
                                               TRUE ~ 'Front Eye')) 

#Creating scEiaD database
scEiaD_pool <- dbConnect(RSQLite::SQLite(), dbname = "~/git/eyeIntegration_app/inst/app/www/2022/scEiaD.sqlite")
#Writing the tables
dbWriteTable(scEiaD_pool, 'scEiaD_CT_table_info', scEiaD_CT_table, row.names = FALSE, overwrite = TRUE)
db_create_index(scEiaD_pool, 'scEiaD_CT_table_info', 'Gene')
dbWriteTable(scEiaD_pool, 'cell_types', scEiaD_CT_table %>% ungroup() %>% dplyr::select(CellType_predict) %>% unique() %>% arrange() %>% filter(!is.na(CellType_predict)), row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'gene_IDs', genes, row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'pseudoBulk',scalePB_out_filter  %>% dplyr::rename(zCount = value), row.names = FALSE, overwrite = TRUE)
db_create_index(scEiaD_pool, 'pseudoBulk', 'Gene')
dbWriteTable(scEiaD_pool, 'scEiaD_meta_counts', meta_counts, row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'ct_site', ct_site,  row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'Date_DB_Created', Sys.Date() %>% as.character() %>% enframe(name = NULL) %>% dplyr::select(DB_Created = value), row.names = FALSE, overwrite=TRUE)

# disconnect pools
dbDisconnect(scEiaD_pool)
dbDisconnect(gene_pool_2023)
