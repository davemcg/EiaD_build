library(RSQLite)
library(tidyverse)
library(pool)
library(dplyr)
library(readr)
library(rtracklayer)
library(dtplyr)


#Initialize gene_pool
gene_pool_2022 <- dbConnect(RSQLite::SQLite(), dbname = "~/git/eyeIntegration_app/inst/app/www/2022/eyeIntegration_2022_human.sqlite")

#Load necessary files
eyeIntegration22 <- data.table::fread("data/eyeIntegration22_meta_2023_03_03.csv.gz") %>% 
  as_tibble() %>% 
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details %>% gsub('Adult Tissue ', '', .))) 
gene_counts <- vroom::vroom("gene_counts/gene_TPM.csv.gz") %>% data.table::as.data.table()
############ load count matrix to ID zero count genes#
count_matrix <- vroom::vroom("gene_counts/gene_counts_matrix.csv.gz") %>% data.table::as.data.table()
genes_above_zero <- count_matrix$gene_id[rowSums(count_matrix[,2:ncol(count_matrix)]) > 0]
#####################################################
mapping_data <- bind_rows(vroom::vroom("mapping_data/recount3_mapping_information.csv.gz"),
                            vroom::vroom("mapping_data/local_mapping_information.csv.gz"),
                          vroom::vroom("mapping_data/gtex_mapping_information.csv.gz")) %>% 
  filter(sra.sample_acc.x %in% eyeIntegration22$sample_accession)
# http://duffel.rail.bio/recount3/human/new_annotations/gene_sums/human.gene_sums.G029.gtf.gz
gene_annotations <- 
  rtracklayer::import("data/human.gene_sums.G029.gtf.gz") %>% 
  as.data.frame()
# http://duffel.rail.bio/recount3/human/new_annotations/exon_sums/human.exon_sums.G029.gtf.gz
exon_annotations <- 
  rtracklayer::import("data/human.exon_sums.G029.gtf.gz") %>% 
  as.data.frame()

#Dropping duplicate colnames prior to left_join
exon_annotations <- exon_annotations %>% 
  select(-c(seqnames, strand, phase, gene_name, tag, start, source, level, end, type, havana_gene, width, score, gene_type)) %>% 
  unique()
gene_data <- gene_annotations %>% left_join(exon_annotations, by="gene_id")
gene_names <- gene_data %>% select(gene_id, gene_name) %>% unique()

#Left_joining gene_name to gene_counts
# using data.table for speed
gc_list <- list()
for (i in gene_counts$sample_accession %>% unique()){
  print(i)
  gc <- gene_counts[sample_accession == i]
  gc[gene_names, on = 'gene_id', gene_name := i.gene_name]

  gc_list[[i]] <- gc %>% collect()
  if (grepl('^\\d+', i)){
    # prepend X to start digit samples
    new_sample <- paste0('X',i)
    print(paste("New name: ", new_sample))
    gc_list[[i]]$sample_accession <- new_sample
  }
}


gene_counts_with_name <- bind_rows(gc_list)

#make gene ID df
gene_IDs <- gene_data %>% dplyr::select(gene_id, gene_name, gene_type) %>% 
  filter(gene_id %in% genes_above_zero) %>% 
  mutate(ID = paste0(gene_name, ' (', gene_id, ')')) %>% select(-gene_name, -gene_id) %>% 
  arrange(ID) %>% 
  unique()

# #make tx ID df
# gene_data <- gene_data %>% mutate(tx_name = paste(gene_data$transcript_name, paste("(", gene_data$transcript_id, ")", sep="")))
# tx_IDs <- gene_data %>% dplyr::select(tx_name, transcript_type) %>% arrange(tx_name) %>% unique()
# names(tx_IDs) <- c("ID", "transcript_type")

#Create lsTPM_gene to obtain mean_rank_deciles for genes by sub_tissue
#names(gene_counts_with_name) <- c("gene_id", "sample_accession", "value", "ID")
lsTPM_gene <- gene_counts_with_name %>% mutate(ID = paste0(gene_name, ' (', gene_id, ')')) %>% select(-gene_name, -gene_id) %>% 
  collect()

#Calculating Mean Rank Decile
mean_rank_decile <- lsTPM_gene %>% 
  left_join(., eyeIntegration22, by = "sample_accession") %>% 
  dplyr::select(Tissue, Sub_Tissue, Source, Perturbation, Age, sample_accession, ID, value ) %>% 
  group_by(Tissue, Sub_Tissue, Source, Perturbation, Age, ID) %>% summarise(meanlsTPM = mean(value)) %>% 
  mutate(Rank = rank(-meanlsTPM, ties.method='first'), Decile = ntile(meanlsTPM, n = 10)) %>% 
  arrange(Tissue, Rank) %>% 
  collect()

#Writing the tables
dbWriteTable(gene_pool_2022, 'mean_rank_decile_gene', mean_rank_decile, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2022,'mean_rank_decile_gene', 'ID')

dbWriteTable(gene_pool_2022, 'mapping_information', mapping_data, row.names = FALSE, overwrite = TRUE)

dbWriteTable(gene_pool_2022, 'metadata', eyeIntegration22, row.names = FALSE, overwrite = TRUE)

dbWriteTable(gene_pool_2022, 'lsTPM_gene', lsTPM_gene, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2022, 'lsTPM_gene', 'ID')

#dbWriteTable(gene_pool_2022, 'tx_IDs', tx_IDs, row.names = FALSE, overwrite = TRUE)
#db_create_index(gene_pool_2022, 'tx_IDs', 'ID')

dbWriteTable(gene_pool_2022, 'gene_IDs', gene_IDs, row.names = FALSE, overwrite = TRUE)
db_create_index(gene_pool_2022, 'gene_IDs', 'ID')

dbWriteTable(gene_pool_2022, 'Date_DB_Created', Sys.Date() %>% as.character() %>% enframe(name = NULL) %>% 
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
dbDisconnect(gene_pool_2022)
