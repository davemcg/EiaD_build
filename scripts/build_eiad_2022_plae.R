library(RSQLite)
library(tidyverse)
library(dbplyr)
library(readr)
library(rtracklayer)
library(dtplyr)
library(Biostrings)


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