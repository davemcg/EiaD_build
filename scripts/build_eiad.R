library(RSQLite)
library(tidyverse)
library(pool)
library(dplyr)
library(readr)
library(rtracklayer)
library(dtplyr)


#Initialize gene_pool
gene_pool_2022 <- dbConnect(RSQLite::SQLite(), dbname = "eyeIntegration_2022_human.sqlite")

#Load necessary files
eyeIntegration22 <- data.table::fread("data/eyeIntegration22_meta_2022_11_06.01.csv") %>% as_tibble() %>% 
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details %>% gsub('Adult Tissue ', '', .))) 
gene_counts <- vroom::vroom("gene_counts/gene_TPM.csv.gz") %>% data.table::as.data.table()
############ load count matrix to ID zero count genes#
count_matrix <- vroom::vroom("gene_counts/gene_counts_matrix.csv.gz") %>% data.table::as.data.table()
genes_above_zero <- count_matrix$gene_id[rowSums(count_matrix[,2:ncol(count_matrix)]) > 0]
#####################################################
mapping_data <- vroom::vroom("mapping_data/mapping_data.csv.gz")
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

#load in scRNA data from plae/scEiaD
pb <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType-Homo_sapiens.pseudoCounts.tsv.gz') 

######################################
# TPM conversion needs the GTF to calculate the gene sizes
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
library(GenomicFeatures) 
# build DB
txdb <- makeTxDbFromGFF("gencode.v35.annotation.gtf.gz",format="gtf")
# aggregate by gene
exons_gene <- exonsBy(txdb,by="gene")
# calc gene sizes
gene_sizes <- as.data.frame(sum(width(reduce(exons_gene))))
# gene name (ENSG) has \.\d+ which are not used in the PB matrix
# so we erase them to facilitate the join later
row.names(gene_sizes) <- gsub('\\.\\d+','',row.names(gene_sizes))
colnames(gene_sizes)[1] <- 'length'
###########################################


# only keep pb genes that also exist gtf (because our PB matrix has mouse/macaque gene names)
pbH <- pb %>% filter(Gene %in% row.names(gene_sizes))
# get length info you calculated in the block above with makeTxDbFromGFF/etc/gene_sizes
## doing this wacky join to ensure the sizes are the in the same order as the matrix
len_info <- pbH$Gene %>% enframe() %>% left_join(gene_sizes %>% as_tibble(rownames = 'value')) 
# remove gene name col and turn into matrix
pbM <- pbH[,2:ncol(pbH)] %>% as.matrix()
# divide each gene by transcript length
TPMpblen <- apply( pbM, 2, function(x){ x / len_info$length } )
# divide again the transcript length and apply 1e6 multiplier
TPMpb <- apply( TPMpblen, 2, function(x) { x / sum(x) * 1E6} )
row.names(TPMpb) <- pbH$Gene

# load in table info
load("/Users/mcgaugheyd/data/EiaD/www/scEiaD_CT_table.Rdata")
names(scEiaD_CT_table) <- c("Gene", "CellType_predict", "organism", "study_accession", "Stage", "Cell # Detected", "Total Cells", "% of Cells Detected", "Meanlog2(Counts+1)")

#Creating scEiaD database
scEiaD_pool <- dbConnect(RSQLite::SQLite(), dbname = "scEiaD.sqlite")

#Writing the tables
dbWriteTable(scEiaD_pool, 'scEiaD_CT_table_info', scEiaD_CT_table, row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'cell_types', scEiaD_CT_table %>% ungroup() %>% select(CellType_predict) %>% unique() %>% arrange() %>% filter(!is.na(CellType_predict)), row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'gene_IDs', scEiaD_CT_table[,1], row.names = FALSE, overwrite = TRUE)
dbWriteTable(scEiaD_pool, 'Date_DB_Created', Sys.Date() %>% as.character() %>% enframe(name = NULL) %>% select(DB_Created = value), row.names = FALSE, overwrite=TRUE)
