# PCA
library(tidyverse)

##################
# Load in metadata
emeta <- data.table::fread('data/eyeIntegration22_meta_2023_03_03.csv.gz') %>% as_tibble()
##################

###################
# Load in counts
# samples are columns
# genes are rows
mat <- vroom::vroom("gene_counts/gene_counts_matrix.csv.gz") %>% data.frame()
###################


######################
# Light processing
# set rownames and remove .digit at the end
row.names(mat) <- mat[,1] %>% gsub('\\.\\d+','',.)
# remove gene column
mat <- mat[,2:ncol(mat)]
##########################


###########################
# Align metadata to counts matrix
gene_count_meta <- colnames(mat) %>% 
  enframe(value = 'run_accession') %>% 
  left_join(emeta %>% as_tibble() %>% 
              # these fields create duplicated entries for the sort-of-odd GTEX metadata
              dplyr::select(-sample_title, -gtex_sra_run_accession) %>% unique())
##############################

##############################
# collapse technical replicates to the sample level
library(DESeq2)
dds <- DESeqDataSetFromMatrix(mat, gene_count_meta, design = ~ Tissue)
dds <- collapseReplicates(dds, groupby = colData(dds)$sample_accession)

# now extract the sample accession level counts
sample_mat <- assay(dds, 'counts')
# merge with TPMpb (pb is pseudoBulk)

# load in scRNA pseudobulk counts data from plae/scEiaD
pb <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType_predict-Homo_sapiens.Staged.pseudoCounts.tsv.gz') 
pb_genes <- pb$Gene
pb_mat <- pb[,-1] %>% as.matrix()
row.names(pb_mat) <- pb_genes

mat_all <- cbind(sample_mat[intersect(row.names(sample_mat) , row.names(pb_mat)),],
                pb_mat[intersect(row.names(sample_mat) , row.names(pb_mat)),])



#####################################
# add in human readable gene names
gene_annotations <-
  rtracklayer::import("data/human.gene_sums.G029.gtf.gz") %>%
  as.data.frame()
row.names(mat_all) <- row.names(mat_all) %>% 
  enframe(value = 'gene_id') %>% 
  left_join(gene_annotations %>% mutate(gene_id = gsub("\\.\\d+$","",gene_id))) %>% 
  mutate(ID = paste0(gene_name, ' (', gene_id, ')')) %>% 
  pull(ID)

########################################################

######################################################
# recreate metadata
## cant use previous gene_count_meta as that was built around run_accession
## now we have sample level counts
## also sc data
meta_mat_all <- colnames(mat_all) %>% 
  enframe(value = 'sample_accession') %>% 
  left_join(gene_count_meta %>% dplyr::select(-run_accession, -name)%>% dplyr::distinct(), by = 'sample_accession') %>% 
  # now label the scRNA pb cell types
  mutate(CellType = case_when(grepl('__', sample_accession) ~ str_extract(sample_accession, '__.*__') %>% gsub('__','',.)),
         Age = case_when(grepl('Matur', sample_accession) ~ 'Adult',
                         grepl("Dev\\.$", sample_accession) ~ 'Fetal',
                         TRUE ~ Age),
         Source = case_when(grepl('__', sample_accession) ~ 'scRNA',
                            TRUE ~ Source))
#########################################
save(meta_mat_all, mat_all, gene_annotations, file = 'data/EiaD_pca_analysis.Rdata')
#write_csv(mat_all %>% as_tibble(rownames = 'Gene'), file = 'data/EiaD_pca_analysis_mat_all.csv.gz')
#write_csv(meta_mat_all, file = 'data/EiaD_pca_analysis_meta_mat_all.csv.gz')
#########################################