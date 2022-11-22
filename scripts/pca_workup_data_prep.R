# PCA
library(tidyverse)

##################
# Load in metadata
emeta <- data.table::fread('data/eyeIntegration22_meta_2022_11_16.01.csv') %>% as_tibble()
##################

###################
# Load in counts
# samples are columns
# genes are rows
mat <- vroom::vroom("gene_counts/gene_counts_matrix.csv.gz") %>% data.frame()
###################


######################
# Light processing
# set rownames
row.names(mat) <- mat[,1]
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
#############################################

#########################
# calculate TPM
# pull counts out
raw_counts <- assay(dds, 'counts')
# get length info ('score' from recount3 gtf)
# http://duffel.rail.bio/recount3/human/new_annotations/gene_sums/human.gene_sums.G029.gtf.gz
gene_annotations <- 
  rtracklayer::import("data/human.gene_sums.G029.gtf.gz") %>% 
  as.data.frame()
len_info <- row.names(raw_counts) %>% enframe() %>% mutate(geneid = str_extract(value, 'ENSG\\d+') %>% gsub(')','',.)) %>% 
  left_join(gene_annotations %>% mutate(geneid = gsub('\\.\\d+','',gene_id))) 
# Divide each gene by transcript length
TPMlen <- apply( raw_counts, 2, function(x){ x / len_info$score } )
# Divide again the transcript length and apply 1e6 multiplier
TPM <- apply( TPMlen, 2, function(x) { x / sum(x) * 1E6} )
#########################################


# merge with TPMpb (pb is pseudoBulk)
## first run a diff script to get scEiaD/plae data
source('scripts/scEiaD_PBcounts_toTPM.R')
# remove .1, etc at the ends of ensgene
row.names(TPM) <- row.names(TPM) %>% gsub("\\.\\d+$","",.)
TPMall <- cbind(TPM[intersect(row.names(TPM) , row.names(TPMpb)),],
                TPMpb[intersect(row.names(TPM) , row.names(TPMpb)),])



#####################################
# add in human readable gene names
row.names(TPMall) <- row.names(TPMall) %>% 
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
meta_TPMall <- colnames(TPMall) %>% 
  enframe(value = 'sample_accession') %>% 
  left_join(gene_count_meta %>% dplyr::select(-run_accession, -name)%>% dplyr::distinct(), by = 'sample_accession') %>% 
  # now label the scRNA pb cell types
  mutate(CellType = gsub("^.*__","",sample_accession)) %>% 
  left_join(sc_tissue_ct %>% dplyr::rename(scTissue = Tissue), by = 'CellType') %>% 
  mutate(Tissue = case_when(is.na(Tissue) ~ scTissue,
                            TRUE ~ Tissue),
         Cohort = case_when(is.na(Cohort) ~ Organ,
                            TRUE ~ Cohort),
         Source = case_when(!is.na(scTissue) ~ 'scRNA',
                            TRUE ~ Source))
#########################################
save(TPM, TPMpb, TPMall, emeta, meta_TPMall, sc_tissue_ct, file = 'data/EiaD_pca_analysis.Rdata')
write_csv(TPMall %>% as_tibble(rownames = 'Gene'), file = 'data/EiaD_pca_analysis_TPMall.csv.gz')
write_csv(meta_TPMall, file = 'data/EiaD_pca_analysis_meta_TPMall.csv.gz')
#########################################