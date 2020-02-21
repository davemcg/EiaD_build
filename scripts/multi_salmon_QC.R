library(tximport)
library(tidyverse)
library(qsmooth)
library(Rtsne)
library(dbscan)
library(edgeR)
library(argparse)
args=commandArgs(trailingOnly = T)
#args=c('sampleTable1218_tissues.tab','ref/gencodeAno_bsc.gtf','/data/swamyvs/autoRNAseq','transcript','testQC.tsv','ref/bad_mapping.txt')
parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--sampleTable', action = 'store', dest = 'sample_metadata_file')
parser$add_argument('--countsFile', action = 'store', dest = 'counts_file')
parser$add_argument('--refGtf', action = 'store', dest = 'gtf_file')
parser$add_argument('--level', action = 'store', dest = 'level')
parser$add_argument('--mappingRates', action = 'store', dest = 'mapping_rate_file')
parser$add_argument('--smoothedTPMs', action = 'store', dest = 'output_file')
parser$add_argument('--fullCorTPMs',action = 'store', dest = 'fullCor_output_file')
list2env(parser$parse_args(), .GlobalEnv)
save.image('testing/msqc.args.Rdata')
#working_dir='/data/swamyvs/autoRNAseq'
# sample_metadata_file='sampleTable_2020_02_18_multi_salmon.tsv'
# dntx_gtf_file= '/data/swamyvs/ocular_transcriptomes_pipeline/data/gtfs/all_tissues.combined_NovelAno.gtf'
# quant_path='/data/swamyvs/autoRNAseq/results/salmon_quant/'
# track_file='/data/swamyvs/ocular_transcriptomes_pipeline/data/misc/TCONS2MSTRG.tsv'
# mapping_rate_file='/data/swamyvs/autoRNAseq/results/mapping_rates.txt'



setwd(working_dir)
process_lib_size_tabs <- function(file){
    df_messy <- read.delim(file, sep = ' ', header = F) 
    sample <- df_messy$V1 %>% as.character %>%  str_split('/') %>% 
        sapply( function (x){idx=(which(grepl('aux_info', x)) -1);return(x[idx])} )
    total_reads <-  df_messy$V3 %>% str_split(',') %>% sapply(function(x)x[1]) %>% as.numeric 
    percent_mapped <- df_messy$V5 %>% str_split(',') %>% sapply(function(x)x[1]) %>% as.numeric %>% {./100}
    df <- tibble(sample_accession=sample, total_reads, mapping_rate=percent_mapped)
    return(df)
}

########
# Qsmooth > remove median counts > remove lowly expressed genes > tSNE > DBSCAN

# load gtf to construct gene <-> transcript_mapping
# load sample tissue classification for qsmooth

gtf <- rtracklayer::readGFF(gtf_file) %>% dplyr::filter(type=='transcript') %>% dplyr::mutate(gene_type = 'protein_coding')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id", "gene_type")]
#colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','Sub_Tissue','origin')
sample_design <- read_tsv(sample_metadata_file) %>% rename(sample_accession=`#sample_accession`)
tpm_df <- readRDS(counts_file) 
id_name=ifelse(level == 'gene', 'gene_name', 'transcript_id')
tpm_df <- tpm_df %>% rename(ID = !!id_name)
tpms <- tpm_df %>% select(-ID ) %>% as.matrix
rownames(tpms) <- tpm_df$ID
mapping_rate <- process_lib_size_tabs(mapping_rate_file)
sample_design <- left_join(sample_design, mapping_rate, by = 'sample_accession') %>% 
    filter(sample_accession %in% colnames(tpms))
#normalize for library size
norm <- DGEList(tpms)
norm <- calcNormFactors(norm)
norm_counts <- norm$counts
#extract scaling factor for each sample and  multiply
correction <- norm$samples %>% data.frame() %>% .[['norm.factors']]
lsTPM_librarySize <- norm_counts %*% diag(correction)
colnames(lsTPM_librarySize) <- colnames(tpms)

#quantile normalize samples
sample_design <- sample_design %>% filter(sample_accession  %in% colnames(lsTPM_librarySize))
qs <- qsmooth(object = lsTPM_librarySize,group_factor = as.factor(sample_design$Tissue))
lstpms_smoothed <- as.data.frame(qsmoothData(qs))
colnames(lstpms_smoothed) <- colnames(lsTPM_librarySize)
tpms_smoothed_filtered <- lstpms_smoothed

# remove batch effect of studies
# and use model.matrix for ~sub_tissue design
# to correct TPM as per limma
# these corrected counts ARE NOT USED BY LIMMA for diff expression
# but are rather for: 
# 	tSNE
# 	PCA
# 	boxplots/heatmaps
#	cor_outlier_identification() below
sample_design$mapping_rate[is.na(sample_design$mapping_rate)] <- 0
design <- model.matrix(~ 0 + as.numeric(sample_design$mapping_rate) + as.factor(sample_design$Sub_Tissue))
tpms_smoothed_filtered_fullCor <- limma::removeBatchEffect(log2(tpms_smoothed_filtered+1), 
                                                           batch = as.factor(sample_design$study_accession),
                                                           design = design)
tpms_smoothed_filtered_fullCor <- (2 ** tpms_smoothed_filtered_fullCor) - 1
tpms_smoothed_filtered_fullCor[tpms_smoothed_filtered_fullCor < 0] <- 0
tpms_smoothed_filtered_fullCor <- tpms_smoothed_filtered_fullCor %>% data.frame()
colnames(tpms_smoothed_filtered_fullCor) <- colnames(tpms_smoothed_filtered)

#trimmed_counts_smoothed_pCor <- tpms_smoothed_filtered_pCor  %>% rownames_to_column('ID') %>% select(-removed)
# this for limma DE
write_csv(tpms_smoothed_filtered %>% rownames_to_column('ID') , path = output_file)
write_csv(tpms_smoothed_filtered_fullCor %>% rownames_to_column('ID') , path = fullCor_output_file)


