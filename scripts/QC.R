library(tximport)
library(tidyverse)
library(qsmooth)
library(Rtsne)
library(dbscan)
library(edgeR)
args=commandArgs(trailingOnly = T)
#args=c('sampleTable1218_tissues.tab','ref/gencodeAno_bsc.gtf','/data/swamyvs/autoRNAseq','transcript','testQC.tsv','ref/bad_mapping.txt')


sample_metadata = args[1]
gtf_file = args[2]
working_dir = args[3]
level = args[4] # transcript or gene level quantification
bad_map_file = args[5]
output_file = args[6]
qc_remove_output_file = args[7]
setwd(working_dir)

# Qsmooth > remove median counts > remove lowly expressed genes > tSNE > DBSCAN

# load gtf to construct gene <-> transcript_mapping
# load sample tissue classification for qsmooth

sample_design <- read_tsv(sample_metadata)
colnames(sample_design)[1] <- 'sample_accession'
gtf <- rtracklayer::readGFF(gtf_file) %>% dplyr::filter(type=='transcript')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id", "gene_type")]
#colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','Sub_Tissue','origin')
removal_log=data.frame()

files0 <-paste0('RE_quant_files/',sample_design$sample_accession, '/quant.sf')
samplenames <- strsplit(files0,'/' )%>% sapply(function(x) x[2])
sample_design <-  filter(sample_design, sample_accession%in%samplenames)

# load data at the transcript level or merge to the gene level
if (level == 'transcript') {
    txi.lsTPMs_tx <- tximport(files=files0,txOut = T, type = "salmon", countsFromAbundance = "lengthScaledTPM")
    txi.lsTPMs <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = "lengthScaledTPM")
    tpms_tx <- as.data.frame(txi.lsTPMs_tx$counts)
    colnames(tpms_tx) <- samplenames
} else {
    txi.lsTPMs <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = "lengthScaledTPM")
}
tpms <- as.data.frame(txi.lsTPMs$counts)
colnames(tpms)  <-  samplenames
## remove samples that failed to meet mapping %
bad_mapping_samples <- read.table(bad_map_file,stringsAsFactors = F) %>% filter(as.numeric(V7) < 35) %>% pull(V2)
if(length(bad_mapping_samples) > 0){
    bm = data.frame(sample=sample_design$sample_accession[sample_design$sample_accession%in%bad_mapping_samples],
                    reason='low salmon mapping rate')
    removal_log = rbind(removal_log, bm)
    }

sample_design <- sample_design[!sample_design$sample_accession%in%bad_mapping_samples,]
samplenames <- sample_design$sample_accession
tpms <- tpms[,samplenames]

tpms[is.na(tpms)] <- 0
keep_genes <- which(rowSums(tpms)>=ncol(tpms))# revove gene w less than an average count of .5
tpms <- tpms[keep_genes,]
sample_medians <- apply(tpms,2,function(x) median(x))
# more genes => lower medians, looks like 3 givees results most similar to david, see comparingQC.r
keep_median <- which(sample_medians>2)# criteria for removing samples with low counts
tpms <- tpms[,keep_median]
if( sum(!keep_median)>0){
    bm= data.frame(colnames(sample=tpms[!keep_median]),reason='low_median_counts')
    removal_log = rbind(removal_log,bm)

}


if (level == 'transcript'){
    #to keep the samples the same between the gene and tx level, filter samples only at the gene level
    #going from cut of 1 > .5 keeps about 40k tx's
    tpms_tx = tpms_tx[,samplenames]
    tpms <- tpms_tx
    tpms[is.na(tpms)] <- 0
    keep_genes <- which(rowSums(tpms)>=ncol(tpms)*1)# revove gene w less than an average count of 1
    tpms <- tpms[keep_genes,keep_median]

}
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
qs <- qsmooth(object = lsTPM_librarySize,groupFactor = as.factor(sample_design$Tissue))
lstpms_smoothed <- as.data.frame(qsmoothData(qs))

colnames(lstpms_smoothed) <- colnames(lsTPM_librarySize)
tpms_smoothed_filtered <- lstpms_smoothed

#cluster with tSNE, then run clustered data throught
# use 3000 most variable protein coding genes
if (level == 'transcript'){
	TPM <- lstpms_smoothed[anno %>% filter(gene_type == 'protein_coding', transcript_id %in% row.names(lstpms_smoothed)) %>% pull(transcript_id) %>% unique(), ]
} else {
	TPM <- lstpms_smoothed[anno %>% filter(gene_type == 'protein_coding', gene_name %in% row.names(lstpms_smoothed)) %>% pull(gene_name) %>% unique(), ]
}
variance <- apply(TPM, 1, var, na.rm=TRUE)
TPM$variance <- variance
high_var_TPM <- TPM %>% arrange(variance) %>% tail(3000) %>% select(-variance)
cor_outlier_identification <- function(df){
    cor_mat <- cor(df, method = 'spearman')# pairwise cor
    avg_cor <- rowMeans(cor_mat)#avg cor per sapmle
    grand_cor <- mean(avg_cor)# global avg cor
    dist <- avg_cor-grand_cor
    D <- dist/median(dist)
    D[is.nan(D)] <- Inf #possile for avg cor and grand cor to be  the same
    names(D) <- colnames(df)
    D
}

scores <- lapply(unique(sample_design$Sub_Tissue), function(x)filter(sample_design, Sub_Tissue == x) %>% 
                     pull(sample_accession) %>% high_var_TPM[,.] %>% cor_outlier_identification ) %>% do.call (c,.)

low_cor <- scores < -17.5

removed <-  names(scores)[low_cor]



if(length(removed) > 0){
    outliers = data.frame(sample = removed,reason = 'outlier')
    removal_log = rbind(removal_log,outliers)
}

trimmed_counts_smoothed <- tpms_smoothed_filtered  %>% rownames_to_column('ID') %>% select( -removed)
write_csv(trimmed_counts_smoothed, path = output_file)
write_tsv(removal_log, path = qc_remove_output_file)
