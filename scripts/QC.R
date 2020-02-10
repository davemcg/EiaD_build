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
mapping_rate_file = args[6]
counts_file = args[7]
output_file = args[8]
qc_remove_output_file = args[9]
cor_scores = args[10]
fullCor_output_file = args[11]
setwd(working_dir)

# Qsmooth > remove median counts > remove lowly expressed genes > tSNE > DBSCAN

# load gtf to construct gene <-> transcript_mapping
# load sample tissue classification for qsmooth
sample_design <- read_tsv(sample_metadata)
colnames(sample_design)[1] <- 'sample_accession'
gtf <- rtracklayer::readGFF(gtf_file) %>% dplyr::filter(type=='transcript') %>% dplyr::mutate(gene_type = 'protein_coding')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id", "gene_type")]
#colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','Sub_Tissue','origin')
removal_log=data.frame()


mapping_rate <- read_delim(mapping_rate_file, delim = ' ', col_names = F) %>%
  mutate(X2 = gsub('%', '', X2) %>% as.numeric()) %>%
  rename(sample_accession = X1, mapping_rate = X2) %>%
  mutate(mapping_rate = case_when(is.na(mapping_rate) ~ 0))

files0 <-paste0('RE_quant_files/',sample_design$sample_accession, '/quant.sf')
samplenames <- strsplit(files0,'/' )%>% sapply(function(x) x[2])
sample_design <-  filter(sample_design, sample_accession%in%samplenames)

# load data at the transcript level or merge to the gene level
if (level == 'transcript') {
    txi.lsTPMs <- tximport(files=files0,txOut = T, type = "salmon", countsFromAbundance = "lengthScaledTPM")
    #txi.lsTPMs <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = "lengthScaledTPM")
    txi.tx.counts <- tximport(files=files0,txOut = T, type = "salmon", countsFromAbundance = 'no')
	counts <- txi.tx.counts$counts 
	colnames(counts) <- samplenames
    counts <- counts %>% as_tibble(rownames = 'ID')
	write_csv(counts, path = counts_file) 
    rm(txi.tx.counts)
    rm(counts) 
	tpms_tx <- as.data.frame(txi.lsTPMs$counts)
    colnames(tpms_tx) <- samplenames
} else {
    txi.lsTPMs <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = "lengthScaledTPM")
    txi.gene.counts <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = 'no')
	counts <- txi.gene.counts$counts 
	colnames(counts) <- samplenames
    counts <- counts %>% as_tibble(rownames = 'ID')
	write_csv(counts, path = counts_file) 
    rm(counts)
    rm(txi.gene.counts) 
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
# add mapping rate
sample_design <- left_join(sample_design, mapping_rate, by = 'sample_accession')

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
#tpms_smoothed_filtered_fullCor <- tpms_smoothed_filtered

# remove study batch effect for counts going to limma
# can't be used a covariate as for some tissues, there is only one study, 
# so it is perfectly confounded and you can't build the model
#tpms_smoothed_filtered_pCor <- limma::removeBatchEffect(log2(tpms_smoothed_filtered+1),
#    batch = as.factor(sample_design$study_accession))
#tpms_smoothed_filtered_pCor <- (2 ** tpms_smoothed_filtered_pCor) - 1
#tpms_smoothed_filtered_pCor[tpms_smoothed_filtered_pCor < 0] <- 0
#tpms_smoothed_filtered_pCor <- tpms_smoothed_filtered_pCor %>% data.frame()
#colnames(tpms_smoothed_filtered_pCor) <- colnames(tpms_smoothed_filtered)

# for cor ID step below only keep protein coding genes (or transcripts)
if (level == 'transcript'){
	TPM <- tpms_smoothed_filtered_fullCor[anno %>% filter(gene_type == 'protein_coding', transcript_id %in% row.names(tpms_smoothed_filtered_fullCor)) %>% pull(transcript_id) %>% unique(), ]
} else {
	TPM <- tpms_smoothed_filtered_fullCor[anno %>% filter(gene_type == 'protein_coding', gene_name %in% row.names(tpms_smoothed_filtered_fullCor)) %>% pull(gene_name) %>% unique(), ]
}
TPM <- TPM %>% data.frame()
colnames(TPM) <- colnames(tpms_smoothed_filtered)
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

scores <- lapply(unique(sample_design$Sub_Tissue), function(x) filter(sample_design, Sub_Tissue == x) %>% 
                     pull(sample_accession) %>% high_var_TPM[,.] %>% cor_outlier_identification ) %>% do.call (c,.)

low_cor <- scores < -17.5

removed <-  names(scores)[low_cor]

if(length(removed) > 0){
    outliers = data.frame(sample = removed,reason = 'outlier')
    removal_log = rbind(removal_log,outliers)
}

#trimmed_counts_smoothed_pCor <- tpms_smoothed_filtered_pCor  %>% rownames_to_column('ID') %>% select(-removed)
# this for limma DE
write_csv(tpms_smoothed_filtered %>% rownames_to_column('ID') %>% select(-removed), path = output_file)
# this for sqlite file for boxplot / tSNE / PCA 
write_csv(tpms_smoothed_filtered_fullCor %>% rownames_to_column('ID') %>% select(-removed), path = fullCor_output_file)
write_tsv(removal_log, path = qc_remove_output_file)
write_tsv(scores %>% as_tibble(rownames = 'sample'), path = cor_scores)
