library(tximport)
library(tidyverse)
library(qsmooth)
library(Rtsne)
library(dbscan)
library(edgeR)
args=commandArgs(trailingOnly = T)
#args=c('sampleTable1218_tissues.tab','ref/gencodeAno_bsc.gtf','/data/swamyvs/autoRNAseq','transcript','testQC.tsv','ref/bad_mapping.txt')
setwd('/data/swamyvs/EiaD_build/')
save(args, file = 'testing/qc_Args.R')


sample_metadata = args[1]
gtf_file = args[2]
working_dir = args[3]
level = args[4] # transcript or gene level quantification
bad_map_file = args[5]
mapping_rate_file = args[6]
output_file = args[7]
qc_remove_output_file = args[8]
cor_scores = args[9]
fullCor_output_file = args[10]
setwd(working_dir)

# Qsmooth > remove median counts > remove lowly expressed genes > tSNE > DBSCAN

# load gtf to construct gene <-> transcript_mapping
# load sample tissue classification for qsmooth
sample_design <- read_tsv(sample_metadata) %>% select(-run_accession) %>% distinct
colnames(sample_design)[1] <- 'sample_accession'
gtf <- rtracklayer::readGFF(gtf_file) %>% dplyr::filter(type=='transcript') %>% dplyr::mutate(gene_type = 'protein_coding')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id", "gene_type")]
#colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','Sub_Tissue','origin')
removal_log=data.frame()


mapping_rate <- read_delim(mapping_rate_file, delim = ' ', col_names = F) %>%
  mutate(X2 = gsub('%', '', X2) %>% as.numeric()) %>%
  rename(sample_accession = X1, mapping_rate = X2) %>%
  mutate(mapping_rate = replace_na(mapping_rate, 0))

files0 <-paste0('RE_quant_files/',sample_design$sample_accession, '/quant.sf')
samplenames <- strsplit(files0,'/' )%>% sapply(function(x) x[2])
sample_design <-  filter(sample_design, sample_accession%in%samplenames) %>% left_join(mapping_rate)

# load data at the transcript level or merge to the gene level
if (level == 'transcript') {
    message('Running transcript-based filtering')
    removed_gene_samples = read_tsv(args[11])
    removed = removed_gene_samples$sample
    txi.tx.counts <- tximport(files=files0,txOut = T, type = "salmon", countsFromAbundance = 'no')
    counts_tx <- txi.tx.counts$counts 
    tpms_tx <- txi.tx.counts$abundance
    colnames(counts_tx) <-colnames(tpms_tx) <-   samplenames
    keep_samples = samplenames[!samplenames %in% removed]
    tpms_tx = tpms_tx[,keep_samples]
    keep_tx <- which(rowSums(tpms_tx)>=ncol(tpms_tx)*1)# revove gene w less than an average count of 1
    tpms_tx <- tpms_tx[keep_tx, ]
    counts_tx  = counts_tx[keep_tx, keep_samples]

    write_csv(counts_tx %>% as.data.frame %>% rownames_to_column('ID'), path = output_file)
    # this for sqlite file for boxplot / tSNE / PCA 
    write_csv(tpms_tx %>% as.data.frame %>% rownames_to_column('ID') , path = fullCor_output_file)
    write_tsv(removal_log, path = qc_remove_output_file)
    write_tsv(tibble(sample=''), path = cor_scores)

} else {
    message('Running Gene-based filtering')
    txi.gene.counts <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = 'no')
    counts <- txi.gene.counts$counts 
    tpms <- txi.gene.counts$abundance
    colnames(counts) <-colnames(tpms) <-   samplenames



    ## remove samples that failed to meet mapping %
    swaroop_samples <- filter(sample_design, study_accession  == 'SRP151763') %>% pull(sample_accession)
    swaroop_tpm <- tpms[,swaroop_samples]
    bad_mapping_samples <- mapping_rate %>% filter(mapping_rate< 60, !sample_accession%in% swaroop_samples) %>% pull(sample_accession)
    if(length(bad_mapping_samples) > 0){
        bm = data.frame(sample=sample_design$sample_accession[sample_design$sample_accession%in%bad_mapping_samples],
                        reason='low salmon mapping rate')
        removal_log = rbind(removal_log, bm)
        }


    sample_design <- sample_design %>% filter(!sample_accession %in% c(bad_mapping_samples, swaroop_samples))
    samplenames <- sample_design$sample_accession
    # add mapping rate

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

    # remove batch effect of studies
    # and use model.matrix for ~sub_tissue design
    # to correct TPM as per limma
    # these corrected counts ARE NOT USED BY LIMMA for diff expression
    # but are rather for: 
    # 	tSNE
    # 	PCA
    # 	boxplots/heatmaps
    #	cor_outlier_identification() below



    # for cor ID step below only keep protein coding genes (or transcripts)

    TPM <- tpms[anno %>% filter(gene_type == 'protein_coding', gene_name %in% row.names(tpms)) %>% pull(gene_name) %>% unique(), ]

    TPM <- TPM %>% data.frame()
    colnames(TPM) <- colnames(tpms)
    variance <- apply(TPM, 1, var, na.rm=TRUE)
    TPM$variance <- variance
    high_var_TPM <- TPM %>% arrange(variance) %>% tail(3000) %>% select(-variance) %>% as_tibble
    cor_outlier_identification <- function(df){
        if( length(df) == 0 |nrow(df) < 3){
            return('')
        }
        cor_mat <- cor(df, method = 'spearman')# pairwise cor
        avg_cor <- rowMeans(cor_mat)#avg cor per sapmle
        grand_cor <- mean(avg_cor)# global avg cor
        dist <- avg_cor-grand_cor
        D <- dist/median(dist)
        D[is.nan(D)] <- Inf #possile for avg cor and grand cor to be  the same
        names(D) <- colnames(df)
        D
    }
    sample_design = sample_design %>% filter(sample_accession %in% colnames(TPM))
    scores <- lapply(unique(sample_design$Sub_Tissue), function(x) {
        print(x)
        filter(sample_design, Sub_Tissue == x) %>% 
                        pull(sample_accession) %>% high_var_TPM[,.] %>% cor_outlier_identification }) %>% do.call (c,.)

    low_cor <- scores < -17.5

    removed <-  names(scores)[low_cor]

    if(length(removed) > 0){
        outliers = data.frame(sample = removed,reason = 'outlier')
        removal_log = rbind(removal_log,outliers)
    }
    save.image('testing/debug_qc.Rdata')
    tpms_complete <- cbind(tpms[, !colnames(tpms) %in% removed], swaroop_tpm[rownames(tpms),])
    counts_complete <- counts[rownames(tpms_complete), colnames(tpms_complete)]
    write_csv(counts_complete %>% as.data.frame %>% rownames_to_column('ID'), path = output_file)
    # this for sqlite file for boxplot / tSNE / PCA 
    write_csv(tpms_complete %>% as.data.frame %>% rownames_to_column('ID') , path = fullCor_output_file)
    write_tsv(removal_log, path = qc_remove_output_file)
    write_tsv(scores %>% as_tibble(rownames = 'sample'), path = cor_scores)
}

