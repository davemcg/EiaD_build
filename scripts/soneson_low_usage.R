#!/
library(tximport)
library(tidyverse)
# sonnesson 2016 inspired ID of low-usage tx
# http://biorxiv.org/content/early/2015/08/24/025387
# working dir, biowulf2
args=commandArgs(trailingOnly=T)
working_dir <- args[1]
gtf_file <- args[2]
#working_dir <- '~/NIH/autoRNAseq'
setwd(working_dir)
# pull in salmon files
files <- list.files(path='quant_files',recursive=TRUE,pattern='quant.sf', full.names = T)
# Gene TX to name conversion
gtf <- rtracklayer::readGFF(gtf_file) %>% filter(type=='transcript')

anno <- gtf %>% select(transcript_id, gene_name)
#save(anno,file = 'ref/txdb.Rdata')
# pull counts
txi <- tximport(files, type = "salmon", tx2gene = anno,  txOut = T)
# get counts
tx_c <- data.frame(txi$counts)

# remove samples with low median counts
sample_medians <- apply(tx_c,2,function(x) median(x))
# remove all with 0 median
samples_to_keep <- which(sample_medians >0)

tx_c <- tx_c[,samples_to_keep] %>%
    mutate(transcript_id =rownames(.)) %>%
    select(transcript_id, everything()) %>%
    left_join(anno, .) %>%
    arrange(gene_name)# 1 is corresponds to the trancript_id column
print(ncol(tx_c))
# get gene name added

#  sum counts by gene for all samples and calculate tx usage ratio
gene_sums <- summarizeToGene(txi,tx2gene = anno) %>%
    .[['counts']] %>%
    .[,samples_to_keep] %>%
    as.data.frame %>%
    mutate(gene_sums = rowSums(.), gene_name=rownames(.)) %>%
    select(gene_name, gene_sums) %>%
    arrange(gene_name)
gene_sums_tx <-left_join(gene_sums, tx_c)
all_ratios <- gene_sums_tx %>%
    select(-transcript_id, -gene_name, -gene_sums) %>%
    {. / gene_sums_tx$gene_sums}
all_ratios[is.nan(as.matrix(all_ratios))] <- 0
# find number of samples for each transcripts which are < 5% of the total
low_usage <- which(rowSums(all_ratios)<=.05)
print(paste(length(low_usage), 'Transcripts Removed'))
tx_to_remove <- gene_sums_tx[low_usage,'transcript_id']
write(tx_to_remove,'tx_for_removal.txt',sep='\n')
