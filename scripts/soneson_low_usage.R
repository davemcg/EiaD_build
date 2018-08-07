library(tximport)
library(data.table)
library(dplyr)
library(readr)
library(rtracklayer)

# sonnesson 2016 inspired ID of low-usage tx
# http://biorxiv.org/content/early/2015/08/24/025387

# working dir, biowulf2

working_dir <- '/data/swamyvs/autoRNAseq'

setwd(working_dir)

# pull in salmon files
files <- list.files(path='quant_files',recursive=TRUE,pattern='quant.sf')

# Gene TX to name conversion
gtf <- readGFF('ref/gencodeAno.gtf',tags = c('gene_id','gene_name','transcript_id'))[,c('gene_id','gene_name','transcript_id')]
anno <-  gtf[complete.cases(gtf),]

# pull counts
txi <- tximport(files, type = "salmon", tx2gene = anno[,2:3], reader = read_tsv, txOut = T)

# get counts
tx_c <- data.frame(txi$counts)

# remove samples with low median counts
sample_medians <- apply(tx_c[,2:ncol(tx_c)],2,function(x) median(x))
# remove all with 0 median
samples_to_keep <- names(sample_medians[sample_medians>0])
tx_temp <- tx_c %>% select_(.dots = samples_to_keep)
tx_c <- cbind(tx_c[,1:2],tx_temp)

# get gene name added
tx_c <- merge(data.frame(anno),tx_c,by.x='Transcript.ID',by.y='row.names')
# group by gene and sum counts by gene
gene_sums <- tx_c %>% group_by(Gene.Name) %>% summarise_each(funs(sum), 3:ncol(tx_c))
# add tx name back
gene_sum_tx <- tx_c %>% select(Transcript.ID, Gene.Name) %>% left_join(., gene_sums)
# matrix divide to get ratio
all_ratios <- tx_c[,3:ncol(tx_c)]/gene_sum_tx[,3:ncol(gene_sum_tx)]
# find number of samples for each transcripts which are < 10% of the total
low_usage <- apply(all_ratios,1, function(x) sum(x<0.05,na.rm=T))
# summary
summary(low_usage)
# num with all samples have < 10% tx usage
total_samples_left_minus_1 <- (ncol(tx_c)-3)
table(low_usage[low_usage>total_samples_left_minus_1])

tx_ids_to_keep <- gene_sum_tx[low_usage<total_samples_left_minus_1+1,] %>% .[['Transcript.ID']]
