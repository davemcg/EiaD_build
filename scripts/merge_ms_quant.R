library(tximport)
library(tidyverse)
library(argparse)
parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--quantPath', action = 'store', dest = 'quant_path')
parser$add_argument('--trackFile', action = 'store', dest = 'track_file')
parser$add_argument('--sampleMetadata', action = 'store', dest = 'sample_metadata_file')
parser$add_argument('--txTPMfile', action = 'store', dest = 'tx_tpm_file')
parser$add_argument('--geneTPMfile', action = 'store', dest = 'gene_tpm_file')
list2env(parser$parse_args(), .GlobalEnv)

read_salmon <- function(path){
    print(path)
    qfiles <- list.files(path,pattern='quant.sf', recursive=T, full.names = T)
    name_idx <- str_split(qfiles[1], '/')[[1]] %>% grep('quant.sf', .) %>% {. -1}
    names <- str_split(qfiles,'/') %>% sapply(function(x) x[name_idx])
    txi <- tximport::tximport(files = qfiles, type = 'salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
    colnames(txi$counts) <- names
    colnames(txi$abundance) <- names
    return(txi)
}



sample_design <- read_tsv(sample_metadata_file) %>% rename(sample_accession=`#sample_accession`)
DNTX_gtf <- rtracklayer::readGFF(dntx_gtf_file)
DNTX_t2g <- DNTX_gtf %>% as_tibble %>% select(transcript_id, gene_name) %>% distinct

ms_subtissues <- unique(sample_design$ms_subtissue)
all_quant_paths <- paste0(quant_path, ms_subtissues)
converter_tab <- read_tsv(track_file)  %>% 
    inner_join(DNTX_t2g)

all_quant <- lapply(all_quant_paths, function(x) read_salmon(x) )

proc_tx_quant <- function(tissue, txi){
    df <- txi$counts 
    ctab <- {converter_tab[,c('transcript_id', tissue)]} %>% filter(!is.na(.[,2]))
    df %>% as.data.frame %>% mutate(!!tissue := rownames(.)) %>% 
        inner_join(ctab) %>% select(-(!!tissue))
}

proc_gene_quant <- function(tissue, txi){
    t2g <-  {converter_tab[,c('transcript_id', tissue, 'gene_name')]} %>% 
        filter(!is.na(.[,2])) %>% 
        {.[,2:3]}
    gene_quant <- summarizeToGene(object = txi, tx2gene= t2g, countsFromAbundance = 'lengthScaledTPM') %>% {.[['counts']]}
    gene_quant %>% as.data.frame %>% mutate(gene_name = rownames(.)) %>% select(gene_name, everything())
    
}

tpms_tx <- lapply(seq_along(ms_subtissues), function(i) proc_tx_quant(ms_subtissues[i], all_quant[[i]])) %>% 
    reduce(full_join)
tpms_gene <- lapply(seq_along(ms_subtissues), function(i) proc_gene_quant(ms_subtissues[i], all_quant[[i]])) %>% 
    reduce(full_join)
tpms_tx[is.na(tpms_tx)]<-0
tpms_gene[is.na(tpms_gene)] <- 0
saveRDS(tpms_tx, tx_tpm_file)
saveRDS(tpms_gene, gene_tpm_file)

