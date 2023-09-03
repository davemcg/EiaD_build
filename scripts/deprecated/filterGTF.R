setwd('/data/swamyvs/autoRNAseq/')
library(rtracklayer)
library(dplyr)
gtf <- readGFF('ref/gencodeAno_bsc_tmp.gtf')%>%filter(gene_type=='protein_coding')
gtf[is.na(gtf)] <- '.'
cols_to_write <- colnames(gtf)[9:24]
gtf$to_write <- sapply(cols_to_write, function(x) paste0(x,' "', gtf[,x],'" '))%>%apply( 1, function(x) paste(x,collapse = ';'))
writeLines('##description: Gencode protein coding genes only
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf','ref/gencodeAno_bsc.gtf')
write.table(gtf[,c(1:8,25)],'ref/gencodeAno_bsc.gtf',sep = '\t', col.names = F,row.names = F, quote = F, append = T)
