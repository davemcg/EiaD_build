setwd('/data/swamyvs/autoRNAseq')
library(rtracklayer)
library(dplyr)
gtf <- readGFF('ref/comb.gtf')
ref <- readGFF('ref/gencodeAno_bsc.gtf')
gtf <- gtf[!is.na(gtf$ref_gene_id),]
gtf[is.na(gtf)] <- '.'
ref_genes <- filter(ref,type=='gene')

ref_genes <- ref_genes[ref_genes$gene_id%in%gtf$ref_gene_id,  ]


gtf <- gtf[,c(1:8,13,10,12,11)]
ref_genes <- ref_genes[,c(1:9,14,11,20)]
cn <- c("seqid", "source", "type", "start", "end" ,"score", "strand", "phase" ,"gene_id", "transcript_id" ,"gene_name", "exon_number"    )
colnames(gtf) <- colnames(ref_genes) <- cn

gtf_final <- rbind(gtf,ref_genes)
gtf_final[is.na(gtf_final)] <- '.'
to_write <- paste( paste0('gene_id ',gtf_final$gene_id),
                       paste0(' transcript_id ','"',gtf_final$transcript_id,'"'),
                       paste0(' exon_numer ','"',gtf_final$exon_number,'"'),
                       paste0(' gene_name ','"',gtf_final$gene_name,'"'),sep = ';')
gtf_final$to_write <- to_write  
gtf_final <- arrange(gtf_final, seqid, start, desc(type))
gtf_final <- arrange(gtf_final, seqid,start, gene_id, transcript_id,type)
write.table(gtf_final[,c(1:8,13)],'ref/combined_final.gtf',sep = '\t',quote = F, col.names = F,row.names = F)


