#setwd('rmats_final')  
library(dplyr)
library(biomaRt)
library(grid)
library(ggplot2
        )
##summary stats
setwd("~/NIH/autoRNAseq/")
all_subtissues <- c("RPE_Stem.Cell.Line","RPE_Cell.Line","Retina_Adult.Tissue","RPE_Fetal.Tissue","ESC_Stem.Cell.Line","Cornea_Adult.Tissue","Cornea_Fetal.Tissue",
                    "Cornea_Cell.Line","Retina_Stem.Cell.Line","RPE_Adult.Tissue")
events <-c("SE.MATS.JC.txt", "RI.MATS.JC.txt",  "MXE.MATS.JC.txt", "A5SS.MATS.JC.txt", "A3SS.MATS.JC.txt") 
subtissue <- all_subtissues[1]

tot_events <- function(subtissue,events){
    all_data <- 
        return(total_events)
}

all_events_summary <- sapply(all_subtissues,function(y) lapply(events,function(x) paste('rmats_final',y,x, sep = '/')%>% 
                                                                   read.table(sep = '\t', header = T, stringsAsFactors = F))%>%sapply( nrow))%>%t()
event_names <- c('skipped_exon','retained_intron','mutually_exclusive_exon','alternative_5_prime_SS','alternative_3_prime_SS')
colnames(all_events_summary) <- event_names
t <- gridExtra::grid.table(all_events_summary)
ggsave('all_event_summary.png')

#####
setwd("~/NIH/autoRNAseq/rmats_final/Retina_Adult.Tissue")
skipped_exon <- read.table("SE.MATS.JC.txt",sep = '\t', header = T,stringsAsFactors = F)
retained_intron <- read.table("RI.MATS.JC.txt",sep = '\t', header = T,stringsAsFactors = F)
counts <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1")
skipped_exon[is.na(skipped_exon)] <- 1
skipped_exon[,counts] <- apply(skipped_exon[,counts],2,function(y) strsplit(y,',')%>%sapply(function(x) as.numeric(x)%>%sum))
total_event_counts <- rowSums(skipped_exon[,counts])
all_novel_SE <- read.table('all.SE.novelevents.txt',header = T,stringsAsFactors = F,sep = '\t')
novel_in_stissue <-semi_join(skipped_exon,all_novel_SE,by=c("chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")) 
# mean <- mean(total_event_counts)
# dev <- sd(total_event_counts)
# hist(total_event_counts)
# max(total_event_counts)
# keep events with an average of 50 counts per group, still a lot left might increase it 
keep <- which(total_event_counts>=100)#trimmed out about half
skipped_exon.filtered <- skipped_exon[keep,]
skipped_exon.filtered$score <- skipped_exon.filtered[,13:30]%>%rowSums()/2
retnet <- read.csv('../../ref/RetNet_genes.csv',stringsAsFactors = F,header = F)$V1

SE_retnet <- skipped_exon[skipped_exon$geneSymbol%in%retnet,]
SE_retnet.fil <-skipped_exon.filtered[skipped_exon.filtered$geneSymbol%in%retnet,]
unique(SE_retnet$geneSymbol)%>%length

mart <- useMart('ENSEMBL_MART_SNP','hsapiens_snp')
snps.chr15 <- getSNPlocs('15')
snps.chr15 <- getBM(attributes = c('refsnp_id','chrom_start','chrom_end'),filters = c('chr_name','start','end'),values =list(chr_name=15,start=65645611,end=65645703),mart = mart)


mart <- useMart('ENSEMBL_MART_ENSEMBL','hsapiens_gene_ensembl')


