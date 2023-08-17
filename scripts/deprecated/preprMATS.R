library(dplyr)
setwd('/data/swamyvs/autoRNAseq')
args = commandArgs(trailingOnly=TRUE)
sample_design <- read.table(args[1],stringsAsFactors = F, header = F, sep = '\t')
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','sub-tissue','origin')
sample_design[!sample_design$tissue%in%c('Retina','RPE','Cornea','ESC'),]$`sub-tissue` <- 'body'
sample_design <-data.frame(sample_design, path=sapply(sample_design$sample_accession,function(x)paste0('STARbams_realigned/',x,'/Aligned.out.bam')),stringsAsFactors = F)

#ignore single ended files for now
sample_design <- sample_design[!is.na(sample_design$path),]
ttypes=unique(sample_design$sub.tissue)[-7]
for(type in ttypes){
  df <- filter(sample_design, sub.tissue==type ,paired=='y')
  writeLines(df$path,paste0('ref/',gsub(' ', '.' ,type),'_PE','.rmats.txt'),sep = ',')
  df <- filter(sample_design, sub.tissue==type ,paired=='n')
  writeLines(df$path,paste0('ref/',gsub(' ', '.' ,type),'_SE','.rmats.txt'),sep = ',')
}

body=read.table('synth_body.txt', stringsAsFactors=F,header=F)
for(i in 1:nrow(body)){
	paste0('STARbams_realigned/',body[1,i],'_PE' ,'/Aligned.out.bam')%>%writeLines('ref/body_PE.rmats.txt',sep=',')
	paste0('STARbams_realigned/',body[1,i],'_SE' ,'/Aligned.out.bam')%>%writeLines('ref/body_SE.rmats.txt',sep=',')
}



