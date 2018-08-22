library(dplyr)
setwd('/data/swamyvs/autoRNAseq')
sample_design <- read.table('sampleTable.txt',stringsAsFactors = F, header = F, sep = '\t')
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','sub-tissue','origin')
sample_design[!sample_design$tissue%in%c('Retina','RPE','Cornea','ESC'),]$`sub-tissue` <- 'body'
sample_design$path <- NA
for( i in 1:nrow(sample_design)){
  id <- sample_design[i,'sample_accession']
  paired <- sample_design[i,'paired']
  pathto <- paste0('tmp/',sample_design[i,'sub-tissue'],'/')
  if(paired=='y'){
    sample_design[i,'path'] <-paste0(paste0(pathto,id,'_1.fastq'),':', paste0(pathto,id,'_2.fastq'))
  }
  # need to figure out how to deal with mixed se/pe files
  # else{
  #   sample_design[i,'path'] <- paste0(pathto,id,'.fastq')
  # }
}
#ignore single ended files for now
sample_design <- sample_design[!is.na(sample_design$path),]
for(type in unique(sample_design$`sub-tissue`)){
  df <- sample_design[sample_design$`sub-tissue`==type,]
  writeLines(df$path,paste0('ref/',gsub(' ', '.' ,type),'.rmats.txt'),sep = ',')
  writeLines(paste0(df$path,'.gz'), paste0('ref/',type,'.tmp'),sep = ',')
}



