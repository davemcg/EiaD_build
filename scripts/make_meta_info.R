library(tidyverse)
library(RSQLite)
library(SRAdb)
args=commandArgs(trailingOnly = T)
sample_metadata = args[1]
gtf_file = args[2]
sqlfile=args[3]
gene_qc_file=args[4]
tx_qc_file=args[5]

sample_design <- read.table(sample_metadata, stringsAsFactors = F, header=F, sep = '\t')
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','sub-tissue','origin')
sra_con <- dbConnect(RSQLite::SQLite(),sqlfile)
que <- "SELECT * FROM sra WHERE sample_accession='BLANK'"
eye_tissue <- c('Retina','RPE','Cornea','ESC')
core_tight <- lapply(sample_design$sample_accession,function(x) dbGetQuery(sra_con, gsub('BLANK',x,que)))%>%do.call(rbind,.)%>%
  distinct%>%select(study_accession ,study_title,study_abstract,sample_accession,sample_attribute)%>%
  left_join(sample_design[,c(1,4,5)],by='sample_accession')
core_tight$origin <-  strsplit(core_tight[,'sub-tissue'],'_')%>%sapply(function(x)ifelse(x[[1]]%in%eye_tissue,x[[2]],'Tissue'))
gene_qc_tpms=read_csv(gene_qc_file)
samples_kept=colnames(gene_qc_tpms)%>%data.frame(sample_accession=.,kept='kept',stringsAsFactors = F)
core_tight=left_join(core_tight,samples_kept, by='sample_accession')
core_tight$kept[is.na(core_tight$kept)]='removed'
save(core_tight,file = 'results/core_tight.Rdata')

gtf <- rtracklayer::readGFF(gtf_file) %>% dplyr::filter(type=='transcript')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id")]%>%filter(gene_name%in%gene_qc_tpms$ID)

genes= anno$gene_name %>% unique
save(genes, file='results/gene_names.Rdata')
tx_qc_tpms=read_csv(tx_qc_file)
tx=gtf[,c("gene_name", "transcript_id")]%>%distinct%>%filter(transcript_id%in%tx_qc_tpms$ID)%>% split(.[,2])%>%
  sapply(function(x) paste0(x[,1],' (',x[,2],')'))
save(tx, file='results/tx_gene_names.Rdata')
