library(tidyverse)
library(RSQLite)

args=commandArgs(trailingOnly = T)
sample_metadata = args[1]
gtf_file = args[2]
sql_file=args[3]
gene_qc_file=args[4]
tx_qc_file=args[5]
setwd(args[6])
metadata_file <- args[7]
tx_file <- args[8]
gene_file <- args[9]

sample_design <- read.table(sample_metadata, stringsAsFactors = F, header=F, sep = '\t')
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','Tissue','Sub_Tissue','Origin')
sra_con <- dbConnect(RSQLite::SQLite(),sql_file)
sra_query <- "SELECT * FROM sra WHERE sample_accession='BLANK'"
eye_tissue <- c('Retina','RPE','Cornea','ESC')
core_tight <- lapply(sample_design$sample_accession,function(x) dbGetQuery(sra_con, gsub('BLANK',x,sra_query)))%>%do.call(rbind,.) %>%
  distinct%>%select(study_accession, study_title, study_abstract, sample_accession, run_accession, sample_attribute) %>%
  left_join(sample_design[,c(1,4,5)],by='sample_accession')
core_tight$Origin <-  strsplit(core_tight[,'Sub_Tissue'],'_') %>% 
  sapply(function(x) ifelse(x[[1]] %in% eye_tissue, x[[2]], 'Tissue'))
# add EMTAB
core_tight <- bind_rows(core_tight, 
                        sample_design %>% 
                          filter(grepl('MTAB', sample_accession)) %>% 
                          select(sample_accession, run_accession, Tissue, Sub_Tissue) %>% 
                          mutate(Origin = 'Adult Tissue',
                                 study_title = 'RNAseq 50 Normal Human Retina',
                                 study_accession = 'E-MTAB-4377',
                                 study_abstract = 'RNA-seq of post-mort retina donor without clinically relevant visual impairment. Ploy-A enriched. 75-nt paired-end. Short time lapse between tissue sampling and cDNA generation.'))
core_tight <- core_tight %>% 
  as.tibble() %>% 
  rowwise() %>% 
  mutate(Tissue = gsub('\\.', ' ', Tissue),
         Sub_Tissue = gsub('\\.', ' ', Sub_Tissue),
         Origin = gsub('\\.', ' ', Origin)) %>% 
  mutate(Tissue = gsub('_', ' - ', Tissue),
         Sub_Tissue = gsub('_', ' - ', Sub_Tissue),
         Origin = gsub('_', ' - ', Origin))


gene_qc_tpms=read_csv(gene_qc_file)
samples_kept=colnames(gene_qc_tpms) %>% 
  data.frame(sample_accession=.,Kept='Kept',stringsAsFactors = F)
core_tight=left_join(core_tight,samples_kept, by='sample_accession')
core_tight$Kept[is.na(core_tight$Kept)]='removed'
save(core_tight,file = metadata_file)

gtf <- rtracklayer::readGFF(gtf_file) %>% 
  dplyr::filter(type=='transcript')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id")] %>% 
  filter(gene_name%in%gene_qc_tpms$ID)

genes= anno$gene_name %>% unique
save(genes, file= gene_file)
tx_qc_tpms=read_csv(tx_qc_file)
tx=gtf[,c("gene_name", "transcript_id")] %>%
  distinct%>%filter(transcript_id %in% tx_qc_tpms$ID)%>% split(.[,2]) %>%
  sapply(function(x) paste0(x[,1],' (',x[,2],')'))
save(tx, file=tx_file)
