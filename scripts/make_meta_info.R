library(tidyverse)
library(RSQLite)

args=commandArgs(trailingOnly = T)
sample_metadata = args[1]
gtf_file = args[2]
sql_file=args[3]
gene_qc_file=args[4]
tx_qc_file=args[5]
samples_remove_gene=args[6]
samples_remove_tx=args[7]
setwd(args[8])
metadata_file <- args[9]
tx_file <- args[10]
gene_file <- args[11]
gene_tx_info <- args[12]

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
# hand annotate retina organoids SRP159246
sample_accession <- c('SRS3729741','SRS3729740','SRS3729739','SRS3729738','SRS3729737','SRS3729736','SRS3729735','SRS3729734','SRS3729733','SRS3729732','SRS3729731','SRS3729730','SRS3729729','SRS3729728','SRS3729727','SRS3729726','SRS3729725','SRS3729724','SRS3729723','SRS3729722','SRS3729721','SRS3729720','SRS3729719','SRS3729718','SRS3729717','SRS3729716','SRS3729715','SRS3729714','SRS3729713','SRS3729712','SRS3729711')
sample_attribute <- c('Day_250_3','Day_250_2','Day_250_1','Day_200_3','Day_200_2','Day_200_1','Day_181_3','Day_181_2','Day_181_1','Day_173_3','Day_173_2','Day_173_1','Day_158_2','Day_158_1','Day_128_3','Day_128_2','Day_128_1','Day_111_3','Day_111_2','Day_111_1','Day_69_3','Day_69_2','Day_69_1','Day_35_3','Day_35_2','Day_35_1','Day_20_2','Day_20_1','Day_10_3','Day_10_2','Day_10_1')
SRP159246 <- cbind(sample_accession, sample_attribute) %>% as_tibble() %>% mutate(study_title = 'Thyroid hormone signaling specifies cone subtypes in human retinal organoids')
core_tight <- left_join(core_tight, SRP159246, by = "sample_accession") %>% 
  mutate(study_title = case_when(is.na(study_title.x) ~ study_title.y, 
                                 TRUE ~ study_title.x),
         sample_attribute = case_when(is.na(sample_attribute.x) ~ sample_attribute.y,
                                      TRUE ~ sample_attribute.x)) %>% 
  select(sample_accession, study_accession, study_title, study_abstract, sample_attribute, 
         run_accession:Kept, -study_title.x, -study_title.y, 
         -sample_attribute.x, -sample_attribute.y) 
save(core_tight,file = metadata_file)

gtf <- rtracklayer::readGFF(gtf_file) %>% 
  dplyr::filter(type=='transcript')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id", "gene_type", "transcript_type")] %>% 
  filter(gene_name%in%gene_qc_tpms$ID)

gene_names = anno$gene_name %>% unique %>% sort()
save(gene_names, file= gene_file)
save(anno, file = gene_tx_info)
tx_qc_tpms=read_csv(tx_qc_file)
geneTX_names=gtf[,c("gene_name", "transcript_id")] %>%
  distinct%>%filter(transcript_id %in% tx_qc_tpms$ID)%>% split(.[,2]) %>%
  sapply(function(x) paste0(x[,1],' (',x[,2],')')) %>%
  sort()
save(geneTX_names, file=tx_file)
