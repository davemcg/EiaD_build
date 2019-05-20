library(tidyverse)
library(RSQLite)
library(data.table)

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

sample_design <- read_tsv(sample_metadata)
colnames(sample_design)[1] <- 'sample_accession'
#colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','Tissue','Sub_Tissue','Origin')
sra_con <- dbConnect(RSQLite::SQLite(),sql_file)
sra_query <- "SELECT * FROM sra WHERE sample_accession='BLANK'"
eye_tissue <- c('Retina','RPE','Cornea','ESC', 'Lens')
core_tight <- lapply(sample_design$sample_accession,function(x) dbGetQuery(sra_con, gsub('BLANK',x,sra_query)))%>%do.call(rbind,.) %>%
  distinct%>%select(study_accession, study_title, study_abstract, sample_accession, run_accession, sample_attribute) %>%
  right_join(sample_design %>% select(sample_accession, Tissue, Sub_Tissue, Origin, Age_Days),by='sample_accession')
#core_tight$Origin <-  strsplit(core_tight[,'Sub_Tissue'],'_') %>% 
#  sapply(function(x) ifelse(x[[1]] %in% eye_tissue, x[[2]], 'Tissue'))
# add EMTAB
core_tight <- core_tight %>% mutate(Origin = case_when(grepl('MTAB', sample_accession) ~ 'Adult Tissue',
                                                       TRUE ~ Origin),
                                    study_title = case_when(grepl('MTAB', sample_accession) ~ 'RNAseq 50 Normal Human Retina',
                                                            TRUE ~ study_title),
                                    study_accession = case_when(grepl('MTAB', sample_accession) ~ 'E-MTAB-4377',
                                                                TRUE ~ study_accession),
                                    study_abstract = case_when(grepl('MTAB', sample_accession) ~ 'RNA-seq of post-mort retina donor without clinically relevant visual impairment. Ploy-A enriched. 75-nt paired-end. Short time lapse between tissue sampling and cDNA generation.',
                                                               TRUE ~ study_abstract))
# prettify formatting
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
# build swaroop AMD project metadata
amd <- read_tsv('~/git/EiaD_build/data/complete_swaroop_metadata_sraAdded.tsv') %>% data.table()
cols <- colnames(amd)[1:27]
amd[, sample_attribute := do.call(paste, 
                                  c(lapply(cols, function(x) paste(x, get(x), sep=": ")),              
                                    sep=" | "))]
amd <- amd %>% mutate(study_title = 'Integrated analysis of genetic variants regulating retinal transcriptome (GREx) identifies genes underlying age-related macular degeneration',
                      study_accession = 'SRP151763',
                      study_abstract = 'Age-related macular degeneration (AMD) is a complex multifactorial disease with at least 34 loci contributing to genetic susceptibility. To gain functional understanding of AMD genetics, we generated transcriptional profiles of retina from 453 individuals including both controls and cases at distinct stages of AMD. We integrated retinal transcriptomes, covering 13,662 protein-coding and 1,462 noncoding genes, with genotypes at over 9 million common single nucleotide polymorphisms (SNPs) for expression quantitative trait loci (eQTL) analysis of a tissue not included in Genotype-Tissue Expression (GTEx) and other large datasets. Cis-eQTL analysis revealed 10,474 genes under genetic regulation, including 4,541 eQTLs detected only in the retina. We then integrated the AMD-genome-wide association studies (GWAS) data with eQTLs and ascertained target genes at six loci. Furthermore, using transcriptome wide association analysis (TWAS), we identified 23 additional AMD-associated genes, including RLBP1, HIC1 and PARP12. Our studies expand the genetic landscape of AMD leading to direct targets for biological evaluation and establish the Genotype-Retina Expression (GREx) database as a resource for post-GWAS interpretation of retina-associated traits including glaucoma and diabetic retinopathy. Overall design: Retinal samples from 523 aged post-mortem human subjects from a spectrum of age-related macular degeneration (AMD) were RNA-seq profiled.') %>% rename(run_accession = run) %>% select(run_accession, study_title, sample_attribute, study_accession, study_abstract)
# build SRP105756
SRP105756 <- read_tsv('~/git/EiaD_build/data/SRP105756.txt') %>% data.table()
cols <- colnames(SRP105756)[9:12]
SRP105756[, sample_attribute := do.call(paste, 
                                  c(lapply(cols, function(x) paste(x, get(x), sep=": ")),              
                                    sep=" | "))]
SRP105756 <- SRP105756 %>% mutate(study_title = 'A gene expression study of the developing human eye',
                                  study_abstract = 'There is a vast range of tools with which to study the cells and tissues of the human body, however the scarcity of human embryonic material as a resource for direct study means that there is still limited data on the specific expression patterns which exist during human development. While embryonic studies are common place in other mammalian organisms and these have contributed greatly to our knowledge, characterising the events which occur during human development specifically is crucial in order to extrapolate the differences which exist between humans and other organisms and understand the human situation more completely. Establishing and maintaining collections of human developmental tissue for direct study requires significant time and resources, and is not a viable option for all. For this reason the creation of an atlas defining the key events which occur during human eye development would be of great value to the field. Through the Human Developmental Biology Resource we have collected human developmental eye tissue from post conception week 4 to 19 and performed gene expression studies using RNA sequencing. This study begins to reveal the spatiotemporal gene expression patterns which occur during normal human eye ontogenesis and represents a benchmark for comparison with development, disease and cellular differentiation studies. Overall design: 32 samples, which include 10 embryonic, 19 foetal and 3 adult samples',
                                  study_accession = 'SRP105756') %>% rename(run_accession = Run) %>% select(run_accession, sample_attribute, study_title, study_accession, study_abstract)
# hand annotate retina organoids SRP159246 and swaroop fetal retina SRP119766
sample_accession <- c('SRS3729741','SRS3729740','SRS3729739','SRS3729738','SRS3729737','SRS3729736','SRS3729735','SRS3729734','SRS3729733','SRS3729732','SRS3729731','SRS3729730','SRS3729729','SRS3729728','SRS3729727','SRS3729726','SRS3729725','SRS3729724','SRS3729723','SRS3729722','SRS3729721','SRS3729720','SRS3729719','SRS3729718','SRS3729717','SRS3729716','SRS3729715','SRS3729714','SRS3729713','SRS3729712','SRS3729711')
sample_attribute <- c('Day_250_3','Day_250_2','Day_250_1','Day_200_3','Day_200_2','Day_200_1','Day_181_3','Day_181_2','Day_181_1','Day_173_3','Day_173_2','Day_173_1','Day_158_2','Day_158_1','Day_128_3','Day_128_2','Day_128_1','Day_111_3','Day_111_2','Day_111_1','Day_69_3','Day_69_2','Day_69_1','Day_35_3','Day_35_2','Day_35_1','Day_20_2','Day_20_1','Day_10_3','Day_10_2','Day_10_1')
SRP159246 <- cbind(sample_accession, sample_attribute) %>% as_tibble() %>% mutate(study_title = 'Thyroid hormone signaling specifies cone subtypes in human retinal organoids')
sample_accession <- c('SRS2582170','SRS2582169','SRS2582168','SRS2582167','SRS2582166','SRS2582165','SRS2582164','SRS2582163','SRS2582162','SRS2582161','SRS2582160','SRS2582159','SRS2582158','SRS2582157','SRS2582156','SRS2582155','SRS2582154','SRS2582153','SRS2582152','SRS2582151','SRS2582171','SRS2582150','SRS2582149','SRS2582148','SRS2582147')
sample_attribute <- c('D132M','D96M','D73M','D132NC','D96NC','D59C-2','D59C','D132P','D96P','D73P','D59P-2','D59P','D136','D132','D125','D115','D107','D105','D94-2','D94','D80','D67','D57','D53','D52/54')
SRP119766 <- cbind(sample_accession, sample_attribute) %>% as_tibble() %>% mutate(study_title = 'Molecular anatomy of the developing human retina')

core_tight <- left_join(core_tight, bind_rows(SRP159246,SRP119766), by = "sample_accession") %>% 
  left_join(., bind_rows(amd, SRP105756), by = "run_accession") %>% 
  mutate(study_title = case_when(is.na(study_title.x) ~ study_title.y, 
                                 TRUE ~ study_title.x),
         sample_attribute = case_when(is.na(sample_attribute.x) ~ sample_attribute.y,
                                      TRUE ~ sample_attribute.x), 
         study_abstract = case_when(is.na(study_abstract.x) ~ study_abstract.y,
                                    TRUE ~ study_abstract.x)) %>% 
  select(sample_accession, study_accession, study_title, study_abstract, sample_attribute, 
         run_accession:Kept, -study_title.x, -study_title.y, -study_abstract.x,
         -study_abstract.y, -sample_attribute.x, -sample_attribute.y) 

# rewrite Origin from Sub_Tissue
#core_tight <- core_tight %>% mutate(Origin = case_when( grepl('Transformed', Sub_Tissue, ignore.case = T) ~ 'Cell Line',
#                                               grepl('Fetal', Sub_Tissue) ~ 'Fetal Tissue',
#                                               grepl('Stem', Sub_Tissue) ~ 'Stem Cell',
#                                               grepl('Organoid', Sub_Tissue) ~ 'Organoid',
#                                               grepl('Cell Line', Sub_Tissue) ~ 'Cell Line',
#                                               TRUE ~ 'Adult Tissue')) 

# update Kept field with specific reasons
# e.g. failed because of poor salmon alignment stats or tsne clustering (likely mislabled tissue/tube swap)
reasons <- read_tsv(samples_remove_gene)
core_tight <- left_join(core_tight, reasons, by = c('sample_accession' = 'sample')) %>%
  mutate(Kept = case_when(!is.na(reason) ~ reason,
                          TRUE ~ Kept)) %>%
  select(-reason)

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
