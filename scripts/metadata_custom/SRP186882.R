# SRP186882
library(reutils)
library(xml2)
library(tidyverse)

srr <- read_csv('data/SRP186882.sra.txt')


SRP186882_list <- list()
for (i in srr %>% filter(Experiment %in% c('SRX9972989','SRX9972990','SRX9972991')) %>% pull(Run) ){
  print(i)
  sra_grab <- efetch(uid = i, db = 'sra', retmode = 'xml')
  full_xml <- read_xml(content(sra_grab, 'text')) %>% as_list()
  exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
  exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
  common_names <- lapply(exp_set_list, names) %>% purrr::reduce(intersect)
  
  res <- lapply(exp_set_list , function(x) x[common_names]) %>% 
    do.call(cbind,.) %>% 
    as_tibble(rownames = 'category') %>% 
    group_by(category) %>% 
    summarise(EXPERIMENT_PACKAGE = paste(EXPERIMENT_PACKAGE %>% unique(), collapse = ', ')) %>% 
    data.frame()
  out <- res[,2] %>% t()
  colnames(out) <- res[,1]
  
  SRP186882_list[[i]] <- out %>% data.frame()
  Sys.sleep(1) 
}

SRP186882_meta <- SRP186882_list %>% bind_rows() %>% 
  mutate(sample_accession = Pool.Member.IDENTIFIERS.PRIMARY_ID, 
         study_accession = STUDY.IDENTIFIERS.PRIMARY_ID,
         Cohort = 'Eye',
         Tissue = 'Cornea', 
         Sub_Tissue = 'Endothelium',
         Source = 'Native',
         Age = 'Adult',
         study_title = STUDY.DESCRIPTOR.STUDY_TITLE, 
         study_abstract = STUDY.DESCRIPTOR.STUDY_ABSTRACT,
         sample_title = glue::glue("{EXPERIMENT.TITLE}; {SAMPLE.TITLE}"),
         sample_attribute = glue::glue("{SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE}"),
         BioSample = Pool.Member.IDENTIFIERS.EXTERNAL_ID,
         #Library_Notes = EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_CONSTRUCTION_PROTOCOL,
         region = NA) %>% select(sample_accession:region) %>% 
  left_join(srr %>% mutate(run_accession = Run) %>% select(run_accession, BioSample) %>% as_tibble()) 

SRP186882_meta <- SRP186882_meta[,colnames(SRP186882_meta)[!grepl('\\.', colnames(SRP186882_meta))]]
