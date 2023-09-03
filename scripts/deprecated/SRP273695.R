# SRP273695
library(reutils)
library(xml2)
library(tidyverse)

srr <- read_csv('data/SRP273695.sra.txt')


SRP273695_list <- list()
for (i in srr %>% pull(Run) ){
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
  
  SRP273695_list[[i]] <- out %>% data.frame()
  Sys.sleep(1) 
}

# as we are updating existing entries!
samples <- c("SRS9761682", "SRS9761686", "SRS9761687", "SRS9761691", "SRS9761692", 
             "SRS9761694", "SRS9761695", "SRS9761697", "SRS9761699", "SRS9761701", 
             "SRS9761702", "SRS9761706", "SRS9761707", "SRS9761683", "SRS9761684", 
             "SRS9761685", "SRS9761688", "SRS9761689", "SRS9761690", "SRS9761693", 
             "SRS9761696", "SRS9761698", "SRS9761700", "SRS9761705", "SRS9761703", 
             "SRS9761704", "SRS9761708", "SRS9761709")

SRP273695_meta <- SRP273695_list %>% bind_rows() %>% 
  filter(EXPERIMENT.DESIGN.SAMPLE_DESCRIPTOR.IDENTIFIERS.PRIMARY_ID %in% samples) %>% 
  mutate(sample_accession = Pool.Member.IDENTIFIERS.PRIMARY_ID, 
         study_accession = STUDY.IDENTIFIERS.PRIMARY_ID,
         Cohort = 'Eye',
         Tissue = case_when(grepl("Retina", SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE) ~ 'Retina',
                            TRUE ~ 'RPE'),
         Source = 'Tissue',
         Age = 'Adult',
         region = case_when(grepl("Periphral", SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE) ~ 'Peripheral',
                            TRUE ~ 'Macula'),
         Sub_Tissue = region,
         study_title = STUDY.DESCRIPTOR.STUDY_TITLE, 
         study_abstract = STUDY.DESCRIPTOR.STUDY_ABSTRACT,
         sample_title = glue::glue("{EXPERIMENT.TITLE}; {SAMPLE.TITLE}"),
         sample_attribute = glue::glue("{SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE}"),
         BioSample = Pool.Member.IDENTIFIERS.EXTERNAL_ID,
         data_location = 'local') %>% 
  select(sample_accession:BioSample) %>% 
  left_join(srr %>% select(run_accession = Run, BioSample)) 


SRP273695_meta <- SRP273695_meta[,colnames(SRP273695_meta)[!grepl('\\.', colnames(SRP273695_meta))]]

