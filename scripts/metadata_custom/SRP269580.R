# SRP269580
library(reutils)
library(xml2)
library(tidyverse)

srr <- read_csv('data/SRP269580.sra.txt')


SRP269580_list <- list()
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
  
  SRP269580_list[[i]] <- out %>% data.frame()
  Sys.sleep(1) 
}

SRP269580_meta <- SRP269580_list %>% bind_rows() %>% 
  filter(!grepl("Chamber", EXPERIMENT.TITLE)) %>% 
  mutate(sample_accession = Pool.Member.IDENTIFIERS.PRIMARY_ID, 
         study_accession = STUDY.IDENTIFIERS.PRIMARY_ID,
         Cohort = 'Eye',
         Tissue = 'RPE', 
         Sub_Tissue = NA, 
         Source = 'iPSC', 
         Source_details = case_when(grepl("Transwell", EXPERIMENT.TITLE) ~ 'Transwell',
                                    grepl("PCL", EXPERIMENT.TITLE ) ~ 'PCL',
                                    grepl("Soy", EXPERIMENT.TITLE ) ~ 'Soy'),
         Age = NA,
         Age_Days = NA, 
         study_title = STUDY.DESCRIPTOR.STUDY_TITLE, 
         study_abstract = STUDY.DESCRIPTOR.STUDY_ABSTRACT,
         sample_attribute = glue::glue("{EXPERIMENT.TITLE}; {SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE}"),
         BioSample = Pool.Member.IDENTIFIERS.EXTERNAL_ID,
         #Library_Notes = EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_CONSTRUCTION_PROTOCOL
         region = NA) %>% 
  select(sample_accession:region) %>% 
  left_join(srr %>% select(run_accession = Run, BioSample))
