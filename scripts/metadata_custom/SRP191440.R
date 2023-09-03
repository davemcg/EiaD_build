# SRP191440
library(xml2)
library(tidyverse)
full_xml <- read_xml('data/SRP191440.xml') %>% as_list()
exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
common_names <- lapply(exp_set_list, names) %>% purrr::reduce(intersect)

res <- lapply(exp_set_list , function(x) x[common_names]) %>%   do.call(rbind,.) %>% as.data.frame()


srr_srx <- read_csv('data/SRP191440.txt')

SRP191440_meta <- res %>% 
  filter(grepl('glucose', SAMPLE.TITLE)) %>% 
  mutate(sample_accession = Pool.Member.IDENTIFIERS.PRIMARY_ID, 
                      study_accession = STUDY.IDENTIFIERS.PRIMARY_ID,
                      Cohort = 'Eye',
               Tissue = 'RPE', 
               Source = 'Primary Culture', 
               Age_Days = NA, 
               study_title = STUDY.DESCRIPTOR.STUDY_TITLE, 
               study_abstract = STUDY.DESCRIPTOR.STUDY_ABSTRACT,
               sample_attribute = glue::glue("{SAMPLE.TITLE}; {SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE}"),
               BioSample = Pool.Member.IDENTIFIERS.EXTERNAL_ID,
               region = NA,
               Library_Notes = EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_CONSTRUCTION_PROTOCOL) %>% 
  select(sample_accession:Library_Notes) %>% 
  left_join(srr_srx %>% select(run_accession = Run, BioSample))
  

write_tsv(SRP191440_meta, 'data/SRP191440_eiad_meta.tsv.gz')

