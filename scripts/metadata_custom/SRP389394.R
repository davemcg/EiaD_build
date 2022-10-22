# SRP389394
library(reutils)
library(xml2)
library(tidyverse)

srr <- read_csv('data/SRP389394.sra.txt')
full_xml <- read_xml('data/SRP389394.xml') %>% as_list()
exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
common_names <- lapply(exp_set_list, names) %>% purrr::reduce(intersect)

SRP389394_xml <- lapply(exp_set_list , function(x) x[common_names]) %>%   do.call(rbind,.) %>% as.data.frame()


SRP389394_meta <- SRP389394_xml %>% mutate(sample_accession = Pool.Member.IDENTIFIERS.PRIMARY_ID, 
                                           study_accession = STUDY.IDENTIFIERS.PRIMARY_ID,
                                           Cohort = 'Eye',
                                           Tissue = 'RPE', 
                                           Sub_Tissue = NA, 
                                           Source = case_when(grepl('ahRPE', SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE) ~ 'Primary culture',
                                                              grepl('ARPE19', SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE) ~ 'Cell Line',
                                                              grepl('iPSC', SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE) ~ 'iPSC'),
                                           Source_details = case_when(grepl('ARPE19', SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE) ~ 'ARPE-19'),
                                           Age = NA,
                                           Age_Days = NA, 
                                           study_title = STUDY.DESCRIPTOR.STUDY_TITLE, 
                                           study_abstract = STUDY.DESCRIPTOR.STUDY_ABSTRACT,
                                           sample_attribute = glue::glue("{SAMPLE.TITLE}; {SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE}"),
                                           BioSample = Pool.Member.IDENTIFIERS.EXTERNAL_ID,
                                           #Library_Notes = EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_CONSTRUCTION_PROTOCOL
                                           region = NA) %>% 
  select(sample_accession:region) %>% 
  left_join(srr %>% select(run_accession = Run, BioSample))
