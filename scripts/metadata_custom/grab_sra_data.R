library(reutils)
library(xml2)
library(tidyverse)

# 2022-10-20
# Importing Metadata
eyeIntegration22 <- read_csv("https://hpc.nih.gov/~parikhpp/EiaD/2022_metadata.csv", col_types = cols(...1 = col_skip()))


sra_meta_list <- list()

for (i in eyeIntegration22 %>% filter(study_accession != 'SRP012682') %>% pull(run_accession) %>% unique() ){
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
 
  sra_meta_list[[i]] <- out %>% data.frame()
  Sys.sleep(1) 
}


gtex_sra <- read_csv('~/Downloads/SraRunTable-3.txt')


gtex_srr <- eyeIntegration22 %>% 
  filter(study_accession == 'SRP012682') %>% 
  mutate(run_accession = gsub('\\.\\d$','',run_accession)) %>%  
  left_join(gtex_sra %>% 
              mutate(run_accession = paste0(gsub('-','.', biospecimen_repository_sample_id))) %>% 
              select(Run, run_accession) %>% 
              unique()) %>% 
  pull(Run)

gtex_sra_meta_list <- list()
for (i in gtex_srr %>% unique() ){
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
  
  gtex_sra_meta_list[[i]] <- out %>% data.frame()
  Sys.sleep(1) 
}


eiad_sra_meta <- bind_rows(sra_meta_list %>% bind_rows(.id = 'run_accession'),
                           gtex_sra_meta_list %>% bind_rows(.id = 'run_accession'))
                           
                           
                           
                           
eiad_sra_meta %>% write_tsv('data/2022-10-20_EiaD_SRA_meta.tsv.gz')
