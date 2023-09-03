library(xml2)
library(tidyverse)
full_xml <- read_xml('~/Downloads/SraExperimentPackage.xml') %>% as_list()
exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
common_names <- lapply(exp_set_list, names) %>% purrr::reduce(intersect)
res <-lapply(exp_set_list , function(x) x[common_names]) %>%   do.call(rbind,.) %>% as.data.frame()
write_tsv(res, 'new_sra_manual_Search.tsv')

