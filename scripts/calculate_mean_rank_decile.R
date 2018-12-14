#calculate mean_rank_deciles for genes by sub_tissue
args = commandArgs(trailingOnly=TRUE)
working_dir <- args[1]
lsTPM_file <- args[2]
metadata_file <- args[3]
output_file <- args[4]

library(tidyverse)

setwd(working_dir)
lsTPM <- read_csv(lsTPM_file)

load(metadata_file)

core_tight$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight$sample_accession)
#core_tight$Sub_Tissue <- gsub('_',' - ',core_tight$Sub_Tissue)

mean_rank_decile <- lsTPM %>% 
  gather(sample_accession, value, -ID) %>% 
  left_join(.,core_tight) %>% 
  arrange(-value, Sub_Tissue) %>% 
  dplyr::select(Sub_Tissue, sample_accession, ID, value ) %>% 
  group_by(Sub_Tissue, ID) %>% summarise(meanlsTPM = mean(value)) %>% 
  mutate(Rank = rank(-meanlsTPM, ties.method='first'), Decile = ntile(meanlsTPM, n = 10)) %>% 
  arrange(Sub_Tissue, Rank)

write_tsv(mean_rank_decile, path = output_file)

