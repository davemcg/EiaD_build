# some code to yank in some of the recount3/snaptron files

# stuff needed (from unifier output)
# junctions.sqlite ## holds the junction locations (chr start stop) for each sample 
# junctions.bgz ## (bgzipped file with the locations....maybe not so needed?)
# junctions.bgz.tbi
# samples.tsv ## maps the rail_id to the run_accession

library(pool)
library(tidyverse)

# connect to the snaptron sqlite file
pool <- dbPool(RSQLite::SQLite(), dbname = '~/data/eiad_rse/junctions.sqlite')
## dangerously yank all info into memory
explosion <- pool %>% tbl('intron') %>% as_tibble() # only has one table 
# mapping info between rail_id (used in sqlite) and the run_accession
## also has junction_count and junction_coverage info
## needed for norm
snap_samp <- read_tsv('~/data/eiad_rse/samples.tsv')
emeta <- data.table::fread('data/eyeIntegration22_meta_2023_03_03.csv.gz') %>% as_tibble()

# CPM
## https://www.reneshbedre.com/blog/expression_units.html

# CPM = (counts * 1e6) / 
#   junction_coverage

big_explosion <- explosion %>% 
  # remove leading comma
  mutate(samples = gsub('^,','',samples)) %>% 
  # make data long
  separate_rows(samples, sep = ',') %>% 
  # split rail_id and counts into separate columns
  separate(samples, c('rail_id','counts'), ':') %>% 
  left_join(snap_samp %>% 
              mutate(rail_id = as.character(rail_id))) %>% 
  left_join(emeta, by = c("sample_id" = "run_accession")) %>% 
  mutate(counts = as.integer(counts),
         junction_coverage = as.integer(junction_coverage),
         CPM = (counts * 1e6) / 
           junction_coverage)

big_explosion_by_tissue <- big_explosion %>%
  group_by(chrom, start, end, 
           length, strand, annotated, donor,
           acceptor, left_annotated, right_annotated, Tissue) %>% 
  summarise(CPM = mean(CPM))

wide <- big_explosion_by_tissue %>% 
ungroup() %>% 
  pivot_wider(names_from = Tissue, values_from = CPM)
