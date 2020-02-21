library(xml2)
library(tidyverse)
setwd('/Volumes/data/autoRNAseq/')
setwd('/data/swamyvs/autoRNAseq/')
xml <- read_xml('ref/swaroop_sra_metadata.xml') %>% as_list
samples <- xml$EXPERIMENT_PACKAGE_SET

res <- lapply(samples, function(x) unlist(x, recursive = T)) %>% do.call(rbind, .) %>% as_tibble()
rownames(res) <- 1:nrow(res)
missing_age_sample <- res[202,]# had to manually look at this, everything else looks is fine
res <- res[-202,]

res <- res %>% select(id=EXPERIMENT.TITLE, sample=EXPERIMENT.DESIGN.SAMPLE_DESCRIPTOR.IDENTIFIERS.PRIMARY_ID, run=RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID)
metadata <- read_csv('ref/swaroop_metadata.csv') %>% 
    mutate(patient_number=str_split(r_id, '_') %>% sapply(function(x) x[[1]]) %>% as.numeric)
res <- res %>% mutate(patient_number= id %>% as.character %>% str_split('_|-') %>% sapply(function(x) x[[3]]) %>% as.numeric)

complete_metadata <- inner_join(res, metadata) %>% mutate(tissue='Retina', 
                                                          subtissue= case_when(mgs_level==1 ~ 'Retina_Adult.Tissue',
                                                                               mgs_level==2 ~ 'Retina_Adult.Tissue_AMD.MGS-2',
                                                                               mgs_level==3 ~ 'Retina_Adult.Tissue_AMD.MGS-3',
                                                                               mgs_level==4 ~ 'Retina_Adult.Tissue_AMD.MGS-4'),
                                                          paired='y',
                                                          origin='Adult.Tissue',
                                                          Age_Days='.',
                                                          Source='swaroop_AMD') 
write_tsv(complete_metadata, 'ref/complete_swaroop_metadata_sraAdded.tsv')
metadata_trimmed <- complete_metadata %>% select(sample, run, paired, tissue, subtissue, origin, Age_Days, Source)
write_tsv(metadata_trimmed, 'ref/sampleTable_swaroopAMD_formatted.tsv')


SRA_search <- read_tsv('~/Downloads/SraRunTable.txt')
search_by_study <- SRA_search %>% group_by(SRA_Study) %>% summarise(n_samples=length(SRA_Sample), n_run=length(Run))
full_xml <- read_xml('~/Downloads/SraExperimentPackage.xml') %>% as_list()
k <- full_xml$EXPERIMENT_PACKAGE_SET
i <- lapply(k, function(x) unlist(x, recursive = T))
k <- lapply(i, names) %>% reduce(intersect)
res <-lapply(i, function(x) x[k]) %>%   do.call(rbind,.) %>% as.data.frame()

