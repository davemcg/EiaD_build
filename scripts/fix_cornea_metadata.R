cornea <- core_tight_2019 %>% filter(Tissue == 'Cornea') %>% 
  mutate(sub = case_when(grepl('Limbus', sample_attribute, ignore.case = T) ~ 'Limbus',
                         grepl('endo', sample_attribute, ignore.case = T) ~ 'Endothelium',
                         grepl('epi', sample_attribute, ignore.case = T) ~ 'Epithelium',
                         grepl('stroma', sample_attribute, ignore.case = T) ~ 'Stroma',
                         grepl('CEC', sample_attribute, ignore.case = T) ~ 'Endothelium',
                         TRUE ~ Origin),
         state = case_when(grepl('ES', sample_attribute) ~ 'Stem Cell',
                           grepl('imortalization: SV', sample_attribute) ~ 'Cell Line',
                           grepl('imortalization: te', sample_attribute) ~ 'Cell Line',
                           grepl('16-18wk', sample_attribute) ~ 'Fetal',
                           TRUE ~ ''))    %>% 
  mutate(nSub_Tissue = case_when(state == '' ~ sub,
                                 TRUE ~ paste0(state, ' ', sub))) %>% 
  select(sample_accession, sample_attribute, Tissue, Sub_Tissue, nSub_Tissue) %>% 
  mutate(tiss = paste0(Tissue, '_', gsub(' ','.',nSub_Tissue)))  %>% 
  select(sample_accession, tiss) %>% 
  data.frame() 
  
  write_tsv(left_join(meta, cornea, by = c("X1" = "sample_accession")) %>% mutate(X5 = case_when(is.na(tiss) ~ X5, TRUE ~ tiss)) %>% select(-tiss) %>% unique(), path = '~/git/eyeIntegration_data_build/sampleTable_2019-01-25_tissues.tab', col_names = F)
