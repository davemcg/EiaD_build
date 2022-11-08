emeta %>% 
  as_tibble() %>% 
  filter(study_accession != 'SRP012682') %>% 
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details),
         Perturbation = gsub("Adult Tissue","", Perturbation)) %>% 
  mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue), 
         Source = case_when(is.na(Source) ~ '', TRUE ~ Source), 
         Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
         Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation)) %>% 
  group_by(Tissue, Source, Sub_Tissue, Age, Perturbation) %>% 
  summarise(Count = n()) %>% 
  ggplot(aes(x ='', y=Count)) + 
  ggh4x::facet_nested(Tissue+Source+Sub_Tissue+Age+Perturbation ~ ., scale = 'free', space=  'free', switch = 'y') + 
  geom_bar(stat = 'identity') +
  xlab('') + 
  coord_flip()  + 
  theme_minimal() + 
  theme(strip.text.y.left = element_text(angle = 0))  + 
  theme(strip.background = element_rect(colour="gray", fill="gray"),
        strip.switch.pad.wrap = margin(t = 0, r = 0, b = 0, l = 0))



e17 %>% 
  filter(study_accession != 'SRP012682') %>% 
  group_by(Tissue, Sub_Tissue) %>% summarise(Count = n()) %>% ggplot(aes(x=Sub_Tissue, y = Count)) + 
  geom_bar(stat= 'identity') + 
  ggforce::facet_col(~Tissue, scales = 'free_y', space= 'free') + 
  cowplot::theme_cowplot() + 
  coord_flip()
