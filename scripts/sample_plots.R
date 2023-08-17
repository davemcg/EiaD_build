svg(filename = '2023_counts.sample_count.01.svg', width = 16, height = 12)
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
  summarise(`Study Count` = length(unique(study_accession)),
            `Sample Count` = n()) %>% 
  mutate(Tissue = case_when(grepl("Trab",Tissue) ~ 'TM', grepl("EyeLid", Tissue) ~ "Eye Lid", TRUE ~ Tissue),
         Source = case_when(grepl("Primary", Source) ~ "P. Culture",
                            TRUE ~ Source)
         ) %>% 
  ggplot(aes(x ='', y=`Sample Count`)) + 
  ggh4x::facet_nested(Tissue+Source+Sub_Tissue+Age+Perturbation ~ ., scale = 'free', space=  'free', switch = 'y',
                      strip = ggh4x::strip_nested(size = "variable"))+ 
  geom_bar(stat = 'identity') +
  xlab('') + 
  coord_flip()  + 
  theme_minimal() + 
  theme(strip.text.y.left = element_text(angle = 0, size = 12), text = element_text(size = 16))  + 
  theme(strip.background = element_rect(colour="gray", fill="gray"),
        strip.switch.pad.wrap = margin(t = 0, r = 0, b = 0, l = 0)) 

dev.off()





svg(filename = '2023_counts.study_count.01.svg', width = 16, height = 12)
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
  summarise(`Study Count` = length(unique(study_accession)),
            `Sample Count` = n()) %>%
  mutate(Tissue = case_when(grepl("Trab",Tissue) ~ 'TM', grepl("EyeLid", Tissue) ~ "Eye Lid", TRUE ~ Tissue),
         Source = case_when(grepl("Primary", Source) ~ "P. Culture",
                            TRUE ~ Source)
  ) %>%
  ggplot(aes(x ='', y=`Study Count`)) +
  ggh4x::facet_nested(Tissue+Source+Sub_Tissue+Age+Perturbation ~ ., scale = 'free', space=  'free', switch = 'y',
                      strip = ggh4x::strip_nested(size = "variable"))+
  geom_bar(stat = 'identity') +
  xlab('') +
  coord_flip()  +
  theme_minimal() +
  theme(strip.text.y.left = element_text(angle = 0, size = 12), text = element_text(size = 16))  +
  theme(strip.background = element_rect(colour="gray", fill="gray"),
        strip.switch.pad.wrap = margin(t = 0, r = 0, b = 0, l = 0))

dev.off()

e17 %>% 
  filter(study_accession != 'SRP012682') %>% 
  group_by(Tissue, Sub_Tissue) %>% summarise(Count = n()) %>% ggplot(aes(x=Sub_Tissue, y = Count)) + 
  geom_bar(stat= 'identity') + 
  ggforce::facet_col(~Tissue, scales = 'free_y', space= 'free') + 
  cowplot::theme_cowplot() + 
  coord_flip()
