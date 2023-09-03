library(tidyverse)
svg(filename = '2023_counts.03.svg', width = 18, height = 12)
emeta <- data.table::fread("~/git/EiaD_build/data/eyeIntegration23_meta_2023_09_01.built.csv.gz")

a <- emeta %>%
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
  ggplot(aes(x ='', y=`Study Count`, label = `Study Count`)) +
  ggh4x::facet_nested(Tissue+Source+Sub_Tissue+Age+Perturbation ~ ., scale = 'free', space=  'free', switch = 'y',
                      strip = ggh4x::strip_nested(size = "variable"))+
  geom_bar(stat = 'identity', fill = "gray20") +
  xlab('') +
  coord_flip()  +
  theme_minimal() + 
  theme(strip.text.y.left = element_text(angle = 0, size = 12), text = element_text(size = 16))  + 
  theme(strip.background = element_rect(colour="gray20", fill="gray20"),
        strip.switch.pad.wrap = margin(t = 0, r = 0, b = 0, l = 0),
        strip.text.y = element_text(color="white")
  ) +
  geom_text(nudge_y = 0.5, size = 5)

b <- emeta %>%
  as_tibble() %>%
  filter(study_accession != 'SRP012682') %>%
  mutate(Perturbation = case_when(grepl('MGS', Source_details) ~ Source_details),
         Perturbation = gsub("Adult Tissue","", Perturbation)) %>%
  mutate(Sub_Tissue = case_when(is.na(Sub_Tissue) ~ '', TRUE ~ Sub_Tissue),
         Source = case_when(is.na(Source) ~ '', TRUE ~ Source),
         Age = case_when(is.na(Age) ~ '', TRUE ~ Age),
         Perturbation = case_when(is.na(Perturbation) ~ '', TRUE ~ Perturbation),
         Sex = Sex_ML) %>%
  group_by(Tissue, Source, Sub_Tissue, Age, Perturbation, Sex) %>%
  summarise(
            `Sample Count` = n()) %>%
  mutate(Tissue = case_when(grepl("Trab",Tissue) ~ 'TM', grepl("EyeLid", Tissue) ~ "Eye Lid", TRUE ~ Tissue),
         Source = case_when(grepl("Primary", Source) ~ "P. Culture",
                            TRUE ~ Source)
  ) %>%
  ggplot(aes(x ='', y=`Sample Count`, label = `Sample Count`, fill = Sex)) +
  ggh4x::facet_nested(Tissue+Source+Sub_Tissue+Age+Perturbation ~ ., scale = 'free', space=  'free', switch = 'y',
                      strip = ggh4x::strip_nested(size = "variable"))+
  geom_bar(stat = 'identity') +
  xlab('') +
  coord_flip()  +
  theme_minimal() +   
  theme(axis.text.y=element_blank(), 
        strip.text.y.left = element_blank(),
        axis.ticks.y=element_blank(), 
        text = element_text(size = 16)
  ) +
  geom_text(position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("violet","lightskyblue"))

cowplot::plot_grid(a,b,nrow=1, rel_widths = c(1,0.6))
dev.off()

system("cp 2023_counts.03.svg ~/git/eyeIntegration_app/inst/app/www/2023/")
