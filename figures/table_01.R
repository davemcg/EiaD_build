load('/Volumes/data/projects/nei/mcgaughey/auto_eyeIntegration/results/core_tight.Rdata')
core_tight_2019 <- core_tight
load('~/git/Human_eyeIntegration_App/www/core_tight.Rdata')
core_tight$Sub_Tissue <- gsub('_',' - ', core_tight$Sub_Tissue)
core_tight <- core_tight %>% select(-run_accession) %>% unique()
core_both <- bind_rows(core_tight_2019 %>% mutate(Version = '2019') %>% select(-run_accession) %>% filter(Kept == 'Kept'), core_tight %>% mutate(Version = '2017')) %>% unique()
core_tight %>% mutate(Tissue = factor(Tissue, levels = sort(unique(core_tight_2019$Tissue))), Sub_Tissue = factor(Sub_Tissue, levels = sort(unique(core_tight_2019$Sub_Tissue)))) %>% group_by(Tissue, Sub_Tissue) %>% summarise(Count=n()) %>% ggplot(aes(x=Sub_Tissue, y=Count, fill=Tissue)) + geom_bar(stat = 'identity') + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3))+ xlab('')
core_tight %>% filter(grepl('Retina|Cornea|Lens|RPE|ESC', Tissue)) %>%  group_by(Tissue, Sub_Tissue) %>% summarise(Count=n()) %>% ggplot(aes(x=Sub_Tissue, y=Count, fill=Tissue)) + geom_bar(stat = 'identity', position = 'dodge', color='grey') + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3))+ xlab('')
core_both %>% filter(grepl('Retina|Cornea|Lens|RPE|ESC', Tissue)) %>% unique() %>% group_by(Tissue, Sub_Tissue, Version) %>% summarise(Count=n()) %>% ggplot(aes(x=Sub_Tissue, y=Count, fill=Tissue, colour = Version)) + geom_bar(stat = 'identity', position = 'dodge') + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3)) + scale_colour_manual(values = c('grey','black')) + xlab('')
